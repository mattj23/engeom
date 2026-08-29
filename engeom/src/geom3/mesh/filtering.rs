//! This module has implementations of different ways of filtering/reducing a mesh

use crate::common::points::{dist, mean_point, triangle_area};
use crate::common::{IndexMask, PCoords};
use crate::geom3::mesh::MeshData3;
use crate::geom3::mesh::algorithms::subsets::{compact_by_masks, compute_unique_point_mask};
use crate::geom3::mesh::nav_structure::MeshNav;
use crate::geom3::mesh::patches::PatchFilter;
use crate::{Mesh3, Point3, SelectOp, Selection, SurfacePoint3, UnitVec3, Vector3};
use crate::{Plane3, Result};
use parry3d_f64::query::PointQuery;
use std::collections::HashMap;
use std::f64::consts::PI;

pub struct TriangleFilter<'a> {
    mesh: &'a Mesh3,
    mask: IndexMask,
}

// Internal implementations
impl TriangleFilter<'_> {
    /// Get the indices of the triangles which would need to be checked for an operation of the
    /// specified type. If the operation is `SelectOp::Add`, then the triangles that are not in the
    /// current selection will be returned. If the operation is `SelectOp::Remove`, or
    /// `SelectOp::Keep` then the triangles that are in the current selection will be returned.
    fn to_check(&self, mode: SelectOp) -> IndexMask {
        match mode {
            // When adding, we want to check all faces that are not currently selected
            SelectOp::Add => {
                let mut check_mask = self.mask.clone();
                check_mask.not_mut();
                check_mask
            }

            // When removing or keeping, we want to check all faces that are currently selected
            SelectOp::Remove | SelectOp::KeepOnly => self.mask.clone(),
        }
    }

    fn mutate_pass_list(mut self, mode: SelectOp, pass_mask: &IndexMask) -> Self {
        match mode {
            SelectOp::Add => self.mask.or_mut(pass_mask).unwrap(),
            SelectOp::Remove => {
                let mut flipped = pass_mask.clone();
                flipped.not_mut();
                self.mask.and_mut(&flipped).unwrap();
            }
            SelectOp::KeepOnly => {
                self.mask.and_mut(pass_mask).unwrap();
            }
        };

        self
    }
}

// Public API
impl TriangleFilter<'_> {
    /// Collect the indices of the triangles that have been filtered
    pub fn collect_indices(self) -> Vec<usize> {
        self.mask.to_indices()
    }

    /// Take the mask of indices that have been filtered
    pub fn take_mask(self) -> IndexMask {
        self.mask
    }

    /// Create a new mesh from the filtered faces, carrying every attribute across.
    ///
    /// returns: `Result<Mesh3>`, failing if the selection is empty, since a mesh needs at least one
    /// face to build an acceleration structure over
    pub fn into_mesh(self) -> Result<Mesh3> {
        self.mesh.extract_subset_faces(&self.mask)
    }

    /// Perform a direct mask operation on the current selection. This will modify the currently
    /// selected faces based on the operation:
    ///
    /// - `SelectOp::Add`: Add the triangles in the mask to the current selection.
    /// - `SelectOp::Remove`: Remove the triangles in the mask from the current selection.
    /// - `SelectOp::KeepOnly`: Keep only the triangles which are in both the current selection
    ///   _and_ the mask, removing all others.
    ///
    /// # Arguments
    ///
    /// * `mask`: The mask of indices to apply the operation to. This mask should have the same
    ///   length as the number of faces in the mesh.
    /// * `op`: The operation to perform on the current selection
    ///
    /// returns: TriangleFilter
    pub fn by_mask(self, mask: &IndexMask, op: SelectOp) -> Result<Self> {
        let new_mask = match op {
            SelectOp::Add => self.mask.or(mask)?,
            SelectOp::Remove => {
                let mut new_mask = self.mask.not();
                new_mask.or_mut(mask)?;
                new_mask.not_mut();
                new_mask
            }
            SelectOp::KeepOnly => self.mask.and(mask)?,
        };

        Ok(TriangleFilter {
            mesh: self.mesh,
            mask: new_mask,
        })
    }

    /// Perform a selection operation with triangles whose vertices are within a certain distance
    /// of a test point. The selection can allow a triangle with _any_ vertex to be included, or
    /// it can require that _all_ vertices of the triangle are within the distance, depending on
    /// the value of `all_vertices`.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to check triangle vertices against
    /// * `max_dist`: the maximum allowable distance between the point and vertices
    /// * `all_vertices`: if `true`, all vertices of the triangle must be within the distance for
    ///   the triangle to be included; if `false`, only one vertex needs to be within the distance.
    /// * `mode`: the type of operation to perform with the triangles that meet the distance
    ///   criterial
    ///
    /// returns: TriangleFilter
    pub fn vertices_near_point(
        self,
        point: &impl PCoords<3>,
        max_dist: f64,
        all_vertices: bool,
        mode: SelectOp,
    ) -> Self {
        let check_mask = self.to_check(mode);
        let mut op_mask = IndexMask::new(self.mesh.faces().len(), false);

        for i in check_mask.iter_true() {
            let face = self.mesh.shape.triangle(i as u32);

            // Check if the triangle is above the plane
            if all_vertices {
                if dist(&face.a, point) <= max_dist
                    && dist(&face.b, point) <= max_dist
                    && dist(&face.c, point) <= max_dist
                {
                    op_mask.set(i, true);
                }
            } else if dist(&face.a, point) <= max_dist
                || dist(&face.b, point) <= max_dist
                || dist(&face.c, point) <= max_dist
            {
                op_mask.set(i, true);
            }
        }

        self.mutate_pass_list(mode, &op_mask)
    }

    /// Perform a selection operation based on the position of triangles relative to a plane. This
    /// function will check the position of each triangle's vertices against the plane and include
    /// it in the operation based on whether any vertex (`all_vertices=false`) or all vertices
    /// (`all_vertices=true`) lie in the positive half-space defined by the plane.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to check against
    /// * `all_vertices`: if `true`, all vertices of the triangle must be above the plane for the
    ///   triangle to be included in the operation; if `false`, only one vertex needs to be above
    ///   the plane.
    /// * `mode`: the type of operation to perform with valid triangles
    ///
    /// returns: TriangleFilter
    /// Select the triangles belonging to connected patches which pass a filter.
    ///
    /// Patch structure is computed among the faces this operation would consider, not the whole
    /// mesh, which is what the other steps do and what makes the builder compose. Under
    /// `SelectOp::KeepOnly` and `SelectOp::Remove` that is the current selection, so a preceding
    /// step which carves the mesh up decides what counts as connected here; under `SelectOp::Add`
    /// it is everything not currently selected.
    ///
    /// Note that this builds a `MeshNav` of its own, so a chain which uses it more than once pays
    /// for the adjacency each time.
    ///
    /// # Arguments
    ///
    /// * `filter`: which patches are worth keeping
    /// * `mode`: how the passing triangles combine with the current selection
    ///
    /// returns: `Result<TriangleFilter>`
    pub fn keep_patches(self, filter: &PatchFilter, mode: SelectOp) -> Result<Self> {
        let check_mask = self.to_check(mode);
        let nav = MeshNav::new(self.mesh);
        let op_mask = nav.patch_mask(filter, Some(&check_mask))?;

        Ok(self.mutate_pass_list(mode, &op_mask))
    }

    pub fn above_plane(self, plane: &Plane3, all_vertices: bool, mode: SelectOp) -> Self {
        let check_mask = self.to_check(mode);
        let mut op_mask = IndexMask::new(self.mesh.faces().len(), false);

        for i in check_mask.iter_true() {
            let face = self.mesh.shape.triangle(i as u32);

            // Check if the triangle is above the plane
            if all_vertices {
                if plane.point_is_positive(&face.a)
                    && plane.point_is_positive(&face.b)
                    && plane.point_is_positive(&face.c)
                {
                    op_mask.set(i, true);
                }
            } else if plane.point_is_positive(&face.a)
                || plane.point_is_positive(&face.b)
                || plane.point_is_positive(&face.c)
            {
                op_mask.set(i, true);
            }
        }

        self.mutate_pass_list(mode, &op_mask)
    }

    /// Select triangles that are facing a certain direction within a specified angle. This
    /// function will check the angle between the normal of each triangle and the specified normal
    /// vector. If the angle is less than the specified angle, the triangle will be included in the
    /// operation.
    ///
    /// The `mode` parameter will determine if the triangles are added to the current selection,
    /// removed from it, or if the selection is modified to retain only the triangles that meet
    /// the direction criteria.
    ///
    /// # Arguments
    ///
    /// * `normal`: the normal vector to check against. This does not need to be normalized.
    /// * `angle`: the angle in radians to check against. If the angle between the triangle's normal
    ///   and the specified normal is less than this angle, the triangle will be included in the
    ///   operation.
    /// * `mode`: what kind of operation is done with triangles that meet the directional criteria.
    ///
    /// returns: TriangleFilter
    pub fn facing(self, normal: &Vector3, angle: f64, mode: SelectOp) -> Self {
        let check_mask = self.to_check(mode);
        let mut op_mask = IndexMask::new(self.mesh.faces().len(), false);

        for i in check_mask.iter_true() {
            let n = self.mesh.shape.triangle(i as u32).normal();
            if let Some(nv) = n {
                if nv.angle(normal) < angle {
                    op_mask.set(i, true);
                } else {
                    op_mask.set(i, false);
                }
            } else {
                op_mask.set(i, false);
            }
        }

        self.mutate_pass_list(mode, &op_mask)
    }

    /// Modify the list of indices to only include triangles that are within a certain distance of
    /// their closest projection onto another mesh. The distance can require that all points of the
    /// triangle are within the tolerance, or just one.
    ///
    /// There are two additional optional tolerances that can be applied.
    ///
    /// 1. A planar tolerance, which checks the distance of the vertex projected onto the plane of
    ///    the reference mesh triangle and looks at how far it is from the projection point. This
    ///    is useful to filter out triangles that go past the edge of the reference mesh.
    /// 2. An angle tolerance, which checks the angle between the normal of the current triangle
    ///    and the normal of the reference triangle. This is useful to filter out triangles that
    ///    are not facing the same direction as the reference mesh.
    ///
    /// # Arguments
    ///
    /// * `other`:
    /// * `all_points`: all points of the triangle must be within the tolerances if this is set to
    ///   `true`, otherwise only one point must be within the tolerances.
    /// * `distance_tol`:
    /// * `planar_tol`:
    /// * `angle_tol`:
    ///
    /// returns: TriangleFilter
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn near_mesh(
        self,
        other: &Mesh3,
        all_points: bool,
        distance_tol: f64,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
        mode: SelectOp,
    ) -> Self {
        let mut check = MeshNearCheck::new(self.mesh, other, distance_tol, planar_tol, angle_tol);
        let to_check = self.to_check(mode);
        let mut passes = IndexMask::new(self.mesh.faces().len(), false);
        for i in to_check.iter_true() {
            let tri = self.mesh.faces()[i];
            let face = self.mesh.shape.triangle(i as u32);

            let keep = if all_points {
                check.near_check(tri[0], face.normal())
                    && check.near_check(tri[1], face.normal())
                    && check.near_check(tri[2], face.normal())
            } else {
                check.near_check(tri[0], face.normal())
                    || check.near_check(tri[1], face.normal())
                    || check.near_check(tri[2], face.normal())
            };

            passes.set(i, keep);
        }

        self.mutate_pass_list(mode, &passes)
    }

    /// Expand/dilate the selection of triangles based on shared vertices with the currently
    /// selected set of triangles. Though it's named `expand`, it can also be used to shrink the
    /// selection by using `SelectOp::Remove` or `SelectOp::KeepOnly`.
    ///
    /// The function will check the triangles that are currently selected and will mutate the
    /// selection based on triangles that share vertices with any of the currently selected
    /// triangles. If a triangle shares a vertex with any of the currently selected triangles,
    /// what happens next depends on the `mode`:
    ///
    /// - `SelectOp::Add`: The triangle will be added to the selection. This will expand the
    ///   selection at the border by a single row of triangles, similar to a dilation operation in
    ///   image processing.
    /// - `SelectOp::Remove`: The triangle will be removed from the selection. This will shrink the
    ///   selection at the border by a single row of triangles, similar to an erosion operation in
    ///   image processing.
    /// - `SelectOp::KeepOnly`: This is a no-op.
    ///
    /// An optional `exclude` mask can be provided to exclude certain triangles from being even
    /// considered for the operation.
    ///
    /// # Arguments
    ///
    /// * `exclude`: An optional mask of indices that should be excluded from the operation. If
    ///   `None`, all triangles will be considered for the operation, if it contains a mask, the
    ///   mask can be thought of as a region which the expansion is not allowed to enter or which
    ///   is not allowed to erode.
    /// * `mode`: The operation to perform on the current selection. This can be `SelectOp::Add`,
    ///   or `SelectOp::Remove`. The `SelectOp::KeepOnly` is a no-op and will not change the
    ///   selection.
    ///
    /// returns: Result<TriangleFilter, Box<dyn Error, Global>>
    pub fn expand(self, exclude: Option<&IndexMask>, mode: SelectOp) -> Result<Self> {
        // Get the mask of indices that we want to check
        let check_mask = if let Some(exclude) = exclude {
            let mut check_base = self.to_check(mode);
            check_base.and_not_mut(exclude)?;
            check_base
        } else {
            self.to_check(mode)
        };

        // Get a mask of the vertices that are currently considered selected
        let vert_mask = match mode {
            // If we're adding new faces, we'll start with the vertices that are part of triangles
            // that are currently selected by the filter
            SelectOp::Add => self.mesh.compute_unique_point_mask(&self.mask),

            // If we're removing or keeping faces, we start with the vertices that are part of
            // triangles that are NOT currently selected by the filter
            SelectOp::Remove | SelectOp::KeepOnly => {
                let mut flipped = self.mask.clone();
                flipped.not_mut();
                self.mesh.compute_unique_point_mask(&flipped)
            }
        }
        .expect("Failed to create point mask from face mask, was the face mask valid?");

        // Now we'll check the triangles in the check mask, and if they contain any of the vertices
        // in the vertex mask, we'll add them to the pass list
        let mut passes = IndexMask::new(self.mesh.faces().len(), false);
        for i in check_mask.iter_true() {
            let t = self.mesh.faces()[i];
            if vert_mask.get(t[0] as usize)
                || vert_mask.get(t[1] as usize)
                || vert_mask.get(t[2] as usize)
            {
                passes.set(i, true);
            }
        }

        Ok(self.mutate_pass_list(mode, &passes))
    }

    /// This is a shorthand for calling `expand` multiple times in a row. It will apply the `expand`
    /// operation `n` times using the specified `mode`.
    ///
    /// # Arguments
    ///
    /// * `n`: the number of times to call `expand`.
    /// * `exclude`: An optional mask of indices that should be excluded from the operation. If
    ///   `None`, all triangles will be considered for the operation, if it contains a mask, the
    ///   mask can be thought of as a region which the expansion is not allowed to enter or which
    ///   is not allowed to erode.
    /// * `mode`: the operation to perform on the current selection. This can be `SelectOp::Add`,
    ///   or `SelectOp::Remove`. `SelectOp::KeepOnly` is a no-op and will not change the
    ///   selection.
    ///
    /// returns: Result<TriangleFilter, Box<dyn Error, Global>>
    pub fn expand_n(self, n: usize, exclude: Option<&IndexMask>, mode: SelectOp) -> Result<Self> {
        let mut filter = self;
        for _ in 0..n {
            filter = filter.expand(exclude, mode)?;
        }
        Ok(filter)
    }

    pub fn faces_overlap(
        self,
        other: &Mesh3,
        angle_tol: f64,
        distance_tol: f64,
        mode: SelectOp,
    ) -> Self {
        // Project every vertex onto the other mesh
        let projected: Vec<Option<Point3>> = self
            .mesh
            .points()
            .iter()
            .map(|v| {
                other
                    .shape
                    .project_local_point_with_max_dist(v, false, distance_tol)
                    .map(|p| p.point)
            })
            .collect();

        let to_check = self.to_check(mode);
        let mut pass_mask = IndexMask::new(self.mesh.faces().len(), false);
        for i in to_check.iter_true() {
            let f = self.mesh.faces()[i];
            let Some(v0) = projected[f[0] as usize] else {
                continue;
            };
            let Some(v1) = projected[f[1] as usize] else {
                continue;
            };
            let Some(v2) = projected[f[2] as usize] else {
                continue;
            };

            let tri = self.mesh.tri_mesh().triangle(i as u32);
            let area_original = tri.area();
            if area_original < 1e-12 {
                continue;
            }

            let Some(face_normal) = tri.normal() else {
                continue;
            };

            // Check that the centroid falls on a triangle of the other mesh with a normal
            // facing the same direction
            let centroid = mean_point(&[tri.a, tri.b, tri.c]);
            let mp = other.surface_closest_to(&centroid);
            if mp.normal().angle(&face_normal) > PI * 0.45 {
                continue;
            }

            // Check that the angle to the centroid is within the angle tolerance.
            //
            // A centroid which already lies on the other surface has no meaningful direction to its
            // own projection: the difference is floating point residue, and normalizing it produces
            // a direction made of rounding noise which lands perpendicular to the normal about as
            // often as not. Zero separation is the strongest agreement there is, so it passes
            // without the test. This mirrors the guard in `Mesh3::measure_point_deviation`.
            let v_to_centroid = mp.point() - centroid;
            if v_to_centroid.norm() > 1e-6 {
                let a_to_centroid = face_normal.angle(&v_to_centroid.normalize());
                if a_to_centroid > angle_tol && a_to_centroid < (PI - angle_tol) {
                    continue;
                }
            }

            // What's the area of the triangle formed by the projected points?
            let area_proj = triangle_area(&v0, &v1, &v2);
            if area_proj < 1e-12 {
                continue;
            }

            if dist(&v0, &tri.a) > distance_tol
                || dist(&v1, &tri.b) > distance_tol
                || dist(&v2, &tri.c) > distance_tol
            {
                continue;
            }

            let e0 = v1 - v0;
            let e1 = v2 - v0;
            let n = e0.cross(&e1).normalize();
            if face_normal.angle(&n) > angle_tol {
                continue;
            }

            pass_mask.set(i, true);
        }

        self.mutate_pass_list(mode, &pass_mask)
    }
}

impl Mesh3 {
    /// Create a new mask with the same length as the number of faces in the mesh, initialized to
    /// the specified value.
    pub fn face_mask(&self, value: bool) -> IndexMask {
        IndexMask::new(self.faces().len(), value)
    }

    /// Create a new mask with the same length as the number of points in the mesh, initialized to
    /// the specified value.
    pub fn point_mask(&self, value: bool) -> IndexMask {
        IndexMask::new(self.points().len(), value)
    }

    /// Start an operation to filter the faces of the mesh. This function will return a filter
    /// handle that can be used to add or remove faces from the selection while maintaining
    /// an immutable reference to the mesh.
    ///
    /// The filter can be started with no faces selected (`Selection::None`), all faces selected
    /// (`Selection::All`), or a specific set of faces selected (`Selection::Indices(Vec<usize>)`).
    /// Each successive filter operation will modify the selection the selected indices.
    ///
    /// # Arguments
    ///
    /// * `start`: The initial selection of faces to start with, either `Selection::None`,
    ///   `Selection::All`, or `Selection::Indices(Vec<usize>)`
    ///
    /// returns: TriangleFilter
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn face_select(&self, start: Selection) -> TriangleFilter<'_> {
        let mask = match start {
            Selection::None => IndexMask::new(self.faces().len(), false),
            Selection::All => IndexMask::new(self.faces().len(), true),
            Selection::Indices(i) => IndexMask::try_from_indices(&i, self.faces().len())
                .expect("Invalid indices for face selection"),
            Selection::Mask(m) => m,
        };

        TriangleFilter { mesh: self, mask }
    }

    /// Extract points and faces from the mesh based on a mask of face indices. This is a step
    /// towards creating a new mesh, but can be used independently.  To directly construct a new
    /// mesh, use `extract_subset_faces` instead, which also carries the attributes across.
    ///
    /// # Arguments
    ///
    /// * `mask`: a mask of face indices that will be used to filter the points and faces. Must
    ///   have the same length as the number of faces in the mesh, or the function will return an
    ///   error.
    ///
    /// returns: `Result<(Vec<Point3>, Vec<[u32; 3]>)>`
    pub fn compute_points_and_faces_from_mask(
        &self,
        mask: &IndexMask,
    ) -> Result<(Vec<Point3>, Vec<[u32; 3]>)> {
        let point_mask = self.compute_unique_point_mask(mask)?;
        compact_by_masks(self.points(), self.faces(), &point_mask, mask)
    }

    /// Create a new mesh from a mask of face indices.
    ///
    /// Only points referenced by a surviving face are kept, so any point the selection orphans is
    /// dropped. The surviving points are renumbered and the faces re-indexed to match.
    ///
    /// **Every attribute is carried through**, in both domains, selected by the same masks the
    /// geometry was. The result is not solid regardless of what this mesh was, since a subset of a
    /// closed surface generally is not closed.
    ///
    /// # Arguments
    ///
    /// * `mask`: a mask of face indices to be part of the new mesh. Must have the same length as
    ///   the number of faces in the mesh, or the function will return an error.
    ///
    /// returns: `Result<Mesh3>`, failing if the mask is the wrong length or selects no faces
    pub fn extract_subset_faces(&self, mask: &IndexMask) -> Result<Self> {
        let point_mask = self.compute_unique_point_mask(mask)?;
        let (points, faces) = compact_by_masks(self.points(), self.faces(), &point_mask, mask)?;
        let point_attrs = self.point_attrs.subset(&point_mask)?;
        let face_attrs = self.face_attrs.subset(mask)?;

        Self::from_data(
            MeshData3::new_with_attrs(points, faces, point_attrs, face_attrs)?,
            false,
        )
    }

    /// Create a new mesh from a list of face indices. The indices correspond with elements in the
    /// `faces()` slice.
    ///
    /// This is `extract_subset_faces` with the selection given as indices rather than a mask, and
    /// behaves identically: orphaned points are dropped, the survivors are renumbered, and every
    /// attribute is carried through. Because the selection becomes a mask, the faces of the result
    /// are in ascending index order regardless of the order the indices were given in, and a
    /// repeated index selects its face once rather than duplicating it.
    ///
    /// # Arguments
    ///
    /// * `indices`: the indices of the faces to keep, each of which must be less than the face
    ///   count
    ///
    /// returns: `Result<Mesh3>`
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use engeom::{Mesh3, Vector3, SelectOp, Selection};
    /// let mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
    /// let indices = mesh.face_select(Selection::None)
    ///     .facing(&Vector3::z(), PI / 2.0, SelectOp::Add)
    ///     .collect_indices();
    /// let new_mesh = mesh.extract_subset_faces_from_indices(&indices).unwrap();
    ///
    /// assert_eq!(new_mesh.faces().len(), 2);
    /// assert_eq!(new_mesh.points().len(), 4);
    /// ```
    pub fn extract_subset_faces_from_indices(&self, indices: &[usize]) -> Result<Self> {
        let mask = IndexMask::try_from_indices(indices, self.faces().len())?;
        self.extract_subset_faces(&mask)
    }

    pub fn compute_unique_point_mask(&self, face_mask: &IndexMask) -> Result<IndexMask> {
        compute_unique_point_mask(self.faces(), face_mask, self.points().len())
    }
}

struct MeshNearCheck<'a> {
    this_mesh: &'a Mesh3,
    ref_mesh: &'a Mesh3,
    checked: HashMap<u32, bool>,
    distance_tol: f64,
    planar_tol: Option<f64>,
    angle_tol: Option<f64>,
}

impl<'a> MeshNearCheck<'a> {
    fn new(
        this_mesh: &'a Mesh3,
        ref_mesh: &'a Mesh3,
        distance_tol: f64,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
    ) -> Self {
        Self {
            this_mesh,
            ref_mesh,
            checked: HashMap::new(),
            distance_tol,
            planar_tol,
            angle_tol,
        }
    }

    fn store_and_return(&mut self, vertex_index: u32, result: bool) -> bool {
        self.checked.insert(vertex_index, result);
        result
    }

    fn near_check(&mut self, vertex_index: u32, face_normal: Option<UnitVec3>) -> bool {
        if let Some(&checked) = self.checked.get(&vertex_index) {
            checked
        } else {
            let p = self.this_mesh.points()[vertex_index as usize];

            let is_ok = if let Some((prj, ri, _loc)) =
                self.ref_mesh.project_with_max_dist(&p, self.distance_tol)
            {
                if self.planar_tol.is_none() && self.angle_tol.is_none() {
                    true
                } else if let Some(rn) = self.ref_mesh.shape.triangle(ri).normal() {
                    // We need to get the normal of the reference triangle
                    let rsp = SurfacePoint3::new(prj.point, rn);

                    let check_planar = if let Some(planar_tol) = self.planar_tol {
                        rsp.planar_distance(&p) <= planar_tol
                    } else {
                        true
                    };

                    let check_angle = if let Some(angle_tol) = self.angle_tol {
                        if let Some(face_normal) = face_normal {
                            face_normal.angle(&rn) <= angle_tol
                        } else {
                            // No face normal, so we can't check the angle, assume it's bad?
                            false
                        }
                    } else {
                        true
                    };

                    check_planar && check_angle
                } else {
                    false
                }
            } else {
                false
            };

            self.store_and_return(vertex_index, is_ok)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Iso3;
    use crate::SelectOp::Add;
    use std::f64::consts::PI;

    #[test]
    fn test_triangles_facing() -> Result<()> {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let selection = mesh
            .face_select(Selection::None)
            .facing(&Vector3::z(), PI / 2.0, Add);

        let new_mesh = selection.into_mesh()?;
        assert_eq!(new_mesh.faces().len(), 2);

        for t in new_mesh.tri_mesh().triangles() {
            let n = t.normal().unwrap();
            assert!(n.dot(&Vector3::z()) > 0.0);
        }

        Ok(())
    }

    /// A body with a small flyer parked well away from it.
    fn body_with_flyer() -> Mesh3 {
        let mut mesh = Mesh3::create_box(20.0, 20.0, 20.0, false);
        let mut flyer = Mesh3::create_box(0.5, 0.5, 0.5, false);
        flyer.transform_in_place(&Iso3::translation(300.0, 0.0, 0.0));
        mesh.append_in_place(&flyer).unwrap();
        mesh
    }

    #[test]
    fn keep_patches_selects_the_surviving_patches() -> Result<()> {
        let mesh = body_with_flyer();
        let kept = mesh
            .face_select(Selection::All)
            .keep_patches(&PatchFilter::keep_largest(), SelectOp::KeepOnly)?
            .into_mesh()?;

        assert_eq!(kept.faces().len(), 12);
        Ok(())
    }

    /// Patch structure is computed among the faces the mode would consider, so a step which cuts
    /// the mesh up first changes what counts as connected. Here only the +z faces of each box
    /// survive the facing step, which leaves them touching at neither an edge nor a vertex.
    #[test]
    fn keep_patches_runs_within_the_current_selection() -> Result<()> {
        let mesh = body_with_flyer();

        let facing_up = mesh
            .face_select(Selection::None)
            .facing(&Vector3::z(), PI / 4.0, Add);
        assert_eq!(facing_up.take_mask().count_true(), 4, "two faces per box");

        // Both boxes contribute an equal-area pair, so the rank cut keeps one pair, not one face.
        let kept = mesh
            .face_select(Selection::None)
            .facing(&Vector3::z(), PI / 4.0, Add)
            .keep_patches(&PatchFilter::keep_largest(), SelectOp::KeepOnly)?
            .take_mask();

        assert_eq!(kept.count_true(), 2);
        Ok(())
    }

    /// `Remove` takes the passing patches out of the selection rather than keeping them.
    #[test]
    fn keep_patches_can_remove_instead() -> Result<()> {
        let mesh = body_with_flyer();
        let remaining = mesh
            .face_select(Selection::All)
            .keep_patches(&PatchFilter::keep_largest(), SelectOp::Remove)?
            .take_mask();

        // The body passed the filter, so it is what gets removed, leaving the flyer.
        assert_eq!(remaining.count_true(), 12);
        for face in 0..12 {
            assert!(!remaining.get(face));
        }
        Ok(())
    }

    /// Two copies of the same mesh in the same place overlap completely. This used to select
    /// nothing, because the direction from a face centroid to its own projection is floating point
    /// residue, and normalizing it produced a direction perpendicular to the face normal.
    #[test]
    fn faces_overlap_accepts_a_coincident_copy() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, true);
        let same = Mesh3::create_box(2.0, 2.0, 2.0, true);

        let selected = mesh
            .face_select(Selection::None)
            .faces_overlap(&same, 0.1, 0.1, Add)
            .collect_indices();

        assert_eq!(selected.len(), mesh.faces().len());
    }

    /// The separated case still has to be rejected, so the guard above did not simply disable the
    /// direction test.
    #[test]
    fn faces_overlap_rejects_a_mesh_which_is_far_away() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, true);
        let mut apart = Mesh3::create_box(2.0, 2.0, 2.0, true);
        apart.transform_in_place(&Iso3::translation(50.0, 0.0, 0.0));

        let selected = mesh
            .face_select(Selection::None)
            .faces_overlap(&apart, 0.1, 0.1, Add)
            .collect_indices();

        assert!(selected.is_empty());
    }
}
