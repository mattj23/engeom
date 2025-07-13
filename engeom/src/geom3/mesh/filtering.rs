//! This module has implementations of different ways of filtering/reducing a mesh

use crate::Result;
use crate::common::{IndexMask, SelectOp, Selection};
use crate::{Mesh, Point3, SurfacePoint3, UnitVec3, Vector3};
use itertools::Itertools;
use std::collections::{HashMap, HashSet};

pub struct TriangleFilter<'a> {
    mesh: &'a Mesh,
    mask: IndexMask,
}

impl TriangleFilter<'_> {
    /// Collect the indices of the triangles that have been filtered
    pub fn collect_indices(self) -> Vec<usize> {
        self.mask.to_indices()
    }

    /// Take the mask of indices that have been filtered
    pub fn take_mask(self) -> IndexMask {
        self.mask
    }

    /// Create a new mesh from the filtered indices
    pub fn create_mesh(self) -> Mesh {
        self.mesh.create_from_mask(&self.mask).unwrap()
    }

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
            SelectOp::Remove | SelectOp::Keep => self.mask.clone(),
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
            SelectOp::Keep => {
                self.mask.and_mut(pass_mask).unwrap();
            }
        };

        self
    }

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
        other: &Mesh,
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

    pub fn expand(self, mode: SelectOp) -> Self {
        // Get the mask of indices that we want to check
        let check_mask = self.to_check(mode);

        // Get a mask of the vertices that are currently considered selected
        let vert_mask = match mode {
            // If we're adding new faces, we'll start with the vertices that are part of triangles
            // that are currently selected by the filter
            SelectOp::Add => self.mesh.unique_vertex_mask(&self.mask),

            // If we're removing or keeping faces, we start with the vertices that are part of
            // triangles that are NOT currently selected by the filter
            SelectOp::Remove | SelectOp::Keep => {
                let mut flipped = self.mask.clone();
                flipped.not_mut();
                self.mesh.unique_vertex_mask(&flipped)
            }
        }
        .expect("Failed to create vertex mask from face mask, was the face mask valid?");

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

        self.mutate_pass_list(mode, &passes)
    }

    pub fn expand_n(self, n: usize, mode: SelectOp) -> Self {
        let mut filter = self;
        for _ in 0..n {
            filter = filter.expand(mode);
        }
        filter
    }
}

impl Mesh {
    /// Create a new mask with the same length as the number of faces in the mesh, initialized to
    /// the specified value.
    pub fn new_face_mask(&self, value: bool) -> IndexMask {
        IndexMask::new(self.faces().len(), value)
    }

    pub fn new_vertex_mask(&self, value: bool) -> IndexMask {
        IndexMask::new(self.vertices().len(), value)
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
    pub fn face_select(&self, start: Selection) -> TriangleFilter {
        let mask = match start {
            Selection::None => IndexMask::new(self.faces().len(), false),
            Selection::All => IndexMask::new(self.faces().len(), true),
            Selection::Indices(i) => IndexMask::try_from_indices(&i, self.faces().len())
                .expect("Invalid indices for face selection"),
            Selection::Mask(m) => m,
        };

        TriangleFilter { mesh: self, mask }
    }

    pub fn create_from_mask(&self, mask: &IndexMask) -> Result<Self> {
        let vertex_mask = self.unique_vertex_mask(mask)?;

        // The map_back array will map the old vertex indices to the new ones
        let mut map_back = vec![u32::MAX; self.vertices().len()];
        let mut new_verts = Vec::new();

        for (new_i, old_i) in vertex_mask.iter_true().enumerate() {
            map_back[old_i] = new_i as u32;
            new_verts.push(self.vertices()[old_i]);
        }

        let mut new_faces = Vec::new();
        for i in mask.iter_true() {
            let t = self.faces()[i];
            new_faces.push([
                map_back[t[0] as usize],
                map_back[t[1] as usize],
                map_back[t[2] as usize],
            ]);
        }

        Ok(Self::new(new_verts, new_faces, false))
    }

    /// Create a new mesh from a list of triangle indices. The indices correspond with elements in
    /// the `triangles()` slice. This function will iterate through the triangle indices,
    /// taking the three vertices associated with each index and marking them for inclusion in the
    /// new mesh. Then it will recreate the triangles, remapping them to the new vertex indices.
    ///
    /// # Arguments
    ///
    /// * `indices`: A slice of usize values that correspond to the indices of the triangles in the
    ///   original mesh. There cannot be any duplicate indices, or the function will return a
    ///   non-manifold mesh.
    ///
    /// returns: Mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use engeom::{Mesh, Vector3, SelectOp, Selection};
    /// let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
    /// let indices = mesh.face_select(Selection::None)
    ///     .facing(&Vector3::z(), PI / 2.0, SelectOp::Add)
    ///     .collect_indices();
    /// let new_mesh = mesh.create_from_indices(&indices);
    ///
    /// assert_eq!(new_mesh.faces().len(), 2);
    /// assert_eq!(new_mesh.vertices().len(), 4);
    /// ```
    pub fn create_from_indices(&self, indices: &[usize]) -> Self {
        let to_keep = self.unique_vertices(indices);
        // The map_back array will map the old vertex indices to the new ones
        let map_back: HashMap<u32, u32> = to_keep
            .iter()
            .enumerate()
            .map(|(i, v)| (*v, i as u32))
            .collect();

        let vertices: Vec<Point3> = to_keep
            .iter()
            .map(|i| self.vertices()[*i as usize])
            .collect();

        let triangles = indices
            .iter()
            .map(|i| {
                let t = self.faces()[*i];
                [map_back[&t[0]], map_back[&t[1]], map_back[&t[2]]]
            })
            .collect_vec();

        Self::new(vertices, triangles, false)
    }

    fn face_mask_matches(&self, face_mask: &IndexMask) -> bool {
        face_mask.len() == self.faces().len()
    }

    fn check_face_mask(&self, face_mask: &IndexMask) -> Result<()> {
        if !self.face_mask_matches(face_mask) {
            Err("Face mask length does not match the number of faces in the mesh".into())
        } else {
            Ok(())
        }
    }

    /// Using a mask of face indices, this function will create a vertex mask that contains only
    /// the vertices that are used in the triangles specified by the face mask.
    ///
    /// # Arguments
    ///
    /// * `face_mask`: a mask of face indices that will be used to filter the vertices. Must have
    ///   the same length as the number of faces in the mesh, or the function will return an error.
    ///
    /// returns: Result<IndexMask, Box<dyn Error, Global>>
    fn unique_vertex_mask(&self, face_mask: &IndexMask) -> Result<IndexMask> {
        self.check_face_mask(face_mask)?;

        let mut vertex_mask = IndexMask::new(self.vertices().len(), false);
        for i in face_mask.iter_true() {
            let t = self.faces()[i];
            vertex_mask.set(t[0] as usize, true);
            vertex_mask.set(t[1] as usize, true);
            vertex_mask.set(t[2] as usize, true);
        }

        Ok(vertex_mask)
    }

    fn unique_vertices(&self, triangle_indices: &[usize]) -> Vec<u32> {
        let mut to_save = HashSet::new();
        for i in triangle_indices {
            let t = self.faces()[*i];
            to_save.insert(t[0]);
            to_save.insert(t[1]);
            to_save.insert(t[2]);
        }

        // Now we can sort them in order
        let mut keep_order = to_save.iter().copied().collect_vec();
        keep_order.sort_unstable();

        keep_order
    }
}

struct MeshNearCheck<'a> {
    this_mesh: &'a Mesh,
    ref_mesh: &'a Mesh,
    checked: HashMap<u32, bool>,
    distance_tol: f64,
    planar_tol: Option<f64>,
    angle_tol: Option<f64>,
}

impl<'a> MeshNearCheck<'a> {
    fn new(
        this_mesh: &'a Mesh,
        ref_mesh: &'a Mesh,
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
            let p = self.this_mesh.vertices()[vertex_index as usize];

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
    use crate::SelectOp::Add;
    use std::f64::consts::PI;

    #[test]
    fn test_triangles_facing() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let selection = mesh
            .face_select(Selection::None)
            .facing(&Vector3::z(), PI / 2.0, Add);

        let new_mesh = selection.create_mesh();
        assert_eq!(new_mesh.faces().len(), 2);

        for t in new_mesh.tri_mesh().triangles() {
            let n = t.normal().unwrap();
            assert!(n.dot(&Vector3::z()) > 0.0);
        }
    }
}
