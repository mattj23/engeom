//! This module contains an abstraction for a mesh of triangles, represented by vertices and their
//! indices into the vertex list.  This abstraction is built around the `TriMesh` type from the
//! `parry3d` crate.

mod collisions;
mod conformal;
mod edges;
pub mod filtering;
pub mod half_edge;
mod measurement;
mod nav_structure;
mod outline;
mod queries;
pub mod sampling;
mod uv_mapping;

use crate::common::{IndexMask, PCoords};
use crate::geom3::IsoExtensions3;
use crate::io::{deflate_bytes, u_bytes_to_mesh};
use crate::na::SVector;
use crate::{Iso3, Point2, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
pub use collisions::MeshCollisionSet;
pub use edges::MeshEdges;
pub use half_edge::HalfEdgeMesh;
pub use nav_structure::MeshNav;
use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::shape::{TriMesh, TriMeshFlags};
use parry3d_f64::{shape, transformation};
pub use uv_mapping::UvMapping;

/// A struct that represents a point on the surface of a mesh, including the index of the face
/// on which it lies, its barycentric coordinates, and the point/normal representation in space.
/// This representation has no link back to the original mesh, so the face index and barycentric
/// coordinates will be invalid if (1) the mesh is modified, or (2) if you attempt to use them on
/// a different mesh.
#[derive(Debug, Clone, Copy)]
pub struct MeshSurfPoint {
    /// The index of the face on which this point lies.
    pub face_index: u32,

    /// The barycentric coordinates of the point on the face.
    pub bc: [f64; 3],

    /// The surface point (point + normal) corresponding to this barycentric coordinate.
    pub sp: SurfacePoint3,
}

impl MeshSurfPoint {
    /// Create a new `MeshSurfPoint` from the given face index, barycentric coordinates, and
    /// surface point.
    pub fn new(face_index: u32, bc: [f64; 3], sp: SurfacePoint3) -> Self {
        Self { face_index, bc, sp }
    }

    /// Get the point in space corresponding to this surface point.
    pub fn point(&self) -> Point3 {
        self.sp.point
    }

    /// Get the normal at this surface point.
    pub fn normal(&self) -> UnitVec3 {
        self.sp.normal
    }

    pub fn transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            face_index: self.face_index,
            bc: self.bc,
            sp: iso * self.sp,
        }
    }
}

impl Default for MeshSurfPoint {
    fn default() -> Self {
        Self {
            face_index: 0,
            bc: [0.0, 0.0, 0.0],
            sp: SurfacePoint3::default(),
        }
    }
}

impl PCoords<3> for MeshSurfPoint {
    fn coords(&self) -> SVector<f64, 3> {
        self.sp.point.coords
    }
}

/// This is a triangle mesh optimized for collision detection and geometric queries. It is built on
/// top of the `parry3d` library's `TriMesh` type, which provides efficient storage and querying of
/// triangle meshes. This mesh has some basic functionality for interrogating its structure, and
/// some very basic functionality for editing.  However, it is not a structure optimized for
/// editing or modification.
#[derive(Clone)]
pub struct Mesh {
    shape: TriMesh,
    is_solid: bool,
    uv: Option<UvMapping>,
}

// ===============================================================================================
// Core access
// ===============================================================================================

impl Mesh {
    /// Get a reference to the AABB of the underlying mesh in the local coordinate system.
    pub fn aabb(&self) -> Aabb {
        self.shape.local_aabb()
    }

    /// Gets a reference to the underlying `TriMesh` object to provide direct access to
    /// the `parry3d` API.
    pub fn tri_mesh(&self) -> &TriMesh {
        &self.shape
    }

    /// Return a flag indicating whether the mesh is considered "solid" or not for the purposes of
    /// distance queries. If a mesh is "solid", then distance queries for points on the inside of
    /// the mesh will return a zero distance.
    pub fn is_solid(&self) -> bool {
        self.is_solid
    }

    /// Get a reference to the vertices of the mesh.
    pub fn vertices(&self) -> &[Point3] {
        self.shape.vertices()
    }

    /// Get a reference to the face indices of the mesh.
    pub fn faces(&self) -> &[[u32; 3]] {
        self.shape.indices()
    }
}

// ===============================================================================================
// General creation methods
// ===============================================================================================
impl Mesh {
    /// Create a new mesh from a list of vertices and a list of triangles.  Additional options can
    /// be set to merge duplicate vertices and delete degenerate triangles.
    ///
    /// # Arguments
    ///
    /// * `vertices`:
    /// * `triangles`:
    /// * `is_solid`:
    /// * `merge_duplicates`:
    /// * `delete_degenerate`:
    /// * `uv`:
    ///
    /// returns: Result<Mesh, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new_with_options(
        vertices: Vec<Point3>,
        triangles: Vec<[u32; 3]>,
        is_solid: bool,
        merge_duplicates: bool,
        delete_degenerate: bool,
        uv: Option<Vec<Point2>>,
    ) -> Result<Self> {
        let mut flags = TriMeshFlags::empty();
        if merge_duplicates {
            flags |= TriMeshFlags::MERGE_DUPLICATE_VERTICES;
            flags |= TriMeshFlags::DELETE_DUPLICATE_TRIANGLES;
        }
        if delete_degenerate {
            flags |= TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES;
            flags |= TriMeshFlags::DELETE_DEGENERATE_TRIANGLES;
        }

        let uv_mapping = if let Some(uv) = uv {
            Some(UvMapping::new(uv, triangles.clone())?)
        } else {
            None
        };

        let shape = TriMesh::with_flags(vertices, triangles, flags)?;
        Ok(Self {
            shape,
            is_solid,
            uv: uv_mapping,
        })
    }

    pub fn new(vertices: Vec<Point3>, triangles: Vec<[u32; 3]>, is_solid: bool) -> Self {
        let shape = TriMesh::new(vertices, triangles).expect("Failed to create TriMesh");
        Self {
            shape,
            is_solid,
            uv: None,
        }
    }
    pub fn new_take_trimesh(shape: TriMesh, is_solid: bool) -> Self {
        Self {
            shape,
            is_solid,
            uv: None,
        }
    }
}

// ===============================================================================================
// Mutation/Transformation
// ===============================================================================================
impl Mesh {
    /// Transform the mesh in place by applying the given transformation to all vertices.
    pub fn transform_by(&mut self, transform: &Iso3) {
        self.shape.transform_vertices(transform);
    }

    /// Create a new mesh by scaling all vertices uniformly.
    ///
    /// # Arguments
    ///
    /// * `scale`: a scale factor to apply to all vertices
    ///
    /// returns: Mesh
    pub fn new_scaled_uniform(&self, scale: f64) -> Self {
        let new_shape = self
            .shape
            .clone()
            .scaled(&Vector3::new(scale, scale, scale));
        Mesh::new_take_trimesh(new_shape, self.is_solid)
    }

    /// Create a new mesh by offsetting each vertex along its smoothed vertex normal.
    ///
    /// The offset is applied as `vertex + offset * normal`, where the normal is the
    /// normalized per-vertex normal computed from adjacent face normals.
    ///
    /// Positive offsets expand the mesh outward; negative offsets shrink it inward.
    /// The original mesh is not modified.
    ///
    /// # Arguments
    ///
    /// * `offset`: The distance to offset each vertex along its normal.
    ///
    /// returns: Mesh
    pub fn new_offset_vertices(&self, offset: f64) -> Self {
        // These are already normalized
        let normals = self.get_vertex_normals();

        let updated = self
            .vertices()
            .iter()
            .zip(normals.iter())
            .map(|(v, n)| v + offset * n)
            .collect();

        Self::new(updated, self.faces().to_vec(), self.is_solid)
    }
}

// ===============================================================================================
// Unsorted
// ===============================================================================================

impl Mesh {
    pub fn calc_edges(&self) -> Result<MeshEdges<'_>> {
        MeshEdges::new(self)
    }

    /// Return a convex hull of the points in the mesh.
    pub fn convex_hull(&self) -> Self {
        let (vertices, faces) = transformation::convex_hull(self.shape.vertices());
        Self::new(vertices, faces, true)
    }

    pub fn append(&mut self, other: &Mesh) -> Result<()> {
        // For now, both meshes must have an empty UV mapping
        if self.uv.is_some() || other.uv.is_some() {
            return Err("Cannot append meshes with UV mappings".into());
        }

        self.shape.append(&other.shape);
        Ok(())
    }

    pub fn uv(&self) -> Option<&UvMapping> {
        self.uv.as_ref()
    }

    pub fn uv_to_3d(&self, uv: &Point2) -> Option<MeshSurfPoint> {
        let (i, bc) = self.uv()?.triangle(uv)?;
        self.at_barycentric(i, bc).ok()
    }

    pub fn project_to_uv(&self, p: &impl PCoords<3>) -> Option<Point2> {
        let uv_map = self.uv()?;
        let mp = self.surf_closest_to(p);
        Some(uv_map.point(mp.face_index, mp.bc))
    }

    pub fn uv_with_tol(
        &self,
        point: &Point3,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(Point2, f64)> {
        if let Some(uv_map) = self.uv() {
            let point = if let Some(transform) = transform {
                transform * point
            } else {
                *point
            };

            if let Some((prj, id, loc)) = self.project_with_tol(&point, max_dist, max_angle, None) {
                let triangle = self.shape.triangle(id);
                if let Some(normal) = triangle.normal() {
                    let uv = uv_map.point(id, loc.barycentric_coordinates().unwrap());
                    // Now find the depth
                    let sp = SurfacePoint3::new(prj.point, normal);
                    Some((uv, sp.scalar_projection(&point)))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
    /// Create a new `MeshNav` structure for this mesh. This structure is used to efficiently
    /// navigate the mesh through edges and faces.  It is recommended to use this if you will be
    /// performing multiple structural queries on the mesh, so that the structure does not need to
    /// be recomputed each time.
    pub fn nav(&self) -> MeshNav<'_> {
        MeshNav::new(self)
    }

    /// Calculates the patches in the mesh. If you are going to be doing multiple queries of the
    /// structure of the mesh, either use the half-edge representation, or generate a `MeshNav`
    /// through the `nav()` method to avoid having to recompute the mesh structure each time.
    ///
    /// # Arguments
    ///
    /// * `mask`:
    ///
    /// returns: Result<Vec<IndexMask, Global>, Box<dyn Error, Global>>
    pub fn get_patches(&self, mask: Option<&IndexMask>) -> Result<Vec<IndexMask>> {
        let nav = self.nav();
        nav.patches(mask)
    }

    /// Gets the boundary points of each patch in the mesh.  This function will return a list of
    /// lists of points, where each list of points is the boundary of a patch.  Note that this
    /// function will not work on non-manifold meshes.
    ///
    /// returns: Result<Vec<Vec<usize, Global>, Global>>
    pub fn get_patch_boundary_points(&self) -> Result<Vec<Vec<Point3>>> {
        let edges = MeshEdges::new(self)?;

        let mut b_loops = Vec::new();
        for b_loop in edges.boundary_loops.iter() {
            b_loops.push(
                b_loop
                    .iter()
                    .map(|vi| self.vertices()[*vi as usize])
                    .collect(),
            );
        }

        Ok(b_loops)
    }

    pub fn get_face_normals(&self) -> Result<Vec<UnitVec3>> {
        let mut result = Vec::new();
        for t in self.shape.triangles() {
            if let Some(n) = t.normal() {
                result.push(n);
            } else {
                return Err("Failed to get normal".into());
            }
        }

        Ok(result)
    }

    /// Calculates and returns the areas of all triangular faces in the shape.
    ///
    /// This function iterates over all the triangles in the shape, computes the area
    /// of each triangle using the `area` method, and collects the results into a vector.
    pub fn get_face_areas(&self) -> Vec<f64> {
        self.shape.triangles().map(|t| t.area()).collect::<Vec<_>>()
    }

    /// Compute smooth per-vertex normals by averaging the normals of all adjacent triangles
    /// weighted by triangle area. At the end of the computation, the normals are normalized to
    /// have unit length.
    ///
    /// Be aware that vertices that are not referenced by any valid triangle keep the zero vector.
    ///
    /// Also, be aware that this may not produce the results you expect for meshes with large flat
    /// surfaces represented by multiple triangles. For example, on a cube mesh, not all corner
    /// vertices will point along the diagonals, since each vertex will have some faces where it
    /// touches two triangles and may have some faces where it touches only one triangle, making
    /// the weights uneven.
    pub fn get_vertex_normals(&self) -> Vec<Vector3> {
        let mut sums: Vec<Vector3> = vec![Vector3::new(0.0, 0.0, 0.0); self.shape.vertices().len()];
        let mut weights = vec![0.0; self.shape.vertices().len()];

        for (face_i, tri) in self.shape.triangles().enumerate() {
            let indices = self.shape.indices()[face_i];
            let a = tri.area();
            if let Some(n) = tri.normal() {
                for i in indices {
                    sums[i as usize] += n.into_inner() * a;
                    weights[i as usize] += a;
                }
            }
        }

        // Normalize the normals
        for i in 0..sums.len() {
            if weights[i] > 0.0 {
                let v = sums[i] / weights[i];
                sums[i] = v.normalize();
            }
        }

        sums
    }
}

// ===============================================================================================
// Shape creation methods
// ===============================================================================================

impl Mesh {
    pub fn create_cone(half_height: f64, radius: f64, steps: usize) -> Self {
        let cone = shape::Cone::new(half_height, radius);
        let (vertices, faces) = cone.to_trimesh(steps as u32);

        Self::new(vertices, faces, true)
    }

    pub fn create_capsule(
        p0: &Point3,
        p1: &Point3,
        radius: f64,
        n_theta: usize,
        n_phi: usize,
    ) -> Self {
        let capsule = shape::Capsule::new(*p0, *p1, radius);
        let (vertices, faces) = capsule.to_trimesh(n_theta as u32, n_phi as u32);

        Self::new(vertices, faces, true)
    }

    /// Create a spherical mesh centered at the origin. The `n_theta` and `n_phi` parameters control
    /// the tessellation density.
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius of the sphere.
    /// * `n_theta` - Number of subdivisions around the polar direction.
    /// * `n_phi` - Number of subdivisions around the azimuthal direction.
    ///
    /// returns: Mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Mesh;
    /// use approx::assert_relative_eq;
    ///
    /// let n_t = 14;
    /// let n_p = 15;
    /// let sphere = Mesh::create_sphere(1.0, n_t, n_p);
    ///
    /// assert_eq!(sphere.vertices().len(), n_t * (n_p - 1) + 2);
    ///
    /// // Verify that the vertices are on the surface of the unit sphere.
    /// for vertex in sphere.vertices() {
    ///     let dist_from_origin = vertex.coords.norm();
    ///     assert_relative_eq!(dist_from_origin, 1.0)
    /// }
    /// ```
    pub fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let sphere = shape::Ball::new(radius);
        let (vertices, faces) = sphere.to_trimesh(n_theta as u32, n_phi as u32);

        Self::new(vertices, faces, true)
    }

    /// Create a box mesh with the given dimensions, centered at the origin.
    ///
    /// # Arguments
    ///
    /// * `length`: the dimension of the box in the x direction
    /// * `width`: the dimension of the box in the y direction
    /// * `height`: the dimension of the box in the z direction
    /// * `is_solid`: whether the box is solid or hollow, used for some specific distance queries
    ///   in the underlying parry library
    ///
    /// returns: Mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Mesh;
    /// use approx::assert_relative_eq;
    /// let mesh = Mesh::create_box(2.0, 4.0, 6.0, false);
    /// assert_relative_eq!(mesh.aabb().maxs.x, 1.0);
    /// assert_relative_eq!(mesh.aabb().maxs.y, 2.0);
    /// assert_relative_eq!(mesh.aabb().maxs.z, 3.0);
    /// assert_relative_eq!(mesh.aabb().mins.x, -1.0);
    /// assert_relative_eq!(mesh.aabb().mins.y, -2.0);
    /// assert_relative_eq!(mesh.aabb().mins.z, -3.0);
    /// ```
    pub fn create_box(length: f64, width: f64, height: f64, is_solid: bool) -> Self {
        let bx = shape::Cuboid::new(Vector3::new(length / 2.0, width / 2.0, height / 2.0));
        let (vertices, triangles) = bx.to_trimesh();
        Self::new(vertices, triangles, is_solid)
    }

    /// Create a cylindrical mesh centered at the origin and aligned with the local `y` axis.
    /// The `radius` controls the cylinder radius, `height` its full height, and `steps`
    /// controls the tessellation density around the circumference.
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius of the cylinder.
    /// * `height` - Full height of the cylinder (along the y-axis).
    /// * `steps` - Number of subdivisions around the cylinder axis.
    ///
    /// returns: Mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Mesh, Point3};
    /// use approx::assert_relative_eq;
    ///
    /// let cyl = Mesh::create_cylinder(1.0, 4.0, 16);
    ///
    /// assert_relative_eq!(cyl.aabb().mins.z, -1.0);
    /// assert_relative_eq!(cyl.aabb().maxs.z,  1.0);
    /// assert_relative_eq!(cyl.aabb().mins.x, -1.0);
    /// assert_relative_eq!(cyl.aabb().maxs.x,  1.0);
    /// assert_relative_eq!(cyl.aabb().mins.y, -2.0);
    /// assert_relative_eq!(cyl.aabb().maxs.y,  2.0);
    ///
    /// for vertex in cyl.vertices() {
    ///     let proj = Point3::new(vertex.x, 0.0, vertex.z);
    ///     assert_relative_eq!(proj.coords.norm(), 1.0);
    /// }
    /// ```
    pub fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        let cyl = shape::Cylinder::new(height / 2.0, radius);
        let (vertices, faces) = cyl.to_trimesh(steps as u32);

        Self::new(vertices, faces, true)
    }

    pub fn create_rect_beam_between(
        p0: &Point3,
        p1: &Point3,
        width: f64,
        height: f64,
        up: &Vector3,
    ) -> Result<Self> {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;
        let box_geom = shape::Cuboid::new(Vector3::new(width / 2.0, height / 2.0, v.norm() / 2.0));

        // I think this is OK?
        let transform = Iso3::try_from_basis_zy(&v, up, Some(pc))?;

        let (vertices, faces) = box_geom.to_trimesh();
        let mut mesh = Self::new(vertices, faces, true);
        mesh.transform_by(&transform);
        Ok(mesh)
    }

    pub fn create_cylinder_between(p0: &Point3, p1: &Point3, radius: f64, steps: usize) -> Self {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;
        let cyl = shape::Cylinder::new(v.norm() / 2.0, radius);

        // I think this is OK?
        let transform = Iso3::try_from_basis_yz(&v, &Vector3::z(), Some(pc))
            .unwrap_or(Iso3::try_from_basis_yx(&v, &Vector3::x(), Some(pc)).unwrap());

        let (vertices, faces) = cyl.to_trimesh(steps as u32);
        let mut mesh = Self::new(vertices, faces, true);
        mesh.transform_by(&transform);
        mesh
    }

    /// Load a Stanford bunny mesh embedded in the binary with 453 vertices and 948 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res3.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res4() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_4.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 1889 vertices and 3851 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res3.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res3() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_3.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 8171 vertices and 16301 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res2.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res2() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_2.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::stanford_bun_4;
    use approx::assert_relative_eq;

    #[test]
    fn vertex_normals_match_vertex_count_and_are_normalized() {
        let mesh = stanford_bun_4();
        let normals = mesh.get_vertex_normals();

        assert_eq!(normals.len(), mesh.vertices().len());

        for normal in normals {
            assert_relative_eq!(normal.norm(), 1.0, epsilon = 1.0e-12);
        }
    }

    #[test]
    fn new_scaled_uniform_scales_spherical_radius() {
        let radius = 1.0;
        let scale = 2.5;
        let mesh = Mesh::create_sphere(radius, 100, 100);
        let scaled = mesh.new_scaled_uniform(scale);

        assert_eq!(mesh.vertices().len(), scaled.vertices().len());

        for vertex in scaled.vertices() {
            assert_relative_eq!(vertex.coords.norm(), radius * scale, epsilon = 1.0e-12);
        }
    }

    #[test]
    fn new_offset_vertices_preserves_spherical_radius() {
        let radius = 1.0;
        let offset = 0.1;
        let mesh = Mesh::create_sphere(radius, 100, 100);
        let offset_mesh = mesh.new_offset_vertices(offset);

        assert_eq!(mesh.vertices().len(), offset_mesh.vertices().len());

        for vertex in offset_mesh.vertices() {
            assert_relative_eq!(vertex.coords.norm(), radius + offset, epsilon = 1.0e-5);
        }
    }
}
