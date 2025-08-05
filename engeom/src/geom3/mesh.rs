//! This module contains an abstraction for a mesh of triangles, represented by vertices and their
//! indices into the vertex list.  This abstraction is built around the `TriMesh` type from the
//! `parry3d` crate.

mod collisions;
mod conformal;
mod edges;
mod faces;
pub mod filtering;
pub mod half_edge;
mod measurement;
mod nav_structure;
mod outline;
mod queries;
pub mod sampling;
mod uv_mapping;

use crate::common::{IndexMask, PCoords};
use crate::geom3::{Aabb3, IsoExtensions3};
use crate::na::SVector;
use crate::{Iso3, Point2, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
pub use collisions::MeshCollisionSet;
pub use edges::MeshEdges;
pub use half_edge::HalfEdgeMesh;
pub use nav_structure::MeshNav;
use parry3d_f64::shape::{TriMesh, TriMeshFlags};
use parry3d_f64::{shape, transformation};
use parry3d_f64::query::PointQueryWithLocation;
pub use uv_mapping::UvMapping;

/// A struct which represents a point on the surface of a mesh, including the index of the face
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

// Core access
impl Mesh {
    /// Get a reference to the AABB of the underlying mesh in the local coordinate system.
    pub fn aabb(&self) -> &Aabb3 {
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

impl Mesh {
    pub fn calc_edges(&self) -> Result<MeshEdges> {
        MeshEdges::new(self)
    }

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

    /// Transform the mesh in place by applying the given transformation to all vertices.
    pub fn transform_by(&mut self, transform: &Iso3) {
        self.shape.transform_vertices(transform);
    }

    pub fn uv_to_3d(&self, uv: &Point2) -> Option<MeshSurfPoint> {
        let (i, bc) = self.uv()?.triangle(uv)?;
        if let Ok(msp) = self.at_barycentric(i, bc) {
            Some(msp)
        } else {
            None
        }
    }

    pub fn project_to_uv(&self, p: &impl PCoords<3>) -> Option<Point2> {
        let Some(uv_map) = self.uv() else {
            return None;
        };

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

    pub fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let sphere = shape::Ball::new(radius);
        let (vertices, faces) = sphere.to_trimesh(n_theta as u32, n_phi as u32);

        Self::new(vertices, faces, true)
    }

    pub fn create_box(length: f64, width: f64, height: f64, is_solid: bool) -> Self {
        let bx = shape::Cuboid::new(Vector3::new(length / 2.0, width / 2.0, height / 2.0));
        let (vertices, triangles) = bx.to_trimesh();
        Self::new(vertices, triangles, is_solid)
    }

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

    /// Create a new `MeshNav` structure for this mesh. This structure is used to efficiently
    /// navigate the mesh through edges and faces.  It is recommended to use this if you will be
    /// performing multiple structural queries on the mesh, so that the structure does not need to
    /// be recomputed each time.
    pub fn nav(&self) -> MeshNav {
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
        let edges = MeshEdges::new(&self)?;

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

    pub fn get_vertex_normals(&self) -> Vec<Vector3> {
        let mut sums: Vec<Vector3> = vec![Vector3::new(0.0, 0.0, 0.0); self.shape.vertices().len()];
        let mut counts = vec![0; self.shape.vertices().len()];

        for (indices, tri) in self.shape.indices().iter().zip(self.shape.triangles()) {
            if let Some(n) = tri.normal() {
                for i in indices {
                    sums[*i as usize] += n.into_inner();
                    counts[*i as usize] += 1;
                }
            }
        }

        // Normalize the normals
        for i in 0..sums.len() {
            if counts[i] > 0 {
                let v = sums[i] / counts[i] as f64;
                sums[i] = v.normalize();
            }
        }

        sums
    }
}
