//! This module contains an abstraction for a mesh of triangles, represented by vertices and their
//! indices into the vertex list.  This abstraction is built around the `TriMesh` type from the
//! `parry3d` crate.

mod patches;
mod serialization;
mod uv_mapping;

use crate::{Iso3, Plane3, Point2, Point3, SurfacePoint3, Result};
use std::f64::consts::PI;

pub use self::serialization::{MeshData, MeshFlatData};
pub use self::uv_mapping::UvMapping;
use crate::common::indices::index_vec;
use crate::common::SurfacePointCollection;
use crate::geom3::points::points_sample_poisson_disk;
use parry3d_f64::query::{PointProjection, PointQueryWithLocation, SplitResult};
use parry3d_f64::shape::{TriMesh, TrianglePointLocation};
use rand::prelude::SliceRandom;
use rand::thread_rng;

#[derive(Clone)]
pub struct Mesh {
    shape: TriMesh,
    is_solid: bool,
    uv: Option<UvMapping>,
}

impl Mesh {
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

    pub fn append(&mut self, other: &Mesh) -> Result<()> {
        // For now, both meshes must have an empty UV mapping
        if self.uv.is_some() || other.uv.is_some() {
            return Err("Cannot append meshes with UV mappings".into());
        }

        Ok(self.shape.append(&other.shape))
    }

    pub fn new_with_uv(
        vertices: Vec<Point3>,
        triangles: Vec<[u32; 3]>,
        is_solid: bool,
        uv: Option<UvMapping>,
    ) -> Self {
        let shape =
            TriMesh::new(vertices, triangles).expect("Failed to create TriMesh with UV mapping");
        Self {
            shape,
            is_solid,
            uv,
        }
    }

    pub fn uv(&self) -> Option<&UvMapping> {
        self.uv.as_ref()
    }

    pub fn vertices(&self) -> &[Point3] {
        self.shape.vertices()
    }

    pub fn triangles(&self) -> &[[u32; 3]] {
        self.shape.indices()
    }

    pub fn transform(&mut self, transform: &Iso3) {
        self.shape.transform_vertices(transform);
    }

    pub fn surf_closest_to(&self, point: &Point3) -> SurfacePoint3 {
        let result = self
            .shape
            .project_local_point_and_get_location(point, self.is_solid);
        let (projection, (tri_id, _location)) = result;
        let triangle = self.shape.triangle(tri_id);
        let normal = triangle.normal().unwrap(); // When could this fail? On a degenerate tri?
        SurfacePoint3::new(projection.point, normal)
    }

    pub fn uv_to_3d(&self, uv: &Point2) -> Option<SurfacePoint3> {
        let (i, bc) = self.uv()?.triangle(uv)?;
        let t = self.shape.triangle(i as u32);
        let coords = t.a.coords * bc[0] + t.b.coords * bc[1] + t.c.coords * bc[2];

        t.normal().map(|n| SurfacePoint3::new(coords.into(), n))
    }

    pub fn point_closest_to(&self, point: &Point3) -> Point3 {
        let (result, _) = self
            .shape
            .project_local_point_and_get_location(point, self.is_solid);
        result.point
    }

    pub fn project_with_max_dist(
        &self,
        point: &Point3,
        max_dist: f64,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        self.shape
            .project_local_point_and_get_location_with_max_dist(point, self.is_solid, max_dist)
            .map(|(prj, (id, loc))| (prj, id, loc))
    }

    /// Given a test point, return its projection onto the mesh *if and only if* it is within the
    /// given distance tolerance from the mesh and the angle between the normal of the triangle and
    /// the +/- vector from the triangle to the point is less than the given angle tolerance.
    ///
    /// When a test point projects onto to the face of a triangle, the vector from the triangle
    /// point to the test point will be parallel to the triangle normal, by definition.  The angle
    /// tolerance will come into effect when the test point projects to an edge or vertex.  This
    /// will happen occasionally when the test point is near an edge with two triangles that reflex
    /// away from the point, and it will happen when the test point is beyond the edge of the mesh.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to project onto the mesh
    /// * `max_dist`: the maximum search distance from the test point to find a projection
    /// * `max_angle`: the max allowable angle deviation between the mesh normal at the projection
    ///   and the vector from the projection to the test point
    /// * `transform`: an optional transform to apply to the test point before projecting it onto
    ///  the mesh
    ///
    /// returns: Option<(PointProjection, u32, TrianglePointLocation)>
    pub fn project_with_tol(
        &self,
        point: &Point3,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        let point = if let Some(transform) = transform {
            transform * point
        } else {
            *point
        };

        let result = self
            .shape
            .project_local_point_and_get_location_with_max_dist(&point, self.is_solid, max_dist);
        if let Some((prj, (id, loc))) = result {
            let local = point - &prj.point;
            let triangle = self.shape.triangle(id);
            if let Some(normal) = triangle.normal() {
                let angle = normal.angle(&local).abs();
                if angle < max_angle || angle > PI - max_angle {
                    Some((prj, id, loc))
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

    /// Return the indices of the points in the given list that project onto the mesh within the
    /// given distance tolerance and angle tolerance.  An optional transform can be provided to
    /// transform the points before projecting them onto the mesh.
    ///
    /// # Arguments
    ///
    /// * `points`:
    /// * `max_dist`:
    /// * `max_angle`:
    /// * `transform`:
    ///
    /// returns: Vec<usize, Global>
    pub fn indices_in_tol(
        &self,
        points: &[Point3],
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Vec<usize> {
        let mut result = Vec::new();
        for (i, point) in points.iter().enumerate() {
            if self
                .project_with_tol(point, max_dist, max_angle, transform)
                .is_some()
            {
                result.push(i);
            }
        }
        result
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
                    let uv = uv_map.point(id as usize, loc.barycentric_coordinates().unwrap());
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

    pub fn tri_mesh(&self) -> &TriMesh {
        &self.shape
    }

    pub fn is_solid(&self) -> bool {
        self.is_solid
    }

    pub fn create_box(width: f64, height: f64, depth: f64, is_solid: bool) -> Self {
        let (vertices, triangles) = box_geom(width, height, depth);
        Self::new(vertices, triangles, is_solid)
    }

    pub fn sample_uniform(&self, n: usize) -> Vec<SurfacePoint3> {
        let mut cumulative_areas = Vec::new();
        let mut total_area = 0.0;
        for tri in self.shape.triangles() {
            total_area += tri.area();
            cumulative_areas.push(total_area);
        }

        let mut result = Vec::new();
        for _ in 0..n {
            let r = rand::random::<f64>() * total_area;
            let tri_id = cumulative_areas
                .binary_search_by(|a| a.partial_cmp(&r).unwrap())
                .unwrap_or_else(|i| i);
            let tri = self.shape.triangle(tri_id as u32);
            let r1 = rand::random::<f64>();
            let r2 = rand::random::<f64>();
            let a = 1.0 - r1.sqrt();
            let b = r1.sqrt() * (1.0 - r2);
            let c = r1.sqrt() * r2;
            let v = tri.a.coords * a + tri.b.coords * b + tri.c.coords * c;
            result.push(SurfacePoint3::new(Point3::from(v), tri.normal().unwrap()));
        }

        result
    }

    pub fn sample_poisson(&self, radius: f64) -> Vec<SurfacePoint3> {
        let starting = self.sample_dense(radius * 0.5);
        // TODO: this can be more efficient without all the copying
        let points = starting.clone_points();
        let mut rng = thread_rng();
        let mut indices = index_vec(None, starting.len());
        indices.shuffle(&mut rng);

        let to_take = points_sample_poisson_disk(&points, radius, &indices);
        to_take.into_iter().map(|i| starting[i]).collect()
    }

    pub fn sample_dense(&self, max_spacing: f64) -> Vec<SurfacePoint3> {
        let mut sampled = Vec::new();
        for face in self.shape.triangles() {
            // Find the angle closest to 90 degrees
            let ua = &face.b - &face.a;
            let va = &face.c - &face.a;

            let ub = &face.a - &face.b;
            let vb = &face.c - &face.b;

            let uc = &face.a - &face.c;
            let vc = &face.b - &face.c;

            let aa = ua.angle(&va).abs() - PI / 2.0;
            let ab = ub.angle(&vb).abs() - PI / 2.0;
            let ac = uc.angle(&vc).abs() - PI / 2.0;

            let (u, v, p) = if aa < ab && aa < ac {
                (ua, va, face.a)
            } else if ab < aa && ab < ac {
                (ub, vb, face.b)
            } else {
                (uc, vc, face.c)
            };

            let nu = u.norm() / max_spacing;
            let nv = v.norm() / max_spacing;

            for ui in 0..nu as usize {
                for vi in 0..nv as usize {
                    let uf = ui as f64 / nu;
                    let vf = vi as f64 / nv;
                    if uf + vf <= 1.0 {
                        let p = p + u * uf + v * vf;
                        let sp = SurfacePoint3::new(p, face.normal().unwrap());
                        sampled.push(sp);
                    }
                }
            }
        }

        sampled
    }

    pub fn get_patches(&self) -> Vec<Vec<usize>> {
        patches::compute_patch_indices(self)
    }

    /// Gets the boundary points of each patch in the mesh.  This function will return a list of
    /// lists of points, where each list of points is the boundary of a patch.  Note that this
    /// function will not work on non-manifold meshes.
    ///
    /// returns: Vec<Vec<usize, Global>, Global>
    pub fn get_patch_boundary_points(&self) -> Vec<Vec<Point3>> {
        let patches = self.get_patches();
        patches
            .iter()
            .map(|patch| patches::compute_boundary_points(self, patch))
            .flatten()
            .collect()
    }

    pub fn split(&self, plane: &Plane3) -> SplitResult<Mesh> {
        let result = self.shape.local_split(&plane.normal, plane.d, 1.0e-6);
        match result {
            SplitResult::Pair(a, b) => {
                let mesh_a = Mesh::new_take_trimesh(a, false);
                let mesh_b = Mesh::new_take_trimesh(b, false);
                SplitResult::Pair(mesh_a, mesh_b)
            }
            SplitResult::Negative => SplitResult::Negative,
            SplitResult::Positive => SplitResult::Positive,
        }
    }
}

fn box_geom(width: f64, height: f64, depth: f64) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let vertices = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(width, 0.0, 0.0),
        Point3::new(0.0, 0.0, depth),
        Point3::new(width, 0.0, depth),
        Point3::new(0.0, height, 0.0),
        Point3::new(width, height, 0.0),
        Point3::new(0.0, height, depth),
        Point3::new(width, height, depth),
    ];

    let triangles = vec![
        [4, 7, 5],
        [4, 6, 7],
        [0, 2, 4],
        [2, 6, 4],
        [0, 1, 2],
        [1, 3, 2],
        [1, 5, 7],
        [1, 7, 3],
        [2, 3, 7],
        [2, 7, 6],
        [0, 4, 1],
        [1, 4, 5],
    ];

    (vertices, triangles)
}
