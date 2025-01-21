//! This module contains an abstraction for mapping triangles in a mesh to a 2D UV space.

use crate::geom2::Point2;
use crate::{Result, To3D};
use parry3d_f64::query::PointQueryWithLocation;
use parry3d_f64::shape::TriMesh;

#[derive(Clone)]
pub struct UvMapping {
    triangles: Vec<[Point2; 3]>,
    tri_map: TriMesh,
}

impl UvMapping {
    pub fn new(triangles: Vec<[Point2; 3]>) -> Result<Self> {
        let tri_map = tri_map_from_triangles(&triangles)?;
        Ok(Self { triangles, tri_map })
    }

    pub fn new_from_vertices(vertices: &[Point2]) -> Result<Self> {
        let mut triangles = Vec::new();
        for tri in vertices.chunks_exact(3) {
            triangles.push([tri[0], tri[1], tri[2]]);
        }
        Self::new(triangles)
    }

    pub fn to_vertices(&self) -> Vec<Point2> {
        let mut vertices = Vec::new();
        for tri in self.triangles() {
            vertices.push(tri[0]);
            vertices.push(tri[1]);
            vertices.push(tri[2]);
        }
        vertices
    }

    pub fn triangles(&self) -> &[[Point2; 3]] {
        &self.triangles
    }

    /// Given a triangle ID and a barycentric coordinate, return the corresponding point in the
    /// 2D UV space.
    ///
    /// # Arguments
    ///
    /// * `tri_id`: The ID of the triangle to map.
    /// * `barycentric`: The barycentric coordinate of the point to map on the triangle
    ///
    /// returns: OPoint<f64, Const<2>>
    pub fn point(&self, tri_id: usize, barycentric: [f64; 3]) -> Point2 {
        let tri = self.triangles()[tri_id];
        let p = tri[0].coords * barycentric[0]
            + tri[1].coords * barycentric[1]
            + tri[2].coords * barycentric[2];
        Point2::from(p)
    }

    /// Given a point in the UV space, return the corresponding triangle ID and barycentric
    /// coordinates of the closest point in the UV map.
    ///
    /// # Arguments
    ///
    /// * `point`: the point in UV space to test
    ///
    /// returns: Option<(usize, [f64; 3])>
    pub fn triangle(&self, point: &Point2) -> Option<(usize, [f64; 3])> {
        let result = self
            .tri_map
            .project_local_point_and_get_location(&point.to_3d(), false);
        let (_, (t_id, loc)) = result;
        Some((t_id as usize, loc.barycentric_coordinates().unwrap()))
    }
}

fn tri_map_from_triangles(tris: &[[Point2; 3]]) -> Result<TriMesh> {
    let mut vertices = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();

    for tri in tris {
        let i = vertices.len() as u32;
        vertices.push(tri[0].to_3d());
        vertices.push(tri[1].to_3d());
        vertices.push(tri[2].to_3d());

        triangles.push([i, i + 1, i + 2]);
    }

    TriMesh::new(vertices, triangles)
        .map_err(|e| format!("Failed to create TriMesh from triangles: {:?}", e).into())
}
