//! This module contains an abstraction for mapping triangles in a mesh to a 2D UV space.

use crate::Result;
use crate::geom2::Point2;
use parry2d_f64::query::PointQueryWithLocation;
use parry2d_f64::shape::TriMesh;

#[derive(Clone)]
pub struct UvMapping {
    tri_map: TriMesh,
}

impl UvMapping {
    pub fn new(vertices: Vec<Point2>, faces: Vec<[u32; 3]>) -> Result<Self> {
        let tri_map = TriMesh::new(vertices, faces)?;
        Ok(Self { tri_map })
    }

    pub fn faces(&self) -> &[[u32; 3]] {
        self.tri_map.indices()
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
        let tri = self.tri_map.triangle(tri_id as u32);
        let p = tri.a.coords * barycentric[0]
            + tri.b.coords * barycentric[1]
            + tri.c.coords * barycentric[2];
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
            .project_local_point_and_get_location(point, false);
        let (_, (t_id, loc)) = result;
        Some((t_id as usize, loc.barycentric_coordinates().unwrap()))
    }
}
