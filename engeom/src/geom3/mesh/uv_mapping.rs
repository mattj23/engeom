//! This module contains an abstraction for mapping triangles in a mesh to a 2D UV space.

use crate::common::points::barycentric;
use crate::raster2::RasterMapping;
use crate::{Point2, Result, Vector2};
use parry2d_f64::bounding_volume::SimdAabb;
use parry2d_f64::math::{Point, SIMD_WIDTH, SimdReal};
use parry2d_f64::na::{SimdBool, SimdValue};
use parry2d_f64::partitioning::{SimdVisitStatus, SimdVisitor};
use parry2d_f64::query::PointQuery;
use parry2d_f64::shape::TriMesh;

/// A `UvMapping` is a structure that represents a two-way mapping between a two-dimensional
/// space (typically referred to as UV to substitute for x and y) and the surface of a 3D mesh.
/// This nomenclature is common in computer graphics, where UV coordinates are used to represent
/// positions on a texture to be applied to the faces of a 3D polygonal model.
///
/// Unlike in normal texture applications, this implementation not only allows for points on the
/// mesh faces to be mapped to the UV plane, but also allows for any arbitrary point in the UV plane
/// to be mapped back to the corresponding triangle in the mesh and its barycentric coordinates.
///
/// Low distortion UV mappings have a number of applications in engineering, as they allow for a
/// transformation between a manifold surface and a 2D plane, where certain types of measurements,
/// transformations, and analyses can be performed more easily.
///
/// They can be chained with a second two-way mapping, transforming between the UV space and a
/// raster space, which can allow certain types of image-processing operations to be performed
/// like filtering, smoothing, convolutions, and as inputs to convolutional neural networks.
#[derive(Clone)]
pub struct UvMapping {
    tri_map: TriMesh,
}

impl UvMapping {
    /// Create a new UV mapping from a set of vertices and faces. If this mapping is used to
    /// provide two-way transformation between a 2D UV space and a 3D mesh, the vertices should
    /// be in the same order as their untransformed counterparts in the 3D mesh, and the list of
    /// faces should be copied exactly from the 3D mesh.
    ///
    /// This mapping uses an internal acceleration structure to allow fast lookups of points in the
    /// UV space and find their corresponding triangle and barycentric coordinates. That structure
    /// maintains its own list of faces, so the faces should be cloned from the original face.
    ///
    /// # Arguments
    ///
    /// * `vertices`: a list of points in the UV space, which correspond with the vertices of the
    ///   3D mesh conformally mapped to the 2D plane.
    /// * `faces`: a list of triangles in the UV space, which should be identical to the values
    ///   and order of the faces in the 3D mesh. Each triangle is represented by a list of three
    ///   vertex indices, which are indices into the `vertices` list.  The face at index `i` in the
    ///   UV mapping should be the same as the face at index `i` in the 3D mesh.
    ///
    /// returns: Result<UvMapping, Box<dyn Error, Global>>
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
    pub fn point(&self, tri_id: u32, barycentric: [f64; 3]) -> Point2 {
        let tri = self.tri_map.triangle(tri_id);
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
    /// returns: Option<(u32, [f64; 3])>
    pub fn triangle(&self, point: &Point2) -> Option<(u32, [f64; 3])> {
        // self.tri_map.qbvh().traverse_depth_first()

        let mut visitor = TriangleVisitor::new(&self.tri_map, *point);
        let _ = self.tri_map.qbvh().traverse_depth_first(&mut visitor);
        visitor.result
    }

    /// Creates a raster mapping for the UV space based on the triangle mesh. The raster mapping
    /// will have its origin at the minimum point of the UV space, and will be padded by the
    /// specified `padding` in pixels. The size of each pixel in the raster mapping is defined by
    /// `px_size`.
    ///
    /// # Arguments
    ///
    /// * `px_size`: the size of each pixel in world units
    /// * `padding`: the number of pixels to pad the raster mapping on each side
    ///
    /// returns: RasterMapping
    pub fn make_raster_mapping(&self, px_size: f64, padding: usize) -> RasterMapping {
        let origin: Point2 = self.tri_map.local_aabb().mins.clone();
        let max = self.tri_map.local_aabb().maxs.clone();
        let size = max - origin;

        let width_px = (size.x / px_size).ceil() as usize + padding * 2;
        let height_px = (size.y / px_size).ceil() as usize + padding * 2;
        let pad = padding as f64 * px_size;
        let padded_origin = origin - Vector2::new(pad, pad);
        RasterMapping::new(padded_origin, (height_px, width_px), px_size, None)
    }
}

struct TriangleVisitor<'a> {
    tri_mesh: &'a TriMesh,
    point: Point2,
    result: Option<(u32, [f64; 3])>,
}

impl<'a> TriangleVisitor<'a> {
    pub fn new(tri_mesh: &'a TriMesh, point: Point2) -> Self {
        Self {
            tri_mesh,
            point,
            result: None,
        }
    }
}

impl SimdVisitor<u32, SimdAabb> for TriangleVisitor<'_> {
    fn visit(
        &mut self,
        bv: &SimdAabb,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus {
        let simd_point: Point<SimdReal> = Point::splat(self.point);
        let mask = bv.contains_local_point(&simd_point);

        if let Some(data) = data {
            let bitmask = mask.bitmask();

            for (ii, data) in data.into_iter().enumerate() {
                if (bitmask & (1 << ii)) != 0 {
                    let Some(index) = data else { continue };
                    let tri = self.tri_mesh.triangle(*index);
                    if tri.contains_local_point(&self.point) {
                        let bc = barycentric(&tri.a, &tri.b, &tri.c, &self.point);
                        self.result = Some((*index, bc));
                        return SimdVisitStatus::ExitEarly;
                    }
                }
            }
        }

        SimdVisitStatus::MaybeContinue(mask)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::mean_point;
    use approx::assert_relative_eq;

    #[test]
    fn uv_mapping_bc_works() {
        let vertices = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];
        let faces = vec![[0, 1, 2]];

        let test_point = mean_point(&vertices);

        let uv_mapping = UvMapping::new(vertices, faces).unwrap();
        let (tri_id, bary) = uv_mapping.triangle(&test_point).unwrap();
        assert_eq!(tri_id, 0);
        assert_relative_eq!(bary[0], 1.0 / 3.0, epsilon = 1e-6);
        assert_relative_eq!(bary[1], 1.0 / 3.0, epsilon = 1e-6);
        assert_relative_eq!(bary[2], 1.0 / 3.0, epsilon = 1e-6);
    }

    #[test]
    fn uv_mapping_miss() {
        let vertices = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];
        let faces = vec![[0, 1, 2]];
        let test_point = Point2::new(2.0, 2.0); // Outside the triangle
        let uv_mapping = UvMapping::new(vertices, faces).unwrap();
        assert!(uv_mapping.triangle(&test_point).is_none());
    }
}
