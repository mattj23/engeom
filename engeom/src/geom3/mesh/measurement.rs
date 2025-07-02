//! This module contains features for taking measurements on meshes

use crate::common::DistMode;
use crate::metrology::Distance3;
use crate::{Mesh, Point3, UnitVec3};

impl Mesh {
    /// Compute the deviation of a point from this mesh (this mesh is considered the reference) and
    /// return it as a Length Measurement object.
    ///
    /// The deviation is the distance from the point to its closest projection onto the mesh using
    /// the specified distance mode.  The direction of the measurement is the direction between the
    /// point and the projection, flipped into the positive half-space of the mesh surface at the
    /// projection point.
    ///
    /// If the distance is less than a very small floating point epsilon, the direction will be
    /// taken directly from the mesh surface normal.
    ///
    /// The first point `.a` of the measurement is the reference point, and the second point `.b`
    /// is the test point.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to measure the deviation from
    /// * `dist_mode`: whether to use the point-to-point distance or the scalar projection distance
    ///   when computing the deviation. This will have an effect near the edges of the mesh, in
    ///   which the `ToPlane` mode will not penalize a point for being off the mesh surface.
    ///
    /// returns: Length<3>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn measure_point_deviation(&self, point: &Point3, dist_mode: DistMode) -> Distance3 {
        let closest = self.surf_closest_to(point);

        // In both cases, the measurement point `b` will remain the test point and `a` will be the
        // where the reference point, what will change is the direction of the measurement

        let d = match dist_mode {
            DistMode::ToPoint => {
                let v = point - closest.point;
                if v.norm() < 1e-6 {
                    closest.normal
                } else if closest.normal.dot(&v) > 0.0 {
                    UnitVec3::new_normalize(v)
                } else {
                    -UnitVec3::new_normalize(v)
                }
            }
            DistMode::ToPlane => closest.normal,
        };

        Distance3::new(closest.point, *point, Some(d))
    }
}
