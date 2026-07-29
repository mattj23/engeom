//! This module contains features for taking measurements on meshes

use crate::common::DistMode;
use crate::metrology::Distance3;
use crate::{Mesh3, Point3, SurfacePoint3, UnitVec3};

impl Mesh3 {
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
        let closest = self.surface_closest_to(point).sp;

        // In both cases, the measurement point `b` will remain the test point and `a` will be the
        // where the reference point, what will change is the direction of the measurement

        let d = deviation_direction(&closest, point, dist_mode);

        Distance3::new(closest.point, *point, Some(d))
    }

    /// Compute the signed deviation of many points from this mesh at once, returning the scalar
    /// values rather than full measurements.
    ///
    /// Each value is what `measure_point_deviation` would report for the same point and mode, and
    /// the two share their sign convention. This exists because the scalar is what a bulk
    /// comparison actually wants, and building a `Distance3` per point to read one number off it
    /// costs two points and a direction each time.
    ///
    /// # Arguments
    ///
    /// * `points`: the test points to measure
    /// * `dist_mode`: whether to use the point-to-point distance or the scalar projection distance,
    ///   which differ near the edges of the mesh
    ///
    /// returns: `Vec<f64>`, one signed deviation per test point, positive on the outside
    pub fn measure_deviations(&self, points: &[Point3], dist_mode: DistMode) -> Vec<f64> {
        points
            .iter()
            .map(|point| {
                let closest = self.surface_closest_to(point).sp;
                let d = deviation_direction(&closest, point, dist_mode);
                d.dot(&(point - closest.point))
            })
            .collect()
    }
}

/// The direction a deviation measurement is taken along, given a test point and its projection.
///
/// A test point sitting on the surface has no meaningful direction to its own projection, so the
/// surface normal is used rather than normalizing what would be floating point residue.
fn deviation_direction(closest: &SurfacePoint3, point: &Point3, dist_mode: DistMode) -> UnitVec3 {
    match dist_mode {
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrology::Measurement;
    use approx::assert_relative_eq;

    /// The bulk and single-point paths share a direction helper, so they must agree exactly rather
    /// than merely closely.
    #[test]
    fn measure_deviations_matches_the_single_point_measurement() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, true);
        let points = vec![
            Point3::new(0.0, 0.0, 3.0),   // clearly outside
            Point3::new(0.0, 0.0, -3.0),  // clearly outside, other side
            Point3::new(0.0, 0.0, 1.0),   // exactly on the surface
            Point3::new(5.0, 5.0, 5.0),   // off the corner, where the modes differ
            Point3::new(0.1, -0.2, 0.05), // inside
        ];

        for mode in [DistMode::ToPoint, DistMode::ToPlane] {
            let bulk = mesh.measure_deviations(&points, mode);
            assert_eq!(bulk.len(), points.len());

            for (i, point) in points.iter().enumerate() {
                let single = mesh.measure_point_deviation(point, mode).value();
                assert_relative_eq!(bulk[i], single, epsilon = 1.0e-15);
            }
        }
    }

    #[test]
    fn measure_deviations_signs_the_outside_positive() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let values = mesh.measure_deviations(
            &[Point3::new(0.0, 0.0, 2.0), Point3::new(0.0, 0.0, 0.5)],
            DistMode::ToPoint,
        );

        assert_relative_eq!(values[0], 1.0, epsilon = 1.0e-12);
        assert!(values[1] < 0.0);
    }

    #[test]
    fn measure_deviations_of_nothing_is_nothing() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        assert!(mesh.measure_deviations(&[], DistMode::ToPoint).is_empty());
    }
}
