//! This module has tools for simulating laser profile triangulation sensors, which work by emitting
//! a laser line into a scene and using a detector to measure the angle of the reflected laser line.
//! These are commonly used sensors in metrology because of the relatively high performance they
//! are capable of achieving at a disproportionately low cost compared to other technologies.
//!
//! Examples include the LMI Gocator 2D profile series, the Micro-Epsilon scanCONTROL series,
//! and the Keyence LJ series.

use crate::na::Translation3;
use crate::sensors::SimulatedPointSensor;
use crate::{Iso3, Mesh, Point3, PointCloud, PointCloudFeatures, UnitVec3, Vector3};
use crate::{Result, SurfacePoint3};
use parry3d_f64::query::{Ray, RayCast};
use std::f64::consts::PI;

/// Represents the geometry of a laser profile line sensor, which emits a laser line into a scene
/// and detects the reflection of that line to triangulate the distance to points on a surface.
#[derive(Debug, Clone)]
pub struct LaserProfile {
    pub emitter_z: f64,
    pub detector_y: f64,
    pub detector_z: f64,
    pub volume_width: f64,
    pub volume_z_min: f64,
    pub volume_z_max: f64,
    pub resolution: usize,
    pub angle_limit: Option<f64>,
}

impl LaserProfile {
    /// Create the base geometry of a laser profile line sensor, which emits a laser line into a
    /// scene and detects the reflection of that line to triangulate the distance to points on a
    /// surface.
    ///
    /// The general coordinate system is specified in X and Z. The center of the detection volume
    /// is at the origin, with the laser line ranging from the -X direction to the +X direction.
    /// The +Z direction points directly up towards the emitter.  The +Y direction is orthogonal to
    /// laser line and is typically the direction which the sensor will be panned.
    ///
    /// The geometry is specified with the following assumptions:
    ///   - The laser line is emitted from a point directly on the +Z axis, with no offset in
    ///     the X or Y direction.
    ///   - The detector is not offset in the X direction, and can be specified with a Y and
    ///     Z offset from the center of the detection volume.
    ///   - The detection volume is trapezoidal, and its flat top and bottom are specified by a
    ///     maximum and minimum Z value.
    ///   - The detection volume's with is specified at Z=0, and is symmetrical around X=0.
    ///
    /// # Arguments
    ///
    /// * `emitter_z`: The Z coordinate of the laser emitter. This is the height from the volume
    ///   center where the laser fans into a triangle.
    /// * `detector_y`: The Y coordinate of the detector's optical center. This is the out-of-plane
    ///   offset from the plane of the laser line.
    /// * `detector_z`: The Z coordinate of the detector's optical center. This is the height from
    ///   the volume center where the detector's optical center is located.
    /// * `volume_width`: The width of the detection volume at Z=0. The volume is assumed to be
    ///   symmetrical around the X axis, ranging from -volume_width/2 to +volume_width/2.
    /// * `volume_z_min`: The minimum Z value of the detection volume. This is the bottom of the
    ///   trapezoidal volume, the farthest distance from the emitter where the sensor will still
    ///   return points.
    /// * `volume_z_max`: The maximum Z value of the detection volume. This is the top of the
    ///   trapezoidal volume, the closest distance to the emitter where the sensor will still
    ///   return points.
    /// * `resolution`: The number of rays to cast across the laser line. This is the number of
    ///   points that will be returned in the point cloud.
    /// * `angle_limit`: An optional angle limit in radians. If specified, the sensor will only
    ///   return a point if the angle between the surface normal at the point and the detector is
    ///   less than this limit.
    ///
    /// returns: Result<LaserProfileSensor, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new(
        emitter_z: f64,
        detector_y: f64,
        detector_z: f64,
        volume_width: f64,
        volume_z_min: f64,
        volume_z_max: f64,
        resolution: usize,
        angle_limit: Option<f64>,
    ) -> Self {
        Self {
            emitter_z,
            detector_y,
            detector_z,
            volume_width,
            volume_z_min,
            volume_z_max,
            resolution,
            angle_limit,
        }
    }
}

fn obstruction_limit(obstruction: Option<&Mesh>, ray: &Ray, iso: &Iso3) -> f64 {
    obstruction
        .map(|ob| {
            ob.tri_mesh()
                .cast_ray(iso, ray, f64::MAX, false)
                .unwrap_or(f64::MAX)
        })
        .unwrap_or(f64::MAX)
}

impl SimulatedPointSensor for LaserProfile {
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>) {
        let limit = self.angle_limit.unwrap_or(PI / 2.0);

        let mut points = Vec::new();
        let mut normals = Vec::new();
        let emitter = Point3::new(0.0, 0.0, self.emitter_z);
        let detector = Point3::new(0.0, self.detector_y, self.detector_z);
        let spacing = self.volume_width / (self.resolution - 1) as f64;

        for i in 0..self.resolution {
            let focus = Point3::new(i as f64 * spacing - self.volume_width / 2.0, 0.0, 0.0);
            let surf_point = SurfacePoint3::new_normalize(emitter, focus - emitter);
            let ray = Ray::from(&surf_point);

            // Calculate the min and max T values for the ray which correspond with the volume's
            // min and max Z values.
            let (min_t, max_t) = min_max_t(&ray, self.volume_z_min, self.volume_z_max);

            // Check if the emitted ray intersects with an obstruction. If it does, the ob_limit
            // will be less than the f64::MAX value
            let ob_limit = obstruction_limit(obstruction, &ray, iso);

            // The range limit is a way of accounting for both
            let range_limit = ob_limit.min(max_t);

            // Check if the emitted ray intersects with the target before the obstruction
            if let Some(ri) =
                target
                    .tri_mesh()
                    .cast_ray_and_get_normal(iso, &ray, range_limit, false)
            {
                // Check that we're at least at the minimum range
                if ri.time_of_impact < min_t {
                    continue;
                }

                // Check the normal to the emitted ray
                let n = ri.normal * -1.0;
                if ray.dir.angle(&n) > limit {
                    continue;
                }

                // Create the witness ray
                let impact = ray.point_at(ri.time_of_impact);
                let witness = Ray::new(emitter, &impact - detector);

                // Check if the witness ray intersects with the obstruction
                let ob_limit = obstruction_limit(obstruction, &witness, iso);
                if ob_limit < 1.0 - 1e-4 {
                    continue;
                }

                // Check if the witness ray intersects with the target before expected
                if let Some(_) = target.tri_mesh().cast_ray(iso, &witness, 1.0 - 1e-4, false) {
                    continue;
                }

                points.push(impact);
                normals.push(UnitVec3::new_normalize(ri.normal));
            }
        }

        // This should be safe because we assembled the points and normals together
        let cloud = PointCloud::try_new(points, Some(normals), None).unwrap();

        (cloud, None)
    }
}

/// Calculate the minimum and maximum time of flight values for a ray that intersects with the
/// volume min and max Z values.  Because the ray will be originating from the +Z direction and
/// pointing in -Z, the minimum time of flight will correspond with the volume's maximum Z value
/// and the maximum time of flight will correspond with the volume's minimum Z value.
///
/// # Arguments
///
/// * `ray`: The ray to calculate the time of flight limits for. MUST BE NORMALIZED.
/// * `volume_z_min`: The minimum Z value of the volume. This is the bottom of the trapezoidal
///   volume, farthest from the emitter.
/// * `volume_z_max`: The maximum Z value of the volume. This is the top of the trapezoidal
///   volume, closest to the emitter.
///
/// returns: (f64, f64)
fn min_max_t(ray: &Ray, volume_z_min: f64, volume_z_max: f64) -> (f64, f64) {
    let to_max = ray.origin.z - volume_z_max;
    let to_min = ray.origin.z - volume_z_min;
    let scale = 1.0 / -ray.dir.z;

    (scale * to_max, scale * to_min)
}

#[derive(Debug, Clone)]
pub struct PanningLaserProfile {
    laser_line: LaserProfile,
    y_step: f64,
    steps: usize,
}

impl PanningLaserProfile {
    pub fn new(laser_line: LaserProfile, y_step: f64, steps: usize) -> Self {
        Self {
            laser_line,
            y_step,
            steps,
        }
    }
}

impl SimulatedPointSensor for PanningLaserProfile {
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>) {
        let pan_vector = Vector3::new(0.0, self.y_step, 0.0);
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..self.steps {
            let shift = Translation3::from(-pan_vector * i as f64);
            let inv = shift.inverse();
            let t = shift * iso;
            let (cloud, _) = self.laser_line.get_points(target, obstruction, &t);
            points.extend(cloud.points().iter().map(|p| inv * p));
            normals.extend(cloud.normals().unwrap());
        }

        let cloud = PointCloud::try_new(points, Some(normals), None).unwrap();

        (cloud, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn time_of_flight_limits() {
        let z_min = -5.0;
        let z_max = 10.0;
        let emitter = Point3::new(0.0, 0.0, 20.0);
        let focus = Point3::new(13.0, 0.0, 0.0);
        let surf_point = SurfacePoint3::new_normalize(emitter, focus - emitter);
        let ray = Ray::from(&surf_point);

        let (min_t, max_t) = min_max_t(&ray, z_min, z_max);

        assert_relative_eq!(z_max, surf_point.at_distance(min_t).z, epsilon = 1e-8);
        assert_relative_eq!(z_min, surf_point.at_distance(max_t).z, epsilon = 1e-8);
    }
}
