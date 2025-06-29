//! This module has tools for simulating laser triangulation sensors, which work by emitting a
//! laser line into a scene and using a detector to measure the angle of the reflected laser line.
//! These are commonly used sensors in metrology because of the relatively high performance they
//! are capable of achieving at a disproportionately low cost compared to other technologies.

use crate::na::Translation3;
use crate::sensors::SimulatedPointSensor;
use crate::{Iso3, Mesh, Point3, PointCloud, PointCloudFeatures, UnitVec3, Vector3};
use parry3d_f64::query::{Ray, RayCast};
use std::f64::consts::PI;

/// Represents a laser line sensor where a line of rays are cast from a single point and must
/// intersect a surface and be witnessed by another point in the sensor.
#[derive(Debug, Clone)]
pub struct LaserLine {
    ray_origin: Point3,
    detect_origin: Point3,
    line_start: Point3,
    line_end: Point3,
    min_range: f64,
    max_range: f64,
    rays: usize,
    angle_limit: Option<f64>,
}

impl LaserLine {
    pub fn new(
        ray_origin: Point3,
        detect_origin: Point3,
        line_start: Point3,
        line_end: Point3,
        min_range: f64,
        max_range: f64,
        rays: usize,
        angle_limit: Option<f64>,
    ) -> crate::Result<Self> {
        if rays < 2 {
            return Err("Number of rays must be at least 2".into());
        }
        Ok(Self {
            ray_origin,
            detect_origin,
            line_start,
            line_end,
            min_range,
            max_range,
            rays,
            angle_limit,
        })
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

impl SimulatedPointSensor for LaserLine {
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>) {
        let v = self.line_end - self.line_start;
        let limit = self.angle_limit.unwrap_or(PI / 2.0);

        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..self.rays {
            let f = (i as f64 / (self.rays - 1) as f64).clamp(0.0, 1.0);
            let n = ((self.line_start + v * f) - self.ray_origin).normalize();

            // This is the emitted ray
            let ray = Ray::new(self.ray_origin, n);

            // Check if the emitted ray intersects with an obstruction. If it does, the ob_limit
            // will be less than the f64::MAX value
            let ob_limit = obstruction_limit(obstruction, &ray, iso);

            let range_limit = ob_limit.min(self.max_range);

            // Check if the emitted ray intersects with the target before the obstruction
            if let Some(ri) =
                target
                    .tri_mesh()
                    .cast_ray_and_get_normal(iso, &ray, range_limit, false)
            {
                // Check that we're at least at the minimum range
                if ri.time_of_impact < self.min_range {
                    continue;
                }

                // Check the normal to the emitted ray
                let n = ri.normal * -1.0;
                if ray.dir.angle(&n) > limit {
                    continue;
                }

                // Create the witness ray
                let impact = ray.point_at(ri.time_of_impact);
                let witness = Ray::new(self.detect_origin, &impact - self.detect_origin);

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

#[derive(Debug, Clone)]
pub struct PanningLaserLine {
    laser_line: LaserLine,
    pan_vector: Vector3,
    steps: usize,
}

impl PanningLaserLine {
    pub fn new(laser_line: LaserLine, pan_vector: Vector3, steps: usize) -> Self {
        Self {
            laser_line,
            pan_vector,
            steps,
        }
    }
}

impl SimulatedPointSensor for PanningLaserLine {
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>) {
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..self.steps {
            let shift = Translation3::from(-self.pan_vector * i as f64);
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
