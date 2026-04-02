//! Pinhole camera model for projecting 3D world points onto a 2D image plane and for
//! back-projecting image points into 3D rays.
//!
//! The camera coordinate system follows the standard computer-vision convention:
//!   - +X points right in the image
//!   - +Y points down in the image
//!   - +Z points out along the optical axis (into the scene)
//!
//! The `Iso3` supplied to every operation is the **camera-to-world** transform (i.e. applying
//! it to a point expressed in camera space yields the corresponding world-space point).  The
//! inverse is applied internally when needed to map world → camera coordinates.

use crate::common::PCoords;
use crate::{Iso3, Point2, Point3, UnitVec3, Vector3};
use parry3d_f64::query::Ray;

/// A pinhole camera defined by its intrinsic parameters.
///
/// # Intrinsic parameters
/// * `fx` / `fy` – focal lengths in pixels along the horizontal and vertical axes
/// * `cx` / `cy` – principal point (optical center) in pixel coordinates
#[derive(Debug, Clone)]
pub struct PinholeCamera {
    pub fx: f64,
    pub fy: f64,
    pub cx: f64,
    pub cy: f64,
}

impl PinholeCamera {
    /// Create a new pinhole camera from its intrinsic parameters.
    ///
    /// # Arguments
    /// * `fx` – horizontal focal length in pixels
    /// * `fy` – vertical focal length in pixels
    /// * `cx` – horizontal coordinate of the principal point in pixels
    /// * `cy` – vertical coordinate of the principal point in pixels
    pub fn new(fx: f64, fy: f64, cx: f64, cy: f64) -> Self {
        Self { fx, fy, cx, cy }
    }

    /// Create a camera with equal horizontal and vertical focal lengths and a principal point
    /// at the center of an image with the given pixel dimensions.
    ///
    /// # Arguments
    /// * `focal_length` – focal length in pixels (applied to both axes)
    /// * `width` – image width in pixels
    /// * `height` – image height in pixels
    pub fn from_focal_length(focal_length: f64, width: u32, height: u32) -> Self {
        Self {
            fx: focal_length,
            fy: focal_length,
            cx: width as f64 / 2.0,
            cy: height as f64 / 2.0,
        }
    }

    /// Project a single world-space point onto the image plane.
    ///
    /// Returns `None` if the point is behind the camera (Z ≤ 0 in camera space).
    ///
    /// # Arguments
    /// * `world_point`: the point to project, expressed in world coordinates
    /// * `iso` camera-to-world transform; its inverse maps world → camera space
    pub fn project(&self, world_point: &impl PCoords<3>, iso: &Iso3) -> Option<Point2> {
        let cam = iso.inverse() * Point3::from(world_point.coords());
        if cam.z <= 0.0 {
            return None;
        }
        let u = self.fx * cam.x / cam.z + self.cx;
        let v = self.fy * cam.y / cam.z + self.cy;
        Some(Point2::new(u, v))
    }

    /// Project a slice of world-space points onto the image plane.
    ///
    /// Points that are behind the camera produce `None` in the returned vector.
    ///
    /// # Arguments
    /// * `world_points`: points expressed in world coordinates
    /// * `iso`: camera-to-world transform
    pub fn project_many(
        &self,
        world_points: &[impl PCoords<3>],
        iso: &Iso3,
    ) -> Vec<Option<Point2>> {
        world_points.iter().map(|p| self.project(p, iso)).collect()
    }

    /// Back-project an image-plane point into a ray in world space.
    ///
    /// The ray originates at the camera center (in world space) and points toward the scene along
    /// the direction corresponding to the given pixel.
    ///
    /// # Arguments
    /// * `image_point`: pixel coordinates `(u, v)` in the image plane
    /// * `iso`: camera-to-world transform; applied to map the camera-space ray into world space
    pub fn back_project(&self, image_point: &Point2, iso: &Iso3) -> Ray {
        let dx = (image_point.x - self.cx) / self.fx;
        let dy = (image_point.y - self.cy) / self.fy;
        // Direction in camera space (unit vector along the ray toward the scene)
        let dir_cam = UnitVec3::new_normalize(Vector3::new(dx, dy, 1.0));
        // Transform origin and direction into world space
        let origin_world = iso * Point3::origin();
        let dir_world = iso * dir_cam;
        Ray::new(origin_world, dir_world.into_inner())
    }

    /// Back-project a slice of image-plane points into rays in world space.
    ///
    /// # Arguments
    /// * `image_points`: pixel coordinates in the image plane
    /// * `iso`: camera-to-world transform
    pub fn back_project_many(&self, image_points: &[Point2], iso: &Iso3) -> Vec<Ray> {
        image_points
            .iter()
            .map(|p| self.back_project(p, iso))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::na::{Translation3, UnitQuaternion};
    use approx::assert_relative_eq;

    fn identity_iso() -> Iso3 {
        Iso3::identity()
    }

    /// Camera placed 10 units behind the world origin along +Z, looking in the +Z direction.
    fn camera_at_z_neg10() -> Iso3 {
        Iso3::from_parts(
            Translation3::new(0.0, 0.0, -10.0),
            UnitQuaternion::identity(),
        )
    }

    #[test]
    fn project_point_on_optical_axis() {
        // A point on the optical axis should land exactly at the principal point.
        let cam = PinholeCamera::new(500.0, 500.0, 320.0, 240.0);
        let iso = identity_iso();
        let world_pt = Point3::new(0.0, 0.0, 5.0);
        let img = cam.project(&world_pt, &iso).unwrap();
        assert_relative_eq!(img.x, 320.0, epsilon = 1e-10);
        assert_relative_eq!(img.y, 240.0, epsilon = 1e-10);
    }

    #[test]
    fn project_known_point() {
        // With fx=fy=1, cx=cy=0: u = X/Z, v = Y/Z
        let cam = PinholeCamera::new(1.0, 1.0, 0.0, 0.0);
        let iso = identity_iso();
        let world_pt = Point3::new(2.0, 3.0, 4.0);
        let img = cam.project(&world_pt, &iso).unwrap();
        assert_relative_eq!(img.x, 2.0 / 4.0, epsilon = 1e-10);
        assert_relative_eq!(img.y, 3.0 / 4.0, epsilon = 1e-10);
    }

    #[test]
    fn project_behind_camera_returns_none() {
        let cam = PinholeCamera::new(500.0, 500.0, 320.0, 240.0);
        let iso = identity_iso();
        // Z is negative → behind the camera
        let world_pt = Point3::new(0.0, 0.0, -1.0);
        assert!(cam.project(&world_pt, &iso).is_none());
    }

    #[test]
    fn back_project_principal_point_along_optical_axis() {
        // Back-projecting the principal point should yield a ray along +Z in world space.
        let cam = PinholeCamera::new(500.0, 500.0, 320.0, 240.0);
        let iso = identity_iso();
        let principal = Point2::new(320.0, 240.0);
        let ray = cam.back_project(&principal, &iso);
        assert_relative_eq!(ray.origin.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(ray.origin.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(ray.origin.z, 0.0, epsilon = 1e-10);
        // Direction should be (0, 0, 1) normalised
        let dir = ray.dir;
        assert_relative_eq!(dir.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(dir.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(dir.z, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn project_then_back_project_roundtrip() {
        // The back-projected ray should pass through the original world point.
        let cam = PinholeCamera::new(800.0, 800.0, 400.0, 300.0);
        let iso = camera_at_z_neg10();
        let world_pts = vec![
            Point3::new(1.0, 2.0, 5.0),
            Point3::new(-3.0, 0.5, 8.0),
            Point3::new(0.0, 0.0, 1.0),
        ];

        for world_pt in &world_pts {
            let img = cam.project(world_pt, &iso).unwrap();
            let ray = cam.back_project(&img, &iso);

            // The world point should lie on the back-projected ray.
            // Compute t: world_pt = ray.origin + t * ray.dir
            let t = (world_pt - ray.origin).dot(&ray.dir);
            let on_ray = ray.origin + ray.dir * t;

            assert_relative_eq!(on_ray.x, world_pt.x, epsilon = 1e-8);
            assert_relative_eq!(on_ray.y, world_pt.y, epsilon = 1e-8);
            assert_relative_eq!(on_ray.z, world_pt.z, epsilon = 1e-8);
        }
    }

    #[test]
    fn camera_origin_is_at_iso_translation() {
        // The ray origin for any image point should be the camera centre in world space.
        let cam = PinholeCamera::new(500.0, 500.0, 320.0, 240.0);
        let iso = camera_at_z_neg10();
        let img_pt = Point2::new(100.0, 150.0);
        let ray = cam.back_project(&img_pt, &iso);
        assert_relative_eq!(ray.origin.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(ray.origin.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(ray.origin.z, -10.0, epsilon = 1e-10);
    }

    #[test]
    fn from_focal_length_sets_principal_point_at_centre() {
        let cam = PinholeCamera::from_focal_length(600.0, 1280, 960);
        assert_relative_eq!(cam.cx, 640.0, epsilon = 1e-10);
        assert_relative_eq!(cam.cy, 480.0, epsilon = 1e-10);
        assert_relative_eq!(cam.fx, 600.0, epsilon = 1e-10);
        assert_relative_eq!(cam.fy, 600.0, epsilon = 1e-10);
    }
}
