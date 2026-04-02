use crate::common::PCoords;
use crate::geom3::iso3::IsoExtensions3;
use crate::geom3::manifold::Manifold1Pos3;
use crate::{Iso3, Plane3, Point3, Result, UnitVec3, Vector3};

/// A flat circle in 3D space, defined by a radius and a world isometry. The circle can be thought
/// of as a circle of radius `r` sitting at the origin of the XY plane, then transformed by the
/// isometry to position and orientation in 3D space.
#[derive(Debug, Clone)]
pub struct Circle3 {
    radius: f64,
    iso: Iso3,
}

impl Circle3 {
    /// Returns the radius of the circle.
    pub fn r(&self) -> f64 {
        self.radius
    }

    /// Returns the isometry that defines the circle's position and orientation in world space.
    pub fn iso(&self) -> &Iso3 {
        &self.iso
    }

    /// Returns the world-space center of the circle.
    pub fn center(&self) -> Point3 {
        self.iso * Point3::origin()
    }

    /// Returns the world-space unit normal of the circle's plane.
    pub fn normal(&self) -> UnitVec3 {
        UnitVec3::new_normalize(self.iso * Vector3::z())
    }

    /// Returns the plane that the circle lies in
    pub fn plane(&self) -> Plane3 {
        self.iso * Plane3::xy()
    }

    /// Returns a new circle transformed by the given isometry, without modifying the original.
    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            radius: self.radius,
            iso: self.iso * iso,
        }
    }

    /// Transforms this circle in place by the given isometry.
    pub fn transform_by(&mut self, iso: &Iso3) {
        self.iso = self.iso * iso;
    }

    /// Create a circle from a radius and an isometry. The circle will lie on the XY plane of the
    /// isometry's local frame, centered at the isometry's origin, with the normal along the
    /// isometry's z-axis.
    pub fn new(radius: f64, iso: Iso3) -> Self {
        Self { radius, iso }
    }

    /// Returns the `Manifold1Pos3` at the given angle (in radians) on the circle. The arc length
    /// `l` is measured from angle 0 (the isometry's local x-axis direction).
    pub fn at_angle(&self, angle: f64) -> Manifold1Pos3 {
        let (sin_a, cos_a) = angle.sin_cos();
        let local_point = Point3::new(self.radius * cos_a, self.radius * sin_a, 0.0);
        let local_tangent = Vector3::new(-sin_a, cos_a, 0.0);
        Manifold1Pos3::new(
            self.radius * angle,
            self.iso * local_point,
            UnitVec3::new_normalize(self.iso.rotation * local_tangent),
        )
    }

    /// Returns the angle (in radians) of the point on the circle closest to `test_point`. The
    /// angle is in the range `(-π, π]` and is measured from the isometry's local x-axis.
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: f64
    pub fn closest_angle(&self, test_point: &impl PCoords<3>) -> f64 {
        let local = self.iso.inverse() * Point3::from(test_point.coords());
        local.y.atan2(local.x)
    }

    /// Returns the `Manifold1Pos3` at the point on the circle closest to `test_point`. The arc
    /// length `l` is measured from angle 0 (the isometry's local x-axis direction).
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: Manifold1Pos3
    pub fn closest_position(&self, test_point: &impl PCoords<3>) -> Manifold1Pos3 {
        self.at_angle(self.closest_angle(test_point))
    }

    /// Create a circle from a center point, a normal direction, and a radius. An arbitrary
    /// orientation around the normal is chosen for the isometry's x/y axes.
    ///
    /// # Arguments
    ///
    /// * `center`: the center point of the circle in world space
    /// * `normal`: the normal direction of the circle's plane
    /// * `radius`: the radius of the circle
    ///
    /// returns: Result<Circle3, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Point3, UnitVec3, Vector3};
    /// use engeom::geom3::Circle3;
    ///
    /// let center = Point3::new(1.0, 2.0, 3.0);
    /// let normal = UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0));
    /// let circle = Circle3::from_point_normal(&center, &normal, 5.0).unwrap();
    ///
    /// assert_relative_eq!(circle.center(), center, epsilon = 1e-12);
    /// assert_relative_eq!(circle.normal().into_inner(), normal.into_inner(), epsilon = 1e-12);
    /// assert_relative_eq!(circle.r(), 5.0, epsilon = 1e-12);
    /// ```
    pub fn from_point_normal(center: &Point3, normal: &UnitVec3, radius: f64) -> Result<Self> {
        // Choose a reference vector not parallel to the normal, so we can form a basis.
        let reference = if normal.z.abs() < 0.9 {
            Vector3::z()
        } else {
            Vector3::x()
        };
        let iso = Iso3::try_from_basis_zy(&normal.into_inner(), &reference, Some(*center))?;
        Ok(Self::new(radius, iso))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::linear_space;
    use crate::geom3::tests::{random_iso3, random_point3};
    use crate::{Curve3, SurfacePoint3};
    use approx::assert_relative_eq;
    use rand::RngExt;
    use std::f64::consts::{FRAC_PI_2, PI};

    /// Build a circle with a non-trivial orientation
    fn tilted_circle() -> Circle3 {
        let center = Point3::new(1.0, 2.0, 3.0);
        let normal = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0));
        Circle3::from_point_normal(&center, &normal, 4.0).unwrap()
    }

    #[test]
    fn at_angle_zero_is_on_x_axis() {
        let circle = tilted_circle();
        let pos = circle.at_angle(0.0);

        // arc length at angle 0 is 0
        assert_relative_eq!(pos.l, 0.0, epsilon = 1e-12);

        // point must lie on the circle
        let d = (pos.point - circle.center()).norm();
        assert_relative_eq!(d, circle.r(), epsilon = 1e-12);
    }

    #[test]
    fn at_angle_arc_length_and_on_circle() {
        let circle = tilted_circle();
        for &angle in &[0.0, 0.5, 1.0, -1.5, PI - 0.001] {
            let pos = circle.at_angle(angle);
            assert_relative_eq!(pos.l, circle.r() * angle, epsilon = 1e-12);

            let d = (pos.point - circle.center()).norm();
            assert_relative_eq!(d, circle.r(), epsilon = 1e-12);
        }
    }

    #[test]
    fn stress_closest_position() -> Result<()> {
        let n = 1000;

        for _ in 0..n {
            let circle = random_circle();

            let points = linear_space(-PI, PI, 1000)
                .iter()
                .map(|t| circle.at_angle(*t).point)
                .collect::<Vec<_>>();
            let curve = Curve3::from_points(&points, 1e-10)?;

            for _ in 0..10 {
                let test_pt = random_point3();
                let expected = curve.at_closest_to_point(&test_pt).point();
                let test_result = circle.closest_position(&test_pt);

                assert_relative_eq!(test_result.point, expected, epsilon = 5e-2);
            }
        }

        Ok(())
    }

    #[test]
    fn stress_closest_angle_roundtrip() {
        let n = 1000;
        for _ in 0..n {
            let circle = random_circle();
            for &angle in &[0.0, 0.5, 1.0, -1.5, PI - 0.001] {
                let pt = circle.at_angle(angle).point;
                assert_relative_eq!(circle.closest_angle(&pt), angle, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn stress_check_tangent() {
        let n = 1000;
        let mut rng = rand::rng();

        for _ in 0..n {
            let circle = random_circle();
            let angle = rng.random_range(-PI..PI);
            let m0 = circle.at_angle(angle);
            let m = circle.at_angle(angle + 1.0f64.to_radians());
            assert!(m0.as_surface_point().scalar_projection(&m) > 0.0);
        }
    }

    fn random_circle() -> Circle3 {
        let mut rng = rand::rng();
        let r = rng.random_range(0.8..5.0);
        let iso = random_iso3();
        Circle3::new(r, iso)
    }
}
