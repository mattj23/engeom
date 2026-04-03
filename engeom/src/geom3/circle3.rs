use crate::common::PCoords;
use crate::geom3::iso3::IsoExtensions3;
use crate::geom3::manifold::Manifold1Pos3;
use crate::{Iso3, Plane3, Point3, Result, UnitVec3, Vector3};
use parry3d_f64::na::UnitQuaternion;
use std::ops;

/// A flat circle in 3D space, defined by a radius and a world isometry. The circle can be thought
/// of as a circle of radius `r` sitting at the origin of the XY plane, then transformed by the
/// isometry to position and orientation in 3D space.
#[derive(Debug, Clone, PartialEq)]
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
        let mut new_circle = self.clone();
        new_circle.transform_by(iso);
        new_circle
    }

    /// Transforms this circle in place by the given isometry.
    pub fn transform_by(&mut self, iso: &Iso3) {
        self.iso = iso * self.iso;
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

    /// Intersects the circle with a plane, returning 0, 1, or 2 intersection angles (in radians).
    ///
    /// Returns an empty vec if the plane does not intersect the circle, a single angle if the
    /// plane is tangent to the circle, or two angles if the plane cuts through it. Returns empty
    /// for the degenerate case where the plane is parallel to (or coincident with) the circle's
    /// plane. Angles are in the range `(-π, π]` and can be passed directly to `at_angle`.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to intersect with
    ///
    /// returns: Vec<f64>
    pub fn intersect_plane(&self, plane: &Plane3) -> Vec<f64> {
        let x_axis = self.iso * Vector3::x();
        let y_axis = self.iso * Vector3::y();
        let a = self.radius * plane.normal.dot(&x_axis);
        let b = self.radius * plane.normal.dot(&y_axis);
        let c = plane.d - plane.normal.dot(&self.center().coords);

        let r_amp = (a * a + b * b).sqrt();
        if r_amp < 1e-10 {
            // Plane normal is (anti)parallel to circle normal — plane is parallel to the circle.
            return vec![];
        }

        let ratio = c / r_amp;
        if ratio > 1.0 + 1e-10 || ratio < -1.0 - 1e-10 {
            return vec![];
        }

        let phi = b.atan2(a);
        let delta = ratio.clamp(-1.0, 1.0).acos();

        if delta < 1e-9 || (std::f64::consts::PI - delta) < 1e-9 {
            vec![phi + delta]
        } else {
            vec![phi + delta, phi - delta]
        }
    }

    /// Rotates the circle's isometry around its normal so that the point currently at `angle`
    /// becomes the new zero angle (aligned with the local x-axis).
    pub fn set_zero_angle(&mut self, angle: f64) {
        let rot = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), angle);
        self.iso.rotation *= rot;
    }

    /// Returns the angle (in radians) of the point on the circle that maximizes the dot product
    /// with the given direction vector. Returns an error if the direction is parallel to the
    /// circle's normal (all points on the circle are equidistant in that direction).
    ///
    /// # Arguments
    ///
    /// * `direction`: a vector in world space; the returned angle maximizes `point · direction`
    ///
    /// returns: Result<f64>
    pub fn max_extent_angle(&self, direction: &Vector3) -> Result<f64> {
        let local = self.iso.inverse().rotation * direction;
        let mag = (local.x * local.x + local.y * local.y).sqrt();
        if mag < 1e-10 {
            return Err("direction is parallel to the circle's normal".into());
        }
        Ok(local.y.atan2(local.x))
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

impl ops::Mul<Circle3> for Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: Circle3) -> Circle3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<&Circle3> for Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: &Circle3) -> Circle3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<Circle3> for &Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: Circle3) -> Circle3 {
        rhs.new_transformed_by(self)
    }
}

impl ops::Mul<&Circle3> for &Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: &Circle3) -> Circle3 {
        rhs.new_transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Curve3;
    use crate::common::linear_space;
    use crate::geom3::tests::RandomGeometry;
    use approx::assert_relative_eq;
    use rand::RngExt;
    use std::f64::consts::PI;

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
        let mut rg = RandomGeometry::new();

        for _ in 0..n {
            let circle = random_circle();

            let points = linear_space(-PI, PI, 1000)
                .iter()
                .map(|t| circle.at_angle(*t).point)
                .collect::<Vec<_>>();
            let curve = Curve3::from_points(&points, 1e-10)?;

            for _ in 0..10 {
                let test_pt = rg.point3(10.0);
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
        let mut rg = RandomGeometry::new();
        let r = rg.sample_f64(0.8, 5.0);
        let iso = rg.iso3(10.0);
        Circle3::new(r, iso)
    }

    #[test]
    fn plane_through_diameter_gives_two_angles() {
        // Circle in XY plane, radius 3 at origin. Plane XZ (y=0) cuts at (±3, 0, 0).
        let circle = Circle3::from_point_normal(
            &Point3::origin(),
            &UnitVec3::new_normalize(Vector3::z()),
            3.0,
        )
        .unwrap();
        let plane = Plane3::xz(); // y=0
        let angles = circle.intersect_plane(&plane);
        assert_eq!(angles.len(), 2);
        // Both angles must produce points on the plane and on the circle
        for &a in &angles {
            let pt = circle.at_angle(a).point;
            assert_relative_eq!(plane.signed_distance_to_point(&pt), 0.0, epsilon = 1e-10);
            assert_relative_eq!((pt - circle.center()).norm(), circle.r(), epsilon = 1e-10);
        }
        // The two points should be antipodal
        let p0 = circle.at_angle(angles[0]).point;
        let p1 = circle.at_angle(angles[1]).point;
        let mid = Point3::from((p0.coords + p1.coords) / 2.0);
        assert_relative_eq!(mid, circle.center(), epsilon = 1e-10);
    }

    #[test]
    fn plane_tangent_gives_one_angle() {
        // Circle in XY plane at origin, radius 2. Plane y=2 is tangent at (0,2,0).
        let circle = Circle3::from_point_normal(
            &Point3::origin(),
            &UnitVec3::new_normalize(Vector3::z()),
            2.0,
        )
        .unwrap();
        let plane = Plane3::new(Vector3::y_axis(), 2.0);
        let angles = circle.intersect_plane(&plane);
        assert_eq!(angles.len(), 1);
        let pt = circle.at_angle(angles[0]).point;
        assert_relative_eq!(pt, Point3::new(0.0, 2.0, 0.0), epsilon = 1e-10);
    }

    #[test]
    fn plane_misses_gives_empty() {
        let circle = Circle3::from_point_normal(
            &Point3::origin(),
            &UnitVec3::new_normalize(Vector3::z()),
            1.0,
        )
        .unwrap();
        let plane = Plane3::new(Vector3::y_axis(), 5.0); // y=5, outside circle
        assert!(circle.intersect_plane(&plane).is_empty());
    }

    #[test]
    fn parallel_plane_gives_empty() {
        // Plane parallel to the circle's own plane
        let circle = Circle3::from_point_normal(
            &Point3::origin(),
            &UnitVec3::new_normalize(Vector3::z()),
            1.0,
        )
        .unwrap();
        let plane = Plane3::new(Vector3::z_axis(), 1.0); // z=1, parallel but offset
        assert!(circle.intersect_plane(&plane).is_empty());
    }

    #[test]
    fn stress_intersect_plane_angles_on_circle_and_plane() {
        // For a random circle and a plane through two known on-circle points, verify the
        // returned angles produce points that lie on both the circle and the plane.
        for _ in 0..500 {
            let circle = random_circle();
            let mut rng = rand::rng();
            let angle = rng.random_range(-PI..PI);
            // Build a plane that passes through a diameter of the circle
            let p1 = circle.at_angle(angle).point;
            let p2 = circle.at_angle(angle + PI).point;
            let some_other = p1 + Vector3::z() * 3.0;
            let plane = Plane3::from((&p1, &p2, &some_other));

            let angles = circle.intersect_plane(&plane);
            assert!(
                !angles.is_empty(),
                "expected intersection with diameter plane"
            );
            for &a in &angles {
                let pt = circle.at_angle(a).point;
                assert_relative_eq!(
                    plane.signed_distance_to_point(&pt).abs(),
                    0.0,
                    epsilon = 1e-8
                );
                assert_relative_eq!((pt - circle.center()).norm(), circle.r(), epsilon = 1e-8);
            }
        }
    }

    #[test]
    fn stress_max_extent_angle() {
        let n = 1000;
        let mut rng = rand::rng();

        for _ in 0..n {
            let circle = random_circle();
            let dir = Vector3::new(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
            );

            // Skip degenerate directions (parallel to normal)
            let Ok(angle) = circle.max_extent_angle(&dir) else {
                continue;
            };

            let best = circle.at_angle(angle).point;
            let best_dot = best.coords.dot(&dir);

            // Sample 5000 points and verify none exceeds the returned dot product
            for &a in linear_space(-PI, PI, 5000).iter() {
                let pt = circle.at_angle(a).point;
                assert!(
                    pt.coords.dot(&dir) <= best_dot + 1e-6,
                    "found point with larger dot product than max_extent_angle result"
                );
            }
        }
    }

    #[test]
    fn stress_set_zero_angle() {
        let n = 1000;
        let mut rng = rand::rng();

        for _ in 0..n {
            let original = random_circle();
            let angle = rng.random_range(-PI..PI);

            // Capture the manifold position at `angle` before re-zeroing
            let expected = original.at_angle(angle);

            let mut circle = original.clone();
            circle.set_zero_angle(angle);

            // After re-zeroing, angle 0 should land on the same world-space point
            let result = circle.at_angle(0.0);

            assert_relative_eq!(result.point, expected.point, epsilon = 1e-12);
            assert_relative_eq!(
                result.direction.into_inner(),
                expected.direction.into_inner(),
                epsilon = 1e-12
            );
        }
    }

    // -------------------------------------------------------------------------
    // Transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn new_transformed_by_identity_preserves_all() {
        let circle = tilted_circle();
        let result = circle.new_transformed_by(&Iso3::identity());
        assert_relative_eq!(result.r(), circle.r(), epsilon = 1e-12);
        assert_relative_eq!(result.center(), circle.center(), epsilon = 1e-12);
        assert_relative_eq!(
            result.normal().into_inner(),
            circle.normal().into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn transform_by_identity_preserves_all() {
        let mut circle = tilted_circle();
        let original_center = circle.center();
        let original_normal = circle.normal().into_inner();
        let original_r = circle.r();
        circle.transform_by(&Iso3::identity());
        assert_relative_eq!(circle.r(), original_r, epsilon = 1e-12);
        assert_relative_eq!(circle.center(), original_center, epsilon = 1e-12);
        assert_relative_eq!(
            circle.normal().into_inner(),
            original_normal,
            epsilon = 1e-12
        );
    }

    #[test]
    fn stress_transform_by() {
        let mut rg = RandomGeometry::new();
        for _ in 0..1000 {
            let original = random_circle();
            let iso = rg.iso3(10.0);

            let moved = original.new_transformed_by(&iso);

            for a in linear_space(-PI, PI, 50).iter() {
                let test_point = moved.at_angle(*a).point;
                let expected = original.at_angle(*a).point;
                assert_relative_eq!(test_point, iso * expected, epsilon = 1e-12);
            }
        }
    }
}
