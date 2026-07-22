use crate::common::PCoords;
use crate::{Iso3, Plane3, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
use std::ops;

mod fitting;

/// A flat circle in 3D space, defined by a center point, a unit normal, and a radius. The circle
/// consists of every point in the plane through `center` perpendicular to `normal` at distance
/// `radius` from `center`.
///
/// This is one of `engeom`'s 3D geometric primitives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle3 {
    pub center: Point3,
    pub normal: UnitVec3,
    pub radius: f64,
}

impl Circle3 {
    /// Create a circle from a center point, a unit normal direction, and a radius.
    ///
    /// # Arguments
    ///
    /// * `center`: the center point of the circle in world space
    /// * `normal`: the normal direction of the circle's plane
    /// * `radius`: the radius of the circle
    ///
    /// returns: Circle3
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
    /// let circle = Circle3::new(center, normal, 5.0);
    ///
    /// assert_relative_eq!(circle.center, center, epsilon = 1e-12);
    /// assert_relative_eq!(circle.normal.into_inner(), normal.into_inner(), epsilon = 1e-12);
    /// assert_relative_eq!(circle.r(), 5.0, epsilon = 1e-12);
    /// ```
    pub fn new(center: Point3, normal: UnitVec3, radius: f64) -> Self {
        Self {
            center,
            normal,
            radius,
        }
    }

    /// Returns the radius of the circle.
    pub fn r(&self) -> f64 {
        self.radius
    }

    /// Computes and returns the plane that the circle lies in.
    pub fn plane(&self) -> Plane3 {
        Plane3::new(self.normal, self.normal.dot(&self.center.coords))
    }

    /// Returns a new circle transformed by the given isometry, without modifying the original.
    pub fn transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            center: iso * self.center,
            normal: UnitVec3::new_normalize(iso.rotation * self.normal.into_inner()),
            radius: self.radius,
        }
    }

    /// Returns a new circle with the normal direction reversed, without modifying the original.
    pub fn normal_reversed(&self) -> Self {
        Self {
            normal: -self.normal,
            ..*self
        }
    }

    /// Returns the point on the circle's perimeter closest to `test_point`, paired with the
    /// tangent direction at that point. The tangent follows the right-hand rule around the
    /// circle's normal (`tangent = normal × (point - center)`).
    ///
    /// Returns `None` if `test_point` lies on the axis through the center along the normal, where
    /// every point on the circle is equally close.
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: Option<SurfacePoint3>
    pub fn closest_point(&self, test_point: &impl PCoords<3>) -> Option<SurfacePoint3> {
        let v = Point3::from(test_point.coords()) - self.center;
        let in_plane = v - self.normal.into_inner() * v.dot(&self.normal);
        let dir = UnitVec3::try_new(in_plane, 1e-10)?;
        let point = self.center + dir.into_inner() * self.radius;
        let tangent = UnitVec3::new_normalize(self.normal.cross(&dir));
        Some(SurfacePoint3::new(point, tangent))
    }

    /// Intersects the circle with a plane, returning 0, 1, or 2 intersection points.
    ///
    /// Returns an empty vec if the plane does not intersect the circle, or for the degenerate
    /// case where the plane is parallel to (or coincident with) the circle's own plane.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to intersect with
    ///
    /// returns: Vec<Point3>
    pub fn intersect_plane(&self, plane: &Plane3) -> Vec<Point3> {
        let Some(line) = self.plane().intersect_plane(plane) else {
            return vec![];
        };

        let closest = line.closest_point(&self.center);
        let d = (closest - self.center).norm();
        if d > self.radius + 1e-10 {
            return vec![];
        }
        if (d - self.radius).abs() < 1e-10 {
            return vec![closest];
        }

        let h = (self.radius * self.radius - d * d).sqrt();
        let dir = line.direction.normalize();
        vec![closest + dir * h, closest - dir * h]
    }

    /// Returns the point on the circle that maximizes the dot product with the given direction
    /// vector. Returns an error if the direction is parallel to the circle's normal (all points
    /// on the circle are equidistant in that direction).
    ///
    /// # Arguments
    ///
    /// * `direction`: a vector in world space; the returned point maximizes `point · direction`
    ///
    /// returns: Result<Point3>
    pub fn max_extent_point(&self, direction: &Vector3) -> Result<Point3> {
        let in_plane = direction - self.normal.into_inner() * direction.dot(&self.normal);
        let dir = UnitVec3::try_new(in_plane, 1e-10)
            .ok_or("direction is parallel to the circle's normal")?;
        Ok(self.center + dir.into_inner() * self.radius)
    }

    /// Create the unique circle passing through three points in 3D space. The circle lies in the
    /// plane defined by the three points, with its normal following the right-hand rule from `p0`
    /// to `p1` to `p2`. Returns an error if the points are collinear (or coincident).
    ///
    /// # Arguments
    ///
    /// * `p0`: the first point
    /// * `p1`: the second point
    /// * `p2`: the third point
    ///
    /// returns: Result<Circle3, Box<dyn Error, Global>>
    pub fn from_3_points(
        p0: &impl PCoords<3>,
        p1: &impl PCoords<3>,
        p2: &impl PCoords<3>,
    ) -> Result<Self> {
        let a = Point3::from(p0.coords());
        let b = Point3::from(p1.coords());
        let c = Point3::from(p2.coords());

        // Circumcenter of the triangle, using the standard vector formula relative to `c`.
        let av = a - c;
        let bv = b - c;
        let axb = av.cross(&bv);
        let denom = 2.0 * axb.norm_squared();
        if denom < 1e-20 {
            return Err("Points are collinear".into());
        }

        let to_center = (av.norm_squared() * bv - bv.norm_squared() * av).cross(&axb) / denom;
        let center = c + to_center;
        let radius = to_center.norm();
        let normal = UnitVec3::new_normalize(axb);
        Ok(Circle3::new(center, normal, radius))
    }
}

impl ops::Mul<Circle3> for Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: Circle3) -> Circle3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<&Circle3> for Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: &Circle3) -> Circle3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<Circle3> for &Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: Circle3) -> Circle3 {
        rhs.transformed_by(self)
    }
}

impl ops::Mul<&Circle3> for &Iso3 {
    type Output = Circle3;
    fn mul(self, rhs: &Circle3) -> Circle3 {
        rhs.transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Curve3;
    use crate::common::linear_space;
    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// Build a circle with a non-trivial orientation
    fn tilted_circle() -> Circle3 {
        let center = Point3::new(1.0, 2.0, 3.0);
        let normal = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0));
        Circle3::new(center, normal, 4.0)
    }

    fn random_circle() -> Circle3 {
        let mut rg = RandomGeometry3::new();
        let r = rg.f64(0.8, 5.0);
        let center = rg.point(10.0);
        let normal = rg.unit_vec();
        Circle3::new(center, normal, r)
    }

    /// Test-only parametrization of a circle's perimeter, used to generate sample points for
    /// verifying `Circle3`'s point-based queries against an independent method. `Circle3` itself
    /// has no angle-based API; this exists purely as test infrastructure.
    fn sample_circle_point(circle: &Circle3, t: f64) -> Point3 {
        let n = circle.normal.into_inner();
        let reference = if n.z.abs() < 0.9 {
            Vector3::z()
        } else {
            Vector3::x()
        };
        let x_axis = reference.cross(&n).normalize();
        let y_axis = n.cross(&x_axis);
        circle.center + x_axis * (circle.radius * t.cos()) + y_axis * (circle.radius * t.sin())
    }

    fn sample_circle_points(circle: &Circle3, n: usize) -> Vec<Point3> {
        linear_space(-PI, PI, n)
            .iter()
            .map(|&t| sample_circle_point(circle, t))
            .collect()
    }

    #[test]
    fn stress_closest_point() -> Result<()> {
        let n = 100;
        let mut rg = RandomGeometry3::new();

        for _ in 0..n {
            let circle = random_circle();
            let points = sample_circle_points(&circle, 1000);
            let curve = Curve3::from_points(&points, 1e-10)?;

            for _ in 0..10 {
                let test_pt = rg.point(10.0);
                let expected = curve.at_closest_to_point(&test_pt).point();
                let test_result = circle.closest_point(&test_pt).unwrap();

                assert_relative_eq!(test_result.point, expected, epsilon = 5e-2);
            }
        }

        Ok(())
    }

    #[test]
    fn closest_point_returns_none_on_axis() {
        let circle = tilted_circle();

        // The center itself is on the axis through the center along the normal.
        assert!(circle.closest_point(&circle.center).is_none());

        // So is any other point directly along the normal from the center.
        let on_axis = circle.center + circle.normal.into_inner() * 3.0;
        assert!(circle.closest_point(&on_axis).is_none());
    }

    #[test]
    fn stress_closest_point_tangent_is_consistent() {
        let n = 1000;
        let mut rg = RandomGeometry3::new();

        for _ in 0..n {
            let circle = random_circle();
            let test_pt = rg.point(10.0);
            let Some(sp) = circle.closest_point(&test_pt) else {
                continue;
            };

            let radial = UnitVec3::new_normalize(sp.point - circle.center);
            assert_relative_eq!(sp.normal.dot(&radial), 0.0, epsilon = 1e-10);
            assert_relative_eq!(sp.normal.dot(&circle.normal), 0.0, epsilon = 1e-10);
            assert_relative_eq!(
                sp.normal.into_inner(),
                circle.normal.cross(&radial),
                epsilon = 1e-10
            );
        }
    }

    #[test]
    fn plane_through_diameter_gives_two_points() {
        // Circle in XY plane, radius 3 at origin. Plane XZ (y=0) cuts at (±3, 0, 0).
        let circle = Circle3::new(Point3::origin(), UnitVec3::new_normalize(Vector3::z()), 3.0);
        let plane = Plane3::xz(); // y=0
        let points = circle.intersect_plane(&plane);
        assert_eq!(points.len(), 2);
        for &pt in &points {
            assert_relative_eq!(plane.signed_distance_to_point(&pt), 0.0, epsilon = 1e-10);
            assert_relative_eq!((pt - circle.center).norm(), circle.r(), epsilon = 1e-10);
        }
        // The two points should be antipodal
        let mid = Point3::from((points[0].coords + points[1].coords) / 2.0);
        assert_relative_eq!(mid, circle.center, epsilon = 1e-10);
    }

    #[test]
    fn plane_tangent_gives_one_point() {
        // Circle in XY plane at origin, radius 2. Plane y=2 is tangent at (0,2,0).
        let circle = Circle3::new(Point3::origin(), UnitVec3::new_normalize(Vector3::z()), 2.0);
        let plane = Plane3::new(Vector3::y_axis(), 2.0);
        let points = circle.intersect_plane(&plane);
        assert_eq!(points.len(), 1);
        assert_relative_eq!(points[0], Point3::new(0.0, 2.0, 0.0), epsilon = 1e-10);
    }

    #[test]
    fn plane_misses_gives_empty() {
        let circle = Circle3::new(Point3::origin(), UnitVec3::new_normalize(Vector3::z()), 1.0);
        let plane = Plane3::new(Vector3::y_axis(), 5.0); // y=5, outside circle
        assert!(circle.intersect_plane(&plane).is_empty());
    }

    #[test]
    fn parallel_plane_gives_empty() {
        // Plane parallel to the circle's own plane
        let circle = Circle3::new(Point3::origin(), UnitVec3::new_normalize(Vector3::z()), 1.0);
        let plane = Plane3::new(Vector3::z_axis(), 1.0); // z=1, parallel but offset
        assert!(circle.intersect_plane(&plane).is_empty());
    }

    #[test]
    fn stress_intersect_plane_points_on_circle_and_plane() {
        // For a random circle and a plane through two known on-circle points, verify the
        // returned points lie on both the circle and the plane.
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let circle = random_circle();
            let t = rg.angle_sym_pi();
            // Build a plane that passes through a diameter of the circle
            let p1 = sample_circle_point(&circle, t);
            let p2 = sample_circle_point(&circle, t + PI);
            let some_other = p1 + Vector3::z() * 3.0;
            let plane = Plane3::from_3_points(&p1, &p2, &some_other).unwrap();

            let points = circle.intersect_plane(&plane);
            assert!(
                !points.is_empty(),
                "expected intersection with diameter plane"
            );
            for &pt in &points {
                assert_relative_eq!(
                    plane.signed_distance_to_point(&pt).abs(),
                    0.0,
                    epsilon = 1e-8
                );
                assert_relative_eq!((pt - circle.center).norm(), circle.r(), epsilon = 1e-8);
            }
        }
    }

    #[test]
    fn stress_max_extent_point() {
        let n = 1000;
        let mut rg = RandomGeometry3::new();

        for _ in 0..n {
            let circle = random_circle();
            let dir = rg.vector(1.0);

            // Skip degenerate directions (parallel to normal)
            let Ok(best) = circle.max_extent_point(&dir) else {
                continue;
            };
            let best_dot = best.coords.dot(&dir);

            // Sample 5000 points and verify none exceeds the returned dot product
            for pt in sample_circle_points(&circle, 5000) {
                assert!(
                    pt.coords.dot(&dir) <= best_dot + 1e-6,
                    "found point with larger dot product than max_extent_point result"
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // Transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn transformed_by_identity_preserves_all() {
        let circle = tilted_circle();
        let result = circle.transformed_by(&Iso3::identity());
        assert_relative_eq!(result.r(), circle.r(), epsilon = 1e-12);
        assert_relative_eq!(result.center, circle.center, epsilon = 1e-12);
        assert_relative_eq!(
            result.normal.into_inner(),
            circle.normal.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn stress_transformed_by() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..1000 {
            let original = random_circle();
            let iso = rg.iso3(10.0);

            let moved = original.transformed_by(&iso);
            assert_relative_eq!(moved.r(), original.r(), epsilon = 1e-9);
            assert_relative_eq!(moved.center, iso * original.center, epsilon = 1e-9);

            for pt in sample_circle_points(&original, 50) {
                let moved_point = iso * pt;
                let back_check = (moved_point - moved.center).norm();
                assert_relative_eq!(back_check, moved.r(), epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn normal_reversed_reverses_normal_preserves_center_and_radius() {
        let circle = tilted_circle();
        let reversed = circle.normal_reversed();
        assert_relative_eq!(reversed.center, circle.center, epsilon = 1e-12);
        assert_relative_eq!(reversed.r(), circle.r(), epsilon = 1e-12);
        assert_relative_eq!(
            reversed.normal.into_inner(),
            -circle.normal.into_inner(),
            epsilon = 1e-12
        );
    }
}
