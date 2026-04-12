use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom3::circle3::Circle3;
use crate::geom3::plane3::Plane3;
use crate::geom3::sphere3::Sphere3;
use crate::{Iso3, Point3, Vector3};
use std::ops;

/// A parameterized line in 3D space: `P(t) = origin + t * direction`.
///
/// The direction is not required to be normalized; use `new_normalize` for unit-speed
/// parameterization where `t` equals unit length.
#[derive(Debug, Clone, PartialEq)]
pub struct Line3 {
    origin: Point3,
    direction: Vector3,
}

impl Line3 {
    pub fn x_axis() -> Self {
        Self::new(Point3::origin(), Vector3::x())
    }
    pub fn y_axis() -> Self {
        Self::new(Point3::origin(), Vector3::y())
    }
    pub fn z_axis() -> Self {
        Self::new(Point3::origin(), Vector3::z())
    }

    /// Create a line from an origin point and a direction vector (stored as-is, not normalized).
    pub fn new(origin: Point3, direction: Vector3) -> Self {
        Self { origin, direction }
    }

    /// Moves the origin of the line by a given amount along the direction of the line. A positive
    /// `delta_t` moves the origin forward along the direction of the line, while a negative
    /// `delta_t` moves the origin backward along the direction of the line. The line is modified
    /// in place.
    ///
    /// # Arguments
    ///
    /// * `delta_t`: the amount to move the origin along the direction of the line. If the direction
    ///   is not of unit length, keep in mind this shift will be proportional to the length of the
    ///   direction vector.
    ///
    /// returns: ()
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn shift_origin(&mut self, delta_t: f64) {
        self.origin += self.direction * delta_t;
    }

    /// Creates a new line with the origin shifted by a given amount along the direction of the
    /// line. The direction of the new line is the same as the original line. The original is left
    /// unchanged.
    ///
    /// # Arguments
    ///
    /// * `delta_t`: the amount to move the origin along the direction of the line. If the direction
    ///   is not of unit length, keep in mind this shift will be proportional to the length of the
    ///   direction vector.
    ///
    /// returns: Line3
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new_shifted_origin(&self, delta_t: f64) -> Self {
        Self {
            origin: self.origin + self.direction * delta_t,
            direction: self.direction,
        }
    }

    /// Create a line from an origin point and a direction vector, normalizing the direction so
    /// that the parameter `t` equals arc length from the origin.
    pub fn new_normalize(origin: Point3, direction: Vector3) -> Self {
        Self::new(origin, direction.normalize())
    }

    /// Create a line through two points. The direction is `p2 - p1` (not normalized).
    pub fn from_points(p1: Point3, p2: Point3) -> Self {
        Self::new(p1, p2 - p1)
    }

    pub fn origin(&self) -> Point3 {
        self.origin
    }

    pub fn direction(&self) -> Vector3 {
        self.direction
    }

    /// Normalizes the direction vector in place so that `t` equals arc length from the origin.
    pub fn normalize(&mut self) {
        self.direction = self.direction.normalize();
    }

    /// Returns a new line with the same origin but a normalized direction, so that `t` equals arc
    /// length from the origin.
    pub fn normalized(&self) -> Self {
        Self::new(self.origin, self.direction.normalize())
    }

    /// Returns the point at parameter `t`: `P(t) = origin + t * direction`.
    pub fn at(&self, t: f64) -> Point3 {
        self.origin + self.direction * t
    }

    /// Returns the parameter `t` such that `P(t)` is the closest point on the line to `point`.
    /// Equivalent to the scalar projection of `(point - origin)` onto `direction`, divided by
    /// `|direction|²`.
    pub fn scalar_project(&self, point: &impl PCoords<3>) -> f64 {
        (Point3::from(point.coords()) - self.origin).dot(&self.direction)
            / self.direction.norm_squared()
    }

    /// Returns the closest point on the line to `point`.
    pub fn closest_point(&self, point: &impl PCoords<3>) -> Point3 {
        self.at(self.scalar_project(point))
    }

    /// Returns the perpendicular distance from `point` to the line.
    pub fn distance_to(&self, point: &impl PCoords<3>) -> f64 {
        let pt = Point3::from(point.coords());
        (pt - self.closest_point(&pt)).norm()
    }

    /// Returns a new line with both origin and direction transformed by the given isometry.
    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        let mut clone = self.clone();
        clone.transform_by(iso);
        clone
    }

    /// Transforms this line in place by the given isometry.
    pub fn transform_by(&mut self, iso: &Iso3) {
        self.origin = iso * self.origin;
        self.direction = iso.rotation * self.direction;
    }

    /// Intersects the line with a plane, returning the parameter `t` at the intersection, or
    /// `None` if the line is parallel to (or lies within) the plane.
    pub fn intersect_plane(&self, plane: &Plane3) -> Option<f64> {
        let denom = plane.normal.dot(&self.direction);
        if denom.abs() < 1e-10 {
            return None;
        }
        Some((plane.d - plane.normal.dot(&self.origin.coords)) / denom)
    }

    /// Intersects the line with a sphere, returning 0, 1, or 2 parameters `t` at which the line
    /// meets the sphere surface.
    pub fn intersect_sphere(&self, sphere: &Sphere3) -> Vec<f64> {
        solve_sphere_quadratic(&self.origin, &self.direction, &sphere.center(), sphere.r())
    }

    /// Intersect the line with the plane of a circle and check if the intersection point is within
    /// the circle. Returns `None` if the line does not intersect the circle.
    ///
    /// # Arguments
    ///
    /// * `circle`: a reference to the circle to intersect with.
    ///
    /// returns: Option<f64>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    /// Projects this line onto a plane. The projected origin is the closest point on the plane to
    /// `self.origin`, and the projected direction is `self.direction` with its normal component
    /// removed. Returns `None` if the line is perpendicular to the plane (the projection
    /// degenerates to a point).
    pub fn project_onto_plane(&self, plane: &Plane3) -> Option<Self> {
        let projected_direction =
            self.direction - plane.normal.into_inner() * plane.normal.dot(&self.direction);
        if projected_direction.norm_squared() < 1e-20 {
            return None;
        }
        Some(Self::new(
            plane.project_point(&self.origin),
            projected_direction,
        ))
    }

    /// Intersect the line with the plane of a circle and check if the intersection point is within
    /// the circle. Returns `None` if the line does not intersect the circle.
    ///
    /// # Arguments
    ///
    /// * `circle`: a reference to the circle to intersect with.
    ///
    /// returns: Option<f64>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn intersect_circle(&self, circle: &Circle3) -> Option<f64> {
        if let Some(t) = self.intersect_plane(&circle.plane()) {
            let p = self.at(t);
            if dist(&p, &circle.center()) <= circle.r() {
                Some(t)
            } else {
                None
            }
        } else {
            None
        }
    }
}

/// Solve `|origin + t*dir - center|² = r²` for t, returning the real roots.
fn solve_sphere_quadratic(origin: &Point3, dir: &Vector3, center: &Point3, r: f64) -> Vec<f64> {
    let d = origin - center;
    let a = dir.norm_squared();
    let b = 2.0 * d.dot(dir);
    let c = d.norm_squared() - r * r;
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -1e-10 {
        vec![]
    } else if discriminant.abs() <= 1e-10 {
        vec![-b / (2.0 * a)]
    } else {
        let sq = discriminant.sqrt();
        vec![(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)]
    }
}

impl ops::Mul<Line3> for Iso3 {
    type Output = Line3;
    fn mul(self, rhs: Line3) -> Line3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<&Line3> for Iso3 {
    type Output = Line3;
    fn mul(self, rhs: &Line3) -> Line3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<Line3> for &Iso3 {
    type Output = Line3;
    fn mul(self, rhs: Line3) -> Line3 {
        rhs.new_transformed_by(self)
    }
}

impl ops::Mul<&Line3> for &Iso3 {
    type Output = Line3;
    fn mul(self, rhs: &Line3) -> Line3 {
        rhs.new_transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::sphere3::Sphere3;
    use crate::geom3::tests::RandomGeometry;
    use crate::{Plane3, UnitVec3};
    use approx::assert_relative_eq;

    fn x_axis_line() -> Line3 {
        Line3::new(Point3::origin(), Vector3::x())
    }

    #[test]
    fn at_zero_is_origin() {
        let line = x_axis_line();
        assert_relative_eq!(line.at(0.0), Point3::origin(), epsilon = 1e-12);
    }

    #[test]
    fn at_one_is_origin_plus_direction() {
        let line = Line3::new(Point3::new(1.0, 2.0, 3.0), Vector3::new(0.0, 1.0, 0.0));
        assert_relative_eq!(line.at(1.0), Point3::new(1.0, 3.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(-1.0), Point3::new(1.0, 1.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn new_normalize_gives_unit_direction() {
        let line = Line3::new_normalize(Point3::origin(), Vector3::new(3.0, 0.0, 0.0));
        assert_relative_eq!(line.direction().norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn from_points_direction_is_difference() {
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(4.0, 0.0, 0.0);
        let line = Line3::from_points(p1, p2);
        assert_relative_eq!(
            line.direction(),
            Vector3::new(3.0, 0.0, 0.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn scalar_project_point_on_line() {
        let line = x_axis_line();
        let pt = Point3::new(5.0, 0.0, 0.0);
        assert_relative_eq!(line.scalar_project(&pt), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_project_perpendicular_offset() {
        // Point directly above the origin — projection is at t=0
        let line = x_axis_line();
        let pt = Point3::new(0.0, 3.0, 0.0);
        assert_relative_eq!(line.scalar_project(&pt), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_project_non_unit_direction() {
        // Direction has length 2; at(1) = (2,0,0), at(2) = (4,0,0)
        let line = Line3::new(Point3::origin(), Vector3::new(2.0, 0.0, 0.0));
        let pt = Point3::new(4.0, 7.0, 0.0);
        assert_relative_eq!(line.scalar_project(&pt), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_project_normalized_equals_arc_length() {
        let line = Line3::new_normalize(Point3::origin(), Vector3::new(1.0, 1.0, 0.0));
        let pt = line.at(3.7);
        assert_relative_eq!(line.scalar_project(&pt), 3.7, epsilon = 1e-12);
    }

    #[test]
    fn closest_point_on_line_is_point_itself() {
        let line = x_axis_line();
        let pt = Point3::new(3.0, 0.0, 0.0);
        assert_relative_eq!(line.closest_point(&pt), pt, epsilon = 1e-12);
    }

    #[test]
    fn closest_point_perpendicular_drop() {
        let line = x_axis_line();
        let pt = Point3::new(4.0, 3.0, 0.0);
        assert_relative_eq!(
            line.closest_point(&pt),
            Point3::new(4.0, 0.0, 0.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn distance_to_point_on_line_is_zero() {
        let line = x_axis_line();
        assert_relative_eq!(line.distance_to(&line.at(5.0)), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn distance_to_known_value() {
        // Line along X, point at (0, 3, 4) distance = 5
        let line = x_axis_line();
        assert_relative_eq!(
            line.distance_to(&Point3::new(0.0, 3.0, 4.0)),
            5.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn stress_closest_point_is_perpendicular() {
        let mut rg = RandomGeometry::new();
        for _ in 0..500 {
            let iso = rg.iso3(10.0);
            let line = Line3::new(iso * Point3::origin(), iso.rotation * Vector3::x());
            let pt = rg.point3(10.0);
            let cp = line.closest_point(&pt);
            // Vector from closest point to test point must be perpendicular to direction
            assert_relative_eq!((pt - cp).dot(&line.direction()), 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn intersect_plane_hits() {
        // Line along Z, plane z=5
        let line = Line3::new(Point3::origin(), Vector3::z());
        let plane = Plane3::new(Vector3::z_axis(), 5.0);
        let t = line.intersect_plane(&plane).unwrap();
        assert_relative_eq!(t, 5.0, epsilon = 1e-12);
        assert_relative_eq!(
            plane.signed_distance_to_point(&line.at(t)),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn intersect_plane_parallel_returns_none() {
        let line = Line3::new(Point3::origin(), Vector3::x());
        let plane = Plane3::new(Vector3::z_axis(), 1.0); // z=1, parallel to XY
        assert!(line.intersect_plane(&plane).is_none());
    }

    #[test]
    fn intersect_plane_oblique() {
        let line = Line3::new(Point3::new(0.0, 0.0, -3.0), Vector3::new(1.0, 1.0, 1.0));
        let plane = Plane3::xy(); // z=0
        let t = line.intersect_plane(&plane).unwrap();
        let pt = line.at(t);
        assert_relative_eq!(pt.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn stress_intersect_plane_result_on_plane() {
        let mut rg = RandomGeometry::new();

        for _ in 0..500 {
            let iso = rg.iso3(10.0);
            let origin = iso * Point3::origin();
            let dir = rg.vector3(2.0);
            let line = Line3::new(origin, dir);
            let plane = Plane3::new(UnitVec3::new_normalize(iso.rotation * Vector3::z()), 2.0);

            if let Some(t) = line.intersect_plane(&plane) {
                assert_relative_eq!(
                    plane.signed_distance_to_point(&line.at(t)).abs(),
                    0.0,
                    epsilon = 1e-10
                );
            }
        }
    }

    #[test]
    fn line_through_sphere_two_intersections() {
        let sphere = Sphere3::new(Point3::origin(), 3.0);
        let line = Line3::new(Point3::new(0.0, 0.0, -10.0), Vector3::z());
        let ts = line.intersect_sphere(&sphere);
        assert_eq!(ts.len(), 2);
        for &t in &ts {
            assert_relative_eq!(
                (line.at(t) - sphere.center()).norm(),
                sphere.r(),
                epsilon = 1e-10
            );
        }
        let t_vals: Vec<f64> = {
            let mut v = ts.clone();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v
        };
        // Origin at z=-10 going +z: hits sphere at z=-3 (t=7) and z=+3 (t=13)
        assert_relative_eq!(t_vals[0], 7.0, epsilon = 1e-10);
        assert_relative_eq!(t_vals[1], 13.0, epsilon = 1e-10);
    }

    #[test]
    fn line_misses_sphere() {
        let sphere = Sphere3::new(Point3::origin(), 1.0);
        let line = Line3::new(Point3::new(5.0, 0.0, 0.0), Vector3::z());
        assert!(line.intersect_sphere(&sphere).is_empty());
    }

    #[test]
    fn line_tangent_to_sphere() {
        let sphere = Sphere3::new(Point3::origin(), 1.0);
        // Line at x=1, y=0 along Z is tangent to sphere
        let line = Line3::new(Point3::new(1.0, 0.0, 0.0), Vector3::z());
        let ts = line.intersect_sphere(&sphere);
        assert_eq!(ts.len(), 1);
        let pt = line.at(ts[0]);
        assert_relative_eq!((pt - sphere.center()).norm(), sphere.r(), epsilon = 1e-10);
        assert_relative_eq!(pt.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn stress_intersect_sphere_points_on_surface() {
        let mut rg = RandomGeometry::new();
        for _ in 0..1000 {
            let sphere = Sphere3::new(rg.point3(2.0), rg.f64(0.1, 3.0));
            let line = Line3::new(rg.point3(3.0), rg.vector3(2.0));

            for &t in &line.intersect_sphere(&sphere) {
                let dist = (line.at(t) - sphere.center()).norm();
                assert_relative_eq!(dist, sphere.r(), epsilon = 1e-8);
            }
        }
    }

    #[test]
    fn project_onto_plane_origin_lies_on_plane() {
        let mut rg = RandomGeometry::new();
        for _ in 0..500 {
            let iso = rg.iso3(10.0);
            let plane = Plane3::xy().new_transformed_by(&iso);
            let line = Line3::new(rg.point3(10.0), rg.vector3(1.0));
            if let Some(proj) = line.project_onto_plane(&plane) {
                assert_relative_eq!(
                    plane.signed_distance_to_point(&proj.origin()).abs(),
                    0.0,
                    epsilon = 1e-10
                );
            }
        }
    }

    #[test]
    fn project_onto_plane_direction_parallel_to_plane() {
        let mut rg = RandomGeometry::new();
        for _ in 0..500 {
            let iso = rg.iso3(10.0);
            let plane = Plane3::xy().new_transformed_by(&iso);
            let line = Line3::new(rg.point3(10.0), rg.vector3(1.0));
            if let Some(proj) = line.project_onto_plane(&plane) {
                // projected direction must be perpendicular to the plane normal
                assert_relative_eq!(plane.normal.dot(&proj.direction()), 0.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn project_onto_plane_perpendicular_line_returns_none() {
        // Line along Z projected onto the XY plane: direction collapses to zero
        let plane = Plane3::xy();
        let line = Line3::new(Point3::new(1.0, 2.0, 5.0), Vector3::z());
        assert!(line.project_onto_plane(&plane).is_none());
    }

    #[test]
    fn project_onto_plane_already_in_plane() {
        // A line lying in the XY plane should be unchanged
        let plane = Plane3::xy();
        let line = Line3::new(Point3::new(1.0, 2.0, 0.0), Vector3::x());
        let proj = line.project_onto_plane(&plane).unwrap();
        assert_relative_eq!(proj.origin(), line.origin(), epsilon = 1e-12);
        assert_relative_eq!(proj.direction(), line.direction(), epsilon = 1e-12);
    }

    #[test]
    fn new_transformed_by_isometry_preserves_point_on_line() {
        let mut rg = RandomGeometry::new();

        for _ in 0..500 {
            let original = Line3::new(rg.point3(10.0), rg.vector3(1.0));
            let iso = rg.iso3(10.0);
            let transformed = original.new_transformed_by(&iso);

            for t in [-3.0, -1.0, -0.001, 0.0, 0.001, 1.0, 3.0] {
                let p0 = original.at(t);
                let p1 = transformed.at(t);
                assert_relative_eq!(iso * p0, p1, epsilon = 1e-12);
            }
        }
    }
}
