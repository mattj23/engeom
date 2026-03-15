use crate::AngleDir::{Ccw, Cw};
use crate::common::points::dist;
use crate::common::{
    ANGLE_TOL, PCoords, angle_in_direction, angle_signed_pi, shortest_angle_between,
};
use crate::geom2::aabb2::arc_aabb2;
use crate::geom2::{Aabb2, BoundaryElement, HasBounds2, ManifoldPosition2, directed_angle, rot90};
use crate::{AngleInterval, Circle2, Point2, UnitVec2};
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Arc2 {
    center: Point2,
    radius: f64,
    angle0: f64,
    angle: f64,
    aabb: Aabb2,
}

impl Arc2 {
    pub fn new(circle: Circle2, angle0: f64, angle: f64) -> Self {
        let aabb = arc_aabb2(&circle, angle0, angle);
        Self {
            center: circle.center,
            radius: circle.r(),
            angle0,
            angle,
            aabb,
        }
    }

    /// Create an arc from a center point, a radius, starting at `angle0` and extending for
    /// `angle` radians.
    ///
    /// # Arguments
    ///
    /// * `center`: The arc center point
    /// * `radius`: The arc radius
    /// * `angle0`: The angle in radians (with respect to the x-axis) at which the arc starts
    /// * `angle`: The angle in radians which the arc sweeps through, beginning at `angle0`. A
    ///   positive value indicates a counter-clockwise sweep, while a negative value indicates a
    ///   clockwise sweep.
    ///
    /// returns: Arc2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn circle_angles(center: Point2, radius: f64, angle0: f64, angle: f64) -> Self {
        let circle = Circle2::from_point(center, radius);
        let aabb = arc_aabb2(&circle, angle0, angle);
        Self {
            center,
            radius,
            angle0,
            angle,
            aabb,
        }
    }

    /// Create an arc from a center point, a radius, a point on the perimeter, and an included
    /// angle starting at the point.
    ///
    /// # Arguments
    ///
    /// * `center`: The arc center point
    /// * `radius`: The arc radius
    /// * `point`: A point on the perimeter of the arc at which the arc starting point should be
    ///   located
    /// * `angle`: The angle in radians which the arc sweeps through, beginning at the point. A
    ///   positive value indicates a counter-clockwise sweep, while a negative value indicates a
    ///   clockwise sweep.
    ///
    /// returns: Arc2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn circle_point_angle(center: Point2, radius: f64, point: Point2, angle: f64) -> Self {
        let circle = Circle2::from_point(center, radius);
        let angle0 = circle.angle_of_point(&point);
        let aabb = arc_aabb2(&circle, angle0, angle);
        Self {
            center,
            radius,
            angle0,
            angle,
            aabb,
        }
    }

    /// Create an arc from three points. The arc will begin at the first point, pass through the
    /// second point, and end at the third point.
    ///
    /// # Arguments
    ///
    /// * `p0`: The starting point of the arc
    /// * `p1`: A point on the arc, between the start and end points
    /// * `p2`: The ending point of the arc
    ///
    /// returns: Arc2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn three_points(p0: Point2, p1: Point2, p2: Point2) -> Self {
        let circle = Circle2::from_3_points(&p0, &p1, &p2).unwrap();
        let angle0 = circle.angle_of_point(&p0);
        let v0 = p0 - circle.center;
        let v2 = p2 - circle.center;

        let det = (p1.x - p0.x) * (p1.y + p0.y)
            + (p2.x - p1.x) * (p2.y + p1.y)
            + (p0.x - p2.x) * (p0.y + p2.y);
        let angle = if det < 0.0 {
            directed_angle(&v0, &v2, Ccw)
        } else {
            -directed_angle(&v0, &v2, Cw)
        };

        let aabb = arc_aabb2(&circle, angle0, angle);
        Self {
            center: circle.center,
            radius: circle.r(),
            angle0,
            angle,
            aabb,
        }
    }

    pub fn angle0(&self) -> f64 {
        self.angle0
    }

    pub fn angle(&self) -> f64 {
        self.angle
    }

    pub fn length(&self) -> f64 {
        self.radius * self.angle.abs()
    }

    pub fn circle(&self) -> Circle2 {
        Circle2::from_point(self.center, self.radius)
    }

    pub fn center(&self) -> Point2 {
        self.center
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn point_at_angle(&self, angle: f64) -> Point2 {
        self.circle().point_at_angle(self.angle0 + angle)
    }

    pub fn point_at_fraction(&self, fraction: f64) -> Point2 {
        self.point_at_angle(self.angle * fraction)
    }

    pub fn point_at_length(&self, length: f64) -> Point2 {
        self.point_at_fraction(length / self.length())
    }

    pub fn a(&self) -> Point2 {
        self.point_at_angle(0.0)
    }

    pub fn b(&self) -> Point2 {
        self.point_at_angle(self.angle)
    }

    pub fn is_ccw(&self) -> bool {
        self.angle > 0.0
    }

    pub fn angle_interval(&self) -> AngleInterval {
        AngleInterval::new(self.angle0, self.angle)
    }

    pub fn is_theta_on_arc(&self, theta: f64) -> bool {
        self.angle_interval().contains(theta)
    }

    pub fn theta_to_fraction(&self, theta: f64) -> f64 {
        let theta = angle_signed_pi(theta);
        if shortest_angle_between(theta, self.angle0).abs() < ANGLE_TOL {
            return 0.0;
        }
        if shortest_angle_between(theta, self.angle0 + self.angle).abs() < ANGLE_TOL {
            return 1.0;
        }

        if self.is_ccw() {
            angle_in_direction(self.angle0, theta, Ccw) / self.angle
        } else {
            angle_in_direction(self.angle0, theta, Cw) / -self.angle
        }
    }

    pub fn at_fraction(&self, fraction: f64) -> Point2 {
        self.point_at_fraction(fraction)
    }

    /// Discretize the arc into a series of points such that the maximum distance between the arc
    /// and the line segments connecting the points is less than the specified tolerance. The first
    /// and last points will be the start and end points of the arc.
    ///
    /// # Arguments
    ///
    /// * `tol`: The maximum allowable distance between the arc and the line segments connecting
    ///   the points.
    ///
    /// returns: Vec<OPoint<f64, Const<2>>, Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn make_points(&self, tol: f64) -> Vec<Point2> {
        // The difference between the radius and r*cos(theta/2) is the error between the apogee of
        // each arc segment (of angle theta) and the chord connecting its endpoints. We want to
        // calculate what the angle theta should be that would meet that tolerance, knowing that we
        // must subdivide the total arc angle so that each segment is no larger than that angle.
        let max_theta = 2.0 * (1.0 - tol / self.radius()).acos();
        let n_segments = (self.angle.abs() / max_theta).ceil().max(1.0) as usize;
        let angle_step = self.angle / n_segments as f64;
        let mut points = Vec::with_capacity(n_segments + 1);

        for i in 0..=n_segments {
            let angle = angle_step * i as f64;
            points.push(self.point_at_angle(angle));
        }
        points
    }
}

impl HasBounds2 for Arc2 {
    fn aabb(&self) -> &Aabb2 {
        &self.aabb
    }
}

impl BoundaryElement for Arc2 {
    fn length(&self) -> f64 {
        Arc2::length(self)
    }

    fn at_length(&self, length: f64) -> ManifoldPosition2 {
        // Get the point at the specified length along the arc
        let point = self.point_at_length(length);

        // The normal direction will be either towards or away from the center depending on whether
        // the arc is going clockwise or counterclockwise
        let normal = UnitVec2::new_normalize((point - self.center) * self.angle.signum());

        // The manifold direction will be the normal direction rotated 90 degrees counter-clockwise
        let direction = rot90(Ccw) * normal;

        ManifoldPosition2::new(length, point, direction, normal)
    }

    fn closest_to_point(&self, point: &impl PCoords<2>) -> ManifoldPosition2 {
        let theta = self.circle().angle_of_point(point);
        let t = if self.is_theta_on_arc(theta) {
            self.theta_to_fraction(theta) * self.length()
        } else {
            let d0 = dist(&self.at_start(), point);
            let d1 = dist(&self.at_end(), point);
            if d0 < d1 { 0.0 } else { self.length() }
        };

        self.at_length(t)
    }

    fn aabb(&self) -> Aabb2 {
        self.aabb
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::mid_point;

    use crate::geom2::tests::Random2;
    use crate::{Arc2, Curve2};
    use approx::assert_relative_eq;

    use std::f64::consts::PI;

    #[test]
    fn three_point_arc_ccw() {
        let p0 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p2 = Point2::new(0.0, -1.0);
        let arc = Arc2::three_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center());
        assert_relative_eq!(1.0, arc.radius());
        assert_relative_eq!(0.0, arc.angle0);
        assert_relative_eq!(3.0 * PI / 2.0, arc.angle);
    }

    #[test]
    fn three_point_arc_cw() {
        let p2 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p0 = Point2::new(0.0, -1.0);
        let arc = Arc2::three_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center());
        assert_relative_eq!(1.0, arc.radius());
        assert_relative_eq!(-PI / 2.0, arc.angle0);
        assert_relative_eq!(-3.0 * PI / 2.0, arc.angle);
    }

    #[test]
    fn make_points_tol() {
        // This test generates points along a quarter-circle arc and checks that the maximum
        // distance from the arc to the line segments connecting the points is within the specified
        // tolerance.
        let arc = Arc2::circle_angles(Point2::origin(), 1.0, 0.0, PI / 2.0);
        let tol = 0.001;
        let points = arc.make_points(tol);
        assert!(points.len() >= 3);

        for i in 0..points.len() - 1 {
            let p0 = points[i];
            let p1 = points[i + 1];
            let mid = mid_point(&p0, &p1);

            let d = arc.circle().distance_to(&mid);
            assert!(d.abs() < tol, "Distance {} exceeds tolerance {}", d, tol);
            assert_relative_eq!(0.0, arc.circle().distance_to(&p0), epsilon = 1e-8);
            assert_relative_eq!(0.0, arc.circle().distance_to(&p1), epsilon = 1e-8);
        }
    }

    #[test]
    fn make_points_ends() {
        // This test checks that the first and last points generated by make_points are exactly
        // the start and end points of the arc.
        let arc = Arc2::circle_angles(Point2::origin(), 1.0, 0.0, PI / 2.0);
        let tol = 0.01;
        let points = arc.make_points(tol);

        assert_relative_eq!(arc.a(), points.first().unwrap(), epsilon = 1e-8);
        assert_relative_eq!(arc.b(), points.last().unwrap(), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_center() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::circle_angles(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(0.0, 2.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(0.0, 1.0), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_start() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::circle_angles(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(2.0, -1.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(1.0, 0.0), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_end() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::circle_angles(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(-2.0, -1.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(-1.0, 0.0), epsilon = 1e-8);
    }

    #[test]
    fn stress_test_closest() {
        let mut rnd = Random2::new();
        let tol = 0.00001;

        for _ in 0..100 {
            let circle = Circle2::from_point(rnd.point(10.0), rnd.positive(5.0) + 0.5);
            let arc = Arc2::new(circle, rnd.angle_sym_pi(), rnd.angle_sym_2pi());

            let points = arc.make_points(tol);
            let curve = Curve2::from_points(&points, tol * 0.01, false).unwrap();
            let arc_len = arc.length() / (points.len() - 1) as f64;

            for _ in 0..1000 {
                let test_point = rnd.point(15.0);
                let expected = curve.at_closest_to_point(&test_point);
                let actual = arc.closest_to_point(&test_point);

                // The difference in distance must be within the original tolerance
                let expected_dist = dist(&expected, &test_point);
                let actual_dist = dist(&actual.point, &test_point);
                assert_relative_eq!(expected_dist, actual_dist, epsilon = tol);

                // The distance between the two points must be less than the distance between
                // curve vertices
                let d = dist(&expected, &actual);
                assert!(
                    d < arc_len,
                    "Distance {} exceeds arc segment length {}",
                    d,
                    arc_len
                );

                let direction_angle = expected.direction().angle(&actual.direction);
                let normal_angle = expected.normal().angle(&actual.normal);

                // Because of the discretization, the directions and normals will have more noise
                // in them, so we use a looser tolerance here
                assert!(
                    direction_angle.abs() < 2e-2,
                    "Direction angle difference too large: {}",
                    direction_angle
                );
                assert!(
                    normal_angle.abs() < 2e-2,
                    "Normal angle difference too large: {}",
                    normal_angle
                );
            }
        }
    }
}
