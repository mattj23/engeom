use crate::AngleDir::{Ccw, Cw};
use crate::common::points::dist;
use crate::common::{
    ANGLE_TOL, PCoords, angle_in_direction, angle_signed_pi, shortest_angle_between,
};
use crate::geom2::aabb2::arc_aabb2;
use crate::geom2::circle2::intersection_line_circle;
use crate::geom2::{Aabb2, BoundaryElement, HasBounds2, ManifoldPosition2, directed_angle, rot90};
use crate::{AngleInterval, Circle2, Iso2, Point2, SurfacePoint2, UnitVec2, Vector2};
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Arc2 {
    pub circle: Circle2,
    pub angle0: f64,
    pub angle: f64,
    aabb: Aabb2,
}

impl Arc2 {
    pub fn new(circle: Circle2, angle0: f64, angle: f64) -> Self {
        let aabb = arc_aabb2(&circle, angle0, angle);
        Self {
            circle,
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
            circle,
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
            circle,
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
            circle,
            angle0,
            angle,
            aabb,
        }
    }

    pub fn length(&self) -> f64 {
        self.circle.ball.radius * self.angle.abs()
    }

    pub fn center(&self) -> Point2 {
        self.circle.center
    }

    pub fn radius(&self) -> f64 {
        self.circle.ball.radius
    }

    pub fn point_at_angle(&self, angle: f64) -> Point2 {
        self.circle.point_at_angle(self.angle0 + angle)
    }

    pub fn point_at_fraction(&self, fraction: f64) -> Point2 {
        self.point_at_angle(self.angle * fraction)
    }

    pub fn point_at_length(&self, length: f64) -> Point2 {
        self.point_at_fraction(length / self.length())
    }

    pub fn start(&self) -> Point2 {
        self.point_at_angle(0.0)
    }

    pub fn end(&self) -> Point2 {
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
        let normal = UnitVec2::new_normalize((point - self.circle.center) * self.angle.signum());

        // The manifold direction will be the normal direction rotated 90 degrees clockwise
        let direction = rot90(Cw) * normal;

        ManifoldPosition2::new(length, point, direction, normal)
    }

    fn closest_to_point(&self, point: &impl PCoords<2>) -> ManifoldPosition2 {
        let theta = self.circle.angle_of_point(point);
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
    use crate::Arc2;
    use crate::geom2::Ray2;
    use approx::assert_relative_eq;
    use imageproc::point::Point;
    use rand::Rng;
    use std::f64::consts::PI;
    use test_case::test_case;

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
}
