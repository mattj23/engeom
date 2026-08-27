use crate::AngleDir::{Ccw, Cw};
use crate::common::consensus::Magsac;
use crate::common::points::dist;
use crate::common::{
    ANGLE_TOL, PCoords, angle_in_direction, angle_signed_pi, angle_to_2pi, linear_space,
    shortest_angle_between,
};
use crate::geom2::aabb2::arc_aabb2;
use crate::geom2::{Aabb2, BoundaryElement2, Manifold1Pos2, directed_angle, rot90};
use crate::{AngleInterval, Circle2, IntervalOps, Point2, Result, UnitVec2};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// A circular arc in 2D space, defined by a center point, a radius, a start angle, and an
/// angular span.
///
/// This is one of `engeom`'s 2D geometric primitives.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Arc2 {
    pub center: Point2,
    pub radius: f64,
    pub angle0: f64,
    pub angle: f64,
}

impl Arc2 {
    /// Create an arc from a circle, a start angle, and a sweep angle.
    ///
    /// # Arguments
    ///
    /// * `circle`: the circle of which the arc is a part
    /// * `angle0`: the angle in radians (with respect to the x-axis) at which the arc starts
    /// * `angle`: the angle in radians which the arc sweeps through, beginning at `angle0`. A
    ///   positive value indicates a counter-clockwise sweep, while a negative value indicates a
    ///   clockwise sweep.
    ///
    /// returns: Arc2
    pub fn from_circle(circle: Circle2, angle0: f64, angle: f64) -> Self {
        Self {
            center: circle.center,
            radius: circle.r(),
            angle0,
            angle,
        }
    }

    pub fn from_ends(
        start: &impl PCoords<2>,
        end: &impl PCoords<2>,
        center: &impl PCoords<2>,
        clockwise: bool,
    ) -> Result<Self> {
        let r = dist(start, center);
        let mismatch = r - dist(start, center);
        if mismatch.abs() > 1e-8 {
            return Err(format!("Arc start and end points do not coincide: {}", mismatch).into());
        };

        let c = Circle2::new(center.coords().x, center.coords().y, r);
        let t0 = c.angle_of_point(start);
        let v0 = start.coords() - c.center.coords();
        let v1 = end.coords() - c.center.coords();
        let t = if clockwise {
            -directed_angle(&v0, &v1, Cw)
        } else {
            directed_angle(&v0, &v1, Ccw)
        };

        // Special case for full circle
        let t = if dist(start, end) < 1e-8 {
            if clockwise { -PI * 2.0 } else { PI * 2.0 }
        } else {
            t
        };

        Ok(Arc2::from_circle(c, t0, t))
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
    pub fn new(center: Point2, radius: f64, angle0: f64, angle: f64) -> Self {
        Self {
            center,
            radius,
            angle0,
            angle,
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
    pub fn from_point_angle(center: Point2, radius: f64, point: Point2, angle: f64) -> Self {
        let circle = Circle2::from_point(center, radius);
        let angle0 = circle.angle_of_point(&point);
        Self {
            center,
            radius,
            angle0,
            angle,
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
    pub fn from_3_points(p0: Point2, p1: Point2, p2: Point2) -> Self {
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

        Self {
            center: circle.center,
            radius: circle.r(),
            angle0,
            angle,
        }
    }

    /// Fit a circular arc to a set of points using MAGSAC++ robust consensus estimation.
    ///
    /// A robust circle is estimated with the same MAGSAC++ consensus fit as
    /// [`Circle2::from_consensus`], rejecting gross outliers. The arc's angular bounds are then set
    /// to the smallest sector (about the fitted circle's center) that contains all of the *inlier*
    /// points, so outliers influence neither the circle nor the arc's extent. The returned arc
    /// always has a non-negative (counter-clockwise) sweep.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to fit the arc to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `min_r`: an optional minimum radius; candidate circles smaller than this are rejected
    /// * `max_r`: an optional maximum radius; candidate circles larger than this are rejected
    /// * `options`: an optional [`Magsac`] configuration to override the iteration count, refinement
    ///   steps, confidence, or RNG seed. Its `sigma_max` field is overridden by the `sigma_max`
    ///   argument.
    ///
    /// returns: Result<Arc2, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point2],
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Self> {
        let min_r = min_r.unwrap_or(0.0);
        let max_r = max_r.unwrap_or(f64::INFINITY);

        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit =
            magsac.fit_filtered::<2, Circle2, _>(points, |c| c.r() >= min_r && c.r() <= max_r)?;
        let circle = fit.model;

        if fit.inliers.is_empty() {
            return Err("Consensus fit produced no inliers to bound the arc".into());
        }

        // Angle of every inlier about the fitted circle's center, normalized to [0, 2*PI).
        let mut angles: Vec<f64> = fit
            .inliers
            .iter()
            .map(|&i| angle_to_2pi(circle.angle_of_point(&points[i])))
            .collect();
        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // The smallest covering sector is the complement of the largest angular gap between
        // consecutive inlier angles. The sorted angles are treated as a cyclic sequence, so the
        // wrap from the last angle back to the first is one of the candidate gaps; the arc starts
        // at the angle just after that largest gap and sweeps counter-clockwise across the rest.
        let n = angles.len();
        let mut max_gap = 0.0;
        let mut start_idx = 0;
        for i in 0..n {
            let next = (i + 1) % n;
            let gap = if next == 0 {
                angles[0] + 2.0 * PI - angles[n - 1]
            } else {
                angles[next] - angles[i]
            };
            if gap > max_gap {
                max_gap = gap;
                start_idx = next;
            }
        }

        let angle0 = angles[start_idx];
        let sweep = 2.0 * PI - max_gap;
        Ok(Arc2::from_circle(circle, angle0, sweep))
    }

    pub fn length(&self) -> f64 {
        self.radius * self.angle.abs()
    }

    pub fn circle(&self) -> Circle2 {
        Circle2::from_point(self.center, self.radius)
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
        AngleInterval::new_start_angle(self.angle0, self.angle)
    }

    pub fn is_theta_on_arc(&self, theta: f64) -> bool {
        self.angle_interval().contains_value(theta)
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

    /// Returns the axis-aligned bounding box of the arc, computed on demand.
    pub fn aabb(&self) -> Aabb2 {
        arc_aabb2(&self.circle(), self.angle0, self.angle)
    }
}

impl BoundaryElement2 for Arc2 {
    fn length(&self) -> f64 {
        Arc2::length(self)
    }

    fn at_length(&self, length: f64) -> Manifold1Pos2 {
        // Get the point at the specified length along the arc
        let point = self.point_at_length(length);

        // The normal direction will be either towards or away from the center depending on whether
        // the arc is going clockwise or counterclockwise
        let normal = UnitVec2::new_normalize((point - self.center) * self.angle.signum());

        // The manifold direction will be the normal direction rotated 90 degrees counter-clockwise
        let direction = rot90(Ccw) * normal;

        Manifold1Pos2::new(length, point, direction, normal)
    }

    ///
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: Manifold1Pos2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn closest_to_point(&self, point: &dyn PCoords<2>) -> Manifold1Pos2 {
        let point = Point2::from(point.coords());
        let theta = self.circle().angle_of_point(&point);
        let t = if self.is_theta_on_arc(theta) {
            self.theta_to_fraction(theta) * self.length()
        } else {
            let d0 = dist(&self.at_start(), &point);
            let d1 = dist(&self.at_end(), &point);
            if d0 < d1 { 0.0 } else { self.length() }
        };

        self.at_length(t)
    }

    fn aabb(&self) -> Aabb2 {
        self.aabb()
    }

    fn to_points(&self, tol: f64) -> Vec<Point2> {
        let theta = 2.0 * ((self.radius - tol) / self.radius).acos();
        let n = (self.angle.abs() / theta).ceil() as usize + 1;
        let mut points = Vec::with_capacity(n);
        for x in linear_space(0.0, 1.0, n).iter() {
            points.push(self.point_at_fraction(*x));
        }

        points
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::mid_point;

    use crate::common::consensus::Magsac;
    use crate::common::random_geometry::RandomGeometry2;
    use crate::{Arc2, Circle2, Curve2, Result};
    use approx::assert_relative_eq;

    use std::f64::consts::PI;

    #[test]
    fn three_point_arc_ccw() {
        let p0 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p2 = Point2::new(0.0, -1.0);
        let arc = Arc2::from_3_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center);
        assert_relative_eq!(1.0, arc.radius);
        assert_relative_eq!(0.0, arc.angle0);
        assert_relative_eq!(3.0 * PI / 2.0, arc.angle);
    }

    #[test]
    fn three_point_arc_cw() {
        let p2 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p0 = Point2::new(0.0, -1.0);
        let arc = Arc2::from_3_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center);
        assert_relative_eq!(1.0, arc.radius);
        assert_relative_eq!(-PI / 2.0, arc.angle0);
        assert_relative_eq!(-3.0 * PI / 2.0, arc.angle);
    }

    #[test]
    fn to_points_tol() {
        // This test generates points along a quarter-circle arc and checks that the maximum
        // distance from the arc to the line segments connecting the points is within the specified
        // tolerance.
        let arc = Arc2::new(Point2::origin(), 1.0, 0.0, PI / 2.0);
        let tol = 0.001;
        let points = arc.to_points(tol);
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
    fn to_points_ends() {
        // This test checks that the first and last points generated by to_points are exactly
        // the start and end points of the arc.
        let arc = Arc2::new(Point2::origin(), 1.0, 0.0, PI / 2.0);
        let tol = 0.01;
        let points = arc.to_points(tol);

        assert_relative_eq!(arc.a(), points.first().unwrap(), epsilon = 1e-8);
        assert_relative_eq!(arc.b(), points.last().unwrap(), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_center() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::new(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(0.0, 2.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(0.0, 1.0), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_start() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::new(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(2.0, -1.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(1.0, 0.0), epsilon = 1e-8);
    }

    #[test]
    fn closest_simple_end() {
        // Arc starts at (1,0), ends at (-1,0), center at (0,0)
        let arc = Arc2::new(Point2::origin(), 1.0, 0.0, PI);
        let test_point = Point2::new(-2.0, -1.0);
        let closest = arc.closest_to_point(&test_point);

        assert_relative_eq!(closest.point, Point2::new(-1.0, 0.0), epsilon = 1e-8);
    }

    #[test]
    fn stress_test_closest() {
        let mut rnd = RandomGeometry2::new();
        let tol = 0.00001;

        for _ in 0..100 {
            let circle = Circle2::from_point(rnd.point(10.0), rnd.positive(5.0) + 0.5);
            let arc = Arc2::from_circle(circle, rnd.angle_sym_pi(), rnd.angle_sym_2pi());

            let points = arc.to_points(tol);
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

    #[test]
    fn arc_from_consensus_bounds_inlier_sector() -> Result<()> {
        use std::f64::consts::TAU;

        // Inliers spread over the CCW sector [0, PI] (an upper half) of a known circle, with a
        // small deterministic radial perturbation.
        let (cx, cy, r) = (2.0, -1.0, 1.3);
        let mut points = Vec::new();
        let inlier_count = 90;
        for i in 0..inlier_count {
            let t = PI * i as f64 / (inlier_count - 1) as f64;
            let rr = r + 0.003 * (5.0 * t).sin();
            points.push(Point2::new(cx + rr * t.cos(), cy + rr * t.sin()));
        }

        // A dense cluster of gross outliers well off the circle (near the empty lower sector).
        for i in 0..40 {
            let t = TAU * i as f64 / 40.0;
            points.push(Point2::new(cx + 0.4 * t.cos(), cy - 4.0 + 0.4 * t.sin()));
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(400),
            refinement_steps: 4,
            confidence: 0.99,
            seed: Some(42),
        };
        let arc = Arc2::from_consensus(&points, 0.02, None, None, Some(magsac))?;

        // The fitted circle is recovered.
        assert_relative_eq!(arc.center.x, cx, epsilon = 5.0e-3);
        assert_relative_eq!(arc.center.y, cy, epsilon = 5.0e-3);
        assert_relative_eq!(arc.radius, r, epsilon = 5.0e-3);

        // The arc spans only the inlier sector [0, PI], not the outlier-adjacent lower half.
        assert!(
            arc.angle > 0.0,
            "sweep should be counter-clockwise (positive)"
        );
        assert_relative_eq!(arc.angle0, 0.0, epsilon = 2.0e-2);
        assert_relative_eq!(arc.angle, PI, epsilon = 2.0e-2);

        Ok(())
    }
}
