mod fitting;
mod queries;

use parry3d_f64::na::{Point, SVector, Unit};

/// A cubic Bézier curve in D-dimensional space, defined by four control points.
///
/// The curve is parameterized by a scalar `t`, conventionally taken in the range `[0, 1]`. At
/// `t = 0` the curve passes through the first control point, and at `t = 1` it passes through the
/// fourth. The two interior control points influence the curve's shape but are not generally
/// interpolated by it.
///
/// The position at parameter `t` is the standard cubic Bernstein polynomial combination of the
/// four control points:
///
/// `B(t) = (1 - t)^3 P0 + 3 (1 - t)^2 t P1 + 3 (1 - t) t^2 P2 + t^3 P3`
#[derive(Debug, Clone, PartialEq)]
pub struct CubicSpline<const D: usize> {
    pub p0: Point<f64, D>,
    pub p1: Point<f64, D>,
    pub p2: Point<f64, D>,
    pub p3: Point<f64, D>,
}

impl<const D: usize> CubicSpline<D> {
    /// Creates a new cubic Bézier curve from its four control points, in order from the start of
    /// the curve to the end.
    pub fn new(p0: Point<f64, D>, p1: Point<f64, D>, p2: Point<f64, D>, p3: Point<f64, D>) -> Self {
        Self { p0, p1, p2, p3 }
    }

    /// Returns the position on the curve at the given parameter `t`.
    ///
    /// At `t = 0` this returns the first control point, and at `t = 1` it returns the fourth.
    /// Values outside `[0, 1]` are accepted and will extrapolate the underlying polynomial.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// assert_relative_eq!(curve.position(0.0), Point2::new(0.0, 0.0));
    /// assert_relative_eq!(curve.position(1.0), Point2::new(3.0, 0.0));
    /// assert_relative_eq!(curve.position(0.5), Point2::new(1.5, 0.75));
    /// ```
    pub fn position(&self, t: f64) -> Point<f64, D> {
        let u = 1.0 - t;
        let b0 = u * u * u;
        let b1 = 3.0 * u * u * t;
        let b2 = 3.0 * u * t * t;
        let b3 = t * t * t;
        Point::from(
            self.p0.coords * b0 + self.p1.coords * b1 + self.p2.coords * b2 + self.p3.coords * b3,
        )
    }

    /// Returns the derivative of the curve at parameter `t` as a vector (the velocity at `t` if
    /// the curve is interpreted as a parametric trajectory).
    ///
    /// The derivative of a cubic Bezier is a quadratic Bezier whose control vectors are scaled
    /// finite differences of the original control points:
    ///
    /// `B'(t) = 3 (1 - t)^2 (P1 - P0) + 6 (1 - t) t (P2 - P1) + 3 t^2 (P3 - P2)`
    ///
    /// At `t = 0` this evaluates to `3 (P1 - P0)`, the initial tangent direction (scaled). At
    /// `t = 1` it evaluates to `3 (P3 - P2)`, the final tangent direction (scaled).
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, Vector2};
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// assert_relative_eq!(curve.derivative(0.0), Vector2::new(3.0, 3.0));
    /// assert_relative_eq!(curve.derivative(1.0), Vector2::new(3.0, -3.0));
    /// ```
    pub fn derivative(&self, t: f64) -> SVector<f64, D> {
        let u = 1.0 - t;
        (self.p1.coords - self.p0.coords) * (3.0 * u * u)
            + (self.p2.coords - self.p1.coords) * (6.0 * u * t)
            + (self.p3.coords - self.p2.coords) * (3.0 * t * t)
    }

    /// Returns the second derivative of the curve at parameter `t` as a vector (the acceleration
    /// at `t` if the curve is interpreted as a parametric trajectory).
    ///
    /// The second derivative of a cubic Bezier is a linear function of `t`:
    ///
    /// `B''(t) = 6 (1 - t) (P0 - 2 P1 + P2) + 6 t (P1 - 2 P2 + P3)`
    ///
    /// At `t = 0` this evaluates to `6 (P0 - 2 P1 + P2)`. At `t = 1` it evaluates to
    /// `6 (P1 - 2 P2 + P3)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, Vector2};
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// assert_relative_eq!(curve.second_derivative(0.0), Vector2::new(0.0, -6.0));
    /// assert_relative_eq!(curve.second_derivative(1.0), Vector2::new(0.0, -6.0));
    /// ```
    pub fn second_derivative(&self, t: f64) -> SVector<f64, D> {
        let u = 1.0 - t;
        (self.p0.coords - self.p1.coords * 2.0 + self.p2.coords) * (6.0 * u)
            + (self.p1.coords - self.p2.coords * 2.0 + self.p3.coords) * (6.0 * t)
    }

    /// Returns the unit tangent vector of the curve at parameter `t`, i.e. the derivative
    /// normalized to unit length.
    ///
    /// The derivative must be non-zero at `t` for the result to be well-defined. The derivative
    /// can vanish at cusps: for example, when `p0 == p1` (the curve "rests" at the start before
    /// moving) or when `p2 == p3` (it "rests" at the end). Callers that need to handle those cases
    /// gracefully should use [`derivative`](Self::derivative) directly and test the magnitude.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// let t0 = curve.tangent(0.0);
    /// assert_relative_eq!(t0.x, 0.5_f64.sqrt());
    /// assert_relative_eq!(t0.y, 0.5_f64.sqrt());
    /// ```
    pub fn tangent(&self, t: f64) -> Unit<SVector<f64, D>> {
        Unit::new_normalize(self.derivative(t))
    }

    /// Returns the curvature magnitude of the curve at parameter `t`.
    ///
    /// Curvature is the magnitude of the rate of change of the unit tangent with respect to arc
    /// length. For a parametric curve in any dimension this can be written as:
    ///
    /// `κ(t) = √(|B'|² |B''|² − (B'·B'')²) / |B'|³`
    ///
    /// (The numerator is the Lagrange-identity form of `|B' × B''|`, which sidesteps the fact that
    /// the cross product is only natively defined in 3D.) The result is always non-negative; its
    /// reciprocal is the radius of the osculating circle at that point.
    ///
    /// Returns `NaN` if `B'(t) == 0` (a cusp), since curvature is undefined there.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A straight cubic has zero curvature everywhere.
    /// let line = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(2.0, 0.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// assert_relative_eq!(line.curvature(0.3), 0.0, epsilon = 1e-12);
    /// ```
    pub fn curvature(&self, t: f64) -> f64 {
        let d1 = self.derivative(t);
        let d2 = self.second_derivative(t);
        let d1_sq = d1.norm_squared();
        let d2_sq = d2.norm_squared();
        let dot = d1.dot(&d2);
        // Lagrange identity: |a × b|² = |a|²|b|² - (a·b)². Clamp at zero to guard against tiny
        // negative values from floating-point cancellation.
        let cross_sq = (d1_sq * d2_sq - dot * dot).max(0.0);
        cross_sq.sqrt() / (d1_sq * d1_sq.sqrt())
    }

    /// Returns the real roots of the derivative of each component of the curve, found via the
    /// quadratic formula.
    ///
    /// For each dimension `d`, the d-th component of `B'(t)` is a quadratic polynomial in `t`
    /// whose roots are the parameter values where that component has a local extremum. The
    /// method returns those roots as a fixed-size `[[f64; 2]; D]`: row `d` holds up to two
    /// roots for dimension `d`, with unused slots filled with `f64::NAN`. When two roots exist
    /// the smaller one is in slot 0.
    ///
    /// Per dimension, with `a = P1 - P0`, `b = P2 - P1`, `c = P3 - P2`, the component derivative
    /// (up to the constant factor of 3) is
    ///
    /// `(a - 2 b + c) t² + 2 (b - a) t + a`.
    ///
    /// Roots are not filtered to `[0, 1]`. Callers interested only in extrema along the curve as
    /// drawn should discard roots outside that range.
    ///
    /// When a component is everywhere stationary (all four control points share the same value
    /// in that dimension) the derivative is identically zero and every `t` is technically a
    /// root; in this case both slots are `NaN` for that dimension.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// let roots = curve.derivative_roots();
    /// // x is monotonic in t, no roots.
    /// assert!(roots[0][0].is_nan() && roots[0][1].is_nan());
    /// // y peaks at the midpoint.
    /// assert_relative_eq!(roots[1][0], 0.5);
    /// assert!(roots[1][1].is_nan());
    /// ```
    pub fn derivative_roots(&self) -> [[f64; 2]; D] {
        std::array::from_fn(|d| {
            let a = self.p1.coords[d] - self.p0.coords[d];
            let b = self.p2.coords[d] - self.p1.coords[d];
            let c = self.p3.coords[d] - self.p2.coords[d];
            solve_quadratic_real_roots(a - 2.0 * b + c, 2.0 * (b - a), a)
        })
    }

    /// Returns an adaptive polyline approximation of the curve such that the linear interpolation
    /// between any two consecutive points deviates from the underlying spline by no more than the
    /// specified `tolerance` (measured as Euclidean distance).
    ///
    /// The returned `Vec` always starts at `p0` and ends at `p3`. Subdivision is adaptive: regions
    /// where the curve is locally close to straight produce widely spaced points, while regions of
    /// high curvature are subdivided more finely.
    ///
    /// The bound is enforced via the convex-hull property of cubic Bezier curves. Splitting is
    /// performed in control-point space using de Casteljau's algorithm; the flatness of each
    /// sub-curve is measured as the maximum perpendicular distance from its interior control
    /// points (`p1`, `p2`) to its chord line (`p0` to `p3`). Because the curve lies entirely
    /// within the convex hull of its control points, this distance is an upper bound on how far
    /// the sub-curve can deviate from its chord.
    ///
    /// # Arguments
    ///
    /// * `tolerance`: the maximum allowed deviation between the polyline and the underlying spline.
    ///   Must be positive. Smaller values produce more points.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    ///
    /// // A straight line collapses to just its endpoints, regardless of tolerance.
    /// let straight = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(2.0, 0.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let pts = straight.polyline(1e-6);
    /// assert_eq!(pts.len(), 2);
    /// ```
    pub fn polyline(&self, tolerance: f64) -> Vec<Point<f64, D>> {
        const MAX_DEPTH: u32 = 32;
        let mut out = Vec::with_capacity(8);
        out.push(self.p0);
        self.flatten_into(tolerance, MAX_DEPTH, &mut out);
        out
    }

    /// Returns the corners (min, max) of the tight axis-aligned bounding box of the curve over
    /// the parameter range `[0, 1]`.
    ///
    /// Per-dimension extrema of a cubic Bezier occur either at the endpoints `t = 0` and
    /// `t = 1` or at interior points where the corresponding component of `B'(t)` vanishes.
    /// Those interior candidates come from [`derivative_roots`](Self::derivative_roots); roots
    /// outside `[0, 1]` and `NaN` slots are ignored.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    ///
    /// let (lo, hi) = curve.compute_bounds();
    /// assert_relative_eq!(lo, Point2::new(0.0, 0.0));
    /// assert_relative_eq!(hi, Point2::new(3.0, 0.75));
    /// ```
    pub fn compute_bounds(&self) -> (Point<f64, D>, Point<f64, D>) {
        let mut lo = SVector::<f64, D>::zeros();
        let mut hi = SVector::<f64, D>::zeros();
        for d in 0..D {
            let a = self.p0.coords[d];
            let b = self.p3.coords[d];
            lo[d] = a.min(b);
            hi[d] = a.max(b);
        }

        let roots = self.derivative_roots();
        for d in 0..D {
            for &t in &roots[d] {
                if !(0.0..=1.0).contains(&t) {
                    continue;
                }
                let u = 1.0 - t;
                let v = u * u * u * self.p0.coords[d]
                    + 3.0 * u * u * t * self.p1.coords[d]
                    + 3.0 * u * t * t * self.p2.coords[d]
                    + t * t * t * self.p3.coords[d];
                if v < lo[d] {
                    lo[d] = v;
                }
                if v > hi[d] {
                    hi[d] = v;
                }
            }
        }

        (Point::from(lo), Point::from(hi))
    }

    fn flatten_into(&self, tolerance: f64, depth_remaining: u32, out: &mut Vec<Point<f64, D>>) {
        if depth_remaining == 0 || self.chord_perp_distance() <= tolerance {
            out.push(self.p3);
        } else {
            let (left, right) = self.split(0.5);
            left.flatten_into(tolerance, depth_remaining - 1, out);
            right.flatten_into(tolerance, depth_remaining - 1, out);
        }
    }

    /// Splits the curve at parameter `t` using de Casteljau's algorithm, returning the left and
    /// right sub-curves as new `CubicSpline` instances. Concatenating the left and right curves
    /// reproduces the original.
    ///
    /// It is the client code's responsibility to either ensure that the parameter `t` is in the
    /// range `[0, 1]` or to correctly handle the result of not doing so.
    pub fn split(&self, t: f64) -> (Self, Self) {
        let u = 1.0 - t;
        let p01 = self.p0.coords * u + self.p1.coords * t;
        let p12 = self.p1.coords * u + self.p2.coords * t;
        let p23 = self.p2.coords * u + self.p3.coords * t;
        let p012 = p01 * u + p12 * t;
        let p123 = p12 * u + p23 * t;
        let p0123 = p012 * u + p123 * t;

        let left = Self::new(
            self.p0,
            Point::from(p01),
            Point::from(p012),
            Point::from(p0123),
        );
        let right = Self::new(
            Point::from(p0123),
            Point::from(p123),
            Point::from(p23),
            self.p3,
        );
        (left, right)
    }

    /// Splits the curve at parameter `t`, returning `None` if `t` is not in `[0, 1]`.
    ///
    /// See [`split`](Self::split) for details.
    pub fn try_split(&self, t: f64) -> Option<(Self, Self)> {
        if (0.0..=1.0).contains(&t) {
            Some(self.split(t))
        } else {
            None
        }
    }

    /// Returns the maximum perpendicular distance from the interior control points (`p1`, `p2`) to
    /// the infinite chord line through `p0` and `p3`. When `p0` and `p3` coincide, falls back to
    /// the maximum Euclidean distance from `p1` and `p2` to that shared point.
    fn chord_perp_distance(&self) -> f64 {
        let chord = self.p3.coords - self.p0.coords;
        let chord_len_sq = chord.norm_squared();

        if chord_len_sq < f64::EPSILON {
            let d1 = (self.p1.coords - self.p0.coords).norm();
            let d2 = (self.p2.coords - self.p0.coords).norm();
            return d1.max(d2);
        }

        let v1 = self.p1.coords - self.p0.coords;
        let v2 = self.p2.coords - self.p0.coords;
        let perp1 = v1 - chord * (chord.dot(&v1) / chord_len_sq);
        let perp2 = v2 - chord * (chord.dot(&v2) / chord_len_sq);
        perp1.norm().max(perp2.norm())
    }
}

/// Returns the real roots of `a t² + b t + c` via the quadratic formula, handling the linear
/// (`a == 0`) and identically-zero (`a == b == 0`) degenerate cases.
///
/// Unused slots are filled with `f64::NAN`. When two real roots exist the smaller one is in
/// slot 0.
fn solve_quadratic_real_roots(a: f64, b: f64, c: f64) -> [f64; 2] {
    let nan = f64::NAN;
    if a == 0.0 {
        if b == 0.0 {
            return [nan, nan];
        }
        return [-c / b, nan];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        [nan, nan]
    } else if disc == 0.0 {
        [-b / (2.0 * a), nan]
    } else {
        let sqrt_disc = disc.sqrt();
        let inv = 1.0 / (2.0 * a);
        let r1 = (-b - sqrt_disc) * inv;
        let r2 = (-b + sqrt_disc) * inv;
        if r1 <= r2 { [r1, r2] } else { [r2, r1] }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Point2, Point3, Vector2};
    use approx::assert_relative_eq;

    fn sample_2d() -> CubicSpline<2> {
        CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(3.0, 0.0),
        )
    }

    #[test]
    fn endpoints_match_control_points() {
        let c = sample_2d();
        assert_relative_eq!(c.position(0.0), c.p0);
        assert_relative_eq!(c.position(1.0), c.p3);
    }

    #[test]
    fn midpoint_matches_bernstein_formula() {
        let c = sample_2d();
        let t = 0.5;
        let u = 1.0 - t;
        let expected = Point2::from(
            c.p0.coords * (u * u * u)
                + c.p1.coords * (3.0 * u * u * t)
                + c.p2.coords * (3.0 * u * t * t)
                + c.p3.coords * (t * t * t),
        );
        assert_relative_eq!(c.position(t), expected);
    }

    #[test]
    fn straight_line_is_linear_in_t() {
        // Control points collinear and evenly spaced collapse the curve to a straight line.
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert_relative_eq!(
                c.position(t),
                Point3::new(3.0 * t, 0.0, 0.0),
                epsilon = 1e-12
            );
        }
    }

    /// Distance from a point to a line segment, used in the polyline-tolerance tests.
    fn dist_point_to_segment<const D: usize>(
        p: Point<f64, D>,
        a: Point<f64, D>,
        b: Point<f64, D>,
    ) -> f64 {
        let ab = b.coords - a.coords;
        let ab_len_sq = ab.norm_squared();
        if ab_len_sq < f64::EPSILON {
            return (p.coords - a.coords).norm();
        }
        let t = (ab.dot(&(p.coords - a.coords)) / ab_len_sq).clamp(0.0, 1.0);
        (p.coords - (a.coords + ab * t)).norm()
    }

    fn max_curve_to_polyline_distance<const D: usize>(
        curve: &CubicSpline<D>,
        polyline: &[Point<f64, D>],
        samples: usize,
    ) -> f64 {
        let mut max_dist: f64 = 0.0;
        for i in 0..=samples {
            let t = i as f64 / samples as f64;
            let pt = curve.position(t);
            let min_to_polyline = polyline
                .windows(2)
                .map(|w| dist_point_to_segment(pt, w[0], w[1]))
                .fold(f64::INFINITY, f64::min);
            if min_to_polyline > max_dist {
                max_dist = min_to_polyline;
            }
        }
        max_dist
    }

    #[test]
    fn derivative_at_endpoints() {
        let c = sample_2d();
        let d0 = c.derivative(0.0);
        let d1 = c.derivative(1.0);
        // d/dt at t=0 should be 3 * (P1 - P0); at t=1, 3 * (P3 - P2).
        assert_relative_eq!(d0, 3.0 * (c.p1.coords - c.p0.coords));
        assert_relative_eq!(d1, 3.0 * (c.p3.coords - c.p2.coords));
    }

    #[test]
    fn derivative_matches_finite_difference() {
        let c = sample_2d();
        let h = 1e-6;
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let approx = (c.position(t + h).coords - c.position(t - h).coords) / (2.0 * h);
            let exact = c.derivative(t);
            assert_relative_eq!(approx, exact, epsilon = 1e-6);
        }
    }

    #[test]
    fn tangent_is_unit_length() {
        let c = sample_2d();
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let tangent = c.tangent(t);
            assert_relative_eq!(tangent.norm(), 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn tangent_at_endpoints_aligns_with_arm_directions() {
        let c = sample_2d();
        // At t=0, the tangent is parallel to P1 - P0.
        let arm_start = (c.p1.coords - c.p0.coords).normalize();
        assert_relative_eq!(c.tangent(0.0).into_inner(), arm_start, epsilon = 1e-12);
        // At t=1, the tangent is parallel to P3 - P2.
        let arm_end = (c.p3.coords - c.p2.coords).normalize();
        assert_relative_eq!(c.tangent(1.0).into_inner(), arm_end, epsilon = 1e-12);
    }

    #[test]
    fn second_derivative_at_endpoints() {
        let c = sample_2d();
        // B''(0) = 6*(P0 - 2 P1 + P2) = 6*((0,0) - (2,2) + (2,1)) = 6*(0,-1) = (0,-6)
        assert_relative_eq!(
            c.second_derivative(0.0),
            Vector2::new(0.0, -6.0),
            epsilon = 1e-12
        );
        // B''(1) = 6*(P1 - 2 P2 + P3) = 6*((1,1) - (4,2) + (3,0)) = 6*(0,-1) = (0,-6)
        assert_relative_eq!(
            c.second_derivative(1.0),
            Vector2::new(0.0, -6.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn second_derivative_matches_finite_difference() {
        let c = sample_2d();
        let h = 1e-4;
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let approx = (c.derivative(t + h) - c.derivative(t - h)) / (2.0 * h);
            assert_relative_eq!(approx, c.second_derivative(t), epsilon = 1e-6);
        }
    }

    #[test]
    fn curvature_at_endpoints_sample_2d() {
        // For sample_2d, by hand:
        // B'(0) = 3*(P1-P0) = (3, 3),  |B'|² = 18
        // B''(0) = 6*(P0 - 2 P1 + P2) = (0, -6),  |B''|² = 36, B'·B'' = -18
        // |B' × B''|² = 18*36 - 18² = 324 → numerator = 18
        // |B'|³ = 18 * √18 = 54 √2
        // κ(0) = 18 / (54 √2) = 1 / (3 √2)
        let c = sample_2d();
        let expected = 1.0 / (3.0 * 2.0_f64.sqrt());
        assert_relative_eq!(c.curvature(0.0), expected, epsilon = 1e-12);
        assert_relative_eq!(c.curvature(1.0), expected, epsilon = 1e-12); // symmetric
    }

    #[test]
    fn curvature_of_straight_line_is_zero() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert_relative_eq!(c.curvature(t), 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn curvature_matches_numerical_estimate() {
        let c = sample_2d();
        let h = 1e-4;
        for i in 1..10 {
            let t = i as f64 / 10.0;
            let d1 = (c.position(t + h).coords - c.position(t - h).coords) / (2.0 * h);
            let d2 = (c.position(t + h).coords - 2.0 * c.position(t).coords
                + c.position(t - h).coords)
                / (h * h);
            let d1_sq = d1.norm_squared();
            let d2_sq = d2.norm_squared();
            let dot = d1.dot(&d2);
            let numerical = ((d1_sq * d2_sq - dot * dot).max(0.0)).sqrt() / (d1_sq * d1_sq.sqrt());
            assert_relative_eq!(c.curvature(t), numerical, epsilon = 1e-5);
        }
    }

    #[test]
    fn derivative_roots_sample_2d() {
        let c = sample_2d();
        let roots = c.derivative_roots();
        // x component is linear in t (constant derivative), no roots.
        assert!(roots[0][0].is_nan());
        assert!(roots[0][1].is_nan());
        // y component peaks at t = 0.5.
        assert_relative_eq!(roots[1][0], 0.5, epsilon = 1e-12);
        assert!(roots[1][1].is_nan());
    }

    #[test]
    fn derivative_roots_straight_line_has_none() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let roots = c.derivative_roots();
        for dim in 0..3 {
            assert!(roots[dim][0].is_nan(), "dim {} slot 0", dim);
            assert!(roots[dim][1].is_nan(), "dim {} slot 1", dim);
        }
    }

    #[test]
    fn derivative_roots_constant_component_returns_nan() {
        // y is identically zero in all control points: A = B = C = 0 in that component.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
        );
        let roots = c.derivative_roots();
        assert!(roots[1][0].is_nan());
        assert!(roots[1][1].is_nan());
    }

    #[test]
    fn derivative_roots_zero_the_component_derivative() {
        // Pick a curve with a genuine quadratic-with-two-roots component and check that the
        // returned roots really do zero that component of B', and are returned sorted.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 3.0),
            Point2::new(2.0, -3.0),
            Point2::new(3.0, 0.0),
        );
        let roots = c.derivative_roots();
        assert!(roots[1][0].is_finite() && roots[1][1].is_finite());
        assert!(roots[1][0] <= roots[1][1]);
        for &t in &roots[1] {
            assert_relative_eq!(c.derivative(t)[1], 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn compute_bounds_sample_2d() {
        // sample_2d: x is monotonic 0→3, y peaks at t=0.5 with value 0.75.
        let c = sample_2d();
        let (lo, hi) = c.compute_bounds();
        assert_relative_eq!(lo, Point2::new(0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(hi, Point2::new(3.0, 0.75), epsilon = 1e-12);
    }

    #[test]
    fn compute_bounds_straight_line() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let (lo, hi) = c.compute_bounds();
        assert_relative_eq!(lo, Point3::new(0.0, 0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(hi, Point3::new(3.0, 0.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn compute_bounds_loop_with_coincident_endpoints() {
        // p0 == p3 at origin; the curve bulges into the upper half-plane between the inner
        // control points. The bounding box must capture that bulge, not just the endpoints.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(-1.0, 2.0),
            Point2::new(0.0, 0.0),
        );
        let (lo, hi) = c.compute_bounds();

        // Independently approximate the bounds with a fine polyline and compare. The polyline
        // tolerance gives an upper bound on how far each polyline vertex can fall short of the
        // true extremum along the curve, so we compare with that as the slack.
        let tol = 1e-6;
        let pts = c.polyline(tol);
        let mut sample_lo = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut sample_hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in &pts {
            for d in 0..2 {
                if p.coords[d] < sample_lo.coords[d] {
                    sample_lo.coords[d] = p.coords[d];
                }
                if p.coords[d] > sample_hi.coords[d] {
                    sample_hi.coords[d] = p.coords[d];
                }
            }
        }
        assert_relative_eq!(lo, sample_lo, epsilon = tol);
        assert_relative_eq!(hi, sample_hi, epsilon = tol);
    }

    #[test]
    fn polyline_endpoints_are_control_endpoints() {
        let c = sample_2d();
        let pts = c.polyline(1e-3);
        assert!(pts.len() >= 2);
        assert_relative_eq!(*pts.first().unwrap(), c.p0);
        assert_relative_eq!(*pts.last().unwrap(), c.p3);
    }

    #[test]
    fn polyline_collapses_straight_line() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let pts = c.polyline(1e-9);
        assert_eq!(pts.len(), 2);
    }

    #[test]
    fn polyline_respects_tolerance() {
        let c = sample_2d();
        for &tol in &[1e-1, 1e-2, 1e-3, 1e-4] {
            let pts = c.polyline(tol);
            let max_err = max_curve_to_polyline_distance(&c, &pts, 2000);
            assert!(
                max_err <= tol + 1e-9,
                "tolerance {} violated: max_err={} with {} points",
                tol,
                max_err,
                pts.len()
            );
        }
    }

    #[test]
    fn polyline_tighter_tolerance_yields_more_points() {
        let c = sample_2d();
        let coarse = c.polyline(1e-2).len();
        let fine = c.polyline(1e-6).len();
        assert!(fine > coarse, "fine={} should be > coarse={}", fine, coarse);
    }

    #[test]
    fn polyline_works_with_coincident_endpoints() {
        // A loop where p0 == p3 — the degenerate-chord branch of chord_perp_distance must trigger
        // subdivision, and the resulting polyline must still bound the spline within tolerance.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(-1.0, 2.0),
            Point2::new(0.0, 0.0),
        );
        let tol = 1e-3;
        let pts = c.polyline(tol);
        assert!(pts.len() > 2);
        let max_err = max_curve_to_polyline_distance(&c, &pts, 2000);
        assert!(max_err <= tol + 1e-9, "max_err={}", max_err);
    }

    #[test]
    fn split_at_half_preserves_curve() {
        let c = sample_2d();
        let (left, right) = c.split(0.5);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            // Left sub-curve at parameter t corresponds to the original at t/2.
            assert_relative_eq!(left.position(t), c.position(0.5 * t), epsilon = 1e-12);
            // Right sub-curve at parameter t corresponds to the original at 0.5 + t/2.
            assert_relative_eq!(
                right.position(t),
                c.position(0.5 + 0.5 * t),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn split_at_arbitrary_t_preserves_curve() {
        let c = sample_2d();
        for &s in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let (left, right) = c.split(s);
            // Left sub-curve maps [0,1] -> original [0, s].
            for i in 0..=10 {
                let t = i as f64 / 10.0;
                assert_relative_eq!(left.position(t), c.position(s * t), epsilon = 1e-12);
            }
            // Right sub-curve maps [0,1] -> original [s, 1].
            for i in 0..=10 {
                let t = i as f64 / 10.0;
                assert_relative_eq!(
                    right.position(t),
                    c.position(s + (1.0 - s) * t),
                    epsilon = 1e-12
                );
            }
        }
    }

    #[test]
    fn try_split_returns_none_outside_range() {
        let c = sample_2d();
        assert!(c.try_split(-0.1).is_none());
        assert!(c.try_split(1.1).is_none());
        assert!(c.try_split(f64::NAN).is_none());
        assert!(c.try_split(0.5).is_some());
    }

    #[test]
    fn bernstein_basis_sums_to_one() {
        // If all four control points are identical, every position must equal that point.
        let p = Point2::new(7.5, -2.25);
        let c = CubicSpline::new(p, p, p, p);
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            assert_relative_eq!(c.position(t), p, epsilon = 1e-12);
        }
    }
}
