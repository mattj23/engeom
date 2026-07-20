use crate::AngleDir::Cw;
use crate::common::cubic_spline::CubicSpline;
use crate::common::solve_quadratic_real_roots;
use crate::geom2::rot90;
use crate::na::SVector;
use crate::{Circle2, UnitVec2};

/// A cubic Bézier curve in 2D space, defined by four control points.
///
/// This is one of `engeom`'s 2D geometric primitives.
///
/// This is the two-dimensional specialization of the dimension-generic
/// [`CubicSpline`](CubicSpline); see that type for the shared constructors and queries (`new`,
/// `position`, `derivative`, `tangent`, `curvature`, `split`, `polyline`, and so on). The methods
/// defined directly on `CubicSpline2` here (`normal`, `curvature_circle`, `find_inflections`) are
/// the ones that only make sense in 2D.
pub type CubicSpline2 = CubicSpline<2>;

impl CubicSpline2 {
    /// Return the unit normal vector of the curve at parameter `t`. The normal vector is the
    /// tangent vector rotated clockwise by 90 degrees.
    pub fn normal(&self, t: f64) -> UnitVec2 {
        rot90(Cw) * self.tangent(t)
    }

    /// Returns the circle of curvature (osculating circle) tangent to the curve at parameter `t`:
    /// the unique circle that matches the curve's position, tangent direction, and curvature there.
    /// Its radius is the reciprocal of the [`curvature`](CubicSpline::curvature), and its center is
    /// the center of curvature, lying on the concave side of the curve.
    ///
    /// Returns `None` where no finite circle exists: where the curve is locally straight (zero
    /// curvature, so the osculating circle degenerates to a line of infinite radius) or where the
    /// curvature is undefined (a cusp, where `B'(t) = 0`).
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A symmetric arch: at the apex t = 0.5 the curvature is 2/3, so the osculating circle has
    /// // radius 3/2 and sits below the apex point (1.5, 0.75) on the concave side.
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let circle = curve.curvature_circle(0.5).unwrap();
    /// assert_relative_eq!(circle.r(), 1.5, epsilon = 1e-12);
    /// assert_relative_eq!(circle.center, Point2::new(1.5, -0.75), epsilon = 1e-12);
    /// ```
    pub fn curvature_circle(&self, t: f64) -> Option<Circle2> {
        let d1 = self.derivative(t);
        let d2 = self.second_derivative(t);
        let speed_sq = d1.norm_squared();

        // The curvature vector is the component of acceleration normal to the velocity, divided by
        // the squared speed. It points from the curve toward the center of curvature and has a
        // magnitude equal to the (unsigned) curvature κ, so the center of curvature is
        // `p + k_vec / |k_vec|²` (= `p + N / κ`) and the radius of curvature is `1 / |k_vec|`.
        let t_hat = d1 / speed_sq.sqrt();
        let k_vec = (d2 - t_hat * d2.dot(&t_hat)) / speed_sq;
        let k_sq = k_vec.norm_squared();

        // A vanishing curvature vector (straight curve) gives an infinite radius, and a cusp makes
        // `t_hat` and hence `k_vec` non-finite; neither yields a circle.
        if !k_sq.is_finite() || k_sq == 0.0 {
            return None;
        }

        let center = self.position(t) + k_vec / k_sq;
        Some(Circle2::from_point(center, 1.0 / k_sq.sqrt()))
    }

    /// Returns the parameter values of any inflection points of the curve, found via the
    /// quadratic formula.
    ///
    /// In 2D, an inflection point is where the signed curvature crosses zero. Equivalently,
    /// the scalar cross product `B'(t) × B''(t)` vanishes. For a cubic Bezier the cubic-degree
    /// term of that cross product cancels identically, leaving a quadratic in `t`, so a cubic
    /// Bezier has at most two inflection points.
    ///
    /// With `a = P1 - P0`, `b = P2 - P1`, `c = P3 - P2`, let
    ///
    /// - `A = a - 2 b + c`
    /// - `Q = 2 (b - a)`
    ///
    /// Then the inflection quadratic is `(A × Q) t² + 2 (A × a) t + (Q × a) = 0`, where `u × v`
    /// here denotes the 2D scalar cross `u.x v.y - u.y v.x`.
    ///
    /// Results are returned as a fixed-size `[f64; 2]`, with unused slots filled with
    /// `f64::NAN` and the smaller root (when two exist) in slot 0. Roots are not filtered to
    /// `[0, 1]`. When the entire curve is a straight line, the quadratic is identically zero
    /// and both slots are `NaN`.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // S-shaped cubic with control points symmetric about t = 0.5: one inflection there.
    /// let s = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, -1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let infl = s.find_inflections();
    /// assert_relative_eq!(infl[0], 0.5);
    /// assert!(infl[1].is_nan());
    /// ```
    pub fn find_inflections(&self) -> [f64; 2] {
        let a = self.p1.coords - self.p0.coords;
        let b = self.p2.coords - self.p1.coords;
        let c = self.p3.coords - self.p2.coords;
        let av = a - 2.0 * b + c;
        let qv = 2.0 * (b - a);

        // 2D scalar cross product: u × v = u.x * v.y - u.y * v.x
        let cross2 = |u: SVector<f64, 2>, v: SVector<f64, 2>| u.x * v.y - u.y * v.x;

        let alpha = cross2(av, qv);
        let beta = 2.0 * cross2(av, a);
        let gamma = cross2(qv, a);

        solve_quadratic_real_roots(alpha, beta, gamma)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point2;
    use approx::assert_relative_eq;

    fn sample_2d() -> CubicSpline2 {
        CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(3.0, 0.0),
        )
    }

    #[test]
    fn curvature_circle_radius_matches_curvature() {
        // The osculating-circle radius must be the reciprocal of the curvature, and the circle
        // must pass through the curve point with that point on its boundary.
        let c = sample_2d();
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let circle = c.curvature_circle(t).expect("finite curvature");
            assert_relative_eq!(circle.r(), 1.0 / c.curvature(t), epsilon = 1e-9);
            // The curve point lies on the circle: its distance to the center equals the radius.
            let p = c.position(t);
            assert_relative_eq!((p - circle.center).norm(), circle.r(), epsilon = 1e-9);
        }
    }

    #[test]
    fn curvature_circle_center_on_concave_side() {
        // sample_2d bulges upward, so the center of curvature is below each curve point.
        let c = sample_2d();
        let t = 0.5;
        let circle = c.curvature_circle(t).expect("finite curvature");
        assert_relative_eq!(circle.r(), 1.5, epsilon = 1e-12);
        assert_relative_eq!(circle.center, Point2::new(1.5, -0.75), epsilon = 1e-12);
    }

    #[test]
    fn curvature_circle_straight_line_is_none() {
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
        );
        assert!(c.curvature_circle(0.3).is_none());
    }

    #[test]
    fn curvature_circle_at_cusp_is_none() {
        // p1 == p2 with mirrored arms gives a cusp at t = 0.5 where B'(t) = 0.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 0.0),
        );
        assert!(c.curvature_circle(0.5).is_none());
    }
}
