use crate::AngleDir::Cw;
use crate::UnitVec2;
use crate::common::cubic_spline::CubicSpline;
use crate::common::solve_quadratic_real_roots;
use crate::geom2::rot90;
use crate::na::SVector;

pub type CubicSpline2 = CubicSpline<2>;

impl CubicSpline2 {
    /// Return the unit normal vector of the curve at parameter `t`. The normal vector is the
    /// tangent vector rotated clockwise by 90 degrees.
    pub fn normal(&self, t: f64) -> UnitVec2 {
        rot90(Cw) * self.tangent(t)
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
