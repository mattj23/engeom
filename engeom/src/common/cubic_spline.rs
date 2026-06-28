//! A dimension-generic cubic Bézier curve and the operations built on top of it.
//!
//! The core type is [`CubicSpline`], a single cubic Bézier segment in `D`-dimensional space
//! defined by four control points. Beyond evaluation ([`position`](CubicSpline::position),
//! [`derivative`](CubicSpline::derivative), [`tangent`](CubicSpline::tangent),
//! [`curvature`](CubicSpline::curvature)), it provides structural analysis
//! ([`find_cusp`](CubicSpline::find_cusp), [`find_curvature_zeros`](CubicSpline::find_curvature_zeros)),
//! de Casteljau [`split`](CubicSpline::split)ting, an adaptive [`polyline`](CubicSpline::polyline)
//! approximation, and a tight axis-aligned [`compute_bounds`](CubicSpline::compute_bounds).
//!
//! The type's surface is split across a few places:
//!
//! - The dimension-generic curve and its operations live in this module.
//! - Spatial queries (closest-point projection and the acceleration structure backing it) are
//!   exposed as [`CubicSplineQueries`].
//! - Curve fitting lives in the `fitting` submodule.
//! - Operations that only make sense in 2D (signed-curvature inflections, the left/right normal)
//!   are implemented on `CubicSpline2` (`CubicSpline<2>`) in the `geom2` module, not here. In
//!   particular, the 2D `find_inflections` is the signed-curvature analog of the generic
//!   [`find_curvature_zeros`](CubicSpline::find_curvature_zeros).

mod fitting;
mod queries;

use crate::common::{Line, solve_quadratic_real_roots};
pub use fitting::{SplineBuildFn, SplineFitResult, fit_spline_to_points};
use parry3d_f64::na::{Point, SVector, Unit};
pub use queries::CubicSplineQueries;
use serde::{Deserialize, Serialize};

/// A value of some generic type `T` associated with a parameter value `t` along a spline.
///
/// This pairs a sample of any per-parameter quantity (a position, a derivative, a curvature, etc.)
/// with the parameter `t` at which it was taken.
#[derive(Debug, Clone, Copy)]
pub struct SplineValue<T> {
    /// The parameter value along the spline at which `value` was evaluated.
    pub t: f64,

    /// The value of type `T` at parameter `t`.
    pub value: T,
}

impl<T> SplineValue<T> {
    /// Creates a new [`SplineValue`] pairing a parameter `t` with its associated `value`.
    pub fn new(t: f64, value: T) -> Self {
        Self { t, value }
    }
}

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
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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

    /// Return the position and derivative direction of the curve at parameter `t` in the form of
    /// a parameterized line.
    pub fn line_at(&self, t: f64) -> Line<D> {
        Line::new(self.position(t), self.derivative(t))
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

    /// Returns the parameter `t` of a cusp if one exists in `[0, 1]`, otherwise `None`.
    ///
    /// A cusp on a parametric curve is a point where the velocity vector vanishes, i.e.
    /// `B'(t) = 0`. Every component of `B'` must be zero simultaneously, so any cusp `t`
    /// is necessarily a root of *every* component's derivative quadratic.
    ///
    /// Candidates are drawn from [`derivative_roots`](Self::derivative_roots) restricted to
    /// `[0, 1]`, and each is verified by checking that the full derivative magnitude at that
    /// `t` is within a tolerance scaled to the control-point spacing. If multiple cusps exist,
    /// the smallest `t` is returned.
    ///
    /// Returns `None` for the fully degenerate curve in which all four control points
    /// coincide; that case is best detected explicitly by the caller.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A regular cubic with no cusp.
    /// let smooth = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// assert!(smooth.find_cusp().is_none());
    ///
    /// // p1 == p2 with mirrored start/end arms gives a cusp at t = 0.5.
    /// let cusped = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(0.0, 0.0),
    /// );
    /// let t = cusped.find_cusp().expect("expected a cusp");
    /// assert_relative_eq!(t, 0.5);
    /// ```
    pub fn find_cusp(&self) -> Option<f64> {
        let scale_sq = (self.p1.coords - self.p0.coords)
            .norm_squared()
            .max((self.p2.coords - self.p1.coords).norm_squared())
            .max((self.p3.coords - self.p2.coords).norm_squared());
        if scale_sq == 0.0 {
            return None;
        }
        let eps_sq = scale_sq * 1e-20;

        let roots = self.derivative_roots();
        let mut best: Option<f64> = None;
        for dim_roots in &roots {
            for &t in dim_roots {
                if !(0.0..=1.0).contains(&t) {
                    continue;
                }
                if self.derivative(t).norm_squared() <= eps_sq && best.map_or(true, |b| t < b) {
                    best = Some(t);
                }
            }
        }
        best
    }

    /// Returns parameter values in `[0, 1]` where the curve's curvature is zero.
    ///
    /// This is the dimension-generic analog of the 2D-specific
    /// [`find_inflections`](CubicSpline::<2>::find_inflections). Instead of locating sign
    /// changes of the signed curvature (which only exists in 2D), this method finds parameters
    /// where `B'(t)` and `B''(t)` are linearly dependent, i.e. where the
    /// Lagrange-identity squared cross product
    ///
    /// `|B'(t)|² |B''(t)|² − (B'(t) · B''(t))²`
    ///
    /// equals zero. In 2D this gives the same in-range parameter values as
    /// `find_inflections`. In higher dimensions, true zero-curvature points are rare and
    /// require the curve to be locally straight or planar.
    ///
    /// Candidates are enumerated from the per-axis-pair "cross" quadratics. For each pair
    /// `(i, j)` with `i < j`, the polynomial `B'_i B''_j − B'_j B''_i` is a quadratic in `t`,
    /// and its roots are necessary but not sufficient conditions for zero curvature. Each
    /// candidate is then verified by checking the Lagrange identity against a tolerance scaled
    /// to the control-polygon edge lengths.
    ///
    /// Returns up to two distinct in-range values, with unused slots filled with `f64::NAN`
    /// and the smaller value in slot 0. A curve that is everywhere straight (control points
    /// collinear) returns `[NaN, NaN]`, because no isolated zero exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point3;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A 2D S-shape embedded in the z = 0 plane: one inflection at t = 0.5.
    /// let s = CubicSpline::new(
    ///     Point3::new(0.0, 0.0, 0.0),
    ///     Point3::new(1.0, 1.0, 0.0),
    ///     Point3::new(2.0, -1.0, 0.0),
    ///     Point3::new(3.0, 0.0, 0.0),
    /// );
    /// let zeros = s.find_curvature_zeros();
    /// assert_relative_eq!(zeros[0], 0.5);
    /// assert!(zeros[1].is_nan());
    /// ```
    pub fn find_curvature_zeros(&self) -> [f64; 2] {
        let r = self.p1.coords - self.p0.coords;
        let bv = self.p2.coords - self.p1.coords;
        let cv = self.p3.coords - self.p2.coords;
        let s = 2.0 * (bv - r);
        let w = r - 2.0 * bv + cv;

        let mut candidates: Vec<f64> = Vec::new();
        for i in 0..D {
            for j in (i + 1)..D {
                let alpha = s[i] * w[j] - s[j] * w[i];
                let beta = 2.0 * (r[i] * w[j] - r[j] * w[i]);
                let gamma = r[i] * s[j] - r[j] * s[i];
                for &t in &solve_quadratic_real_roots(alpha, beta, gamma) {
                    if (0.0..=1.0).contains(&t) {
                        candidates.push(t);
                    }
                }
            }
        }
        candidates.sort_by(|a, b| a.partial_cmp(b).unwrap());
        candidates.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

        let scale_sq = r.norm_squared().max(s.norm_squared()).max(w.norm_squared());
        if scale_sq == 0.0 {
            return [f64::NAN; 2];
        }
        let eps_sq = scale_sq * scale_sq * 1e-12;

        let mut result = [f64::NAN; 2];
        let mut idx = 0;
        for t in candidates {
            if idx >= 2 {
                break;
            }
            let d1 = self.derivative(t);
            let d2 = self.second_derivative(t);
            let lagrange = d1.norm_squared() * d2.norm_squared() - d1.dot(&d2).powi(2);
            if lagrange.max(0.0) <= eps_sq {
                result[idx] = t;
                idx += 1;
            }
        }
        result
    }

    /// Finds the point of maximum curvature on the curve over the parameter range `[0, 1]` and
    /// returns it as a [`SplineValue<f64>`]: the parameter `t` at which the maximum occurs paired
    /// with the curvature magnitude there as its `value`. Recover the position on the curve with
    /// [`position`](Self::position) if needed.
    ///
    /// Unlike the per-component extrema in [`derivative_roots`](Self::derivative_roots) or the
    /// curvature *zeros* in [`find_curvature_zeros`](Self::find_curvature_zeros), the maximum of
    /// the curvature function `κ(t)` has no clean closed form: `κ(t)` is a rational function whose
    /// stationary points are the roots of a high-degree polynomial. This method instead brackets
    /// the global maximum with a dense uniform scan and then refines it locally with a
    /// golden-section search, which converges to the parameter to near machine precision provided
    /// the scan is fine enough to isolate the peak.
    ///
    /// Parameters where the curvature is undefined (cusps, where `B'(t) = 0` and
    /// [`curvature`](Self::curvature) returns `NaN`) are treated as having no curvature for the
    /// purposes of the search, so the result is the largest *finite* curvature found. As the
    /// curvature grows without bound approaching a cusp, the reported maximum for a cusped curve
    /// will sit just to one side of the cusp rather than exactly on it.
    ///
    /// For a fully degenerate curve whose curvature is `NaN` everywhere (e.g. all four control
    /// points coincident), the returned `t` is `0.0` and the curvature value is `NaN`.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A symmetric arch: curvature peaks at the apex, t = 0.5.
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let max = curve.find_max_curvature();
    /// // The parameter localizes to about sqrt(machine epsilon) at a flat peak; the curvature
    /// // value itself is far more accurate.
    /// assert_relative_eq!(max.t, 0.5, epsilon = 1e-6);
    /// assert_relative_eq!(max.value, 2.0 / 3.0, epsilon = 1e-9);
    /// ```
    pub fn find_max_curvature(&self) -> SplineValue<f64> {
        // Curvature undefined (cusp) is treated as "no curvature" so the scan ignores it and
        // converges to the largest finite value instead of chasing the NaN.
        let kappa = |t: f64| {
            let k = self.curvature(t);
            if k.is_finite() { k } else { f64::NEG_INFINITY }
        };

        // TODO: I think this approach is conservative and probably a candidate for improvement

        // Coarse uniform scan to bracket the global maximum. The resolution must be fine enough
        // to land inside the basin of the tallest peak; a cubic has at most a handful of
        // curvature maxima, so a few hundred samples is comfortably sufficient.
        const SAMPLES: usize = 256;
        let step = 1.0 / SAMPLES as f64;
        let mut best_t = 0.0;
        let mut best_k = kappa(0.0);
        for i in 1..=SAMPLES {
            let t = i as f64 / SAMPLES as f64;
            let k = kappa(t);
            if k > best_k {
                best_k = k;
                best_t = t;
            }
        }

        // No finite curvature anywhere (e.g. fully coincident control points): nothing to refine.
        if !best_k.is_finite() {
            return SplineValue::new(0.0, self.curvature(0.0));
        }

        // Golden-section refinement on the bracket spanning the samples on either side of the
        // best one. This assumes a single maximum within the bracket, which the dense scan above
        // is responsible for guaranteeing.
        let a = (best_t - step).max(0.0);
        let b = (best_t + step).min(1.0);
        let t = golden_section_max(kappa, a, b);

        // The refined interior point can, in degenerate cases, score below the bracketing scan
        // sample; keep whichever parameter actually yields the larger curvature.
        let t = if kappa(t) >= best_k { t } else { best_t };
        SplineValue::new(t, self.curvature(t))
    }

    /// Returns every local maximum of the curvature over the parameter range `[0, 1]`, each as a
    /// [`SplineValue<f64>`] pairing the parameter `t` where the maximum occurs with the curvature
    /// magnitude there. The results are ordered by ascending `t`.
    ///
    /// Where [`find_max_curvature`](Self::find_max_curvature) returns only the single global peak,
    /// this returns all of them, including local maxima at the domain endpoints `t = 0` and
    /// `t = 1` when the curvature there exceeds its inward neighbor.
    ///
    /// The same coarse-scan-then-golden-section strategy as `find_max_curvature` is used: the
    /// curvature is sampled on a dense uniform grid, each grid point that is at least as large as
    /// both of its neighbors is taken as a bracket for a single peak, and that peak is refined with
    /// a golden-section search. Parameters where the curvature is undefined (cusps, where
    /// `B'(t) = 0`) are treated as having no curvature, so they neither register as maxima nor
    /// interrupt the detection of nearby finite peaks.
    ///
    /// A curve with no curvature variation (a straight line, or fully coincident control points)
    /// has no isolated maximum and yields an empty vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A symmetric arch has a single curvature peak at the apex, t = 0.5.
    /// let curve = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 1.0),
    ///     Point2::new(2.0, 1.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let maxima = curve.find_curvature_maxima();
    /// assert_eq!(maxima.len(), 1);
    /// assert_relative_eq!(maxima[0].t, 0.5, epsilon = 1e-6);
    /// assert_relative_eq!(maxima[0].value, 2.0 / 3.0, epsilon = 1e-9);
    /// ```
    pub fn find_curvature_maxima(&self) -> Vec<SplineValue<f64>> {
        // Curvature undefined (cusp) is treated as "no curvature" so cusps neither read as maxima
        // nor block detection of finite peaks beside them, matching `find_max_curvature`.
        let kappa = |t: f64| {
            let k = self.curvature(t);
            if k.is_finite() { k } else { f64::NEG_INFINITY }
        };

        // Coarse uniform scan: a cubic has only a handful of curvature maxima, so a few hundred
        // samples is comfortably fine enough to isolate each peak in its own bracket.
        const SAMPLES: usize = 256;
        let step = 1.0 / SAMPLES as f64;
        let samples: Vec<f64> = (0..=SAMPLES).map(|i| kappa(i as f64 * step)).collect();

        let mut result = Vec::new();
        for i in 0..=SAMPLES {
            let here = samples[i];
            if !here.is_finite() {
                continue;
            }

            // A maximum must not be lower than either neighbor. The asymmetric comparison
            // (`>` left, `>=` right) counts a flat-topped run of equal samples exactly once, at its
            // left edge. Out-of-domain neighbors are vacuously satisfied so the endpoints can
            // register, but we additionally require a strict rise over at least one real neighbor,
            // so a perfectly flat curvature profile (a straight line, κ ≡ 0) reports nothing rather
            // than spuriously flagging an endpoint.
            let has_left = i > 0;
            let has_right = i < SAMPLES;
            let left_ok = !has_left || here > samples[i - 1];
            let right_ok = !has_right || here >= samples[i + 1];
            let strictly_above_neighbor =
                (has_left && here > samples[i - 1]) || (has_right && here > samples[i + 1]);
            if !(left_ok && right_ok && strictly_above_neighbor) {
                continue;
            }

            let t_i = i as f64 * step;
            let a = (t_i - step).max(0.0);
            let b = (t_i + step).min(1.0);
            let t = golden_section_max(kappa, a, b);

            // The refined interior point can score below the bracketing sample in degenerate cases;
            // keep whichever parameter actually yields the larger curvature.
            let t = if kappa(t) >= here { t } else { t_i };
            result.push(SplineValue::new(t, self.curvature(t)));
        }

        result
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
        // TODO: unify this with whatever ends up happening with AABBs
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

    /// Returns the arc length of the curve over the parameter range `[t0, t1]`, computed by
    /// numerical integration of the speed `|B'(t)|` with composite Gauss-Legendre quadrature.
    ///
    /// The arc length of a parametric curve is the integral of its speed:
    ///
    /// `L = ∫ |B'(t)| dt`
    ///
    /// For a cubic Bezier the speed is the square root of a quartic polynomial in `t`, which has
    /// no elementary antiderivative, so the integral is evaluated numerically rather than in
    /// closed form. The interval `[t0, t1]` is split into [`ARC_LENGTH_PANELS`] equal panels and
    /// each panel is integrated with the fixed 10-point Gauss-Legendre rule in
    /// [`GAUSS_LEGENDRE_10`]. A 10-point rule is exact for polynomials up to degree 19, so for the
    /// smooth speed function of a single cubic segment the per-panel error is negligible and the
    /// result is accurate to near machine precision.
    ///
    /// The parameters are not restricted to `[0, 1]`; a wider range integrates over the
    /// extrapolated polynomial. The integral is signed by direction, so when `t1 < t0` the result
    /// is the negative of the length traversed; pass the parameters in increasing order for a
    /// conventional non-negative length.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// // A straight cubic spanning (0,0) to (3,0): its arc length is just the chord length.
    /// let line = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(2.0, 0.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// assert_relative_eq!(line.arc_length_between(0.0, 1.0), 3.0, epsilon = 1e-9);
    /// ```
    pub fn arc_length_between(&self, t0: f64, t1: f64) -> f64 {
        let h = (t1 - t0) / ARC_LENGTH_PANELS as f64;
        let mut total = 0.0;
        for i in 0..ARC_LENGTH_PANELS {
            let a = t0 + i as f64 * h;
            let mid = a + 0.5 * h;
            let half = 0.5 * h;
            let mut panel = 0.0;
            for &(x, w) in &GAUSS_LEGENDRE_10 {
                let t = mid + half * x;
                panel += w * self.derivative(t).norm();
            }
            total += panel * half;
        }
        total
    }

    /// Returns the total arc length of the curve over the parameter range `[0, 1]`.
    ///
    /// This is a convenience wrapper over [`arc_length_between`](Self::arc_length_between) for the
    /// full curve; see that method for the integration details.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    /// use approx::assert_relative_eq;
    ///
    /// let line = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(2.0, 0.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// assert_relative_eq!(line.arc_length(), 3.0, epsilon = 1e-9);
    /// ```
    pub fn arc_length(&self) -> f64 {
        self.arc_length_between(0.0, 1.0)
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

    /// Consumes the curve and builds its [`CubicSplineQueries`] acceleration structure, ready for
    /// repeated spatial queries such as closest-point projection.
    ///
    /// This is a convenience equivalent to `CubicSplineQueries::from(self)`; prefer whichever reads
    /// better at the call site. Building the structure is a one-time cost paid up front so that
    /// individual queries are cheap, so reuse the returned value across many queries rather than
    /// rebuilding it per query.
    pub fn into_query(self) -> CubicSplineQueries<D> {
        CubicSplineQueries::new(self)
    }
}

/// Number of equal panels the arc-length integral is split into before the Gauss-Legendre rule is
/// applied to each. Composing several panels keeps the per-panel error of the fixed-order rule
/// negligible even for the curviest single cubic segments.
const ARC_LENGTH_PANELS: usize = 16;

/// Nodes and weights of the 10-point Gauss-Legendre quadrature rule on the reference interval
/// `[-1, 1]`, stored as `(node, weight)` pairs. The rule integrates polynomials up to degree 19
/// exactly; used by [`CubicSpline::arc_length_between`] to integrate the curve's speed.
const GAUSS_LEGENDRE_10: [(f64, f64); 10] = [
    (-0.973_906_528_517_171_7, 0.066_671_344_308_688_1),
    (-0.865_063_366_688_984_5, 0.149_451_349_150_580_6),
    (-0.679_409_568_299_024_4, 0.219_086_362_515_982_0),
    (-0.433_395_394_129_247_2, 0.269_266_719_309_996_3),
    (-0.148_874_338_981_631_2, 0.295_524_224_714_752_9),
    (0.148_874_338_981_631_2, 0.295_524_224_714_752_9),
    (0.433_395_394_129_247_2, 0.269_266_719_309_996_3),
    (0.679_409_568_299_024_4, 0.219_086_362_515_982_0),
    (0.865_063_366_688_984_5, 0.149_451_349_150_580_6),
    (0.973_906_528_517_171_7, 0.066_671_344_308_688_1),
];

/// Golden-section search for the location of the maximum of `f` over the bracket `[a, b]`,
/// assuming a single maximum within it. Returns the parameter at which the maximum is attained.
///
/// The caller is responsible for supplying a bracket that contains exactly one maximum (typically
/// the span between the samples on either side of a coarse-scan peak); the search converges to that
/// peak to near machine precision.
fn golden_section_max(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    const INV_PHI: f64 = 0.618_033_988_749_894_8;
    let mut c = b - (b - a) * INV_PHI;
    let mut d = a + (b - a) * INV_PHI;
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..100 {
        if (b - a) <= 1e-15 {
            break;
        }
        if fc >= fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) * INV_PHI;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) * INV_PHI;
            fd = f(d);
        }
    }
    0.5 * (a + b)
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
    fn find_cusp_smooth_curve_returns_none() {
        assert!(sample_2d().find_cusp().is_none());
    }

    #[test]
    fn find_cusp_p1_eq_p2_loop_at_half() {
        // p0 == p3 and p1 == p2 with mirrored arms: B'(0.5) = 0 exactly.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 0.0),
        );
        let t = c.find_cusp().expect("expected a cusp at 0.5");
        assert_relative_eq!(t, 0.5, epsilon = 1e-12);
        // Sanity: derivative at the reported cusp parameter is (numerically) zero.
        assert_relative_eq!(c.derivative(t).norm(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn find_cusp_at_start_when_p0_eq_p1() {
        // p0 == p1 forces B'(0) = 0 and so a cusp at t = 0.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
        );
        let t = c.find_cusp().expect("expected a cusp at 0");
        assert_relative_eq!(t, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn find_cusp_at_end_when_p2_eq_p3() {
        // p2 == p3 forces B'(1) = 0 and so a cusp at t = 1.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 0.0),
        );
        let t = c.find_cusp().expect("expected a cusp at 1");
        assert_relative_eq!(t, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn find_inflections_symmetric_s_curve() {
        // Anti-symmetric in y about t = 0.5: a single inflection right at the midpoint.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, -1.0),
            Point2::new(3.0, 0.0),
        );
        let infl = c.find_inflections();
        assert_relative_eq!(infl[0], 0.5, epsilon = 1e-12);
        assert!(infl[1].is_nan());
    }

    #[test]
    fn find_inflections_sample_2d_has_none() {
        // sample_2d bulges in a single direction across t in [0, 1]: no inflections.
        let c = sample_2d();
        let infl = c.find_inflections();
        assert!(infl[0].is_nan());
        assert!(infl[1].is_nan());
    }

    #[test]
    fn find_inflections_straight_line_returns_nan() {
        // Collinear control points: cross product is identically zero.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
        );
        let infl = c.find_inflections();
        assert!(infl[0].is_nan());
        assert!(infl[1].is_nan());
    }

    #[test]
    fn find_inflections_signed_curvature_changes_sign() {
        // Build a curve with one inflection in (0, 1) and verify that the signed curvature
        // (B' × B'') flips sign across the reported parameter.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(2.0, -1.0),
            Point2::new(3.0, 1.0),
        );
        let infl = c.find_inflections();
        let signed = |t: f64| {
            let d1 = c.derivative(t);
            let d2 = c.second_derivative(t);
            d1.x * d2.y - d1.y * d2.x
        };
        let mut found_in_range = false;
        for &t in &infl {
            if (0.0..=1.0).contains(&t) {
                found_in_range = true;
                assert_relative_eq!(signed(t), 0.0, epsilon = 1e-10);
                let h = 1e-3;
                assert!(
                    signed(t - h) * signed(t + h) < 0.0,
                    "no sign change around t = {}",
                    t
                );
            }
        }
        assert!(found_in_range, "expected an inflection in [0, 1]");
    }

    #[test]
    fn find_cusp_all_coincident_returns_none() {
        let p = Point2::new(1.5, 2.5);
        let c = CubicSpline::new(p, p, p, p);
        assert!(c.find_cusp().is_none());
    }

    #[test]
    fn find_curvature_zeros_3d_planar_s_curve() {
        // 2D S-shape lifted into the xy-plane of 3D: inflection at t = 0.5.
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, -1.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let zeros = c.find_curvature_zeros();
        assert_relative_eq!(zeros[0], 0.5, epsilon = 1e-12);
        assert!(zeros[1].is_nan());
    }

    #[test]
    fn find_curvature_zeros_3d_straight_line_returns_nan() {
        // Curvature is identically zero, so there is no isolated zero to report.
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let zeros = c.find_curvature_zeros();
        assert!(zeros[0].is_nan());
        assert!(zeros[1].is_nan());
    }

    #[test]
    fn find_curvature_zeros_2d_matches_find_inflections_in_range() {
        // For 2D curves, find_curvature_zeros must agree with find_inflections on any roots
        // that lie inside [0, 1].
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, -1.0),
            Point2::new(3.0, 0.0),
        );
        let generic = c.find_curvature_zeros();
        let specific = c.find_inflections();
        assert_relative_eq!(generic[0], specific[0], epsilon = 1e-12);
        assert!(generic[1].is_nan() && specific[1].is_nan());
    }

    #[test]
    fn find_curvature_zeros_generic_3d_curve_has_none() {
        // A non-planar 3D curve whose per-pair candidate roots disagree across pairs: the
        // verification step must reject them all.
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.5),
            Point3::new(2.0, 2.0, 1.5),
            Point3::new(3.0, 0.0, 2.0),
        );
        let zeros = c.find_curvature_zeros();
        assert!(zeros[0].is_nan());
        assert!(zeros[1].is_nan());
    }

    #[test]
    fn find_curvature_zeros_sample_2d_has_none() {
        let zeros = sample_2d().find_curvature_zeros();
        assert!(zeros[0].is_nan());
        assert!(zeros[1].is_nan());
    }

    #[test]
    fn find_max_curvature_symmetric_arch() {
        // sample_2d is symmetric about t = 0.5 and bulges in one direction; curvature peaks at
        // the apex. By hand at t = 0.5: B' = (3, 0), B'' = (0, -6), κ = 18 / 27 = 2/3.
        let c = sample_2d();
        let max = c.find_max_curvature();
        assert_relative_eq!(max.t, 0.5, epsilon = 1e-6);
        assert_relative_eq!(c.position(max.t), Point2::new(1.5, 0.75), epsilon = 1e-6);
        assert_relative_eq!(max.value, 2.0 / 3.0, epsilon = 1e-9);
    }

    #[test]
    fn find_max_curvature_agrees_with_dense_sampling() {
        // Cross-check the analytic refinement against a brute-force scan on an asymmetric curve.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 2.0),
            Point2::new(2.5, 1.5),
            Point2::new(3.0, 0.0),
        );
        let max = c.find_max_curvature();

        let mut brute_t = 0.0;
        let mut brute_k = f64::NEG_INFINITY;
        for i in 0..=100_000 {
            let s = i as f64 / 100_000.0;
            let ks = c.curvature(s);
            if ks.is_finite() && ks > brute_k {
                brute_k = ks;
                brute_t = s;
            }
        }
        assert_relative_eq!(max.t, brute_t, epsilon = 1e-4);
        assert!(
            max.value >= brute_k - 1e-9,
            "refined {} < brute {}",
            max.value,
            brute_k
        );
    }

    #[test]
    fn find_max_curvature_straight_line_is_zero() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let max = c.find_max_curvature();
        assert_relative_eq!(max.value, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn find_max_curvature_all_coincident_is_nan() {
        let p = Point2::new(1.5, 2.5);
        let c = CubicSpline::new(p, p, p, p);
        let max = c.find_max_curvature();
        assert_eq!(max.t, 0.0);
        assert!(max.value.is_nan());
    }

    #[test]
    fn find_curvature_maxima_single_peak() {
        // sample_2d has one curvature peak, at the apex t = 0.5 with κ = 2/3.
        let c = sample_2d();
        let maxima = c.find_curvature_maxima();
        assert_eq!(maxima.len(), 1);
        assert_relative_eq!(maxima[0].t, 0.5, epsilon = 1e-6);
        assert_relative_eq!(maxima[0].value, 2.0 / 3.0, epsilon = 1e-9);
    }

    #[test]
    fn find_curvature_maxima_two_peaks_on_s_curve() {
        // A symmetric S-curve has an inflection (zero curvature) at t = 0.5 and a curvature peak
        // on each lobe, mirrored about the midpoint.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, -1.0),
            Point2::new(3.0, 0.0),
        );
        let maxima = c.find_curvature_maxima();
        assert_eq!(maxima.len(), 2);
        // Ordered by ascending t, mirrored about 0.5, with equal curvature by symmetry.
        assert!(maxima[0].t < 0.5 && maxima[1].t > 0.5);
        assert_relative_eq!(maxima[0].t, 1.0 - maxima[1].t, epsilon = 1e-6);
        assert_relative_eq!(maxima[0].value, maxima[1].value, epsilon = 1e-9);
    }

    #[test]
    fn find_curvature_maxima_are_local_maxima_and_include_global() {
        // Every reported maximum must actually be a local maximum of the curvature, and the global
        // maximum from find_max_curvature must appear among them.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 2.0),
            Point2::new(2.5, 1.5),
            Point2::new(3.0, 0.0),
        );
        let maxima = c.find_curvature_maxima();
        assert!(!maxima.is_empty());

        // Ascending in t.
        for w in maxima.windows(2) {
            assert!(w[0].t < w[1].t);
        }

        // Each is at least as curved as nearby parameters (clamped to the domain).
        let h = 1e-4;
        for m in &maxima {
            let lo = (m.t - h).max(0.0);
            let hi = (m.t + h).min(1.0);
            assert!(m.value >= c.curvature(lo) - 1e-6);
            assert!(m.value >= c.curvature(hi) - 1e-6);
        }

        // The global peak is one of the reported maxima.
        let global = c.find_max_curvature();
        let best = maxima
            .iter()
            .map(|m| m.value)
            .fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(best, global.value, epsilon = 1e-9);
    }

    #[test]
    fn find_curvature_maxima_straight_line_is_empty() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        assert!(c.find_curvature_maxima().is_empty());
    }

    #[test]
    fn find_curvature_maxima_all_coincident_is_empty() {
        let p = Point2::new(1.5, 2.5);
        let c = CubicSpline::new(p, p, p, p);
        assert!(c.find_curvature_maxima().is_empty());
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
        // A loop where p0 == p3: the degenerate-chord branch of chord_perp_distance must trigger
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

    /// Reference arc length via a dense composite midpoint integration of the speed, used to
    /// cross-check the Gauss-Legendre quadrature.
    fn arc_length_reference<const D: usize>(c: &CubicSpline<D>, t0: f64, t1: f64) -> f64 {
        const N: usize = 2_000_000;
        let h = (t1 - t0) / N as f64;
        let mut total = 0.0;
        for i in 0..N {
            let t = t0 + (i as f64 + 0.5) * h;
            total += c.derivative(t).norm();
        }
        total * h
    }

    #[test]
    fn arc_length_straight_line_is_chord_length() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        assert_relative_eq!(c.arc_length(), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn arc_length_matches_dense_integration() {
        let c = sample_2d();
        let reference = arc_length_reference(&c, 0.0, 1.0);
        assert_relative_eq!(c.arc_length(), reference, epsilon = 1e-9);
    }

    #[test]
    fn arc_length_between_is_additive() {
        let c = sample_2d();
        let whole = c.arc_length();
        let split = c.arc_length_between(0.0, 0.37) + c.arc_length_between(0.37, 1.0);
        assert_relative_eq!(whole, split, epsilon = 1e-12);
    }

    #[test]
    fn arc_length_between_subrange_matches_reference() {
        let c = sample_2d();
        let reference = arc_length_reference(&c, 0.2, 0.8);
        assert_relative_eq!(c.arc_length_between(0.2, 0.8), reference, epsilon = 1e-9);
    }

    #[test]
    fn arc_length_between_reversed_is_negated() {
        let c = sample_2d();
        assert_relative_eq!(
            c.arc_length_between(0.8, 0.2),
            -c.arc_length_between(0.2, 0.8),
            epsilon = 1e-12
        );
    }

    #[test]
    fn arc_length_matches_fine_polyline() {
        // The arc length must agree with the length of a finely subdivided polyline approximation.
        let c = sample_2d();
        let pts = c.polyline(1e-7);
        let poly_len: f64 = pts
            .windows(2)
            .map(|w| (w[1].coords - w[0].coords).norm())
            .sum();
        assert_relative_eq!(c.arc_length(), poly_len, epsilon = 1e-5);
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
