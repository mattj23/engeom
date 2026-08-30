//! This module generalizes fitting of a cubic Bézier spline to sample points.
//!
//! The mechanism mirrors the boundary fitting in `geom2::boundary2`: the caller supplies a builder
//! closure that turns a parameter vector into a [`CubicSpline`], an initial guess for those
//! parameters, and the sample points to fit against. A Levenberg-Marquardt minimization then drives
//! the parameters so that the spline best fits the points in a least-squares sense, where each
//! residual is the distance from a point to its closest point on the spline.
//!
//! Writing the builder so that the parameterization inherently contains the desired constraints is
//! the easiest way to get a well-behaved fit. Because the residuals are closest-point distances,
//! there is no inherent pressure preventing the curve from sliding or growing past the points
//! unless the parameterization forbids it.
//!
//! # Jacobian
//!
//! The Jacobian is analytic with respect to the control points and uses the envelope theorem. At
//! the closest point `p(t*)` to a sample `x`, the closest-point parameter `t*` is stationary (or
//! clamped to an end of the curve), so the derivative of the distance `|p(t*) - x|` with respect to
//! control point `q_j` is just `B_j(t*) * u`, where `B_j` is the cubic Bernstein weight and `u` is
//! the unit vector from `x` to `p(t*)`. This does not require the points to be projected again. The
//! derivative with respect to the builder's parameters is then obtained through the chain rule,
//! using forward differences for the control points. This is inexpensive because it requires no
//! projections. The result is a damped tangent-distance minimization in the taxonomy of Wang,
//! Pottmann, and Liu (2006), with the damping supplied by Levenberg-Marquardt.

use super::{CubicSpline, SplineValue, bernstein_basis};
use crate::Result;
use crate::na::{DMatrix, DVector, Dyn, Matrix, Owned, Point, U1, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

/// A builder closure that turns a parameter vector into a [`CubicSpline`] of dimension `D`. The
/// closure may return an `Err` if the parameters describe invalid geometry, but it must succeed on
/// the initial guess or the fitting routine will abort before it starts.
pub type SplineBuildFn<const D: usize> = Box<dyn Fn(&DVector<f64>) -> Result<CubicSpline<D>>>;

/// Forward-difference step used to differentiate the builder's control points with respect to its
/// parameters.
const DELTA: f64 = 1e-6;

/// Below this closest-point distance, the direction from a sample point to the curve is undefined,
/// so the point's Jacobian row is zeroed. The residual itself is already zero to this precision,
/// so the point contributes nothing to the step either way.
const MIN_DIST: f64 = 1e-12;

/// A Jacobian column whose norm is below this fraction of the largest column norm is numerically
/// zero and is set to exactly zero. The Levenberg-Marquardt solver scales each parameter by its
/// column norm and treats only an exactly zero column as unscaled; a roundoff-level column would
/// otherwise get a scale of ~1e-16 and produce an enormous, meaningless step along it. This
/// happens whenever the current spline is insensitive to a parameter, such as with a collinear
/// initial guess, where moving a control point along the line changes no distance.
const ZERO_COLUMN_RATIO: f64 = 1e-12;

/// The result of a successful spline fit: the fitted spline, the optimized parameters that produce
/// it, and the per-point residuals (unweighted closest-point distances) at those parameters.
#[derive(Clone)]
pub struct SplineFitResult<const D: usize> {
    pub spline: CubicSpline<D>,
    pub params: DVector<f64>,
    pub residuals: DVector<f64>,
}

impl<const D: usize> SplineFitResult<D> {
    pub fn new(spline: CubicSpline<D>, params: DVector<f64>, residuals: DVector<f64>) -> Self {
        Self {
            spline,
            params,
            residuals,
        }
    }
}

/// Given a set of points, a function which takes a parameter vector and produces a [`CubicSpline`],
/// and an initial guess, this function will attempt to perform a Levenberg-Marquardt best fit of
/// the parameters to produce a `CubicSpline` that minimizes the residuals of the points projected
/// onto the curve.
///
/// You will need to provide a `DVector` with a reasonable initial guess, and a [`SplineBuildFn`]
/// that accepts a `DVector` of that size and produces a spline. The function is allowed to return
/// an `Err(Box<dyn Error>)` if the parameters produce invalid geometry, but the initial guess
/// values provided may not fail or this function will exit before attempting the minimization.
///
/// Because the residuals are the distances from the points to the nearest point on the curve, be
/// aware that there is no inherent pressure preventing the curve geometry from sliding or growing
/// beyond the sample points _unless_ it exists within the nature of the parameterization. For this
/// reason it is best to write the builder function so that any constraints are inherently contained
/// within it.
///
/// This is the unweighted form of [`fit_spline_to_points_weighted`].
///
/// # Arguments
///
/// * `points`: A slice of `Point<f64, D>` that the spline will be fit to. The residuals are
///   calculated as the distance from each point to the closest point on the curve, so the points
///   will "pull" the closest area of the curve towards themselves, but won't have any effect on the
///   curve farther away.
/// * `builder`: A function that takes a `DVector` of parameters and produces a `CubicSpline`.
/// * `initial`: An initial guess for the parameters. The `builder` function must return an
///   `Ok(CubicSpline)` with this input for the fitting algorithm to initialize. The algorithm will
///   use the size of the initial guess to determine the number of parameters, so make sure it
///   matches what the `builder` function is expecting.
///
/// returns: Result<SplineFitResult<D>>
///
/// # Examples
///
/// ```
/// // Fit a cubic spline to points sampled along a known curve, recovering its inner control
/// // points while holding the endpoints fixed.
/// // ------------------------------------------------------------------------------------------
/// use engeom::{Point2, DVector};
/// use engeom::common::cubic_spline::{CubicSpline, SplineBuildFn, fit_spline_to_points};
/// use approx::assert_relative_eq;
///
/// // The curve we are trying to recover, and a set of points sampled along it.
/// let truth = CubicSpline::new(
///     Point2::new(0.0, 0.0),
///     Point2::new(1.0, 2.0),
///     Point2::new(2.0, 2.0),
///     Point2::new(3.0, 0.0),
/// );
/// let points: Vec<Point2> = (0..=20)
///     .map(|i| truth.position(i as f64 / 20.0))
///     .collect();
///
/// // The endpoints are held fixed; only the two interior control points are free, encoded as
/// // x,y pairs in the parameter vector.
/// let builder: SplineBuildFn<2> = Box::new(|p: &DVector| {
///     Ok(CubicSpline::new(
///         Point2::new(0.0, 0.0),
///         Point2::new(p[0], p[1]),
///         Point2::new(p[2], p[3]),
///         Point2::new(3.0, 0.0),
///     ))
/// });
///
/// let initial = DVector::from(vec![1.0, 0.0, 2.0, 0.0]);
/// let result = fit_spline_to_points(&points, &builder, initial).unwrap();
///
/// let expected = DVector::from(vec![1.0, 2.0, 2.0, 2.0]);
/// assert_relative_eq!(result.params, expected, epsilon = 1.0e-6);
///
/// // The fitted spline is returned directly, there is no need to call the builder again.
/// assert_relative_eq!(result.spline.p1, truth.p1, epsilon = 1.0e-6);
/// ```
pub fn fit_spline_to_points<const D: usize>(
    points: &[Point<f64, D>],
    builder: &SplineBuildFn<D>,
    initial: DVector<f64>,
) -> Result<SplineFitResult<D>> {
    fit_spline_to_points_weighted(points, None, builder, initial)
}

/// The weighted form of [`fit_spline_to_points`]. Each point's residual and Jacobian row are scaled
/// by a fixed per-point weight, so points with larger weights pull harder on the curve and points
/// with zero weight are ignored entirely. Weights remain fixed throughout the solve, making this
/// suitable as the refinement step in an outer robust loop.
///
/// # Arguments
///
/// * `points`: The points the spline will be fit to; see [`fit_spline_to_points`].
/// * `weights`: If `Some`, a slice the same length as `points` of non-negative weights. `None` is
///   equivalent to all weights being one.
/// * `builder`: A function that takes a `DVector` of parameters and produces a `CubicSpline`.
/// * `initial`: An initial guess for the parameters; see [`fit_spline_to_points`].
///
/// returns: `Result<SplineFitResult<D>>`, whose `residuals` are the unweighted closest-point
/// distances.
///
/// # Failure
///
/// Returns an error if the number of weights differs from the number of points, the builder rejects
/// the initial parameters, or the solver does not converge successfully.
pub fn fit_spline_to_points_weighted<const D: usize>(
    points: &[Point<f64, D>],
    weights: Option<&[f64]>,
    builder: &SplineBuildFn<D>,
    initial: DVector<f64>,
) -> Result<SplineFitResult<D>> {
    let weights = match weights {
        Some(w) if w.len() != points.len() => {
            return Err(format!(
                "expected {} weights for {} points, got {}",
                points.len(),
                points.len(),
                w.len()
            )
            .into());
        }
        Some(w) => DVector::from_column_slice(w),
        None => DVector::from_element(points.len(), 1.0),
    };

    let problem = SplineFit::try_new(points, weights, builder, initial)?;
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let state = result.state.unwrap();
        let residuals = DVector::from_iterator(points.len(), state.feet.iter().map(|f| f.value));
        Ok(SplineFitResult::new(state.spline, result.params, residuals))
    } else {
        Err(format!("Fitting failed: {:?}", report.termination).into())
    }
}

/// The spline produced by the current parameters together with the closest-point projection of
/// every sample point onto it.
struct FitState<const D: usize> {
    spline: CubicSpline<D>,
    /// For each sample point, the parameter of its closest point on the spline and the distance to
    /// that point.
    feet: Vec<SplineValue<f64>>,
}

impl<const D: usize> FitState<D> {
    fn new(spline: CubicSpline<D>, points: &[Point<f64, D>]) -> Self {
        let queries = spline.into_query();
        let feet = points.iter().map(|p| queries.project_point(p)).collect();
        Self {
            spline: queries.into_spline(),
            feet,
        }
    }
}

/// The Levenberg-Marquardt problem backing [`fit_spline_to_points_weighted`]. Each residual is the
/// weighted distance from a sample point to its closest point on the spline built from the current
/// parameters. See the module documentation for the Jacobian.
struct SplineFit<'a, const D: usize> {
    points: &'a [Point<f64, D>],
    weights: DVector<f64>,
    builder: &'a SplineBuildFn<D>,
    params: DVector<f64>,
    /// `None` when the builder rejected the current parameters.
    state: Option<FitState<D>>,
}

impl<'a, const D: usize> SplineFit<'a, D> {
    fn try_new(
        points: &'a [Point<f64, D>],
        weights: DVector<f64>,
        builder: &'a SplineBuildFn<D>,
        initial: DVector<f64>,
    ) -> Result<Self> {
        // Check to make sure that the initial value doesn't fail
        let _ = builder(&initial)?;

        let mut problem = SplineFit {
            points,
            weights,
            builder,
            params: DVector::zeros(initial.len()),
            state: None,
        };

        problem.set_params(&initial);

        Ok(problem)
    }

    /// The derivative of the `4 * D` control-point coordinates with respect to the builder's
    /// parameters, computed using forward differences. A parameter that the builder rejects when
    /// disturbed receives a zero column.
    fn control_jacobian(&self, spline: &CubicSpline<D>) -> DMatrix<f64> {
        let base = control_coords(spline);
        let mut dctrl = DMatrix::zeros(4 * D, self.params.len());

        for k in 0..self.params.len() {
            let mut params = self.params.clone();
            params[k] += DELTA;
            let Ok(disturbed) = (self.builder)(&params) else {
                continue;
            };
            let c = control_coords(&disturbed);
            for r in 0..4 * D {
                dctrl[(r, k)] = (c[r] - base[r]) / DELTA;
            }
        }

        dctrl
    }
}

/// Sets every column of `jac` whose norm is below [`ZERO_COLUMN_RATIO`] times the largest column
/// norm to exactly zero.
fn zero_negligible_columns(jac: &mut DMatrix<f64>) {
    let norms: Vec<f64> = jac.column_iter().map(|c| c.norm()).collect();
    let largest = norms.iter().copied().fold(0.0, f64::max);
    for (k, norm) in norms.iter().enumerate() {
        if *norm < ZERO_COLUMN_RATIO * largest {
            jac.column_mut(k).fill(0.0);
        }
    }
}

/// The control point coordinates of a spline flattened in `p0, p1, p2, p3` order, `D` values each.
fn control_coords<const D: usize>(spline: &CubicSpline<D>) -> Vec<f64> {
    [&spline.p0, &spline.p1, &spline.p2, &spline.p3]
        .iter()
        .flat_map(|p| p.coords.iter().copied())
        .collect()
}

impl<const D: usize> LeastSquaresProblem<f64, Dyn, Dyn> for SplineFit<'_, D> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params = x.clone();
        self.state = (self.builder)(&self.params)
            .ok()
            .map(|spline| FitState::new(spline, self.points));
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let state = self.state.as_ref()?;
        Some(DVector::from_iterator(
            self.points.len(),
            state
                .feet
                .iter()
                .zip(self.weights.iter())
                .map(|(f, w)| f.value * w),
        ))
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let state = self.state.as_ref()?;
        let n_params = self.params.len();
        let dctrl = self.control_jacobian(&state.spline);

        let mut jac =
            Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(self.points.len(), n_params);

        for (i, (point, foot)) in self.points.iter().zip(state.feet.iter()).enumerate() {
            if foot.value < MIN_DIST {
                continue;
            }

            // Unit vector from the sample point to its closest point on the curve
            let u = (state.spline.position(foot.t) - point) / foot.value;
            let b = bernstein_basis(foot.t);
            let w = self.weights[i];

            for k in 0..n_params {
                let mut sum = 0.0;
                for (j, bj) in b.iter().enumerate() {
                    for (a, ua) in u.iter().enumerate() {
                        sum += bj * ua * dctrl[(j * D + a, k)];
                    }
                }
                jac[(i, k)] = w * sum;
            }
        }

        zero_negligible_columns(&mut jac);
        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::{RandomGeometry2, RandomGeometry3};
    use crate::{Point2, Point3, Vector2};
    use approx::assert_relative_eq;

    fn truth() -> CubicSpline<2> {
        CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(2.0, 2.0),
            Point2::new(3.0, 0.0),
        )
    }

    /// Builder with the endpoints fixed and the two interior control points free.
    fn inner_builder() -> SplineBuildFn<2> {
        Box::new(|p: &DVector<f64>| {
            Ok(CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(p[0], p[1]),
                Point2::new(p[2], p[3]),
                Point2::new(3.0, 0.0),
            ))
        })
    }

    /// A nonlinear builder with fixed endpoints and interior control points defined by an angle and
    /// an arm length from each end. This exercises the chain rule through the builder.
    fn polar_builder() -> SplineBuildFn<2> {
        Box::new(|p: &DVector<f64>| {
            let p0 = Point2::new(0.0, 0.0);
            let p3 = Point2::new(3.0, 0.0);
            let p1 = p0 + Vector2::new(p[0].cos(), p[0].sin()) * p[1];
            let p2 = p3 + Vector2::new(p[2].cos(), p[2].sin()) * p[3];
            Ok(CubicSpline::new(p0, p1, p2, p3))
        })
    }

    /// A finite-difference Jacobian computed directly by disturbing each parameter, rebuilding the
    /// spline, and projecting every point again.
    fn numeric_jacobian<const D: usize>(fit: &SplineFit<D>) -> DMatrix<f64> {
        let base = fit.residuals().unwrap();
        let mut jac = DMatrix::zeros(fit.points.len(), fit.params.len());
        let h = 1e-7;
        for k in 0..fit.params.len() {
            let mut params = fit.params.clone();
            params[k] += h;
            let spline = (fit.builder)(&params).unwrap();
            let state = FitState::new(spline, fit.points);
            for i in 0..fit.points.len() {
                jac[(i, k)] = (state.feet[i].value * fit.weights[i] - base[i]) / h;
            }
        }
        jac
    }

    fn assert_jacobians_agree<const D: usize>(fit: &SplineFit<D>) {
        let analytic = fit.jacobian().unwrap();
        let numeric = numeric_jacobian(fit);
        for i in 0..analytic.nrows() {
            for k in 0..analytic.ncols() {
                assert_relative_eq!(
                    analytic[(i, k)],
                    numeric[(i, k)],
                    epsilon = 1e-5,
                    max_relative = 1e-4
                );
            }
        }
    }

    #[test]
    fn recover_inner_control_points() {
        let curve = truth();
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();

        let initial = DVector::from(vec![1.0, 0.0, 2.0, 0.0]);
        let result = fit_spline_to_points(&points, &inner_builder(), initial).unwrap();

        let expected = DVector::from(vec![1.0, 2.0, 2.0, 2.0]);
        assert_relative_eq!(result.params, expected, epsilon = 1e-6);
        assert_relative_eq!(result.spline.p1, curve.p1, epsilon = 1e-6);
        assert_relative_eq!(result.spline.p2, curve.p2, epsilon = 1e-6);
    }

    #[test]
    fn residuals_are_small_after_fit() {
        let curve = truth();
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();

        let initial = DVector::from(vec![1.2, 0.5, 1.8, 0.5]);
        let result = fit_spline_to_points(&points, &inner_builder(), initial).unwrap();

        for r in result.residuals.iter() {
            assert!(r.abs() < 1e-6, "residual {} too large", r);
        }
    }

    #[test]
    fn initial_guess_that_fails_is_an_error() {
        let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
        let builder: SplineBuildFn<2> =
            Box::new(|_: &DVector<f64>| Err("always fails".to_string().into()));
        let initial = DVector::from(vec![0.0, 0.0, 0.0, 0.0]);
        assert!(fit_spline_to_points(&points, &builder, initial).is_err());
    }

    #[test]
    fn wrong_weight_count_is_an_error() {
        let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
        let initial = DVector::from(vec![0.0, 0.0, 0.0, 0.0]);
        let result =
            fit_spline_to_points_weighted(&points, Some(&[1.0]), &inner_builder(), initial);
        assert!(result.is_err());
    }

    #[test]
    fn analytic_jacobian_matches_finite_differences_2d() {
        let mut rg = RandomGeometry2::from_seed(11);
        for _ in 0..20 {
            // Use random noisy samples around the reference curve, kept off the curve so the
            // offset direction is well defined, and random parameters away from the optimum.
            let curve = truth();
            let points: Vec<Point2> = (0..=15)
                .map(|i| curve.position(i as f64 / 15.0) + rg.gaussian_vector::<2>(0.1))
                .collect();
            let weights: Vec<f64> = (0..points.len()).map(|_| rg.f64(0.2, 2.0)).collect();
            let params = DVector::from(vec![
                rg.f64(0.3, 1.3),
                rg.f64(1.0, 3.0),
                rg.f64(1.8, 2.8),
                rg.f64(1.0, 3.0),
            ]);

            let builder = polar_builder();
            let fit = SplineFit::try_new(
                &points,
                DVector::from_column_slice(&weights),
                &builder,
                params,
            )
            .unwrap();
            assert_jacobians_agree(&fit);
        }
    }

    #[test]
    fn analytic_jacobian_matches_finite_differences_3d() {
        let mut rg = RandomGeometry3::from_seed(7);
        let builder: SplineBuildFn<3> = Box::new(|p: &DVector<f64>| {
            Ok(CubicSpline::new(
                Point3::new(p[0], p[1], p[2]),
                Point3::new(p[3], p[4], p[5]),
                Point3::new(p[6], p[7], p[8]),
                Point3::new(p[9], p[10], p[11]),
            ))
        });

        for _ in 0..20 {
            let curve =
                CubicSpline::new(rg.point(2.0), rg.point(2.0), rg.point(2.0), rg.point(2.0));
            let points: Vec<Point3> = (0..=15)
                .map(|i| curve.position(i as f64 / 15.0) + rg.gaussian_vector::<3>(0.1))
                .collect();
            let params = DVector::from_iterator(
                12,
                control_coords(&curve).iter().map(|c| c + rg.f64_sym(0.2)),
            );

            let fit = SplineFit::try_new(
                &points,
                DVector::from_element(points.len(), 1.0),
                &builder,
                params,
            )
            .unwrap();
            assert_jacobians_agree(&fit);
        }
    }

    #[test]
    fn collinear_initial_guess_converges() {
        // A straight-line initial guess is degenerate: moving a control point along the line
        // changes no distance, so the Jacobian columns for the x coordinates contain only
        // roundoff. Those columns must be treated as exactly zero, or the solver takes a giant
        // step along them.
        let curve = truth();
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let initial = DVector::from(vec![1.0, 0.0, 2.0, 0.0]);
        let result = fit_spline_to_points(&points, &inner_builder(), initial).unwrap();
        assert_relative_eq!(result.spline.p1, curve.p1, epsilon = 1e-6);
        assert_relative_eq!(result.spline.p2, curve.p2, epsilon = 1e-6);
    }

    #[test]
    fn zero_weight_points_are_ignored() {
        let curve = truth();
        let mut points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let mut weights = vec![1.0; points.len()];

        // A gross outlier that would drag the fit away from the reference curve.
        points.push(Point2::new(1.5, 5.0));
        weights.push(0.0);

        let initial = DVector::from(vec![1.2, 0.5, 1.8, 0.5]);
        let result =
            fit_spline_to_points_weighted(&points, Some(&weights), &inner_builder(), initial)
                .unwrap();

        let expected = DVector::from(vec![1.0, 2.0, 2.0, 2.0]);
        assert_relative_eq!(result.params, expected, epsilon = 1e-6);
    }
}
