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

use super::{CubicSpline, CubicSplineQueries};
use crate::Result;
use crate::na::{DVector, Dyn, Matrix, Owned, Point, U1, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

/// A builder closure that turns a parameter vector into a [`CubicSpline`] of dimension `D`. The
/// closure may return an `Err` if the parameters describe invalid geometry, but it must succeed on
/// the initial guess or the fitting routine will abort before it starts.
pub type SplineBuildFn<const D: usize> = Box<dyn Fn(&DVector<f64>) -> Result<CubicSpline<D>>>;

const DELTA: f64 = 1e-6;

/// The result of a successful spline fit: the optimized parameters and the per-point residuals
/// (closest-point distances) at those parameters.
#[derive(Clone)]
pub struct SplineFitResult {
    pub params: DVector<f64>,
    pub residuals: DVector<f64>,
}

impl SplineFitResult {
    pub fn new(params: DVector<f64>, residuals: DVector<f64>) -> Self {
        Self { params, residuals }
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
/// returns: Result<SplineFitResult>
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
/// ```
pub fn fit_spline_to_points<const D: usize>(
    points: &[Point<f64, D>],
    builder: &SplineBuildFn<D>,
    initial: DVector<f64>,
) -> Result<SplineFitResult> {
    let problem = SplineFit::try_new(points, builder, initial)?;
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let residuals = result.residuals.unwrap();
        Ok(SplineFitResult::new(result.params, residuals))
    } else {
        Err(format!("Fitting failed: {:?}", report.termination).into())
    }
}

/// The Levenberg-Marquardt problem backing [`fit_spline_to_points`]. Each residual is the distance
/// from a sample point to its closest point on the spline built from the current parameters, and
/// the jacobian is computed by forward finite differences over the parameters.
struct SplineFit<'a, const D: usize> {
    points: &'a [Point<f64, D>],
    params: DVector<f64>,
    builder: &'a SplineBuildFn<D>,
    residuals: Option<DVector<f64>>,
}

impl<'a, const D: usize> SplineFit<'a, D> {
    fn try_new(
        points: &'a [Point<f64, D>],
        builder: &'a SplineBuildFn<D>,
        initial: DVector<f64>,
    ) -> Result<Self> {
        // Check to make sure that the initial value doesn't fail
        let _ = builder(&initial)?;

        let mut problem = SplineFit {
            points,
            params: DVector::zeros(initial.len()),
            builder,
            residuals: None,
        };

        problem.set_params(&initial);

        Ok(problem)
    }

    /// Projects every sample point onto `queries` and collects the closest-point distances.
    fn distances(&self, queries: &CubicSplineQueries<D>) -> DVector<f64> {
        let mut res = DVector::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = queries.project_point(&self.points[i]).distance;
        }
        res
    }
}

impl<const D: usize> LeastSquaresProblem<f64, Dyn, Dyn> for SplineFit<'_, D> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params = x.clone();
        self.residuals = (self.builder)(&self.params)
            .ok()
            .map(|spline| self.distances(&spline.into_query()));
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        self.residuals.clone()
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let residuals = self.residuals.as_ref()?;

        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            residuals.len(),
            self.params.len(),
        );

        for k in 0..self.params.len() {
            let mut params = self.params.clone();
            params[k] += DELTA;
            let Ok(disturbed) = (self.builder)(&params) else {
                continue;
            };
            let queries = disturbed.into_query();

            for i in 0..residuals.len() {
                let d = queries.project_point(&self.points[i]).distance;
                jac[(i, k)] = (d - residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point2;
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

    #[test]
    fn recover_inner_control_points() {
        let curve = truth();
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();

        let initial = DVector::from(vec![1.0, 0.0, 2.0, 0.0]);
        let result = fit_spline_to_points(&points, &inner_builder(), initial).unwrap();

        let expected = DVector::from(vec![1.0, 2.0, 2.0, 2.0]);
        assert_relative_eq!(result.params, expected, epsilon = 1e-6);
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
}
