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
use crate::common::PCoords;
use crate::common::points::dist;
use crate::common::svd_basis::SvdBasis;
use crate::na::{DMatrix, DVector, Dyn, Matrix, Owned, Point, SVector, U1, Unit, Vector};
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

// ================================================================================================
// Named fitting constructors
// ================================================================================================

impl<const D: usize> CubicSpline<D> {
    /// Fit a cubic Bézier curve to a set of points, holding the two endpoints fixed and solving for
    /// the two interior control points.
    ///
    /// Fixing the endpoints is what makes this problem well-posed: the residuals are closest-point
    /// distances to the curve, which do not change when the curve extends beyond the points, so a
    /// fit with free endpoints has nothing to stop the ends from sliding. Fixing the endpoints
    /// leaves the two interior control points as the parameters to solve for.
    ///
    /// The fit starts from the best of three closed-form seeds and then refines against the true
    /// closest-point distances with a weighted Levenberg-Marquardt minimization
    /// ([`fit_spline_to_points_weighted`]). The seeds are the straight chord between the endpoints,
    /// and two linear least-squares solutions using a chord-length parameterization of the points.
    /// One solution orders the points by their projection onto the chord, while the other uses the
    /// input order. Projection ordering handles unordered input for any curve that does not double
    /// back along its chord. Input ordering handles ordered points from curves that do double back,
    /// while the chord provides a fallback for all other cases. The seed with the lowest weighted
    /// RMS closest-point distance is refined, so the caller does not need to determine which case
    /// applies.
    ///
    /// The fit determines the *curve* accurately, but its control points are only as well determined
    /// as the shape allows. To first order, shifting both interior control points along the curve
    /// reparameterizes it without changing its shape. The objective is therefore very flat in that
    /// direction, and control points fitted to exact data can differ from the generating control
    /// points by a fraction of a percent even when the curves agree much more closely. Compare the
    /// curves, rather than their control points, when checking a fit.
    ///
    /// # Arguments
    ///
    /// * `points`: The points to fit the curve to, in any order. At least two are required.
    /// * `p0`: The start point of the curve, held fixed.
    /// * `p3`: The end point of the curve, held fixed. Must be distinct from `p0`.
    /// * `weights`: If `Some`, a slice the same length as `points` of non-negative weights that
    ///   scale each point's residual. Points with zero weight have no effect on the fit. `None`
    ///   weights all points equally.
    ///
    /// returns: A `Result` containing the fitted curve.
    ///
    /// # Failure
    ///
    /// Returns an error if there are fewer than two points, the endpoints coincide, the weight
    /// count does not match the point count, or the minimization fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    ///
    /// let truth = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(2.0, 2.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let points: Vec<Point2> = (0..=20).map(|i| truth.position(i as f64 / 20.0)).collect();
    ///
    /// let fitted = CubicSpline::from_fit_with_ends(&points, &truth.p0, &truth.p3, None).unwrap();
    ///
    /// // The fitted curve passes through the sampled points.
    /// let queries = fitted.into_query();
    /// for p in &points {
    ///     assert!(queries.project_point(p).value < 1e-6);
    /// }
    /// ```
    pub fn from_fit_with_ends(
        points: &[impl PCoords<D>],
        p0: &Point<f64, D>,
        p3: &Point<f64, D>,
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        check_fit_inputs(points.len(), weights)?;
        let chord = p3 - p0;
        if chord.norm() <= f64::EPSILON {
            return Err("the endpoints of the spline must be distinct".into());
        }

        let points: Vec<Point<f64, D>> = points.iter().map(|p| Point::from(p.coords())).collect();
        fit_with_ends(&points, p0, p3, weights)
    }

    /// Fit a cubic Bézier curve to a set of points whose endpoints are unknown, assuming that the
    /// curve runs from one end of the points' principal axis to the other.
    ///
    /// The principal axis is the direction of greatest variance of the points (the first vector of
    /// an [`SvdBasis`]). The points are ordered by their projection onto it and the two extreme
    /// points anchor the curve's ends: each endpoint is free to slide in the plane through its
    /// extreme point perpendicular to the axis, which absorbs the perpendicular noise of a single
    /// sample while still pinning the ends along the axis (the direction in which a closest-point
    /// fit could otherwise slide freely). The interior control points are seeded and refined as
    /// in [`Self::from_fit_with_ends`]. This approach covers lines, arcs of up to about a half turn,
    /// S-curves, and any other shape that is single-valued along its own principal axis, without
    /// requiring input from the caller beyond the points.
    ///
    /// The assumption is checked before fitting. When the curve doubles back along the axis (a
    /// hairpin, a loop, an arc well past a half turn) the ordering by projection interleaves the
    /// two branches, which shows up as perpendicular jumps between consecutive points that are
    /// comparable to the perpendicular extent of the whole set. The fit is refused with an error
    /// when the upper quartile of that jump ratio exceeds [`MAX_JUMP_RATIO`]. The check is
    /// scale-free and needs no noise estimate. On synthetic sweeps, it rejects hairpins, hooks,
    /// loops, and arcs of about 250 degrees or more while accepting every single-valued shape
    /// tested. Two limitations are worth knowing:
    ///
    /// * It is a check for gross violations. A curve that only slightly overhangs the ends of its
    ///   axis, such as an arc between about 190 and 240 degrees, passes and is fitted between its
    ///   extreme points rather than its true ends.
    /// * Noise-dominated straight data with only a couple of dozen points can be falsely rejected,
    ///   since the jumps between neighbouring noisy samples are then a large fraction of the
    ///   (purely noise) perpendicular extent. Fifty or more points are enough in practice.
    ///
    /// Use [`Self::from_fit_with_ends`] when the endpoints are known.
    ///
    /// # Arguments
    ///
    /// * `points`: The points to fit the curve to, in any order. At least two with positive
    ///   weight are required.
    /// * `weights`: If `Some`, a slice the same length as `points` of non-negative weights that
    ///   scale each point's residual. Points with zero weight take no part in the fit, the
    ///   principal axis, or the choice of endpoints. `None` weights all points equally.
    ///
    /// returns: A `Result` containing the fitted curve.
    ///
    /// # Failure
    ///
    /// Returns an error if there are too few points, the weight count does not match the point
    /// count, the points double back along their principal axis, or the minimization fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::common::cubic_spline::CubicSpline;
    ///
    /// // Points along an S-curve, in scrambled order and with no endpoints given.
    /// let truth = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.5, 0.0),
    ///     Point2::new(1.5, 2.0),
    ///     Point2::new(3.0, 2.0),
    /// );
    /// let points: Vec<Point2> = (0..=40)
    ///     .map(|i| truth.position(((i * 17) % 41) as f64 / 40.0))
    ///     .collect();
    ///
    /// let fitted = CubicSpline::from_fit_principal_axis(&points, None).unwrap();
    /// let queries = fitted.into_query();
    /// for p in &points {
    ///     assert!(queries.project_point(p).value < 1e-6);
    /// }
    /// ```
    pub fn from_fit_principal_axis(
        points: &[impl PCoords<D>],
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        check_fit_inputs(points.len(), weights)?;
        let points: Vec<Point<f64, D>> = points.iter().map(|p| Point::from(p.coords())).collect();
        let ones = vec![1.0; points.len()];
        let w = weights.unwrap_or(&ones);

        let basis = SvdBasis::from_points(&points, weights)
            .ok_or("could not compute a principal axis for the points")?;
        let local: Vec<Point<f64, D>> = points.iter().map(|p| basis.point_to_basis(p)).collect();
        let mut order: Vec<usize> = (0..points.len()).filter(|&i| w[i] > 0.0).collect();
        if order.len() < 2 {
            return Err("at least two points with positive weight are required".into());
        }
        order.sort_by(|&a, &b| local[a][0].total_cmp(&local[b][0]));

        if let Some(ratio) = jump_ratio(&local, &order)
            && ratio > MAX_JUMP_RATIO
        {
            return Err(format!(
                "the points do not run monotonically along their principal axis (jump ratio {:.2} \
                 exceeds {}); supply the endpoints with from_fit_with_ends instead",
                ratio, MAX_JUMP_RATIO
            )
            .into());
        }

        let first = points[order[0]];
        let last = points[*order.last().unwrap()];
        let initial = seed_with_ends(&points, w, &first, &last);

        // Parameters: p1 (D), p2 (D), then the perpendicular offsets of each end (D - 1 each).
        let perp: Vec<SVector<f64, D>> = basis.basis[1..].to_vec();
        let builder: SplineBuildFn<D> = Box::new(move |x: &DVector<f64>| {
            let p1 = Point::from(SVector::<f64, D>::from_iterator(x.iter().take(D).copied()));
            let p2 = Point::from(SVector::<f64, D>::from_iterator(
                x.iter().skip(D).take(D).copied(),
            ));
            let mut p0 = first;
            let mut p3 = last;
            for (k, v) in perp.iter().enumerate() {
                p0 += v * x[2 * D + k];
                p3 += v * x[3 * D - 1 + k];
            }
            Ok(CubicSpline::new(p0, p1, p2, p3))
        });
        let mut x0 = DVector::zeros(4 * D - 2);
        for k in 0..D {
            x0[k] = initial.p1[k];
            x0[D + k] = initial.p2[k];
        }
        Ok(fit_spline_to_points_weighted(&points, weights, &builder, x0)?.spline)
    }
}

/// The largest upper-quartile jump ratio (see [`jump_ratio`]) for which the points are taken to
/// run monotonically along their principal axis in [`CubicSpline::from_fit_principal_axis`]. On
/// synthetic sweeps, hairpins, hooks, loops, and arcs of 270 degrees or more score 0.66 and above,
/// single-valued shapes score 0.37 and below, and noise-dominated straight lines of fifty or more
/// points score 0.38 and below.
pub const MAX_JUMP_RATIO: f64 = 0.5;

/// The upper quartile of the perpendicular jumps between consecutive points in `order`, divided
/// by the perpendicular extent of the whole set. `local` holds the points in a frame whose first
/// coordinate is the principal axis. A single-valued curve produces small jumps (noise plus slope
/// times spacing), while interleaved branches produce jumps comparable to the extent. The upper
/// quartile is used instead of the median because the branches of a smooth hairpin converge toward
/// its turn, making many of its jumps small. Returns `None` when the perpendicular extent is
/// negligible, as it is for a straight line.
fn jump_ratio<const D: usize>(local: &[Point<f64, D>], order: &[usize]) -> Option<f64> {
    let mut lo = [f64::INFINITY; D];
    let mut hi = [f64::NEG_INFINITY; D];
    for &i in order {
        for k in 0..D {
            lo[k] = lo[k].min(local[i][k]);
            hi[k] = hi[k].max(local[i][k]);
        }
    }
    let extent: f64 = (1..D).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    if extent <= 1e-12 * (hi[0] - lo[0]).max(f64::MIN_POSITIVE) {
        return None;
    }
    let mut jumps: Vec<f64> = order
        .windows(2)
        .map(|pair| {
            (1..D)
                .map(|k| (local[pair[1]][k] - local[pair[0]][k]).powi(2))
                .sum::<f64>()
                .sqrt()
                / extent
        })
        .collect();
    jumps.sort_by(f64::total_cmp);
    Some(jumps[3 * jumps.len() / 4])
}

/// The implementation behind [`CubicSpline::from_fit_with_ends`] for points already converted to
/// the spline's dimension.
fn fit_with_ends<const D: usize>(
    points: &[Point<f64, D>],
    p0: &Point<f64, D>,
    p3: &Point<f64, D>,
    weights: Option<&[f64]>,
) -> Result<CubicSpline<D>> {
    let ones = vec![1.0; points.len()];
    let w = weights.unwrap_or(&ones);
    let initial = seed_with_ends(points, w, p0, p3);

    let (p0, p3) = (*p0, *p3);
    let builder: SplineBuildFn<D> = Box::new(move |x: &DVector<f64>| {
        let p1 = Point::from(SVector::<f64, D>::from_iterator(x.iter().take(D).copied()));
        let p2 = Point::from(SVector::<f64, D>::from_iterator(
            x.iter().skip(D).take(D).copied(),
        ));
        Ok(CubicSpline::new(p0, p1, p2, p3))
    });
    let x0 = DVector::from_iterator(
        2 * D,
        initial
            .p1
            .coords
            .iter()
            .chain(initial.p2.coords.iter())
            .copied(),
    );

    Ok(fit_spline_to_points_weighted(points, weights, &builder, x0)?.spline)
}

/// The best closed-form seed for a fit with fixed endpoints. Candidates include the straight chord
/// and the linear least-squares interior points under a chord-length parameterization for each
/// candidate ordering.
fn seed_with_ends<const D: usize>(
    points: &[Point<f64, D>],
    w: &[f64],
    p0: &Point<f64, D>,
    p3: &Point<f64, D>,
) -> CubicSpline<D> {
    let chord = p3 - p0;
    // List candidate seeds from cheapest to most expensive so ties favor the simplest one.
    let mut seeds = vec![CubicSpline::new(
        *p0,
        p0 + chord / 3.0,
        p3 - chord / 3.0,
        *p3,
    )];
    for order in candidate_orders(points, p0, p3) {
        if let Some(params) = chord_length_params(points, &order, w, p0, p3) {
            seeds.extend(seed_interior_points(points, &params, w, p0, p3));
        }
    }
    best_seed(seeds, points, w)
}

impl<const D: usize> CubicSpline<D> {
    /// Fit a cubic Bézier curve to a set of points, holding the two endpoints and the tangent
    /// directions at each end fixed and solving only for the lengths of the two tangent arms (the
    /// distances from `p0` to `p1` and from `p3` to `p2`).
    ///
    /// This is the most constrained and best-conditioned named fit. It has two scalar unknowns,
    /// both of which must be positive, and its geometry is a Hermite-style blend between two known
    /// end conditions. It is a natural choice for a curve that must leave one known feature (a
    /// line, an arc, or another curve) tangentially and arrive at another feature in the same way.
    /// It is also the single-segment core of Schneider's curve-fitting algorithm (Graphics Gems,
    /// 1990).
    ///
    /// Both tangents point in the direction of travel along the curve, from `p0` toward `p3`:
    /// `tangent0` is the direction the curve leaves `p0` (so `p1 = p0 + a0 * tangent0`), and
    /// `tangent3` is the direction the curve is traveling as it arrives at `p3` (so
    /// `p2 = p3 - a1 * tangent3`). The arm lengths `a0` and `a1` are held positive during the fit;
    /// a curve that would need to leave an endpoint against its given tangent cannot be produced.
    ///
    /// The fit starts from the best of several closed-form seeds and then refines against the true
    /// closest-point distances with a weighted Levenberg-Marquardt minimization
    /// ([`fit_spline_to_points_weighted`]). The seeds are Schneider's linear least-squares
    /// solution for the arm lengths under a chord-length parameterization of the points. It is
    /// computed once with the points ordered by projection onto the chord and once in the input
    /// order. A third seed uses one-third of the chord length for each arm as a fallback. See
    /// [`Self::from_fit_with_ends`] for why both point orderings are tried. The seed with the lowest
    /// weighted RMS closest-point distance is refined.
    ///
    /// Unlike [`Self::from_fit_with_ends`], the endpoints may coincide: a closed loop leaving and
    /// returning to the same point with different tangents is a valid curve.
    ///
    /// The fit determines the *curve* accurately, but its control points are only as well determined
    /// as the shape allows. To first order, shifting both interior control points along the curve
    /// reparameterizes it without changing its shape. The objective is therefore very flat in that
    /// direction, and control points fitted to exact data can differ from the generating control
    /// points by a fraction of a percent even when the curves agree much more closely. Compare the
    /// curves, rather than their control points, when checking a fit.
    ///
    /// # Arguments
    ///
    /// * `points`: The points to fit the curve to, in any order. At least two are required.
    /// * `p0`: The start point of the curve, held fixed.
    /// * `tangent0`: The unit direction the curve leaves `p0` in, held fixed.
    /// * `p3`: The end point of the curve, held fixed.
    /// * `tangent3`: The unit direction the curve is traveling in as it arrives at `p3`, held
    ///   fixed.
    /// * `weights`: If `Some`, a slice the same length as `points` of non-negative weights that
    ///   scale each point's residual. Points with zero weight have no effect on the fit. `None`
    ///   weights all points equally.
    ///
    /// returns: A `Result` containing the fitted curve.
    ///
    /// # Failure
    ///
    /// Returns an error if there are fewer than two points, the weight count does not match the
    /// point count, or the minimization fails.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, Vector2, UnitVec2};
    /// use engeom::common::cubic_spline::CubicSpline;
    ///
    /// let truth = CubicSpline::new(
    ///     Point2::new(0.0, 0.0),
    ///     Point2::new(1.0, 2.0),
    ///     Point2::new(2.0, 2.0),
    ///     Point2::new(3.0, 0.0),
    /// );
    /// let points: Vec<Point2> = (0..=20).map(|i| truth.position(i as f64 / 20.0)).collect();
    ///
    /// // Both tangents point in the direction of travel: out of p0 and into p3.
    /// let t0 = UnitVec2::new_normalize(Vector2::new(1.0, 2.0));
    /// let t3 = UnitVec2::new_normalize(Vector2::new(1.0, -2.0));
    /// let fitted =
    ///     CubicSpline::from_fit_hermite(&points, &truth.p0, &t0, &truth.p3, &t3, None).unwrap();
    ///
    /// // The fitted curve passes through the sampled points.
    /// let queries = fitted.into_query();
    /// for p in &points {
    ///     assert!(queries.project_point(p).value < 1e-6);
    /// }
    /// ```
    pub fn from_fit_hermite(
        points: &[impl PCoords<D>],
        p0: &Point<f64, D>,
        tangent0: &Unit<SVector<f64, D>>,
        p3: &Point<f64, D>,
        tangent3: &Unit<SVector<f64, D>>,
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        check_fit_inputs(points.len(), weights)?;

        let points: Vec<Point<f64, D>> = points.iter().map(|p| Point::from(p.coords())).collect();
        let ones = vec![1.0; points.len()];
        let w = weights.unwrap_or(&ones);

        let t0 = tangent0.into_inner();
        let t3 = tangent3.into_inner();
        let (q0, q3) = (*p0, *p3);
        let hermite = move |a0: f64, a1: f64| CubicSpline::new(q0, q0 + t0 * a0, q3 - t3 * a1, q3);

        // Use one-third of the chord as the fallback arm length. When the endpoints coincide, use
        // one-third of the farthest point's distance from p0 instead.
        let reach = points
            .iter()
            .map(|p| dist(p, p0))
            .fold(dist(p0, p3), f64::max);
        let mut seeds = vec![hermite(reach / 3.0, reach / 3.0)];
        for order in candidate_orders(&points, p0, p3) {
            if let Some(params) = chord_length_params(&points, &order, w, p0, p3)
                && let Some((a0, a1)) = seed_arm_lengths(&points, &params, w, p0, &t0, p3, &t3)
            {
                seeds.push(hermite(a0, a1));
            }
        }
        let initial = best_seed(seeds, &points, w);

        let builder: SplineBuildFn<D> = Box::new(move |x: &DVector<f64>| {
            if x[0] <= 0.0 || x[1] <= 0.0 {
                return Err("tangent arm lengths must be positive".into());
            }
            Ok(hermite(x[0], x[1]))
        });
        let x0 = DVector::from(vec![
            dist(&initial.p1, &initial.p0),
            dist(&initial.p2, &initial.p3),
        ]);

        Ok(fit_spline_to_points_weighted(&points, weights, &builder, x0)?.spline)
    }
}

/// Input validation shared by the named fitting constructors.
fn check_fit_inputs(n_points: usize, weights: Option<&[f64]>) -> Result<()> {
    if n_points < 2 {
        return Err("at least two points are required to fit a cubic spline".into());
    }
    if let Some(w) = weights
        && w.len() != n_points
    {
        return Err(format!(
            "expected {} weights for {} points, got {}",
            n_points,
            n_points,
            w.len()
        )
        .into());
    }
    Ok(())
}

/// The two candidate orderings of the points used to seed a fit: sorted by scalar projection onto
/// the chord from `p0` to `p3`, and in their given order. If the chord has zero length, projection
/// ordering degenerates to the given order.
fn candidate_orders<const D: usize>(
    points: &[Point<f64, D>],
    p0: &Point<f64, D>,
    p3: &Point<f64, D>,
) -> [Vec<usize>; 2] {
    let chord = p3 - p0;
    let mut by_projection: Vec<usize> = (0..points.len()).collect();
    by_projection.sort_by(|&a, &b| {
        let pa = (points[a] - p0).dot(&chord);
        let pb = (points[b] - p0).dot(&chord);
        pa.total_cmp(&pb)
    });
    let as_given: Vec<usize> = (0..points.len()).collect();
    [by_projection, as_given]
}

/// Chord-length parameters for the points with positive weight, taken in the given `order`. Each
/// point's parameter is its cumulative distance along the polyline from `p0`, through the ordered
/// points, to `p3`, expressed as a fraction of the total length. Returns `(index, t)` pairs, or
/// `None` if the polyline has zero length. Points with zero weight are omitted from both the
/// polyline and the result, so they have no effect on any seed built from these parameters.
fn chord_length_params<const D: usize>(
    points: &[Point<f64, D>],
    order: &[usize],
    weights: &[f64],
    p0: &Point<f64, D>,
    p3: &Point<f64, D>,
) -> Option<Vec<(usize, f64)>> {
    let mut cumulative = Vec::with_capacity(order.len());
    let mut prev = *p0;
    let mut total = 0.0;
    for &i in order.iter().filter(|&&i| weights[i] > 0.0) {
        total += dist(&prev, &points[i]);
        cumulative.push((i, total));
        prev = points[i];
    }
    total += dist(&prev, p3);
    if total <= f64::EPSILON {
        return None;
    }
    Some(
        cumulative
            .into_iter()
            .map(|(i, s)| (i, s / total))
            .collect(),
    )
}

/// Picks the seed with the lowest weighted RMS closest-point distance. Ties favor the earliest
/// seed.
fn best_seed<const D: usize>(
    seeds: Vec<CubicSpline<D>>,
    points: &[Point<f64, D>],
    weights: &[f64],
) -> CubicSpline<D> {
    seeds
        .into_iter()
        .map(|s| (weighted_rms(&s, points, weights), s))
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, s)| s)
        .expect("at least one seed is always supplied")
}

/// The weighted root-mean-square closest-point distance from `points` to `spline`.
fn weighted_rms<const D: usize>(
    spline: &CubicSpline<D>,
    points: &[Point<f64, D>],
    weights: &[f64],
) -> f64 {
    let queries = spline.clone().into_query();
    let (sum, wsum) = points
        .iter()
        .zip(weights)
        .fold((0.0, 0.0), |(sum, wsum), (p, w)| {
            let d = queries.project_point(p).value;
            (sum + w * d * d, wsum + w)
        });
    if wsum > 0.0 { (sum / wsum).sqrt() } else { 0.0 }
}

/// The linear least-squares seed for the interior control points of a curve with fixed endpoints,
/// given chord-length parameters for the points (see [`chord_length_params`]). With the parameters
/// fixed, the interior control points are the solution to a weighted 2 × 2 normal system in the
/// Bernstein basis, shared across all coordinates.
///
/// Returns `None` if the system is singular (fewer than two distinct parameters carrying weight).
fn seed_interior_points<const D: usize>(
    points: &[Point<f64, D>],
    params: &[(usize, f64)],
    weights: &[f64],
    p0: &Point<f64, D>,
    p3: &Point<f64, D>,
) -> Option<CubicSpline<D>> {
    let (mut a11, mut a12, mut a22) = (0.0, 0.0, 0.0);
    let mut r1 = SVector::<f64, D>::zeros();
    let mut r2 = SVector::<f64, D>::zeros();
    for &(i, t) in params {
        let [b0, b1, b2, b3] = bernstein_basis(t);
        let w = weights[i];
        let rhs = points[i].coords - p0.coords * b0 - p3.coords * b3;
        a11 += w * b1 * b1;
        a12 += w * b1 * b2;
        a22 += w * b2 * b2;
        r1 += rhs * (w * b1);
        r2 += rhs * (w * b2);
    }

    let det = a11 * a22 - a12 * a12;
    if det <= 1e-12 * a11 * a22 {
        return None;
    }
    let p1 = (r1 * a22 - r2 * a12) / det;
    let p2 = (r2 * a11 - r1 * a12) / det;
    Some(CubicSpline::new(*p0, Point::from(p1), Point::from(p2), *p3))
}

/// Schneider's linear least-squares seed for the two tangent arm lengths of a curve with fixed
/// endpoints and end tangents, given chord-length parameters for the points. With the parameters
/// fixed, the curve is linear in the arm lengths: `B(t) = (b0 + b1) p0 + (b2 + b3) p3 + a0 b1 t0 -
/// a1 b2 t3`, so the arm lengths solve a weighted 2 × 2 normal system.
///
/// Returns `None` if the system is singular or either arm length is non-positive, since
/// the fit requires the curve to leave each endpoint along its given tangent.
fn seed_arm_lengths<const D: usize>(
    points: &[Point<f64, D>],
    params: &[(usize, f64)],
    weights: &[f64],
    p0: &Point<f64, D>,
    t0: &SVector<f64, D>,
    p3: &Point<f64, D>,
    t3: &SVector<f64, D>,
) -> Option<(f64, f64)> {
    let (mut a11, mut a12, mut a22) = (0.0, 0.0, 0.0);
    let (mut r1, mut r2) = (0.0, 0.0);
    for &(i, t) in params {
        let [b0, b1, b2, b3] = bernstein_basis(t);
        let w = weights[i];
        let va = t0 * b1;
        let vb = -t3 * b2;
        let rhs = points[i].coords - p0.coords * (b0 + b1) - p3.coords * (b2 + b3);
        a11 += w * va.dot(&va);
        a12 += w * va.dot(&vb);
        a22 += w * vb.dot(&vb);
        r1 += w * va.dot(&rhs);
        r2 += w * vb.dot(&rhs);
    }

    let det = a11 * a22 - a12 * a12;
    if det <= 1e-12 * a11 * a22 {
        return None;
    }
    let a0 = (r1 * a22 - r2 * a12) / det;
    let a1 = (r2 * a11 - r1 * a12) / det;
    (a0 > 0.0 && a1 > 0.0).then_some((a0, a1))
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

    // --------------------------------------------------------------------------------------------
    // from_fit_with_ends
    // --------------------------------------------------------------------------------------------

    /// The maximum distance from samples along `truth` to the `fitted` curve.
    fn max_curve_deviation<const D: usize>(truth: &CubicSpline<D>, fitted: &CubicSpline<D>) -> f64 {
        let queries = fitted.clone().into_query();
        (0..=200)
            .map(|i| {
                queries
                    .project_point(&truth.position(i as f64 / 200.0))
                    .value
            })
            .fold(0.0, f64::max)
    }

    /// Asserts that `fitted` reproduces the shape of `truth` to within 1e-6. Control points are
    /// deliberately not compared: shifting both interior control points along the curve is a
    /// near-reparameterization that barely changes the shape, so they are only weakly determined.
    fn assert_same_curve<const D: usize>(truth: &CubicSpline<D>, fitted: &CubicSpline<D>) {
        let dev = max_curve_deviation(truth, fitted);
        assert!(dev < 1e-6, "curve deviation {} too large", dev);
    }

    fn shuffled<T: Clone>(items: &[T], rg: &mut RandomGeometry2) -> Vec<T> {
        let mut out = items.to_vec();
        for i in (1..out.len()).rev() {
            let j = rg.f64(0.0, (i + 1) as f64).floor() as usize;
            out.swap(i, j.min(i));
        }
        out
    }

    #[test]
    fn with_ends_recovers_exact_curve_2d() {
        let curve = truth();
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let fitted = CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn with_ends_recovers_exact_curve_3d() {
        let curve = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 1.0),
            Point3::new(2.0, 2.0, -1.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let points: Vec<Point3> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let fitted = CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn with_ends_noisy_unordered_points() {
        let mut rg = RandomGeometry2::from_seed(3);
        let curve = truth();
        let sigma = 0.01;
        for _ in 0..10 {
            let points: Vec<Point2> = (0..=40)
                .map(|i| curve.position(i as f64 / 40.0) + rg.gaussian_vector::<2>(sigma))
                .collect();
            let points = shuffled(&points, &mut rg);
            let fitted =
                CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, None).unwrap();
            let dev = max_curve_deviation(&curve, &fitted);
            assert!(dev < 2.0 * sigma, "deviation {} too large", dev);
        }
    }

    #[test]
    fn with_ends_hook_curve_ordered_and_shuffled() {
        // A hook that doubles back along its own chord, so ordering by chord projection is wrong.
        let curve = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(4.0, -1.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 1.0),
        );
        let points: Vec<Point2> = (0..=40).map(|i| curve.position(i as f64 / 40.0)).collect();
        let fitted = CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, None).unwrap();
        assert_same_curve(&curve, &fitted);

        let mut rg = RandomGeometry2::from_seed(5);
        let points = shuffled(&points, &mut rg);
        let fitted = CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn with_ends_zero_weight_outlier_ignored() {
        let curve = truth();
        let mut points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let mut weights = vec![1.0; points.len()];
        points.push(Point2::new(1.5, 5.0));
        weights.push(0.0);
        let fitted =
            CubicSpline::from_fit_with_ends(&points, &curve.p0, &curve.p3, Some(&weights)).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn with_ends_input_errors() {
        let p0 = Point2::new(0.0, 0.0);
        let p3 = Point2::new(3.0, 0.0);
        let one = vec![Point2::new(1.0, 1.0)];
        assert!(CubicSpline::from_fit_with_ends(&one, &p0, &p3, None).is_err());

        let two = vec![Point2::new(1.0, 1.0), Point2::new(2.0, 1.0)];
        assert!(CubicSpline::from_fit_with_ends(&two, &p0, &p0, None).is_err());
        assert!(CubicSpline::from_fit_with_ends(&two, &p0, &p3, Some(&[1.0])).is_err());
    }

    // --------------------------------------------------------------------------------------------
    // from_fit_hermite
    // --------------------------------------------------------------------------------------------

    /// The end tangents of a curve in the direction of travel, as `from_fit_hermite` expects them.
    fn end_tangents<const D: usize>(
        c: &CubicSpline<D>,
    ) -> (Unit<SVector<f64, D>>, Unit<SVector<f64, D>>) {
        (
            Unit::new_normalize(c.p1 - c.p0),
            Unit::new_normalize(c.p3 - c.p2),
        )
    }

    #[test]
    fn hermite_recovers_exact_curve_2d() {
        let curve = truth();
        let (t0, t3) = end_tangents(&curve);
        let points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn hermite_recovers_exact_curve_3d() {
        let curve = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 1.0),
            Point3::new(2.0, 2.0, -1.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let (t0, t3) = end_tangents(&curve);
        let points: Vec<Point3> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn hermite_noisy_unordered_points() {
        let mut rg = RandomGeometry2::from_seed(13);
        let curve = truth();
        let (t0, t3) = end_tangents(&curve);
        let sigma = 0.01;
        for _ in 0..10 {
            let points: Vec<Point2> = (0..=40)
                .map(|i| curve.position(i as f64 / 40.0) + rg.gaussian_vector::<2>(sigma))
                .collect();
            let points = shuffled(&points, &mut rg);
            let fitted =
                CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None)
                    .unwrap();
            let dev = max_curve_deviation(&curve, &fitted);
            assert!(dev < 2.0 * sigma, "deviation {} too large", dev);
        }
    }

    #[test]
    fn hermite_hook_curve_ordered_and_shuffled() {
        let curve = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(4.0, -1.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 1.0),
        );
        let (t0, t3) = end_tangents(&curve);
        let points: Vec<Point2> = (0..=40).map(|i| curve.position(i as f64 / 40.0)).collect();
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None).unwrap();
        assert_same_curve(&curve, &fitted);

        let mut rg = RandomGeometry2::from_seed(5);
        let points = shuffled(&points, &mut rg);
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn hermite_closed_loop_with_coincident_endpoints() {
        // A teardrop that leaves and returns to the same point with different tangents.
        let curve = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 2.0),
            Point2::new(3.0, -2.0),
            Point2::new(0.0, 0.0),
        );
        let (t0, t3) = end_tangents(&curve);
        let points: Vec<Point2> = (0..=40).map(|i| curve.position(i as f64 / 40.0)).collect();
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, None).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn hermite_collinear_degenerate_case() {
        // Points on a straight line with tangents along it: every arm length fits perfectly, so
        // the fit must succeed and stay on the line.
        let p0 = Point2::new(0.0, 0.0);
        let p3 = Point2::new(3.0, 0.0);
        let t = Unit::new_normalize(Vector2::new(1.0, 0.0));
        let points: Vec<Point2> = (0..=20)
            .map(|i| Point2::new(3.0 * i as f64 / 20.0, 0.0))
            .collect();
        let fitted = CubicSpline::from_fit_hermite(&points, &p0, &t, &p3, &t, None).unwrap();
        assert!(fitted.p1.y.abs() < 1e-12 && fitted.p2.y.abs() < 1e-12);
        assert!(fitted.p1.x > 0.0 && fitted.p2.x < 3.0);
        let queries = fitted.into_query();
        for p in &points {
            assert!(queries.project_point(p).value < 1e-9);
        }
    }

    #[test]
    fn hermite_zero_weight_outlier_ignored() {
        let curve = truth();
        let (t0, t3) = end_tangents(&curve);
        let mut points: Vec<Point2> = (0..=20).map(|i| curve.position(i as f64 / 20.0)).collect();
        let mut weights = vec![1.0; points.len()];
        points.push(Point2::new(1.5, 5.0));
        weights.push(0.0);
        let fitted =
            CubicSpline::from_fit_hermite(&points, &curve.p0, &t0, &curve.p3, &t3, Some(&weights))
                .unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn hermite_input_errors() {
        let p0 = Point2::new(0.0, 0.0);
        let p3 = Point2::new(3.0, 0.0);
        let t = Unit::new_normalize(Vector2::new(1.0, 1.0));
        let one = vec![Point2::new(1.0, 1.0)];
        assert!(CubicSpline::from_fit_hermite(&one, &p0, &t, &p3, &t, None).is_err());
        let two = vec![Point2::new(1.0, 1.0), Point2::new(2.0, 1.0)];
        assert!(CubicSpline::from_fit_hermite(&two, &p0, &t, &p3, &t, Some(&[1.0])).is_err());
    }

    // --------------------------------------------------------------------------------------------
    // from_fit_principal_axis
    // --------------------------------------------------------------------------------------------

    fn samples<const D: usize>(c: &CubicSpline<D>, n: usize) -> Vec<Point<f64, D>> {
        (0..n)
            .map(|i| c.position(i as f64 / (n - 1) as f64))
            .collect()
    }

    fn gentle_s() -> CubicSpline<2> {
        CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.5, 0.0),
            Point2::new(1.5, 2.0),
            Point2::new(3.0, 2.0),
        )
    }

    #[test]
    fn principal_axis_exact_shapes() {
        let mut rg = RandomGeometry2::from_seed(21);
        let shapes = [
            truth(),
            gentle_s(),
            // Steep S-curve.
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(2.0, 0.0),
                Point2::new(1.0, 3.0),
                Point2::new(3.0, 3.0),
            ),
            // Wide U, single-valued along its principal axis.
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(4.0, 0.0),
                Point2::new(4.0, 3.0),
                Point2::new(0.0, 3.0),
            ),
            // Straight line.
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.5),
                Point2::new(2.0, 1.0),
                Point2::new(3.0, 1.5),
            ),
        ];
        for curve in shapes {
            let points = shuffled(&samples(&curve, 41), &mut rg);
            let fitted = CubicSpline::from_fit_principal_axis(&points, None).unwrap();
            assert_same_curve(&curve, &fitted);
        }
    }

    #[test]
    fn principal_axis_noisy_shapes() {
        let mut rg = RandomGeometry2::from_seed(22);
        let sigma = 0.01;
        for curve in [truth(), gentle_s()] {
            for _ in 0..5 {
                let points: Vec<Point2> = samples(&curve, 60)
                    .iter()
                    .map(|p| p + rg.gaussian_vector::<2>(sigma))
                    .collect();
                let points = shuffled(&points, &mut rg);
                let fitted = CubicSpline::from_fit_principal_axis(&points, None).unwrap();
                let dev = max_curve_deviation(&curve, &fitted);
                assert!(dev < 4.0 * sigma, "deviation {} too large", dev);
            }
        }
    }

    #[test]
    fn principal_axis_noisy_line_accepted() {
        let mut rg = RandomGeometry2::from_seed(23);
        let points: Vec<Point2> = (0..100)
            .map(|i| Point2::new(3.0 * i as f64 / 99.0, 0.0) + rg.gaussian_vector::<2>(0.3))
            .collect();
        assert!(CubicSpline::from_fit_principal_axis(&points, None).is_ok());
    }

    #[test]
    fn principal_axis_rejects_doubling_back() {
        let mut rg = RandomGeometry2::from_seed(24);
        let shapes = [
            // hairpin
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(4.0, 0.0),
                Point2::new(4.0, 1.0),
                Point2::new(0.0, 1.0),
            ),
            // hook
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(4.0, -1.0),
                Point2::new(4.0, 4.0),
                Point2::new(0.0, 1.0),
            ),
            // loop
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(3.0, 2.0),
                Point2::new(3.0, -2.0),
                Point2::new(0.0, 0.0),
            ),
        ];
        for curve in shapes {
            let points = shuffled(&samples(&curve, 60), &mut rg);
            assert!(CubicSpline::from_fit_principal_axis(&points, None).is_err());
        }

        // A 300 degree arc
        let half = 150.0f64.to_radians();
        let arc: Vec<Point2> = (0..60)
            .map(|i| {
                let a = -half + 2.0 * half * i as f64 / 59.0;
                Point2::new(a.cos(), a.sin())
            })
            .collect();
        assert!(CubicSpline::from_fit_principal_axis(&arc, None).is_err());
    }

    #[test]
    fn principal_axis_zero_weight_points_excluded() {
        // A distant zero-weight point must not become an endpoint or tilt the axis.
        let curve = truth();
        let mut points = samples(&curve, 41);
        let mut weights = vec![1.0; points.len()];
        points.push(Point2::new(20.0, 20.0));
        weights.push(0.0);
        let fitted = CubicSpline::from_fit_principal_axis(&points, Some(&weights)).unwrap();
        assert_same_curve(&curve, &fitted);
    }

    #[test]
    fn principal_axis_3d() {
        let mut rg = RandomGeometry2::from_seed(25);
        let curve = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 1.0),
            Point3::new(2.0, 2.0, -1.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let points = shuffled(&samples(&curve, 41), &mut rg);
        let fitted = CubicSpline::from_fit_principal_axis(&points, None).unwrap();
        assert_same_curve(&curve, &fitted);

        let hairpin = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 1.0),
            Point3::new(4.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 0.0),
        );
        let points = shuffled(&samples(&hairpin, 60), &mut rg);
        assert!(CubicSpline::from_fit_principal_axis(&points, None).is_err());
    }

    #[test]
    fn principal_axis_input_errors() {
        let one = vec![Point2::new(1.0, 1.0)];
        assert!(CubicSpline::from_fit_principal_axis(&one, None).is_err());
        let two = vec![Point2::new(1.0, 1.0), Point2::new(2.0, 1.0)];
        assert!(CubicSpline::from_fit_principal_axis(&two, Some(&[1.0])).is_err());
        assert!(CubicSpline::from_fit_principal_axis(&two, Some(&[1.0, 0.0])).is_err());
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
