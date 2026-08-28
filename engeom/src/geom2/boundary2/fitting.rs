//! This module generalizes fitting of boundary geometry to sample points

use crate::common::consensus::weights::{MagsacWeight, estimate_sigma_max};
use crate::common::points::dist;
use crate::common::{PCoords, RefinementHalt, SPCoords, SolveQuality, TerminationReason};
use crate::geom2::{Boundary2, Manifold1Pos2};
use crate::na::{DVector, Dyn, Matrix, Owned, U1, Vector};
use crate::{Point2, Result, SurfacePoint2, VecDot};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

pub type BndBuildFn = Box<dyn Fn(&DVector<f64>) -> Result<Boundary2>>;
const DELTA: f64 = 1e-6;

/// The residual of a single sample: its distance to the closest position on the boundary, signed
/// by which side of the boundary's normal it fell on.
///
/// The sign carries no extra information for the fit itself. Because the objective squares the
/// residual, and because the jacobian picks up the same sign that the residual does, both `J^T J`
/// and `J^T r` are unchanged by it: a signed and an unsigned fit take identical Levenberg-Marquardt
/// steps and converge to the same parameters.
///
/// It is worth having for two other reasons. The jacobian here is a forward finite difference, and
/// the derivative of an absolute value is wrong for any sample sitting within `DELTA` of the
/// boundary, where the disturbed and undisturbed distances straddle zero. And a robust noise
/// estimate by median absolute deviation assumes residuals that are centered on zero when the fit
/// is good; unsigned distances are all positive, which would have it measure the spread about the
/// typical distance instead and report a scale that is too small.
fn signed_residual(m: &Manifold1Pos2, p: &impl PCoords<2>) -> f64 {
    dist(m, p) * m.scalar_projection(p).signum()
}

/// The residual degrees of freedom for the MAGSAC++ weight function. A boundary residual is a
/// full point-to-manifold distance in the plane, so it follows a chi distribution with two
/// degrees of freedom.
const RESIDUAL_DOF: usize = 2;

/// The default Levenberg-Marquardt evaluation budget, matching the `levenberg-marquardt` crate's
/// own default so that the options struct changes nothing unless a caller asks it to.
const DEFAULT_PATIENCE: usize = 100;

/// Options controlling a boundary fit.
///
/// The default is a plain unweighted least-squares fit with no robust refinement, which is what
/// this module did before the option existed. Robustness is opt-in because it is not free: each
/// refinement round is a whole extra Levenberg-Marquardt solve, and every step of every solve
/// rebuilds the boundary `n_params` times to take a finite-difference jacobian. Call
/// [`BoundaryFitOptions::robust`] when the data is expected to carry outliers and the cost is
/// worth paying.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryFitOptions {
    /// The number of iteratively reweighted refinement rounds to perform after the initial
    /// unweighted solve. Zero, the default, disables robust weighting entirely.
    pub refinement_steps: usize,

    /// The MAGSAC++ upper noise bound, estimated from the initial residuals via the median
    /// absolute deviation when `None`.
    ///
    /// This is in the units of the geometry: a residual is a signed distance from a sample to the
    /// boundary, and nothing normalizes it. Supplying it explicitly is worth doing when the
    /// measurement noise is known, since the estimate is taken from a fit that the outliers have
    /// already influenced.
    pub sigma_max: Option<f64>,

    /// The Levenberg-Marquardt evaluation budget, as a multiplier on the parameter count. Must be
    /// greater than zero.
    pub patience: usize,
}

impl Default for BoundaryFitOptions {
    fn default() -> Self {
        Self {
            refinement_steps: 0,
            sigma_max: None,
            patience: DEFAULT_PATIENCE,
        }
    }
}

impl BoundaryFitOptions {
    /// Options with robust refinement switched on, using the default number of rounds and a noise
    /// bound estimated from the data.
    pub fn robust() -> Self {
        Self {
            refinement_steps: 4,
            ..Default::default()
        }
    }

    fn validate(&self) -> Result<()> {
        if self.patience == 0 {
            return Err("patience must be greater than zero".into());
        }
        if let Some(s) = self.sigma_max
            && (!s.is_finite() || s <= 0.0)
        {
            return Err(
                format!("sigma_max is {s}, but must be finite and strictly positive").into(),
            );
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct BoundaryFitResult {
    pub params: DVector<f64>,

    /// The geometric residuals of the fit, signed distances in the units of the model. These are
    /// the plain residuals, not the weighted values the solver minimized, so they describe the
    /// geometry rather than the fitting machinery.
    pub residuals: DVector<f64>,

    solves: Vec<TerminationReason>,
    halt: Option<RefinementHalt>,
}

impl BoundaryFitResult {
    pub(crate) fn new(
        params: DVector<f64>,
        residuals: DVector<f64>,
        solves: Vec<TerminationReason>,
        halt: Option<RefinementHalt>,
    ) -> Self {
        Self {
            params,
            residuals,
            solves,
            halt,
        }
    }

    /// How every solve which contributed to this fit terminated, the initial one first.
    pub fn solves(&self) -> &[TerminationReason] {
        &self.solves
    }

    /// The number of robust refinement rounds which completed, which is zero for a plain fit.
    pub fn refinement_rounds(&self) -> usize {
        self.solves.len().saturating_sub(1)
    }

    /// The worst quality among the solves which contributed to this fit.
    pub fn quality(&self) -> SolveQuality {
        self.solves
            .iter()
            .map(SolveQuality::from_termination)
            .fold(SolveQuality::Converged, |acc, q| acc.worse_of(q))
    }

    /// Whether every solve which contributed to this fit converged.
    pub fn converged(&self) -> bool {
        self.quality() == SolveQuality::Converged
    }

    /// Why robust refinement stopped early, if it did.
    pub fn halt(&self) -> Option<&RefinementHalt> {
        self.halt.as_ref()
    }
}

// =============================================================================================
// Fitting a boundary to discrete points
// =============================================================================================

/// Given a set of points, a function which takes a parameter vector and produces a [`Boundary2`],
/// and an initial guess, this function will attempt to perform a Levenberg-Marquardt best fit of
/// the parameters to produce a `Boundary2` that minimizes the residuals of the points projected
/// onto the boundary.
///
/// You will need to provide a `DVector` with a reasonable initial guess, and a [`BndBuildFn`] that
/// accepts a `DVector` of that size and produces a boundary. The function is allowed to return a
/// `Err(Box<dyn Error>)` if the parameters produce invalid geometry, but the initial guess values
/// provided may not fail or this function will exit before attempting the minimization.
///
/// Because the residuals are the values of the points projected onto the nearest point on the
/// boundary, be aware that there is no inherent pressure preventing the boundary geometry from
/// growing beyond the sample points _unless_ it exists within the nature of the geometry. For
/// example, if the points lie on a line and the boundary geometry is a segment, there is nothing
/// that will constrain the segment from growing beyond the ends of the points.
///
/// For this reason, it is best to write the builder function so that any constraints are inherently
/// contained within it. For example, write the builder function so that the segment ends are
/// only parameterized in one direction.
///
/// Refer to the examples for more information.
///
/// # Arguments
///
/// * `points`: A slice of `Point2` that the boundary will be fit to. The residuals are calculated
///   as the distance from each point to the closest point on the boundary, so the points will
///   "pull" the closest area of the boundary towards themselves, but won't have any effect on
///   the boundary farther away.
/// * `builder`: A function that takes a `DVector` of parameters and produces a `Boundary2`.
/// * `initial`: An initial guess for the parameters. The `builder` function must return an
///   `Ok(Boundary2)` with this input for the fitting algorithm to initialize. The algorithm will
///   use the size of the initial guess to determine the number of parameters, so make sure it
///   matches what the `builder` function is expecting.
/// * `ignore_ends`: if `true`, points that project onto the ends of an open boundary will have
///   residuals of zero
///
/// returns: Result<BoundaryFitResult>
///
/// # Examples
///
/// ```
/// // This example will fit a boundary consisting of three line segments to a bunch of points
/// // arranged in a triangle.
/// // ------------------------------------------------------------------------------------------
/// use engeom::{Point2, DVector};
/// use engeom::common::{to_points, fill_gaps};
/// use engeom::geom2::{
///     BndBuildFn, BoundaryData2, BoundaryEditor, BoundaryFitOptions, fit_boundary_to_points,
/// };
/// use approx::assert_relative_eq;
///
/// // We'll create the initial corners and then use the `fill_gaps` helper to generate points
/// // between them. Normally this would be your observed data.
/// let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
/// let points = fill_gaps(&corners, 0.1);
///
/// // Here we define the function which creates the boundary from six parameters. This is an
/// // extremely simple parameterization that just encodes the three corners as x,y pairs
/// let builder: BndBuildFn = Box::new(|params: &DVector| {
///    let mut bdata = BoundaryData2::new_open(Point2::new(params[0], params[1]));
///    bdata.add_seg_xy(params[2], params[3]);
///    bdata.add_seg_xy(params[4], params[5]);
///    bdata.add_seg_xy(params[0], params[1]);
///    bdata.try_to_boundary()
/// });
///
/// // We'll provide an initial guess which has the corners in roughly the correct areas, just way
/// // too large for the actual data and not aligned with the edges. Then we'll run the fitting
/// // algorithm.
/// let initial = DVector::from(vec![0.0, 0.0, 4.0, 0.0, 1.0, 7.0]);
/// let result = fit_boundary_to_points(
///     &points, &builder, initial, false, &BoundaryFitOptions::default(),
/// )
/// .unwrap();
///
/// // Finally we'll verify that the corners match the ones we originally provided.
/// let expected = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
/// assert_relative_eq!(result.params, expected, epsilon = 1.0e-6);
///
/// ```
pub fn fit_boundary_to_points(
    points: &[Point2],
    builder: &BndBuildFn,
    initial: DVector<f64>,
    ignore_ends: bool,
    opts: &BoundaryFitOptions,
) -> Result<BoundaryFitResult> {
    let fitting = BoundaryToPoints::new(points, ignore_ends);
    solve(&fitting, builder, initial, opts)
}

struct BoundaryToPoints<'a> {
    points: &'a [Point2],
    ignore_ends: bool,
}

impl<'a> BoundaryToPoints<'a> {
    fn new(points: &'a [Point2], ignore_ends: bool) -> Self {
        BoundaryToPoints {
            points,
            ignore_ends,
        }
    }
}

impl BoundaryFittable for BoundaryToPoints<'_> {
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>) {
        let bounds = f64::EPSILON..(boundary.length() - f64::EPSILON);
        let mut res = DVector::zeros(self.points.len());
        let mut weights = DVector::zeros(self.points.len());
        weights.fill(1.0);

        for i in 0..self.points.len() {
            let (_, m) = boundary.at_closest_to_point(&self.points[i]);
            if self.ignore_ends {
                weights[i] = if bounds.contains(&m.l) { 1.0 } else { 0.0 };
            }
            res[i] = signed_residual(&m, &self.points[i]);
        }

        (res, weights)
    }

    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64 {
        let (_, m) = boundary.at_closest_to_point(&self.points[sample_i]);
        signed_residual(&m, &self.points[sample_i])
    }
}

// =============================================================================================
// Fitting a boundary to discrete surface points
// =============================================================================================

/// Given a set of surface points, a function which takes a parameter vector and produces a
/// [`Boundary2`], and an initial guess, this function will attempt to perform a Levenberg-Marquardt
/// best fit of the parameters to produce a `Boundary2` that minimizes the residuals of the surface
/// points projected onto the boundary and weighted by their dot product with the boundary normal
/// at the projection site.
///
/// You will need to provide a `DVector` with a reasonable initial guess, and a [`BndBuildFn`] that
/// accepts a `DVector` of that size and produces a boundary. The function is allowed to return a
/// `Err(Box<dyn Error>)` if the parameters produce invalid geometry, but the initial guess values
/// provided may not fail or this function will exit before attempting the minimization.
///
/// Refer to the examples for more information.
///
/// # Arguments
///
/// * `points`: a slice of `SurfacePoint2` entities, representing points and their surface normals
/// * `builder`: A function that takes a `DVector` of parameters and produces a `Boundary2`.
/// * `initial`: An initial guess for the parameters. The `builder` function must return an
///   `Ok(Boundary2)` with this input for the fitting algorithm to initialize. The algorithm will
///   use the size of the initial guess to determine the number of parameters, so make sure it
///   matches what the `builder` function is expecting.
/// * `weight_mode`: the mode by which the surface normal dot products produce a weight. If you
///   use `AsIs`, a point can have a negative weight if the normals face in opposite direction. If
///   you use `Abs`, it will de-weight surfaces pointing orthogonally to each other, but won't care
///   about the direction. If you use `ClampPos`, it will only consider points with normals facing
///   in the same direction. I recommend using `ClampPos` if you know the boundary and sample
///   normals should be facing the same direction, and `Abs` if you just want to de-weight
///   orthogonal normals but aren't sure if they're facing the same way.
/// * `ignore_ends`: if `true`, points that project onto the ends of an open boundary will have
///   residuals of zero
///
/// returns: Result<BoundaryFitResult>
///
/// # Examples
///
/// ```
/// // This example shows how a very simple boundary fitting can use the surface normals to reject
/// // samples facing the wrong direction
/// // ------------------------------------------------------------------------------------------
/// use engeom::{VecDot, Vector2, SurfacePoint2, DVector, Point2};
/// use engeom::common::{to_points, fill_gaps, linear_space};
/// use engeom::geom2::{
///     BndBuildFn, BoundaryData2, BoundaryEditor, BoundaryFitOptions,
///     fit_boundary_to_surface_points,
/// };
/// use approx::assert_relative_eq;
///
/// // We'll create ten points facing in +Y at Y=0
/// let good = linear_space(0.0, 10.0, 10)
///     .iter()
///     .map(|x| SurfacePoint2::new(Point2::new(*x, 0.0), Vector2::y_axis()))
///     .collect::<Vec<_>>();
///
/// // We'll create ten bad points facing in +x at Y=1
/// let bad = linear_space(0.0, 10.0, 10)
///     .iter()
///     .map(|x| SurfacePoint2::new(Point2::new(*x, 1.0), Vector2::x_axis()))
///     .collect::<Vec<_>>();
///
/// // We'll combine the good and the bad points together
/// let samples = [good, bad].concat();
///
/// // Our builder function will create a single line segment from X=0 to X=10. Because `engeom`'s
/// // boundaries follows the anti-clockwise winding order convention, the segment's surface normal
/// // will be facing in -Y. If we reversed the order it would face in +Y.
/// let builder: BndBuildFn = Box::new(|params: &DVector| {
///     let mut bdata = BoundaryData2::new_open([0.0, params[0]].into());
///     bdata.add_seg_xy(10.0, params[1]);
///     bdata.try_to_boundary()
/// });
///
/// // We'll create our initial guess right near the bad points.
/// let initial = DVector::from(vec![1.0, 1.0]);
///
/// // When we do the fitting, we'll use the `VecDot::Abs` mode, which will not care that the
/// // boundary normal is facing in the opposite direction of the good points, but will de-weight
/// // the bad points because they are orthogonal.
/// let result =
///     fit_boundary_to_surface_points(
///         &samples, &builder, initial, VecDot::Abs, false, &BoundaryFitOptions::default(),
///     )
///     .unwrap();
///
/// // As we can see, the fit ended at the position of the good points.
/// let expected = DVector::from(vec![0.0, 0.0]);
/// assert_relative_eq!(result.params, expected, epsilon = 1.0e-6);
/// ```
pub fn fit_boundary_to_surface_points(
    points: &[SurfacePoint2],
    builder: &BndBuildFn,
    initial: DVector<f64>,
    weight_mode: VecDot,
    ignore_ends: bool,
    opts: &BoundaryFitOptions,
) -> Result<BoundaryFitResult> {
    let fitting = BoundaryToSurfacePoints::new(points, weight_mode, ignore_ends);
    solve(&fitting, builder, initial, opts)
}

struct BoundaryToSurfacePoints<'a> {
    points: &'a [SurfacePoint2],
    weight_mode: VecDot,
    ignore_ends: bool,
}

impl<'a> BoundaryToSurfacePoints<'a> {
    fn new(points: &'a [SurfacePoint2], weight_mode: VecDot, ignore_ends: bool) -> Self {
        BoundaryToSurfacePoints {
            points,
            weight_mode,
            ignore_ends,
        }
    }
}

impl BoundaryFittable for BoundaryToSurfacePoints<'_> {
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>) {
        let bounds = f64::EPSILON..(boundary.length() - f64::EPSILON);
        let mut res = DVector::zeros(self.points.len());
        let mut weights = DVector::zeros(self.points.len());
        weights.fill(1.0);

        for i in 0..self.points.len() {
            let (_, m) = boundary.at_closest_to_point(&self.points[i]);
            let dot = self.points[i].normal.dot(&m.normal);
            weights[i] = match self.weight_mode {
                VecDot::AsIs => dot,
                VecDot::Abs => dot.abs(),
                VecDot::ClampPos => dot.max(0.0),
            };

            if self.ignore_ends && !bounds.contains(&m.l) {
                weights[i] = 0.0;
            }

            res[i] = signed_residual(&m, &self.points[i]);
        }

        (res, weights)
    }

    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64 {
        let (_, m) = boundary.at_closest_to_point(&self.points[sample_i]);
        signed_residual(&m, &self.points[sample_i])
    }
}

// =============================================================================================
// `BoundaryFit` and `BoundaryFittable` together offer a generic mechanism for performing the
// fitting of boundaries to a set of samples of a finite count.
// =============================================================================================

// =============================================================================================
// The solve
// =============================================================================================

/// Runs the fit: one unweighted Levenberg-Marquardt solve, optionally followed by rounds of
/// iteratively reweighted least squares using MAGSAC++ weights.
///
/// The structure mirrors `geom2::align2::points_to_surface2`, for the same reasons. An initial
/// unweighted solve comes first because the noise scale has to be estimated from residuals, and
/// there are none until something has been fitted. Within each refinement round the weights are
/// held fixed, so that the jacobian stays consistent with the residual it differentiates.
///
/// An `Err` is reserved for having no answer at all: rejected options, a builder which fails on
/// the initial parameters, or an initial solve which broke down. A refinement round which breaks
/// down is rolled back to the previous round's parameters and reported on the result.
fn solve(
    fitting: &dyn BoundaryFittable,
    builder: &BndBuildFn,
    initial: DVector<f64>,
    opts: &BoundaryFitOptions,
) -> Result<BoundaryFitResult> {
    opts.validate()?;

    let lm = LevenbergMarquardt::new().with_patience(opts.patience);
    let n_params = initial.len();

    let problem = BoundaryFit::try_new(fitting, builder, initial)?;
    let (mut problem, termination) = run(&lm, problem);
    if !termination.was_successful() {
        return Err(format!("Fitting failed: {termination:?}").into());
    }

    let mut solves = vec![termination];
    let mut halt = None;

    if opts.refinement_steps > 0 {
        match resolve_sigma_max(opts, &problem) {
            None => halt = Some(RefinementHalt::NoNoiseEstimate),
            Some(sigma_max) => {
                let weighting = MagsacWeight::new(sigma_max, RESIDUAL_DOF);

                for _ in 0..opts.refinement_steps {
                    let weighted = problem.count_if_reweighted(&weighting);
                    if weighted < n_params {
                        halt = Some(RefinementHalt::Underdetermined {
                            weighted,
                            params: n_params,
                        });
                        break;
                    }

                    let last_good = problem.params.clone();
                    problem.apply_magsac_weights(&weighting);

                    let (next, termination) = run(&lm, problem);
                    problem = next;

                    if !SolveQuality::from_termination(&termination).is_usable() {
                        problem.restore(&last_good);
                        halt = Some(RefinementHalt::SolveFailed(termination));
                        break;
                    }
                    solves.push(termination);
                }
            }
        }
    }

    let residuals = problem
        .residuals
        .clone()
        .ok_or("the fitted parameters do not produce a valid boundary")?;

    Ok(BoundaryFitResult::new(
        problem.params,
        residuals,
        solves,
        halt,
    ))
}

fn run<'a>(
    lm: &LevenbergMarquardt<f64>,
    problem: BoundaryFit<'a>,
) -> (BoundaryFit<'a>, TerminationReason) {
    let (result, report) = lm.minimize(problem);
    (result, report.termination)
}

fn resolve_sigma_max(opts: &BoundaryFitOptions, problem: &BoundaryFit<'_>) -> Option<f64> {
    match opts.sigma_max {
        Some(s) => Some(s),
        None => estimate_sigma_max(problem.residuals.as_ref()?.as_slice()),
    }
}

pub trait BoundaryFittable {
    /// Given a boundary, this should return two equally sized `DVector`s of the residuals and
    /// the weights. The weights are not calculated independently during the jabcobian, but simply
    /// are reused from this step. Return the _un-weighted_ residuals, and the weights that will
    /// be used to scale them. For any element `i` in the residual vector, the value of
    /// `residual_only(i, ...)` should be identical.
    ///
    /// Return in the order: `(residuals, weights)`
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>);

    /// Given a boundary, this should return the residual for sample `sample_i`. This step is used
    /// for the numerical jacobians. It should produce values identical to the ones from
    /// `residuals_and_weights`
    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64;
}

struct BoundaryFit<'a> {
    fitting: &'a dyn BoundaryFittable,
    params: DVector<f64>,
    builder: &'a BndBuildFn,
    current: Option<Boundary2>,
    residuals: Option<DVector<f64>>,

    /// The geometric weight of each sample, which depends on where the sample projects and so is
    /// recomputed whenever the parameters move.
    weights: Option<DVector<f64>>,

    /// The robust weight of each sample, held fixed for the whole of a solve and refreshed only
    /// between them. Reweighting inside a solve would leave the jacobian differentiating a
    /// different objective than the residual reports.
    magsac_weights: DVector<f64>,
}

impl<'a> BoundaryFit<'a> {
    fn try_new(
        fitting: &'a dyn BoundaryFittable,
        builder: &'a BndBuildFn,
        initial: DVector<f64>,
    ) -> Result<Self> {
        // Check to make sure that the initial value doesn't fail
        let _ = builder(&initial)?;

        let mut problem = BoundaryFit {
            fitting,
            params: DVector::zeros(initial.len()),
            builder,
            current: None,
            residuals: None,
            weights: None,
            magsac_weights: DVector::zeros(0),
        };

        problem.set_params(&initial);

        // Sized once the first evaluation has told us how many samples there are, and left at one
        // apiece so that the opening solve is unweighted.
        let n = problem.residuals.as_ref().map_or(0, |r| r.len());
        problem.magsac_weights = DVector::from_element(n, 1.0);

        Ok(problem)
    }

    /// The factor applied to both the residual and the jacobian row of a sample.
    ///
    /// The geometric weight enters linearly and the robust weight as a square root, which is not
    /// an oversight. A residual scaled by `w` contributes `w^2 r^2` to the objective, so the
    /// square root is what makes a MAGSAC weight mean what it says. The geometric weight predates
    /// that and is left alone: it is a caller-facing knob whose established behavior is a scale on
    /// the residual, and quietly changing it would move every existing fit.
    fn scale(&self, weights: &DVector<f64>, i: usize) -> f64 {
        weights[i] * self.magsac_weights[i].sqrt()
    }

    /// Recomputes the robust weights from the current residuals.
    fn apply_magsac_weights(&mut self, weighting: &MagsacWeight) {
        if let Some(residuals) = &self.residuals {
            for i in 0..residuals.len() {
                self.magsac_weights[i] = weighting.weight(residuals[i].abs());
            }
        }
    }

    /// How many samples would still carry weight after reweighting. A round which would leave
    /// fewer of them than there are parameters is rank-deficient and must not be run.
    fn count_if_reweighted(&self, weighting: &MagsacWeight) -> usize {
        let (Some(residuals), Some(weights)) = (&self.residuals, &self.weights) else {
            return 0;
        };
        (0..residuals.len())
            .filter(|&i| weights[i].abs() * weighting.weight(residuals[i].abs()) > 0.0)
            .count()
    }

    /// Puts the parameters back to a previous state and rebuilds the fit to match.
    fn restore(&mut self, params: &DVector<f64>) {
        let params = params.clone();
        self.set_params(&params);
    }
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for BoundaryFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params = x.clone();

        if let Ok(boundary) = (self.builder)(&self.params) {
            let (residuals, weights) = self.fitting.residuals_and_weights(&boundary);
            self.residuals = Some(residuals);
            self.weights = Some(weights);
            self.current = Some(boundary);
        } else {
            self.residuals = None;
            self.current = None;
        }
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let Some(residuals) = &self.residuals else {
            return None;
        };
        let Some(weights) = &self.weights else {
            return None;
        };
        let mut res = DVector::zeros(residuals.len());
        for i in 0..residuals.len() {
            res[i] = residuals[i] * self.scale(weights, i);
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let Some(residuals) = &self.residuals else {
            return None;
        };
        let Some(weights) = &self.weights else {
            return None;
        };

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

            for i in 0..residuals.len() {
                let d = self.fitting.residual_only(i, &disturbed);
                jac[(i, k)] = self.scale(weights, i) * (d - residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector2;
    use crate::common::linear_space;
    use crate::common::points::{fill_gaps, to_points};
    use crate::common::random_geometry::RandomGeometry2;
    use crate::geom2::{BoundaryData2, BoundaryEditor};
    use approx::assert_relative_eq;

    /// A closed counter-clockwise unit square from (0,0) to (1,1). Boundary normals point to the
    /// right of the tangent, which for counter-clockwise winding means outward.
    fn ccw_square() -> Boundary2 {
        let mut bdata = BoundaryData2::new_closed();
        let mut cursor = bdata.get_cursor(None);
        cursor.add_seg_xy(0.0, 0.0);
        cursor.add_seg_xy(1.0, 0.0);
        cursor.add_seg_xy(1.0, 1.0);
        cursor.add_seg_xy(0.0, 1.0);
        bdata.try_to_boundary().unwrap()
    }

    #[test]
    fn the_residual_is_signed_by_which_side_of_the_boundary_a_sample_falls_on() {
        let square = ccw_square();

        let outside = Point2::new(0.5, -0.25);
        let inside = Point2::new(0.5, 0.25);

        let (_, m_out) = square.at_closest_to_point(&outside);
        let (_, m_in) = square.at_closest_to_point(&inside);

        let r_out = signed_residual(&m_out, &outside);
        let r_in = signed_residual(&m_in, &inside);

        assert!(
            r_out > 0.0,
            "a point outside should be positive, got {r_out}"
        );
        assert!(r_in < 0.0, "a point inside should be negative, got {r_in}");

        // The magnitude is still the plain distance to the boundary.
        assert_relative_eq!(r_out.abs(), 0.25, epsilon = 1e-12);
        assert_relative_eq!(r_in.abs(), 0.25, epsilon = 1e-12);
    }

    #[test]
    fn signing_leaves_the_size_of_the_residual_alone() {
        // The objective squares the residual, so signing it cannot move the minimum. This pins
        // that down directly: the magnitude at every sample matches the unsigned distance the
        // residual used to be.
        let square = ccw_square();
        let samples = fill_gaps(
            &to_points(&[
                [-0.4, -0.4],
                [1.4, -0.4],
                [1.4, 1.4],
                [-0.4, 1.4],
                [-0.4, -0.4],
            ]),
            0.1,
        );

        for p in &samples {
            let (_, m) = square.at_closest_to_point(p);
            assert_relative_eq!(signed_residual(&m, p).abs(), dist(&m, p), epsilon = 1e-12);
        }
    }

    // =========================================================================================
    // Robust refinement
    // =========================================================================================

    /// The triangle fixture of `simple_triangle`, as (points, builder, initial).
    fn triangle_case() -> (Vec<Point2>, BndBuildFn, DVector<f64>) {
        let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
        let points = fill_gaps(&corners, 0.1);

        let builder: BndBuildFn = Box::new(|params: &DVector<f64>| {
            let mut bdata = BoundaryData2::new_open(Point2::new(params[0], params[1]));
            let mut cursor = bdata.get_cursor(None);
            cursor.add_seg_xy(params[2], params[3]);
            cursor.add_seg_xy(params[4], params[5]);
            cursor.add_seg_xy(params[0], params[1]);
            bdata.try_to_boundary()
        });

        let initial = DVector::from(vec![0.9, 0.9, 3.1, 2.1, 2.1, 4.1]);
        (points, builder, initial)
    }

    #[test]
    fn gross_outliers_are_rejected() {
        // Clean samples carrying a little measurement noise, plus a tenth of them thrown well
        // clear of the triangle. The plain fit has no defense against the strays; the robust one
        // should weight them out and land close to the truth.
        //
        // The noise on the inliers is not decoration. MAGSAC separates outliers from inliers by
        // a noise scale estimated from the residuals, so there has to be a visible inlier scale
        // for it to find. Perfectly clean data plus large strays is the pathological case: the
        // opening unweighted solve is dragged, and the only spread left to measure is the damage
        // itself. That is the same trap `MultiAlignOptions2::max_distance` exists to avoid.
        let (mut points, builder, initial) = triangle_case();

        let mut rg = RandomGeometry2::from_seed(0xb0_11d);
        for p in points.iter_mut() {
            *p += Vector2::new(rg.gaussian_f64(0.0, 0.01), rg.gaussian_f64(0.0, 0.01));
        }
        for (k, p) in points.iter_mut().enumerate() {
            if k.is_multiple_of(10) {
                *p += Vector2::new(0.0, 0.5);
            }
        }

        let truth = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
        let err = |r: &BoundaryFitResult| (&r.params - &truth).norm();

        let naive = fit_boundary_to_points(
            &points,
            &builder,
            initial.clone(),
            false,
            &Default::default(),
        )
        .unwrap();
        let robust = fit_boundary_to_points(
            &points,
            &builder,
            initial.clone(),
            false,
            &BoundaryFitOptions::robust(),
        )
        .unwrap();

        assert!(
            err(&robust) < 0.25 * err(&naive),
            "robust weighting should have suppressed the outliers: naive {}, robust {}",
            err(&naive),
            err(&robust)
        );
        assert_eq!(robust.refinement_rounds(), 4);
        assert!(robust.halt().is_none());

        // A noise bound supplied outright does as well, and does not depend on the estimate being
        // taken from a solve the outliers have already influenced.
        let explicit = fit_boundary_to_points(
            &points,
            &builder,
            initial,
            false,
            &BoundaryFitOptions {
                sigma_max: Some(0.02),
                ..BoundaryFitOptions::robust()
            },
        )
        .unwrap();
        assert!(err(&explicit) < 0.25 * err(&naive));
    }

    #[test]
    fn refinement_is_off_by_default() {
        let (points, builder, initial) = triangle_case();
        let result =
            fit_boundary_to_points(&points, &builder, initial, false, &Default::default()).unwrap();

        assert_eq!(result.refinement_rounds(), 0);
        assert_eq!(result.solves().len(), 1);
        assert!(result.halt().is_none());
    }

    #[test]
    fn refinement_on_essentially_exact_data_leaves_the_fit_alone() {
        // The fixture samples sit on the triangle outright, so a converged fit leaves residuals
        // at the level of floating point noise. Reweighting has nothing to find, and the rounds
        // must not disturb an answer that is already right.
        //
        // Worth knowing: the noise estimate does not bail out here. It only rejects a spread that
        // is zero or non-finite, and residuals of order 1e-16 are neither, so the rounds run on
        // numerical noise rather than reporting `RefinementHalt::NoNoiseEstimate`. That is shared
        // behavior with every alignment solver, not something specific to boundary fitting.
        let (points, builder, initial) = triangle_case();
        let result = fit_boundary_to_points(
            &points,
            &builder,
            initial,
            false,
            &BoundaryFitOptions::robust(),
        )
        .unwrap();

        let truth = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
        assert_relative_eq!(result.params, truth, epsilon = 1.0e-5);
        assert!(
            result.halt().is_none(),
            "unexpected halt: {:?}",
            result.halt()
        );
    }

    #[test]
    fn invalid_options_are_rejected() {
        let (points, builder, initial) = triangle_case();

        let bad_patience = BoundaryFitOptions {
            patience: 0,
            ..Default::default()
        };
        assert!(
            fit_boundary_to_points(&points, &builder, initial.clone(), false, &bad_patience)
                .is_err()
        );

        let bad_sigma = BoundaryFitOptions {
            sigma_max: Some(-1.0),
            ..BoundaryFitOptions::robust()
        };
        assert!(fit_boundary_to_points(&points, &builder, initial, false, &bad_sigma).is_err());
    }

    /// This and `surface_normal_line` assert converged parameters to 1e-6, so between them they
    /// are the regression net for the claim that signing the residual does not move a fit.
    #[test]
    fn simple_triangle() {
        let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
        let points = fill_gaps(&corners, 0.1);

        let builder: BndBuildFn = Box::new(|params: &DVector<f64>| {
            let mut bdata = BoundaryData2::new_open(Point2::new(params[0], params[1]));
            let mut cursor = bdata.get_cursor(None);
            cursor.add_seg_xy(params[2], params[3]);
            cursor.add_seg_xy(params[4], params[5]);
            cursor.add_seg_xy(params[0], params[1]);
            bdata.try_to_boundary()
        });

        let initial = DVector::from(vec![0.0, 0.0, 4.0, 0.0, 1.0, 7.0]);
        let result =
            fit_boundary_to_points(&points, &builder, initial, false, &Default::default()).unwrap();
        let expected = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
        assert_relative_eq!(result.params, expected, epsilon = 1.0e-6);
    }

    #[test]
    fn surface_normal_line() {
        let good = linear_space(0.0, 10.0, 10)
            .iter()
            .map(|x| SurfacePoint2::new(Point2::new(*x, 0.0), Vector2::y_axis()))
            .collect::<Vec<_>>();
        let bad = linear_space(0.0, 10.0, 10)
            .iter()
            .map(|x| SurfacePoint2::new(Point2::new(*x, 1.0), Vector2::x_axis()))
            .collect::<Vec<_>>();

        let samples = [good, bad].concat();

        let builder: BndBuildFn = Box::new(|params: &DVector<f64>| {
            let mut bdata = BoundaryData2::new_open(Point2::new(0.0, params[0]));
            let mut cursor = bdata.get_cursor(None);
            cursor.add_seg_xy(10.0, params[1]);
            bdata.try_to_boundary()
        });
        let initial = DVector::from(vec![1.0, 1.0]);

        let result = fit_boundary_to_surface_points(
            &samples,
            &builder,
            initial,
            VecDot::Abs,
            false,
            &Default::default(),
        )
        .unwrap();
        let expected = DVector::from(vec![0.0, 0.0]);
        assert_relative_eq!(result.params, expected, epsilon = 1.0e-6);
    }
}
