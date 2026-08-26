//! A simultaneous alignment of several curves to each other in one combined Levenberg-Marquardt
//! minimization.
//!
//! This is a bundle adjustment rather than a pose graph optimization: it carries a single
//! transformation for each curve except one, which is held fixed, and solves for all of them at
//! once against the raw correspondences. It was built to register metrology-quality profiles of
//! objects with unambiguous morphology, and it works best on low-noise curves which have already
//! been pre-aligned close to each other with substantial overlap.
//!
//! This is the 2D counterpart of `geom3::align3::multi_mesh`, and is deliberately structured
//! identically to it.
//!
//! # How a correspondence constrains two bodies
//!
//! Every alignment point belongs to one curve and is matched against another, and *both* of those
//! curves are moving. A single residual therefore contributes two blocks to its jacobian row: the
//! forward derivative with respect to the test curve's parameters, and the reverse derivative with
//! respect to the reference curve's. See [`point_surf_jacobian2`] and [`point_surf_jacobian2_rev`].
//!
//! # Coordinate frames
//!
//! Everything is measured in world coordinates. A test point is moved by its own curve's
//! transform, the query is pushed into the reference curve's frame only long enough to find the
//! closest point on it, and that match is brought back out to the world before the residual is
//! taken.
//!
//! This matters because each curve's [`AlignValues2`] describes its transform *to the world*, so
//! the residual and both jacobian blocks have to live there too for the derivatives to be
//! consistent with it.
//!
//! # The distance gate is not optional
//!
//! [`MultiAlignOptions2::max_distance`] must be supplied, and there is no `Default` for the
//! options precisely so that it cannot be forgotten. See its documentation for the measurements
//! behind that, which were taken on the 3D side but describe a property of the bundle rather than
//! of the dimension.
//!
//! # A known mismatch with the robust weighting, inherited from the single-body solver
//!
//! [`RESIDUAL_DOF`] is 2 here, matching `points_to_surface2`, on the grounds that the residual is
//! a Euclidean point-to-curve distance in the plane. That describes the arithmetic rather than the
//! distribution. What sets the degrees of freedom is how many independent dimensions of noise
//! survive into the residual, and against a curve the answer is usually one: noise which slides a
//! test point along the local tangent does not change its distance to the curve at all, because
//! the closest point slides along to follow it. Only the normal component contributes. The
//! residual regains its second dimension only where the projection clamps, at a vertex or at the
//! end of an open curve.
//!
//! The practical effect is that MAGSAC++ treats residuals as drawn from a wider distribution than
//! they are, so robust refinement down-weights outliers less aggressively than it should. The
//! geometry is unaffected. Fixing it means changing the weighting, and `MagsacWeight` will not
//! accept a dof of 1 as written, so it is recorded here rather than tuned around. The same is true
//! of the 3D path, and any fix should be judged against both.

use crate::Result;
use crate::common::align::{RefinementHalt, SolveQuality, TerminationReason};
use crate::common::consensus::weights::MagsacWeight;
use crate::common::points::dist;
use crate::geom2::align2::curve::{AlignmentCurve, CurveSurfPoint};
use crate::geom2::align2::jacobian::{point_surf_jacobian2, point_surf_jacobian2_rev};
use crate::geom2::align2::{AlignValues2, Dof3, MultiAlignParams2};
use crate::geom2::{Alignment2, MultiAlignOutcome2, Point2, SurfacePoint2, Vector2};
use crate::na::{DVector, Dyn, Matrix, Owned, U1, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use rayon::prelude::*;

/// The residual degrees of freedom for the MAGSAC++ weight function. The residual is a full
/// Euclidean point-to-curve distance in the plane, so it follows a chi distribution with two
/// degrees of freedom. See the module documentation for why this overstates the true figure.
const RESIDUAL_DOF: usize = 2;

/// The scale factor that turns a median absolute deviation into a consistent estimate of the
/// standard deviation of normally distributed data.
const MAD_TO_SIGMA: f64 = 1.4826;

// ================================================================================================
// Options
// ================================================================================================

/// Options controlling a simultaneous multi-curve alignment.
///
/// Construct with [`MultiAlignOptions2::new`], which requires the correspondence distance gate,
/// and override the remaining fields you care about. The robust-estimation fields mirror
/// `AlignOptions2` and behave identically; the two gates at the top are specific to multi-body
/// work.
///
/// There is deliberately no `Default`. See [`MultiAlignOptions2::max_distance`] for why the gate
/// cannot be left to a default.
#[derive(Clone, Copy, Debug)]
pub struct MultiAlignOptions2 {
    /// A hard cutoff on correspondence distance, in the units of the geometry. A point whose
    /// closest match on the reference curve is further than this contributes nothing.
    ///
    /// This is required rather than optional, and that is a deliberate reversal of the earlier
    /// advice to leave it off and let robust weighting do the work. On real multi-scan data that
    /// advice is wrong, and expensively so.
    ///
    /// Bodies overlap only partially, so a point in a region the reference never saw has no
    /// meaningful match at any distance. The alignment opens with an *unweighted* least squares
    /// solve, in order to have residuals from which to estimate a noise scale, and unweighted least
    /// squares has no defense whatsoever against a gross outlier. A handful of correspondences
    /// sampled across non-overlapping regions is enough to drag the whole bundle.
    ///
    /// Measured on the 3D path, over seventeen partially overlapping laser scans of a part 127 mm
    /// across, starting from an already converged alignment: with no gate the initial solve moved
    /// the group by 0.83 mm and left ten times as many points beyond a tenth of a millimetre of the
    /// reference. Any gate at all fixed it, taking the drift to 0.09 mm. A gate of 1.0 mm and one
    /// of 0.5 mm gave identical results, which says the offending correspondences were more than a
    /// millimetre out rather than marginal. Robust refinement afterwards did not rescue it, because
    /// the noise scale is estimated by a median absolute deviation from that same bad solve: the
    /// estimate describes the well-fitted core, and MAGSAC then treats the displaced body's
    /// correspondences as outliers and stops pulling on them, entrenching the error it inherited.
    ///
    /// Curves are, if anything, more exposed to this than meshes. A profile sliced out of a scan
    /// carries fewer correspondences than a mesh does, so a given number of bad ones is a larger
    /// share of the total.
    ///
    /// Choose it from the geometry rather than from the expected residual. It only has to be
    /// tighter than the scale of a spurious match, so err generous: a value far larger than the
    /// alignment error still removes the correspondences that matter.
    ///
    /// Note that, like `ignore_off_target` in the single-body case, this is a hard gate
    /// re-evaluated as the points move, so it makes the objective piecewise.
    pub max_distance: f64,

    /// An optional cutoff on the angle, in radians, between a test point's normal and its match's
    /// normal. A correspondence exceeding it contributes nothing.
    ///
    /// This suppresses matches onto the far side of a thin wall, where the geometry is close but
    /// facing the wrong way. `None` accepts any orientation.
    pub max_normal_angle: Option<f64>,

    /// The number of iteratively reweighted refinement rounds to perform after the initial
    /// unweighted solve. Zero disables robust weighting entirely.
    pub refinement_steps: usize,

    /// The MAGSAC++ upper noise bound, estimated from the initial residuals via the median
    /// absolute deviation when `None`. See `AlignOptions2::sigma_max` for the units.
    pub sigma_max: Option<f64>,

    /// The Levenberg-Marquardt evaluation budget, as a multiplier on the parameter count.
    ///
    /// A multi-body solve carries `3 * (n - 1)` parameters, so this budget stretches a long way
    /// with many curves. A solve which exhausts it is not a failure: the alignments are kept and
    /// the outcome reports [`crate::common::SolveQuality::Unconverged`]. That is a common outcome
    /// here, since correspondences flip between every overlapping pair.
    ///
    /// Must be greater than zero.
    pub patience: usize,

    /// An optional degree-of-freedom constraint applied to every non-static curve.
    pub dof: Option<Dof3>,
}

impl MultiAlignOptions2 {
    /// Options with the required correspondence distance gate and defaults for everything else.
    ///
    /// `max_distance` has no safe default and so has to be supplied; see its documentation for the
    /// measurements behind that. Everything else follows the single-body defaults: four rounds of
    /// robust refinement with the noise scale estimated from the data, no normal gate, no degree of
    /// freedom locks.
    pub fn new(max_distance: f64) -> Self {
        Self {
            max_distance,
            max_normal_angle: None,
            refinement_steps: 4,
            sigma_max: None,
            patience: 100,
            dof: None,
        }
    }
}

// ================================================================================================
// Alignment points
// ================================================================================================

/// A single correspondence in a multi-curve adjustment: a sample point on one curve which is being
/// matched against another curve.
#[derive(Clone, Debug)]
pub struct MulCurveAlignPoint {
    /// The index of the curve this point was sampled from.
    pub curve_i: usize,

    /// The sample point, in the coordinates of its own curve.
    pub cp: CurveSurfPoint,

    /// The index of the curve this point is being matched against.
    pub ref_i: usize,

    /// A base weight for the correspondence, applied on top of any robust weighting.
    pub weight: f64,
}

impl MulCurveAlignPoint {
    pub fn new(curve_i: usize, cp: CurveSurfPoint, ref_i: usize, weight: f64) -> Self {
        Self {
            curve_i,
            cp,
            ref_i,
            weight,
        }
    }
}

// ================================================================================================
// Entry points
// ================================================================================================

/// Runs a simultaneous alignment over an explicit set of correspondences.
///
/// Use this when the correspondences have already been chosen, whether by sampling them yourself
/// or by pruning an existing set.
///
/// # Arguments
///
/// * `curves`: the bodies taking part, each with its own optional uncertainty, initial pose, and
///   weight providers
/// * `points`: the correspondences, each naming the curve it was sampled from and the curve it is
///   matched against
/// * `static_i`: the index of the curve to hold fixed, which becomes the frame every result is
///   expressed in
/// * `opts`: the solver options, which must carry a correspondence distance gate
///
/// # Failure
///
/// As with the single-body solver, `Err` is reserved for the case where there is no answer at
/// all: rejected arguments, or an initial solve that broke down. A solve which merely exhausts
/// its budget is reported on the outcome and its alignments kept.
pub fn multi_curve_adjustment_with_points(
    curves: &[AlignmentCurve],
    points: Vec<MulCurveAlignPoint>,
    static_i: usize,
    opts: &MultiAlignOptions2,
) -> Result<MultiAlignOutcome2> {
    validate(curves, &points, static_i, opts)?;
    solve_bundle(
        MultiCurveProblem::new(curves, points, static_i, opts)?,
        opts,
    )
}

// ================================================================================================
// The solve
// ================================================================================================

fn validate(
    curves: &[AlignmentCurve],
    points: &[MulCurveAlignPoint],
    static_i: usize,
    opts: &MultiAlignOptions2,
) -> Result<()> {
    if opts.patience == 0 {
        return Err("patience must be greater than zero".into());
    }
    if let Some(s) = opts.sigma_max
        && (!s.is_finite() || s <= 0.0)
    {
        return Err(format!("sigma_max is {s}, but must be finite and strictly positive").into());
    }
    if !opts.max_distance.is_finite() || opts.max_distance <= 0.0 {
        return Err(format!(
            "max_distance is {}, but must be finite and strictly positive",
            opts.max_distance
        )
        .into());
    }
    if static_i >= curves.len() {
        return Err(format!(
            "the static curve index {} is out of range for {} curves",
            static_i,
            curves.len()
        )
        .into());
    }
    if let Some(p) = points
        .iter()
        .find(|p| p.curve_i >= curves.len() || p.ref_i >= curves.len())
    {
        return Err(format!(
            "an alignment point references curve {} against curve {}, but there are only {} curves",
            p.curve_i,
            p.ref_i,
            curves.len()
        )
        .into());
    }

    // Uncertainty is interpolated between the two vertices of the edge a sample lands on, so a
    // short slice would panic partway through a solve rather than here.
    for (i, c) in curves.iter().enumerate() {
        if let Some(u) = c.uncertainty
            && u.len() != c.curve.count()
        {
            return Err(format!(
                "curve {} has {} vertices but {} uncertainty values",
                i,
                c.curve.count(),
                u.len()
            )
            .into());
        }
    }

    Ok(())
}

fn solve_bundle(
    problem: MultiCurveProblem<'_>,
    opts: &MultiAlignOptions2,
) -> Result<MultiAlignOutcome2> {
    let lm = LevenbergMarquardt::<f64>::new().with_patience(opts.patience);
    let n_params = problem.params.param_count();

    let (mut problem, termination) = run(&lm, problem);
    if !SolveQuality::from_termination(&termination).is_usable() {
        return Err(format!("Failed to align curves to each other: {termination:?}").into());
    }

    let mut solves = vec![termination];
    let mut halt = None;

    if opts.refinement_steps > 0 {
        match resolve_sigma_max(opts, &problem) {
            None => halt = Some(RefinementHalt::NoNoiseEstimate),
            Some(sigma_max) => {
                let weighting = MagsacWeight::new(sigma_max, RESIDUAL_DOF);

                for _ in 0..opts.refinement_steps {
                    problem.refresh_inv_sigma();

                    let weighted = problem.count_if_reweighted(&weighting);
                    if weighted < n_params {
                        halt = Some(RefinementHalt::Underdetermined {
                            weighted,
                            params: n_params,
                        });
                        break;
                    }

                    let last_good = problem.params.storage().clone();
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

    Ok(MultiAlignOutcome2::new(problem.alignments(), solves, halt))
}

fn run<'a>(
    lm: &LevenbergMarquardt<f64>,
    problem: MultiCurveProblem<'a>,
) -> (MultiCurveProblem<'a>, TerminationReason) {
    let (result, report) = lm.minimize(problem);
    (result, report.termination)
}

fn resolve_sigma_max(opts: &MultiAlignOptions2, problem: &MultiCurveProblem<'_>) -> Option<f64> {
    match opts.sigma_max {
        Some(s) => Some(s),
        None => estimate_sigma_max(&problem.normalized_residuals()),
    }
}

/// Estimates a MAGSAC++ `sigma_max` from residuals via the median absolute deviation, which is
/// insensitive to the gross outliers the robust weighting exists to suppress.
fn estimate_sigma_max(residuals: &[f64]) -> Option<f64> {
    let center = median(residuals)?;
    let deviations: Vec<f64> = residuals.iter().map(|r| (r - center).abs()).collect();
    let sigma = MAD_TO_SIGMA * median(&deviations)?;
    (sigma.is_finite() && sigma > 0.0).then_some(sigma)
}

fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let n = sorted.len();
    Some(if n.is_multiple_of(2) {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    } else {
        sorted[n / 2]
    })
}

// ================================================================================================
// The problem
// ================================================================================================

/// The per-correspondence state produced by one pass of `move_points`.
struct Moved {
    p_world: Point2,
    c_world: SurfacePoint2,
    residual: f64,
    target_weight: f64,
    target_sigma: f64,
}

struct MultiCurveProblem<'a> {
    curves: &'a [AlignmentCurve<'a>],
    handles: Vec<MulCurveAlignPoint>,
    params: MultiAlignParams2,

    /// The alignment values of every body for the current parameters, cached so the residual and
    /// jacobian loops don't recompute them per correspondence.
    values: Vec<AlignValues2>,

    /// The test point of each correspondence in world coordinates.
    moved: Vec<Point2>,

    /// The matching surface point of each correspondence, also in world coordinates.
    closest: Vec<SurfacePoint2>,

    /// The signed geometric distance of each correspondence, in model units.
    residuals: Vec<f64>,

    /// Each test point's own measurement standard deviation, interpolated from its curve once at
    /// construction because the point never moves within its own curve.
    test_sigma: Vec<f64>,

    /// The reference curve's standard deviation at each current match.
    target_sigma: Vec<f64>,

    /// The reciprocal of the combined test-and-target standard deviation, or 1.0 where there is no
    /// uncertainty at all. Held fixed across a solve and refreshed between them, for the same
    /// reason as in the single-body case.
    inv_sigma: Vec<f64>,

    /// The weight contributed by the correspondence itself: the point's base weight and the two
    /// geometric gates.
    target_weights: Vec<f64>,

    /// The weight contributed by MAGSAC++ reweighting, held fixed across a solve.
    magsac_weights: Vec<f64>,

    max_distance: f64,
    min_normal_dot: Option<f64>,
}

impl<'a> MultiCurveProblem<'a> {
    fn new(
        curves: &'a [AlignmentCurve],
        handles: Vec<MulCurveAlignPoint>,
        static_i: usize,
        opts: &MultiAlignOptions2,
    ) -> Result<Self> {
        let centers = curves
            .iter()
            .map(|c| c.curve.aabb().center())
            .collect::<Vec<_>>();
        let initial = curves.iter().map(|c| c.transform()).collect::<Vec<_>>();
        let params = MultiAlignParams2::from_centers(static_i, &centers, Some(&initial), opts.dof)?;

        let n = handles.len();

        // A test point is fixed within its own curve, so its uncertainty never changes and is
        // resolved once here rather than on every step.
        let test_sigma = handles
            .iter()
            .map(|h| curves[h.curve_i].sigma_at(&h.cp))
            .collect();

        let values = params.compute_all_values();
        let mut item = Self {
            curves,
            handles,
            params,
            values,
            moved: vec![Point2::origin(); n],
            closest: vec![SurfacePoint2::new(Point2::origin(), Vector2::y_axis()); n],
            residuals: vec![0.0; n],
            test_sigma,
            target_sigma: vec![0.0; n],
            inv_sigma: vec![1.0; n],
            target_weights: vec![1.0; n],
            magsac_weights: vec![1.0; n],
            max_distance: opts.max_distance,
            min_normal_dot: opts.max_normal_angle.map(f64::cos),
        };

        item.move_points();
        item.refresh_inv_sigma();
        Ok(item)
    }

    /// Moves every test point into world coordinates, finds its match on the reference curve,
    /// brings that match back out to the world, and recomputes the residual and the geometric
    /// weights.
    fn move_points(&mut self) {
        self.values = self.params.compute_all_values();

        let values = &self.values;
        let curves = self.curves;
        let max_distance = self.max_distance;
        let min_normal_dot = self.min_normal_dot;

        let collected: Vec<Moved> = self
            .handles
            .par_iter()
            .map(|h| {
                let t_test = values[h.curve_i].transform;
                let t_ref = values[h.ref_i].transform;

                // Into the world, then into the reference curve's frame just long enough to find
                // the closest point on it.
                let p_world = t_test * h.cp.sp.point;
                let query = t_ref.inverse_transform_point(&p_world);
                let station = curves[h.ref_i].curve.at_closest_to_point(&query);
                let cp = CurveSurfPoint::from_station(&station);

                // ...and back out to the world, where the residual is measured.
                let c_world = cp.sp.transformed_by(&t_ref);

                let d = dist(&p_world, &c_world.point);
                let residual = d * c_world.scalar_projection(&p_world).signum();

                let mut weight = h.weight;
                if d > max_distance {
                    weight = 0.0;
                }
                if let Some(min_dot) = min_normal_dot {
                    let n_test = t_test.rotation * h.cp.sp.normal;
                    if n_test.dot(&c_world.normal) < min_dot {
                        weight = 0.0;
                    }
                }

                Moved {
                    p_world,
                    c_world,
                    residual,
                    target_weight: weight,
                    target_sigma: curves[h.ref_i].sigma_at(&cp),
                }
            })
            .collect();

        for (i, m) in collected.into_iter().enumerate() {
            self.moved[i] = m.p_world;
            self.closest[i] = m.c_world;
            self.residuals[i] = m.residual;
            self.target_weights[i] = m.target_weight;
            self.target_sigma[i] = m.target_sigma;
        }
    }

    /// Puts the parameters back to a previous state and rebuilds the correspondences to match.
    fn restore(&mut self, storage: &DVector<f64>) {
        self.params.set_storage(storage);
        self.move_points();
    }

    /// Recombines the test and target standard deviations in quadrature.
    fn refresh_inv_sigma(&mut self) {
        for i in 0..self.handles.len() {
            let t = self.test_sigma[i];
            let r = self.target_sigma[i];
            let combined = (t * t + r * r).sqrt();
            self.inv_sigma[i] = if combined > 0.0 { 1.0 / combined } else { 1.0 };
        }
    }

    fn normalized_residuals(&self) -> Vec<f64> {
        self.residuals
            .iter()
            .zip(self.inv_sigma.iter())
            .map(|(r, inv)| r * inv)
            .collect()
    }

    /// The factor applied to both the residual and the jacobian row of a correspondence.
    fn scale(&self, i: usize) -> f64 {
        (self.target_weights[i] * self.magsac_weights[i]).sqrt() * self.inv_sigma[i]
    }

    fn apply_magsac_weights(&mut self, weighting: &MagsacWeight) {
        for i in 0..self.handles.len() {
            let r = (self.residuals[i] * self.inv_sigma[i]).abs();
            self.magsac_weights[i] = weighting.weight(r);
        }
    }

    fn count_if_reweighted(&self, weighting: &MagsacWeight) -> usize {
        (0..self.handles.len())
            .filter(|&i| {
                let r = (self.residuals[i] * self.inv_sigma[i]).abs();
                self.target_weights[i] * weighting.weight(r) > 0.0
            })
            .count()
    }

    /// Builds one [`Alignment2`] per curve, with that curve's own correspondence residuals.
    fn alignments(&self) -> Vec<Alignment2> {
        let mut grouped = vec![Vec::new(); self.curves.len()];
        for (i, h) in self.handles.iter().enumerate() {
            grouped[h.curve_i].push(self.residuals[i]);
        }

        self.values
            .iter()
            .enumerate()
            .zip(grouped)
            .map(|((i, v), residuals)| {
                let body = self.params.body(i);
                Alignment2::new(v.transform, v.align, body.local, body.offset, residuals)
            })
            .collect()
    }
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for MultiCurveProblem<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params.set_storage(x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.storage().clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.handles.len());
        for i in 0..self.handles.len() {
            res[i] = self.residuals[i] * self.scale(i);
        }
        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            self.handles.len(),
            self.params.param_count(),
        );

        for (i, h) in self.handles.iter().enumerate() {
            let scale = self.scale(i);
            let p = &self.moved[i];
            let c = &self.closest[i];

            // The correspondence constrains the test curve...
            let fwd = point_surf_jacobian2(p, c, &self.values[h.curve_i]) * scale;
            self.params.add_jacobian_block(&mut jac, i, h.curve_i, &fwd);

            // ...and the reference curve, which is also free to move.
            let rev = point_surf_jacobian2_rev(p, c, &self.values[h.ref_i]) * scale;
            self.params.add_jacobian_block(&mut jac, i, h.ref_i, &rev);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Curve2, Iso2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// A correspondence gate wider than the test rectangles themselves, so that the synthetic
    /// cases exercise the solver rather than the gate. Cases which are about the gate set their
    /// own.
    const WIDE_GATE: f64 = 20.0;

    fn test_opts() -> MultiAlignOptions2 {
        MultiAlignOptions2::new(WIDE_GATE)
    }

    /// A closed counter-clockwise rectangle 10 by 6, centered on the origin, so its normals point
    /// outward. Its corners sit at arc lengths 0, 10, 16 and 26 out of a perimeter of 32.
    fn rect_curve() -> Curve2 {
        let points =
            [(-5.0, -3.0), (5.0, -3.0), (5.0, 3.0), (-5.0, 3.0)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    /// Correspondences from every curve onto the one before it, sampled at a uniform arc length
    /// spacing. The half-spacing offset keeps samples off the corners, where the normal of a
    /// polyline is ambiguous.
    fn chain_points(curves: &[Curve2], spacing: f64) -> Vec<MulCurveAlignPoint> {
        let mut points = Vec::new();
        for (i, curve) in curves.iter().enumerate().skip(1) {
            let n = (curve.length() / spacing).floor() as usize;
            for k in 0..n {
                let l = (k as f64 + 0.5) * spacing;
                if let Some(station) = curve.at_length(l) {
                    let cp = CurveSurfPoint::from_station(&station);
                    points.push(MulCurveAlignPoint::new(i, cp, i - 1, 1.0));
                }
            }
        }
        points
    }

    fn alignment_curves<'a>(curves: &'a [Curve2], initial: &'a [Iso2]) -> Vec<AlignmentCurve<'a>> {
        curves
            .iter()
            .zip(initial.iter())
            .map(|(c, t)| AlignmentCurve::new(c, None, Some(t), None))
            .collect()
    }

    // ============================================================================================
    // Recovery
    // ============================================================================================

    #[test]
    fn a_disturbed_curve_is_brought_back_onto_its_twin() -> Result<()> {
        // Two identical rectangles, one of them displaced. The adjustment should undo the
        // displacement, since the two curves coincide when aligned.
        let curves = vec![rect_curve(), rect_curve()];
        let disturb = Iso2::new(Vector2::new(0.2, -0.15), PI / 90.0);
        let initial = vec![Iso2::identity(), disturb];

        let acs = alignment_curves(&curves, &initial);
        let points = chain_points(&curves, 0.4);

        let outcome = multi_curve_adjustment_with_points(&acs, points, 0, &test_opts())?;

        assert_eq!(outcome.len(), 2);
        // Curve 0 is static and must not have moved at all.
        assert_relative_eq!(
            outcome.alignment(0).full_transform().to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-12
        );
        // Curve 1 should have come back to the identity, undoing its initial displacement.
        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-6
        );

        Ok(())
    }

    #[test]
    fn three_curves_settle_together() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve(), rect_curve()];
        let initial = vec![
            Iso2::identity(),
            Iso2::new(Vector2::new(0.15, -0.1), PI / 120.0),
            Iso2::new(Vector2::new(-0.1, 0.12), -PI / 100.0),
        ];

        let acs = alignment_curves(&curves, &initial);
        let points = chain_points(&curves, 0.5);

        let outcome = multi_curve_adjustment_with_points(&acs, points, 0, &test_opts())?;

        for i in 0..3 {
            assert_relative_eq!(
                outcome.alignment(i).full_transform().to_matrix(),
                Iso2::identity().to_matrix(),
                epsilon = 1e-5
            );
        }

        Ok(())
    }

    #[test]
    fn the_static_curve_never_moves() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::translation(0.3, 0.0), Iso2::translation(0.0, 0.25)];
        let acs = alignment_curves(&curves, &initial);

        // Make curve 1 the static one this time, so the check is not trivially about index zero.
        let outcome =
            multi_curve_adjustment_with_points(&acs, chain_points(&curves, 0.5), 1, &test_opts())?;

        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            initial[1].to_matrix(),
            epsilon = 1e-12
        );

        Ok(())
    }

    #[test]
    fn residuals_are_grouped_by_the_curve_they_came_from() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(); 3];
        let acs = alignment_curves(&curves, &initial);
        let points = chain_points(&curves, 0.6);

        let from_1 = points.iter().filter(|p| p.curve_i == 1).count();
        let from_2 = points.iter().filter(|p| p.curve_i == 2).count();

        let outcome = multi_curve_adjustment_with_points(&acs, points, 0, &test_opts())?;

        // Curve 0 sourced no correspondences, so it has no residuals of its own.
        assert_eq!(outcome.alignment(0).residuals().len(), 0);
        assert_eq!(outcome.alignment(1).residuals().len(), from_1);
        assert_eq!(outcome.alignment(2).residuals().len(), from_2);

        Ok(())
    }

    // ============================================================================================
    // Robustness and reporting
    // ============================================================================================

    #[test]
    fn locked_dof_are_honored() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(), Iso2::translation(0.3, 0.0)];
        let acs = alignment_curves(&curves, &initial);

        let opts = MultiAlignOptions2 {
            dof: Some(Dof3::new(false, true, true)),
            ..test_opts()
        };
        let outcome =
            multi_curve_adjustment_with_points(&acs, chain_points(&curves, 0.5), 0, &opts)?;

        // tx is locked, so the x displacement cannot be undone and must survive untouched.
        assert_relative_eq!(
            outcome.alignment(1).full_transform().translation.vector.x,
            0.3,
            epsilon = 1e-9
        );

        Ok(())
    }

    #[test]
    fn exhausted_budget_is_reported_not_raised() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![
            Iso2::identity(),
            Iso2::new(Vector2::new(1.5, -1.0), PI / 12.0),
        ];
        let acs = alignment_curves(&curves, &initial);

        let opts = MultiAlignOptions2 {
            patience: 1,
            ..test_opts()
        };
        let outcome =
            multi_curve_adjustment_with_points(&acs, chain_points(&curves, 0.5), 0, &opts)?;

        assert_eq!(outcome.quality(), SolveQuality::Unconverged);
        assert!(outcome.solves().contains(&TerminationReason::LostPatience));
        // The alignments are still real geometry rather than something degenerate.
        for a in outcome.alignments() {
            assert!(a.full_transform().to_matrix().iter().all(|v| v.is_finite()));
        }

        Ok(())
    }

    #[test]
    fn gross_outliers_are_rejected() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let disturb = Iso2::new(Vector2::new(0.2, -0.15), PI / 120.0);
        let initial = vec![Iso2::identity(), disturb];
        let acs = alignment_curves(&curves, &initial);

        // Corrupt a tenth of the correspondences by throwing their sample points off the curve.
        let mut points = chain_points(&curves, 0.4);
        for (k, p) in points.iter_mut().enumerate() {
            if k.is_multiple_of(10) {
                p.cp.sp.point += Vector2::new(0.0, 5.0);
            }
        }

        let naive = multi_curve_adjustment_with_points(
            &acs,
            points.clone(),
            0,
            &MultiAlignOptions2 {
                refinement_steps: 0,
                ..test_opts()
            },
        )?;
        let robust = multi_curve_adjustment_with_points(&acs, points, 0, &test_opts())?;

        let error = |o: &MultiAlignOutcome2| {
            (o.alignment(1).full_transform().to_matrix() - Iso2::identity().to_matrix()).norm()
        };

        assert!(
            error(&robust) < error(&naive) * 0.5,
            "robust weighting should have suppressed the outliers: naive {}, robust {}",
            error(&naive),
            error(&robust)
        );

        Ok(())
    }

    /// The initial unweighted solve is dragged by correspondences into regions the other curve
    /// never covered, which is why `max_distance` is required rather than optional.
    ///
    /// This is the unit-scale form of a failure found on real 3D data: seventeen partially
    /// overlapping laser scans, starting from an already correct alignment, were moved most of a
    /// millimetre by a handful of such correspondences, and the robust refinement that followed
    /// entrenched the error rather than undoing it. See the module documentation.
    #[test]
    fn without_a_distance_gate_the_unweighted_solve_is_dragged() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(), Iso2::identity()];
        let acs = alignment_curves(&curves, &initial);

        // Both bodies start on the answer, so any motion at all is damage. A twentieth of the
        // correspondences are thrown well clear of the curve, standing in for points sampled over
        // a region the two curves never both saw.
        let mut points = chain_points(&curves, 0.4);
        for (k, p) in points.iter_mut().enumerate() {
            if k.is_multiple_of(20) {
                p.cp.sp.point += Vector2::new(0.0, 8.0);
            }
        }

        let unweighted = |gate: f64| MultiAlignOptions2 {
            refinement_steps: 0,
            ..MultiAlignOptions2::new(gate)
        };
        let drift = |o: &MultiAlignOutcome2| {
            (o.alignment(1).full_transform().to_matrix() - Iso2::identity().to_matrix()).norm()
        };

        let ungated =
            multi_curve_adjustment_with_points(&acs, points.clone(), 0, &unweighted(WIDE_GATE))?;
        let gated = multi_curve_adjustment_with_points(&acs, points, 0, &unweighted(1.0))?;

        assert!(
            drift(&ungated) > 1e-3,
            "the ungated solve should have been pulled off the answer, but moved only {}",
            drift(&ungated)
        );
        assert!(
            drift(&gated) < drift(&ungated) * 0.1,
            "the gate should have held the solve in place: ungated {}, gated {}",
            drift(&ungated),
            drift(&gated)
        );

        Ok(())
    }

    #[test]
    fn the_distance_gate_excludes_far_correspondences() -> Result<()> {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(), Iso2::identity()];
        let acs = alignment_curves(&curves, &initial);

        let mut points = chain_points(&curves, 0.5);
        let far = points.len() / 2;
        points[far].cp.sp.point += Vector2::new(0.0, 50.0);

        // With a tight gate the distant point contributes nothing, so the fit is unchanged from
        // the identity it already sits at.
        let opts = MultiAlignOptions2 {
            max_distance: 1.0,
            refinement_steps: 0,
            ..test_opts()
        };
        let outcome = multi_curve_adjustment_with_points(&acs, points, 0, &opts)?;

        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-6
        );

        Ok(())
    }

    #[test]
    fn the_normal_gate_excludes_disagreeing_correspondences() -> Result<()> {
        // A point whose normal faces the opposite way from its match is suppressed, which is what
        // keeps a match from landing on the far side of a thin wall.
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(), Iso2::identity()];
        let acs = alignment_curves(&curves, &initial);

        // One correspondence, sampled mid-way along the bottom edge and then flipped so its
        // normal points inward while its match still points outward.
        let station = curves[1].at_length(5.0).unwrap();
        let mut cp = CurveSurfPoint::from_station(&station);
        cp.sp.normal = -cp.sp.normal;

        let opts = MultiAlignOptions2 {
            max_normal_angle: Some(PI / 4.0),
            refinement_steps: 0,
            ..test_opts()
        };
        let points = vec![MulCurveAlignPoint::new(1, cp, 0, 1.0)];
        let problem = MultiCurveProblem::new(&acs, points, 0, &opts)?;

        assert_eq!(problem.target_weights[0], 0.0);

        Ok(())
    }

    // ============================================================================================
    // The jacobian
    // ============================================================================================

    #[test]
    fn jacobian_matches_finite_differences() -> Result<()> {
        // The recovery tests above are insensitive to jacobian errors: a wrong jacobian usually
        // costs Levenberg-Marquardt iterations rather than moving the minimum, so the solve still
        // lands in the right place. This checks the derivative itself.
        //
        // Three bodies are needed, not two. With two, the only reference is the static body, whose
        // jacobian block is dropped anyway, so the reverse derivative would never be exercised.
        // Here body 1 is a reference for body 2's correspondences *and* is free to move, which is
        // the configuration that puts the reverse block to work.
        let curves = vec![rect_curve(), rect_curve(), rect_curve()];
        let initial = vec![
            Iso2::identity(),
            Iso2::new(Vector2::new(0.21, -0.13), 0.03),
            Iso2::new(Vector2::new(-0.17, 0.11), -0.04),
        ];
        let acs = alignment_curves(&curves, &initial);

        // Sample points well away from the corners of the rectangle, so every projection lands in
        // the interior of an edge. The analytic jacobian holds the correspondence fixed, which is
        // right for a locally straight target but not across a corner, and near a corner the
        // finite difference would pick up a correspondence flip the analytic form does not model.
        let probes = [2.0, 5.0, 8.0, 12.0, 14.0, 18.0, 21.0, 24.0, 28.0, 30.0];
        let mut points = Vec::new();
        for body in [1usize, 2] {
            for l in probes {
                let cp = CurveSurfPoint::from_station(&curves[body].at_length(l).unwrap());
                points.push(MulCurveAlignPoint::new(body, cp, body - 1, 1.0));
            }
        }

        let opts = MultiAlignOptions2 {
            refinement_steps: 0,
            ..test_opts()
        };
        let mut problem = MultiCurveProblem::new(&acs, points, 0, &opts)?;

        let x0 = problem.params.storage().clone();
        let analytic = problem.jacobian().unwrap();
        assert_eq!(
            analytic.ncols(),
            6,
            "two free bodies, three parameters each"
        );

        let eps = 1e-7;
        for k in 0..x0.len() {
            let mut lo = x0.clone();
            lo[k] -= eps;
            problem.set_params(&lo);
            let r_lo = problem.residuals().unwrap();

            let mut hi = x0.clone();
            hi[k] += eps;
            problem.set_params(&hi);
            let r_hi = problem.residuals().unwrap();

            for row in 0..r_lo.len() {
                let numeric = (r_hi[row] - r_lo[row]) / (2.0 * eps);
                assert_relative_eq!(analytic[(row, k)], numeric, epsilon = 1e-5);
            }
        }

        Ok(())
    }

    #[test]
    fn a_reference_body_receives_jacobian_columns() -> Result<()> {
        // A narrower statement of the same thing, which fails loudly if the reverse block is ever
        // dropped: correspondences sourced from body 2 must produce nonzero derivatives in body
        // 1's columns, because body 1 is what they are matched against.
        let curves = vec![rect_curve(), rect_curve(), rect_curve()];
        let initial = vec![
            Iso2::identity(),
            Iso2::translation(0.2, 0.0),
            Iso2::translation(0.0, 0.2),
        ];
        let acs = alignment_curves(&curves, &initial);

        let cp = CurveSurfPoint::from_station(&curves[2].at_length(5.0).unwrap());
        let points = vec![MulCurveAlignPoint::new(2, cp, 1, 1.0)];

        let problem = MultiCurveProblem::new(&acs, points, 0, &test_opts())?;
        let jac = problem.jacobian().unwrap();

        let body1 = problem.params.column_offset(1).unwrap();
        let body2 = problem.params.column_offset(2).unwrap();

        let block_norm =
            |start: usize| -> f64 { (0..3).map(|k| jac[(0, start + k)].abs()).sum::<f64>() };

        assert!(
            block_norm(body2) > 1e-9,
            "the test body should have a nonzero jacobian block"
        );
        assert!(
            block_norm(body1) > 1e-9,
            "the reference body should have a nonzero jacobian block too"
        );

        Ok(())
    }

    #[test]
    fn a_self_correspondence_cancels_to_a_zero_row() -> Result<()> {
        // A point matched against its own curve carries no information: moving that curve moves
        // the test point and its match together, so the distance between them cannot change. The
        // two jacobian blocks land in the same columns and must cancel, which they do because the
        // reverse form is the negation of the forward one when both share a body.
        //
        // This is why the blocks are accumulated rather than assigned. Overwriting would keep only
        // the reverse block and report a spurious sensitivity for a correspondence that has none.
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(), Iso2::translation(0.3, -0.2)];
        let acs = alignment_curves(&curves, &initial);

        let cp = CurveSurfPoint::from_station(&curves[1].at_length(5.0).unwrap());
        let points = vec![MulCurveAlignPoint::new(1, cp, 1, 1.0)];

        let problem = MultiCurveProblem::new(&acs, points, 0, &test_opts())?;
        let jac = problem.jacobian().unwrap();

        for k in 0..jac.ncols() {
            assert_relative_eq!(jac[(0, k)], 0.0, epsilon = 1e-12);
        }

        Ok(())
    }

    // ============================================================================================
    // Validation
    // ============================================================================================

    #[test]
    fn invalid_options_and_indices_are_rejected() {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(); 2];
        let acs = alignment_curves(&curves, &initial);
        let points = chain_points(&curves, 1.0);

        let bad_patience = MultiAlignOptions2 {
            patience: 0,
            ..test_opts()
        };
        assert!(
            multi_curve_adjustment_with_points(&acs, points.clone(), 0, &bad_patience).is_err()
        );

        let bad_sigma = MultiAlignOptions2 {
            sigma_max: Some(-1.0),
            ..test_opts()
        };
        assert!(multi_curve_adjustment_with_points(&acs, points.clone(), 0, &bad_sigma).is_err());

        let bad_distance = MultiAlignOptions2 {
            max_distance: 0.0,
            ..test_opts()
        };
        assert!(
            multi_curve_adjustment_with_points(&acs, points.clone(), 0, &bad_distance).is_err()
        );

        // A static index past the end of the curve list.
        assert!(multi_curve_adjustment_with_points(&acs, points, 5, &test_opts()).is_err());
    }

    #[test]
    fn out_of_range_alignment_points_are_rejected() {
        let curves = vec![rect_curve(), rect_curve()];
        let initial = vec![Iso2::identity(); 2];
        let acs = alignment_curves(&curves, &initial);

        let mut points = chain_points(&curves, 1.0);
        points[0].ref_i = 7;

        let err = multi_curve_adjustment_with_points(&acs, points, 0, &test_opts())
            .unwrap_err()
            .to_string();
        assert!(err.contains("only 2 curves"), "unexpected message: {err}");
    }

    #[test]
    fn a_mismatched_uncertainty_slice_is_rejected() {
        // Uncertainty is interpolated between the vertices of the edge a sample lands on, so a
        // short slice would panic partway through a solve. It has to be caught up front.
        let curves = vec![rect_curve(), rect_curve()];
        let initial = [Iso2::identity(); 2];
        let short = [0.01; 2];

        let acs = vec![
            AlignmentCurve::new(&curves[0], None, Some(&initial[0]), None),
            AlignmentCurve::new(&curves[1], Some(&short), Some(&initial[1]), None),
        ];

        let err =
            multi_curve_adjustment_with_points(&acs, chain_points(&curves, 1.0), 0, &test_opts())
                .unwrap_err()
                .to_string();
        assert!(
            err.contains("uncertainty values"),
            "unexpected message: {err}"
        );
    }
}
