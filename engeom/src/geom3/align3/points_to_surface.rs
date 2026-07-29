use crate::common::SPCoords;
use crate::common::align::{RefinementHalt, SolveQuality, TerminationReason};
use crate::common::consensus::weights::MagsacWeight;
use crate::common::dist;
use crate::geom3::align3::jacobian::{copy_jacobian, point_surf_jacobian};
use crate::geom3::align3::{
    AlignOptions3, AlignParams3, AlignSurfMatch3, AlignValues3, SurfaceTarget3,
};
use crate::geom3::{AlignOutcome3, Alignment3};
use crate::na::{Dyn, Matrix, Owned, U1, U6, Vector};
use crate::{Point3, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

/// The residual degrees of freedom for the MAGSAC++ weight function. The residual here is a full
/// Euclidean point-to-target distance in space, so it follows a chi distribution with three
/// degrees of freedom. (`MagsacWeight` requires at least 2; a point-to-plane residual would be a
/// one-dimensional projection and would need the weight function extended.)
const RESIDUAL_DOF: usize = 3;

/// The scale factor that turns a median absolute deviation into a consistent estimate of the
/// standard deviation of normally distributed data.
const MAD_TO_SIGMA: f64 = 1.4826;

/// The number of free parameters in a 3D alignment (tx, ty, tz, rx, ry, rz).
const N_PARAMS: usize = 6;

/// Performs a Levenberg-Marquardt minimization to align a set of 3D points to a surface target.
///
/// The points are repeatedly projected onto their closest position on the target as the solver
/// moves them, so the correspondences are re-established at every step rather than being fixed up
/// front.
///
/// By default this is a robust alignment: an initial unweighted solve is followed by
/// [`AlignOptions3::refinement_steps`] rounds of iteratively reweighted least squares using
/// MAGSAC++ weights. See [`AlignOptions3`] for the weighting and per-point uncertainty controls.
///
/// # Failure
///
/// An `Err` means there is no usable answer at all: the arguments were rejected, or the initial
/// solve broke down numerically. Everything softer than that is reported on the returned
/// [`AlignOutcome3`] rather than raised as an error, because in each case there is still a real
/// alignment to hand back:
///
/// * A solve which exhausts its evaluation budget leaves behind the best parameters it found, so
///   the alignment is kept and the outcome reports [`crate::common::SolveQuality::Unconverged`].
/// * A refinement round which breaks down is rolled back, and the alignment from the last round
///   that succeeded is returned with [`RefinementHalt::SolveFailed`] recorded on the outcome.
///
/// Callers who do not care about any of this can finish with `?.into_alignment()`.
///
/// # Arguments
///
/// * `points`: the 3D points to be aligned, in their own local coordinate system
/// * `target`: an entity that implements the [`SurfaceTarget3`] trait, such as a `Mesh3`, which
///   will be the stationary target of the points.
/// * `params`: the alignment parameters, see [`AlignParams3`] for details
/// * `opts`: options controlling weighting and filtering, see [`AlignOptions3`]
///
/// returns: Result<AlignOutcome<Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn points_to_surface3(
    points: &[Point3],
    target: &impl SurfaceTarget3,
    params: AlignParams3,
    opts: &AlignOptions3,
) -> Result<AlignOutcome3> {
    validate(points, opts)?;

    let lm = LevenbergMarquardt::<f64>::new().with_patience(opts.patience);

    // The first solve is unweighted. Its residuals are what the sigma_max estimate is drawn from,
    // since residuals taken before any alignment would describe the initial misplacement rather
    // than the noise.
    let (mut problem, termination) =
        solve(&lm, PointsToSurface3::new(points, target, params, opts));

    // This is the one place a solve failure is fatal. There is no earlier result to fall back to,
    // so there is nothing to report.
    if !SolveQuality::from_termination(&termination).is_usable() {
        return Err(format!("Failed to align points to surface: {termination:?}").into());
    }

    let mut solves = vec![termination];
    let mut halt = None;

    if opts.refinement_steps > 0 {
        match resolve_sigma_max(opts, &problem) {
            None => halt = Some(RefinementHalt::NoNoiseEstimate),
            Some(sigma_max) => {
                let weighting = MagsacWeight::new(sigma_max, RESIDUAL_DOF);

                for _ in 0..opts.refinement_steps {
                    // The matches have moved since the last refresh, so any target-side
                    // uncertainty is re-read before it feeds the reweighting.
                    problem.refresh_inv_sigma();

                    // If reweighting would leave the problem underdetermined, keep the last good
                    // result rather than solving something rank-deficient.
                    let weighted = problem.count_if_reweighted(&weighting);
                    if weighted < N_PARAMS {
                        halt = Some(RefinementHalt::Underdetermined {
                            weighted,
                            params: N_PARAMS,
                        });
                        break;
                    }

                    // Refinement improves on an alignment which is already usable, so a round
                    // which breaks down costs nothing as long as its parameters can be discarded.
                    // `AlignParams3` is a pair of isometries and a 6-vector, so keeping a copy of
                    // the last good state is cheaper than the solve it protects.
                    let last_good = problem.params.clone();

                    problem.apply_magsac_weights(&weighting);
                    let (next, termination) = solve(&lm, problem);
                    problem = next;

                    if !SolveQuality::from_termination(&termination).is_usable() {
                        problem.restore(last_good);
                        halt = Some(RefinementHalt::SolveFailed(termination));
                        break;
                    }

                    solves.push(termination);
                }
            }
        }
    }

    let c = problem.params.compute_values();
    let alignment = Alignment3::new(
        c.transform,
        c.align,
        problem.params.local,
        problem.params.offset,
        // The geometric residuals, not the weighted/normalized ones the solver minimizes. These
        // are what a caller wants to inspect: real deviations in the units of the geometry.
        problem.residuals.clone(),
    );

    Ok(AlignOutcome3::new(alignment, solves, halt))
}

/// Checks the arguments which the solver cannot recover from. These are caller mistakes rather
/// than difficult data, so they are raised as errors instead of being reported on the outcome.
fn validate(points: &[Point3], opts: &AlignOptions3) -> Result<()> {
    if opts.patience == 0 {
        return Err("patience must be greater than zero".into());
    }

    if let Some(s) = opts.sigma_max
        && (!s.is_finite() || s <= 0.0)
    {
        return Err(format!("sigma_max is {s}, but must be finite and strictly positive").into());
    }

    if let Some(sigma) = opts.point_sigma {
        if sigma.len() != points.len() {
            return Err(format!(
                "point_sigma has {} entries but there are {} points",
                sigma.len(),
                points.len()
            )
            .into());
        }
        if let Some(i) = sigma.iter().position(|s| !s.is_finite() || *s <= 0.0) {
            return Err(format!(
                "point_sigma[{}] is {}, but every entry must be finite and strictly positive",
                i, sigma[i]
            )
            .into());
        }
    }

    Ok(())
}

/// Runs a single Levenberg-Marquardt solve, returning the problem in whatever state the solver
/// left it along with how it terminated.
///
/// Classifying the termination is deliberately left to the caller: whether a given outcome is
/// fatal depends on whether there is an earlier result to fall back to.
fn solve<'a, T: SurfaceTarget3>(
    lm: &LevenbergMarquardt<f64>,
    problem: PointsToSurface3<'a, T>,
) -> (PointsToSurface3<'a, T>, TerminationReason) {
    let (result, report) = lm.minimize(problem);
    (result, report.termination)
}

/// Determines the MAGSAC++ noise bound, either from the explicit option or by estimating it from
/// the current residuals. Returns `None` when the residual spread is too degenerate to estimate
/// from, in which case robust refinement is skipped.
///
/// An explicit `sigma_max` is validated up front by `validate`, so it is trusted here.
fn resolve_sigma_max<T: SurfaceTarget3>(
    opts: &AlignOptions3,
    problem: &PointsToSurface3<'_, T>,
) -> Option<f64> {
    match opts.sigma_max {
        Some(s) => Some(s),
        None => estimate_sigma_max(&problem.normalized_residuals()),
    }
}

/// Estimates a MAGSAC++ `sigma_max` from a set of residuals via the median absolute deviation.
///
/// MAD is used rather than the standard deviation because it is insensitive to the gross outliers
/// the robust weighting exists to suppress: contaminating up to half the data cannot move it
/// arbitrarily, whereas a single distant point can dominate a standard deviation.
///
/// Returns `None` when the spread is zero or non-finite, which happens when the fit is already
/// essentially exact and there is nothing to reweight.
fn estimate_sigma_max(residuals: &[f64]) -> Option<f64> {
    let center = median(residuals)?;
    let deviations: Vec<f64> = residuals.iter().map(|r| (r - center).abs()).collect();
    let sigma = MAD_TO_SIGMA * median(&deviations)?;

    if sigma.is_finite() && sigma > 0.0 {
        Some(sigma)
    } else {
        None
    }
}

/// The median of a slice of finite values, or `None` if the slice is empty.
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

struct PointsToSurface3<'a, T: SurfaceTarget3> {
    points: &'a [Point3],
    target: &'a T,
    params: AlignParams3,

    /// The alignment values for the current parameters. Cached here so that `jacobian()` doesn't
    /// have to recompute what `move_points` already worked out.
    current: AlignValues3,

    /// The test points after being moved by the current transform, in the target's space.
    moved: Vec<Point3>,

    /// The closest match on the target for each moved point.
    closest: Vec<AlignSurfMatch3>,

    /// The signed geometric distance from each moved point to its match, in model units.
    residuals: Vec<f64>,

    /// Each test point's own measurement standard deviation, or 0.0 when none was supplied.
    test_sigma: Vec<f64>,

    /// The reciprocal of each point's combined test-and-target standard deviation, or 1.0 where
    /// there is no uncertainty information at all. Stored inverted so the hot path multiplies.
    ///
    /// Like `magsac_weights`, this is held fixed across a single Levenberg-Marquardt solve and
    /// refreshed between them. The target's contribution depends on which match a point found, so
    /// recomputing it inside the solve would make the true residual derivative pick up a
    /// `d(sigma)/d(params)` term that the analytic jacobian does not model.
    inv_sigma: Vec<f64>,

    /// The weight contributed by the target itself, recomputed whenever the points move.
    target_weights: Vec<f64>,

    /// The weight contributed by MAGSAC++ reweighting. All 1.0 until the first refinement round,
    /// and held fixed across a single Levenberg-Marquardt solve so the jacobian stays consistent
    /// with the residual it differentiates.
    magsac_weights: Vec<f64>,

    ignore_off_target: bool,
}

impl<'a, T: SurfaceTarget3> PointsToSurface3<'a, T> {
    fn new(
        points: &'a [Point3],
        target: &'a T,
        params: AlignParams3,
        opts: &AlignOptions3,
    ) -> Self {
        let current = params.compute_values();
        let test_sigma = match opts.point_sigma {
            Some(s) => s.to_vec(),
            None => vec![0.0; points.len()],
        };

        let mut x = Self {
            points,
            target,
            params,
            current,
            moved: vec![Point3::origin(); points.len()],
            closest: vec![AlignSurfMatch3::default(); points.len()],
            residuals: vec![0.0; points.len()],
            test_sigma,
            inv_sigma: vec![1.0; points.len()],
            target_weights: vec![1.0; points.len()],
            magsac_weights: vec![1.0; points.len()],
            ignore_off_target: opts.ignore_off_target,
        };

        x.move_points();

        // The target's uncertainty isn't known until the points have found their matches, so the
        // combined scale is resolved once here rather than in the field initializer above.
        x.refresh_inv_sigma();
        x
    }

    /// Moves the points by the current transform, finds the closest position on the target to
    /// each, and recomputes the geometric residuals and the target-contributed weights.
    ///
    /// This deliberately does not touch `magsac_weights`, which must stay fixed for the duration
    /// of a single solve.
    fn move_points(&mut self) {
        self.current = self.params.compute_values();
        let transform = self.current.transform;

        for i in 0..self.points.len() {
            let m = transform * self.points[i];
            let c = self.target.find_align_match(&m);

            // The residual is the distance between the test point and the closest point on the
            // target, adjusted for the direction of the scalar projection.
            self.residuals[i] = dist(&m, &c.point) * c.scalar_projection(&m).signum();

            self.target_weights[i] = if self.ignore_off_target {
                c.weight * f64::from(c.is_on)
            } else {
                c.weight
            };

            self.moved[i] = m;
            self.closest[i] = c;
        }
    }

    /// Puts the problem back to a previous set of parameters, discarding whatever the solver left
    /// behind, and recomputes the matches and residuals to match.
    ///
    /// Used to roll back a refinement round which broke down. The stale weights are left alone
    /// because a rolled-back round is the end of the refinement loop, and the reported residuals
    /// are the geometric ones which `move_points` has just rebuilt.
    fn restore(&mut self, params: AlignParams3) {
        self.params = params;
        self.move_points();
    }

    /// Recomputes the combined uncertainty scale from the test points' own standard deviations
    /// and whatever uncertainty the target reported at each current match.
    ///
    /// The two are independent measurements of position, so the variance of their difference is
    /// the sum of their variances and the standard deviations combine in quadrature. Where
    /// neither side reports any uncertainty the scale falls back to 1.0, leaving the residual as
    /// a plain geometric distance.
    fn refresh_inv_sigma(&mut self) {
        for i in 0..self.points.len() {
            let test = self.test_sigma[i];
            let target = self.closest[i].sigma;
            let combined = (test * test + target * target).sqrt();

            self.inv_sigma[i] = if combined > 0.0 { 1.0 / combined } else { 1.0 };
        }
    }

    /// The residuals expressed as dimensionless multiples of each point's own standard deviation.
    /// Identical to the geometric residuals when no uncertainty was supplied on either side.
    fn normalized_residuals(&self) -> Vec<f64> {
        self.residuals
            .iter()
            .zip(self.inv_sigma.iter())
            .map(|(r, inv)| r * inv)
            .collect()
    }

    /// The factor applied to both the residual and the jacobian row for point `i`.
    ///
    /// This folds together two separate things. The `inv_sigma` term scales the residual into
    /// units of the point's own noise, and being part of the residual definition it must appear
    /// in the derivative too. The `sqrt` of the weight is the standard way to express a weighted
    /// least-squares problem to a solver that minimizes a plain sum of squares: the solver sees
    /// `(sqrt(w) * r)^2`, which is the intended `w * r^2`.
    fn scale(&self, i: usize) -> f64 {
        (self.target_weights[i] * self.magsac_weights[i]).sqrt() * self.inv_sigma[i]
    }

    /// Recomputes the MAGSAC++ weights from the current residuals.
    fn apply_magsac_weights(&mut self, weighting: &MagsacWeight) {
        for i in 0..self.points.len() {
            let r = (self.residuals[i] * self.inv_sigma[i]).abs();
            self.magsac_weights[i] = weighting.weight(r);
        }
    }

    /// How many points would still carry nonzero weight if `weighting` were applied now. Used to
    /// avoid stepping into a rank-deficient solve.
    fn count_if_reweighted(&self, weighting: &MagsacWeight) -> usize {
        (0..self.points.len())
            .filter(|&i| {
                let r = (self.residuals[i] * self.inv_sigma[i]).abs();
                self.target_weights[i] * weighting.weight(r) > 0.0
            })
            .count()
    }
}

impl<T: SurfaceTarget3> LeastSquaresProblem<f64, Dyn, U6> for PointsToSurface3<'_, T> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params.set_storage(*x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params.storage()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.residuals[i] * self.scale(i);
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surf_jacobian(p, c, &self.current) * self.scale(i);
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{clone_points, mean_point, transform_points};
    use crate::na::{Translation3, UnitQuaternion};
    use crate::tests::engine_blade;
    use crate::{Iso3, Mesh3, SelectOp, Selection, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn box_mesh() -> Mesh3 {
        Mesh3::create_box(10.0, 5.0, 2.0, false)
    }

    fn box_points(mesh: &Mesh3) -> Vec<Point3> {
        clone_points(&mesh.sample_poisson(0.5, None))
    }

    fn small_disturbance() -> Iso3 {
        Iso3::from_parts(
            Translation3::new(0.3, -0.2, 0.1),
            UnitQuaternion::from_euler_angles(PI / 60.0, -PI / 90.0, PI / 45.0),
        )
    }

    // ============================================================================================
    // Recovery: the behavior that existed before the robust path was introduced must survive it.
    // ============================================================================================

    #[test]
    fn simple_box_disturbed() -> Result<()> {
        // This test is to verify that a simple test against a box that doesn't have large rotations
        // produces a result that is roughly the inverse of the disturbance
        let mesh = box_mesh();
        let points = clone_points(&mesh.sample_poisson(0.1, None));
        let disturb = Iso3::from_parts(
            Translation3::new(3.0, 2.0, 1.0),
            UnitQuaternion::from_euler_angles(PI / 8.0, PI / 12.0, PI / 16.0),
        );

        let params = AlignParams3::from_origin(None);
        let to_align = transform_points(&points, &disturb);
        let outcome = points_to_surface3(&to_align, &mesh, params, &AlignOptions3::default())?;

        assert_relative_eq!(
            disturb.inverse(),
            outcome.alignment().full_transform(),
            epsilon = 1e-8
        );
        Ok(())
    }

    #[test]
    fn blade_example() -> Result<()> {
        let mesh = engine_blade();
        let mask = mesh
            .face_select(Selection::None)
            .facing(&Vector3::y(), PI / 4.0, SelectOp::Add)
            .take_mask();
        let expected_points = clone_points(&mesh.sample_poisson(2.0, Some(&mask)));

        let disturb = Iso3::from_parts(
            Translation3::new(-100.0, 150.0, 0.0),
            UnitQuaternion::new(Vector3::new(1.0, 1.0, 1.0).normalize() * PI / 6.0),
        );

        let to_align = transform_points(&expected_points, &disturb);

        let params = AlignParams3::from_center(mean_point(&to_align), None);
        let outcome = points_to_surface3(&to_align, &mesh, params, &AlignOptions3::default())?;

        let aligned = transform_points(&to_align, outcome.alignment().full_transform());

        let max_deviation = aligned
            .iter()
            .zip(expected_points.iter())
            .map(|(a, e)| (a - e).norm())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        assert!(
            max_deviation < 1e-6,
            "Max deviation is too high: {}",
            max_deviation
        );

        Ok(())
    }

    #[test]
    fn locked_tx_is_not_recovered() -> Result<()> {
        use crate::geom3::align3::Dof6;

        let mesh = box_mesh();
        let points = box_points(&mesh);

        // A pure x-translation, which the alignment is forbidden from undoing.
        let disturb = Iso3::translation(0.4, 0.0, 0.0);
        let to_align = transform_points(&points, &disturb);

        let dof = Dof6::new(false, true, true, true, true, true);
        let params = AlignParams3::from_origin(Some(dof));
        let outcome = points_to_surface3(&to_align, &mesh, params, &AlignOptions3::default())?;

        // With `local` and `offset` both identity, the full transform is exactly the alignment
        // transform, so a locked tx must leave the x translation at precisely zero.
        assert_eq!(
            outcome.alignment().full_transform().translation.vector.x,
            0.0
        );

        Ok(())
    }

    // ============================================================================================
    // Robust weighting and uncertainty
    // ============================================================================================

    #[test]
    fn gross_outliers_are_rejected() -> Result<()> {
        let mesh = box_mesh();
        let expected = box_points(&mesh);
        let disturb = small_disturbance();
        let mut to_align = transform_points(&expected, &disturb);

        // Throw roughly 10% of the points well away from the target, spread through the set so
        // the survivors still constrain all six degrees of freedom.
        let n = to_align.len();
        for k in 0..n / 10 {
            to_align[k * 10] += Vector3::new(0.0, 0.0, 3.0);
        }

        let params = AlignParams3::from_center(mean_point(&to_align), None);

        let naive = points_to_surface3(
            &to_align,
            &mesh,
            params.clone(),
            &AlignOptions3::no_refine(),
        )?;
        let robust = points_to_surface3(&to_align, &mesh, params, &AlignOptions3::default())?;

        // Measure against the uncorrupted points only; the outliers are not supposed to land
        // anywhere in particular.
        let clean: Vec<usize> = (0..n).filter(|i| !i.is_multiple_of(10)).collect();
        let dev_of = |a: &AlignOutcome3| {
            let t = transform_points(&to_align, a.alignment().full_transform());
            clean
                .iter()
                .map(|&i| (t[i] - expected[i]).norm())
                .fold(0.0_f64, f64::max)
        };

        let naive_dev = dev_of(&naive);
        let robust_dev = dev_of(&robust);

        assert!(
            robust_dev < naive_dev * 0.25,
            "the robust fit should have rejected the outliers: \
             naive deviation {naive_dev}, robust deviation {robust_dev}"
        );

        Ok(())
    }

    #[test]
    fn high_sigma_point_has_less_influence() -> Result<()> {
        // Isolate the per-point uncertainty mechanism from MAGSAC by disabling refinement, so the
        // only thing distinguishing the two alignments below is the residual normalization.
        let mesh = box_mesh();
        let expected = box_points(&mesh);
        let mut to_align = expected.clone();

        // Displace a single point off the target. With uniform uncertainty it drags the fit; when
        // it is declared far noisier than its neighbors it should barely register.
        let bad = to_align.len() / 2;
        to_align[bad] += Vector3::new(0.0, 0.0, 1.0);

        let params = AlignParams3::from_center(mean_point(&to_align), None);

        let uniform = points_to_surface3(
            &to_align,
            &mesh,
            params.clone(),
            &AlignOptions3::no_refine(),
        )?;

        let mut sigma = vec![0.01; to_align.len()];
        sigma[bad] = 100.0;
        let weighted = points_to_surface3(
            &to_align,
            &mesh,
            params,
            &AlignOptions3 {
                refinement_steps: 0,
                point_sigma: Some(&sigma),
                ..Default::default()
            },
        )?;

        // The remaining points already sit exactly on the target, so a fit that correctly ignores
        // the noisy point should leave them alone.
        let clean: Vec<usize> = (0..expected.len()).filter(|i| *i != bad).collect();
        let dev_of = |a: &AlignOutcome3| {
            let t = transform_points(&to_align, a.alignment().full_transform());
            clean
                .iter()
                .map(|&i| (t[i] - expected[i]).norm())
                .fold(0.0_f64, f64::max)
        };

        let uniform_dev = dev_of(&uniform);
        let weighted_dev = dev_of(&weighted);

        assert!(
            weighted_dev < uniform_dev * 0.5,
            "declaring the point noisy should have cut its influence: \
             uniform deviation {uniform_dev}, weighted deviation {weighted_dev}"
        );

        Ok(())
    }

    /// Wraps a target so that every match it produces reports a fixed measurement uncertainty.
    struct WithSigma<'a>(&'a Mesh3, f64);

    impl SurfaceTarget3 for WithSigma<'_> {
        fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
            self.0.find_align_match(p).with_sigma(self.1)
        }
    }

    #[test]
    fn target_sigma_combines_in_quadrature() -> Result<()> {
        // Two runs that should be mathematically indistinguishable: splitting the uncertainty
        // across the test points and the target, versus putting the quadrature sum entirely on
        // the test points. If the combination rule is wrong in any way other than a global scale,
        // these diverge.
        //
        // The per-point uncertainty deliberately varies, because a uniform sigma is only a global
        // scale on the objective and would leave the minimizer unchanged.
        let mesh = box_mesh();
        let mut to_align = box_points(&mesh);
        let bad = to_align.len() / 2;
        to_align[bad] += Vector3::new(0.0, 0.0, 1.0);
        to_align[1] += Vector3::new(0.0, 0.5, 0.0);

        let target_sigma = 0.04;
        let split: Vec<f64> = (0..to_align.len())
            .map(|i| if i == bad { 0.5 } else { 0.01 })
            .collect();
        let combined: Vec<f64> = split
            .iter()
            .map(|s| (s * s + target_sigma * target_sigma).sqrt())
            .collect();

        let params = AlignParams3::from_center(mean_point(&to_align), None);

        fn opts(sigma: &[f64]) -> AlignOptions3<'_> {
            AlignOptions3 {
                refinement_steps: 0,
                point_sigma: Some(sigma),
                ..Default::default()
            }
        }

        // Uncertainty split between the test points and the target...
        let a = points_to_surface3(
            &to_align,
            &WithSigma(&mesh, target_sigma),
            params.clone(),
            &opts(&split),
        )?;

        // ...versus the same total placed entirely on the test points.
        let b = points_to_surface3(
            &to_align,
            &WithSigma(&mesh, 0.0),
            params.clone(),
            &opts(&combined),
        )?;

        assert_relative_eq!(
            a.alignment().full_transform().to_matrix(),
            b.alignment().full_transform().to_matrix(),
            epsilon = 1e-9
        );

        // Guard against the comparison passing because the target sigma did nothing at all.
        let c = points_to_surface3(&to_align, &WithSigma(&mesh, 0.0), params, &opts(&split))?;
        let delta = (c.alignment().full_transform().to_matrix()
            - a.alignment().full_transform().to_matrix())
        .norm();
        assert!(
            delta > 1e-6,
            "target sigma had no effect on the alignment (delta {delta})"
        );

        Ok(())
    }

    #[test]
    fn mesh_reports_its_own_point_stdev() -> Result<()> {
        // A mesh carrying per-vertex uncertainty must surface it on the match, interpolated to
        // wherever the projection landed rather than snapped to a vertex.
        let mut mesh = box_mesh();
        mesh.set_point_stdev(Some(vec![0.25; mesh.point_count()]))?;

        let m = mesh.find_align_match(&Point3::new(1.0, 1.0, 5.0));
        assert_relative_eq!(m.sigma, 0.25, epsilon = 1e-12);

        // ...and a mesh without the attribute reports no uncertainty at all.
        let plain = box_mesh();
        assert_eq!(
            plain.find_align_match(&Point3::new(1.0, 1.0, 5.0)).sigma,
            0.0
        );

        Ok(())
    }

    #[test]
    fn reported_residuals_are_geometric() -> Result<()> {
        // The residuals on the result are real distances in model units, not the internally
        // weighted and normalized values the solver minimizes.
        let mesh = box_mesh();
        let mut to_align = box_points(&mesh);
        let bad = to_align.len() / 2;
        to_align[bad] += Vector3::new(0.0, 0.0, 0.4);

        let sigma = vec![0.01; to_align.len()];
        let opts = AlignOptions3 {
            refinement_steps: 0,
            point_sigma: Some(&sigma),
            ..Default::default()
        };

        let params = AlignParams3::from_center(mean_point(&to_align), None);
        let outcome = points_to_surface3(&to_align, &mesh, params, &opts)?;

        // With sigma = 0.01 throughout, a normalized residual would be 100x the geometric one, so
        // this assertion fails loudly if the wrong vector is reported.
        let largest = outcome
            .alignment()
            .residuals()
            .iter()
            .map(|r| r.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            largest < 1.0,
            "residuals look normalized rather than geometric, largest was {largest}"
        );

        Ok(())
    }

    // ============================================================================================
    // Reporting: what survives as a result, and what is genuinely an error
    // ============================================================================================

    #[test]
    fn exhausted_budget_is_reported_not_raised() -> Result<()> {
        // A patience of 1 allows only `1 * (6 + 1)` function evaluations. Paired with a
        // disturbance large enough to need more than that, the solver runs out of budget. It still
        // leaves behind the best parameters it found, so this must come back as a usable (if
        // unproven) alignment rather than an error.
        //
        // The disturbance has to be large: `small_disturbance` converges by `xtol` within the
        // budget, which is a genuine convergence and correctly reported as such.
        let mesh = box_mesh();
        let disturb = Iso3::from_parts(
            Translation3::new(3.0, 2.0, 1.0),
            UnitQuaternion::from_euler_angles(PI / 8.0, PI / 12.0, PI / 16.0),
        );
        let to_align = transform_points(&box_points(&mesh), &disturb);

        let opts = AlignOptions3 {
            patience: 1,
            ..Default::default()
        };

        let outcome = points_to_surface3(&to_align, &mesh, AlignParams3::from_origin(None), &opts)?;

        assert_eq!(outcome.quality(), SolveQuality::Unconverged);
        assert!(!outcome.converged());
        assert!(
            outcome.solves().contains(&TerminationReason::LostPatience),
            "expected a LostPatience among {:?}",
            outcome.solves()
        );

        let t = outcome.alignment().full_transform().to_matrix();
        assert!(t.iter().all(|v| v.is_finite()));

        Ok(())
    }

    #[test]
    fn converged_alignment_reports_converged() -> Result<()> {
        let mesh = box_mesh();
        let to_align = transform_points(&box_points(&mesh), &small_disturbance());

        let opts = AlignOptions3 {
            sigma_max: Some(0.5),
            refinement_steps: 4,
            ..Default::default()
        };

        let outcome = points_to_surface3(&to_align, &mesh, AlignParams3::from_origin(None), &opts)?;

        assert!(outcome.converged());
        assert_eq!(outcome.quality(), SolveQuality::Converged);
        assert_eq!(outcome.refinement_rounds(), 4);
        assert_eq!(outcome.solves().len(), 5);
        assert_eq!(outcome.halt(), None);

        Ok(())
    }

    #[test]
    fn degenerate_noise_estimate_halts_refinement() -> Result<()> {
        // When every residual is identical the median absolute deviation collapses and there is no
        // noise scale to estimate from. Refinement is skipped, which is a normal outcome and must
        // be visible rather than silent.
        //
        // This uses a target that reports the query point back as its own match, because points
        // sampled from a real mesh do not quite reach a hard zero: projecting them back onto the
        // mesh leaves residuals at the 1e-15 level, whose spread is small but nonzero, and
        // refinement runs against that. `estimate_sigma_max` is unit-tested separately for the
        // degenerate case; this test covers the plumbing that turns its `None` into a halt.
        struct ExactTarget;

        impl SurfaceTarget3 for ExactTarget {
            fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
                AlignSurfMatch3::new(*p, Vector3::z_axis(), true, 1.0)
            }
        }

        let mesh = box_mesh();
        let to_align = box_points(&mesh);

        let outcome = points_to_surface3(
            &to_align,
            &ExactTarget,
            AlignParams3::from_origin(None),
            &AlignOptions3::default(),
        )?;

        assert_eq!(outcome.halt(), Some(&RefinementHalt::NoNoiseEstimate));
        assert_eq!(outcome.refinement_rounds(), 0);

        Ok(())
    }

    #[test]
    fn underdetermined_reweighting_halts_refinement() -> Result<()> {
        // Scaling the points about their centroid leaves a deviation no rigid transform can
        // remove, so every residual stays large. Paired with a tiny noise bound, MAGSAC drives
        // every weight to zero and the next solve would be rank-deficient.
        let mesh = box_mesh();
        let expected = box_points(&mesh);
        let center = mean_point(&expected);
        let to_align: Vec<Point3> = expected
            .iter()
            .map(|p| center + (p - center) * 1.5)
            .collect();

        let opts = AlignOptions3 {
            sigma_max: Some(1e-6),
            ..Default::default()
        };

        let outcome = points_to_surface3(
            &to_align,
            &mesh,
            AlignParams3::from_center(center, None),
            &opts,
        )?;

        assert_eq!(
            outcome.halt(),
            Some(&RefinementHalt::Underdetermined {
                weighted: 0,
                params: N_PARAMS
            })
        );
        assert_eq!(outcome.refinement_rounds(), 0);
        assert!(
            outcome
                .alignment()
                .residuals()
                .iter()
                .all(|r| r.is_finite())
        );

        Ok(())
    }

    #[test]
    fn failed_initial_solve_is_an_error() {
        // With no points there are no residuals, so the solver has nothing to work with and there
        // is no earlier result to fall back to. This is the one case that must still be an `Err`.
        let mesh = box_mesh();
        let result = points_to_surface3(
            &[],
            &mesh,
            AlignParams3::from_origin(None),
            &AlignOptions3::default(),
        );

        let err = result.unwrap_err().to_string();
        assert!(err.contains("Failed to align"), "unexpected message: {err}");
    }

    #[test]
    fn patience_is_validated() {
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let opts = AlignOptions3 {
            patience: 0,
            ..Default::default()
        };

        let err = points_to_surface3(&points, &mesh, AlignParams3::from_origin(None), &opts)
            .unwrap_err()
            .to_string();
        assert!(err.contains("patience"), "unexpected message: {err}");
    }

    #[test]
    fn sigma_max_is_validated() {
        let mesh = box_mesh();
        let points = box_points(&mesh);

        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let opts = AlignOptions3 {
                sigma_max: Some(bad),
                ..Default::default()
            };

            let result = points_to_surface3(&points, &mesh, AlignParams3::from_origin(None), &opts);
            assert!(result.is_err(), "sigma_max {bad} should have been rejected");
        }
    }

    #[test]
    fn point_sigma_is_validated() {
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let short = vec![0.1; points.len() - 1];
        let err = points_to_surface3(
            &points,
            &mesh,
            AlignParams3::from_origin(None),
            &AlignOptions3 {
                point_sigma: Some(&short),
                ..Default::default()
            },
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("entries"), "unexpected message: {err}");

        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut sigma = vec![0.1; points.len()];
            sigma[3] = bad;

            let result = points_to_surface3(
                &points,
                &mesh,
                AlignParams3::from_origin(None),
                &AlignOptions3 {
                    point_sigma: Some(&sigma),
                    ..Default::default()
                },
            );
            assert!(
                result.is_err(),
                "sigma value {bad} should have been rejected"
            );
        }
    }

    // ============================================================================================
    // Supporting pieces
    // ============================================================================================

    #[test]
    fn median_handles_both_parities() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), Some(2.5));
        assert_eq!(median(&[]), None);
    }

    #[test]
    fn sigma_estimate_rejects_degenerate_spread() {
        assert_eq!(estimate_sigma_max(&[2.0; 10]), None);
        assert_eq!(estimate_sigma_max(&[]), None);
    }
}
