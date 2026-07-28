use crate::Result;
use crate::common::consensus::weights::MagsacWeight;
use crate::common::dist;
use crate::geom2::Alignment2;
use crate::geom2::Point2;
use crate::geom2::align2::jacobian::{copy_jacobian, point_surf_jacobian2};
use crate::geom2::align2::{
    AlignOptions2, AlignParams2, AlignSurfMatch2, AlignValues2, SurfaceTarget2,
};
use crate::na::{Dyn, Matrix, Owned, U1, U3, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

/// The residual degrees of freedom for the MAGSAC++ weight function. The residual here is a full
/// Euclidean point-to-target distance in the plane, so it follows a chi distribution with two
/// degrees of freedom. (`MagsacWeight` requires at least 2; a point-to-plane residual would be a
/// one-dimensional projection and would need the weight function extended.)
const RESIDUAL_DOF: usize = 2;

/// The scale factor that turns a median absolute deviation into a consistent estimate of the
/// standard deviation of normally distributed data.
const MAD_TO_SIGMA: f64 = 1.4826;

/// Performs a Levenberg-Marquardt minimization to align a set of 2D points to a surface target.
///
/// The points are repeatedly projected onto their closest position on the target as the solver
/// moves them, so the correspondences are re-established at every step rather than being fixed up
/// front.
///
/// By default this is a robust alignment: an initial unweighted solve is followed by
/// [`AlignOptions2::refinement_steps`] rounds of iteratively reweighted least squares using
/// MAGSAC++ weights. See [`AlignOptions2`] for the weighting and per-point uncertainty controls.
///
/// # Arguments
///
/// * `points`: the 2D points to be aligned, in their own local coordinate system
/// * `target`: an entity that implements the [`SurfaceTarget2`] trait, such as a `Curve2` or a
///   `Boundary2`, which will be the stationary target of the points.
/// * `params`: the alignment parameters, see [`AlignParams2`] for details
/// * `opts`: options controlling weighting and filtering, see [`AlignOptions2`]
///
/// returns: Result<Alignment<Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
pub fn points_to_surface2(
    points: &[Point2],
    target: &impl SurfaceTarget2,
    params: AlignParams2,
    opts: &AlignOptions2,
) -> Result<Alignment2> {
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

    // The first solve is unweighted. Its residuals are what the sigma_max estimate is drawn from,
    // since residuals taken before any alignment would describe the initial misplacement rather
    // than the noise.
    let mut problem = solve(PointsToSurface2::new(points, target, params, opts))?;

    if opts.refinement_steps > 0
        && let Some(sigma_max) = resolve_sigma_max(opts, &problem)
    {
        let weighting = MagsacWeight::new(sigma_max, RESIDUAL_DOF);
        for _ in 0..opts.refinement_steps {
            // The matches have moved since the last refresh, so any target-side uncertainty is
            // re-read before it feeds the reweighting.
            problem.refresh_inv_sigma();

            // If reweighting would leave the problem underdetermined, keep the last good result
            // rather than solving something rank-deficient.
            if problem.count_if_reweighted(&weighting) < N_PARAMS {
                break;
            }
            problem.apply_magsac_weights(&weighting);
            problem = solve(problem)?;
        }
    }

    let c = problem.params.current_values();
    Ok(Alignment2::new(
        c.transform,
        c.align,
        problem.params.local,
        problem.params.offset,
        // The geometric residuals, not the weighted/normalized ones the solver minimizes. These
        // are what a caller wants to inspect: real deviations in the units of the geometry.
        problem.residuals.clone(),
    ))
}

/// The number of free parameters in a 2D alignment (tx, ty, rz).
const N_PARAMS: usize = 3;

fn solve<T: SurfaceTarget2>(problem: PointsToSurface2<'_, T>) -> Result<PointsToSurface2<'_, T>> {
    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    if report.termination.was_successful() {
        Ok(result)
    } else {
        Err(format!(
            "Failed to align points to surface: {:?}",
            report.termination
        )
        .into())
    }
}

/// Determines the MAGSAC++ noise bound, either from the explicit option or by estimating it from
/// the current residuals. Returns `None` when no usable estimate exists, in which case robust
/// refinement is skipped.
fn resolve_sigma_max<T: SurfaceTarget2>(
    opts: &AlignOptions2,
    problem: &PointsToSurface2<'_, T>,
) -> Option<f64> {
    match opts.sigma_max {
        Some(s) if s.is_finite() && s > 0.0 => Some(s),
        Some(_) => None,
        None => {
            let normalized = problem.normalized_residuals();
            estimate_sigma_max(&normalized)
        }
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

struct PointsToSurface2<'a, T: SurfaceTarget2> {
    points: &'a [Point2],
    target: &'a T,
    params: AlignParams2,

    /// The alignment values for the current parameters. Cached here so that `jacobian()` doesn't
    /// have to recompute what `move_points` already worked out.
    current: AlignValues2,

    /// The test points after being moved by the current transform, in the target's space.
    moved: Vec<Point2>,

    /// The closest match on the target for each moved point.
    closest: Vec<AlignSurfMatch2>,

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

    ignore_off: bool,
}

impl<'a, T: SurfaceTarget2> PointsToSurface2<'a, T> {
    fn new(
        points: &'a [Point2],
        target: &'a T,
        params: AlignParams2,
        opts: &AlignOptions2,
    ) -> Self {
        let current = params.current_values();
        let test_sigma = match opts.point_sigma {
            Some(s) => s.to_vec(),
            None => vec![0.0; points.len()],
        };

        let mut x = Self {
            points,
            target,
            params,
            current,
            moved: vec![Point2::origin(); points.len()],
            closest: vec![AlignSurfMatch2::default(); points.len()],
            residuals: vec![0.0; points.len()],
            test_sigma,
            inv_sigma: vec![1.0; points.len()],
            target_weights: vec![1.0; points.len()],
            magsac_weights: vec![1.0; points.len()],
            ignore_off: opts.ignore_off,
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
        self.current = self.params.current_values();
        let transform = self.current.transform;

        for i in 0..self.points.len() {
            let m = transform * self.points[i];
            let c = self.target.align_surf_closest_to(&m);

            // The residual is the distance between the test point and the closest point on the
            // target, adjusted for the direction of the scalar projection.
            self.residuals[i] = dist(&m, &c.point) * c.dn(&m).signum();

            self.target_weights[i] = if self.ignore_off {
                c.weight * f64::from(c.is_on)
            } else {
                c.weight
            };

            self.moved[i] = m;
            self.closest[i] = c;
        }
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

impl<T: SurfaceTarget2> LeastSquaresProblem<f64, Dyn, U3> for PointsToSurface2<'_, T> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U3>;
    type ParameterStorage = Owned<f64, U3>;

    fn set_params(&mut self, x: &Vector<f64, U3, Self::ParameterStorage>) {
        self.params.set_storage(*x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        self.params.get_storage()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.residuals[i] * self.scale(i);
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U3, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U3, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surf_jacobian2(p, c, &self.current) * self.scale(i);
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{mean_point, transform_points};
    use crate::geom2::align2::Dof3;
    use crate::geom2::{Boundary2, BoundaryData2, BoundaryEditor};
    use crate::{Curve2, Iso2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn to_pts(v: &[(f64, f64)]) -> Vec<Point2> {
        v.iter().map(|(x, y)| Point2::new(*x, *y)).collect()
    }

    fn closed_curve(p: &[(f64, f64)]) -> Curve2 {
        Curve2::from_points(&to_pts(p), 1e-8, true).unwrap()
    }

    fn rect_curve() -> Curve2 {
        closed_curve(&[(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (0.0, 1.0)])
    }

    fn rect_points() -> Vec<Point2> {
        to_pts(&[
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 0.0),
            (5.0, 0.25),
            (5.0, 0.75),
            (1.0, 1.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (0.0, 0.25),
            (0.0, 0.75),
        ])
    }

    /// A closed CCW "stadium" boundary: two horizontal segments joined by semicircular arcs at
    /// each end. Exercises a `Boundary2` target containing genuine `Arc2` elements rather than
    /// only segments.
    fn stadium_boundary() -> Boundary2 {
        let mut data = BoundaryData2::new_closed();
        let mut cursor = data.get_cursor(None);
        // Bottom edge, left to right
        cursor.add_seg_xy(2.0, 0.0);
        // Right cap, counter-clockwise around (2.0, 0.5)
        cursor.add_arc_xy(2.0, 0.5, 2.0, 1.0, false);
        // Top edge, right to left
        cursor.add_seg_xy(0.0, 1.0);
        // Left cap, counter-clockwise around (0.0, 0.5)
        cursor.add_arc_xy(0.0, 0.5, 0.0, 0.0, false);
        data.try_to_boundary().unwrap()
    }

    fn small_disturbance() -> Iso2 {
        Iso2::new(crate::Vector2::new(0.05, -0.05), 10.0 * PI / 180.0)
    }

    /// The largest distance between the aligned points and where they should have landed.
    fn max_deviation(moved: &[Point2], expected: &[Point2], result: &Alignment2) -> f64 {
        transform_points(moved, result.full())
            .iter()
            .zip(expected.iter())
            .map(|(a, e)| (a - e).norm())
            .fold(0.0_f64, f64::max)
    }

    // ============================================================================================
    // The unweighted behavior from step 4 must survive the introduction of the robust path.
    // ============================================================================================

    #[test]
    fn curve_recovers_disturbance() -> Result<()> {
        let curve = rect_curve();
        let disturb = small_disturbance();
        let moved = transform_points(&rect_points(), &disturb);

        let params = AlignParams2::new_at_origin(None);
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-10
        );
        Ok(())
    }

    #[test]
    fn curve_recovers_disturbance_at_center() -> Result<()> {
        // Same as above, but rotating about the centroid of the moved points rather than the
        // world origin. The recovered full transform must be identical either way.
        let curve = rect_curve();
        let disturb = small_disturbance();
        let moved = transform_points(&rect_points(), &disturb);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-10
        );
        Ok(())
    }

    #[test]
    fn boundary_with_arcs_recovers_disturbance() -> Result<()> {
        let boundary = stadium_boundary();
        let points = boundary.to_points(1e-4)?;

        let disturb = Iso2::new(crate::Vector2::new(0.03, -0.02), 5.0 * PI / 180.0);
        let moved = transform_points(&points, &disturb);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);
        let result = points_to_surface2(&moved, &boundary, params, &AlignOptions2::default())?;

        // The sampled points sit on the theoretical boundary, so the alignment should recover the
        // disturbance almost exactly. The tolerance is looser than the curve case because the
        // points are a chordal approximation of the arcs, not exact positions on them.
        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-6
        );
        Ok(())
    }

    #[test]
    fn locked_tx_is_not_recovered() -> Result<()> {
        let curve = rect_curve();
        let points = rect_points();

        // A pure x-translation, which the alignment is forbidden from undoing.
        let disturb = Iso2::translation(0.3, 0.0);
        let moved = transform_points(&points, &disturb);

        let dof = Dof3::new(false, true, true);
        let params = AlignParams2::new_at_origin(Some(dof));
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        // With `local` and `offset` both identity, the full transform is exactly the alignment
        // transform, so a locked tx must leave the x translation at precisely zero.
        assert_eq!(result.full().translation.vector.x, 0.0);

        // ...and the disturbance genuinely was not recovered.
        let max_dev = max_deviation(&moved, &points, &result);
        assert!(
            max_dev > 0.1,
            "locked tx should have prevented recovery, max deviation was {max_dev}"
        );

        Ok(())
    }

    // ============================================================================================
    // Robust weighting
    // ============================================================================================

    #[test]
    fn clean_data_recovery_survives_irls() -> Result<()> {
        // On noiseless data the automatic sigma_max estimate is degenerate (the MAD collapses to
        // zero once the fit is exact) and refinement is skipped, so an explicit sigma_max is used
        // here to force every reweighting round to actually run. Robust weighting must not
        // degrade a case that plain least squares already solves exactly.
        let curve = rect_curve();
        let disturb = small_disturbance();
        let moved = transform_points(&rect_points(), &disturb);

        let opts = AlignOptions2 {
            sigma_max: Some(0.5),
            refinement_steps: 4,
            ..Default::default()
        };

        let params = AlignParams2::new_at_origin(None);
        let result = points_to_surface2(&moved, &curve, params, &opts)?;

        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-10
        );
        Ok(())
    }

    #[test]
    fn gross_outliers_are_rejected() -> Result<()> {
        let curve = rect_curve();
        let expected = rect_points();
        let disturb = small_disturbance();
        let mut moved = transform_points(&expected, &disturb);

        // Corrupt 20% of the points by throwing them well away from the target. Their indices are
        // spread out so the surviving points still constrain all three degrees of freedom.
        moved[2] += crate::Vector2::new(0.0, -1.5);
        moved[7] += crate::Vector2::new(0.0, 2.0);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);

        // Plain least squares is dragged off by the outliers.
        let naive = points_to_surface2(
            &moved,
            &curve,
            params.clone(),
            &AlignOptions2 {
                refinement_steps: 0,
                ..Default::default()
            },
        )?;

        // The robust path, with sigma_max estimated automatically from the initial residuals.
        let robust = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        // Measure against the eight uncorrupted points only; the two outliers are not supposed to
        // land anywhere in particular.
        let clean: Vec<usize> = (0..expected.len()).filter(|i| *i != 2 && *i != 7).collect();
        let dev_of = |a: &Alignment2| {
            let t = transform_points(&moved, a.full());
            clean
                .iter()
                .map(|&i| (t[i] - expected[i]).norm())
                .fold(0.0_f64, f64::max)
        };

        let naive_dev = dev_of(&naive);
        let robust_dev = dev_of(&robust);

        assert!(
            naive_dev > 0.05,
            "the outliers should have pulled the unweighted fit off, but deviation was {naive_dev}"
        );
        assert!(
            robust_dev < 1e-6,
            "the robust fit should have rejected the outliers, but deviation was {robust_dev}"
        );

        Ok(())
    }

    #[test]
    fn high_sigma_point_has_less_influence() -> Result<()> {
        // Isolate the per-point uncertainty mechanism from MAGSAC by disabling refinement, so the
        // only thing distinguishing the two alignments below is the residual normalization.
        let curve = rect_curve();
        let expected = rect_points();
        let mut moved = expected.clone();

        // Displace a single point off the target. With uniform uncertainty it drags the fit; when
        // it is declared far noisier than its neighbors it should barely register.
        let bad = 4;
        moved[bad] += crate::Vector2::new(0.4, 0.0);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);

        let uniform = points_to_surface2(
            &moved,
            &curve,
            params.clone(),
            &AlignOptions2 {
                refinement_steps: 0,
                ..Default::default()
            },
        )?;

        let mut sigma = vec![0.01; moved.len()];
        sigma[bad] = 10.0;
        let weighted = points_to_surface2(
            &moved,
            &curve,
            params,
            &AlignOptions2 {
                refinement_steps: 0,
                point_sigma: Some(&sigma),
                ..Default::default()
            },
        )?;

        // The remaining points already sit exactly on the target, so a fit that correctly ignores
        // the noisy point should leave them alone.
        let clean: Vec<usize> = (0..expected.len()).filter(|i| *i != bad).collect();
        let dev_of = |a: &Alignment2| {
            let t = transform_points(&moved, a.full());
            clean
                .iter()
                .map(|&i| (t[i] - expected[i]).norm())
                .fold(0.0_f64, f64::max)
        };

        let uniform_dev = dev_of(&uniform);
        let weighted_dev = dev_of(&weighted);

        assert!(
            weighted_dev < uniform_dev * 0.1,
            "declaring the point noisy should have cut its influence sharply: \
             uniform deviation {uniform_dev}, weighted deviation {weighted_dev}"
        );

        Ok(())
    }

    /// Wraps a target so that every match it produces reports a fixed measurement uncertainty.
    /// Neither `Curve2` nor `Boundary2` carries per-vertex uncertainty today, so this stands in
    /// for a target built from measured data.
    struct WithSigma<'a>(&'a Curve2, f64);

    impl SurfaceTarget2 for WithSigma<'_> {
        fn align_surf_closest_to(&self, p: &Point2) -> AlignSurfMatch2 {
            self.0.align_surf_closest_to(p).with_sigma(self.1)
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
        // scale on the objective and would leave the minimizer unchanged, making the comparison
        // vacuous.
        let curve = rect_curve();
        let expected = rect_points();
        let mut moved = expected.clone();
        moved[4] += crate::Vector2::new(0.4, 0.0);
        moved[1] += crate::Vector2::new(0.0, -0.2);

        let target_sigma = 0.04;
        let split: Vec<f64> = (0..moved.len())
            .map(|i| if i == 4 { 0.5 } else { 0.01 })
            .collect();
        let combined: Vec<f64> = split
            .iter()
            .map(|s| (s * s + target_sigma * target_sigma).sqrt())
            .collect();

        let params = AlignParams2::new_at_center(mean_point(&moved), None);

        fn opts(sigma: &[f64]) -> AlignOptions2<'_> {
            AlignOptions2 {
                refinement_steps: 0,
                point_sigma: Some(sigma),
                ..Default::default()
            }
        }

        // Uncertainty split between the test points and the target...
        let a = points_to_surface2(
            &moved,
            &WithSigma(&curve, target_sigma),
            params.clone(),
            &opts(&split),
        )?;

        // ...versus the same total placed entirely on the test points.
        let b = points_to_surface2(
            &moved,
            &WithSigma(&curve, 0.0),
            params.clone(),
            &opts(&combined),
        )?;

        assert_relative_eq!(a.full().to_matrix(), b.full().to_matrix(), epsilon = 1e-9);

        // Guard against the comparison passing because the target sigma did nothing at all: with
        // only the test-side sigma, the result must differ from both of the above.
        let c = points_to_surface2(&moved, &WithSigma(&curve, 0.0), params, &opts(&split))?;
        let delta = (c.full().to_matrix() - a.full().to_matrix()).norm();
        assert!(
            delta > 1e-6,
            "target sigma had no effect on the alignment (delta {delta})"
        );

        Ok(())
    }

    #[test]
    fn target_sigma_alone_is_enough() -> Result<()> {
        // Target uncertainty must work without any test-side `point_sigma` at all, and must not
        // divide by zero when a match reports no uncertainty.
        let curve = rect_curve();
        let moved = rect_points();
        let params = AlignParams2::new_at_center(mean_point(&moved), None);

        let with = points_to_surface2(
            &moved,
            &WithSigma(&curve, 0.05),
            params.clone(),
            &AlignOptions2::default(),
        )?;
        let without = points_to_surface2(
            &moved,
            &WithSigma(&curve, 0.0),
            params,
            &AlignOptions2::default(),
        )?;

        for a in [&with, &without] {
            for r in a.residuals() {
                assert!(r.is_finite(), "residual was not finite: {r}");
            }
        }

        // A uniform sigma is a global scale on the objective, so the minimizer is unchanged.
        assert_relative_eq!(
            with.full().to_matrix(),
            without.full().to_matrix(),
            epsilon = 1e-9
        );

        Ok(())
    }

    #[test]
    fn reported_residuals_are_geometric() -> Result<()> {
        // The residuals on the result are real distances in model units, not the internally
        // weighted and normalized values the solver minimizes.
        let curve = rect_curve();
        let expected = rect_points();
        let mut moved = expected.clone();
        moved[4] += crate::Vector2::new(0.4, 0.0);

        let sigma = vec![0.01; moved.len()];
        let opts = AlignOptions2 {
            refinement_steps: 0,
            point_sigma: Some(&sigma),
            ..Default::default()
        };

        let params = AlignParams2::new_at_center(mean_point(&moved), None);
        let result = points_to_surface2(&moved, &curve, params, &opts)?;

        // With sigma = 0.01 throughout, a normalized residual would be 100x the geometric one, so
        // this assertion fails loudly if the wrong vector is reported.
        let largest = result
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

    #[test]
    fn point_sigma_length_is_validated() {
        let curve = rect_curve();
        let points = rect_points();
        let sigma = vec![0.1; points.len() - 1];

        let opts = AlignOptions2 {
            point_sigma: Some(&sigma),
            ..Default::default()
        };

        let err = points_to_surface2(&points, &curve, AlignParams2::new_at_origin(None), &opts)
            .unwrap_err()
            .to_string();
        assert!(err.contains("entries"), "unexpected message: {err}");
    }

    #[test]
    fn point_sigma_values_are_validated() {
        let curve = rect_curve();
        let points = rect_points();

        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut sigma = vec![0.1; points.len()];
            sigma[3] = bad;

            let opts = AlignOptions2 {
                point_sigma: Some(&sigma),
                ..Default::default()
            };

            let result =
                points_to_surface2(&points, &curve, AlignParams2::new_at_origin(None), &opts);
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
    fn sigma_estimate_is_robust_to_outliers() {
        // Eleven values with unit spread, plus one enormous outlier. A standard deviation would
        // be dominated by the outlier; the MAD-based estimate should barely notice it.
        let mut values: Vec<f64> = (-5..=5).map(|i| i as f64).collect();
        let clean = estimate_sigma_max(&values).unwrap();

        values.push(10_000.0);
        let contaminated = estimate_sigma_max(&values).unwrap();

        assert!(
            (contaminated - clean).abs() < 0.5 * clean,
            "outlier moved the estimate too far: {clean} -> {contaminated}"
        );
    }

    #[test]
    fn sigma_estimate_rejects_degenerate_spread() {
        // Every residual identical means no spread to estimate from.
        assert_eq!(estimate_sigma_max(&[2.0; 10]), None);
        assert_eq!(estimate_sigma_max(&[]), None);
    }
}
