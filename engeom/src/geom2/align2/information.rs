//! Tools for asking how well a set of points constrains a 2D alignment, and for choosing a
//! subset of them that preserves that conditioning.
//!
//! This is the 2D counterpart of `geom3::align3::information`, and is deliberately structured
//! identically to it.
//!
//! # What this is for
//!
//! A points-to-surface alignment is a least-squares problem, and every test point contributes one
//! row `j_i` to its jacobian. The **information matrix** is the sum of the rank-1 contributions of
//! those rows,
//!
//! $$ H = \sum_i w_i j_i j_i^T $$
//!
//! and it is the object that answers both questions this module exists for. Its inverse is (up to
//! a noise scale) the covariance of the fitted parameters, so it says how well each degree of
//! freedom is pinned down; and because it is a *sum* over points, the contribution of any single
//! point is a rank-1 update that can be added or removed in closed form.
//!
//! Two distinct uses follow:
//!
//! - **Diagnosis.** Is this alignment well conditioned, and if not, along what motion is it free
//!   to slide? See [`AlignInformation2::marginal_precision`] and
//!   [`AlignInformation2::weak_directions`].
//! - **Pruning.** Which points are actually carrying the alignment, and what is the smallest
//!   subset that would constrain it just as well? See [`AlignInformation2::leverage`] and
//!   [`AlignInformation2::select_d_optimal`].
//!
//! # Why scoring points individually is not enough
//!
//! The obvious approach is to score each point on its own and keep the best `n`. That
//! double-counts redundancy: two points sitting on the same flat stretch each score well, but once
//! either is in the set the other adds almost nothing. Greedy D-optimal selection avoids this for
//! free, because a candidate is scored by its *marginal* gain against everything already chosen,
//! and that gain collapses as soon as something equivalent has been picked.
//!
//! # A caveat about units
//!
//! The translation columns of `H` are dimensionless (a dot product of two unit vectors) while the
//! rotation column carries units of length, because a rotation partial scales with the distance
//! from the local origin. `H` is therefore heterogeneous, and two consequences follow.
//!
//! [`AlignInformation2::marginal_precision`] is safe: each entry concerns a single parameter, so
//! its units are internally consistent. [`AlignInformation2::weak_directions`] mixes translation
//! and rotation in one vector and is only meaningful when the two are comparably scaled, which in
//! practice means the local origin should sit near the middle of the point set (as
//! `AlignParams2::from_center` arranges). Placing the origin far away inflates the rotation term
//! and will make rotation look artificially well constrained.

use crate::geom2::align2::jacobian::point_surf_jacobian2;
use crate::geom2::align2::{AlignParams2, AlignStorage2, Dof3, SurfaceTarget2};
use crate::na::{DMatrix, DVector, Matrix3};
use crate::{Point2, Result};

/// The number of alignment parameters (tx, ty, rz).
const N_PARAMS: usize = 3;

/// The relative size of the regularizing ridge used to start greedy selection. See
/// [`AlignInformation2::select_d_optimal`] for why one is needed.
const DEFAULT_RIDGE_REL: f64 = 1e-9;

/// The information content of a set of test points with respect to a 2D alignment.
///
/// Built with [`AlignInformation2::from_points`], which projects the points onto a target and
/// takes the jacobian row of each correspondence, or with [`AlignInformation2::from_rows`] when
/// the rows have already been computed. See the module documentation for what it is good for.
#[derive(Clone, Debug)]
pub struct AlignInformation2 {
    /// The jacobian row for each point. Locked degrees of freedom are already zero in these,
    /// because `point_surf_jacobian2` zeroes their columns.
    rows: Vec<AlignStorage2>,

    /// The weight of each point's contribution to the information matrix.
    weights: Vec<f64>,

    /// Which degrees of freedom are free to move.
    dof: Dof3,

    /// The indices, within the three parameters, of the degrees of freedom which are free. All of
    /// the linear algebra happens in this subspace, because `H` is identically zero along a locked
    /// direction and would not be invertible in the full space.
    free: Vec<usize>,

    /// The information matrix in the free subspace, `k` by `k` where `k = free.len()`.
    h: DMatrix<f64>,
}

impl AlignInformation2 {
    /// Builds the information content of a set of points against a target, in the pose described
    /// by `params`.
    ///
    /// Each point is moved by the current transform, projected onto the target, and turned into a
    /// jacobian row exactly as [`crate::geom2::align2::points_to_surface2`] would do. A point is
    /// weighted by the target's own `weight` for the correspondence, divided by the variance the
    /// target reports there, so that a match on a noisier part of the target carries
    /// proportionally less information.
    ///
    /// This deliberately ignores robust weighting. The question being asked is a geometric one
    /// ("how well does this arrangement of points pin down the pose"), and the answer should not
    /// depend on which points happen to be outliers in some particular fit.
    ///
    /// # Arguments
    ///
    /// * `points`: the test points, in their own local coordinate system
    /// * `target`: the stationary entity the points would be aligned to
    /// * `params`: the alignment parameterization, which fixes the local origin the rotation
    ///   partial is taken about and which degrees of freedom are free
    pub fn from_points(
        points: &[Point2],
        target: &impl SurfaceTarget2,
        params: &AlignParams2,
    ) -> Self {
        let current = params.compute_values();
        let mut rows = Vec::with_capacity(points.len());
        let mut weights = Vec::with_capacity(points.len());

        for p in points {
            let m = current.transform * p;
            let c = target.find_align_match(&m);

            rows.push(point_surf_jacobian2(&m, &c, &current));
            weights.push(if c.sigma > 0.0 {
                c.weight / (c.sigma * c.sigma)
            } else {
                c.weight
            });
        }

        Self::assemble(rows, weights, params.dof)
    }

    /// Builds the information content directly from jacobian rows which have already been
    /// computed, for callers assembling a problem this module doesn't know how to build (a
    /// multi-body adjustment, for instance).
    ///
    /// Rows are expected to already have their locked columns zeroed, which is what
    /// `point_surf_jacobian2` produces when given a constrained [`Dof3`].
    ///
    /// Returns an error if the two slices disagree in length, or if any weight is negative or
    /// non-finite.
    pub fn from_rows(rows: Vec<AlignStorage2>, weights: Vec<f64>, dof: Dof3) -> Result<Self> {
        if rows.len() != weights.len() {
            return Err(format!(
                "there are {} jacobian rows but {} weights",
                rows.len(),
                weights.len()
            )
            .into());
        }

        if let Some(i) = weights.iter().position(|w| !w.is_finite() || *w < 0.0) {
            return Err(format!(
                "weights[{}] is {}, but every weight must be finite and non-negative",
                i, weights[i]
            )
            .into());
        }

        Ok(Self::assemble(rows, weights, dof))
    }

    fn assemble(rows: Vec<AlignStorage2>, weights: Vec<f64>, dof: Dof3) -> Self {
        let free = free_indices(&dof);
        let k = free.len();

        let mut h = DMatrix::zeros(k, k);
        for (row, &w) in rows.iter().zip(weights.iter()) {
            let j = restrict(row, &free);
            // The rank-1 contribution of this point, `w * j * j^T`. This uses the general rank-1
            // update rather than the symmetric one, because `syger` writes only the lower
            // triangle and every consumer here (the inverse, the quadratic forms) needs the whole
            // matrix populated.
            h.ger(w, &j, &j, 1.0);
        }

        Self {
            rows,
            weights,
            dof,
            free,
            h,
        }
    }

    /// The number of points this was built from.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Whether there are no points at all.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// The degrees of freedom that were free when this was built.
    pub fn dof(&self) -> Dof3 {
        self.dof
    }

    /// How many degrees of freedom are free, which is the rank the information matrix would have
    /// if the points constrained the pose completely.
    pub fn free_dof_count(&self) -> usize {
        self.free.len()
    }

    /// The information matrix, expanded back into the full three parameters with zeros in the rows
    /// and columns of any locked degree of freedom.
    ///
    /// Note the unit caveat in the module documentation: the translation block is dimensionless,
    /// the rotation entry carries units of length squared, and the cross terms fall in between.
    pub fn information(&self) -> Matrix3<f64> {
        let mut full = Matrix3::zeros();
        for (a, &i) in self.free.iter().enumerate() {
            for (b, &j) in self.free.iter().enumerate() {
                full[(i, j)] = self.h[(a, b)];
            }
        }
        full
    }

    /// The marginal precision of each degree of freedom, `1 / [H^-1]_kk`, or `None` for a degree
    /// of freedom that was locked.
    ///
    /// This is the precision of a parameter *accounting for its correlation with the others*,
    /// which is the quantity that matters and is not the same as the diagonal of `H`. The
    /// distinction is what a perturbation experiment would measure: displace the points along one
    /// degree of freedom, let the other two re-optimize to absorb as much of it as they can, and
    /// the residual sum of squares that survives is `delta^2` times this value. A degree of
    /// freedom that looks well constrained on the diagonal of `H` can still be nearly free if
    /// another degree of freedom can imitate it.
    ///
    /// Larger is better constrained. The units differ between translation and rotation entries.
    ///
    /// Returns an error if the points do not constrain every free degree of freedom, since then
    /// `H` is singular and no finite precision exists. [`AlignInformation2::weak_directions`]
    /// still works in that case and will show which motions are free.
    pub fn marginal_precision(&self) -> Result<[Option<f64>; N_PARAMS]> {
        let inverse = self.h.clone().try_inverse().ok_or_else(|| {
            format!(
                "the information matrix is singular: the points do not constrain all {} free \
                 degrees of freedom. Use `weak_directions` to see which motions are free.",
                self.free.len()
            )
        })?;

        let mut out = [None; N_PARAMS];
        for (a, &i) in self.free.iter().enumerate() {
            out[i] = Some(1.0 / inverse[(a, a)]);
        }
        Ok(out)
    }

    /// The eigen-decomposition of the information matrix, as `(eigenvalue, direction)` pairs
    /// sorted from least to most constrained.
    ///
    /// The first entry is the motion the point set does least to prevent. This is usually more
    /// informative than the three per-axis numbers, because the weak direction of an alignment is
    /// rarely axis-aligned: points on a circle can spin about its center, points on a straight
    /// line slide along it, and a shallow arc is weakest in a coupled turn-and-shift. An
    /// eigenvalue at or near zero means that motion is unconstrained.
    ///
    /// Directions are unit vectors over the three parameters, with zeros in any locked slot. See
    /// the module documentation for why the relative scaling of the translation and rotation
    /// components depends on where the local origin was placed.
    pub fn weak_directions(&self) -> Vec<(f64, AlignStorage2)> {
        if self.free.is_empty() {
            return Vec::new();
        }

        let eigen = self.h.clone().symmetric_eigen();
        let mut pairs: Vec<(f64, AlignStorage2)> = eigen
            .eigenvalues
            .iter()
            .enumerate()
            .map(|(a, &value)| {
                let column = eigen.eigenvectors.column(a);
                let mut direction = AlignStorage2::zeros();
                for (b, &i) in self.free.iter().enumerate() {
                    direction[i] = column[b];
                }
                (value, direction)
            })
            .collect();

        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        pairs
    }

    /// The leverage of each point, `w_i * j_i^T H^-1 j_i`.
    ///
    /// Leverage is the classical measure of how much an observation influences its own fitted
    /// value. It is self-normalizing: the values sum to the number of free degrees of freedom, so
    /// the average is `free_dof_count / len()` and anything well above that is a point the
    /// alignment leans on.
    ///
    /// This scores each point *independently*, which makes it a good diagnostic but a poor basis
    /// for pruning, since a pair of interchangeable points both score high. Use
    /// [`AlignInformation2::select_d_optimal`] to choose a subset.
    ///
    /// Returns an error if the information matrix is singular, for the same reason as
    /// [`AlignInformation2::marginal_precision`].
    pub fn leverage(&self) -> Result<Vec<f64>> {
        let inverse = self.h.clone().try_inverse().ok_or_else(|| {
            format!(
                "the information matrix is singular: the points do not constrain all {} free \
                 degrees of freedom, so leverage is undefined.",
                self.free.len()
            )
        })?;

        Ok(self
            .rows
            .iter()
            .zip(self.weights.iter())
            .map(|(row, &w)| {
                let j = restrict(row, &self.free);
                w * (j.transpose() * &inverse * &j)[(0, 0)]
            })
            .collect())
    }

    /// Greedily chooses a subset of `count` points which preserves the conditioning of the
    /// alignment as well as a subset of that size can.
    ///
    /// The criterion is D-optimality: maximize the determinant of the information matrix of the
    /// chosen subset. Adding a point multiplies that determinant by `1 + w_i j_i^T H^-1 j_i`, so
    /// each round picks the point with the largest such marginal gain and folds it in with a
    /// Sherman-Morrison rank-1 update of the running inverse. Because a candidate is measured
    /// against what has already been chosen rather than in isolation, redundant points de-select
    /// themselves: the second of two interchangeable points has almost nothing left to add.
    ///
    /// Cost is `O(count * len())` with a fixed `k` by `k` update, where `k` is at most three.
    ///
    /// The running inverse has to start somewhere invertible, so it begins at a small multiple of
    /// the identity scaled to the data (see `ridge_rel`). This is the standard Bayesian form of
    /// D-optimality, and its influence is negligible once `k` points have been selected.
    ///
    /// Returns indices into the original point set, in the order they were chosen, so truncating
    /// the result yields the best smaller subset. Fewer than `count` indices come back if the
    /// candidates run out or if no remaining point can add anything.
    ///
    /// # Arguments
    ///
    /// * `count`: how many points to choose. Values above the number of available points, or
    ///   below the free degree-of-freedom count, are permitted but the latter cannot produce a
    ///   well-conditioned subset.
    /// * `ridge_rel`: the size of the initial regularizer relative to the mean diagonal of the
    ///   full information matrix. `None` uses a sensible default; there is rarely a reason to
    ///   set it.
    pub fn select_d_optimal(&self, count: usize, ridge_rel: Option<f64>) -> Vec<usize> {
        let k = self.free.len();
        if k == 0 || count == 0 || self.rows.is_empty() {
            return Vec::new();
        }

        // Scale the ridge to the problem so it stays negligible regardless of the units in play.
        // A completely empty information matrix (every weight zero) falls back to unity.
        let scale = self.h.trace() / k as f64;
        let ridge = ridge_rel.unwrap_or(DEFAULT_RIDGE_REL) * if scale > 0.0 { scale } else { 1.0 };

        // The inverse of `ridge * I` is `I / ridge`.
        let mut inverse = DMatrix::identity(k, k) / ridge;

        // Precompute the restricted rows once; the inner loop touches all of them every round.
        let restricted: Vec<DVector<f64>> =
            self.rows.iter().map(|r| restrict(r, &self.free)).collect();

        let mut chosen = Vec::with_capacity(count.min(self.rows.len()));
        let mut taken = vec![false; self.rows.len()];

        for _ in 0..count.min(self.rows.len()) {
            let mut best = None;
            let mut best_gain = 0.0;

            for i in 0..self.rows.len() {
                if taken[i] || self.weights[i] <= 0.0 {
                    continue;
                }
                let j = &restricted[i];
                let gain = self.weights[i] * (j.transpose() * &inverse * j)[(0, 0)];
                if gain > best_gain {
                    best_gain = gain;
                    best = Some(i);
                }
            }

            let Some(i) = best else { break };

            // Sherman-Morrison for `H + w j j^T`:
            //     (H + w j j^T)^-1 = H^-1 - (w H^-1 j j^T H^-1) / (1 + w j^T H^-1 j)
            let j = &restricted[i];
            let hj = &inverse * j;
            let denom = 1.0 + best_gain;
            inverse -= (&hj * hj.transpose()) * (self.weights[i] / denom);

            taken[i] = true;
            chosen.push(i);
        }

        chosen
    }
}

/// The indices of the free degrees of freedom, in parameter order.
fn free_indices(dof: &Dof3) -> Vec<usize> {
    [dof.tx, dof.ty, dof.rz]
        .iter()
        .enumerate()
        .filter_map(|(i, &free)| free.then_some(i))
        .collect()
}

/// Pulls the free components out of a full three-parameter jacobian row.
fn restrict(row: &AlignStorage2, free: &[usize]) -> DVector<f64> {
    DVector::from_iterator(free.len(), free.iter().map(|&i| row[i]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::align2::{AlignOrigin2, AlignSurfMatch2};
    use crate::{Curve2, Iso2, UnitVec2, Vector2};
    use approx::assert_relative_eq;

    /// A flat target along `y = 0`, so every match reports the same upward normal. Its
    /// information content is known in closed form, which is what makes it useful here.
    struct GroundLine;

    impl SurfaceTarget2 for GroundLine {
        fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
            AlignSurfMatch2::new(
                Point2::new(p.x, 0.0),
                UnitVec2::new_unchecked(Vector2::y()),
                true,
                1.0,
            )
        }
    }

    /// Points spread along `y = 1`, one unit above `GroundLine`.
    fn line_points() -> Vec<Point2> {
        (-6..=6).map(|i| Point2::new(i as f64, 1.0)).collect()
    }

    /// A closed counter-clockwise rectangle 10 by 5, centered on the origin. Being closed, it
    /// constrains every degree of freedom.
    fn rect_curve() -> Curve2 {
        let points =
            [(-5.0, -2.5), (5.0, -2.5), (5.0, 2.5), (-5.0, 2.5)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    /// Points walked around the rectangle at a uniform arc length spacing. The perimeter is 30,
    /// so a spacing of 0.3 gives 100 of them, and the first twenty all sit on the bottom edge.
    fn rect_points(curve: &Curve2) -> Vec<Point2> {
        let spacing = 0.3;
        let n = (curve.length() / spacing).floor() as usize;
        (0..n)
            .filter_map(|k| {
                curve
                    .at_length((k as f64 + 0.5) * spacing)
                    .map(|s| s.point())
            })
            .collect()
    }

    // ============================================================================================
    // Conditioning against geometry whose answer is known in advance
    // ============================================================================================

    #[test]
    fn a_straight_line_leaves_one_motion_unconstrained() {
        // Points against a straight line can only ever resist motion across it. Sliding along the
        // line moves every point parallel to the surface and changes no residual, so exactly one
        // of the three degrees of freedom is free. Rotation is not free, because turning the set
        // lifts points away from the line in proportion to how far out they sit.
        let info = AlignInformation2::from_points(
            &line_points(),
            &GroundLine,
            &AlignParams2::from_origin(None),
        );

        let weak = info.weak_directions();
        assert_eq!(weak.len(), 3);

        assert!(
            weak[0].0.abs() < 1e-9,
            "expected an unconstrained motion, got {}",
            weak[0].0
        );
        assert!(
            weak[1].0 > 1.0,
            "expected the remaining motions to be constrained, got {}",
            weak[1].0
        );

        // ...and that motion is sliding along the line, with no ty or rz component at all.
        let direction = weak[0].1;
        assert_relative_eq!(direction[0].abs(), 1.0, epsilon = 1e-9);
        assert!(direction[1].abs() < 1e-9);
        assert!(direction[2].abs() < 1e-9);
    }

    #[test]
    fn a_singular_problem_reports_rather_than_lying() {
        let info = AlignInformation2::from_points(
            &line_points(),
            &GroundLine,
            &AlignParams2::from_origin(None),
        );

        let err = info.marginal_precision().unwrap_err().to_string();
        assert!(err.contains("singular"), "unexpected message: {err}");
        assert!(info.leverage().is_err());
    }

    #[test]
    fn locking_the_free_motion_makes_a_line_well_conditioned() {
        // The same points, but with the motion they cannot see locked out. What remains is fully
        // determined, so the analysis should now succeed.
        let dof = Dof3::new(false, true, true);
        let params = AlignParams2::from_origin(Some(dof));
        let info = AlignInformation2::from_points(&line_points(), &GroundLine, &params);

        assert_eq!(info.free_dof_count(), 2);

        let precision = info.marginal_precision().unwrap();
        assert_eq!(precision[0], None, "tx was locked");
        for i in [1usize, 2] {
            assert!(
                precision[i].unwrap() > 0.0,
                "parameter {i} should be constrained"
            );
        }
    }

    #[test]
    fn a_closed_curve_constrains_everything() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let params = AlignParams2::from_origin(None);
        let info = AlignInformation2::from_points(&points, &curve, &params);

        assert_eq!(info.free_dof_count(), 3);
        let precision = info.marginal_precision().unwrap();
        for (i, p) in precision.iter().enumerate() {
            assert!(
                p.unwrap() > 0.0,
                "parameter {i} should be constrained on a closed curve"
            );
        }
    }

    // ============================================================================================
    // Leverage
    // ============================================================================================

    #[test]
    fn leverage_sums_to_the_free_dof_count() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 3.0, epsilon = 1e-9);
    }

    #[test]
    fn leverage_sums_to_the_free_dof_count_when_constrained() {
        // The identity `sum(h_i) == k` has to hold on the free subspace too, not just at k = 3.
        let curve = rect_curve();
        let points = rect_points(&curve);
        let dof = Dof3::new(false, true, true);
        let info =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(Some(dof)));

        assert_eq!(info.free_dof_count(), 2);
        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 2.0, epsilon = 1e-9);
    }

    // ============================================================================================
    // Greedy selection
    // ============================================================================================

    #[test]
    fn selection_prefers_isolated_points_over_a_redundant_cluster() {
        // This is the property that separates greedy selection from independent scoring. A tight
        // cluster of near-identical points all contribute the same row, so after the first one is
        // taken the rest add nothing; the isolated points elsewhere on the rectangle each bring
        // something new and should be chosen first.
        let curve = rect_curve();

        let mut points = Vec::new();

        // Twenty near-duplicates in one spot on the bottom edge.
        for i in 0..20 {
            points.push(Point2::new(1.0 + i as f64 * 1e-4, -2.5));
        }
        let cluster = 0..20;

        // Four isolated points spread over the other three edges.
        let isolated = [
            Point2::new(5.0, 0.0),
            Point2::new(-5.0, 1.0),
            Point2::new(0.0, 2.5),
            Point2::new(-3.0, 2.5),
        ];
        points.extend_from_slice(&isolated);

        let info =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));
        // Four picks, which is how many distinct positions are on offer once the cluster has
        // contributed the one row it has to give. Twenty of the twenty-four candidates are
        // cluster members, so picking arbitrarily would take about three of them.
        let picked = info.select_d_optimal(4, None);
        assert_eq!(picked.len(), 4);

        let from_cluster = picked.iter().filter(|&&i| cluster.contains(&i)).count();
        assert_eq!(
            from_cluster, 1,
            "expected exactly one cluster member among {picked:?}"
        );

        // Asked for more than the geometry can offer, it does go back to the cluster: with three
        // degrees of freedom and every distinct position already taken, a near-duplicate is the
        // best remaining candidate. That is the selection working, not failing.
        let more = info.select_d_optimal(5, None);
        assert_eq!(
            more.iter().filter(|&&i| cluster.contains(&i)).count(),
            2,
            "expected the fifth pick to fall back to the cluster: {more:?}"
        );
    }

    #[test]
    fn selection_is_a_prefix_so_it_can_be_truncated() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let long = info.select_d_optimal(20, None);
        let short = info.select_d_optimal(8, None);

        assert_eq!(&long[..8], &short[..]);
    }

    #[test]
    fn a_selected_subset_stays_well_conditioned() {
        // The point of pruning: a small chosen subset should hold the alignment up nearly as well
        // as the whole set, and far better than the same number of arbitrary points.
        let curve = rect_curve();
        let points = rect_points(&curve);
        let params = AlignParams2::from_origin(None);
        let info = AlignInformation2::from_points(&points, &curve, &params);

        let n = 20;
        let picked = info.select_d_optimal(n, None);
        let chosen: Vec<Point2> = picked.iter().map(|&i| points[i]).collect();
        let arbitrary: Vec<Point2> = points.iter().take(n).copied().collect();

        // The comparison is on the smallest eigenvalue, which measures how close the subset comes
        // to leaving some motion unconstrained. It is used rather than `marginal_precision`
        // because an arbitrary subset can be outright singular (the first `n` samples walk along
        // a single edge of the rectangle), and there is no finite precision to compare then.
        let weakest = |pts: &[Point2]| -> f64 {
            AlignInformation2::from_points(pts, &curve, &params).weak_directions()[0].0
        };

        let chosen_weakest = weakest(&chosen);
        let arbitrary_weakest = weakest(&arbitrary);

        assert!(
            chosen_weakest > arbitrary_weakest * 10.0,
            "the selected subset should be far better conditioned than an arbitrary one: \
             selected {chosen_weakest}, arbitrary {arbitrary_weakest}"
        );

        // ...and unlike the arbitrary subset, it should constrain the pose outright.
        assert!(
            AlignInformation2::from_points(&chosen, &curve, &params)
                .marginal_precision()
                .is_ok(),
            "the selected subset should leave no motion unconstrained"
        );
    }

    #[test]
    fn selection_respects_the_available_count() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        assert!(info.select_d_optimal(0, None).is_empty());
        assert_eq!(
            info.select_d_optimal(points.len() + 50, None).len(),
            points.len()
        );
    }

    #[test]
    fn zero_weight_points_are_never_selected() {
        let rows = vec![
            AlignStorage2::new(1.0, 0.0, 0.0),
            AlignStorage2::new(0.0, 1.0, 0.0),
            AlignStorage2::new(0.0, 0.0, 1.0),
        ];
        let weights = vec![1.0, 0.0, 1.0];

        let info = AlignInformation2::from_rows(rows, weights, Dof3::all()).unwrap();
        let picked = info.select_d_optimal(3, None);

        assert!(!picked.contains(&1), "a zero-weight point was selected");
        assert_eq!(picked.len(), 2);
    }

    // ============================================================================================
    // Construction and validation
    // ============================================================================================

    #[test]
    fn information_matrix_is_expanded_with_zeros_for_locked_dof() {
        let rows = vec![AlignStorage2::new(1.0, 2.0, 0.0)];
        let dof = Dof3::new(true, true, false);
        let info = AlignInformation2::from_rows(rows, vec![1.0], dof).unwrap();

        let h = info.information();
        // The free block is the outer product of (1, 2) with itself.
        assert_relative_eq!(h[(0, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(h[(0, 1)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(h[(1, 1)], 4.0, epsilon = 1e-12);
        // ...and every locked row and column is zero.
        for j in 0..3 {
            assert_eq!(h[(2, j)], 0.0);
            assert_eq!(h[(j, 2)], 0.0);
        }
    }

    #[test]
    fn mismatched_rows_and_weights_are_rejected() {
        let rows = vec![AlignStorage2::zeros(); 3];
        let err = AlignInformation2::from_rows(rows, vec![1.0; 2], Dof3::all())
            .unwrap_err()
            .to_string();
        assert!(err.contains("weights"), "unexpected message: {err}");
    }

    #[test]
    fn invalid_weights_are_rejected() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let rows = vec![AlignStorage2::zeros(); 2];
            let result = AlignInformation2::from_rows(rows, vec![1.0, bad], Dof3::all());
            assert!(result.is_err(), "weight {bad} should have been rejected");
        }
    }

    #[test]
    fn a_fully_locked_problem_degrades_gracefully() {
        let dof = Dof3::new(false, false, false);
        let info = AlignInformation2::from_rows(vec![AlignStorage2::zeros(); 4], vec![1.0; 4], dof)
            .unwrap();

        assert_eq!(info.free_dof_count(), 0);
        assert!(info.weak_directions().is_empty());
        assert!(info.select_d_optimal(4, None).is_empty());
        assert_eq!(info.information(), Matrix3::zeros());
    }

    #[test]
    fn the_local_origin_is_respected() {
        // The rotation partial is taken about the local origin, so moving it changes the rotation
        // entry of the information matrix. This guards against the origin being quietly ignored.
        let curve = rect_curve();
        let points = rect_points(&curve);

        let at_origin =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));
        let offset = AlignInformation2::from_points(
            &points,
            &curve,
            &AlignParams2::from_center(Point2::new(50.0, 0.0), None),
        );

        let delta = (at_origin.information() - offset.information()).norm();
        assert!(
            delta > 1.0,
            "the local origin had no effect (delta {delta})"
        );
    }

    #[test]
    fn a_transformed_pose_is_analyzed_where_it_currently_sits() {
        // The analysis describes the pose described by `params`, not the pose the points were
        // authored in, so a working offset has to be honored.
        let curve = rect_curve();
        let points = rect_points(&curve);

        let here =
            AlignInformation2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let displaced = AlignParams2::new(
            AlignOrigin2::Origin,
            Some(Iso2::translation(0.0, 20.0)),
            None,
        );
        let there = AlignInformation2::from_points(&points, &curve, &displaced);

        let delta = (here.information() - there.information()).norm();
        assert!(
            delta > 1e-6,
            "the working offset had no effect (delta {delta})"
        );
    }
}
