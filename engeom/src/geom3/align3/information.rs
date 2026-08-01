//! Tools for asking how well a set of points constrains a 3D alignment, and for choosing a
//! subset of them that preserves that conditioning.
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
//!   to slide? See [`AlignInformation3::marginal_precision`] and
//!   [`AlignInformation3::weak_directions`].
//! - **Pruning.** Which points are actually carrying the alignment, and what is the smallest
//!   subset that would constrain it just as well? See [`AlignInformation3::leverage`] and
//!   [`AlignInformation3::select_d_optimal`].
//!
//! Pruning is the reason this exists. A simultaneous alignment of fifteen or twenty scans is
//! expensive in proportion to the number of alignment points, and most of those points are
//! redundant with their neighbors.
//!
//! # Why scoring points individually is not enough
//!
//! The obvious approach is to score each point on its own and keep the best `n`. That
//! double-counts redundancy: two points sitting on the same flat patch each score well, but once
//! either is in the set the other adds almost nothing. Greedy D-optimal selection avoids this for
//! free, because a candidate is scored by its *marginal* gain against everything already chosen,
//! and that gain collapses as soon as something equivalent has been picked.
//!
//! # A caveat about units
//!
//! The translation columns of `H` are dimensionless (a dot product of two unit vectors) while the
//! rotation columns carry units of length, because a rotation partial scales with the distance
//! from the local origin. `H` is therefore heterogeneous, and two consequences follow.
//!
//! [`AlignInformation3::marginal_precision`] is safe: each entry concerns a single parameter, so
//! its units are internally consistent. [`AlignInformation3::weak_directions`] mixes translation
//! and rotation in one vector and is only meaningful when the two are comparably scaled, which in
//! practice means the local origin should sit near the middle of the point set (as
//! `AlignParams3::from_center` arranges). Placing the origin far away inflates the rotation block
//! and will make rotations look artificially well constrained.

use crate::geom3::align3::jacobian::point_surf_jacobian;
use crate::geom3::align3::{AlignParams3, AlignStorage3, Dof6, SurfaceTarget3};
use crate::na::{DMatrix, DVector, Matrix6};
use crate::{Point3, Result};

/// The number of alignment parameters (tx, ty, tz, rx, ry, rz).
const N_PARAMS: usize = 6;

/// The relative size of the regularizing ridge used to start greedy selection. See
/// [`AlignInformation3::select_d_optimal`] for why one is needed.
const DEFAULT_RIDGE_REL: f64 = 1e-9;

/// The information content of a set of test points with respect to a 3D alignment.
///
/// Built with [`AlignInformation3::from_points`], which projects the points onto a target and
/// takes the jacobian row of each correspondence, or with [`AlignInformation3::from_rows`] when
/// the rows have already been computed. See the module documentation for what it is good for.
#[derive(Clone, Debug)]
pub struct AlignInformation3 {
    /// The jacobian row for each point. Locked degrees of freedom are already zero in these,
    /// because `point_surf_jacobian` zeroes their columns.
    rows: Vec<AlignStorage3>,

    /// The weight of each point's contribution to the information matrix.
    weights: Vec<f64>,

    /// Which degrees of freedom are free to move.
    dof: Dof6,

    /// The indices, within the six parameters, of the degrees of freedom which are free. All of
    /// the linear algebra happens in this subspace, because `H` is identically zero along a locked
    /// direction and would not be invertible in the full space.
    free: Vec<usize>,

    /// The information matrix in the free subspace, `k` by `k` where `k = free.len()`.
    h: DMatrix<f64>,
}

impl AlignInformation3 {
    /// Builds the information content of a set of points against a target, in the pose described
    /// by `params`.
    ///
    /// Each point is moved by the current transform, projected onto the target, and turned into a
    /// jacobian row exactly as [`crate::geom3::align3::points_to_surface3`] would do. A point is
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
    ///   partials are taken about and which degrees of freedom are free
    pub fn from_points(
        points: &[Point3],
        target: &impl SurfaceTarget3,
        params: &AlignParams3,
    ) -> Self {
        let current = params.compute_values();
        let mut rows = Vec::with_capacity(points.len());
        let mut weights = Vec::with_capacity(points.len());

        for p in points {
            let m = current.transform * p;
            let c = target.find_align_match(&m);

            rows.push(point_surf_jacobian(&m, &c, &current));
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
    /// `point_surf_jacobian` produces when given a constrained [`Dof6`].
    ///
    /// Returns an error if the two slices disagree in length, or if any weight is negative or
    /// non-finite.
    pub fn from_rows(rows: Vec<AlignStorage3>, weights: Vec<f64>, dof: Dof6) -> Result<Self> {
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

    fn assemble(rows: Vec<AlignStorage3>, weights: Vec<f64>, dof: Dof6) -> Self {
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
    pub fn dof(&self) -> Dof6 {
        self.dof
    }

    /// How many degrees of freedom are free, which is the rank the information matrix would have
    /// if the points constrained the pose completely.
    pub fn free_dof_count(&self) -> usize {
        self.free.len()
    }

    /// The information matrix, expanded back into the full six parameters with zeros in the rows
    /// and columns of any locked degree of freedom.
    ///
    /// Note the unit caveat in the module documentation: the translation block is dimensionless,
    /// the rotation block carries units of length squared, and the cross terms fall in between.
    pub fn information(&self) -> Matrix6<f64> {
        let mut full = Matrix6::zeros();
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
    /// degree of freedom, let the other five re-optimize to absorb as much of it as they can, and
    /// the residual sum of squares that survives is `delta^2` times this value. A degree of
    /// freedom that looks well constrained on the diagonal of `H` can still be nearly free if
    /// another degree of freedom can imitate it.
    ///
    /// Larger is better constrained. The units differ between translation and rotation entries.
    ///
    /// Returns an error if the points do not constrain every free degree of freedom, since then
    /// `H` is singular and no finite precision exists. [`AlignInformation3::weak_directions`]
    /// still works in that case and will show which motions are unconstrained.
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
    /// informative than the six per-axis numbers, because the weak direction of an alignment is
    /// rarely axis-aligned: points on a cylinder slide along a helix, points on a plane slide in
    /// any in-plane direction and spin about the normal, and a shallow dish is weakest in a
    /// coupled tilt-and-shift. An eigenvalue at or near zero means that motion is unconstrained.
    ///
    /// Directions are unit vectors over the six parameters, with zeros in any locked slot. See the
    /// module documentation for why the relative scaling of the translation and rotation
    /// components depends on where the local origin was placed.
    pub fn weak_directions(&self) -> Vec<(f64, AlignStorage3)> {
        if self.free.is_empty() {
            return Vec::new();
        }

        let eigen = self.h.clone().symmetric_eigen();
        let mut pairs: Vec<(f64, AlignStorage3)> = eigen
            .eigenvalues
            .iter()
            .enumerate()
            .map(|(a, &value)| {
                let column = eigen.eigenvectors.column(a);
                let mut direction = AlignStorage3::zeros();
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
    /// [`AlignInformation3::select_d_optimal`] to choose a subset.
    ///
    /// Returns an error if the information matrix is singular, for the same reason as
    /// [`AlignInformation3::marginal_precision`].
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
    /// Cost is `O(count * len())` with a fixed `k` by `k` update, where `k` is at most six.
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
fn free_indices(dof: &Dof6) -> Vec<usize> {
    [dof.tx, dof.ty, dof.tz, dof.rx, dof.ry, dof.rz]
        .iter()
        .enumerate()
        .filter_map(|(i, &free)| free.then_some(i))
        .collect()
}

/// Pulls the free components out of a full six-parameter jacobian row.
fn restrict(row: &AlignStorage3, free: &[usize]) -> DVector<f64> {
    DVector::from_iterator(free.len(), free.iter().map(|&i| row[i]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::align3::AlignSurfMatch3;
    use crate::{Iso3, Mesh3, UnitVec3, Vector3};
    use approx::assert_relative_eq;

    /// A flat target in the z = 0 plane, so every match reports the same upward normal. Its
    /// information content is known in closed form, which is what makes it useful here.
    struct GroundPlane;

    impl SurfaceTarget3 for GroundPlane {
        fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
            AlignSurfMatch3::new(
                Point3::new(p.x, p.y, 0.0),
                UnitVec3::new_unchecked(Vector3::z()),
                true,
                1.0,
            )
        }
    }

    /// Points scattered over a patch of the z = 1 plane, one unit above `GroundPlane`.
    fn plane_points() -> Vec<Point3> {
        let mut points = Vec::new();
        for i in -3..=3 {
            for j in -3..=3 {
                points.push(Point3::new(i as f64, j as f64, 1.0));
            }
        }
        points
    }

    fn box_mesh() -> Mesh3 {
        Mesh3::create_box(10.0, 5.0, 2.0, false)
    }

    fn box_points(mesh: &Mesh3) -> Vec<Point3> {
        mesh.sample_poisson(0.5, None)
            .iter()
            .map(|p| p.point())
            .collect()
    }

    // ============================================================================================
    // Conditioning against geometry whose answer is known in advance
    // ============================================================================================

    #[test]
    fn a_plane_leaves_three_motions_unconstrained() {
        // Points on a plane can only ever resist motion along its normal. Sliding in x or y, or
        // spinning about z, moves every point along the surface and changes no residual, so
        // exactly three of the six degrees of freedom are free.
        let info = AlignInformation3::from_points(
            &plane_points(),
            &GroundPlane,
            &AlignParams3::from_origin(None),
        );

        let weak = info.weak_directions();
        assert_eq!(weak.len(), 6);

        // Three eigenvalues at zero, three well away from it.
        let unconstrained: Vec<f64> = weak.iter().take(3).map(|(v, _)| *v).collect();
        for v in &unconstrained {
            assert!(v.abs() < 1e-9, "expected an unconstrained motion, got {v}");
        }
        assert!(
            weak[3].0 > 1.0,
            "expected the remaining motions to be constrained, got {}",
            weak[3].0
        );

        // ...and those three motions live entirely in the tx / ty / rz subspace. Any individual
        // eigenvector may mix them, so the test is that they have no tz / rx / ry component.
        for (_, direction) in weak.iter().take(3) {
            for i in [2usize, 3, 4] {
                assert!(
                    direction[i].abs() < 1e-9,
                    "unconstrained direction leaked into parameter {i}: {direction:?}"
                );
            }
        }
    }

    #[test]
    fn a_singular_problem_reports_rather_than_lying() {
        let info = AlignInformation3::from_points(
            &plane_points(),
            &GroundPlane,
            &AlignParams3::from_origin(None),
        );

        let err = info.marginal_precision().unwrap_err().to_string();
        assert!(err.contains("singular"), "unexpected message: {err}");
        assert!(info.leverage().is_err());
    }

    #[test]
    fn locking_the_free_motions_makes_a_plane_well_conditioned() {
        // The same points, but with the three motions they cannot see locked out. What remains is
        // fully determined, so the analysis should now succeed.
        let dof = Dof6::new(false, false, true, true, true, false);
        let params = AlignParams3::from_origin(Some(dof));
        let info = AlignInformation3::from_points(&plane_points(), &GroundPlane, &params);

        assert_eq!(info.free_dof_count(), 3);

        let precision = info.marginal_precision().unwrap();
        assert_eq!(precision[0], None, "tx was locked");
        assert_eq!(precision[1], None, "ty was locked");
        assert_eq!(precision[5], None, "rz was locked");
        for i in [2usize, 3, 4] {
            assert!(
                precision[i].unwrap() > 0.0,
                "parameter {i} should be constrained"
            );
        }
    }

    #[test]
    fn a_box_constrains_everything() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let params = AlignParams3::from_origin(None);
        let info = AlignInformation3::from_points(&points, &mesh, &params);

        assert_eq!(info.free_dof_count(), 6);
        let precision = info.marginal_precision().unwrap();
        for (i, p) in precision.iter().enumerate() {
            assert!(
                p.unwrap() > 0.0,
                "parameter {i} should be constrained on a closed box"
            );
        }
    }

    // ============================================================================================
    // Leverage
    // ============================================================================================

    #[test]
    fn leverage_sums_to_the_free_dof_count() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 6.0, epsilon = 1e-9);
    }

    #[test]
    fn leverage_sums_to_the_free_dof_count_when_constrained() {
        // The identity `sum(h_i) == k` has to hold on the free subspace too, not just at k = 6.
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let dof = Dof6::new(false, true, true, true, false, true);
        let info =
            AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(Some(dof)));

        assert_eq!(info.free_dof_count(), 4);
        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 4.0, epsilon = 1e-9);
    }

    // ============================================================================================
    // Greedy selection
    // ============================================================================================

    #[test]
    fn selection_prefers_isolated_points_over_a_redundant_cluster() {
        // This is the property that separates greedy selection from independent scoring. A tight
        // cluster of near-identical points all contribute the same row, so after the first one is
        // taken the rest add nothing; the isolated points elsewhere on the box each bring
        // something new and should be chosen first.
        let mesh = box_mesh();

        let mut points = Vec::new();

        // Twenty near-duplicates in one spot on the top face.
        for i in 0..20 {
            points.push(Point3::new(1.0 + i as f64 * 1e-4, 1.0, 1.0));
        }
        let cluster = 0..20;

        // Six isolated points spread over other faces.
        let isolated = [
            Point3::new(-4.0, -2.0, 1.0),
            Point3::new(4.0, 2.0, -1.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(-5.0, 1.0, 0.5),
            Point3::new(0.0, 2.5, -0.5),
            Point3::new(2.0, -2.5, 0.5),
        ];
        points.extend_from_slice(&isolated);

        let info = AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));
        let picked = info.select_d_optimal(7, None);

        assert_eq!(picked.len(), 7);

        // Every isolated point should be taken before a second member of the cluster is.
        let from_cluster = picked.iter().filter(|&&i| cluster.contains(&i)).count();
        assert_eq!(
            from_cluster, 1,
            "expected exactly one cluster member among {picked:?}"
        );
    }

    #[test]
    fn selection_is_a_prefix_so_it_can_be_truncated() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let long = info.select_d_optimal(20, None);
        let short = info.select_d_optimal(8, None);

        assert_eq!(&long[..8], &short[..]);
    }

    #[test]
    fn a_selected_subset_stays_well_conditioned() {
        // The point of pruning: a small chosen subset should hold the alignment up nearly as well
        // as the whole set, and far better than the same number of arbitrary points.
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let params = AlignParams3::from_origin(None);
        let info = AlignInformation3::from_points(&points, &mesh, &params);

        let n = 40;
        let picked = info.select_d_optimal(n, None);
        let chosen: Vec<Point3> = picked.iter().map(|&i| points[i]).collect();
        let arbitrary: Vec<Point3> = points.iter().take(n).copied().collect();

        // The comparison is on the smallest eigenvalue, which measures how close the subset comes
        // to leaving some motion unconstrained. It is used rather than `marginal_precision`
        // because an arbitrary subset can be outright singular (the first `n` Poisson samples
        // tend to land on a single face of the box), and there is no finite precision to compare
        // in that case.
        let weakest = |pts: &[Point3]| -> f64 {
            AlignInformation3::from_points(pts, &mesh, &params).weak_directions()[0].0
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
            AlignInformation3::from_points(&chosen, &mesh, &params)
                .marginal_precision()
                .is_ok(),
            "the selected subset should leave no motion unconstrained"
        );
    }

    #[test]
    fn selection_respects_the_available_count() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        assert!(info.select_d_optimal(0, None).is_empty());
        assert_eq!(
            info.select_d_optimal(points.len() + 50, None).len(),
            points.len()
        );
    }

    #[test]
    fn zero_weight_points_are_never_selected() {
        let rows = vec![
            AlignStorage3::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            AlignStorage3::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            AlignStorage3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        ];
        let weights = vec![1.0, 0.0, 1.0];
        let dof = Dof6::new(true, true, true, false, false, false);

        let info = AlignInformation3::from_rows(rows, weights, dof).unwrap();
        let picked = info.select_d_optimal(3, None);

        assert!(!picked.contains(&1), "a zero-weight point was selected");
        assert_eq!(picked.len(), 2);
    }

    // ============================================================================================
    // Construction and validation
    // ============================================================================================

    #[test]
    fn information_matrix_is_expanded_with_zeros_for_locked_dof() {
        let rows = vec![AlignStorage3::new(1.0, 2.0, 3.0, 0.0, 0.0, 0.0)];
        let dof = Dof6::new(true, true, true, false, false, false);
        let info = AlignInformation3::from_rows(rows, vec![1.0], dof).unwrap();

        let h = info.information();
        // The free block is the outer product of (1, 2, 3) with itself.
        assert_relative_eq!(h[(0, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(h[(0, 1)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(h[(2, 2)], 9.0, epsilon = 1e-12);
        // ...and every locked row and column is zero.
        for i in 3..6 {
            for j in 0..6 {
                assert_eq!(h[(i, j)], 0.0);
                assert_eq!(h[(j, i)], 0.0);
            }
        }
    }

    #[test]
    fn mismatched_rows_and_weights_are_rejected() {
        let rows = vec![AlignStorage3::zeros(); 3];
        let err = AlignInformation3::from_rows(rows, vec![1.0; 2], Dof6::all())
            .unwrap_err()
            .to_string();
        assert!(err.contains("weights"), "unexpected message: {err}");
    }

    #[test]
    fn invalid_weights_are_rejected() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let rows = vec![AlignStorage3::zeros(); 2];
            let result = AlignInformation3::from_rows(rows, vec![1.0, bad], Dof6::all());
            assert!(result.is_err(), "weight {bad} should have been rejected");
        }
    }

    #[test]
    fn a_fully_locked_problem_degrades_gracefully() {
        let dof = Dof6::new(false, false, false, false, false, false);
        let info = AlignInformation3::from_rows(vec![AlignStorage3::zeros(); 4], vec![1.0; 4], dof)
            .unwrap();

        assert_eq!(info.free_dof_count(), 0);
        assert!(info.weak_directions().is_empty());
        assert!(info.select_d_optimal(4, None).is_empty());
        assert_eq!(info.information(), Matrix6::zeros());
    }

    #[test]
    fn the_local_origin_is_respected() {
        // Rotation partials are taken about the local origin, so moving it changes the rotation
        // block of the information matrix. This guards against the origin being quietly ignored.
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let at_origin =
            AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));
        let offset = AlignInformation3::from_points(
            &points,
            &mesh,
            &AlignParams3::from_center(Point3::new(50.0, 0.0, 0.0), None),
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
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let here = AlignInformation3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let displaced = AlignParams3::new(
            crate::geom3::align3::AlignOrigin3::Origin,
            Some(Iso3::translation(0.0, 0.0, 20.0)),
            None,
        );
        let there = AlignInformation3::from_points(&points, &mesh, &displaced);

        let delta = (here.information() - there.information()).norm();
        assert!(
            delta > 1e-6,
            "the working offset had no effect (delta {delta})"
        );
    }
}
