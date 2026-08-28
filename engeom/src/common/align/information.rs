//! The dimension-generic core of the alignment information analysis.
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
//! Nothing about that depends on the dimension: the analysis runs on jacobian rows of `N`
//! parameters, whatever those parameters mean. [`AlignInformation`] holds it once, generic over
//! `N` and the degree-of-freedom lock type `D`, whose [`AlignDofs`] implementation is the only
//! dimension-specific ingredient. The public faces are the aliases
//! `geom2::align2::AlignInformation2` and `geom3::align3::AlignInformation3`, whose modules
//! supply the dimension-specific `from_points` constructors, geometric intuition, and unit
//! caveats.

use crate::Result;
use crate::na::{DMatrix, DVector, SMatrix, SVector};

/// The relative size of the regularizing ridge used to start greedy selection. See
/// [`AlignInformation::select_d_optimal`] for why one is needed.
const DEFAULT_RIDGE_REL: f64 = 1e-9;

/// The degree-of-freedom locks of a dimension, viewed as a subset of its `N` parameters.
///
/// Implemented by `Dof3` (2D: `tx`, `ty`, `rz`) and `Dof6` (3D: `tx`, `ty`, `tz`, `rx`, `ry`,
/// `rz`).
pub trait AlignDofs<const N: usize>: Copy {
    /// The indices of the free degrees of freedom, in parameter order.
    fn free_indices(&self) -> Vec<usize>;
}

/// The information content of a set of test points with respect to an alignment.
///
/// Built through the dimension-specific `from_points` constructors on the `AlignInformation2` and
/// `AlignInformation3` aliases, which project points onto a target and take the jacobian row of
/// each correspondence, or with [`AlignInformation::from_rows`] when the rows have already been
/// computed. See the alias modules' documentation for what it is good for.
#[derive(Clone, Debug)]
pub struct AlignInformation<D: AlignDofs<N>, const N: usize> {
    /// The jacobian row for each point. Locked degrees of freedom are already zero in these,
    /// because the jacobian functions zero their columns.
    rows: Vec<SVector<f64, N>>,

    /// The weight of each point's contribution to the information matrix.
    weights: Vec<f64>,

    /// Which degrees of freedom are free to move.
    dof: D,

    /// The indices, within the `N` parameters, of the degrees of freedom which are free. All of
    /// the linear algebra happens in this subspace, because `H` is identically zero along a locked
    /// direction and would not be invertible in the full space.
    free: Vec<usize>,

    /// The information matrix in the free subspace, `k` by `k` where `k = free.len()`.
    h: DMatrix<f64>,
}

impl<D: AlignDofs<N>, const N: usize> AlignInformation<D, N> {
    /// Builds the information content directly from jacobian rows which have already been
    /// computed, for callers assembling a problem this module doesn't know how to build (a
    /// multi-body adjustment, for instance).
    ///
    /// Rows are expected to already have their locked columns zeroed, which is what the jacobian
    /// functions produce when given a constrained degree-of-freedom set.
    ///
    /// Returns an error if the two slices disagree in length, or if any weight is negative or
    /// non-finite.
    pub fn from_rows(rows: Vec<SVector<f64, N>>, weights: Vec<f64>, dof: D) -> Result<Self> {
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

    pub(crate) fn assemble(rows: Vec<SVector<f64, N>>, weights: Vec<f64>, dof: D) -> Self {
        let free = dof.free_indices();
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
    pub fn dof(&self) -> D {
        self.dof
    }

    /// How many degrees of freedom are free, which is the rank the information matrix would have
    /// if the points constrained the pose completely.
    pub fn free_dof_count(&self) -> usize {
        self.free.len()
    }

    /// The information matrix, expanded back into the full `N` parameters with zeros in the rows
    /// and columns of any locked degree of freedom.
    ///
    /// Note the unit caveat in the alias module's documentation: the translation block is
    /// dimensionless, the rotation entries carry units of length squared, and the cross terms
    /// fall in between.
    pub fn information(&self) -> SMatrix<f64, N, N> {
        let mut full = SMatrix::zeros();
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
    /// degree of freedom, let the others re-optimize to absorb as much of it as they can, and
    /// the residual sum of squares that survives is `delta^2` times this value. A degree of
    /// freedom that looks well constrained on the diagonal of `H` can still be nearly free if
    /// another degree of freedom can imitate it.
    ///
    /// Larger is better constrained. The units differ between translation and rotation entries.
    ///
    /// Returns an error if the points do not constrain every free degree of freedom, since then
    /// `H` is singular and no finite precision exists. [`AlignInformation::weak_directions`]
    /// still works in that case and will show which motions are free.
    pub fn marginal_precision(&self) -> Result<[Option<f64>; N]> {
        let inverse = self.h.clone().try_inverse().ok_or_else(|| {
            format!(
                "the information matrix is singular: the points do not constrain all {} free \
                 degrees of freedom. Use `weak_directions` to see which motions are free.",
                self.free.len()
            )
        })?;

        let mut out = [None; N];
        for (a, &i) in self.free.iter().enumerate() {
            out[i] = Some(1.0 / inverse[(a, a)]);
        }
        Ok(out)
    }

    /// The eigen-decomposition of the information matrix, as `(eigenvalue, direction)` pairs
    /// sorted from least to most constrained.
    ///
    /// The first entry is the motion the point set does least to prevent. This is usually more
    /// informative than the per-axis numbers, because the weak direction of an alignment is
    /// rarely axis-aligned; see the alias modules' documentation for geometric examples in each
    /// dimension. An eigenvalue at or near zero means that motion is unconstrained.
    ///
    /// Directions are unit vectors over the `N` parameters, with zeros in any locked slot. See
    /// the alias module's documentation for why the relative scaling of the translation and
    /// rotation components depends on where the local origin was placed.
    pub fn weak_directions(&self) -> Vec<(f64, SVector<f64, N>)> {
        if self.free.is_empty() {
            return Vec::new();
        }

        let eigen = self.h.clone().symmetric_eigen();
        let mut pairs: Vec<(f64, SVector<f64, N>)> = eigen
            .eigenvalues
            .iter()
            .enumerate()
            .map(|(a, &value)| {
                let column = eigen.eigenvectors.column(a);
                let mut direction = SVector::zeros();
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
    /// [`AlignInformation::select_d_optimal`] to choose a subset.
    ///
    /// Returns an error if the information matrix is singular, for the same reason as
    /// [`AlignInformation::marginal_precision`].
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
    /// Cost is `O(count * len())` with a fixed `k` by `k` update, where `k` is at most `N`.
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

/// Pulls the free components out of a full `N`-parameter jacobian row.
fn restrict<const N: usize>(row: &SVector<f64, N>, free: &[usize]) -> DVector<f64> {
    DVector::from_iterator(free.len(), free.iter().map(|&i| row[i]))
}
