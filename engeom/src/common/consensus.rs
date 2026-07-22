//! This module holds common tools used for consensus based robust estimation. My initial goal
//! (as of July '26) is to rework what are currently a bunch of scattered RANSAC creation methods
//! for some of the geometric primitives into something unified and based on a more up-to-date
//! approach. I also want to make sure that consensus based creation methods are provided for all
//! geometric primitives.
//!
//! I would like very much to remove the requirement for the caller to provide a fixed inlier
//! threshold, since I see that being the main obstacle to quick and efficient usage of the tools.
//! To that end, I think the MAGSAC/MAGSAC++ algorithms are probably the most suitable, as well as
//! being at or near the current state-of-the-art.
//!
//! I investigated the [Inlier](https://github.com/soraxas/inlier) library, which in turn is a Rust
//! port of the Python/C++ [SupeRANSAC](https://github.com/danini/superansac) library. It is less
//! of a building block, however, and more of a feature-complete tool built around several end
//! applications (homographies, fundamental/essential matrix estimation, etc.). It has a lot of
//! internal modularity to tune itself to different applications, which makes it hard to follow and
//! is overkill for what I need here. Lastly, it has a lot of dependencies that have version
//! mismatches with ones this library relies on, which gives me pause.
//!
//! # Marginalizing Sample Consensus
//!
//! [MAGSAC++, a fast, reliable, and accurate robust estimator](https://openaccess.thecvf.com/content_CVPR_2020/papers/Barath_MAGSAC_a_Fast_Reliable_and_Accurate_Robust_Estimator_CVPR_2020_paper.pdf)
//!
//! MAGSAC, and the authors' improvements in MAGSAC++, attempts to remove the inlier threshold by
//! substituting the hard threshold for a noise parameter σ and then marginalizing it out as a
//! nuisance parameter.  The MAGSAC++ improvements, in my understanding, make some assumptions
//! about distribution of noise in the inliers, and then use that to come up with a clever mechanism
//! of analytically integrating over σ.
//!
//! In practice, it seems that the differences from the original RANSAC happen in two places:
//!   - It uses an iteratively reweighted least squares refinement step when fitting each random
//!     model (which is the clever trick for analytically marginalizing σ mentioned above)
//!   - It uses a loss/quality function based on the same σ-marginalized interval as the weight
//!     function in order to determine the best model
//!
//! Like RANSAC, it can take advantage of better sampling and termination methodologies, but that's
//! separate from the method itself.
//!
//! The caveat is that MAGSAC++ still requires a `σ_max` parameter...representing the absolute
//! upper limit of expected noise...which feels a lot like we've come back to the requirement for an
//! inlier threshold. There's good news and bad news here; the good news is that MAGSAC++'s
//! performance is less sensitive to the value of `σ_max` than RANSAC is to the inlier threshold,
//! so long as you don't pick a `σ_max` that is too small (in which case they're about equally bad).
//! The bad news is that MAGSAC++'s performance will still degrade as `σ_max` gets too large for
//! the actual amount of noise present in the data, first by a loss of accuracy in the output, then
//! by a loss of robustness (it will start picking worse models).
//!
//! That said, the `σ_max` requirement in MAGSAC++ is no worse than the threshold requirement for
//! RANSAC, and it seems promising that the quality of the models produced by it will be superior
//! to RANSAC, so it seems worth implementing.
//!
//! # Geometric Entities
//!
//! The goal of this module is to create a unified framework that can be used to implement
//! consensus creation methods for the following geometric entities.
//!
//! Two-dimensional entities:
//!   - `Line2`
//!   - `Circle2`
//!   - `CubicSpline2`
//!
//! Three-dimensional entities:
//!   - `Line3`
//!   - `Circle3`
//!   - `CubicSpline3`
//!   - `Plane3`
//!   - `Sphere3`
//!   - `Cylinder3`
//!   - `Cone3`
//!
//! Fitting will need to be against true geometric residuals, for which we'll use the
//! Levenberg-Marquardt solver this library already relies on.  There's some amount of
//! generalization that can be done here.
//!
//! # API
//!
//! The framework has three pieces:
//!   - [`ConsensusModel`], a dimension-generic trait implemented by each geometric primitive. It
//!     knows how to build a candidate from a minimal random sample, evaluate a point's geometric
//!     residual, and perform one weighted least-squares refinement of itself.
//!   - [`Magsac`], the estimator configuration and driver. Its [`Magsac::fit`] /
//!     [`Magsac::fit_filtered`] methods run the sample-and-refine loop over any `ConsensusModel`.
//!   - [`ConsensusFit`], the result, holding the estimated model, the inlier indices, and the
//!     model's marginalized quality score.
//!
//! Primitive-specific trait implementations live next to their types (for example the `Circle2`
//! implementation is in `geom2/circle2/consensus.rs`), reusing that primitive's existing
//! constructors, distance methods, and analytic-jacobian Levenberg-Marquardt problems.

mod weights;

use crate::Result;
use crate::na::Point;
use rand::SeedableRng;
use rand::distr::{Distribution, Uniform};
use rand::prelude::StdRng;
use weights::MagsacWeight;

/// A geometric primitive that can be estimated from a set of points by the MAGSAC++ consensus
/// driver in [`Magsac`].
///
/// Implementors describe three things: how to build a candidate model from a minimal random sample,
/// how far an arbitrary point lies from the model, and how to perform a single weighted
/// least-squares refinement from a current estimate. The driver composes these into the full
/// marginalizing-sample-consensus loop.
pub trait ConsensusModel<const D: usize>: Sized {
    /// The minimum number of points needed to instantiate a candidate model (for example, three
    /// for a circle or two for a line).
    const SAMPLE_SIZE: usize;

    /// The degrees of freedom of the residual distribution, used by the MAGSAC++ weighting. For
    /// every primitive here the residual is a Euclidean point-to-model distance in `D`-dimensional
    /// space, so this is the ambient dimension `D` (which must be at least 2).
    const RESIDUAL_DOF: usize = D;

    /// Build a candidate model from a minimal sample of exactly [`Self::SAMPLE_SIZE`] points.
    /// Returns `None` if the sample is degenerate (for example, collinear points for a circle).
    fn from_sample(sample: &[Point<f64, D>]) -> Option<Self>;

    /// The geometric residual (signed or unsigned distance) from `point` to the model surface. Only
    /// its magnitude is used for scoring; the sign, where meaningful, keeps the residual smooth for
    /// the least-squares refinement.
    fn residual(&self, point: &Point<f64, D>) -> f64;

    /// Perform one weighted least-squares refit from `initial`, weighting each point by the matching
    /// entry of `weights`. Implemented with an analytic-jacobian Levenberg-Marquardt problem
    /// specific to the primitive. Returns `None` if the solve fails.
    fn refine_weighted(points: &[Point<f64, D>], weights: &[f64], initial: &Self) -> Option<Self>;
}

/// Configuration and driver for MAGSAC++ robust estimation of a [`ConsensusModel`].
///
/// Construct one with [`Magsac::new`], adjust any fields, then call [`Magsac::fit`] or
/// [`Magsac::fit_filtered`].
#[derive(Debug, Clone)]
pub struct Magsac {
    /// The upper bound on the expected inlier noise. This replaces RANSAC's hard inlier threshold;
    /// the estimator is far less sensitive to it, so long as it is not chosen smaller than the
    /// actual noise.
    pub sigma_max: f64,

    /// The maximum number of minimal-sample iterations. Defaults to 500 when `None`. The loop can
    /// terminate earlier once the adaptive confidence criterion is met.
    pub max_iterations: Option<usize>,

    /// The number of iteratively reweighted least-squares refinement steps applied to each
    /// candidate model (the σ-consensus local optimization).
    pub refinement_steps: usize,

    /// The probability that at least one all-inlier sample is drawn, used for adaptive termination.
    pub confidence: f64,

    /// An optional fixed RNG seed for reproducible sampling. When `None`, a random seed is used.
    pub seed: Option<u64>,
}

impl Magsac {
    /// Create a configuration with the given noise bound and sensible defaults (500 iterations, 4
    /// refinement steps, 0.99 confidence, random seed).
    pub fn new(sigma_max: f64) -> Self {
        Self {
            sigma_max,
            max_iterations: None,
            refinement_steps: 4,
            confidence: 0.99,
            seed: None,
        }
    }

    /// Estimate a model from `points`, accepting any candidate the primitive can build.
    pub fn fit<const D: usize, M: ConsensusModel<D>>(
        &self,
        points: &[Point<f64, D>],
    ) -> Result<ConsensusFit<M>> {
        self.fit_filtered(points, |_| true)
    }

    /// Estimate a model from `points`, rejecting any candidate for which `accept` returns `false`.
    /// This is how caller constraints (such as radius limits on a circle) are enforced during both
    /// sampling and refinement.
    pub fn fit_filtered<const D: usize, M, F>(
        &self,
        points: &[Point<f64, D>],
        accept: F,
    ) -> Result<ConsensusFit<M>>
    where
        M: ConsensusModel<D>,
        F: Fn(&M) -> bool,
    {
        if self.sigma_max <= 0.0 || !self.sigma_max.is_finite() {
            return Err("sigma_max must be a positive, finite value".into());
        }
        let n = points.len();
        if n < M::SAMPLE_SIZE {
            return Err("Not enough points to draw a minimal sample".into());
        }

        let weighting = MagsacWeight::new(self.sigma_max, M::RESIDUAL_DOF);
        let cutoff = weighting.cutoff();

        let seed = self.seed.unwrap_or_else(rand::random);
        let mut rng = StdRng::seed_from_u64(seed);
        let sampler = Uniform::new(0, n)?;

        let hard_cap = self.max_iterations.unwrap_or(500).max(1);
        let mut dynamic_cap = hard_cap;

        let mut sample = Vec::with_capacity(M::SAMPLE_SIZE);
        let mut weight_buf = vec![0.0; n];
        let mut best: Option<(M, f64)> = None;

        let mut iter = 0;
        while iter < hard_cap.min(dynamic_cap) {
            iter += 1;

            if !draw_sample(points, M::SAMPLE_SIZE, &sampler, &mut rng, &mut sample) {
                continue;
            }
            let Some(mut model) = M::from_sample(&sample) else {
                continue;
            };
            if !accept(&model) {
                continue;
            }

            // σ-consensus local optimization: iteratively reweight and refit.
            refine(
                &mut model,
                points,
                &weighting,
                &mut weight_buf,
                self.refinement_steps,
                &accept,
            );
            let score = fill_weights(&model, points, &weighting, &mut weight_buf);

            let improved = match &best {
                None => true,
                Some((_, best_score)) => score > *best_score,
            };
            if improved {
                // Update the adaptive iteration budget from the new best inlier ratio.
                let inliers = count_inliers(&model, points, cutoff);
                let ratio = inliers as f64 / n as f64;
                dynamic_cap = required_iterations(ratio, M::SAMPLE_SIZE, self.confidence);
                best = Some((model, score));
            }
        }

        let (mut model, _) = best.ok_or("Consensus estimation failed to find any valid model")?;

        // Final polish over all points, then collect the inlier set and final score.
        refine(
            &mut model,
            points,
            &weighting,
            &mut weight_buf,
            self.refinement_steps,
            &accept,
        );
        let score = fill_weights(&model, points, &weighting, &mut weight_buf);
        let inliers = (0..n)
            .filter(|&i| model.residual(&points[i]).abs() < cutoff)
            .collect();

        Ok(ConsensusFit {
            model,
            inliers,
            score,
        })
    }
}

/// The result of a successful [`Magsac`] estimation.
#[derive(Debug, Clone)]
pub struct ConsensusFit<M> {
    /// The estimated model.
    pub model: M,

    /// The indices (into the input slice) of the points classified as inliers, i.e. whose residual
    /// magnitude is below the MAGSAC++ cutoff `k · σ_max`.
    pub inliers: Vec<usize>,

    /// The marginalized quality score of the model (the sum of the MAGSAC++ point weights). Larger
    /// is better; it is only meaningful for comparing models fit with the same configuration.
    pub score: f64,
}

/// Draw `size` distinct point indices uniformly and collect the corresponding points into `out`.
/// Returns `false` if distinct indices could not be drawn in a reasonable number of attempts.
fn draw_sample<const D: usize>(
    points: &[Point<f64, D>],
    size: usize,
    sampler: &Uniform<usize>,
    rng: &mut StdRng,
    out: &mut Vec<Point<f64, D>>,
) -> bool {
    out.clear();
    let mut chosen: [usize; 8] = [usize::MAX; 8];
    let mut count = 0;
    let mut attempts = 0;
    while count < size {
        let i = sampler.sample(rng);
        if !chosen[..count].contains(&i) {
            chosen[count] = i;
            out.push(points[i]);
            count += 1;
        }
        attempts += 1;
        if attempts > size * 32 {
            return false;
        }
    }
    true
}

/// Run up to `steps` iteratively reweighted refinements of `model` in place, stopping early if a
/// refit fails or produces a rejected model.
fn refine<const D: usize, M: ConsensusModel<D>>(
    model: &mut M,
    points: &[Point<f64, D>],
    weighting: &MagsacWeight,
    weight_buf: &mut [f64],
    steps: usize,
    accept: &impl Fn(&M) -> bool,
) {
    for _ in 0..steps {
        fill_weights(model, points, weighting, weight_buf);
        match M::refine_weighted(points, weight_buf, model) {
            Some(next) if accept(&next) => *model = next,
            _ => break,
        }
    }
}

/// Fill `weight_buf` with the MAGSAC++ weight of every point under `model` and return their sum
/// (the model's quality score).
fn fill_weights<const D: usize, M: ConsensusModel<D>>(
    model: &M,
    points: &[Point<f64, D>],
    weighting: &MagsacWeight,
    weight_buf: &mut [f64],
) -> f64 {
    let mut score = 0.0;
    for (i, p) in points.iter().enumerate() {
        let r = model.residual(p).abs();
        let w = weighting.weight(r);
        weight_buf[i] = w;
        score += w;
    }
    score
}

/// Count points whose residual magnitude falls below the MAGSAC++ cutoff.
fn count_inliers<const D: usize, M: ConsensusModel<D>>(
    model: &M,
    points: &[Point<f64, D>],
    cutoff: f64,
) -> usize {
    points
        .iter()
        .filter(|p| model.residual(p).abs() < cutoff)
        .count()
}

/// The number of RANSAC iterations needed to draw an all-inlier sample with the given confidence,
/// for an inlier ratio `w` and minimal sample size `m`: `log(1 − confidence) / log(1 − wᵐ)`.
fn required_iterations(inlier_ratio: f64, sample_size: usize, confidence: f64) -> usize {
    if inlier_ratio <= 0.0 {
        return usize::MAX;
    }
    let w_m = inlier_ratio.powi(sample_size as i32);
    if w_m >= 1.0 {
        return 1;
    }
    let num = (1.0 - confidence).ln();
    let den = (1.0 - w_m).ln();
    (num / den).ceil().max(1.0) as usize
}
