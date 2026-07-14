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