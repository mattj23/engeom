//! MAGSAC++ consensus fitting for [`Circle2`].
//!
//! This wires `Circle2` into the generic consensus framework in [`crate::common::consensus`]. The
//! minimal sample is three points (via [`Circle2::from_3_points`]), the residual is the signed
//! distance to the perimeter (via [`Circle2::distance_to`]), and the weighted refinement reuses the
//! analytic-jacobian [`CircleFit`] Levenberg-Marquardt problem with fixed per-point weights.

use super::{Circle2, CircleFit};
use crate::Point2;
use crate::Result;
use crate::common::consensus::{ConsensusModel, Magsac};
use levenberg_marquardt::LevenbergMarquardt;

impl ConsensusModel<2> for Circle2 {
    const SAMPLE_SIZE: usize = 3;

    fn from_sample(sample: &[Point2]) -> Option<Self> {
        Circle2::from_3_points(&sample[0], &sample[1], &sample[2]).ok()
    }

    fn residual(&self, point: &Point2) -> f64 {
        self.distance_to(point)
    }

    fn refine_weighted(points: &[Point2], weights: &[f64], initial: &Circle2) -> Option<Circle2> {
        let problem = CircleFit::with_weights(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then_some(result.circle)
    }
}

impl Circle2 {
    /// Fit a circle to a set of points using MAGSAC++ robust consensus estimation.
    ///
    /// Unlike a fixed-threshold RANSAC, this takes an upper bound on the inlier noise (`sigma_max`)
    /// rather than a hard inlier/outlier threshold, and refines each candidate with
    /// noise-marginalized iteratively reweighted least squares. It is substantially less sensitive
    /// to `sigma_max` than RANSAC is to its threshold, as long as `sigma_max` is not chosen smaller
    /// than the actual noise.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to fit the circle to
    /// * `sigma_max`: the upper bound on the expected inlier noise (in the same units as the points)
    /// * `min_r`: an optional minimum radius; candidate circles smaller than this are rejected
    /// * `max_r`: an optional maximum radius; candidate circles larger than this are rejected
    /// * `options`: an optional [`Magsac`] configuration to override iteration count, refinement
    ///   steps, confidence, or the RNG seed. Its `sigma_max` field is overridden by the `sigma_max`
    ///   argument.
    ///
    /// returns: Result<Circle2, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point2],
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Circle2> {
        let min_r = min_r.unwrap_or(0.0);
        let max_r = max_r.unwrap_or(f64::INFINITY);

        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit =
            magsac.fit_filtered::<2, Circle2, _>(points, |c| c.r() >= min_r && c.r() <= max_r)?;
        Ok(fit.model)
    }
}

#[cfg(test)]
mod tests {
    use crate::common::consensus::Magsac;
    use crate::{Circle2, Point2, Result};
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    #[test]
    fn circle_from_consensus_rejects_outliers() -> Result<()> {
        let (cx, cy, r) = (2.0, 3.0, 1.2);

        // Inliers spread around the circle with a small deterministic radial perturbation.
        let mut points = Vec::new();
        let inlier_count = 120;
        for i in 0..inlier_count {
            let t = TAU * i as f64 / inlier_count as f64;
            let rr = r + 0.003 * (7.0 * t).sin();
            points.push(Point2::new(cx + rr * t.cos(), cy + rr * t.sin()));
        }

        // A dense cluster of gross outliers well off the circle.
        for i in 0..50 {
            let t = TAU * i as f64 / 50.0;
            points.push(Point2::new(
                cx + 5.0 + 0.5 * t.cos(),
                cy - 4.0 + 0.5 * t.sin(),
            ));
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(400),
            refinement_steps: 4,
            confidence: 0.99,
            seed: Some(42),
        };
        let fit = magsac.fit_filtered::<2, Circle2, _>(&points, |_| true)?;

        assert_relative_eq!(fit.model.center.x, cx, epsilon = 5.0e-3);
        assert_relative_eq!(fit.model.center.y, cy, epsilon = 5.0e-3);
        assert_relative_eq!(fit.model.r(), r, epsilon = 5.0e-3);

        // No outlier is classified as an inlier, and nearly all true inliers are recovered.
        assert!(fit.inliers.iter().all(|&i| i < inlier_count));
        assert!(fit.inliers.len() > inlier_count * 9 / 10);

        Ok(())
    }

    #[test]
    fn from_consensus_convenience_matches() -> Result<()> {
        let (cx, cy, r) = (-1.5, 0.75, 0.9);
        let n = 80;
        let points: Vec<Point2> = (0..n)
            .map(|i| {
                let t = TAU * i as f64 / n as f64;
                Point2::new(cx + r * t.cos(), cy + r * t.sin())
            })
            .collect();

        let circle = Circle2::from_consensus(&points, 0.01, None, None, None)?;
        assert_relative_eq!(circle.center.x, cx, epsilon = 5.0e-3);
        assert_relative_eq!(circle.center.y, cy, epsilon = 5.0e-3);
        assert_relative_eq!(circle.r(), r, epsilon = 5.0e-3);
        Ok(())
    }
}
