//! Circle fitting for [`Circle2`], both ordinary least squares and MAGSAC++ robust consensus.
//!
//! Both entry points share a single refinement engine, [`CircleFit`], a weighted least-squares
//! problem over a circle's three degrees of freedom (center and radius) solved with
//! Levenberg-Marquardt against the true geometric radial residuals:
//!   - [`Circle2::from_fit`] seeds it with a closed-form algebraic (Kåsa-style) estimate and refines
//!     once. It is fast but not robust to gross outliers.
//!   - [`Circle2::from_consensus`] drives it through the MAGSAC++ sample-and-refine loop, where each
//!     minimal sample is three points (via [`Circle2::from_3_points`]) and the same [`CircleFit`]
//!     performs the iteratively reweighted σ-consensus refinement of every candidate.

use super::Circle2;
use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::{Point2, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry2d_f64::na::{Dyn, Matrix, Matrix3, Owned, U1, U3, Vector, Vector3};

impl ConsensusModel<2> for Circle2 {
    type Point = Point2;
    const SAMPLE_SIZE: usize = 3;

    fn from_sample(sample: &[Point2]) -> Option<Self> {
        Circle2::from_3_points(&sample[0], &sample[1], &sample[2]).ok()
    }

    fn residual(&self, point: &Point2) -> f64 {
        self.distance_to(point)
    }

    fn refine_weighted(points: &[Point2], weights: &[f64], initial: &Circle2) -> Option<Circle2> {
        CircleFit::refine(points, weights, initial)
    }
}

/// Closed-form algebraic circle fit (Kåsa-style) by weighted linear least squares over already
/// converted points and weights. This is the shared core of [`Circle2::from_fit_algebraic`].
///
/// Writing the circle equation as `x² + y² = a·x + b·y + c`, each point contributes a linear
/// equation in the unknowns `s = [a, b, c]`, where `a = 2cx`, `b = 2cy`, and `c = r² − cx² − cy²`.
/// Solving the 3×3 weighted normal equations gives the center directly and the radius from
/// `r² = c + cx² + cy²`. Returns `None` if the normal matrix is singular (fewer than three points,
/// or collinear inputs).
fn algebraic_circle_fit(points: &[Point2], weights: &[f64]) -> Option<Circle2> {
    let mut m = Matrix3::<f64>::zeros();
    let mut v = Vector3::<f64>::zeros();
    for (p, &w) in points.iter().zip(weights) {
        let row = Vector3::new(p.x, p.y, 1.0);
        let target = p.x * p.x + p.y * p.y;
        m += w * row * row.transpose();
        v += w * target * row;
    }

    let s = m.lu().solve(&v)?;
    let cx = 0.5 * s[0];
    let cy = 0.5 * s[1];
    let r_squared = s[2] + cx * cx + cy * cy;
    if r_squared <= 0.0 {
        return None;
    }
    Some(Circle2::new(cx, cy, r_squared.sqrt()))
}

impl Circle2 {
    /// Fit a circle to a set of points using only the closed-form algebraic (Kåsa-style) weighted
    /// least-squares method, without the geometric Levenberg-Marquardt refinement that
    /// [`Circle2::from_fit`] applies on top of it.
    ///
    /// Writing the circle equation as `x² + y² = a·x + b·y + c`, each point contributes a linear
    /// equation in the unknowns `s = [a, b, c]` (with `a = 2cx`, `b = 2cy`, `c = r² − cx² − cy²`);
    /// solving the 3×3 weighted normal equations gives the center and radius directly. The algebraic
    /// error minimized is not the true geometric radial distance, so the result is slightly biased
    /// (most noticeably for partial arcs), but it is fast, closed-form, and needs no initial guess,
    /// which makes it the seed for [`Circle2::from_fit`].
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of at least three non-collinear coordinates to fit the circle to
    /// * `weights`: if `Some`, a slice the same length as `points` giving the weight to multiply each
    ///   point's contribution by; if `None`, all points are weighted equally.
    ///
    /// returns: Result<Circle2, Box<dyn Error, Global>>
    pub fn from_fit_algebraic(points: &[impl PCoords<2>], weights: Option<&[f64]>) -> Result<Self> {
        let pts: Vec<Point2> = points.iter().map(|p| Point2::from(p.coords())).collect();
        if pts.len() < 3 {
            return Err("At least three points are required to fit a circle".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; pts.len()];
                &ones
            }
        };

        algebraic_circle_fit(&pts, weights)
            .ok_or_else(|| "Failed to fit circle: points are collinear or degenerate".into())
    }

    /// Fit a circle to a set of points by ordinary least squares. A closed-form algebraic
    /// (Kåsa-style) estimate from [`Circle2::from_fit_algebraic`] provides the initial guess, which
    /// is then refined against the true geometric radial residuals with the same weighted
    /// [`CircleFit`] Levenberg-Marquardt engine used by the consensus fit. Optional weights may be
    /// provided in a slice of `f64` with the same number of elements as `points`, where the weight
    /// `i` corresponds with the point `i`.
    ///
    /// This is not robust to gross outliers; for that, use [`Circle2::from_consensus`].
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of at least three non-collinear coordinates to fit the circle to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Circle2, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// use approx::assert_relative_eq;
    ///
    /// let points = vec![
    ///     Point2::new(-1.0, 0.0),
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(0.0, -1.0),
    /// ];
    ///
    /// let circle = Circle2::from_fit(&points, None).unwrap();
    /// assert_relative_eq!(circle.x(), 0.0);
    /// assert_relative_eq!(circle.y(), 0.0);
    /// assert_relative_eq!(circle.r(), 1.0);
    /// ```
    pub fn from_fit(points: &[impl PCoords<2>], weights: Option<&[f64]>) -> Result<Self> {
        // The algebraic fit validates the point count and produces the initial guess.
        let guess = Circle2::from_fit_algebraic(points, weights)?;

        let pts: Vec<Point2> = points.iter().map(|p| Point2::from(p.coords())).collect();
        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; pts.len()];
                &ones
            }
        };

        // Refine the algebraic guess against true geometric residuals with the weighted LM engine.
        CircleFit::refine(&pts, weights, &guess).ok_or_else(|| "Failed to refine circle fit".into())
    }

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

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

/// A weighted least-squares problem for refining a circle, shared by both [`Circle2::from_fit`] and
/// the consensus [`ConsensusModel`] implementation. The three parameters are the center coordinates
/// (`[0..2]`) and radius (`[2]`); the jacobian of the weighted point-to-perimeter distances is
/// analytic.
struct CircleFit<'a> {
    /// The points to be fit to the circle.
    points: &'a [Point2],

    /// The parameters being fit
    x: Vector3<f64>,

    /// The current active circle
    circle: Circle2,

    /// The active base residuals
    base_residuals: Residuals,

    /// The per-point weights, held fixed across the solve
    weights: Residuals,
}

impl<'a> CircleFit<'a> {
    /// Create a circle fit with fixed, externally supplied per-point weights. `weights` must have
    /// the same length as `points`.
    fn new(points: &'a [Point2], weights: &[f64], initial: &Circle2) -> Self {
        let x = Vector3::new(initial.center.x, initial.center.y, initial.r());
        let circle = *initial;

        let mut base_residuals = Residuals::zeros(points.len());
        compute_residuals_mut(points, &circle, &mut base_residuals);

        Self {
            points,
            x,
            circle,
            base_residuals,
            weights: Residuals::from_column_slice(weights),
        }
    }

    /// Refine `initial` against the weighted geometric radial residuals of `points` with a single
    /// Levenberg-Marquardt solve, returning the optimized circle or `None` if the solve fails. This
    /// is the shared entry point for both the least-squares and consensus fits.
    fn refine(points: &[Point2], weights: &[f64], initial: &Circle2) -> Option<Circle2> {
        let problem = CircleFit::new(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then_some(result.circle)
    }
}

fn compute_residuals_mut(points: &[Point2], circle: &Circle2, residuals: &mut Residuals) {
    for (i, p) in points.iter().enumerate() {
        residuals[i] = circle.distance_to(p)
    }
}

impl LeastSquaresProblem<f64, Dyn, U3> for CircleFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U3>;
    type ParameterStorage = Owned<f64, U3>;

    fn set_params(&mut self, x: &Vector<f64, U3, Self::ParameterStorage>) {
        self.x = *x;
        self.circle = Circle2::new(x[0], x[1], x[2]);
        compute_residuals_mut(self.points, &self.circle, &mut self.base_residuals);
    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        self.x
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Residuals::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.base_residuals[i] * self.weights[i];
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U3, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U3, Self::JacobianStorage>::zeros(self.points.len());
        for (i, p) in self.points.iter().enumerate() {
            // Find the vector from the center of the circle to the point
            let v = p - self.circle.center;

            // Normalize it
            let n = v.normalize();

            // Fill in the jacobian for this row
            jac[(i, 0)] = -n.x * self.weights[i];
            jac[(i, 1)] = -n.y * self.weights[i];
            jac[(i, 2)] = -self.weights[i];
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use crate::common::consensus::Magsac;
    use crate::common::random_geometry::RandomGeometry2;
    use crate::{Circle2, Point2, Result};
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    /// Generate `n` points around a circle's perimeter with independent radial Gaussian noise.
    fn circle_noise(rg: &mut RandomGeometry2, c: &Circle2, n: usize, sigma: f64) -> Vec<Point2> {
        (0..n)
            .map(|i| {
                let t = TAU * i as f64 / n as f64;
                let r = c.r() + rg.gaussian_f64(0.0, sigma);
                Point2::new(c.x() + r * t.cos(), c.y() + r * t.sin())
            })
            .collect()
    }

    // from_fit tests ─────────────────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_clean_circle() -> Result<()> {
        let mut rg = RandomGeometry2::from_seed(3);
        let expected = Circle2::new(2.0, 3.0, 1.0);
        let points = circle_noise(&mut rg, &expected, 60, 0.0);

        let fit = Circle2::from_fit(&points, None)?;
        assert_relative_eq!(fit.center, expected.center, epsilon = 1e-9);
        assert_relative_eq!(fit.r(), expected.r(), epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn from_fit_with_noise_is_close() -> Result<()> {
        let mut rg = RandomGeometry2::from_seed(8);
        let expected = Circle2::new(2.0, 3.0, 1.0);
        let points = circle_noise(&mut rg, &expected, 500, 0.01);

        let fit = Circle2::from_fit(&points, None)?;
        assert_relative_eq!(fit.center, expected.center, epsilon = 3e-3);
        assert_relative_eq!(fit.r(), expected.r(), epsilon = 3e-3);
        Ok(())
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() -> Result<()> {
        let mut rg = RandomGeometry2::from_seed(21);
        let expected = Circle2::new(-1.0, 0.5, 2.0);
        let points = circle_noise(&mut rg, &expected, 40, 0.02);

        let unweighted = Circle2::from_fit(&points, None)?;
        let weights = vec![1.0; points.len()];
        let weighted = Circle2::from_fit(&points, Some(&weights))?;
        assert_relative_eq!(weighted.center, unweighted.center, epsilon = 1e-9);
        assert_relative_eq!(weighted.r(), unweighted.r(), epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn from_fit_too_few_points_is_error() {
        let points = [Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
        assert!(Circle2::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_collinear_points_is_error() {
        let points = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
        ];
        assert!(Circle2::from_fit(&points, None).is_err());
    }

    // from_consensus tests ───────────────────────────────────────────────────

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
