//! Sphere fitting for [`Sphere3`], both ordinary least squares and MAGSAC++ robust consensus.
//!
//! Both entry points share a single refinement engine, [`Sphere3Fit`], a weighted least-squares
//! problem over a sphere's four degrees of freedom (center and radius) solved with
//! Levenberg-Marquardt against the true geometric radial residuals:
//!   - [`Sphere3::from_fit`] seeds it with a closed-form algebraic (Kåsa-style) estimate and
//!     refines once. It is fast but not robust to gross outliers.
//!   - [`Sphere3::from_consensus`] drives it through the MAGSAC++ sample-and-refine loop, where each
//!     minimal sample is four points (via [`Sphere3::from_4_points`]) and the same [`Sphere3Fit`]
//!     performs the iteratively reweighted σ-consensus refinement of every candidate.

use super::Sphere3;
use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::{Point3, Result, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Matrix4, Owned, U1, U4, Vector, Vector4};

/// The signed radial distance from a point to a sphere's surface: `|point − center| − radius`. It
/// is positive outside the sphere and negative inside; only its magnitude is used for scoring, but
/// keeping the sign makes the residual smooth through the surface for the least-squares refinement.
fn point_sphere_distance(center: &Point3, radius: f64, point: &Point3) -> f64 {
    (point - center).norm() - radius
}

/// Closed-form algebraic sphere fit (Kåsa-style) by weighted linear least squares.
///
/// Writing the sphere equation as `x² + y² + z² = a·x + b·y + c·z + d`, each point contributes a
/// linear equation in the unknowns `s = [a, b, c, d]`, where `a = 2cx`, `b = 2cy`, `c = 2cz`, and
/// `d = r² − |center|²`. Solving the 4×4 weighted normal equations gives the center directly and
/// the radius from `r² = d + |center|²`. Returns `None` if the normal matrix is singular (fewer
/// than four points, or coplanar/collinear inputs).
fn algebraic_sphere_fit(points: &[Point3], weights: &[f64]) -> Option<Sphere3> {
    let mut m = Matrix4::zeros();
    let mut v = Vector4::zeros();
    for (p, &w) in points.iter().zip(weights) {
        let row = Vector4::new(p.x, p.y, p.z, 1.0);
        let target = p.coords.norm_squared();
        m += w * row * row.transpose();
        v += w * target * row;
    }

    let s = m.lu().solve(&v)?;
    let center = Point3::new(0.5 * s[0], 0.5 * s[1], 0.5 * s[2]);
    let r_squared = s[3] + center.coords.norm_squared();
    if r_squared <= 0.0 {
        return None;
    }
    Some(Sphere3::new(&center, r_squared.sqrt()))
}

impl Sphere3 {
    /// Fit a sphere to a set of points by ordinary least squares. A closed-form algebraic
    /// (Kåsa-style) estimate provides the initial guess, which is then refined against the true
    /// geometric radial residuals with the same weighted [`Sphere3Fit`] Levenberg-Marquardt engine
    /// used by the consensus fit. Optional weights may be provided in a slice of `f64` with the same
    /// number of elements as `points`, where the weight `i` corresponds with the point `i`.
    ///
    /// This is not robust to gross outliers; for that, use [`Sphere3::from_consensus`].
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of at least four non-coplanar coordinates to fit the sphere to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Sphere3, Box<dyn Error, Global>>
    pub fn from_fit(points: &[impl PCoords<3>], weights: Option<&[f64]>) -> Result<Self> {
        let pts: Vec<Point3> = points.iter().map(|p| Point3::from(p.coords())).collect();
        if pts.len() < 4 {
            return Err("At least four points are required to fit a sphere".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; pts.len()];
                &ones
            }
        };

        let guess = algebraic_sphere_fit(&pts, weights)
            .ok_or("Failed to fit sphere: points are coplanar or degenerate")?;

        // Refine the algebraic guess against true geometric residuals with the weighted LM engine.
        Sphere3Fit::refine(&pts, weights, &guess)
            .ok_or_else(|| "Failed to refine sphere fit".into())
    }

    /// Fit a sphere to a set of 3D points using MAGSAC++ robust consensus estimation.
    ///
    /// Unlike an ordinary least-squares fit ([`Sphere3::from_fit`]), this rejects gross outliers by
    /// taking an upper bound on the inlier noise (`sigma_max`) rather than a hard inlier/outlier
    /// threshold, and refines each candidate with noise-marginalized iteratively reweighted least
    /// squares. It is substantially less sensitive to `sigma_max` than RANSAC is to its threshold,
    /// as long as `sigma_max` is not chosen smaller than the actual noise.
    ///
    /// # Arguments
    ///
    /// * `points`: the 3D points to fit the sphere to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `min_r`: an optional minimum radius; candidate spheres smaller than this are rejected
    /// * `max_r`: an optional maximum radius; candidate spheres larger than this are rejected
    /// * `options`: an optional [`Magsac`] configuration; its `sigma_max` is overridden by the
    ///   `sigma_max` argument
    ///
    /// returns: Result<Sphere3, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point3],
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Sphere3> {
        let min_r = min_r.unwrap_or(0.0);
        let max_r = max_r.unwrap_or(f64::INFINITY);

        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit =
            magsac.fit_filtered::<3, Sphere3, _>(points, |s| s.r() >= min_r && s.r() <= max_r)?;
        Ok(fit.model)
    }
}

impl ConsensusModel<3> for Sphere3 {
    type Point = Point3;
    const SAMPLE_SIZE: usize = 4;

    fn from_sample(sample: &[Point3]) -> Option<Self> {
        // `from_4_points` already rejects coplanar (and coincident) samples.
        Sphere3::from_4_points(&sample[0], &sample[1], &sample[2], &sample[3]).ok()
    }

    fn residual(&self, point: &Point3) -> f64 {
        point_sphere_distance(&self.center, self.radius, point)
    }

    fn refine_weighted(points: &[Point3], weights: &[f64], initial: &Sphere3) -> Option<Sphere3> {
        Sphere3Fit::refine(points, weights, initial)
    }
}

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

/// A weighted least-squares problem for refining a sphere, shared by both [`Sphere3::from_fit`] and
/// the consensus [`ConsensusModel`] implementation. The four parameters are the absolute center
/// coordinates (`[0..3]`) and radius (`[3]`); a sphere has no orientation, so no relative
/// parameterization is needed. The jacobian of the weighted point-to-surface distances is analytic.
struct Sphere3Fit<'a> {
    points: &'a [Point3],
    weights: &'a [f64],
    params: Vector4<f64>,
    residuals: Residuals,
}

impl<'a> Sphere3Fit<'a> {
    fn new(points: &'a [Point3], weights: &'a [f64], base: &Sphere3) -> Self {
        let params = Vector4::new(base.center.x, base.center.y, base.center.z, base.radius);
        let mut problem = Self {
            points,
            weights,
            params,
            residuals: Residuals::zeros(points.len()),
        };
        problem.set_params(&params);
        problem
    }

    /// Refine `initial` against the weighted geometric radial residuals of `points` with a single
    /// Levenberg-Marquardt solve, returning the optimized sphere or `None` if the solve fails. This
    /// is the shared entry point for both the least-squares and consensus fits.
    fn refine(points: &[Point3], weights: &[f64], initial: &Sphere3) -> Option<Sphere3> {
        let problem = Sphere3Fit::new(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then(|| {
            let (center, radius) = result.sphere_params(&result.params);
            Sphere3::new(&center, radius)
        })
    }

    /// Reconstruct the `(center, radius)` described by a parameter vector, clamping the radius to a
    /// small positive floor.
    fn sphere_params(&self, p: &Vector4<f64>) -> (Point3, f64) {
        (Point3::new(p[0], p[1], p[2]), p[3].max(1e-9))
    }
}

impl LeastSquaresProblem<f64, Dyn, U4> for Sphere3Fit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U4>;
    type ParameterStorage = Owned<f64, U4>;

    fn set_params(&mut self, x: &Vector<f64, U4, Self::ParameterStorage>) {
        self.params = *x;
        let (center, radius) = self.sphere_params(&self.params);
        for i in 0..self.points.len() {
            self.residuals[i] =
                self.weights[i] * point_sphere_distance(&center, radius, &self.points[i]);
        }
    }

    fn params(&self) -> Vector<f64, U4, Self::ParameterStorage> {
        self.params
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        Some(self.residuals.clone())
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U4, Self::JacobianStorage>> {
        let (center, _radius) = self.sphere_params(&self.params);
        let mut jac = Matrix::<f64, Dyn, U4, Self::JacobianStorage>::zeros(self.points.len());

        for (i, point) in self.points.iter().enumerate() {
            let v = point - center;
            let rho = v.norm();
            // Unit vector from the center toward the point; the derivative of the radial distance
            // with respect to the center is its negation. Degenerate at the center (rho ≈ 0), where
            // the gradient is undefined; use zero there.
            let u = if rho > 1e-12 {
                v / rho
            } else {
                Vector3::zeros()
            };
            let w = self.weights[i];
            jac[(i, 0)] = -w * u.x;
            jac[(i, 1)] = -w * u.y;
            jac[(i, 2)] = -w * u.z;
            jac[(i, 3)] = -w;
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use crate::common::consensus::{ConsensusModel, Magsac};
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom3::Sphere3;
    use crate::{Point3, Result};
    use approx::assert_relative_eq;

    fn offset() -> Sphere3 {
        Sphere3::new(&Point3::new(1.0, -2.0, 3.0), 2.5)
    }

    /// Generate `n` points scattered over a sphere's surface with isotropic Gaussian noise `sigma`.
    fn sphere_noise(
        rg: &mut RandomGeometry3,
        sphere: &Sphere3,
        n: usize,
        sigma: f64,
    ) -> Vec<Point3> {
        (0..n)
            .map(|_| {
                let dir = rg.unit_vec();
                sphere.center + dir.into_inner() * sphere.r() + rg.gaussian_vector::<3>(sigma)
            })
            .collect()
    }

    fn assert_matches(result: &Sphere3, expected: &Sphere3) {
        assert_relative_eq!(result.center, expected.center, epsilon = 5.0e-3);
        assert_relative_eq!(result.r(), expected.r(), epsilon = 5.0e-3);
    }

    // from_fit tests ─────────────────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_clean_sphere() {
        let mut rg = RandomGeometry3::from_seed(3);
        let expected = offset();
        let points = sphere_noise(&mut rg, &expected, 50, 0.0);

        let fit = Sphere3::from_fit(&points, None).unwrap();
        assert_relative_eq!(fit.center, expected.center, epsilon = 1e-9);
        assert_relative_eq!(fit.r(), expected.r(), epsilon = 1e-9);
        for p in &points {
            assert_relative_eq!((p - fit.center).norm(), fit.r(), epsilon = 1e-8);
        }
    }

    #[test]
    fn from_fit_with_noise_is_close() {
        let mut rg = RandomGeometry3::from_seed(8);
        let expected = Sphere3::new(&Point3::new(-1.0, 4.0, 2.0), 3.0);
        let points = sphere_noise(&mut rg, &expected, 400, 0.01);

        let fit = Sphere3::from_fit(&points, None).unwrap();
        assert_relative_eq!(fit.center, expected.center, epsilon = 5e-3);
        assert_relative_eq!(fit.r(), expected.r(), epsilon = 5e-3);
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() {
        let mut rg = RandomGeometry3::from_seed(21);
        let expected = offset();
        let points = sphere_noise(&mut rg, &expected, 40, 0.02);

        let unweighted = Sphere3::from_fit(&points, None).unwrap();
        let weights = vec![1.0; points.len()];
        let weighted = Sphere3::from_fit(&points, Some(&weights)).unwrap();
        assert_relative_eq!(weighted.center, unweighted.center, epsilon = 1e-9);
        assert_relative_eq!(weighted.r(), unweighted.r(), epsilon = 1e-9);
    }

    #[test]
    fn from_fit_too_few_points_is_error() {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        assert!(Sphere3::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_coplanar_points_is_error() {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
        ];
        assert!(Sphere3::from_fit(&points, None).is_err());
    }

    // from_consensus tests ───────────────────────────────────────────────────

    #[test]
    fn from_consensus_recovers_clean_sphere() -> Result<()> {
        let mut rg = RandomGeometry3::from_seed(11);
        let expected = offset();
        let points = sphere_noise(&mut rg, &expected, 120, 0.0);

        let result = Sphere3::from_consensus(&points, 0.01, None, None, None)?;
        assert_matches(&result, &expected);
        Ok(())
    }

    #[test]
    fn from_consensus_rejects_outliers() -> Result<()> {
        let mut rg = RandomGeometry3::from_seed(29);
        let expected = offset();
        let inliers = sphere_noise(&mut rg, &expected, 150, 0.01);
        let inlier_count = inliers.len();
        let mut points = inliers.clone();

        // A dense cluster of gross outliers, well off the sphere's surface.
        let center = expected.center + crate::Vector3::new(6.0, -5.0, 4.0);
        for _ in 0..50 {
            points.push(center + rg.gaussian_vector::<3>(1.0));
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(500),
            refinement_steps: 5,
            confidence: 0.99,
            seed: Some(17),
        };
        let fit = magsac.fit_filtered::<3, Sphere3, _>(&points, |_| true)?;

        assert_matches(&fit.model, &expected);
        // Every inlier should lie very close to the recovered sphere.
        for i in &inliers {
            assert!(fit.model.residual(i).abs() < 0.01 * 6.0);
        }
        // No outlier should be classified as an inlier.
        assert!(fit.inliers.iter().all(|&i| i < inlier_count));
        assert!(fit.inliers.len() > inlier_count * 9 / 10);
        Ok(())
    }
}
