//! Circle fitting for [`Circle3`], both ordinary least squares and MAGSAC++ robust consensus.
//!
//! All entry points share a single refinement engine, [`Circle3Fit`], a weighted least-squares
//! problem over a circle's six degrees of freedom (center, plane orientation, and radius) solved
//! with Levenberg-Marquardt against the true geometric point-to-circle distances:
//!   - [`Circle3::from_fit`] seeds it with a closed-form estimate (best-fit plane via [`SvdBasis3`]
//!     plus an in-plane algebraic circle fit) and refines once. It is fast but not robust to gross
//!     outliers.
//!   - [`Circle3::from_consensus`] drives the same [`Circle3Fit`] through the MAGSAC++
//!     sample-and-refine loop in the native 3D dimension: every candidate is a full 3D circle,
//!     residuals are true point-to-circle distances, and the refinement adjusts all six degrees of
//!     freedom, so it is robust to out-of-plane outliers at the cost of a heavier per-iteration
//!     refinement.
//!   - [`Circle3::from_consensus_planar`] instead reduces the problem to 2D by projecting the points
//!     onto the best-fit plane (via an [`SvdBasis3`]), running the [`Circle2`] consensus fit there,
//!     and lifting the result back into 3D. It is fast and effective when the points are genuinely
//!     near-planar, but its plane estimate is corrupted by gross out-of-plane outliers.

use super::Circle3;
use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::geom3::{IsoExtensions3, SvdBasis3};
use crate::{Circle2, Iso3, Point2, Point3, Result, UnitVec3, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Matrix3, Owned, U1, U6, UnitQuaternion, Vector, Vector6};

/// The Euclidean distance from a point to the nearest point on a 3D circle defined by its center,
/// unit normal, and radius. It combines the out-of-plane offset with the in-plane radial deviation.
fn point_circle_distance(center: &Point3, normal: &UnitVec3, radius: f64, point: &Point3) -> f64 {
    let v = point - center;
    let h = v.dot(&normal.into_inner());
    let in_plane = v - normal.into_inner() * h;
    let rho = in_plane.norm();
    ((rho - radius).powi(2) + h * h).sqrt()
}

impl ConsensusModel<3> for Circle3 {
    type Point = Point3;
    const SAMPLE_SIZE: usize = 3;

    fn from_sample(sample: &[Point3]) -> Option<Self> {
        Circle3::from_3_points(&sample[0], &sample[1], &sample[2]).ok()
    }

    fn residual(&self, point: &Point3) -> f64 {
        point_circle_distance(&self.center, &self.normal, self.r(), point)
    }

    fn refine_weighted(points: &[Point3], weights: &[f64], initial: &Circle3) -> Option<Circle3> {
        Circle3Fit::refine(points, weights, initial)
    }
}

/// Closed-form circle estimate for a set of 3D points, used as the initial guess for the
/// least-squares refinement. The best-fit plane is found via an [`SvdBasis3`], the points are
/// projected into that plane, and an in-plane weighted algebraic (Kåsa-style) circle fit gives the
/// center and radius, which are then lifted back into 3D with the plane's normal. Returns `None` if
/// the plane cannot be estimated (fewer than three points, or collinear inputs) or the in-plane
/// algebraic system is singular.
fn algebraic_circle_fit(points: &[Point3], weights: &[f64]) -> Option<Circle3> {
    let basis = SvdBasis3::from_points(points, Some(weights))?;

    // `to_local` maps world points into the basis frame, whose third axis is the plane normal;
    // dropping the local z coordinate projects each point onto the plane.
    let to_local = Iso3::from(&basis);
    let to_world = to_local.inverse();

    // Weighted Kåsa algebraic circle fit in the plane: for the circle `x² + y² = a·x + b·y + c`,
    // each point gives a linear equation in `s = [a, b, c]`, with `a = 2cx`, `b = 2cy`, and
    // `c = r² − cx² − cy²`.
    let mut m = Matrix3::zeros();
    let mut v = Vector3::zeros();
    for (p, &w) in points.iter().zip(weights) {
        let local = to_local.transform_point(p);
        let row = Vector3::new(local.x, local.y, 1.0);
        let target = local.x * local.x + local.y * local.y;
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

    let center = to_world.transform_point(&Point3::new(cx, cy, 0.0));
    let normal = UnitVec3::new_normalize(to_world.rotation * Vector3::z());
    Some(Circle3::new(center, normal, r_squared.sqrt()))
}

impl Circle3 {
    /// Fit a circle to a set of 3D points by ordinary least squares. A closed-form estimate (the
    /// best-fit plane via an [`SvdBasis3`] combined with an in-plane algebraic circle fit) provides
    /// the initial guess, which is then refined against the true geometric point-to-circle distances
    /// with the same weighted [`Circle3Fit`] Levenberg-Marquardt engine used by the consensus fit.
    /// Optional weights may be provided in a slice of `f64` with the same number of elements as
    /// `points`, where the weight `i` corresponds with the point `i`.
    ///
    /// This is not robust to gross outliers; for that, use [`Circle3::from_consensus`].
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of at least three non-collinear coordinates to fit the circle to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Circle3, Box<dyn Error, Global>>
    pub fn from_fit(points: &[impl PCoords<3>], weights: Option<&[f64]>) -> Result<Self> {
        let pts: Vec<Point3> = points.iter().map(|p| Point3::from(p.coords())).collect();
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

        let guess = algebraic_circle_fit(&pts, weights)
            .ok_or("Failed to fit circle: points are collinear or degenerate")?;

        // Refine the algebraic guess against true geometric residuals with the weighted LM engine.
        Circle3Fit::refine(&pts, weights, &guess)
            .ok_or_else(|| "Failed to refine circle fit".into())
    }

    /// Fit a circle to a set of 3D points by first projecting them onto their best-fit plane and
    /// running the [`Circle2`] MAGSAC++ consensus fit in 2D, then lifting the result back into 3D.
    ///
    /// This is fast and works well when the points genuinely lie on or very near a common plane,
    /// such as when working with points generated by cross-section. Because the plane is estimated
    /// from all points (including outliers) via an [`SvdBasis3`], gross out-of-plane outliers will
    /// degrade the result; use [`Circle3::from_consensus`] in that case.
    ///
    /// # Arguments
    ///
    /// * `points`: the 3D points to fit the circle to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `min_r`: an optional minimum radius; candidate circles smaller than this are rejected
    /// * `max_r`: an optional maximum radius; candidate circles larger than this are rejected
    /// * `options`: an optional [`Magsac`] configuration; its `sigma_max` is overridden by the
    ///   `sigma_max` argument
    ///
    /// returns: Result<Circle3, Box<dyn Error, Global>>
    pub fn from_consensus_planar(
        points: &[Point3],
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Circle3> {
        let basis = SvdBasis3::from_points(points, None)
            .ok_or("Failed to compute an SVD basis from the given points")?;

        // `to_local` maps world points into the basis frame, where the third axis is the plane
        // normal; dropping the local z coordinate projects each point onto the plane.
        let to_local = Iso3::from(&basis);
        let to_world = to_local.inverse();

        let points_2d: Vec<Point2> = points
            .iter()
            .map(|p| {
                let local = to_local.transform_point(p);
                Point2::new(local.x, local.y)
            })
            .collect();

        let circle = Circle2::from_consensus(&points_2d, sigma_max, min_r, max_r, options)?;

        let center = to_world.transform_point(&Point3::new(circle.x(), circle.y(), 0.0));
        let normal = UnitVec3::new_normalize(to_world.rotation * Vector3::z());
        Ok(Circle3::new(center, normal, circle.r()))
    }

    /// Fit a circle to a set of 3D points using MAGSAC++ robust consensus estimation in the native
    /// 3D dimension. Each candidate is a full 3D circle and residuals are true point-to-circle
    /// distances, so this is robust to out-of-plane outliers (unlike
    /// [`Circle3::from_consensus_planar`]).
    ///
    /// # Arguments
    ///
    /// * `points`: the 3D points to fit the circle to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `min_r`: an optional minimum radius; candidate circles smaller than this are rejected
    /// * `max_r`: an optional maximum radius; candidate circles larger than this are rejected
    /// * `options`: an optional [`Magsac`] configuration; its `sigma_max` is overridden by the
    ///   `sigma_max` argument
    ///
    /// returns: Result<Circle3, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point3],
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Circle3> {
        let min_r = min_r.unwrap_or(0.0);
        let max_r = max_r.unwrap_or(f64::INFINITY);

        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit =
            magsac.fit_filtered::<3, Circle3, _>(points, |c| c.r() >= min_r && c.r() <= max_r)?;
        Ok(fit.model)
    }
}

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

/// A weighted least-squares problem for refining a 3D circle. The six parameters are an offset
/// applied relative to a fixed base circle: a translation of the center (`[0..3]`), two rotations
/// of the normal about the base circle's local x and y axes (`[3..5]`), and a change in radius
/// (`[5]`). Parameterizing relative to the base avoids any orientation singularity, and the
/// jacobian of the point-to-circle distances is computed by forward finite differences.
struct Circle3Fit<'a> {
    points: &'a [Point3],
    weights: &'a [f64],

    /// Fixed base circle that the parameters are measured relative to.
    base_center: Point3,
    base_normal: UnitVec3,
    base_r: f64,

    /// World-space local x and y axes of the base circle, used as the rotation axes for the normal.
    axis_x: UnitVec3,
    axis_y: UnitVec3,

    params: Vector6<f64>,
    residuals: Residuals,
}

impl<'a> Circle3Fit<'a> {
    fn new(points: &'a [Point3], weights: &'a [f64], base: &Circle3) -> Self {
        let arbitrary = Iso3::from_z_arbitrary_xy(&base.normal, None);
        let axis_x = arbitrary.x();
        let axis_y = arbitrary.y();

        let mut problem = Self {
            points,
            weights,
            base_center: base.center,
            base_normal: base.normal,
            base_r: base.r(),
            axis_x,
            axis_y,
            params: Vector6::zeros(),
            residuals: Residuals::zeros(points.len()),
        };
        problem.set_params(&Vector6::zeros());
        problem
    }

    /// Refine `initial` against the weighted geometric point-to-circle distances of `points` with a
    /// single Levenberg-Marquardt solve, returning the optimized circle or `None` if the solve
    /// fails. This is the shared entry point for both the least-squares and consensus fits.
    fn refine(points: &[Point3], weights: &[f64], initial: &Circle3) -> Option<Circle3> {
        let problem = Circle3Fit::new(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then(|| {
            let (center, normal, radius) = result.circle_params(&result.params);
            Circle3::new(center, normal, radius)
        })
    }

    /// Reconstruct the `(center, normal, radius)` of the circle described by an offset from the base.
    fn circle_params(&self, p: &Vector6<f64>) -> (Point3, UnitVec3, f64) {
        let center = self.base_center + Vector3::new(p[0], p[1], p[2]);
        let rot_x = UnitQuaternion::from_axis_angle(&self.axis_x, p[3]);
        let rot_y = UnitQuaternion::from_axis_angle(&self.axis_y, p[4]);
        let normal = UnitVec3::new_normalize(rot_y * (rot_x * self.base_normal.into_inner()));
        let radius = (self.base_r + p[5]).max(1e-9);
        (center, normal, radius)
    }
}

impl LeastSquaresProblem<f64, Dyn, U6> for Circle3Fit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params = *x;
        let (center, normal, radius) = self.circle_params(&self.params);
        for i in 0..self.points.len() {
            self.residuals[i] =
                self.weights[i] * point_circle_distance(&center, &normal, radius, &self.points[i]);
        }
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        Some(self.residuals.clone())
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        const DELTA: f64 = 1e-7;
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());

        for k in 0..6 {
            let mut p = self.params;
            p[k] += DELTA;
            let (center, normal, radius) = self.circle_params(&p);
            for (i, point) in self.points.iter().enumerate() {
                let d = self.weights[i] * point_circle_distance(&center, &normal, radius, point);
                jac[(i, k)] = (d - self.residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use crate::common::consensus::{ConsensusModel, Magsac};
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom3::Circle3;
    use crate::{Point3, Result, UnitVec3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    fn tilted() -> Circle3 {
        let center = Point3::new(1.0, 2.0, 3.0);
        let normal = UnitVec3::new_normalize(Vector3::new(0.3, -0.5, 1.0));
        Circle3::new(center, normal, 2.5)
    }

    /// Test-only parametrization of a circle's perimeter; `Circle3` has no angle-based API of its
    /// own, so this exists purely to generate sample points for these tests.
    fn sample_circle_point(circle: &Circle3, t: f64) -> Point3 {
        let n = circle.normal.into_inner();
        let reference = if n.z.abs() < 0.9 {
            Vector3::z()
        } else {
            Vector3::x()
        };
        let x_axis = reference.cross(&n).normalize();
        let y_axis = n.cross(&x_axis);
        circle.center + x_axis * (circle.r() * t.cos()) + y_axis * (circle.r() * t.sin())
    }

    fn sample_circle(circle: &Circle3, n: usize) -> Vec<Point3> {
        (0..n)
            .map(|i| sample_circle_point(circle, TAU * i as f64 / n as f64))
            .collect()
    }

    fn assert_matches(result: &Circle3, expected: &Circle3) {
        assert_relative_eq!(result.center, expected.center, epsilon = 5.0e-3);
        assert_relative_eq!(result.r(), expected.r(), epsilon = 5.0e-3);
        // The normal direction of a circle is arbitrary, so compare absolute alignment.
        let dot = result
            .normal
            .into_inner()
            .dot(&expected.normal.into_inner())
            .abs();
        assert_relative_eq!(dot, 1.0, epsilon = 1.0e-3);
    }

    // from_fit tests ─────────────────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_clean_circle() -> Result<()> {
        let expected = tilted();
        let points = sample_circle(&expected, 60);

        let fit = Circle3::from_fit(&points, None)?;
        assert_relative_eq!(fit.center, expected.center, epsilon = 1e-8);
        assert_relative_eq!(fit.r(), expected.r(), epsilon = 1e-8);
        // Every sampled point should lie on the recovered circle.
        for p in &points {
            assert!(fit.residual(p).abs() < 1e-7);
        }
        Ok(())
    }

    #[test]
    fn from_fit_with_noise_is_close() -> Result<()> {
        let mut rg = RandomGeometry3::from_seed(8);
        let expected = tilted();
        let points: Vec<Point3> = sample_circle(&expected, 400)
            .into_iter()
            .map(|p| p + rg.gaussian_vector::<3>(0.01))
            .collect();

        let fit = Circle3::from_fit(&points, None)?;
        assert_matches(&fit, &expected);
        Ok(())
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() -> Result<()> {
        let mut rg = RandomGeometry3::from_seed(21);
        let expected = tilted();
        let points: Vec<Point3> = sample_circle(&expected, 40)
            .into_iter()
            .map(|p| p + rg.gaussian_vector::<3>(0.02))
            .collect();

        let unweighted = Circle3::from_fit(&points, None)?;
        let weights = vec![1.0; points.len()];
        let weighted = Circle3::from_fit(&points, Some(&weights))?;
        assert_relative_eq!(weighted.center, unweighted.center, epsilon = 1e-9);
        assert_relative_eq!(weighted.r(), unweighted.r(), epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn from_fit_too_few_points_is_error() {
        let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        assert!(Circle3::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_collinear_points_is_error() {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        assert!(Circle3::from_fit(&points, None).is_err());
    }

    // from_consensus / construction tests ─────────────────────────────────────

    #[test]
    fn from_3_points_matches_known_circle() -> Result<()> {
        let expected = tilted();
        let p0 = sample_circle_point(&expected, 0.0);
        let p1 = sample_circle_point(&expected, 2.0);
        let p2 = sample_circle_point(&expected, -1.5);
        let circle = Circle3::from_3_points(&p0, &p1, &p2)?;
        assert_matches(&circle, &expected);
        Ok(())
    }

    #[test]
    fn planar_recovers_circle() -> Result<()> {
        let expected = tilted();
        let points = sample_circle(&expected, 100);
        let result = Circle3::from_consensus_planar(&points, 0.01, None, None, None)?;
        assert_matches(&result, &expected);
        Ok(())
    }

    #[test]
    fn native_recovers_circle_with_outliers() -> Result<()> {
        let expected = tilted();
        let mut points = sample_circle(&expected, 120);
        let inlier_count = points.len();

        // A scattered cloud of gross outliers, well off the circle and its plane.
        for i in 0..40 {
            let f = i as f64;
            points.push(
                expected.center
                    + Vector3::new(
                        5.0 + 2.0 * (1.3 * f).sin(),
                        -4.0 + 2.0 * (0.7 * f).cos(),
                        3.0 + (2.1 * f).sin(),
                    ),
            );
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(500),
            refinement_steps: 5,
            confidence: 0.99,
            seed: Some(17),
        };
        let fit = magsac.fit_filtered::<3, Circle3, _>(&points, |_| true)?;

        assert_matches(&fit.model, &expected);
        // No outlier should be classified as an inlier.
        assert!(fit.inliers.iter().all(|&i| i < inlier_count));
        assert!(fit.inliers.len() > inlier_count * 9 / 10);
        Ok(())
    }
}
