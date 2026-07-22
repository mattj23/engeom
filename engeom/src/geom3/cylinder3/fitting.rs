//! Cylinder fitting for [`Cylinder3`] by ordinary least squares against geometric radial residuals.
//!
//! There are two constructors, both refining through the same shared [`CylinderFit`]
//! Levenberg-Marquardt engine (which works over the infinite cylinder's five degrees of freedom:
//! the axis line and the radius):
//!   - [`Cylinder3::from_fit`] takes bare positions plus a caller-supplied initial guess. Positions
//!     alone don't determine a cheap closed-form axis, so an initial `Cylinder3` is required to seed
//!     the refinement.
//!   - [`Cylinder3::from_fit_oriented`] takes oriented [`SurfacePoint3`] samples and self-bootstraps
//!     the axis from the normals (every surface normal is perpendicular to a cylinder's axis, so the
//!     axis is the eigenvector of the normal scatter matrix `Σ wᵢ nᵢ nᵢᵀ` for its smallest
//!     eigenvalue), then estimates a point on the axis and the radius with an in-plane algebraic
//!     circle fit before refining.
//!
//! The normals are used *only* to bootstrap the axis in [`Cylinder3::from_fit_oriented`]; the shared
//! engine and the residual are position-only. The radial residual `dist(point, axis) − radius` does
//! not depend on which side the surface faces, and the scatter matrix `nᵢ nᵢᵀ` is sign-invariant, so
//! inner- and outer-facing cylinders (inward vs. outward normals) fit identically.
//!
//! The least-squares refinement estimates the *infinite* cylinder (axis line + radius); the returned
//! [`Cylinder3`]'s `center` and `length` are then set to span the axial extent of the input points,
//! exactly as [`Segment`](crate::common::Segment) bounds a fitted line.

use super::Cylinder3;
use crate::common::PCoords;
use crate::geom3::IsoExtensions3;
use crate::{Iso3, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Matrix3, Owned, U1, U5, UnitQuaternion, Vector, Vector5};

/// The signed radial distance from a point to an infinite cylinder's surface: the perpendicular
/// distance from the point to the axis line, minus the radius. It is positive outside the cylinder
/// and negative inside, and its sign does not depend on which way the surface faces, so inner- and
/// outer-facing cylinders produce identical residuals.
fn point_cylinder_distance(
    center: &Point3,
    direction: &UnitVec3,
    radius: f64,
    point: &Point3,
) -> f64 {
    let v = point - center;
    let axial = v.dot(&direction.into_inner());
    let radial = v - direction.into_inner() * axial;
    radial.norm() - radius
}

/// Estimate the cylinder axis direction from oriented surface points. On a cylinder every surface
/// normal is perpendicular to the axis, so the axis is the direction least aligned with the normals:
/// the eigenvector of the weighted normal scatter matrix `M = Σ wᵢ nᵢ nᵢᵀ` for its smallest
/// eigenvalue. `nᵢ nᵢᵀ` is sign-invariant, so inconsistently oriented (inward vs. outward) normals
/// give the same result.
///
/// Returns `None` if the normals do not constrain an axis: for a near-planar patch the normals are
/// nearly parallel, leaving the axis direction ambiguous. This is detected when the second-smallest
/// eigenvalue is not clearly separated from zero relative to the largest.
fn axis_from_normals(points: &[SurfacePoint3], weights: &[f64]) -> Option<UnitVec3> {
    let mut m = Matrix3::<f64>::zeros();
    for (sp, &w) in points.iter().zip(weights) {
        let n = sp.normal.into_inner();
        m += w * n * n.transpose();
    }

    let eigen = m.symmetric_eigen();
    let ev = eigen.eigenvalues;
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| ev[a].partial_cmp(&ev[b]).unwrap());
    let (smallest, second, largest) = (order[0], order[1], order[2]);

    // The axis is well-defined only when the normals span a genuine band around it, which requires
    // the second-smallest eigenvalue to be a meaningful fraction of the largest. When the sample is
    // nearly planar the two smallest eigenvalues are both ~0 and the axis is ambiguous.
    if ev[largest] <= 0.0 || ev[second] / ev[largest] < 1e-6 {
        return None;
    }

    let axis = eigen.eigenvectors.column(smallest).into_owned();
    UnitVec3::try_new(axis, 1e-9)
}

/// Given an axis direction, estimate a point on the axis and the radius by projecting the points into
/// the plane perpendicular to the axis and running a weighted algebraic (Kåsa-style) circle fit.
/// Returns `None` if the projected points are collinear or otherwise fail to determine a circle.
fn axis_point_and_radius(
    direction: &UnitVec3,
    points: &[Point3],
    weights: &[f64],
) -> Option<(Point3, f64)> {
    let arbitrary = Iso3::from_z_arbitrary_xy(direction, None);
    let u = arbitrary.x().into_inner();
    let v = arbitrary.y().into_inner();

    // Use the weighted centroid as the in-plane origin to keep the algebraic system well-scaled.
    let mut wsum = 0.0;
    let mut acc = Vector3::zeros();
    for (p, &w) in points.iter().zip(weights) {
        acc += w * p.coords;
        wsum += w;
    }
    if wsum <= 0.0 {
        return None;
    }
    let origin = Point3::from(acc / wsum);

    // Weighted Kåsa fit in the `(u, v)` cross-section: for `x² + y² = a·x + b·y + c`, each projected
    // point gives a linear equation in `s = [a, b, c]`, with `a = 2cx`, `b = 2cy`, `c = r² − cx² −
    // cy²`.
    let mut m = Matrix3::<f64>::zeros();
    let mut b = Vector3::zeros();
    for (p, &w) in points.iter().zip(weights) {
        let d = p - origin;
        let a = d.dot(&u);
        let bb = d.dot(&v);
        let row = Vector3::new(a, bb, 1.0);
        let target = a * a + bb * bb;
        m += w * row * row.transpose();
        b += w * target * row;
    }

    let s = m.lu().solve(&b)?;
    let cu = 0.5 * s[0];
    let cv = 0.5 * s[1];
    let r_squared = s[2] + cu * cu + cv * cv;
    if r_squared <= 0.0 {
        return None;
    }

    let axis_point = origin + u * cu + v * cv;
    Some((axis_point, r_squared.sqrt()))
}

/// Convert the infinite cylinder produced by the least-squares refinement into a finite one by
/// setting `center` and `length` to span the axial extent of `points` along the refined axis. The
/// axis line and radius are left unchanged.
fn bound_axially(cyl: &Cylinder3, points: &[Point3]) -> Cylinder3 {
    let axis = cyl.axis();
    let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in points {
        let t = axis.scalar_project(p);
        t_min = t_min.min(t);
        t_max = t_max.max(t);
    }
    let center = axis.at(t_min);
    Cylinder3::new(center, cyl.direction, cyl.radius, (t_max - t_min).max(0.0))
}

impl Cylinder3 {
    /// Fit a cylinder to a set of points by ordinary least squares, starting from a caller-supplied
    /// initial guess. The infinite cylinder (axis line and radius) is refined against the true
    /// geometric radial residuals with the shared weighted [`CylinderFit`] Levenberg-Marquardt
    /// engine, then bounded to the axial extent of the points.
    ///
    /// Unlike most primitives, this requires an initial `guess`: bare positions do not admit a cheap
    /// closed-form axis estimate. If you have surface normals, use [`Cylinder3::from_fit_oriented`],
    /// which bootstraps the axis for you.
    ///
    /// The normals are irrelevant here; the residual is `dist(point, axis) − radius`, so this fits
    /// inner- and outer-facing cylinders identically. It is not robust to gross outliers.
    ///
    /// # Arguments
    ///
    /// * `points`: the coordinates to fit the cylinder to (at least two, though a well-posed fit
    ///   needs enough points to constrain all five degrees of freedom)
    /// * `guess`: an initial cylinder whose axis and radius seed the refinement
    /// * `weights`: if `Some`, a slice the same length as `points` giving the weight to multiply each
    ///   point residual by; if `None`, all points are weighted equally
    ///
    /// returns: Result<Cylinder3, Box<dyn Error, Global>>
    pub fn from_fit(
        points: &[impl PCoords<3>],
        guess: &Cylinder3,
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        let pts: Vec<Point3> = points.iter().map(|p| Point3::from(p.coords())).collect();
        if pts.len() < 2 {
            return Err("At least two points are required to fit a cylinder".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; pts.len()];
                &ones
            }
        };

        let refined =
            CylinderFit::refine(&pts, weights, guess).ok_or("Failed to refine cylinder fit")?;
        Ok(bound_axially(&refined, &pts))
    }

    /// Fit a cylinder to a set of oriented surface points by ordinary least squares, bootstrapping
    /// the initial axis from the normals. Every surface normal on a cylinder is perpendicular to the
    /// axis, so the axis direction is estimated from the normal scatter matrix (see the module
    /// documentation); a point on the axis and the radius then come from an in-plane algebraic circle
    /// fit. The resulting guess is refined against the true geometric radial residuals with the same
    /// weighted [`CylinderFit`] engine used by [`Cylinder3::from_fit`], then bounded to the axial
    /// extent of the points.
    ///
    /// The normals only seed the axis estimate; the residual is the position-only
    /// `dist(point, axis) − radius`, and the axis bootstrap is sign-invariant, so inner- and
    /// outer-facing cylinders fit identically. It is not robust to gross outliers.
    ///
    /// # Arguments
    ///
    /// * `points`: the oriented surface points to fit the cylinder to (at least three)
    /// * `weights`: if `Some`, a slice the same length as `points` giving the weight to multiply each
    ///   point residual by; if `None`, all points are weighted equally
    ///
    /// returns: Result<Cylinder3, Box<dyn Error, Global>>
    pub fn from_fit_oriented(points: &[SurfacePoint3], weights: Option<&[f64]>) -> Result<Self> {
        if points.len() < 3 {
            return Err("At least three oriented points are required to fit a cylinder".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; points.len()];
                &ones
            }
        };

        let direction = axis_from_normals(points, weights)
            .ok_or("Normals do not constrain a cylinder axis (near-planar sample)")?;

        let pts: Vec<Point3> = points.iter().map(|sp| sp.point).collect();
        let (axis_point, radius) = axis_point_and_radius(&direction, &pts, weights)
            .ok_or("Failed to estimate cylinder radius from points")?;
        let guess = Cylinder3::new(axis_point, direction, radius, 0.0);

        let refined =
            CylinderFit::refine(&pts, weights, &guess).ok_or("Failed to refine cylinder fit")?;
        Ok(bound_axially(&refined, &pts))
    }
}

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

/// A weighted least-squares problem for refining an infinite cylinder (axis line and radius). The
/// five parameters are an offset applied relative to a fixed base cylinder: a translation of the
/// axis point within the plane perpendicular to the axis (`[0..2]`, along the base's perpendicular
/// axes — translation *along* the axis is a gauge freedom and is deliberately excluded), two
/// rotations of the axis direction about those perpendicular axes (`[2..4]`), and a change in radius
/// (`[4]`). Parameterizing relative to the base avoids any orientation singularity, and the jacobian
/// of the point-to-surface distances is computed by forward finite differences.
struct CylinderFit<'a> {
    points: &'a [Point3],
    weights: &'a [f64],

    /// Fixed base point on the axis that the parameters are measured relative to.
    base_center: Point3,
    base_direction: UnitVec3,
    base_r: f64,

    /// World-space axes spanning the plane perpendicular to the base axis, used both as the
    /// translation directions and as the rotation axes for the direction.
    axis_x: UnitVec3,
    axis_y: UnitVec3,

    params: Vector5<f64>,
    residuals: Residuals,
}

impl<'a> CylinderFit<'a> {
    fn new(points: &'a [Point3], weights: &'a [f64], base: &Cylinder3) -> Self {
        let arbitrary = Iso3::from_z_arbitrary_xy(&base.direction, None);
        let axis_x = arbitrary.x();
        let axis_y = arbitrary.y();

        let mut problem = Self {
            points,
            weights,
            base_center: base.center,
            base_direction: base.direction,
            base_r: base.r(),
            axis_x,
            axis_y,
            params: Vector5::zeros(),
            residuals: Residuals::zeros(points.len()),
        };
        problem.set_params(&Vector5::zeros());
        problem
    }

    /// Refine `initial` against the weighted geometric radial residuals of `points` with a single
    /// Levenberg-Marquardt solve, returning the optimized infinite cylinder (with an unbounded axis;
    /// its `center` lies on the axis and its `length` is unset) or `None` if the solve fails. This is
    /// the shared entry point for both fitting constructors and, later, the consensus refinement.
    fn refine(points: &[Point3], weights: &[f64], initial: &Cylinder3) -> Option<Cylinder3> {
        let problem = CylinderFit::new(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then(|| {
            let (center, direction, radius) = result.cylinder_params(&result.params);
            Cylinder3::new(center, direction, radius, 0.0)
        })
    }

    /// Reconstruct the `(axis point, direction, radius)` of the cylinder described by an offset from
    /// the base.
    fn cylinder_params(&self, p: &Vector5<f64>) -> (Point3, UnitVec3, f64) {
        let center =
            self.base_center + self.axis_x.into_inner() * p[0] + self.axis_y.into_inner() * p[1];
        let rot_x = UnitQuaternion::from_axis_angle(&self.axis_x, p[2]);
        let rot_y = UnitQuaternion::from_axis_angle(&self.axis_y, p[3]);
        let direction = UnitVec3::new_normalize(rot_y * (rot_x * self.base_direction.into_inner()));
        let radius = (self.base_r + p[4]).max(1e-9);
        (center, direction, radius)
    }
}

impl LeastSquaresProblem<f64, Dyn, U5> for CylinderFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U5>;
    type ParameterStorage = Owned<f64, U5>;

    fn set_params(&mut self, x: &Vector<f64, U5, Self::ParameterStorage>) {
        self.params = *x;
        let (center, direction, radius) = self.cylinder_params(&self.params);
        for i in 0..self.points.len() {
            self.residuals[i] = self.weights[i]
                * point_cylinder_distance(&center, &direction, radius, &self.points[i]);
        }
    }

    fn params(&self) -> Vector<f64, U5, Self::ParameterStorage> {
        self.params
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        Some(self.residuals.clone())
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U5, Self::JacobianStorage>> {
        const DELTA: f64 = 1e-7;
        let mut jac = Matrix::<f64, Dyn, U5, Self::JacobianStorage>::zeros(self.points.len());

        for k in 0..5 {
            let mut p = self.params;
            p[k] += DELTA;
            let (center, direction, radius) = self.cylinder_params(&p);
            for (i, point) in self.points.iter().enumerate() {
                let d =
                    self.weights[i] * point_cylinder_distance(&center, &direction, radius, point);
                jac[(i, k)] = (d - self.residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom3::Cylinder3;
    use crate::{Point3, UnitVec3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    fn tilted() -> Cylinder3 {
        let center = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::new(0.3, -0.5, 1.0));
        Cylinder3::new(center, direction, 1.5, 6.0)
    }

    /// Deterministically sample `n` oriented surface points spanning the full length of `cyl`. The
    /// axial positions include both endpoints (so the recovered length is exact for clean data) and
    /// the angles are spread by the golden angle so the normals cover a full band around the axis.
    /// When `outward` is false the normals point toward the axis (inner-facing surface).
    fn sample_surface(cyl: &Cylinder3, n: usize, outward: bool) -> Vec<SurfacePoint3> {
        let arbitrary = Iso3::from_z_arbitrary_xy(&cyl.direction, None);
        let u = arbitrary.x().into_inner();
        let v = arbitrary.y().into_inner();
        let axis = cyl.axis();
        (0..n)
            .map(|i| {
                let f = i as f64 / (n - 1) as f64;
                let t = cyl.length * f;
                let ang = TAU * i as f64 * 0.6180339887;
                let radial = u * ang.cos() + v * ang.sin();
                let point = axis.at(t) + radial * cyl.r();
                let normal = if outward { radial } else { -radial };
                SurfacePoint3::new_normalize(point, normal)
            })
            .collect()
    }

    /// Add isotropic Gaussian noise to the positions of a set of surface points, leaving the normals
    /// unchanged.
    fn with_noise(
        rg: &mut RandomGeometry3,
        points: &[SurfacePoint3],
        sigma: f64,
    ) -> Vec<SurfacePoint3> {
        points
            .iter()
            .map(|sp| SurfacePoint3::new(sp.point + rg.gaussian_vector::<3>(sigma), sp.normal))
            .collect()
    }

    /// Assert two cylinders describe the same infinite geometry (coincident axis line, parallel
    /// direction, equal radius) and axial extent, tolerant of axis sign and of which end is `a()`.
    fn assert_cyl_matches(result: &Cylinder3, expected: &Cylinder3, eps: f64) {
        assert_relative_eq!(result.r(), expected.r(), epsilon = eps);
        assert_relative_eq!(result.length, expected.length, epsilon = eps);
        // Direction parallel or anti-parallel.
        let dot = result
            .direction
            .into_inner()
            .dot(&expected.direction.into_inner())
            .abs();
        assert_relative_eq!(dot, 1.0, epsilon = eps);
        // The two axis lines must coincide: each cylinder's base point lies on the other's axis.
        assert!(expected.axis().distance_to(&result.center) < eps);
        assert!(result.axis().distance_to(&expected.center) < eps);
    }

    // from_fit_oriented tests ─────────────────────────────────────────────────

    #[test]
    fn from_fit_oriented_recovers_clean_cylinder() {
        let expected = tilted();
        let points = sample_surface(&expected, 60, true);
        let fit = Cylinder3::from_fit_oriented(&points, None).unwrap();
        assert_cyl_matches(&fit, &expected, 1e-6);
    }

    #[test]
    fn from_fit_oriented_inner_facing_matches_outer_facing() {
        // Inner- and outer-facing surfaces differ only in normal sign; the fit must be identical.
        let expected = tilted();
        let outer =
            Cylinder3::from_fit_oriented(&sample_surface(&expected, 60, true), None).unwrap();
        let inner =
            Cylinder3::from_fit_oriented(&sample_surface(&expected, 60, false), None).unwrap();

        assert_cyl_matches(&inner, &expected, 1e-6);
        assert_relative_eq!(inner.center, outer.center, epsilon = 1e-9);
        assert_relative_eq!(
            inner.direction.into_inner(),
            outer.direction.into_inner(),
            epsilon = 1e-9
        );
        assert_relative_eq!(inner.r(), outer.r(), epsilon = 1e-9);
        assert_relative_eq!(inner.length, outer.length, epsilon = 1e-9);
    }

    #[test]
    fn from_fit_oriented_with_noise_is_close() {
        let mut rg = RandomGeometry3::from_seed(8);
        let expected = tilted();
        let clean = sample_surface(&expected, 400, true);
        let points = with_noise(&mut rg, &clean, 0.01);

        let fit = Cylinder3::from_fit_oriented(&points, None).unwrap();
        assert_cyl_matches(&fit, &expected, 0.05);
    }

    #[test]
    fn from_fit_oriented_uniform_weights_match_unweighted() {
        let mut rg = RandomGeometry3::from_seed(21);
        let expected = tilted();
        let points = with_noise(&mut rg, &sample_surface(&expected, 40, true), 0.02);

        let unweighted = Cylinder3::from_fit_oriented(&points, None).unwrap();
        let weights = vec![1.0; points.len()];
        let weighted = Cylinder3::from_fit_oriented(&points, Some(&weights)).unwrap();
        assert_relative_eq!(weighted.center, unweighted.center, epsilon = 1e-9);
        assert_relative_eq!(
            weighted.direction.into_inner(),
            unweighted.direction.into_inner(),
            epsilon = 1e-9
        );
        assert_relative_eq!(weighted.r(), unweighted.r(), epsilon = 1e-9);
    }

    #[test]
    fn from_fit_oriented_too_few_points_is_error() {
        let expected = tilted();
        let points = sample_surface(&expected, 2, true);
        assert!(Cylinder3::from_fit_oriented(&points, None).is_err());
    }

    #[test]
    fn from_fit_oriented_parallel_normals_is_error() {
        // A flat patch: every point shares the same normal (all on one axis-parallel line of the
        // surface). Parallel normals leave the axis direction ambiguous, so the bootstrap must
        // reject them rather than return a garbage cylinder.
        let expected = tilted();
        let arbitrary = Iso3::from_z_arbitrary_xy(&expected.direction, None);
        let u = arbitrary.x().into_inner();
        let axis = expected.axis();
        let points: Vec<SurfacePoint3> = (0..30)
            .map(|i| {
                let t = expected.length * i as f64 / 29.0;
                // Constant radial direction `u` -> identical normals for every sample.
                SurfacePoint3::new_normalize(axis.at(t) + u * expected.r(), u)
            })
            .collect();
        assert!(Cylinder3::from_fit_oriented(&points, None).is_err());
    }

    // from_fit (guess-based) tests ─────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_clean_cylinder_from_guess() {
        let expected = tilted();
        let surface = sample_surface(&expected, 60, true);
        let points: Vec<Point3> = surface.iter().map(|sp| sp.point).collect();

        // A perturbed but nearby initial guess.
        let guess = Cylinder3::new(
            expected.center + Vector3::new(0.2, -0.15, 0.1),
            UnitVec3::new_normalize(
                expected.direction.into_inner() + Vector3::new(0.05, 0.05, 0.0),
            ),
            expected.r() + 0.3,
            expected.length,
        );

        let fit = Cylinder3::from_fit(&points, &guess, None).unwrap();
        assert_cyl_matches(&fit, &expected, 1e-6);
    }

    #[test]
    fn from_fit_with_noise_is_close() {
        let mut rg = RandomGeometry3::from_seed(13);
        let expected = tilted();
        let surface = with_noise(&mut rg, &sample_surface(&expected, 400, true), 0.01);
        let points: Vec<Point3> = surface.iter().map(|sp| sp.point).collect();

        let guess = Cylinder3::new(
            expected.center + Vector3::new(0.1, 0.1, -0.1),
            expected.direction,
            expected.r() + 0.2,
            expected.length,
        );
        let fit = Cylinder3::from_fit(&points, &guess, None).unwrap();
        assert_cyl_matches(&fit, &expected, 0.05);
    }

    #[test]
    fn from_fit_too_few_points_is_error() {
        let guess = tilted();
        let points = [Point3::new(0.0, 0.0, 0.0)];
        assert!(Cylinder3::from_fit(&points, &guess, None).is_err());
    }
}
