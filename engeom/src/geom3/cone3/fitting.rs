//! Cone fitting for [`Cone3`] by ordinary least squares and MAGSAC++ robust consensus.
//!
//! There are three constructors, all refining through the same shared [`ConeFit`]
//! Levenberg-Marquardt engine (which works over the infinite cone's six degrees of freedom: the apex
//! position, the axis direction, and the half-angle):
//!   - [`Cone3::from_fit`] takes bare positions plus a caller-supplied initial guess. Positions
//!     alone don't determine a cheap closed-form cone, so an initial `Cone3` is required to seed the
//!     refinement.
//!   - [`Cone3::from_fit_oriented`] takes oriented [`SurfacePoint3`] samples and self-bootstraps the
//!     whole cone from the point/normal pairs. Every tangent plane of a cone passes through the apex,
//!     so the apex is the least-squares solution of `nᵢ·a = nᵢ·pᵢ`; the unit directions from the
//!     apex to the points, `uᵢ = normalize(pᵢ − a)`, then lie on a plane whose normal is the axis and
//!     whose offset is `cos α`, giving the axis direction and half-angle in closed form.
//!   - [`Cone3::from_consensus`] takes oriented [`SurfacePoint3`] samples and runs the MAGSAC++
//!     robust estimator, rejecting gross outliers. Each minimal sample is four oriented points (the
//!     same minimum [`Cone3::from_fit_oriented`] requires), reusing the same oriented bootstrap.
//!     Scoring and refinement use only the point positions, so the normals matter only for building
//!     candidates.
//!
//! The normals are used *only* to bootstrap the apex in [`Cone3::from_fit_oriented`] and
//! [`Cone3::from_consensus`]; the shared engine and the residual are position-only. The residual
//! (see [`point_cone_distance`]) treats the "infinite" cone as a single nappe: it opens without
//! bound in the positive axis direction from the apex, but does not extend backward through it, so a
//! point behind the apex is measured against the apex itself rather than the mirror-image nappe that
//! a naive doubly-infinite cone would have. It does not depend on which side the surface faces, and
//! the apex system `nᵢ·a = nᵢ·pᵢ` is invariant under flipping any normal, so inner- and outer-facing
//! cones (inward vs. outward normals, even mixed) fit identically.
//!
//! The least-squares refinement estimates the *infinite* cone (apex + axis + half-angle); the
//! returned [`Cone3`]'s `height` and base `radius` are then set to span the axial extent of the input
//! points from the apex, exactly as [`Cylinder3`](crate::geom3::Cylinder3) bounds a fitted axis.

use super::Cone3;
use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::geom3::{IsoExtensions3, SvdBasis3};
use crate::{Iso3, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Matrix3, Owned, U1, U6, UnitQuaternion, Vector, Vector6};
use std::f64::consts::FRAC_PI_2;

/// Small tolerance for degeneracy checks and non-zero direction extraction.
const EPSILON: f64 = 1e-10;

/// The signed distance from a point to a single-nappe infinite cone's lateral surface: the cone
/// opens without bound in the positive `direction` from `apex`, but does not extend backward through
/// it. With `q` the signed axial coordinate of the point from the apex and `ρ` its radial distance
/// from the axis, the slant surface at half-angle `α` has perpendicular distance `ρ·cos α − q·sin α`
/// (positive outside the cone, negative inside) *when the point's perpendicular foot on the slant
/// line falls at or past the apex* (`t = q·cos α + ρ·sin α ≥ 0`). Otherwise the closest feature is
/// the apex itself, and the distance is the plain Euclidean distance to it, `|point − apex|`. Without
/// this apex clamp, a point well behind the apex could register a small (or even negative) distance
/// against the mirror-image nappe that a naive doubly-infinite cone would have; the clamp keeps the
/// residual growing without bound behind the apex instead. The result is exact near the surface -
/// where cone fits operate - and its sign does not depend on which way the surface faces, so inner-
/// and outer-facing cones produce identical residuals.
fn point_cone_distance(
    apex: &Point3,
    direction: &UnitVec3,
    half_angle: f64,
    point: &Point3,
) -> f64 {
    let v = point - apex;
    let q = v.dot(&direction.into_inner());
    let rho = (v - direction.into_inner() * q).norm();

    let (s, c) = half_angle.sin_cos();
    let t = q * c + rho * s;
    if t < 0.0 { v.norm() } else { rho * c - q * s }
}

/// Estimate the cone apex from oriented surface points. Every tangent plane of a cone passes through
/// the apex, so for each surface point `n·(a − p) = 0`, i.e. `n·a = n·p`; the apex is the
/// least-squares solution of that linear system, `a = (Σ wᵢ nᵢnᵢᵀ)⁻¹ (Σ wᵢ (nᵢ·pᵢ) nᵢ)`. Flipping
/// any normal negates both sides of its equation and leaves the system unchanged, so this is
/// invariant to the surface facing (including mixed inward/outward normals).
///
/// Returns `None` if the normals do not span 3D (the `Σ nᵢnᵢᵀ` matrix is rank-deficient), which
/// happens for a near-cylindrical sample where the normals are all perpendicular to a common axis
/// and the apex recedes to infinity.
fn apex_from_normals(points: &[SurfacePoint3], weights: &[f64]) -> Option<Point3> {
    let mut m = Matrix3::<f64>::zeros();
    let mut b = Vector3::zeros();
    for (sp, &w) in points.iter().zip(weights) {
        let n = sp.normal.into_inner();
        m += w * n * n.transpose();
        b += w * n.dot(&sp.point.coords) * n;
    }

    // Rank guard: the apex is determined only if the normals span 3D. For a (near-)cylindrical
    // sample the smallest eigenvalue of the scatter matrix collapses toward zero.
    let ev = m.symmetric_eigen().eigenvalues;
    let max = ev.max();
    if max <= 0.0 || ev.min() / max < 1e-6 {
        return None;
    }

    let a = m.lu().solve(&b)?;
    Some(Point3::from(a))
}

/// Given the apex, estimate the axis direction and half-angle. The unit directions from the apex to
/// the points, `uᵢ = normalize(pᵢ − a)`, all satisfy `uᵢ·d = cos α`, so they lie on a plane with
/// normal `d` (the axis) at offset `cos α`. Fitting that plane with an [`SvdBasis3`] recovers the
/// axis robustly even from a partial azimuthal arc, and the half-angle comes from the centroid's
/// projection onto the axis (`cos α = mean(uᵢ)·d`, exact regardless of azimuthal coverage). The axis
/// is oriented to point from the apex toward the points. Uses only positions, so it is facing-
/// invariant.
///
/// Returns `None` if the directions are degenerate or the recovered half-angle is not a proper cone
/// angle (too close to `0`, a cylinder, or `π/2`, a flat disc).
fn axis_and_half_angle(
    apex: &Point3,
    points: &[Point3],
    weights: &[f64],
) -> Option<(UnitVec3, f64)> {
    let mut dirs = Vec::with_capacity(points.len());
    let mut dir_weights = Vec::with_capacity(points.len());
    for (p, &w) in points.iter().zip(weights) {
        if let Some(u) = UnitVec3::try_new(p - apex, EPSILON) {
            dirs.push(Point3::from(u.into_inner()));
            dir_weights.push(w);
        }
    }
    if dirs.len() < 3 {
        return None;
    }

    let basis = SvdBasis3::from_points(&dirs, Some(&dir_weights))?;
    let mut axis = basis.smallest().into_inner();

    // `mean(uᵢ)·d = cos α`; orient the axis to point from the apex toward the points.
    let mut cos_alpha = basis.center.coords.dot(&axis);
    if cos_alpha < 0.0 {
        axis = -axis;
        cos_alpha = -cos_alpha;
    }

    let half_angle = cos_alpha.clamp(-1.0, 1.0).acos();
    if !(EPSILON..=FRAC_PI_2 - EPSILON).contains(&half_angle) {
        return None;
    }

    Some((UnitVec3::new_normalize(axis), half_angle))
}

/// Build an unbounded initial cone from oriented points: the apex from the tangent-plane system,
/// then the axis direction and half-angle from the plane the apex-relative unit directions lie on
/// (see [`apex_from_normals`] and [`axis_and_half_angle`]). This is the shared bootstrap used by both
/// [`Cone3::from_fit_oriented`] and the consensus minimal sample. The half-angle is encoded in the
/// returned cone's `height`/`radius` as `(cos α, sin α)`; it is meant to be refined and then bounded.
/// Returns `None` if the sample does not constrain a cone (near-cylindrical normals, degenerate
/// directions, or an improper half-angle).
fn bootstrap_cone(points: &[SurfacePoint3], weights: &[f64]) -> Option<Cone3> {
    let apex = apex_from_normals(points, weights)?;
    let pts: Vec<Point3> = points.iter().map(|sp| sp.point).collect();
    let (direction, half_angle) = axis_and_half_angle(&apex, &pts, weights)?;
    Some(Cone3::new(
        apex,
        direction,
        half_angle.cos(),
        half_angle.sin(),
    ))
}

/// Convert the infinite cone produced by the least-squares refinement into a finite one by setting
/// `height` and base `radius` to span the axial extent of `points` from the apex along the refined
/// axis. The apex, axis direction, and half-angle are left unchanged.
fn bound_axially(cone: &Cone3, points: &[Point3]) -> Cone3 {
    let axis = cone.axis();
    let mut q_max = 0.0_f64;
    for p in points {
        q_max = q_max.max(axis.scalar_project(p));
    }
    let height = q_max.max(0.0);
    let radius = height * cone.half_angle().tan();
    Cone3::new(cone.tip, cone.direction, height, radius)
}

impl Cone3 {
    /// Fit a cone to a set of points by ordinary least squares, starting from a caller-supplied
    /// initial guess. The infinite cone (apex, axis, and half-angle) is refined against the true
    /// geometric surface residuals with the shared weighted [`ConeFit`] Levenberg-Marquardt engine,
    /// then bounded to the axial extent of the points.
    ///
    /// Unlike most primitives, this requires an initial `guess`: bare positions do not admit a cheap
    /// closed-form cone estimate. If you have surface normals, use [`Cone3::from_fit_oriented`],
    /// which bootstraps the cone for you.
    ///
    /// The normals are irrelevant here; the residual is the orthogonal distance to the lateral
    /// surface, so this fits inner- and outer-facing cones identically. It is not robust to gross
    /// outliers.
    ///
    /// # Arguments
    ///
    /// * `points`: the coordinates to fit the cone to (at least two, though a well-posed fit needs
    ///   enough points to constrain all six degrees of freedom)
    /// * `guess`: an initial cone whose apex, axis, and half-angle seed the refinement
    /// * `weights`: if `Some`, a slice the same length as `points` giving the weight to multiply each
    ///   point residual by; if `None`, all points are weighted equally
    ///
    /// returns: Result<Cone3, Box<dyn Error, Global>>
    pub fn from_fit(
        points: &[impl PCoords<3>],
        guess: &Cone3,
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        let pts: Vec<Point3> = points.iter().map(|p| Point3::from(p.coords())).collect();
        if pts.len() < 2 {
            return Err("At least two points are required to fit a cone".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; pts.len()];
                &ones
            }
        };

        let refined = ConeFit::refine(&pts, weights, guess).ok_or("Failed to refine cone fit")?;
        Ok(bound_axially(&refined, &pts))
    }

    /// Fit a cone to a set of oriented surface points by ordinary least squares, bootstrapping the
    /// entire initial cone from the point/normal pairs. The apex is found from the tangent-plane
    /// system `nᵢ·a = nᵢ·pᵢ`, and the axis direction and half-angle from the plane the apex-relative
    /// unit directions lie on (see the module documentation). The resulting guess is refined against
    /// the true geometric surface residuals with the same weighted [`ConeFit`] engine used by
    /// [`Cone3::from_fit`], then bounded to the axial extent of the points.
    ///
    /// The normals only seed the apex estimate; the residual is position-only and every bootstrap
    /// step is facing-invariant, so inner- and outer-facing cones (inward vs. outward normals, even
    /// mixed) fit identically. It is not robust to gross outliers.
    ///
    /// # Arguments
    ///
    /// * `points`: the oriented surface points to fit the cone to (at least four, spread around the
    ///   cone so the normals span 3D)
    /// * `weights`: if `Some`, a slice the same length as `points` giving the weight to multiply each
    ///   point residual by; if `None`, all points are weighted equally
    ///
    /// returns: Result<Cone3, Box<dyn Error, Global>>
    pub fn from_fit_oriented(points: &[SurfacePoint3], weights: Option<&[f64]>) -> Result<Self> {
        if points.len() < 4 {
            return Err("At least four oriented points are required to fit a cone".into());
        }

        let ones;
        let weights = match weights {
            Some(w) => w,
            None => {
                ones = vec![1.0; points.len()];
                &ones
            }
        };

        let guess = bootstrap_cone(points, weights).ok_or(
            "Failed to estimate an initial cone from the oriented points \
             (near-cylindrical normals or degenerate sample)",
        )?;

        let pts: Vec<Point3> = points.iter().map(|sp| sp.point).collect();
        let refined = ConeFit::refine(&pts, weights, &guess).ok_or("Failed to refine cone fit")?;
        Ok(bound_axially(&refined, &pts))
    }

    /// Fit a cone to a set of oriented surface points using MAGSAC++ robust consensus estimation.
    /// Unlike [`Cone3::from_fit_oriented`], this rejects gross outliers by taking an upper bound on
    /// the inlier noise (`sigma_max`) rather than a hard inlier/outlier threshold.
    ///
    /// Each minimal sample is four oriented points, the same minimum [`Cone3::from_fit_oriented`]
    /// requires: the apex comes from the tangent-plane system and the axis/half-angle from the plane
    /// the apex-relative unit directions lie on. Candidate scoring and refinement use only the point
    /// positions (the orthogonal distance to the lateral surface), so the normals are used solely to
    /// build candidates; consequently inner- and outer-facing surfaces (and mixed normals) are
    /// handled identically. The returned cone is bounded to the axial extent of the inlier points.
    ///
    /// Because the infinite cone's `radius` field merely encodes the half-angle until it is bounded,
    /// candidates are filtered by half-angle rather than by radius.
    ///
    /// # Arguments
    ///
    /// * `points`: the oriented surface points to fit the cone to (at least four)
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `min_half_angle`: an optional minimum half-angle, in radians; candidate cones narrower than
    ///   this are rejected
    /// * `max_half_angle`: an optional maximum half-angle, in radians; candidate cones wider than
    ///   this are rejected
    /// * `options`: an optional [`Magsac`] configuration; its `sigma_max` is overridden by the
    ///   `sigma_max` argument
    ///
    /// returns: Result<Cone3, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[SurfacePoint3],
        sigma_max: f64,
        min_half_angle: Option<f64>,
        max_half_angle: Option<f64>,
        options: Option<Magsac>,
    ) -> Result<Cone3> {
        let min_half_angle = min_half_angle.unwrap_or(EPSILON);
        let max_half_angle = max_half_angle.unwrap_or(FRAC_PI_2 - EPSILON);

        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit = magsac.fit_filtered::<3, Cone3, _>(points, |c| {
            let half_angle = c.half_angle();
            half_angle >= min_half_angle && half_angle <= max_half_angle
        })?;

        // Bound the infinite consensus cone to the axial extent of its inliers (falling back to all
        // points in the unlikely event the inlier set is empty).
        let bound_pts: Vec<Point3> = if fit.inliers.is_empty() {
            points.iter().map(|sp| sp.point).collect()
        } else {
            fit.inliers.iter().map(|&i| points[i].point).collect()
        };
        Ok(bound_axially(&fit.model, &bound_pts))
    }
}

impl ConsensusModel<3> for Cone3 {
    /// A cone is estimated from oriented points: the normals are needed to bootstrap the apex.
    type Point = SurfacePoint3;

    /// Four oriented points, the same minimum [`Cone3::from_fit_oriented`] requires: the smallest
    /// sample that reliably fixes the apex (from the normals) and the axis/half-angle (from the plane
    /// the apex-relative directions lie on).
    const SAMPLE_SIZE: usize = 4;

    fn from_sample(sample: &[SurfacePoint3]) -> Option<Self> {
        let weights = vec![1.0; sample.len()];
        bootstrap_cone(sample, &weights)
    }

    fn residual(&self, point: &SurfacePoint3) -> f64 {
        // Position-only surface distance; the normal plays no part in scoring.
        point_cone_distance(&self.tip, &self.direction, self.half_angle(), &point.point)
    }

    fn refine_weighted(
        points: &[SurfacePoint3],
        weights: &[f64],
        initial: &Cone3,
    ) -> Option<Cone3> {
        // Strip the normals and hand the bare positions to the shared position-only LM engine.
        let pts: Vec<Point3> = points.iter().map(|sp| sp.point).collect();
        ConeFit::refine(&pts, weights, initial)
    }
}

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

/// A weighted least-squares problem for refining an infinite cone (apex, axis, and half-angle). The
/// six parameters are an offset applied relative to a fixed base cone: a translation of the apex
/// (`[0..3]`), two rotations of the axis direction about the base's perpendicular axes (`[3..5]`),
/// and a change in half-angle (`[5]`). Parameterizing relative to the base avoids any orientation
/// singularity, and the jacobian of the surface distances is computed analytically (see
/// [`ConeFit::jacobian`]).
struct ConeFit<'a> {
    points: &'a [Point3],
    weights: &'a [f64],

    /// Fixed base apex that the parameters are measured relative to.
    base_apex: Point3,
    base_direction: UnitVec3,
    base_half_angle: f64,

    /// World-space axes spanning the plane perpendicular to the base axis, used as the rotation axes
    /// for the direction.
    axis_x: UnitVec3,
    axis_y: UnitVec3,

    params: Vector6<f64>,
    residuals: Residuals,
}

impl<'a> ConeFit<'a> {
    fn new(points: &'a [Point3], weights: &'a [f64], base: &Cone3) -> Self {
        let arbitrary = Iso3::from_z_arbitrary_xy(&base.direction, None);
        let axis_x = arbitrary.x();
        let axis_y = arbitrary.y();

        let mut problem = Self {
            points,
            weights,
            base_apex: base.tip,
            base_direction: base.direction,
            base_half_angle: base.half_angle(),
            axis_x,
            axis_y,
            params: Vector6::zeros(),
            residuals: Residuals::zeros(points.len()),
        };
        problem.set_params(&Vector6::zeros());
        problem
    }

    /// Refine `initial` against the weighted geometric surface residuals of `points` with a single
    /// Levenberg-Marquardt solve, returning the optimized infinite cone (with its extent unset - the
    /// `height`/`radius` merely encode the half-angle) or `None` if the solve fails. This is the
    /// shared entry point for both fitting constructors and, later, the consensus refinement.
    fn refine(points: &[Point3], weights: &[f64], initial: &Cone3) -> Option<Cone3> {
        let problem = ConeFit::new(points, weights, initial);
        let (result, report) = LevenbergMarquardt::new().minimize(problem);
        report.termination.was_successful().then(|| {
            let (apex, direction, half_angle) = result.cone_params(&result.params);
            Cone3::new(apex, direction, half_angle.cos(), half_angle.sin())
        })
    }

    /// Reconstruct the `(apex, direction, half_angle)` of the cone described by an offset from the
    /// base. The half-angle is clamped to a proper open cone range.
    fn cone_params(&self, p: &Vector6<f64>) -> (Point3, UnitVec3, f64) {
        let apex = self.base_apex + Vector3::new(p[0], p[1], p[2]);
        let rot_x = UnitQuaternion::from_axis_angle(&self.axis_x, p[3]);
        let rot_y = UnitQuaternion::from_axis_angle(&self.axis_y, p[4]);
        let direction = UnitVec3::new_normalize(rot_y * (rot_x * self.base_direction.into_inner()));
        let half_angle = (self.base_half_angle + p[5]).clamp(EPSILON, FRAC_PI_2 - EPSILON);
        (apex, direction, half_angle)
    }
}

impl LeastSquaresProblem<f64, Dyn, U6> for ConeFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params = *x;
        let (apex, direction, half_angle) = self.cone_params(&self.params);
        for i in 0..self.points.len() {
            self.residuals[i] = self.weights[i]
                * point_cone_distance(&apex, &direction, half_angle, &self.points[i]);
        }
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        Some(self.residuals.clone())
    }

    /// Analytic jacobian of the weighted surface residuals with respect to the six offset
    /// parameters. Writing `v = point - apex`, `q = v·direction`, `radial = v - direction*q` with
    /// `rho = |radial|` and `u = radial/rho`: when the perpendicular foot on the slant line falls at
    /// or past the apex (see [`point_cone_distance`]), the residual is
    /// `D = rho·cos(alpha) - q·sin(alpha)`, whose gradient with respect to the apex (as a free point)
    /// is `sin(alpha)*direction - cos(alpha)*u`, and with respect to the direction (treated as a free
    /// vector) is `-(sin(alpha)*v + cos(alpha)*q*u)`. The latter is chained through `dn/dp3` and
    /// `dn/dp4`, the derivatives of the direction with respect to its two rotation parameters;
    /// because those rotations are about the fixed world-space `axis_x`/`axis_y`,
    /// `d/dtheta [R(theta, a) * v] = a x (R(theta, a) * v)` gives those derivatives in closed form at
    /// any iterate, not just at the base cone. The derivative with respect to the half-angle is
    /// `-(rho·sin(alpha) + q·cos(alpha))`. Otherwise (the foot falls behind the apex) the residual is
    /// `D = |v|`, which depends only on the apex, with gradient `-v/|v|`; the direction and
    /// half-angle columns are zero there.
    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let (apex, direction, half_angle) = self.cone_params(&self.params);
        let n = direction.into_inner();
        let (s, c) = half_angle.sin_cos();
        let axis_x = self.axis_x.into_inner();
        let axis_y = self.axis_y.into_inner();

        let rot_x = UnitQuaternion::from_axis_angle(&self.axis_x, self.params[3]);
        let n1 = rot_x * self.base_direction.into_inner();
        let rot_y = UnitQuaternion::from_axis_angle(&self.axis_y, self.params[4]);
        let dn_dp3 = rot_y * axis_x.cross(&n1);
        let dn_dp4 = axis_y.cross(&n);

        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());

        for (i, point) in self.points.iter().enumerate() {
            let w = self.weights[i];
            let v = point - apex;
            let q = v.dot(&n);
            let radial = v - n * q;
            let rho = radial.norm();

            // The radial direction is undefined right on the axis; leave that row's translation and
            // rotation sensitivity at zero rather than dividing by zero.
            let u = if rho > 1e-12 {
                radial / rho
            } else {
                Vector3::zeros()
            };

            // When the perpendicular foot on the slant line falls behind the apex, the residual is
            // the plain distance to the apex instead (see `point_cone_distance`); the direction and
            // half-angle play no part there, only the apex position does.
            let t = q * c + rho * s;
            let (grad_apex, grad_n, grad_half_angle) = if t < 0.0 {
                let d = v.norm();
                let grad_apex = if d > 1e-12 { -v / d } else { Vector3::zeros() };
                (grad_apex, Vector3::zeros(), 0.0)
            } else {
                (n * s - u * c, -(v * s + u * (c * q)), -(rho * s + q * c))
            };

            jac[(i, 0)] = w * grad_apex.x;
            jac[(i, 1)] = w * grad_apex.y;
            jac[(i, 2)] = w * grad_apex.z;
            jac[(i, 3)] = w * grad_n.dot(&dn_dp3);
            jac[(i, 4)] = w * grad_n.dot(&dn_dp4);
            jac[(i, 5)] = w * grad_half_angle;
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom3::Cone3;
    use crate::{Point3, UnitVec3, Vector3};
    use approx::assert_relative_eq;
    use levenberg_marquardt::LeastSquaresProblem;
    use parry3d_f64::na::Vector6;
    use std::f64::consts::TAU;

    fn tilted() -> Cone3 {
        let tip = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::new(0.3, -0.5, 1.0));
        Cone3::new(tip, direction, 6.0, 3.0)
    }

    /// Deterministically sample `n` oriented surface points on `cone`, spread over its full height
    /// (excluding the apex itself) and around the full azimuth (via the golden angle) so the normals
    /// span 3D. When `outward` is false the normals point into the solid (inner-facing surface).
    fn sample_surface(cone: &Cone3, n: usize, outward: bool) -> Vec<SurfacePoint3> {
        let arbitrary = Iso3::from_z_arbitrary_xy(&cone.direction, None);
        let u = arbitrary.x().into_inner();
        let v = arbitrary.y().into_inner();
        let d = cone.direction.into_inner();
        let alpha = cone.half_angle();
        (0..n)
            .map(|i| {
                let f = (i + 1) as f64 / n as f64; // avoid f == 0 (the apex)
                let q = cone.height * f;
                let ang = TAU * i as f64 * 0.6180339887;
                let radial = u * ang.cos() + v * ang.sin();
                let point = cone.tip + d * q + radial * (q * alpha.tan());
                // Outward normal on a cone: cos α along the radial, sin α back toward the apex.
                let normal =
                    (radial * alpha.cos() - d * alpha.sin()) * if outward { 1.0 } else { -1.0 };
                SurfacePoint3::new_normalize(point, normal)
            })
            .collect()
    }

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

    #[test]
    fn point_cone_distance_on_axis_behind_apex_is_distance_to_apex() {
        // A point on the axis, 10 units behind the apex: no slant line projects there (t < 0 for
        // every azimuth), so the closest feature is the apex itself and the distance must equal the
        // plain Euclidean distance to it, not the unclamped `rho*cos(alpha) - q*sin(alpha)` line
        // formula (which would badly underestimate it: at rho = 0, q = -10, alpha = 30 degrees, that
        // formula gives `10*sin(30deg) = 5`, half the true distance of `10`).
        let apex = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0));
        let half_angle = std::f64::consts::FRAC_PI_6;
        let behind = apex - direction.into_inner() * 10.0;

        let d = point_cone_distance(&apex, &direction, half_angle, &behind);
        assert_relative_eq!(d, 10.0, epsilon = 1e-9);
    }

    #[test]
    fn point_cone_distance_far_off_axis_behind_apex_uses_apex_too() {
        // Far off-axis but still behind the apex (t < 0): still clamps to the apex distance rather
        // than the unbounded line formula.
        let apex = Point3::new(0.0, 0.0, 0.0);
        let direction = UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0));
        let half_angle = std::f64::consts::FRAC_PI_4; // 45 degrees, tan = 1
        // q = -1, rho = 0.5: t = q*cos+rho*sin = -0.707 + 0.354 = -0.354 < 0, so still clamped.
        let point = Point3::new(0.5, 0.0, -1.0);

        let d = point_cone_distance(&apex, &direction, half_angle, &point);
        let expected = (point - apex).norm();
        assert_relative_eq!(d, expected, epsilon = 1e-9);
    }

    #[test]
    fn point_cone_distance_on_surface_ahead_of_apex_matches_unclamped_formula() {
        // Ahead of the apex and exactly on the lateral surface: the perpendicular foot is on the
        // ray (t >= 0), so the residual should be (near) zero, matching the unclamped case.
        let cone = tilted();
        let points = sample_surface(&cone, 20, true);
        for sp in &points {
            let d = point_cone_distance(&cone.tip, &cone.direction, cone.half_angle(), &sp.point);
            assert_relative_eq!(d, 0.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn jacobian_matches_finite_differences() {
        let cone = tilted();
        let surface = sample_surface(&cone, 30, true);
        let points: Vec<Point3> = surface.iter().map(|sp| sp.point).collect();
        let weights = vec![1.0; points.len()];

        // Perturb away from the base cone so the jacobian is checked away from `p = 0`, where the
        // fixed-axis rotation derivatives are exercised at a non-trivial angle.
        let mut problem = ConeFit::new(&points, &weights, &cone);
        let probe = Vector6::new(0.05, -0.03, 0.02, 0.2, -0.15, 0.05);
        problem.set_params(&probe);

        let analytic = problem.jacobian().expect("analytic jacobian");
        let base_residuals = problem.residuals().expect("residuals");

        const DELTA: f64 = 1e-6;
        for k in 0..6 {
            let mut perturbed = probe;
            perturbed[k] += DELTA;
            problem.set_params(&perturbed);
            let perturbed_residuals = problem.residuals().expect("residuals");
            problem.set_params(&probe);

            for i in 0..points.len() {
                let numeric = (perturbed_residuals[i] - base_residuals[i]) / DELTA;
                assert_relative_eq!(analytic[(i, k)], numeric, epsilon = 1e-4);
            }
        }
    }

    #[test]
    fn jacobian_matches_finite_differences_behind_apex() {
        // Points solidly behind the apex (large negative `t` even after perturbing the params by
        // `DELTA`), exercising the apex-clamp branch of both the residual and its jacobian.
        let cone = tilted();
        let d = cone.direction.into_inner();
        let arbitrary = Iso3::from_z_arbitrary_xy(&cone.direction, None);
        let u = arbitrary.x().into_inner();
        let points: Vec<Point3> = (0..10)
            .map(|i| cone.tip - d * (5.0 + i as f64) + u * (0.1 * i as f64))
            .collect();
        let weights = vec![1.0; points.len()];

        let mut problem = ConeFit::new(&points, &weights, &cone);
        let probe = Vector6::new(0.05, -0.03, 0.02, 0.2, -0.15, 0.05);
        problem.set_params(&probe);

        let analytic = problem.jacobian().expect("analytic jacobian");
        let base_residuals = problem.residuals().expect("residuals");

        const DELTA: f64 = 1e-6;
        for k in 0..6 {
            let mut perturbed = probe;
            perturbed[k] += DELTA;
            problem.set_params(&perturbed);
            let perturbed_residuals = problem.residuals().expect("residuals");
            problem.set_params(&probe);

            for i in 0..points.len() {
                let numeric = (perturbed_residuals[i] - base_residuals[i]) / DELTA;
                assert_relative_eq!(analytic[(i, k)], numeric, epsilon = 1e-4);
            }
        }
    }

    /// Assert two cones describe the same geometry (apex, axis direction, half-angle) and extent.
    fn assert_cone_matches(result: &Cone3, expected: &Cone3, eps: f64) {
        assert_relative_eq!(result.tip, expected.tip, epsilon = eps);
        assert_relative_eq!(result.r(), expected.r(), epsilon = eps);
        assert_relative_eq!(result.height, expected.height, epsilon = eps);
        assert_relative_eq!(result.half_angle(), expected.half_angle(), epsilon = eps);
        // The axis points from tip toward base, so its sign is determined (not just parallel).
        let dot = result
            .direction
            .into_inner()
            .dot(&expected.direction.into_inner());
        assert_relative_eq!(dot, 1.0, epsilon = eps);
    }

    // from_fit_oriented tests ─────────────────────────────────────────────────

    #[test]
    fn from_fit_oriented_recovers_clean_cone() {
        let expected = tilted();
        let points = sample_surface(&expected, 60, true);
        let fit = Cone3::from_fit_oriented(&points, None).unwrap();
        assert_cone_matches(&fit, &expected, 1e-6);
    }

    #[test]
    fn from_fit_oriented_inner_facing_matches_outer_facing() {
        let expected = tilted();
        let outer = Cone3::from_fit_oriented(&sample_surface(&expected, 60, true), None).unwrap();
        let inner = Cone3::from_fit_oriented(&sample_surface(&expected, 60, false), None).unwrap();

        assert_cone_matches(&inner, &expected, 1e-6);
        assert_relative_eq!(inner.tip, outer.tip, epsilon = 1e-9);
        assert_relative_eq!(
            inner.direction.into_inner(),
            outer.direction.into_inner(),
            epsilon = 1e-9
        );
        assert_relative_eq!(inner.half_angle(), outer.half_angle(), epsilon = 1e-9);
        assert_relative_eq!(inner.height, outer.height, epsilon = 1e-9);
        assert_relative_eq!(inner.r(), outer.r(), epsilon = 1e-9);
    }

    #[test]
    fn from_fit_oriented_mixed_normal_facing_still_fits() {
        // The apex system is invariant under flipping any single normal, so an arbitrary mix of
        // inward- and outward-facing normals must still recover the cone.
        let expected = tilted();
        let mut points = sample_surface(&expected, 60, true);
        for (i, sp) in points.iter_mut().enumerate() {
            if i % 2 == 0 {
                *sp = SurfacePoint3::new(sp.point, -sp.normal);
            }
        }
        let fit = Cone3::from_fit_oriented(&points, None).unwrap();
        assert_cone_matches(&fit, &expected, 1e-6);
    }

    #[test]
    fn from_fit_oriented_with_noise_is_close() {
        let mut rg = RandomGeometry3::from_seed(8);
        let expected = tilted();
        let points = with_noise(&mut rg, &sample_surface(&expected, 500, true), 0.01);
        let fit = Cone3::from_fit_oriented(&points, None).unwrap();
        assert_cone_matches(&fit, &expected, 0.1);
    }

    #[test]
    fn from_fit_oriented_uniform_weights_match_unweighted() {
        let mut rg = RandomGeometry3::from_seed(21);
        let expected = tilted();
        let points = with_noise(&mut rg, &sample_surface(&expected, 60, true), 0.02);

        let unweighted = Cone3::from_fit_oriented(&points, None).unwrap();
        let weights = vec![1.0; points.len()];
        let weighted = Cone3::from_fit_oriented(&points, Some(&weights)).unwrap();
        assert_relative_eq!(weighted.tip, unweighted.tip, epsilon = 1e-9);
        assert_relative_eq!(
            weighted.direction.into_inner(),
            unweighted.direction.into_inner(),
            epsilon = 1e-9
        );
        assert_relative_eq!(
            weighted.half_angle(),
            unweighted.half_angle(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn from_fit_oriented_too_few_points_is_error() {
        let expected = tilted();
        let points = sample_surface(&expected, 3, true);
        assert!(Cone3::from_fit_oriented(&points, None).is_err());
    }

    #[test]
    fn from_fit_oriented_cylindrical_sample_is_error() {
        // Points whose normals are all perpendicular to a common axis (a cylinder) do not determine
        // a cone apex; the bootstrap must reject them rather than return a garbage cone.
        let axis = UnitVec3::new_normalize(Vector3::new(0.3, -0.5, 1.0));
        let arbitrary = Iso3::from_z_arbitrary_xy(&axis, None);
        let u = arbitrary.x().into_inner();
        let v = arbitrary.y().into_inner();
        let base = Point3::new(1.0, 2.0, 3.0);
        let points: Vec<SurfacePoint3> = (0..40)
            .map(|i| {
                let t = i as f64 * 0.25;
                let ang = TAU * i as f64 * 0.6180339887;
                let radial = u * ang.cos() + v * ang.sin();
                let point = base + axis.into_inner() * t + radial * 2.0;
                // Cylindrical (radial) normals -> all perpendicular to `axis`.
                SurfacePoint3::new_normalize(point, radial)
            })
            .collect();
        assert!(Cone3::from_fit_oriented(&points, None).is_err());
    }

    // from_fit (guess-based) tests ─────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_clean_cone_from_guess() {
        let expected = tilted();
        let surface = sample_surface(&expected, 60, true);
        let points: Vec<Point3> = surface.iter().map(|sp| sp.point).collect();

        let guess = Cone3::new(
            expected.tip + Vector3::new(0.2, -0.15, 0.1),
            UnitVec3::new_normalize(
                expected.direction.into_inner() + Vector3::new(0.05, 0.05, 0.0),
            ),
            expected.height + 0.5,
            expected.r() + 0.3,
        );

        let fit = Cone3::from_fit(&points, &guess, None).unwrap();
        assert_cone_matches(&fit, &expected, 1e-6);
    }

    #[test]
    fn from_fit_with_noise_is_close() {
        let mut rg = RandomGeometry3::from_seed(13);
        let expected = tilted();
        let surface = with_noise(&mut rg, &sample_surface(&expected, 500, true), 0.01);
        let points: Vec<Point3> = surface.iter().map(|sp| sp.point).collect();

        let guess = Cone3::new(
            expected.tip + Vector3::new(0.1, 0.1, -0.1),
            expected.direction,
            expected.height + 0.3,
            expected.r() + 0.2,
        );
        let fit = Cone3::from_fit(&points, &guess, None).unwrap();
        assert_cone_matches(&fit, &expected, 0.1);
    }

    #[test]
    fn from_fit_too_few_points_is_error() {
        let guess = tilted();
        let points = [Point3::new(0.0, 0.0, 0.0)];
        assert!(Cone3::from_fit(&points, &guess, None).is_err());
    }

    // from_consensus tests ───────────────────────────────────────────────────

    #[test]
    fn from_consensus_recovers_clean_cone() {
        let expected = tilted();
        let points = sample_surface(&expected, 120, true);
        let options = Magsac {
            sigma_max: 0.01,
            max_iterations: Some(500),
            refinement_steps: 5,
            confidence: 0.99,
            seed: Some(17),
        };
        let fit = Cone3::from_consensus(&points, 0.01, None, None, Some(options)).unwrap();
        assert_cone_matches(&fit, &expected, 5e-3);
    }

    #[test]
    fn from_consensus_rejects_outliers() {
        let mut rg = RandomGeometry3::from_seed(29);
        let expected = tilted();
        let inliers = with_noise(&mut rg, &sample_surface(&expected, 150, true), 0.01);
        let inlier_count = inliers.len();
        let mut points = inliers.clone();

        // Gross outliers spread along the same axial range as the inliers but well outside the true
        // lateral surface radially. A fixed-offset gaussian cluster (as used for the cylinder's
        // equivalent test) risks landing close enough to the surface for some draws to slip under
        // sigma_max by chance, since the cone's residual gap at a single point is easily comparable
        // to the noise scale; anchoring the radial offset to a large multiple of the true radius at
        // each axial position keeps every outlier at a guaranteed distance from the surface instead.
        let arbitrary = Iso3::from_z_arbitrary_xy(&expected.direction, None);
        let u = arbitrary.x().into_inner();
        let v = arbitrary.y().into_inner();
        let d = expected.direction.into_inner();
        for i in 0..50 {
            let f = i as f64 / 49.0;
            let q = expected.height * f;
            let r_at_q = q * expected.half_angle().tan();
            let ang = TAU * i as f64 * 0.6180339887 + 1.0;
            let radial = u * ang.cos() + v * ang.sin();
            let point = expected.tip + d * q + radial * (5.0 * r_at_q + 1.0);
            points.push(SurfacePoint3::new_normalize(
                point,
                rg.unit_vec().into_inner(),
            ));
        }

        let magsac = Magsac {
            sigma_max: 0.03,
            max_iterations: Some(500),
            refinement_steps: 5,
            confidence: 0.99,
            seed: Some(17),
        };
        let fit = magsac
            .fit_filtered::<3, Cone3, _>(&points, |_| true)
            .unwrap();

        // The apex, axis, and half-angle match (height/radius are not checked: the raw consensus
        // model's extent fields merely encode the half-angle until bounded).
        assert_relative_eq!(
            fit.model.half_angle(),
            expected.half_angle(),
            epsilon = 1e-2
        );
        let dot = fit
            .model
            .direction
            .into_inner()
            .dot(&expected.direction.into_inner());
        assert_relative_eq!(dot, 1.0, epsilon = 1e-2);
        assert!((fit.model.tip - expected.tip).norm() < 0.05);

        // No outlier is classified as an inlier, and nearly all true inliers are recovered.
        assert!(fit.inliers.iter().all(|&i| i < inlier_count));
        assert!(fit.inliers.len() > inlier_count * 9 / 10);
    }

    #[test]
    fn from_consensus_too_few_points_is_error() {
        // The minimal sample is four oriented points, so three cannot seed a candidate.
        let expected = tilted();
        let points = sample_surface(&expected, 3, true);
        assert!(Cone3::from_consensus(&points, 0.01, None, None, None).is_err());
    }
}
