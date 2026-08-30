//! Finding the plane that cuts "across" a swept feature of a mesh.
//!
//! A pipe, a conduit, a cone, a machined groove, or a tapered slab is a section swept along some
//! path. At a chosen point on such a feature, one plane produces a section that goes across the
//! sweep rather than obliquely through it. This plane produces a circle on a cylinder or cone,
//! follows the cutter profile on a groove, and makes equal angles with the two faces of a wedge.
//! It is the 3D counterpart of the airfoil thickness line drawn normal to the medial axis.
//!
//! # The criterion
//!
//! Cut the mesh with a trial plane through the point with unit normal `n`. Each cut triangle `i`
//! contributes a segment of length `L_i` with face normal `N_i`, a tilt `phi_i` of that face out
//! of the plane (`sin(phi_i) = N_i . n`), and an in-plane outward normal `m_i`, the normalized
//! in-plane part of `N_i`. If the plane is now slid a distance `dz` along `n`, the segment from
//! face `i` moves outward by `-dz tan(phi_i)` along `m_i`: `-tan(phi)` is the rate at which the
//! section boundary grows as the plane advances through the sweep.
//!
//! That growth field is the sum of the genuine change in the section's shape and a rigid in-plane
//! translation of the whole section, and the translation is precisely what appears when `n` is
//! not the sweep direction: an oblique plane sees the section slide sideways as it advances. The
//! across plane is the one where the translation component vanishes. Projecting the growth field
//! onto rigid translations in the arclength-weighted least-squares sense makes that condition
//!
//! ```text
//! r(n) = sum_i L_i tan(phi_i) m_i = 0
//! ```
//!
//! two equations in the two degrees of freedom of `n`. On a cylinder, every `tan(phi)` is zero.
//! On a cone, `tan(phi)` is uniform and the outward normals of a closed loop sum to zero, wherever
//! the point sits. On a wedge, the two converging faces have equal tilts with opposite `m`, and
//! any rotation of the plane breaks the balance. On a one-sided taper (one face flat, one sloped)
//! the root is the plane normal to the bisector of the two faces, which is the medial-axis
//! behavior intended here. Bending of the sweep path contributes nothing in-plane, so
//! curved pipes and grooves need no special treatment.
//!
//! `tan(phi)` is the weight applied up to a knee angle (70 degrees by default), and from there
//! the weight is tapered to zero at a right angle as `sin(phi) cos(phi)`, scaled to meet `tan`
//! at the knee. Below the knee the criterion is untouched. The taper matters because the section
//! of a real part grazes faces nearly parallel to the plane (a platform beside the feature, a
//! fillet, an end cap) that the sweep model says nothing about; under an unbounded `tan`, one
//! such face outweighs the whole section, and its outward direction `m` flips as it passes
//! through the parallel position. Tapering to zero there keeps the residual bounded and
//! continuous. Three alternatives were tried and rejected. Weighting by `cos^2(phi)` throughout,
//! which makes `n` a principal axis of the section's normal tensor, decays past 45 degrees of
//! tilt, making the far roots attractive; on a cone, all three principal axes are roots.
//! Weighting by `sin(phi)`, bounded by one, gives a wide cone a pull of `2 cos(a) sin(theta)`
//! for a tilt `theta` of the plane, which vanishes as the half-angle `a` grows, and its basin
//! turned out to be small and erratic where the basin for `tan` reaches 45 degrees on a cylinder.
//! A knee at 60 degrees also puts the faces of a wide cone into the taper, where the ordering of
//! the weights reverses and sends the solver in the wrong direction.
//!
//! # The band
//!
//! The balance is not evaluated on the section curve itself but on the band of surface within a
//! half-width `h` of the plane. Each face is weighted by its area inside the band times
//! `cos(phi)`. For a face that crosses the whole band that weight is `2h L_i`, so the criterion
//! is the one above; only its smoothness differs. A segment's length is a piecewise
//! linear function of the plane's orientation, with a slope jump every time the plane crosses a
//! mesh vertex. A small, steep face, such as one caused by roughness in a scanned surface, creates
//! a large jump.
//! Those jumps create spurious local minima in the squared residual near the root, where a
//! Levenberg-Marquardt solve stalls. The band area is the integral of the cut length across the
//! band. It is piecewise quadratic with a continuous derivative, so the solve sees a smooth
//! function.
//! The half-width defaults to a hundredth of the section's length, a few vertex spacings on a
//! typical mesh; it sets the scale over which the roughness is averaged and has no effect on the
//! root of a clean sweep.
//!
//! The coverage matrix `M = sum_i L_i m_i m_i^T` says whether the two equations are independent:
//! it is full rank whenever the section's outward normals span the plane. A closed loop, a half
//! loop, and a U-shaped groove profile all meet this condition. The matrix has rank one for a
//! single flat face, which has no sweep direction to find.
//!
//! # Related work
//!
//! The closest published relative is the rotational symmetry axis of Tagliasacchi, Zhang, and
//! Cohen-Or (2009), a cutting plane whose normal minimizes the variance of the angles to nearby
//! point normals. Variance of tilt is not balance of tilt: the flat side faces of a wedge pull
//! the variance minimum away from the equal-angle plane. The slippage analysis of Gelfand and
//! Guibas (2004) recovers the extrusion direction as the null space of the normal tensor, which
//! is exact only for a constant section. Minimum section area is not usable as a definition
//! either: on a cone it is an area maximum for half-angles over 30 degrees, and with the point on
//! the surface the perpendicular plane is not an area extremum at all.

use super::Mesh3;
use super::section::{TriIntr, chain_segments};
use crate::common::PCoords;
use crate::geom3::IsoExtensions3;
use crate::{IndexMask, Iso3, Plane3, Point3, Result, UnitVec3, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Matrix, Matrix2, Matrix3, Owned, U2, UnitQuaternion, Vector, Vector2};
use parry3d_f64::partitioning::TraversalAction;
use std::f64::consts::FRAC_PI_2;

/// The step in radians used for the central-difference Jacobian. The band residual has a
/// continuous derivative, so the step only has to be small against the scale on which that
/// derivative varies: the band's half-width divided by the section's extent.
const DELTA: f64 = 1e-4;

/// Below this coverage, the two balance equations are not independent, so the direction cannot
/// be solved for.
const MIN_COVERAGE: f64 = 1e-3;

/// Below this sensitivity, the residual does not respond to some tilt of the plane, so the solve
/// has not determined the direction even if the residual is zero.
const MIN_SENSITIVITY: f64 = 1e-3;

/// The default band half-width as a fraction of the section's length at the guess.
const DEFAULT_BAND_FRACTION: f64 = 0.01;

/// A band face takes part in the balance if it lies within this many half-widths of the nearest
/// run of the central section. A face tilted steeply out of the plane extends farther from the
/// section curve than the band is thick; this margin allows for that extension.
const SCOPE_REACH: f64 = 4.0;

/// Options for [`Mesh3::transverse_plane`].
#[derive(Debug, Clone, Copy)]
pub struct TransverseOptions {
    /// The face tilt out of the plane, in radians, beyond which a segment's weight is tapered
    /// from `tan(phi)` to zero at a right angle, keeping the residual continuous as a face passes
    /// through the parallel position. Below this angle, the weight increases with tilt, which the
    /// solver relies on, so it should stay above the tilts the feature's own faces can reach.
    /// Defaults to 70 degrees. Must lie strictly between zero and a right angle.
    pub taper_tilt: f64,

    /// The half-width of the surface band around the plane over which the balance is evaluated,
    /// in the mesh's units. `None` uses a hundredth of the section's length at the guess. See
    /// the module documentation for what the band does.
    pub band: Option<f64>,

    /// The maximum number of residual evaluations before the solve gives up. Defaults to 100.
    pub max_evaluations: usize,

    /// The change in the plane normal, in radians, below which the solve is considered converged.
    /// Defaults to `1e-9`.
    pub tol: f64,
}

impl Default for TransverseOptions {
    fn default() -> Self {
        Self {
            taper_tilt: 70.0_f64.to_radians(),
            band: None,
            max_evaluations: 100,
            tol: 1e-9,
        }
    }
}

impl TransverseOptions {
    fn validate(&self) -> Result<()> {
        if !(self.taper_tilt > 0.0 && self.taper_tilt < FRAC_PI_2) {
            return Err(format!(
                "taper_tilt must lie strictly between 0 and pi/2 radians, got {}",
                self.taper_tilt
            )
            .into());
        }
        if let Some(band) = self.band
            && !(band > 0.0 && band.is_finite())
        {
            return Err(format!("band must be positive and finite, got {band}").into());
        }
        if self.max_evaluations == 0 {
            return Err("max_evaluations must be at least 1".into());
        }
        if !(self.tol > 0.0 && self.tol.is_finite()) {
            return Err(format!("tol must be positive and finite, got {}", self.tol).into());
        }
        Ok(())
    }
}

/// The result of [`Mesh3::transverse_plane`]: the plane found, along with diagnostics showing how
/// well the balance condition was met and how well the plane was determined.
///
/// The condition is that the section's tilt-weighted outward normals balance,
/// `sum_i L_i sin(phi_i) m_i = 0`, where for each cut triangle `L_i` is its segment's length,
/// `phi_i` the tilt of the face out of the plane, and `m_i` the in-plane outward normal, with
/// the weight tapered to zero for faces beyond [`TransverseOptions::taper_tilt`]. Sliding the
/// plane along its normal moves each segment outward at the rate `-tan(phi_i)`, so up to a
/// `cos(phi_i)` weighting this is the statement that the section does not translate sideways as
/// the plane advances along the sweep, which is what makes the plane go across the feature
/// rather than obliquely through it.
#[derive(Debug, Clone)]
pub struct TransversePlane {
    /// The plane through the query point whose section goes across the feature. Its normal is
    /// oriented to point the same way as the guess.
    pub plane: Plane3,

    /// The number of times the section was evaluated, including the Jacobian evaluations that
    /// measure the sensitivity at the solution.
    pub evaluations: usize,

    /// The magnitude of the balance residual at the solution, divided by the section's total
    /// length. Dimensionless; it is the length-weighted mean of the tapered `tan(phi)` left
    /// unbalanced, so a well-converged solve on a clean sweep sits near zero.
    pub residual: f64,

    /// The smaller eigenvalue of the coverage matrix divided by the section's total length.
    /// Ranges from `0` (all outward normals parallel, the direction is undetermined) to `1/2`
    /// (a full closed loop, or any section whose normals are spread evenly around the plane).
    pub coverage: f64,

    /// The smaller singular value of the residual's Jacobian with respect to the plane's tilt at
    /// the solution: how much the normalized residual changes per radian of the tilt it is least
    /// sensitive to. For a constant-section sweep it equals the coverage; a changing section
    /// moves it away from that.
    pub sensitivity: f64,

    /// The number of band faces that took part in the final evaluation.
    pub face_count: usize,

    /// The band half-width that was used.
    pub band: f64,
}

impl Mesh3 {
    /// Finds the plane through `point` whose section cuts across the swept feature the point lies
    /// on: the plane that gives a circle on a cylinder or cone, follows the cutter profile on a
    /// groove, and makes equal angles with the two faces of a wedge. The criterion—the balance of
    /// the faces' tilts around the section—is described with [`TransversePlane`].
    ///
    /// Only the connected run of the section that passes nearest the point takes part. A plane
    /// through one feature usually cuts others elsewhere on the mesh (a plane through a pipe
    /// elbow cuts the pipe twice), and those have no bearing on the section at the point. They
    /// also cannot simply be added: the far side of a circumferential groove is the same feature
    /// swept the opposite way, and summing its imbalance with the near side's cancels the two
    /// where the plane rotates about the shaft axis, leaving a line of false roots. Balancing
    /// each run separately would be sound, and is the extension to make if fragmented sections
    /// (a scanned pipe with holes) turn out to need every piece.
    ///
    /// The plane is found by a Levenberg-Marquardt solve over the two degrees of freedom of its
    /// normal, starting from `guess`. Each evaluation makes one pass over the band of faces around
    /// the trial plane. The solve converges to the nearest balanced plane, and a sweep has more
    /// than one: a plane containing a cylinder's axis is also balanced. The guess therefore has to
    /// be closer to the across plane than to any other, and the trial plane it defines has to go
    /// across the feature in the first place: on a cylinder anything within 45 degrees of the
    /// axis converges, but on a cone of half-angle `a` a plane more than `90 - a` degrees from
    /// the axis no longer cuts a closed section (the conic becomes a hyperbola), and a plane
    /// whose section runs out the end of a finite feature has the same problem. A guess taken
    /// from the geometry (a fitted axis, a principal curvature direction) is normally far inside
    /// those limits.
    ///
    /// # Arguments
    ///
    /// * `point`: the point the plane passes through. It may be on the surface or inside the
    ///   feature (on a pipe's axis, say); the balance condition does not depend on where in the
    ///   section it sits.
    /// * `guess`: an initial direction for the plane normal, roughly along the sweep.
    /// * `faces`: an optional mask limiting the cut to a subset of the faces, which must be the
    ///   same length as the mesh has faces. Use it to isolate the feature from the rest of the
    ///   part: the plane across a groove must not be balanced against the part around it.
    /// * `options`: solve options, or `None` for [`TransverseOptions::default`].
    ///
    /// returns: Result<TransversePlane, Box<dyn Error, Global>>
    ///
    /// # Failure
    ///
    /// An error is returned if the options are invalid, if the plane through the guess misses
    /// the selected faces, if the section's outward normals do not span the plane (a single flat
    /// face, for instance) at either the start or the end of the solve, if the solve does not
    /// converge within the allowed evaluations, or if the residual at the solution is insensitive
    /// to some tilt of the plane.
    pub fn transverse_plane(
        &self,
        point: &impl PCoords<3>,
        guess: &UnitVec3,
        faces: Option<&IndexMask>,
        options: Option<TransverseOptions>,
    ) -> Result<TransversePlane> {
        let options = options.unwrap_or_default();
        options.validate()?;
        let point = Point3::from(point.coords());

        let problem = TransverseFit::new(self, point, guess, faces, options)?;
        let (result, report) = LevenbergMarquardt::new()
            .with_patience(options.max_evaluations)
            .with_xtol(options.tol)
            .minimize(problem);

        if !report.termination.was_successful() {
            return Err(format!(
                "the transverse plane solve did not converge: {:?}",
                report.termination
            )
            .into());
        }

        let balance = result
            .balance
            .as_ref()
            .ok_or("the transverse plane solve ended without a section")?;
        let (e1, e2) = result.in_plane_axes();
        let coverage = balance.coverage(&e1, &e2);
        if coverage < MIN_COVERAGE {
            return Err(underdetermined(coverage));
        }

        let jacobian = result
            .jacobian()
            .ok_or("the section could not be evaluated around the solution")?;
        let sensitivity = jacobian
            .singular_values()
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        if sensitivity < MIN_SENSITIVITY {
            return Err(format!(
                "the balance residual is insensitive to a tilt of the plane at the solution \
                 (sensitivity {sensitivity:.2e}), so the transverse direction is undetermined"
            )
            .into());
        }

        let mut normal = result.normal();
        if normal.dot(guess) < 0.0 {
            normal = -normal;
        }

        Ok(TransversePlane {
            plane: Plane3::from_point_normal(&point, &normal),
            evaluations: report.number_of_evaluations + 4,
            residual: balance.residual().norm(),
            coverage,
            sensitivity,
            face_count: balance.count,
            band: result.band,
        })
    }
}

fn underdetermined(coverage: f64) -> Box<dyn std::error::Error> {
    format!(
        "the section's outward normals do not span the plane (coverage {coverage:.2e}), so the \
         transverse direction is undetermined"
    )
    .into()
}

/// The accumulated balance quantities of one band.
struct Balance {
    /// `sum_i (A_i cos(phi_i) / 2h) w(phi_i) m_i`, lying in the plane.
    r: Vector3,
    /// `sum_i (A_i cos(phi_i) / 2h) m_i m_i^T`, the coverage matrix, lying in the plane.
    m: Matrix3<f64>,
    /// `sum_i A_i cos(phi_i) / 2h`, the band's measure of section length.
    sum_l: f64,
    /// The number of faces that took part.
    count: usize,
}

impl Balance {
    /// The balance residual normalized by the section length.
    fn residual(&self) -> Vector3 {
        self.r / self.sum_l
    }

    /// The smaller eigenvalue of the coverage matrix in the given in-plane basis, divided by the
    /// section length.
    fn coverage(&self, e1: &Vector3, e2: &Vector3) -> f64 {
        let a = e1.dot(&(self.m * e1));
        let b = e1.dot(&(self.m * e2));
        let d = e2.dot(&(self.m * e2));
        let mean = 0.5 * (a + d);
        let half_gap = (0.25 * (a - d).powi(2) + b * b).sqrt();
        (mean - half_gap) / self.sum_l
    }
}

/// Evaluates the balance of the band of `mesh` within `band` of the plane through `point` with
/// normal `n`. Returns `None` if the plane cuts none of the selected faces.
fn balance(
    mesh: &Mesh3,
    point: &Point3,
    n: &UnitVec3,
    faces: Option<&IndexMask>,
    band: f64,
    options: &TransverseOptions,
) -> Result<Option<Balance>> {
    let plane = Plane3::from_point_normal(point, n);

    // The section itself defines the scope: the run of the section curve nearest the point.
    let segments = mesh.section_segments(&plane, faces)?;
    if segments.is_empty() {
        return Ok(None);
    }
    let scope: Vec<&TriIntr> = nearest_chain(&segments, point)
        .into_iter()
        .map(|i| &segments[i])
        .collect();
    let reach = SCOPE_REACH * band;

    let n = n.into_inner();
    let cos_knee = options.taper_tilt.cos();
    let mut result = Balance {
        r: Vector3::zeros(),
        m: Matrix3::zeros(),
        sum_l: 0.0,
        count: 0,
    };

    for face in band_faces(mesh, &plane, band, faces) {
        let Some(face_n) = mesh.shape.triangle(face).normal() else {
            continue;
        };
        let (area, centroid) = clipped_area(mesh, face, &plane, band);
        if area <= 0.0 {
            continue;
        }
        if scope
            .iter()
            .all(|seg| segment_distance(seg, &centroid) > reach)
        {
            continue;
        }

        // s = sin(phi), and the in-plane part of the face normal is cos(phi) m.
        let s = face_n.dot(&n);
        let in_plane = face_n.into_inner() - n * s;
        let c = in_plane.norm();
        if c <= f64::EPSILON {
            // A face parallel to the plane has no in-plane outward direction to balance, and its
            // weight is zero anyway.
            continue;
        }
        let m = in_plane / c;
        // tan(phi) up to the knee, then sin(phi) cos(phi) scaled to meet it there.
        let weight = if c >= cos_knee {
            s / c
        } else {
            s * c / (cos_knee * cos_knee)
        };
        // The band area of a face crossing the whole band is 2h L / cos(phi).
        let length = area * c / (2.0 * band);

        result.r += length * weight * m;
        result.m += length * m * m.transpose();
        result.sum_l += length;
        result.count += 1;
    }

    Ok((result.count > 0).then_some(result))
}

/// Uses the BVH to find mesh faces that may have area within `band` of the plane.
fn band_faces(mesh: &Mesh3, plane: &Plane3, band: f64, faces: Option<&IndexMask>) -> Vec<u32> {
    let mut candidates = Vec::new();
    mesh.shape.bvh().traverse(|node| {
        let (lo, hi) = node
            .aabb()
            .vertices()
            .iter()
            .map(|v| plane.signed_distance_to_point(v))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                (lo.min(d), hi.max(d))
            });
        if lo > band || hi < -band {
            return TraversalAction::Prune;
        }
        if let Some(index) = node.leaf_data()
            && faces.is_none_or(|mask| mask.get(index as usize))
        {
            candidates.push(index);
        }
        TraversalAction::Continue
    });
    candidates
}

/// The area of the face within `band` of the plane and the centroid of that area.
fn clipped_area(mesh: &Mesh3, face: u32, plane: &Plane3, band: f64) -> (f64, Point3) {
    let tri = mesh.shape.triangle(face);
    let mut polygon = vec![tri.a, tri.b, tri.c];
    polygon = clip_polygon(&polygon, |p| band - plane.signed_distance_to_point(p));
    polygon = clip_polygon(&polygon, |p| band + plane.signed_distance_to_point(p));
    if polygon.len() < 3 {
        return (0.0, Point3::origin());
    }

    // Fan triangulation from the first vertex; the polygon is convex.
    let mut area_vec = Vector3::zeros();
    let mut centroid = Vector3::zeros();
    for i in 1..polygon.len() - 1 {
        let cross = (polygon[i] - polygon[0]).cross(&(polygon[i + 1] - polygon[0]));
        let a = 0.5 * cross.norm();
        area_vec += cross;
        centroid += a * (polygon[0].coords + polygon[i].coords + polygon[i + 1].coords) / 3.0;
    }
    let area = 0.5 * area_vec.norm();
    if area <= 0.0 {
        return (0.0, Point3::origin());
    }
    (area, Point3::from(centroid / area))
}

/// Sutherland-Hodgman clipping of a convex polygon against the half-space where `inside` is
/// non-negative; `inside` must be affine in the point.
fn clip_polygon(polygon: &[Point3], inside: impl Fn(&Point3) -> f64) -> Vec<Point3> {
    let mut out = Vec::with_capacity(polygon.len() + 2);
    for i in 0..polygon.len() {
        let a = polygon[i];
        let b = polygon[(i + 1) % polygon.len()];
        let da = inside(&a);
        let db = inside(&b);
        if da >= 0.0 {
            out.push(a);
        }
        if (da >= 0.0) != (db >= 0.0) {
            let t = da / (da - db);
            out.push(a + (b - a) * t);
        }
    }
    out
}

/// The indices of the segments in the connected run which passes nearest the point.
fn nearest_chain(segments: &[TriIntr], point: &Point3) -> Vec<usize> {
    let chains = chain_segments(segments);
    let mut best = (f64::INFINITY, 0);
    for (k, chain) in chains.iter().enumerate() {
        for &i in chain.iter() {
            let d = segment_distance(&segments[i], point);
            if d < best.0 {
                best = (d, k);
            }
        }
    }
    chains.into_iter().nth(best.1).unwrap_or_default()
}

fn segment_distance(seg: &TriIntr, point: &Point3) -> f64 {
    let ab = seg.b - seg.a;
    let len2 = ab.norm_squared();
    let t = if len2 > 0.0 {
        ((point - seg.a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (point - (seg.a + ab * t)).norm()
}

/// The least-squares problem between two rotation angles about the in-plane axes of the initial
/// guess and the two in-plane components of the normalized balance residual.
struct TransverseFit<'a> {
    mesh: &'a Mesh3,
    point: Point3,
    faces: Option<&'a IndexMask>,
    options: TransverseOptions,

    base: UnitVec3,
    axis_x: UnitVec3,
    axis_y: UnitVec3,
    band: f64,

    params: Vector2<f64>,
    balance: Option<Balance>,
}

impl<'a> TransverseFit<'a> {
    fn new(
        mesh: &'a Mesh3,
        point: Point3,
        guess: &UnitVec3,
        faces: Option<&'a IndexMask>,
        options: TransverseOptions,
    ) -> Result<Self> {
        let band = match options.band {
            Some(band) => band,
            None => {
                let plane = Plane3::from_point_normal(&point, guess);
                let segments = mesh.section_segments(&plane, faces)?;
                let length: f64 = segments.iter().map(|seg| seg.length()).sum();
                if length <= 0.0 {
                    return Err(
                        "the plane through the guess does not cut any of the selected faces".into(),
                    );
                }
                DEFAULT_BAND_FRACTION * length
            }
        };

        let frame = Iso3::from_z_arbitrary_xy(guess, None);
        let mut problem = Self {
            mesh,
            point,
            faces,
            options,
            base: *guess,
            axis_x: frame.x(),
            axis_y: frame.y(),
            band,
            params: Vector2::zeros(),
            balance: None,
        };

        let balance = problem
            .evaluate(&problem.params)?
            .ok_or("the plane through the guess does not cut any of the selected faces")?;
        let (e1, e2) = problem.in_plane_axes();
        let coverage = balance.coverage(&e1, &e2);
        if coverage < MIN_COVERAGE {
            return Err(underdetermined(coverage));
        }
        problem.balance = Some(balance);

        Ok(problem)
    }

    /// The rotation that carries the base direction and its in-plane axes to the current
    /// parameters.
    fn rotation_at(&self, p: &Vector2<f64>) -> UnitQuaternion<f64> {
        let rot_x = UnitQuaternion::from_axis_angle(&self.axis_x, p[0]);
        let rot_y = UnitQuaternion::from_axis_angle(&self.axis_y, p[1]);
        rot_y * rot_x
    }

    fn normal_at(&self, p: &Vector2<f64>) -> UnitVec3 {
        UnitVec3::new_normalize(self.rotation_at(p) * self.base.into_inner())
    }

    fn normal(&self) -> UnitVec3 {
        self.normal_at(&self.params)
    }

    /// The in-plane axes of the current plane: the base axes carried along by the rotation.
    fn in_plane_axes(&self) -> (Vector3, Vector3) {
        let rot = self.rotation_at(&self.params);
        (
            rot * self.axis_x.into_inner(),
            rot * self.axis_y.into_inner(),
        )
    }

    fn evaluate(&self, p: &Vector2<f64>) -> Result<Option<Balance>> {
        balance(
            self.mesh,
            &self.point,
            &self.normal_at(p),
            self.faces,
            self.band,
            &self.options,
        )
    }

    /// The two in-plane components of the normalized residual at the given parameters.
    fn residual_at(&self, p: &Vector2<f64>) -> Option<Vector2<f64>> {
        let balance = self.evaluate(p).ok()??;
        let rot = self.rotation_at(p);
        let r = balance.residual();
        Some(Vector2::new(
            r.dot(&(rot * self.axis_x.into_inner())),
            r.dot(&(rot * self.axis_y.into_inner())),
        ))
    }
}

impl LeastSquaresProblem<f64, U2, U2> for TransverseFit<'_> {
    type ResidualStorage = Owned<f64, U2>;
    type JacobianStorage = Owned<f64, U2, U2>;
    type ParameterStorage = Owned<f64, U2>;

    fn set_params(&mut self, x: &Vector<f64, U2, Self::ParameterStorage>) {
        self.params = *x;
        self.balance = self.evaluate(x).ok().flatten();
    }

    fn params(&self) -> Vector<f64, U2, Self::ParameterStorage> {
        self.params
    }

    fn residuals(&self) -> Option<Vector<f64, U2, Self::ResidualStorage>> {
        let balance = self.balance.as_ref()?;
        let (e1, e2) = self.in_plane_axes();
        let r = balance.residual();
        Some(Vector2::new(r.dot(&e1), r.dot(&e2)))
    }

    fn jacobian(&self) -> Option<Matrix<f64, U2, U2, Self::JacobianStorage>> {
        let mut jac = Matrix2::zeros();
        for k in 0..2 {
            let mut plus = self.params;
            plus[k] += DELTA;
            let mut minus = self.params;
            minus[k] -= DELTA;
            let d = (self.residual_at(&plus)? - self.residual_at(&minus)?) / (2.0 * DELTA);
            jac.set_column(k, &d);
        }
        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::arc_segments_for_tol;
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom2::{BoundaryData2, BoundaryEditor};
    use crate::geom3::{ExtrudedBoundary3, RevolvedBoundary3};
    use crate::{Iso3, Vector3};
    use std::f64::consts::TAU;

    /// A test case: a mesh, a point on it, the expected transverse normal, and a face mask.
    struct Case {
        mesh: Mesh3,
        point: Point3,
        expected: UnitVec3,
        mask: Option<IndexMask>,
    }

    impl Case {
        fn new(mesh: Mesh3, point: Point3, expected: Vector3) -> Self {
            Self {
                mesh,
                point,
                expected: UnitVec3::new_normalize(expected),
                mask: None,
            }
        }

        /// Keep only the faces whose center satisfies the predicate.
        fn mask_centers(mut self, keep: impl Fn(&Point3) -> bool) -> Self {
            let centers = self.mesh.compute_face_centers().unwrap();
            let mut mask = self.mesh.face_mask(false);
            for (i, c) in centers.iter().enumerate() {
                mask.set(i, keep(c));
            }
            self.mask = Some(mask);
            self
        }

        /// Keep only the faces whose normal satisfies the predicate.
        fn mask_normals(mut self, keep: impl Fn(&UnitVec3) -> bool) -> Self {
            let normals = self.mesh.compute_face_normals().unwrap();
            let mut mask = self.mesh.face_mask(false);
            for (i, n) in normals.iter().enumerate() {
                mask.set(i, keep(n));
            }
            self.mask = Some(mask);
            self
        }

        /// Move the whole case by a rigid transform.
        fn transformed(self, iso: &Iso3) -> Self {
            let mut mesh = self.mesh;
            mesh.transform_in_place(iso);
            Self {
                mesh,
                point: iso * self.point,
                expected: iso * self.expected,
                mask: self.mask,
            }
        }

        /// The expected normal rotated by `angle` about a perpendicular axis chosen by `rng`.
        fn perturbed_guess(&self, angle: f64, rng: &mut RandomGeometry3) -> UnitVec3 {
            let perp = UnitVec3::new_normalize(self.expected.cross(&rng.unit_vec()));
            UnitQuaternion::from_axis_angle(&perp, angle) * self.expected
        }

        fn solve(
            &self,
            guess: &UnitVec3,
            options: Option<TransverseOptions>,
        ) -> Result<TransversePlane> {
            self.mesh
                .transverse_plane(&self.point, guess, self.mask.as_ref(), options)
        }

        /// Solve from a guess `angle` radians off the answer and return the angular error.
        fn error_from(&self, angle: f64, rng: &mut RandomGeometry3) -> f64 {
            let guess = self.perturbed_guess(angle, rng);
            let result = self.solve(&guess, None).unwrap();
            let actual = result.plane.normal;
            assert!(
                result.coverage > 0.05,
                "coverage {} is too low to trust",
                result.coverage
            );
            actual.dot(&self.expected).abs().min(1.0).acos()
        }
    }

    /// Runs the case as built and again under a random rigid transform. Both runs start from a
    /// guess 25 degrees off the answer and must recover the answer to within `tol` radians.
    fn check(case: Case, tol: f64) {
        let mut rng = RandomGeometry3::from_seed(0x7a3e);
        let e0 = case.error_from(25.0_f64.to_radians(), &mut rng);
        assert!(
            e0 < tol,
            "untransformed: error {e0:.3e} rad exceeds {tol:.1e}"
        );

        let moved = case.transformed(&rng.iso3(10.0));
        let e1 = moved.error_from(25.0_f64.to_radians(), &mut rng);
        assert!(
            e1 < tol,
            "transformed: error {e1:.3e} rad exceeds {tol:.1e}"
        );
    }

    // ---------------------------------------------------------------------------------------
    //  Shapes
    // ---------------------------------------------------------------------------------------

    /// A cylinder along z with its caps masked away, and a point on the wall.
    fn cylinder_wall(point: Point3) -> Case {
        let mesh = Mesh3::create_cylinder(1.0, 4.0, 1e-3).unwrap();
        Case::new(mesh, point, Vector3::z()).mask_normals(|n| n.z.abs() < 0.5)
    }

    /// A cone along z with the given half-angle, base cap masked away, and a point on the axis.
    fn cone(half_angle: f64) -> Case {
        let height = 2.0;
        let radius = height * half_angle.tan();
        let mesh = Mesh3::create_cone(radius, height, 1e-3).unwrap();
        Case::new(mesh, Point3::origin(), Vector3::z()).mask_normals(|n| n.z > -0.5)
    }

    /// A box spanning 4 units along x (from 0 to 4), 2 along y, and 2 along z. Its top face drops
    /// by `k_top` per unit of x, while its bottom face rises by `k_bottom` per unit of x. Scaling
    /// y linearly in x keeps every face planar, so the result is a true wedge with flat sides.
    fn tapered_box(k_top: f64, k_bottom: f64) -> Mesh3 {
        let box_mesh = Mesh3::create_box(4.0, 2.0, 2.0, false);
        let points: Vec<Point3> = box_mesh
            .points()
            .iter()
            .map(|p| {
                let x = p.x + 2.0;
                let top = 1.0 - k_top * x;
                let bottom = -1.0 + k_bottom * x;
                let y = bottom + 0.5 * (p.y + 1.0) * (top - bottom);
                Point3::new(x, y, p.z)
            })
            .collect();
        Mesh3::new(points, box_mesh.faces().to_vec(), false)
    }

    fn extruded(profile: &[(f64, f64)], closed: bool, length: f64) -> Mesh3 {
        let mut data = if closed {
            BoundaryData2::new_closed()
        } else {
            BoundaryData2::new_open_xy(profile[0].0, profile[0].1)
        };
        let mut cursor = data.get_cursor(None);
        for &(x, y) in profile.iter().skip(if closed { 0 } else { 1 }) {
            cursor.add_seg_xy(x, y);
        }
        let b = data.try_to_boundary().unwrap();
        ExtrudedBoundary3::new(b, Iso3::identity(), length)
            .to_mesh(1e-3)
            .unwrap()
    }

    fn revolved(profile: &[(f64, f64)], closed: bool, tol: f64) -> Mesh3 {
        let mut data = if closed {
            BoundaryData2::new_closed()
        } else {
            BoundaryData2::new_open_xy(profile[0].0, profile[0].1)
        };
        let mut cursor = data.get_cursor(None);
        for &(x, y) in profile.iter().skip(if closed { 0 } else { 1 }) {
            cursor.add_seg_xy(x, y);
        }
        let b = data.try_to_boundary().unwrap();
        RevolvedBoundary3::new(b, Iso3::identity(), TAU)
            .to_mesh(tol)
            .unwrap()
    }

    /// The mid-angle of the first ring of a full revolution tessellated at `tol` with largest
    /// radius `r_max`. The faces of a ring are all parallel to the chord across it, so the plane
    /// normal to the chord through a point at the ring's mid-angle is the exact across plane.
    fn ring_mid_angle(r_max: f64, tol: f64) -> f64 {
        let n = arc_segments_for_tol(r_max, TAU, tol).unwrap();
        0.5 * TAU / n as f64
    }

    /// `Iso3::from_ry(t)` carries a profile point `(x, y, 0)` to `(x cos t, y, -x sin t)`.
    fn revolved_point(x: f64, y: f64, t: f64) -> Point3 {
        Point3::new(x * t.cos(), y, -x * t.sin())
    }

    /// The direction of travel of a revolved point at angle `t`.
    fn revolved_tangent(t: f64) -> Vector3 {
        Vector3::new(-t.sin(), 0.0, -t.cos())
    }

    // ---------------------------------------------------------------------------------------
    //  Tests
    // ---------------------------------------------------------------------------------------

    #[test]
    fn cylinder_from_a_point_on_the_wall() {
        check(cylinder_wall(Point3::new(1.0, 0.0, 0.3)), 1e-7);
    }

    #[test]
    fn cylinder_from_a_point_on_the_axis() {
        check(cylinder_wall(Point3::new(0.0, 0.0, -0.4)), 1e-7);
    }

    #[test]
    fn cylinder_without_a_mask_still_finds_the_axis() {
        // The caps are parallel to the across plane, so a plane at any tilt does not intersect
        // them. Leaving them in therefore does not affect the answer.
        let mesh = Mesh3::create_cylinder(1.0, 4.0, 1e-3).unwrap();
        check(
            Case::new(mesh, Point3::new(1.0, 0.0, 0.0), Vector3::z()),
            1e-7,
        );
    }

    #[test]
    fn cones_of_every_half_angle() {
        // At 35.26 degrees, the normal tensor of a cone is isotropic, so a plain eigenvector
        // approach has nothing to select. Past 45 degrees, a naive fixed-point iteration on the
        // translation diverges. Neither poses a difficulty for the root solve.
        for deg in [10.0, 30.0, 35.264, 45.0, 60.0] {
            let case = cone((deg as f64).to_radians());
            check(case, 1e-7);
        }
    }

    #[test]
    fn cone_from_a_point_on_the_surface() {
        let half_angle = 20.0_f64.to_radians();
        let height = 2.0;
        let radius = height * half_angle.tan();
        let mesh = Mesh3::create_cone(radius, height, 1e-3).unwrap();
        // Radius at z = 0 is half the base radius.
        let case = Case::new(mesh, Point3::new(0.5 * radius, 0.0, 0.0), Vector3::z())
            .mask_normals(|n| n.z > -0.5);
        check(case, 1e-7);
    }

    #[test]
    fn symmetric_wedge_balances_its_two_faces() {
        // Top and bottom converge along +x by the same amount; the sweep direction is x. The end
        // faces are masked away so they cannot be cut by a badly tilted trial plane.
        let mesh = tapered_box(0.125, 0.125);
        let case = Case::new(mesh, Point3::new(2.0, 0.75, 0.3), Vector3::x())
            .mask_normals(|n| n.x.abs() < 0.9);
        check(case, 1e-7);
    }

    #[test]
    fn one_sided_wedge_finds_the_bisector() {
        // Flat bottom at y = -1, top sloping from y = 1 down to y = 0.5 over 4 units of x. The
        // across plane is normal to the bisector of the two faces, the medial-axis direction.
        let mesh = tapered_box(0.125, 0.0);
        let top = Vector3::new(4.0, -0.5, 0.0).normalize();
        let bisector = top + Vector3::x();
        let case = Case::new(mesh, Point3::new(2.0, -1.0, 0.3), bisector)
            .mask_normals(|n| n.x.abs() < 0.9);
        check(case, 1e-7);
    }

    #[test]
    fn torus_section_is_normal_to_the_spine() {
        let tol = 1e-3;
        let r_tube = 0.5;
        let r_spine = 2.0;
        // A polygonal tube profile centered at radius r_spine on the axis's equator.
        let profile: Vec<(f64, f64)> = (0..32)
            .map(|i| {
                let a = TAU * i as f64 / 32.0;
                (r_spine + r_tube * a.cos(), r_tube * a.sin())
            })
            .collect();
        let mesh = revolved(&profile, true, tol);
        let t = ring_mid_angle(r_spine + r_tube, tol);
        let point = revolved_point(r_spine + r_tube, 0.0, t);
        let case = Case::new(mesh, point, revolved_tangent(t));
        check(case, 1e-7);
    }

    #[test]
    fn straight_groove_with_and_without_a_mask() {
        // A slot in a flat plate, extruded along z.
        let profile = [
            (-1.5, 0.0),
            (-0.5, 0.0),
            (-0.5, -0.5),
            (0.5, -0.5),
            (0.5, 0.0),
            (1.5, 0.0),
        ];
        let mesh = extruded(&profile, false, 2.0);
        let point = Point3::new(0.1, -0.5, 1.0);
        check(Case::new(mesh.clone(), point, Vector3::z()), 1e-7);

        let case = Case::new(mesh, point, Vector3::z()).mask_centers(|c| c.x.abs() < 0.6);
        assert!(case.mask.as_ref().unwrap().count_true() > 0);
        check(case, 1e-7);
    }

    #[test]
    fn circumferential_groove_on_a_shaft() {
        let tol = 1e-3;
        // Radial x, axial y: a shaft of radius 2 with a notch down to radius 1.5 across |y| < 0.25.
        let profile = [
            (2.0, -1.0),
            (2.0, -0.25),
            (1.5, -0.25),
            (1.5, 0.25),
            (2.0, 0.25),
            (2.0, 1.0),
        ];
        let mesh = revolved(&profile, false, tol);
        let t = ring_mid_angle(2.0, tol);
        let point = revolved_point(1.5, 0.0, t);
        let case = Case::new(mesh, point, revolved_tangent(t))
            .mask_centers(|c| (c.x * c.x + c.z * c.z).sqrt() < 1.99 || c.y.abs() < 0.26);
        assert!(case.mask.as_ref().unwrap().count_true() > 0);
        check(case, 1e-7);
    }

    /// A section broken into pieces by holes in the mesh: only the piece nearest the point takes
    /// part, and it still determines the axis.
    #[test]
    fn a_fragmented_section_uses_the_run_nearest_the_point() {
        let case = cylinder_wall(Point3::new(1.0, 0.0, 0.0)).mask_centers(|c| c.x.abs() > 0.3);
        let mut rng = RandomGeometry3::from_seed(5);
        let guess = case.perturbed_guess(25.0_f64.to_radians(), &mut rng);
        let result = case.solve(&guess, None).unwrap();
        let err = result
            .plane
            .normal
            .dot(&case.expected)
            .abs()
            .min(1.0)
            .acos();
        assert!(err < 1e-7, "error {err:.3e} rad");
        // One arc of roughly 145 degrees, so less than half coverage.
        assert!(
            result.coverage > 0.2 && result.coverage < 0.5,
            "{}",
            result.coverage
        );
    }

    /// Two parallel cylinders side by side: the plane cuts both, and only the loop around the
    /// point's own cylinder takes part.
    #[test]
    fn a_neighbouring_feature_is_left_out() {
        let single = cylinder_wall(Point3::new(1.0, 0.0, 0.0));
        let alone = single.solve(&single.expected, None).unwrap();

        let mut mesh = Mesh3::create_cylinder(1.0, 4.0, 1e-3).unwrap();
        let mut other = Mesh3::create_cylinder(1.0, 4.0, 1e-3).unwrap();
        other.transform_in_place(&Iso3::translation(3.0, 0.0, 0.0));
        mesh.append_in_place(&other).unwrap();
        let case = Case::new(mesh, Point3::new(1.0, 0.0, 0.0), Vector3::z())
            .mask_normals(|n| n.z.abs() < 0.5);
        let mut rng = RandomGeometry3::from_seed(21);
        let guess = case.perturbed_guess(25.0_f64.to_radians(), &mut rng);
        let result = case.solve(&guess, None).unwrap();
        let err = result
            .plane
            .normal
            .dot(&case.expected)
            .abs()
            .min(1.0)
            .acos();
        assert!(err < 1e-7, "error {err:.3e} rad");
        assert_eq!(result.face_count, alone.face_count);
    }

    #[test]
    fn half_cylinder_still_determines_the_axis() {
        let case = cylinder_wall(Point3::new(1.0, 0.0, 0.0)).mask_centers(|c| c.x > 0.05);
        let mut rng = RandomGeometry3::from_seed(3);
        let guess = case.perturbed_guess(25.0_f64.to_radians(), &mut rng);
        let result = case.solve(&guess, None).unwrap();
        let err = result
            .plane
            .normal
            .dot(&case.expected)
            .abs()
            .min(1.0)
            .acos();
        assert!(err < 1e-7, "error {err:.3e} rad");
        // Half a loop has less coverage than a full one but is still well determined.
        assert!(
            result.coverage > 0.2 && result.coverage < 0.5,
            "{}",
            result.coverage
        );
    }

    #[test]
    fn a_full_loop_has_half_coverage() {
        let case = cylinder_wall(Point3::new(0.0, 0.0, 0.0));
        let result = case.solve(&case.expected, None).unwrap();
        assert!((result.coverage - 0.5).abs() < 1e-3, "{}", result.coverage);
        assert!(
            (result.sensitivity - 0.5).abs() < 1e-3,
            "{}",
            result.sensitivity
        );
        assert!(result.residual < 1e-12);
    }

    #[test]
    fn a_single_flat_face_is_underdetermined() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let case =
            Case::new(mesh, Point3::new(0.0, 0.0, 1.0), Vector3::x()).mask_normals(|n| n.z > 0.5);
        let err = case
            .solve(&Vector3::x_axis(), None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("undetermined"), "{err}");
    }

    #[test]
    fn a_plane_missing_the_faces_is_an_error() {
        let case = cylinder_wall(Point3::new(0.0, 0.0, 10.0));
        let err = case
            .solve(&Vector3::z_axis(), None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("does not cut"), "{err}");
    }

    #[test]
    fn invalid_options_are_rejected() {
        let case = cylinder_wall(Point3::origin());
        let bad = TransverseOptions {
            taper_tilt: FRAC_PI_2,
            ..Default::default()
        };
        assert!(case.solve(&Vector3::z_axis(), Some(bad)).is_err());
    }

    /// The longitudinal plane of a cylinder is also balanced, but its section is two parallel
    /// lines whose normals do not span the plane, so it is refused rather than returned.
    #[test]
    fn longitudinal_plane_of_a_cylinder_is_refused() {
        let case = cylinder_wall(Point3::origin());
        let err = case
            .solve(&Vector3::x_axis(), None)
            .unwrap_err()
            .to_string();
        assert!(err.contains("undetermined"), "{err}");
    }

    /// Characterizes how far a guess can be from the axis and still converge to it.
    #[test]
    fn basin_of_the_across_plane_reaches_at_least_45_degrees() {
        let mut rng = RandomGeometry3::from_seed(99);
        for deg in [35.0, 45.0] {
            let case = cylinder_wall(Point3::new(1.0, 0.0, 0.0));
            let err = case.error_from((deg as f64).to_radians(), &mut rng);
            assert!(err < 1e-6, "cylinder, guess {deg} deg off: error {err:.3e}");

            let case = cone(20.0_f64.to_radians());
            let err = case.error_from((deg as f64).to_radians(), &mut rng);
            assert!(err < 1e-6, "cone, guess {deg} deg off: error {err:.3e}");
        }
    }
}
