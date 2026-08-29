//! A simultaneous alignment of several meshes to each other in one combined Levenberg-Marquardt
//! minimization.
//!
//! This is a bundle adjustment rather than a pose graph optimization: it carries a single
//! transformation for each mesh except one, which is held fixed, and solves for all of them at
//! once against the raw correspondences. It was built to register metrology-quality scans of
//! objects with unambiguous morphology, and it works best on low-noise meshes which have already
//! been pre-aligned close to each other with substantial overlap.
//!
//! # How a correspondence constrains two bodies
//!
//! Every alignment point belongs to one mesh and is matched against another, and *both* of those
//! meshes are moving. A single residual therefore contributes two blocks to its jacobian row: the
//! forward derivative with respect to the test mesh's parameters, and the reverse derivative with
//! respect to the reference mesh's. See [`point_surf_jacobian`] and [`point_surf_jacobian_rev`].
//!
//! # Coordinate frames
//!
//! Everything is measured in world coordinates. A test point is moved by its own mesh's transform,
//! the query is pushed into the reference mesh's frame only long enough to find the closest
//! surface point, and that match is brought back out to the world before the residual is taken.
//!
//! This matters because each mesh's [`AlignValues3`] describes its transform *to the world*, so
//! the residual and both jacobian blocks have to live there too for the derivatives to be
//! consistent with it. (The previous generation of this module measured in the reference mesh's
//! local frame while passing world-frame parameters to the jacobian, which did not agree.)
//!
//! # The distance gate is not optional
//!
//! [`MultiOptions3::max_distance`] must be supplied, and there is no `Default` for the options
//! precisely so that it cannot be forgotten.
//!
//! Partially overlapping meshes produce correspondences into regions the other mesh never saw, and
//! the alignment opens with an unweighted least squares solve which has no defense against them.
//! On real scan data a handful of such correspondences moved a seventeen mesh bundle by most of a
//! millimetre away from an already correct answer, and the robust refinement that follows made it
//! worse rather than better: the noise scale is estimated from that spoiled solve, so the displaced
//! mesh looks like an outlier and gets weighted out of its own correction. The measurements are on
//! [`MultiOptions3::max_distance`].

use crate::Result;
use crate::common::align::{RefinementHalt, SolveQuality, TerminationReason};
use crate::common::consensus::weights::{MagsacWeight, estimate_sigma_max};
use crate::common::points::dist;
use crate::geom3::align3::jacobian::{point_surf_jacobian, point_surf_jacobian_rev};
use crate::geom3::align3::mesh::{AlignMesh, generate_alignment_points, interpolated_stdev};
use crate::geom3::align3::{AlignValues3, Dof6, GAPParams, MultiParams3};
use crate::geom3::mesh::MeshSurfPoint;
use crate::geom3::{Align3, MultiOutcome3, Point3, SurfacePoint3};
use crate::na::{DMatrix, DVector, Dyn, Matrix, Owned, U1, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use rayon::prelude::*;

/// The residual degrees of freedom for the MAGSAC++ weight function. The residual is a full
/// Euclidean point-to-surface distance in space, so it follows a chi distribution with three
/// degrees of freedom.
const RESIDUAL_DOF: usize = 3;

// ================================================================================================
// Options
// ================================================================================================

/// Options controlling a simultaneous multi-mesh alignment.
///
/// Construct with [`MultiOptions3::new`], which requires the correspondence distance gate,
/// and override the remaining fields you care about. The robust-estimation fields mirror
/// `AlignOptions3` and behave identically; the two gates at the top are specific to multi-mesh
/// work.
///
/// There is deliberately no `Default`. See [`MultiOptions3::max_distance`] for why the gate
/// cannot be left to a default.
#[derive(Clone, Copy, Debug)]
pub struct MultiOptions3 {
    /// A hard cutoff on correspondence distance, in the units of the geometry. A point whose
    /// closest match on the reference mesh is further than this contributes nothing.
    ///
    /// This is required rather than optional, and that is a deliberate reversal of the earlier
    /// advice to leave it off and let robust weighting do the work. On real multi-scan data that
    /// advice is wrong, and expensively so.
    ///
    /// Meshes overlap only partially, so a point in a region the reference never saw has no
    /// meaningful match at any distance. The alignment opens with an *unweighted* least squares
    /// solve, in order to have residuals from which to estimate a noise scale, and unweighted least
    /// squares has no defense whatsoever against a gross outlier. A handful of correspondences
    /// sampled across non-overlapping regions is enough to drag the whole bundle.
    ///
    /// Measured on seventeen partially overlapping laser scans of a part 127 mm across, starting
    /// from an already converged alignment: with no gate the initial solve moved the group by
    /// 0.83 mm and left ten times as many points beyond a tenth of a millimetre of the reference.
    /// Any gate at all fixed it, taking the drift to 0.09 mm. A gate of 1.0 mm and one of 0.5 mm
    /// gave identical results, which says the offending correspondences were more than a
    /// millimetre out rather than marginal. Robust refinement afterwards did not rescue it, because
    /// the noise scale is estimated by a median absolute deviation from that same bad solve: the
    /// estimate describes the well-fitted core, and MAGSAC then treats the displaced mesh's
    /// correspondences as outliers and stops pulling on them, entrenching the error it inherited.
    ///
    /// Choose it from the geometry rather than from the expected residual. It only has to be
    /// tighter than the scale of a spurious match, so err generous: a value far larger than the
    /// alignment error still removes the correspondences that matter.
    ///
    /// Note that, like `ignore_off_target` in the single-body case, this is a hard gate
    /// re-evaluated as the points move, so it makes the objective piecewise.
    pub max_distance: f64,

    /// An optional cutoff on the angle, in radians, between a test point's normal and its match's
    /// normal. A correspondence exceeding it contributes nothing.
    ///
    /// This suppresses matches onto the far side of a thin wall, where the geometry is close but
    /// facing the wrong way. `None` accepts any orientation.
    pub max_normal_angle: Option<f64>,

    /// The number of iteratively reweighted refinement rounds to perform after the initial
    /// unweighted solve. Zero disables robust weighting entirely.
    pub refinement_steps: usize,

    /// The MAGSAC++ upper noise bound, estimated from the initial residuals via the median
    /// absolute deviation when `None`. See `AlignOptions3::sigma_max` for the units.
    pub sigma_max: Option<f64>,

    /// The Levenberg-Marquardt evaluation budget, as a multiplier on the parameter count.
    ///
    /// A multi-body solve carries `6 * (n - 1)` parameters, so this budget stretches a long way
    /// with many meshes. A solve which exhausts it is not a failure: the alignments are kept and
    /// the outcome reports [`crate::common::SolveQuality::Unconverged`]. That is a common outcome
    /// here, since correspondences flip between every overlapping pair.
    ///
    /// Must be greater than zero.
    pub patience: usize,

    /// An optional degree-of-freedom constraint applied to every non-static mesh.
    pub dof: Option<Dof6>,
}

impl MultiOptions3 {
    /// Options with the required correspondence distance gate and defaults for everything else.
    ///
    /// `max_distance` has no safe default and so has to be supplied; see its documentation for the
    /// measurements behind that. Everything else follows the single-body defaults: four rounds of
    /// robust refinement with the noise scale estimated from the data, no normal gate, no degree of
    /// freedom locks.
    pub fn new(max_distance: f64) -> Self {
        Self {
            max_distance,
            max_normal_angle: None,
            refinement_steps: 4,
            sigma_max: None,
            patience: 100,
            dof: None,
        }
    }
}

// ================================================================================================
// Alignment points
// ================================================================================================

/// A single correspondence in a multi-mesh adjustment: a sample point on one mesh which is being
/// matched against the surface of another.
#[derive(Clone, Debug)]
pub struct MeshAlignPoint {
    /// The index of the mesh this point was sampled from.
    pub mesh_i: usize,

    /// The sample point, in the coordinates of its own mesh.
    pub mp: MeshSurfPoint,

    /// The index of the mesh this point is being matched against.
    pub ref_i: usize,

    /// A base weight for the correspondence, applied on top of any robust weighting.
    pub weight: f64,
}

impl MeshAlignPoint {
    /// Creates a correspondence from a sample `mp` on mesh `mesh_i`, matched against reference
    /// mesh `ref_i` with the given weight.
    pub fn new(mesh_i: usize, mp: MeshSurfPoint, ref_i: usize, weight: f64) -> Self {
        Self {
            mesh_i,
            mp,
            ref_i,
            weight,
        }
    }
}

// ================================================================================================
// Entry points
// ================================================================================================

/// Runs a simultaneous alignment over an explicit set of correspondences.
///
/// Use this when the correspondences have already been chosen, whether by
/// [`multi_mesh_adjustment`]'s own sampling or by pruning them with
/// [`crate::geom3::align3::AlignInfo3`].
///
/// # Failure
///
/// As with the single-body solver, `Err` is reserved for the case where there is no answer at
/// all: rejected arguments, or an initial solve that broke down. A solve which merely exhausts
/// its budget is reported on the outcome and its alignments kept.
pub fn multi_mesh_adjustment_with_points(
    meshes: &[AlignMesh],
    points: Vec<MeshAlignPoint>,
    static_i: usize,
    opts: &MultiOptions3,
) -> Result<MultiOutcome3> {
    validate(meshes, &points, static_i, opts)?;
    solve_bundle(MultiMeshProblem::new(meshes, points, static_i, opts)?, opts)
}

/// Samples correspondences between every pair of meshes and runs a simultaneous alignment.
///
/// The static mesh is chosen automatically as the one the others reference most broadly; see
/// `correspondence_matrix` for how that is scored.
pub fn multi_mesh_adjustment(
    meshes: &[AlignMesh],
    opts: &MultiOptions3,
    sample_opts: &GAPParams,
) -> Result<MultiOutcome3> {
    let reference_order = reference_priority(meshes, sample_opts);
    let static_i = reference_order[0];

    // Build the work list so that each unordered pair of meshes produces correspondences in one
    // direction only, never both.
    let mut work_list = Vec::new();
    let mut meshes_to_test = (0..meshes.len()).collect::<Vec<_>>();
    for ref_i in reference_order {
        meshes_to_test.retain(|j| *j != ref_i);
        for &mesh_i in meshes_to_test.iter() {
            work_list.push((mesh_i, ref_i));
        }
    }

    let points = work_list
        .par_iter()
        .map(|&(mesh_i, ref_i)| {
            let t = meshes[ref_i]
                .transform()
                .inv_mul(&meshes[mesh_i].transform());
            let samples =
                generate_alignment_points(meshes[mesh_i].mesh, meshes[ref_i].mesh, &t, sample_opts);

            samples
                .into_iter()
                .map(|mp| {
                    let weight = meshes[mesh_i].weights.map_or(1.0, |providers| {
                        providers.iter().map(|p| p.weight(&mp)).product()
                    });
                    MeshAlignPoint::new(mesh_i, mp, ref_i, weight)
                })
                .collect::<Vec<_>>()
        })
        .flatten()
        .collect::<Vec<_>>();

    multi_mesh_adjustment_with_points(meshes, points, static_i, opts)
}

// ================================================================================================
// The solve
// ================================================================================================

fn validate(
    meshes: &[AlignMesh],
    points: &[MeshAlignPoint],
    static_i: usize,
    opts: &MultiOptions3,
) -> Result<()> {
    if opts.patience == 0 {
        return Err("patience must be greater than zero".into());
    }
    if let Some(s) = opts.sigma_max
        && (!s.is_finite() || s <= 0.0)
    {
        return Err(format!("sigma_max is {s}, but must be finite and strictly positive").into());
    }
    if !opts.max_distance.is_finite() || opts.max_distance <= 0.0 {
        return Err(format!(
            "max_distance is {}, but must be finite and strictly positive",
            opts.max_distance
        )
        .into());
    }
    if static_i >= meshes.len() {
        return Err(format!(
            "the static mesh index {} is out of range for {} meshes",
            static_i,
            meshes.len()
        )
        .into());
    }
    if let Some(p) = points
        .iter()
        .find(|p| p.mesh_i >= meshes.len() || p.ref_i >= meshes.len())
    {
        return Err(format!(
            "an alignment point references mesh {} against mesh {}, but there are only {} meshes",
            p.mesh_i,
            p.ref_i,
            meshes.len()
        )
        .into());
    }

    Ok(())
}

fn solve_bundle(problem: MultiMeshProblem<'_>, opts: &MultiOptions3) -> Result<MultiOutcome3> {
    let lm = LevenbergMarquardt::<f64>::new().with_patience(opts.patience);
    let n_params = problem.params.param_count();

    let (mut problem, termination) = run(&lm, problem);
    if !SolveQuality::from_termination(&termination).is_usable() {
        return Err(format!("Failed to align meshes to each other: {termination:?}").into());
    }

    let mut solves = vec![termination];
    let mut halt = None;

    if opts.refinement_steps > 0 {
        match resolve_sigma_max(opts, &problem) {
            None => halt = Some(RefinementHalt::NoNoiseEstimate),
            Some(sigma_max) => {
                let weighting = MagsacWeight::new(sigma_max, RESIDUAL_DOF);

                for _ in 0..opts.refinement_steps {
                    problem.refresh_inv_sigma();

                    let weighted = problem.count_if_reweighted(&weighting);
                    if weighted < n_params {
                        halt = Some(RefinementHalt::Underdetermined {
                            weighted,
                            params: n_params,
                        });
                        break;
                    }

                    let last_good = problem.params.storage().clone();
                    problem.apply_magsac_weights(&weighting);

                    let (next, termination) = run(&lm, problem);
                    problem = next;

                    if !SolveQuality::from_termination(&termination).is_usable() {
                        problem.restore(&last_good);
                        halt = Some(RefinementHalt::SolveFailed(termination));
                        break;
                    }
                    solves.push(termination);
                }
            }
        }
    }

    Ok(MultiOutcome3::new(problem.alignments(), solves, halt))
}

fn run<'a>(
    lm: &LevenbergMarquardt<f64>,
    problem: MultiMeshProblem<'a>,
) -> (MultiMeshProblem<'a>, TerminationReason) {
    let (result, report) = lm.minimize(problem);
    (result, report.termination)
}

fn resolve_sigma_max(opts: &MultiOptions3, problem: &MultiMeshProblem<'_>) -> Option<f64> {
    match opts.sigma_max {
        Some(s) => Some(s),
        None => estimate_sigma_max(&problem.normalized_residuals()),
    }
}

// ================================================================================================
// The problem
// ================================================================================================

/// The per-correspondence state produced by one pass of `move_points`.
struct Moved {
    p_world: Point3,
    c_world: SurfacePoint3,
    residual: f64,
    target_weight: f64,
    target_sigma: f64,
}

struct MultiMeshProblem<'a> {
    meshes: &'a [AlignMesh<'a>],
    handles: Vec<MeshAlignPoint>,
    params: MultiParams3,

    /// The alignment values of every body for the current parameters, cached so the residual and
    /// jacobian loops don't recompute them per correspondence.
    values: Vec<AlignValues3>,

    /// The test point of each correspondence in world coordinates.
    moved: Vec<Point3>,

    /// The matching surface point of each correspondence, also in world coordinates.
    closest: Vec<SurfacePoint3>,

    /// The signed geometric distance of each correspondence, in model units.
    residuals: Vec<f64>,

    /// Each test point's own measurement standard deviation, interpolated from its mesh once at
    /// construction because the point never moves within its own mesh.
    test_sigma: Vec<f64>,

    /// The reference mesh's standard deviation at each current match.
    target_sigma: Vec<f64>,

    /// The reciprocal of the combined test-and-target standard deviation, or 1.0 where there is no
    /// uncertainty at all. Held fixed across a solve and refreshed between them, for the same
    /// reason as in the single-body case.
    inv_sigma: Vec<f64>,

    /// The weight contributed by the correspondence itself: the point's base weight and the two
    /// geometric gates.
    target_weights: Vec<f64>,

    /// The weight contributed by MAGSAC++ reweighting, held fixed across a solve.
    magsac_weights: Vec<f64>,

    max_distance: f64,
    min_normal_dot: Option<f64>,
}

impl<'a> MultiMeshProblem<'a> {
    fn new(
        meshes: &'a [AlignMesh],
        handles: Vec<MeshAlignPoint>,
        static_i: usize,
        opts: &MultiOptions3,
    ) -> Result<Self> {
        let centers = meshes
            .iter()
            .map(|m| m.mesh.aabb().center())
            .collect::<Vec<_>>();
        let initial = meshes.iter().map(|m| m.transform()).collect::<Vec<_>>();
        let params = MultiParams3::from_centers(static_i, &centers, Some(&initial), opts.dof)?;

        let n = handles.len();

        // A test point is fixed within its own mesh, so its uncertainty never changes and is
        // resolved once here rather than on every step.
        let test_sigma = handles
            .iter()
            .map(|h| sigma_at(&meshes[h.mesh_i], &h.mp))
            .collect();

        let values = params.compute_all_values();
        let mut item = Self {
            meshes,
            handles,
            params,
            values,
            moved: vec![Point3::origin(); n],
            closest: vec![SurfacePoint3::new(Point3::origin(), crate::Vector3::z_axis()); n],
            residuals: vec![0.0; n],
            test_sigma,
            target_sigma: vec![0.0; n],
            inv_sigma: vec![1.0; n],
            target_weights: vec![1.0; n],
            magsac_weights: vec![1.0; n],
            max_distance: opts.max_distance,
            min_normal_dot: opts.max_normal_angle.map(f64::cos),
        };

        item.move_points();
        item.refresh_inv_sigma();
        Ok(item)
    }

    /// Moves every test point into world coordinates, finds its match on the reference mesh,
    /// brings that match back out to the world, and recomputes the residual and the geometric
    /// weights.
    fn move_points(&mut self) {
        self.values = self.params.compute_all_values();

        let values = &self.values;
        let meshes = self.meshes;
        let max_distance = self.max_distance;
        let min_normal_dot = self.min_normal_dot;

        let collected: Vec<Moved> = self
            .handles
            .par_iter()
            .map(|h| {
                let t_test = values[h.mesh_i].transform;
                let t_ref = values[h.ref_i].transform;

                // Into the world, then into the reference mesh's frame just long enough to find
                // the closest surface point.
                let p_world = t_test * h.mp.point();
                let query = t_ref.inverse_transform_point(&p_world);
                let mp = meshes[h.ref_i].mesh.surface_closest_to(&query);

                // ...and back out to the world, where the residual is measured.
                let c_world = SurfacePoint3::new(t_ref * mp.point(), t_ref.rotation * mp.normal());

                let d = dist(&p_world, &c_world.point);
                let residual = d * c_world.scalar_projection(&p_world).signum();

                let mut weight = h.weight;
                if d > max_distance {
                    weight = 0.0;
                }
                if let Some(min_dot) = min_normal_dot {
                    let n_test = t_test.rotation * h.mp.normal();
                    if n_test.dot(&c_world.normal) < min_dot {
                        weight = 0.0;
                    }
                }

                Moved {
                    p_world,
                    c_world,
                    residual,
                    target_weight: weight,
                    target_sigma: sigma_at(&meshes[h.ref_i], &mp),
                }
            })
            .collect();

        for (i, m) in collected.into_iter().enumerate() {
            self.moved[i] = m.p_world;
            self.closest[i] = m.c_world;
            self.residuals[i] = m.residual;
            self.target_weights[i] = m.target_weight;
            self.target_sigma[i] = m.target_sigma;
        }
    }

    /// Puts the parameters back to a previous state and rebuilds the correspondences to match.
    fn restore(&mut self, storage: &DVector<f64>) {
        self.params.set_storage(storage);
        self.move_points();
    }

    /// Recombines the test and target standard deviations in quadrature.
    fn refresh_inv_sigma(&mut self) {
        for i in 0..self.handles.len() {
            let t = self.test_sigma[i];
            let r = self.target_sigma[i];
            let combined = (t * t + r * r).sqrt();
            self.inv_sigma[i] = if combined > 0.0 { 1.0 / combined } else { 1.0 };
        }
    }

    fn normalized_residuals(&self) -> Vec<f64> {
        self.residuals
            .iter()
            .zip(self.inv_sigma.iter())
            .map(|(r, inv)| r * inv)
            .collect()
    }

    /// The factor applied to both the residual and the jacobian row of a correspondence.
    fn scale(&self, i: usize) -> f64 {
        (self.target_weights[i] * self.magsac_weights[i]).sqrt() * self.inv_sigma[i]
    }

    fn apply_magsac_weights(&mut self, weighting: &MagsacWeight) {
        for i in 0..self.handles.len() {
            let r = (self.residuals[i] * self.inv_sigma[i]).abs();
            self.magsac_weights[i] = weighting.weight(r);
        }
    }

    fn count_if_reweighted(&self, weighting: &MagsacWeight) -> usize {
        (0..self.handles.len())
            .filter(|&i| {
                let r = (self.residuals[i] * self.inv_sigma[i]).abs();
                self.target_weights[i] * weighting.weight(r) > 0.0
            })
            .count()
    }

    /// Builds one [`Align3`] per mesh, with that mesh's own correspondence residuals.
    fn alignments(&self) -> Vec<Align3> {
        let mut grouped = vec![Vec::new(); self.meshes.len()];
        for (i, h) in self.handles.iter().enumerate() {
            grouped[h.mesh_i].push(self.residuals[i]);
        }

        self.values
            .iter()
            .enumerate()
            .zip(grouped)
            .map(|((i, v), residuals)| {
                let body = self.params.body(i);
                Align3::new(v.transform, v.align, body.local, body.offset, residuals)
            })
            .collect()
    }
}

/// The measurement uncertainty of a mesh at a surface point, preferring an explicit override on
/// the [`AlignMesh`] over the mesh's own `point_stdev` attribute. Zero when neither is present.
fn sigma_at(mesh: &AlignMesh, mp: &MeshSurfPoint) -> f64 {
    let stdev = mesh.uncertainty.or_else(|| mesh.mesh.point_stdev());
    stdev.map_or(0.0, |s| interpolated_stdev(mesh.mesh, mp, s))
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for MultiMeshProblem<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params.set_storage(x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.storage().clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.handles.len());
        for i in 0..self.handles.len() {
            res[i] = self.residuals[i] * self.scale(i);
        }
        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            self.handles.len(),
            self.params.param_count(),
        );

        for (i, h) in self.handles.iter().enumerate() {
            let scale = self.scale(i);
            let p = &self.moved[i];
            let c = &self.closest[i];

            // The correspondence constrains the test mesh...
            let fwd = point_surf_jacobian(p, c, &self.values[h.mesh_i]) * scale;
            self.params.add_jacobian_block(&mut jac, i, h.mesh_i, &fwd);

            // ...and the reference mesh, which is also free to move.
            let rev = point_surf_jacobian_rev(p, c, &self.values[h.ref_i]) * scale;
            self.params.add_jacobian_block(&mut jac, i, h.ref_i, &rev);
        }

        Some(jac)
    }
}

// ================================================================================================
// Choosing the static mesh
// ================================================================================================

/// Orders the meshes by how broadly the others reference them, most-referenced first. The head of
/// the returned order is the best candidate for the static mesh. The scoring lives in
/// [`crate::common::align::reference`].
fn reference_priority(meshes: &[AlignMesh], params: &GAPParams) -> Vec<usize> {
    crate::common::align::reference::reference_priority(correspondence_matrix(meshes, params))
}

/// The pairwise correspondence counts of the meshes: each `i, j` entry is the number of sample
/// points in mesh `j` which are a good match for mesh `i`.
fn correspondence_matrix(meshes: &[AlignMesh], params: &GAPParams) -> DMatrix<f64> {
    crate::common::align::reference::correspondence_matrix(meshes.len(), |i, j| {
        let t = meshes[i].transform().inv_mul(&meshes[j].transform());
        generate_alignment_points(meshes[j].mesh, meshes[i].mesh, &t, params).len() as f64
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Iso3;
    use crate::{Mesh3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// A correspondence gate wider than the test boxes themselves, so that the synthetic cases
    /// exercise the solver rather than the gate. Cases which are about the gate set their own.
    const WIDE_GATE: f64 = 20.0;

    fn test_opts() -> MultiOptions3 {
        MultiOptions3::new(WIDE_GATE)
    }

    fn box_mesh() -> Mesh3 {
        Mesh3::create_box(10.0, 6.0, 4.0, false)
    }

    /// Correspondences from every mesh onto the one before it, sampled directly from the mesh
    /// surface. This bypasses `generate_alignment_points` so the tests exercise the solver rather
    /// than the candidate filtering.
    fn chain_points(meshes: &[Mesh3], spacing: f64) -> Vec<MeshAlignPoint> {
        let mut points = Vec::new();
        for (i, mesh) in meshes.iter().enumerate().skip(1) {
            for mp in mesh.sample_surface_poisson(spacing, None) {
                points.push(MeshAlignPoint::new(i, mp, i - 1, 1.0));
            }
        }
        points
    }

    fn alignment_meshes<'a>(meshes: &'a [Mesh3], initial: &'a [Iso3]) -> Vec<AlignMesh<'a>> {
        meshes
            .iter()
            .zip(initial.iter())
            .map(|(m, t)| AlignMesh::new(m, None, Some(t), None))
            .collect()
    }

    // ============================================================================================
    // Recovery
    // ============================================================================================

    #[test]
    fn a_disturbed_mesh_is_brought_back_onto_its_twin() -> Result<()> {
        // Two identical boxes, one of them displaced. The adjustment should undo the displacement
        // exactly, since the two surfaces coincide when aligned.
        let meshes = vec![box_mesh(), box_mesh()];
        let disturb = Iso3::new(Vector3::new(0.2, -0.15, 0.1), Vector3::z() * (PI / 90.0));
        let initial = vec![Iso3::identity(), disturb];

        let ams = alignment_meshes(&meshes, &initial);
        let points = chain_points(&meshes, 0.4);

        let outcome = multi_mesh_adjustment_with_points(&ams, points, 0, &test_opts())?;

        assert_eq!(outcome.len(), 2);
        // Mesh 0 is static and must not have moved at all.
        assert_relative_eq!(
            outcome.alignment(0).full_transform().to_matrix(),
            Iso3::identity().to_matrix(),
            epsilon = 1e-12
        );
        // Mesh 1 should have come back to the identity, undoing its initial displacement.
        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            Iso3::identity().to_matrix(),
            epsilon = 1e-6
        );

        Ok(())
    }

    #[test]
    fn three_meshes_settle_together() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh(), box_mesh()];
        let initial = vec![
            Iso3::identity(),
            Iso3::new(Vector3::new(0.15, 0.0, -0.1), Vector3::y() * (PI / 120.0)),
            Iso3::new(Vector3::new(-0.1, 0.12, 0.05), Vector3::x() * (PI / 100.0)),
        ];

        let ams = alignment_meshes(&meshes, &initial);
        let points = chain_points(&meshes, 0.5);

        let outcome = multi_mesh_adjustment_with_points(&ams, points, 0, &test_opts())?;

        for i in 0..3 {
            assert_relative_eq!(
                outcome.alignment(i).full_transform().to_matrix(),
                Iso3::identity().to_matrix(),
                epsilon = 1e-5
            );
        }

        Ok(())
    }

    #[test]
    fn the_static_mesh_never_moves() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![
            Iso3::translation(0.3, 0.0, 0.0),
            Iso3::translation(0.0, 0.25, 0.0),
        ];
        let ams = alignment_meshes(&meshes, &initial);

        // Make mesh 1 the static one this time, so the check is not trivially about index zero.
        let outcome =
            multi_mesh_adjustment_with_points(&ams, chain_points(&meshes, 0.5), 1, &test_opts())?;

        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            initial[1].to_matrix(),
            epsilon = 1e-12
        );

        Ok(())
    }

    #[test]
    fn residuals_are_grouped_by_the_mesh_they_came_from() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(); 3];
        let ams = alignment_meshes(&meshes, &initial);
        let points = chain_points(&meshes, 0.6);

        let from_1 = points.iter().filter(|p| p.mesh_i == 1).count();
        let from_2 = points.iter().filter(|p| p.mesh_i == 2).count();

        let outcome = multi_mesh_adjustment_with_points(&ams, points, 0, &test_opts())?;

        // Mesh 0 sourced no correspondences, so it has no residuals of its own.
        assert_eq!(outcome.alignment(0).residuals().len(), 0);
        assert_eq!(outcome.alignment(1).residuals().len(), from_1);
        assert_eq!(outcome.alignment(2).residuals().len(), from_2);

        Ok(())
    }

    // ============================================================================================
    // Robustness and reporting
    // ============================================================================================

    #[test]
    fn locked_dof_are_honored() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(), Iso3::translation(0.3, 0.0, 0.0)];
        let ams = alignment_meshes(&meshes, &initial);

        let opts = MultiOptions3 {
            dof: Some(Dof6::new(false, true, true, true, true, true)),
            ..test_opts()
        };
        let outcome =
            multi_mesh_adjustment_with_points(&ams, chain_points(&meshes, 0.5), 0, &opts)?;

        // tx is locked, so the x displacement cannot be undone and must survive untouched.
        assert_relative_eq!(
            outcome.alignment(1).full_transform().translation.vector.x,
            0.3,
            epsilon = 1e-9
        );

        Ok(())
    }

    #[test]
    fn exhausted_budget_is_reported_not_raised() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![
            Iso3::identity(),
            Iso3::new(Vector3::new(1.5, -1.0, 0.8), Vector3::z() * (PI / 12.0)),
        ];
        let ams = alignment_meshes(&meshes, &initial);

        let opts = MultiOptions3 {
            patience: 1,
            ..test_opts()
        };
        let outcome =
            multi_mesh_adjustment_with_points(&ams, chain_points(&meshes, 0.5), 0, &opts)?;

        assert_eq!(outcome.quality(), SolveQuality::Unconverged);
        assert!(outcome.solves().contains(&TerminationReason::LostPatience));
        // The alignments are still real geometry rather than something degenerate.
        for a in outcome.alignments() {
            assert!(a.full_transform().to_matrix().iter().all(|v| v.is_finite()));
        }

        Ok(())
    }

    #[test]
    fn gross_outliers_are_rejected() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let disturb = Iso3::new(Vector3::new(0.2, -0.15, 0.1), Vector3::z() * (PI / 120.0));
        let initial = vec![Iso3::identity(), disturb];
        let ams = alignment_meshes(&meshes, &initial);

        // Corrupt a tenth of the correspondences by throwing their sample points off the surface.
        let mut points = chain_points(&meshes, 0.4);
        for (k, p) in points.iter_mut().enumerate() {
            if k.is_multiple_of(10) {
                p.mp.sp.point += Vector3::new(0.0, 0.0, 5.0);
            }
        }

        let naive = multi_mesh_adjustment_with_points(
            &ams,
            points.clone(),
            0,
            &MultiOptions3 {
                refinement_steps: 0,
                ..test_opts()
            },
        )?;
        let robust = multi_mesh_adjustment_with_points(&ams, points, 0, &test_opts())?;

        let error = |o: &MultiOutcome3| {
            (o.alignment(1).full_transform().to_matrix() - Iso3::identity().to_matrix()).norm()
        };

        assert!(
            error(&robust) < error(&naive) * 0.5,
            "robust weighting should have suppressed the outliers: \
             naive {}, robust {}",
            error(&naive),
            error(&robust)
        );

        Ok(())
    }

    /// The initial unweighted solve is dragged by correspondences into regions the other mesh
    /// never covered, which is why `max_distance` is required rather than optional.
    ///
    /// This is the unit-scale form of a failure found on real data: seventeen partially
    /// overlapping laser scans, starting from an already correct alignment, were moved most of a
    /// millimetre by a handful of such correspondences, and the robust refinement that followed
    /// entrenched the error rather than undoing it. See the module documentation.
    #[test]
    fn without_a_distance_gate_the_unweighted_solve_is_dragged() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(), Iso3::identity()];
        let ams = alignment_meshes(&meshes, &initial);

        // Both bodies start on the answer, so any motion at all is damage. A twentieth of the
        // correspondences are thrown well clear of the surface, standing in for points sampled
        // over a region the two meshes never both saw.
        let mut points = chain_points(&meshes, 0.4);
        for (k, p) in points.iter_mut().enumerate() {
            if k.is_multiple_of(20) {
                p.mp.sp.point += Vector3::new(0.0, 0.0, 8.0);
            }
        }

        let unweighted = |gate: f64| MultiOptions3 {
            refinement_steps: 0,
            ..MultiOptions3::new(gate)
        };
        let drift = |o: &MultiOutcome3| {
            (o.alignment(1).full_transform().to_matrix() - Iso3::identity().to_matrix()).norm()
        };

        let ungated =
            multi_mesh_adjustment_with_points(&ams, points.clone(), 0, &unweighted(WIDE_GATE))?;
        let gated = multi_mesh_adjustment_with_points(&ams, points, 0, &unweighted(1.0))?;

        assert!(
            drift(&ungated) > 1e-3,
            "the ungated solve should have been pulled off the answer, but moved only {}",
            drift(&ungated)
        );
        assert!(
            drift(&gated) < drift(&ungated) * 0.1,
            "the gate should have held the solve in place: ungated {}, gated {}",
            drift(&ungated),
            drift(&gated)
        );

        Ok(())
    }

    #[test]
    fn the_distance_gate_excludes_far_correspondences() -> Result<()> {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(), Iso3::identity()];
        let ams = alignment_meshes(&meshes, &initial);

        let mut points = chain_points(&meshes, 0.5);
        let far = points.len() / 2;
        points[far].mp.sp.point += Vector3::new(0.0, 0.0, 50.0);

        // With a tight gate the distant point contributes nothing, so the fit is unchanged from
        // the identity it already sits at.
        let opts = MultiOptions3 {
            max_distance: 1.0,
            refinement_steps: 0,
            ..test_opts()
        };
        let outcome = multi_mesh_adjustment_with_points(&ams, points, 0, &opts)?;

        assert_relative_eq!(
            outcome.alignment(1).full_transform().to_matrix(),
            Iso3::identity().to_matrix(),
            epsilon = 1e-6
        );

        Ok(())
    }

    // ============================================================================================
    // The jacobian
    // ============================================================================================

    #[test]
    fn jacobian_matches_finite_differences() -> Result<()> {
        // The recovery tests above are insensitive to jacobian errors: a wrong jacobian usually
        // costs Levenberg-Marquardt iterations rather than moving the minimum, so the solve still
        // lands in the right place. This checks the derivative itself.
        //
        // Three bodies are needed, not two. With two, the only reference is the static body, whose
        // jacobian block is dropped anyway, so the reverse derivative would never be exercised.
        // Here body 1 is a reference for body 2's correspondences *and* is free to move, which is
        // the configuration that makes the reverse block load-bearing.
        let meshes = vec![box_mesh(), box_mesh(), box_mesh()];
        let initial = vec![
            Iso3::identity(),
            Iso3::new(Vector3::new(0.21, -0.13, 0.09), Vector3::z() * 0.03),
            Iso3::new(Vector3::new(-0.17, 0.11, -0.07), Vector3::x() * 0.04),
        ];
        let ams = alignment_meshes(&meshes, &initial);

        // Sample points well inside the top face of each body, so every projection lands in the
        // interior of a face. The analytic jacobian holds the correspondence fixed, which is exact
        // for a locally planar target but not across an edge, and near an edge the finite
        // difference would pick up a correspondence flip that the analytic form does not model.
        let mut points = Vec::new();
        for body in [1usize, 2] {
            for i in -2..=2 {
                for j in -1..=1 {
                    let probe = Point3::new(i as f64 * 1.5, j as f64 * 1.5, 2.0);
                    let mp = meshes[body].surface_closest_to(&probe);
                    points.push(MeshAlignPoint::new(body, mp, body - 1, 1.0));
                }
            }
        }

        let opts = MultiOptions3 {
            refinement_steps: 0,
            ..test_opts()
        };
        let mut problem = MultiMeshProblem::new(&ams, points, 0, &opts)?;

        let x0 = problem.params.storage().clone();
        let analytic = problem.jacobian().unwrap();
        assert_eq!(analytic.ncols(), 12, "two free bodies, six parameters each");

        let eps = 1e-7;
        for k in 0..x0.len() {
            let mut lo = x0.clone();
            lo[k] -= eps;
            problem.set_params(&lo);
            let r_lo = problem.residuals().unwrap();

            let mut hi = x0.clone();
            hi[k] += eps;
            problem.set_params(&hi);
            let r_hi = problem.residuals().unwrap();

            for row in 0..r_lo.len() {
                let numeric = (r_hi[row] - r_lo[row]) / (2.0 * eps);
                assert_relative_eq!(analytic[(row, k)], numeric, epsilon = 1e-5);
            }
        }

        Ok(())
    }

    #[test]
    fn a_reference_body_receives_jacobian_columns() -> Result<()> {
        // A narrower statement of the same thing, which fails loudly if the reverse block is ever
        // dropped: correspondences sourced from body 2 must produce nonzero derivatives in body
        // 1's columns, because body 1 is what they are matched against.
        let meshes = vec![box_mesh(), box_mesh(), box_mesh()];
        let initial = vec![
            Iso3::identity(),
            Iso3::translation(0.2, 0.0, 0.0),
            Iso3::translation(0.0, 0.2, 0.0),
        ];
        let ams = alignment_meshes(&meshes, &initial);

        let mp = meshes[2].surface_closest_to(&Point3::new(0.0, 0.0, 2.0));
        let points = vec![MeshAlignPoint::new(2, mp, 1, 1.0)];

        let problem = MultiMeshProblem::new(&ams, points, 0, &test_opts())?;
        let jac = problem.jacobian().unwrap();

        let body1 = problem.params.column_offset(1).unwrap();
        let body2 = problem.params.column_offset(2).unwrap();

        let block_norm =
            |start: usize| -> f64 { (0..6).map(|k| jac[(0, start + k)].abs()).sum::<f64>() };

        assert!(
            block_norm(body2) > 1e-9,
            "the test body should have a nonzero jacobian block"
        );
        assert!(
            block_norm(body1) > 1e-9,
            "the reference body should have a nonzero jacobian block too"
        );

        Ok(())
    }

    #[test]
    fn a_self_correspondence_cancels_to_a_zero_row() -> Result<()> {
        // A point matched against its own mesh carries no information: moving that mesh moves the
        // test point and its match together, so the distance between them cannot change. The two
        // jacobian blocks land in the same columns and must cancel exactly, which they do because
        // the reverse form is the negation of the forward one when both share a body.
        //
        // This is why the blocks are accumulated rather than assigned. Overwriting would keep only
        // the reverse block and report a spurious sensitivity for a correspondence that has none.
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(), Iso3::translation(0.3, -0.2, 0.1)];
        let ams = alignment_meshes(&meshes, &initial);

        let mp = meshes[1].surface_closest_to(&Point3::new(1.0, 1.0, 2.0));
        let points = vec![MeshAlignPoint::new(1, mp, 1, 1.0)];

        let problem = MultiMeshProblem::new(&ams, points, 0, &test_opts())?;
        let jac = problem.jacobian().unwrap();

        for k in 0..jac.ncols() {
            assert_relative_eq!(jac[(0, k)], 0.0, epsilon = 1e-12);
        }

        Ok(())
    }

    // ============================================================================================
    // Validation
    // ============================================================================================

    #[test]
    fn invalid_options_and_indices_are_rejected() {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(); 2];
        let ams = alignment_meshes(&meshes, &initial);
        let points = chain_points(&meshes, 1.0);

        let bad_patience = MultiOptions3 {
            patience: 0,
            ..test_opts()
        };
        assert!(multi_mesh_adjustment_with_points(&ams, points.clone(), 0, &bad_patience).is_err());

        let bad_sigma = MultiOptions3 {
            sigma_max: Some(-1.0),
            ..test_opts()
        };
        assert!(multi_mesh_adjustment_with_points(&ams, points.clone(), 0, &bad_sigma).is_err());

        let bad_distance = MultiOptions3 {
            max_distance: 0.0,
            ..test_opts()
        };
        assert!(multi_mesh_adjustment_with_points(&ams, points.clone(), 0, &bad_distance).is_err());

        // A static index past the end of the mesh list.
        assert!(multi_mesh_adjustment_with_points(&ams, points, 5, &test_opts()).is_err());
    }

    #[test]
    fn out_of_range_alignment_points_are_rejected() {
        let meshes = vec![box_mesh(), box_mesh()];
        let initial = vec![Iso3::identity(); 2];
        let ams = alignment_meshes(&meshes, &initial);

        let mut points = chain_points(&meshes, 1.0);
        points[0].ref_i = 7;

        let err = multi_mesh_adjustment_with_points(&ams, points, 0, &test_opts())
            .unwrap_err()
            .to_string();
        assert!(err.contains("only 2 meshes"), "unexpected message: {err}");
    }
}
