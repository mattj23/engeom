//! This module holds the facade that mesh editing operations are driven through.
//!
//! Editing a mesh means working on its topology, which is what a half-edge structure is for, while
//! `Mesh3` is the type with the acceleration structure needed for geometric queries. Converting
//! between them is not free: each direction rebuilds the whole structure, and going back to `Mesh3`
//! rebuilds its BVH.
//!
//! So an editor is built once, holds the half-edge mesh across a whole chain of operations, and
//! converts back once at the end. Doing the same work as a sequence of `Mesh3 -> Mesh3` methods
//! would pay that cost at every step.
//!
//! The editor also keeps a borrow of the mesh it started from, for accuracy checks against the
//! original shape. Note that guaranteed decimation no longer needs it: it carries its own
//! per-vertex error volume on the half-edge mesh, so repeated runs are anchored to the as-measured
//! surface rather than to an intermediate. Operations added later which do need a reference should
//! use this borrow rather than the current state, for the same reason.
//!
//! <div class="warning">
//!
//! **What that error volume survives is not settled.** It carries across repeated
//! `decimate_guaranteed` calls, which is the case it was built for. It demonstrably does not
//! survive [`MeshEditor::decimate_best_effort`], which is why that combination is refused rather
//! than merely documented. Whether it survives [`MeshEditor::smooth`] is an open question and
//! probably it does not, since smoothing moves vertices without widening the radii attached to
//! them. The guard and the reasoning both live on `HalfEdgeMesh3::error_volume_stale`.
//!
//! </div>

use crate::geom3::half_edge3::{
    BestEffortOpts, DecimateOpts, DecimateReport, HalfEdgeMesh3, RepairOpts, RepairReport,
};
use crate::{Mesh3, Result};

/// An account of everything an editing session changed.
///
/// Marked `#[non_exhaustive]`, since operations added later will report into it.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EditReport {
    /// What the ingest had to repair to make the mesh usable, if it was built with a repair.
    pub repair: RepairReport,

    /// How many smoothing passes were applied.
    pub smoothing_passes: usize,

    /// How many edge collapses decimation performed, summed over every run.
    ///
    /// Pooled across both decimation paths, since this counts work done rather than what was
    /// promised about it.
    pub decimate_collapses: usize,
}

impl std::fmt::Display for EditReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ingest: {}", self.repair)?;
        if self.smoothing_passes > 0 {
            write!(f, "; {} smoothing passes", self.smoothing_passes)?;
        }
        if self.decimate_collapses > 0 {
            write!(f, "; {} collapses", self.decimate_collapses)?;
        }
        Ok(())
    }
}

/// A handle for performing a sequence of topology edits on a mesh.
///
/// Build one from a `Mesh3`, chain operations on it, then convert back with
/// [`MeshEditor::into_mesh`]. See the module documentation for why this is a persistent handle
/// rather than a set of methods on `Mesh3`.
///
/// ```
/// use engeom::Mesh3;
/// use engeom::geom3::mesh::MeshEditor;
/// use engeom::geom3::half_edge3::RepairOpts;
///
/// let mesh = Mesh3::stanford_bunny_res4();
/// let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default())?;
/// editor.smooth(2)?;
/// let cleaned = editor.into_mesh(false, true)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct MeshEditor<'a> {
    half_edge: HalfEdgeMesh3,
    original: &'a Mesh3,
    report: EditReport,
}

impl<'a> MeshEditor<'a> {
    /// Start an editing session, failing on the first element the half-edge structure will not
    /// accept.
    ///
    /// This is the strict path, appropriate when the mesh is supposed to be clean. Use
    /// [`MeshEditor::with_repair`] on measured data, which is routinely non-manifold.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh to edit, borrowed for the life of the session
    ///
    /// returns: `Result<MeshEditor>`
    pub fn new(mesh: &'a Mesh3) -> Result<Self> {
        Ok(Self {
            half_edge: HalfEdgeMesh3::try_from(mesh)?,
            original: mesh,
            report: EditReport::default(),
        })
    }

    /// Start an editing session, repairing whatever topology the half-edge structure will not
    /// accept.
    ///
    /// What the repair had to change is recorded in the session's report, under
    /// `EditReport::repair`. Check it rather than assuming the mesh came through untouched.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh to edit, borrowed for the life of the session
    /// * `opts`: which repairs to attempt
    ///
    /// returns: `Result<MeshEditor>`
    pub fn with_repair(mesh: &'a Mesh3, opts: &RepairOpts) -> Result<Self> {
        let (half_edge, repair) = HalfEdgeMesh3::from_mesh_repaired(mesh, opts)?;

        Ok(Self {
            half_edge,
            original: mesh,
            report: EditReport {
                repair,
                ..Default::default()
            },
        })
    }

    /// The mesh this session started from, which is what operations measure deviation against.
    pub fn original(&self) -> &'a Mesh3 {
        self.original
    }

    /// What this session has changed so far.
    pub fn report(&self) -> &EditReport {
        &self.report
    }

    /// Borrow the underlying half-edge mesh, for operations this facade does not wrap.
    pub fn half_edge(&self) -> &HalfEdgeMesh3 {
        &self.half_edge
    }

    /// Mutably borrow the underlying half-edge mesh. See [`MeshEditor::half_edge`].
    ///
    /// Anything done through this is invisible to the session's report, so a caller which edits
    /// here is responsible for describing what it did.
    pub fn half_edge_mut(&mut self) -> &mut HalfEdgeMesh3 {
        &mut self.half_edge
    }

    /// The current vertex count, including any marked deleted but not yet collected.
    pub fn vertex_count(&self) -> usize {
        self.half_edge.vertex_count()
    }

    /// The current face count, including any marked deleted but not yet collected.
    pub fn face_count(&self) -> usize {
        self.half_edge.face_count()
    }

    /// Smooth the mesh, fitting a plane to each vertex's one-ring and moving the vertex onto the
    /// average height of that ring.
    ///
    /// This moves measured points and is not bounded by any tolerance, so it changes the geometry
    /// rather than only the topology.
    ///
    /// # Arguments
    ///
    /// * `iterations`: how many passes to apply. Zero is a no-op.
    ///
    /// returns: `Result<&mut Self>`, to allow chaining
    pub fn smooth(&mut self, iterations: usize) -> Result<&mut Self> {
        for _ in 0..iterations {
            self.half_edge.neighborhood_smooth()?;
        }

        self.report.smoothing_passes += iterations;
        Ok(self)
    }

    /// Decimate the mesh, keeping both surfaces within the tolerance of each other.
    ///
    /// The bound holds over the surfaces in both directions, and it holds across a chain rather
    /// than per call: the error volume accumulates on the half-edge mesh, so several runs at a
    /// given tolerance land inside that tolerance rather than inside the sum of them. See the
    /// `half_edge3::decimate` module documentation for the mechanism.
    ///
    /// # Errors
    ///
    /// Refuses if [`MeshEditor::decimate_best_effort`] has already run in this session. That path
    /// leaves the error volume describing less than what has happened to the surface, so a bound
    /// computed from it afterwards would not be one. Convert with [`MeshEditor::to_mesh`] and start
    /// a fresh session if the best-effort result is deliberately the new reference. The refusal
    /// comes from [`HalfEdgeMesh3::decimate_guaranteed`], which is where the error volume lives.
    ///
    /// # Arguments
    ///
    /// * `opts`: how decimation should behave
    ///
    /// returns: `Result<DecimateReport>`
    pub fn decimate_guaranteed(&mut self, opts: &DecimateOpts) -> Result<DecimateReport> {
        let report = self.half_edge.decimate_guaranteed(opts)?;
        self.report.decimate_collapses += report.collapses;
        Ok(report)
    }

    /// Decimate the mesh by estimated deviation rather than by a guaranteed bound.
    ///
    /// Faster and considerably more aggressive than
    /// [`decimate_guaranteed`](MeshEditor::decimate_guaranteed), and the tolerance it takes is an
    /// estimate which the result routinely exceeds. Use it for display meshes and for geometry
    /// which scaffolds a measurement rather than carrying one.
    ///
    /// **This does not accumulate a bound the way the guaranteed path does.** It holds every
    /// *original vertex* inside the tolerance and says nothing about the surface between them, so
    /// a chain which mixes the two is only as strong as the weaker one. Once this has run, the
    /// error volume a later `decimate_guaranteed` reasons from no longer describes everything that
    /// has happened to the surface. See the `half_edge3::decimate::best_effort` module
    /// documentation for what the estimate is and the ways it can be wrong.
    ///
    /// # Arguments
    ///
    /// * `opts`: how decimation should behave
    ///
    /// returns: `Result<DecimateReport>`
    pub fn decimate_best_effort(&mut self, opts: &BestEffortOpts) -> Result<DecimateReport> {
        let report = self.half_edge.decimate_best_effort(opts)?;
        self.report.decimate_collapses += report.collapses;
        Ok(report)
    }

    /// Build a `Mesh3` from the session's current state, leaving the session usable.
    ///
    /// This compacts the half-edge structure first, discarding elements marked deleted, so any
    /// handle or index taken from it beforehand is invalid afterwards.
    ///
    /// # Arguments
    ///
    /// * `is_solid`: whether the result should answer distance queries as a solid
    /// * `allow_attribute_loss`: whether to accept losing the original's attributes
    ///
    /// returns: `Result<Mesh3>`, failing if the original carried attributes and the loss was not
    /// accepted
    pub fn to_mesh(&mut self, is_solid: bool, allow_attribute_loss: bool) -> Result<Mesh3> {
        // The half-edge structure has no index mapping back to the original's points or faces, so
        // there is nothing to carry attributes on. This is the same bargain `convex_hull` makes.
        self.original
            .check_attribute_loss("a mesh editing session", allow_attribute_loss)?;

        self.half_edge.to_mesh(is_solid)
    }

    /// Consume the session and build a `Mesh3` from its final state.
    ///
    /// See [`MeshEditor::to_mesh`] for the arguments.
    pub fn into_mesh(mut self, is_solid: bool, allow_attribute_loss: bool) -> Result<Mesh3> {
        self.to_mesh(is_solid, allow_attribute_loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::dist;

    #[test]
    fn a_clean_mesh_round_trips_unchanged() {
        let mesh = crate::tests::engine_blade();
        let editor = MeshEditor::new(&mesh).unwrap();

        assert_eq!(editor.face_count(), mesh.faces().len());
        assert!(editor.report().repair.is_clean());

        let out = editor.into_mesh(false, false).unwrap();
        assert_eq!(out.faces().len(), mesh.faces().len());
        assert_eq!(out.points().len(), mesh.points().len());
    }

    /// The strict constructor refuses what the repairing one accepts, and the repair is reported
    /// rather than hidden.
    #[test]
    fn the_repairing_constructor_takes_a_dirty_mesh() {
        let mesh = Mesh3::stanford_bunny_res4();

        let editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();
        let repair = editor.report().repair;

        assert!(!repair.is_clean());
        assert_eq!(repair.faces_rejected_by_ingest, 0);
        assert_eq!(editor.face_count(), 868);
    }

    #[test]
    fn smoothing_is_counted_and_moves_points() {
        let mesh = Mesh3::stanford_bunny_res4();
        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();

        let before = editor.half_edge().clone_vertices().unwrap();
        editor.smooth(2).unwrap();
        let after = editor.half_edge().clone_vertices().unwrap();

        assert_eq!(editor.report().smoothing_passes, 2);
        assert_eq!(before.len(), after.len());

        let moved = before
            .iter()
            .zip(after.iter())
            .filter(|(a, b)| dist(*a, *b) > 1e-12)
            .count();
        assert!(moved > 0, "smoothing should have moved some points");
    }

    #[test]
    fn zero_smoothing_passes_change_nothing() {
        let mesh = crate::tests::engine_blade();
        let mut editor = MeshEditor::new(&mesh).unwrap();

        let before = editor.half_edge().clone_vertices().unwrap();
        editor.smooth(0).unwrap();
        let after = editor.half_edge().clone_vertices().unwrap();

        assert_eq!(editor.report().smoothing_passes, 0);
        assert_eq!(before, after);
    }

    #[test]
    fn operations_chain() {
        let mesh = Mesh3::stanford_bunny_res4();
        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();

        editor.smooth(1).unwrap().smooth(1).unwrap();

        assert_eq!(editor.report().smoothing_passes, 2);
    }

    /// Both decimation paths are reachable from the editor and report into the same counter.
    ///
    /// The best-effort path is the more aggressive of the two at a given tolerance, which is the
    /// only reason to reach for it, so that ordering is asserted rather than just the plumbing.
    #[test]
    fn both_decimation_paths_are_reachable() {
        const TOL: f64 = 0.01;

        let mesh = Mesh3::stanford_bunny_res4();

        let mut guaranteed = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();
        let before = guaranteed.face_count();
        let a = guaranteed
            .decimate_guaranteed(&DecimateOpts::to_tolerance(TOL))
            .unwrap();

        let mut best = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();
        let b = best
            .decimate_best_effort(&BestEffortOpts::to_tolerance(TOL))
            .unwrap();

        for (report, editor) in [(a, &guaranteed), (b, &best)] {
            assert!(
                report.collapses > 0,
                "nothing collapsed at a tolerance of {TOL}"
            );
            assert!(editor.face_count() < before);
            assert_eq!(editor.report().decimate_collapses, report.collapses);
        }

        assert!(
            best.face_count() < guaranteed.face_count(),
            "best effort should be the more aggressive path: {} against {}",
            best.face_count(),
            guaranteed.face_count()
        );
    }

    /// A bound cannot be computed from an error volume that has stopped describing the surface.
    ///
    /// The refusal is one-directional on purpose: best effort makes no claim, so it has nothing to
    /// lose by following a guaranteed pass, and only the reverse order is unsound.
    #[test]
    fn a_guarantee_cannot_follow_a_best_effort_pass() {
        let mesh = Mesh3::stanford_bunny_res4();

        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();
        editor
            .decimate_best_effort(&BestEffortOpts::to_tolerance(0.01))
            .unwrap();

        let err = editor
            .decimate_guaranteed(&DecimateOpts::to_tolerance(0.01))
            .expect_err("a guaranteed pass should refuse to follow a best-effort one");
        assert!(
            err.to_string().contains("best-effort"),
            "the refusal should say why: {err}"
        );

        // The other order is fine, and so is repeating either path on its own.
        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();
        editor
            .decimate_guaranteed(&DecimateOpts::to_tolerance(0.01))
            .unwrap();
        editor
            .decimate_guaranteed(&DecimateOpts::to_tolerance(0.01))
            .unwrap();
        editor
            .decimate_best_effort(&BestEffortOpts::to_tolerance(0.01))
            .unwrap();
    }

    /// The original is held for the whole session, so a later operation still measures against the
    /// mesh that was handed in rather than against what earlier operations left behind.
    #[test]
    fn the_original_is_the_mesh_the_session_started_from() {
        let mesh = Mesh3::stanford_bunny_res4();
        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();

        editor.smooth(3).unwrap();

        assert_eq!(editor.original().faces().len(), mesh.faces().len());
        assert_eq!(editor.original().points().len(), mesh.points().len());
    }

    /// Attributes cannot survive the half-edge round trip, so losing them has to be opted into,
    /// matching what `convex_hull` already does.
    #[test]
    fn attributes_must_be_given_up_explicitly() {
        let mut mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let labels: Vec<u32> = (0..mesh.faces().len() as u32).collect();
        mesh.set_face_labels(Some(labels)).unwrap();

        let mut editor = MeshEditor::new(&mesh).unwrap();
        assert!(
            editor.to_mesh(false, false).is_err(),
            "dropping attributes should require an opt-in"
        );

        let out = editor.to_mesh(false, true).unwrap();
        assert!(out.face_labels().is_none());
    }

    /// A mesh with no attributes has nothing to lose, so the opt-in is not required.
    #[test]
    fn a_mesh_without_attributes_needs_no_opt_in() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let editor = MeshEditor::new(&mesh).unwrap();

        assert!(editor.into_mesh(false, false).is_ok());
    }

    /// `to_mesh` leaves the session usable, so the state can be inspected part way through a chain.
    #[test]
    fn to_mesh_leaves_the_session_usable() {
        let mesh = Mesh3::stanford_bunny_res4();
        let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default()).unwrap();

        let first = editor.to_mesh(false, true).unwrap();
        editor.smooth(1).unwrap();
        let second = editor.to_mesh(false, true).unwrap();

        assert_eq!(first.faces().len(), second.faces().len());
        assert_eq!(editor.report().smoothing_passes, 1);
    }

    #[test]
    fn is_solid_is_carried_onto_the_result() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);

        let solid = MeshEditor::new(&mesh)
            .unwrap()
            .into_mesh(true, false)
            .unwrap();
        assert!(solid.is_solid());

        let hollow = MeshEditor::new(&mesh)
            .unwrap()
            .into_mesh(false, false)
            .unwrap();
        assert!(!hollow.is_solid());
    }
}
