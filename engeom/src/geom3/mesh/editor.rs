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
//! The editor also keeps a borrow of the mesh it started from, and this is for doing fast accuracy
//! checks against the original shape. Operations which need to know how far the surface has moved
//! _must_ measure against that original, so across a chain such as
//! `repair -> decimate -> smooth -> decimate` the deviation is always taken from the as-measured
//! surface and never from an already-degraded intermediate, otherwise error will compound.

use crate::geom3::mesh::half_edge::{HalfEdgeMesh, RepairOpts, RepairReport};
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
}

impl std::fmt::Display for EditReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ingest: {}", self.repair)?;
        if self.smoothing_passes > 0 {
            write!(f, "; {} smoothing passes", self.smoothing_passes)?;
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
/// use engeom::geom3::mesh::half_edge::RepairOpts;
///
/// let mesh = Mesh3::stanford_bunny_res4();
/// let mut editor = MeshEditor::with_repair(&mesh, &RepairOpts::default())?;
/// editor.smooth(2)?;
/// let cleaned = editor.into_mesh(false, true)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct MeshEditor<'a> {
    half_edge: HalfEdgeMesh,
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
            half_edge: HalfEdgeMesh::try_from(mesh)?,
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
        let (half_edge, repair) = HalfEdgeMesh::from_mesh_repaired(mesh, opts)?;

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
    pub fn half_edge(&self) -> &HalfEdgeMesh {
        &self.half_edge
    }

    /// Mutably borrow the underlying half-edge mesh. See [`MeshEditor::half_edge`].
    ///
    /// Anything done through this is invisible to the session's report, so a caller which edits
    /// here is responsible for describing what it did.
    pub fn half_edge_mut(&mut self) -> &mut HalfEdgeMesh {
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

        self.half_edge.garbage_collect()?;

        let vertices = self.half_edge.clone_vertices()?;
        let faces = self.half_edge.clone_faces()?;

        if faces.is_empty() {
            return Err("The editing session left no faces, and a mesh needs at least one".into());
        }

        Ok(Mesh3::new(vertices, faces, is_solid))
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
