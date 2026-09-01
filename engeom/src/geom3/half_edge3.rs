//! This module provides a half-edge mesh structure, built on the `alum` library.
//!
//! It sits beside [`mesh`](crate::geom3::mesh) rather than inside it. The two are alternative
//! representations of the same thing, a connectivity-first one here and a triangle-soup one there,
//! and despite their interconnectedness, neither is technically part of the other, or part of
//! a common larger hierarchy.
//!
//! Conversion runs both ways through `TryFrom`, and [`MeshEditor`](crate::geom3::mesh::MeshEditor)
//! is the handle that drives an edit session from the `Mesh3` side.
//!
//! The main type is [`HalfEdgeMesh3`], which owns an `alum` polygon mesh rather than being an alias
//! for one like it was in earlier versions of `engeom`.  We went away from that for two reasons:
//!
//! The first is that Rust's orphan rule forbids an inherent `impl` on a foreign type, so while
//! this was an alias every operation had to arrive as a separate extension trait which the caller
//! then had to import, which was cumbersome.
//!
//! The second is the public surface; an alias makes every `alum` method, associated type, and
//! error variant part of this crate's public API, so any change upstream is a potential breaking
//! change here. This hasn't been a heavily developed part of the library, and during the time
//! between when I originally authored it and revisiting it now, `alum` went from 0.6 to 0.7 and
//! added some new requirements that broke the build.  By wrapping the type I'm putting a seam here
//! that can absorb some of that churn.
//!
//! That said, because this is one of the least trodden parts of `engeom`, I don't have a good feel
//! for how tooling around it is going to work in practice, so [`HalfEdgeMesh3::as_alum`],
//! [`HalfEdgeMesh3::as_alum_mut`], and [`HalfEdgeMesh3::into_alum`] hand out the underlying mesh for
//! anything this wrapper does not cover, including the whole of `alum`'s decimation and subdivision
//! machinery.
//!
//! [`NaAdaptor`] is the glue that lets `alum` work in terms of `nalgebra` types.

mod decimate;
#[cfg(test)]
pub(crate) mod difference;
mod repair;
mod smoothing;

use crate::{Mesh3, Point3, Result, Vector3};
use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use decimate::{
    BestEffortDecimator, BestEffortOpts, DecimateOpts, DecimateReport, DecimateStats,
    DecimateTarget, ErrorMethod, Placement, QuadricKind, ToleranceVolumeDecimator,
};
pub use repair::{RepairOpts, RepairReport, Repaired, repair_buffers};
use std::collections::HashMap;
use std::error::Error;

/// The `alum` mesh type that [`HalfEdgeMesh3`] wraps.
///
/// This is public so the escape-hatch accessors have a nameable return type. Naming it in your own
/// code is the explicit decision to depend on `alum`'s API, with the stability that implies; the
/// wrapper's own methods are this crate's to keep stable.
pub type AlumMesh = alum::PolyMeshT<3, NaAdaptor>;

/// A half-edge mesh, which stores connectivity explicitly rather than as a soup of triangles.
///
/// This is the representation to reach for when an operation needs to navigate topology, such as
/// walking the ring of faces around a vertex or collapsing an edge. `Mesh3` remains the right type
/// for geometric queries, since it carries a BVH and this does not.
///
/// Build one from a `Mesh3` with [`HalfEdgeMesh3::from_mesh_repaired`] for measured data, or with
/// `TryFrom` when the input is known to be clean. Neither direction carries attributes, `is_solid`,
/// or the UV mapping, because this structure has nowhere to put them.
pub struct HalfEdgeMesh3 {
    inner: AlumMesh,

    /// Per-vertex scalars which have to outlive a single operation, keyed by name.
    ///
    /// `alum`'s properties are anonymous: creating one hands back a fresh buffer with no way to
    /// ask the mesh whether an equivalent already exists. Anything which needs to accumulate
    /// across successive calls therefore has to keep its own handle, and this is where those
    /// handles live. Decimation's error volume is the motivating case, since it is only a bound on
    /// deviation from the true original if a second run continues from what the first one left
    /// rather than resetting and measuring against that run's own output.
    ///
    /// The properties are registered with the mesh, so they are remapped by
    /// [`HalfEdgeMesh3::garbage_collect`] like any other.
    vertex_scalars: HashMap<String, alum::VProperty<f64>>,

    /// Whether an operation has run whose effect the decimation error volume does not account for.
    ///
    /// [`HalfEdgeMesh3::decimate_guaranteed`] refuses once this is set, because the volume it
    /// reasons from would no longer contain the original surface and the bound it returned would
    /// not be a bound.
    ///
    /// This lives here rather than on [`MeshEditor`](crate::geom3::mesh::MeshEditor) because the
    /// error volume lives here: whether it still describes the surface is a fact about the mesh,
    /// not about a session wrapping it. The editor reports it, it does not own it.
    ///
    // TODO: state accumulation across a chain has not been thought through formally, and this flag
    // is a guard on the one case that is definitely wrong rather than a model of the problem.
    // `decimate_best_effort` sets it because it makes no claim about the surface between the
    // original vertices. Two other operations look like the same defect and currently do not:
    //
    // - `neighborhood_smooth` moves vertices while their error radii stay put, so the volume stops
    //   containing the original surface just as surely. This is the one worth fixing first, and
    //   fixing it means widening rather than refusing: the pass is two-phase and already holds
    //   every vertex displacement, and growing each radius by its own vertex's displacement
    //   satisfies both invariants exactly. A point x with weights w maps to y with
    //   |y - x| <= sum(w_k * d_k) while e_new(y) = e_old(x) + sum(w_k * d_k), so the new ball
    //   contains the old one and vice versa.
    // - `as_alum_mut` hands out arbitrary mutation and is invisible here by construction.
    //
    // Guarding only one of the three implies the other two are safe, which is not established.
    // What is needed is a statement of what each operation promises about the error volume, and
    // then either a guard per operation or a way for an operation to seed the volume to cover what
    // it did. Gueziec sanctions the latter: the nesting proof's base case needs a *valid* starting
    // volume, not a zero one, so a seeded volume is as good as an accumulated one.
    error_volume_stale: bool,
}

impl Default for HalfEdgeMesh3 {
    fn default() -> Self {
        Self::new()
    }
}

impl HalfEdgeMesh3 {
    /// Create an empty half-edge mesh.
    pub fn new() -> Self {
        Self {
            inner: AlumMesh::new(),
            vertex_scalars: HashMap::new(),
            error_volume_stale: false,
        }
    }

    /// Fetch a named per-vertex scalar property, creating it if this is the first request.
    ///
    /// An existing property is returned exactly as it was, carrying whatever previous
    /// operations left in it. A new one is filled with `default`.
    ///
    /// # Arguments
    ///
    /// * `name`: the key the property is remembered under, which callers should namespace
    /// * `default`: the value a newly created property takes, and the value `alum` gives to
    ///   vertices added later
    ///
    /// returns: `Result<alum::VProperty<f64>>`
    pub fn vertex_prop_f64(&mut self, name: &str, default: f64) -> Result<alum::VProperty<f64>> {
        if let Some(p) = self.vertex_scalars.get(name) {
            return Ok(p.clone());
        }

        let prop = self.inner.create_vertex_prop(default);
        self.vertex_scalars.insert(name.to_string(), prop.clone());
        Ok(prop)
    }

    /// Borrow the underlying `alum` mesh.
    ///
    /// This is the escape hatch for anything the wrapper does not expose, `alum`'s decimation and
    /// subdivision traits in particular. Using it couples your code to `alum`'s API rather than to
    /// this crate's.
    pub fn as_alum(&self) -> &AlumMesh {
        &self.inner
    }

    /// Mutably borrow the underlying `alum` mesh. See [`HalfEdgeMesh3::as_alum`].
    pub fn as_alum_mut(&mut self) -> &mut AlumMesh {
        &mut self.inner
    }

    /// Consume the wrapper and return the underlying `alum` mesh.
    pub fn into_alum(self) -> AlumMesh {
        self.inner
    }

    /// Wrap an `alum` mesh which was built elsewhere.
    ///
    /// Any named vertex scalars the mesh carried are lost, since `alum` gives no way to recover a
    /// handle to an anonymous property. Round-tripping a mesh through `alum` therefore resets a
    /// decimation error volume to nothing, which reads as "this mesh is its own reference".
    pub fn from_alum(inner: AlumMesh) -> Self {
        Self {
            inner,
            vertex_scalars: HashMap::new(),
            error_volume_stale: false,
        }
    }

    /// Whether the decimation error volume still accounts for everything done to this mesh.
    ///
    /// False once [`HalfEdgeMesh3::decimate_best_effort`] has run, after which
    /// [`HalfEdgeMesh3::decimate_guaranteed`] refuses. See the field of the same name for what is
    /// and is not covered by this.
    pub fn error_volume_is_current(&self) -> bool {
        !self.error_volume_stale
    }

    /// Record that an operation has invalidated the decimation error volume.
    ///
    /// Call this after mutating the mesh through [`HalfEdgeMesh3::as_alum_mut`] in a way that moves
    /// the surface, so that a later guaranteed decimation refuses rather than returning a bound it
    /// cannot support.
    pub fn invalidate_error_volume(&mut self) {
        self.error_volume_stale = true;
    }

    pub(crate) fn set_error_volume_stale(&mut self) {
        self.error_volume_stale = true;
    }

    pub(crate) fn error_volume_stale(&self) -> bool {
        self.error_volume_stale
    }

    /// Build a `Mesh3` from the current state.
    ///
    /// This compacts the structure first, discarding elements marked deleted, so any handle or
    /// index taken from it beforehand is invalid afterwards.
    ///
    /// Nothing is carried but geometry. A half-edge mesh has nowhere to hold attributes, `is_solid`
    /// or the UV mapping, so a caller which started from a `Mesh3` carrying any of those is losing
    /// them here; [`MeshEditor::to_mesh`](crate::geom3::mesh::MeshEditor::to_mesh) is the entry
    /// point which checks for that, because it is the one which still has the original to check
    /// against.
    ///
    /// # Arguments
    ///
    /// * `is_solid`: whether the result should answer distance queries as a solid
    ///
    /// returns: `Result<Mesh3>`, failing if the structure has no faces left
    pub fn to_mesh(&mut self, is_solid: bool) -> Result<Mesh3> {
        self.garbage_collect()?;

        let vertices = self.clone_vertices()?;
        let faces = self.clone_faces()?;

        if faces.is_empty() {
            return Err(
                "The half-edge mesh has no faces left, and a mesh needs at least one".into(),
            );
        }

        Ok(Mesh3::new(vertices, faces, is_solid))
    }

    /// The number of vertices, including any which are deleted but not yet collected.
    pub fn vertex_count(&self) -> usize {
        self.inner.num_vertices()
    }

    /// The number of faces, including any which are deleted but not yet collected.
    pub fn face_count(&self) -> usize {
        self.inner.num_faces()
    }

    /// The number of edges, including any which are deleted but not yet collected.
    pub fn edge_count(&self) -> usize {
        self.inner.num_edges()
    }

    /// Compact the structure, discarding elements which were marked deleted.
    ///
    /// Element indices are not stable across this, so anything holding handles or indices into the
    /// mesh must re-derive them afterwards.
    pub fn garbage_collect(&mut self) -> Result<()> {
        self.inner
            .garbage_collection()
            .map_err(|e| format!("Garbage collection failed: {:?}", e).into())
    }

    /// Copy the vertex positions out into a plain buffer, in vertex order.
    pub fn clone_vertices(&self) -> Result<Vec<Point3>> {
        let points = self.inner.points();
        let points = points.try_borrow().map_err(|_| "Failed to borrow points")?;
        Ok(points.iter().map(|v| Point3::from(*v)).collect())
    }

    /// Copy the faces out into a plain buffer as triples of vertex indices, triangulating anything
    /// with more than three sides.
    pub fn clone_faces(&self) -> Result<Vec<[u32; 3]>> {
        let f_status = self.inner.face_status_prop();
        let f_status = f_status
            .try_borrow()
            .map_err(|_| "Failed to borrow face status")?;
        Ok(self
            .inner
            .triangulated_vertices(&f_status)
            .map(|f| [f[0].index(), f[1].index(), f[2].index()])
            .collect())
    }
}

impl TryFrom<&HalfEdgeMesh3> for Mesh3 {
    type Error = Box<dyn Error>;

    fn try_from(value: &HalfEdgeMesh3) -> std::result::Result<Self, Self::Error> {
        let vertices = value.clone_vertices()?;
        let faces = value.clone_faces()?;
        Ok(Mesh3::new(vertices, faces, false))
    }
}

impl TryFrom<&Mesh3> for HalfEdgeMesh3 {
    type Error = Box<dyn Error>;

    /// Build a half-edge mesh from a `Mesh3`, failing on the first element the structure will not
    /// accept.
    ///
    /// This is the strict path, appropriate when the input is supposed to be clean. Use
    /// [`HalfEdgeMesh3::from_mesh_repaired`] on measured data, which is routinely non-manifold.
    ///
    /// Note that this is a test of what can be inserted, not of manifoldness: `alum` refuses an
    /// edge shared by three faces and refuses a fold, but accepts a bowtie whose fans are both
    /// open.
    fn try_from(value: &Mesh3) -> std::result::Result<Self, Self::Error> {
        let mut result = HalfEdgeMesh3::new();
        let mut indices = Vec::new();
        for v in value.points() {
            let handle = result
                .inner
                .add_vertex(v.coords)
                .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
            indices.push(handle);
        }

        for f in value.faces() {
            result
                .inner
                .add_tri_face(
                    indices[f[0] as usize],
                    indices[f[1] as usize],
                    indices[f[2] as usize],
                )
                .map_err(|e| format!("Failed to add face: {:?}", e))?;
        }

        Ok(result)
    }
}

/// Teaches `alum` to work in terms of `nalgebra`'s vector and scalar types.
pub struct NaAdaptor {}

impl Adaptor<3> for NaAdaptor {
    type Vector = Vector3;
    type Scalar = f64;

    fn vector(coords: [Self::Scalar; 3]) -> Self::Vector {
        Vector3::new(coords[0], coords[1], coords[2])
    }

    fn zero_vector() -> Self::Vector {
        Vector3::zeros()
    }

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
    }
}

impl VectorLengthAdaptor<3> for NaAdaptor {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        v.norm()
    }

    fn vector_length_squared(v: Self::Vector) -> Self::Scalar {
        v.norm_squared()
    }
}

impl VectorNormalizeAdaptor<3> for NaAdaptor {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for NaAdaptor {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(&b)
    }
}

impl VectorAngleAdaptor for NaAdaptor {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle(&b)
    }
}

impl CrossProductAdaptor for NaAdaptor {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(&b)
    }
}

impl FloatScalarAdaptor<3> for NaAdaptor {
    fn scalarf32(val: f32) -> Self::Scalar {
        val as f64
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_mesh_has_nothing_in_it() {
        let he = HalfEdgeMesh3::new();
        assert_eq!(he.vertex_count(), 0);
        assert_eq!(he.face_count(), 0);
        assert_eq!(he.edge_count(), 0);
        assert!(he.clone_vertices().unwrap().is_empty());
        assert!(he.clone_faces().unwrap().is_empty());
    }

    #[test]
    fn a_clean_mesh_round_trips() {
        let mesh = crate::tests::engine_blade();
        let he = HalfEdgeMesh3::try_from(&mesh).unwrap();

        assert_eq!(he.vertex_count(), mesh.points().len());
        assert_eq!(he.face_count(), mesh.faces().len());

        let back = Mesh3::try_from(&he).unwrap();
        assert_eq!(back.points().len(), mesh.points().len());
        assert_eq!(back.faces().len(), mesh.faces().len());
    }

    /// The escape hatch has to actually reach `alum`'s own traits, since that is the whole reason
    /// it exists. This is the path the decimator will take.
    /// The staleness guard is on the mesh, so the raw half-edge path is protected too rather than
    /// only the editor facade that wraps it.
    #[test]
    fn a_guarantee_cannot_follow_a_best_effort_pass() {
        use crate::geom3::half_edge3::{BestEffortOpts, DecimateOpts};

        let mesh = Mesh3::stanford_bunny_res4();
        let (mut he, _) = HalfEdgeMesh3::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();

        assert!(he.error_volume_is_current());
        he.decimate_best_effort(&BestEffortOpts::to_tolerance(0.01))
            .unwrap();
        assert!(!he.error_volume_is_current());

        let err = he
            .decimate_guaranteed(&DecimateOpts::to_tolerance(0.01))
            .expect_err("a guaranteed pass should refuse to follow a best-effort one");
        assert!(
            err.to_string().contains("best-effort"),
            "the refusal should say why: {err}"
        );
    }

    /// Rejected options must not mark the mesh, since nothing was done to it.
    #[test]
    fn a_refused_best_effort_call_leaves_the_volume_alone() {
        use crate::geom3::half_edge3::BestEffortOpts;

        let mesh = Mesh3::stanford_bunny_res4();
        let (mut he, _) = HalfEdgeMesh3::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();

        assert!(
            he.decimate_best_effort(&BestEffortOpts::to_tolerance(-1.0))
                .is_err()
        );
        assert!(he.error_volume_is_current());
    }

    #[test]
    fn the_escape_hatch_reaches_alum_decimation() {
        use alum::decimate::HasDecimation;
        use alum::decimate::quadric::{QuadricDecimater, QuadricType};

        let mesh = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&mesh).unwrap();
        let before = he.face_count();

        let mut decimator = QuadricDecimater::new(he.as_alum(), QuadricType::Triangle).unwrap();
        let collapses = he
            .as_alum_mut()
            .decimate_to_face_count(&mut decimator, before / 2)
            .unwrap();

        assert!(collapses > 0, "decimation should have collapsed something");

        he.garbage_collect().unwrap();
        assert!(
            he.face_count() < before,
            "face count should have dropped from {}",
            before
        );
    }

    #[test]
    fn wrapping_and_unwrapping_preserves_the_mesh() {
        let mesh = crate::tests::stanford_bun_4();
        let he = HalfEdgeMesh3::from_mesh_repaired(&mesh, &RepairOpts::default())
            .unwrap()
            .0;
        let faces = he.face_count();

        let rewrapped = HalfEdgeMesh3::from_alum(he.into_alum());
        assert_eq!(rewrapped.face_count(), faces);
    }
}
