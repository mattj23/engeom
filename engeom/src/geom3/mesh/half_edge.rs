//! This module provides a half-edge mesh structure, built on the `alum` library.
//!
//! The main type is [`HalfEdgeMesh`], which owns an `alum` polygon mesh rather than being an alias
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
//! for how tooling around it is going to work in practice, so [`HalfEdgeMesh::as_alum`],
//! [`HalfEdgeMesh::as_alum_mut`], and [`HalfEdgeMesh::into_alum`] hand out the underlying mesh for
//! anything this wrapper does not cover, including the whole of `alum`'s decimation and subdivision
//! machinery.
//!
//! [`NaAdaptor`] is the glue that lets `alum` work in terms of `nalgebra` types.

mod repair;
mod smoothing;

use crate::{Mesh3, Point3, Result, Vector3};
use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use repair::{RepairOpts, RepairReport, Repaired, repair_buffers};
use std::error::Error;

/// The `alum` mesh type that [`HalfEdgeMesh`] wraps.
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
/// Build one from a `Mesh3` with [`HalfEdgeMesh::from_mesh_repaired`] for measured data, or with
/// `TryFrom` when the input is known to be clean. Neither direction carries attributes, `is_solid`,
/// or the UV mapping, because this structure has nowhere to put them.
pub struct HalfEdgeMesh {
    inner: AlumMesh,
}

impl Default for HalfEdgeMesh {
    fn default() -> Self {
        Self::new()
    }
}

impl HalfEdgeMesh {
    /// Create an empty half-edge mesh.
    pub fn new() -> Self {
        Self {
            inner: AlumMesh::new(),
        }
    }

    /// Borrow the underlying `alum` mesh.
    ///
    /// This is the escape hatch for anything the wrapper does not expose, `alum`'s decimation and
    /// subdivision traits in particular. Using it couples your code to `alum`'s API rather than to
    /// this crate's.
    pub fn as_alum(&self) -> &AlumMesh {
        &self.inner
    }

    /// Mutably borrow the underlying `alum` mesh. See [`HalfEdgeMesh::as_alum`].
    pub fn as_alum_mut(&mut self) -> &mut AlumMesh {
        &mut self.inner
    }

    /// Consume the wrapper and return the underlying `alum` mesh.
    pub fn into_alum(self) -> AlumMesh {
        self.inner
    }

    /// Wrap an `alum` mesh which was built elsewhere.
    pub fn from_alum(inner: AlumMesh) -> Self {
        Self { inner }
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

impl TryFrom<&HalfEdgeMesh> for Mesh3 {
    type Error = Box<dyn Error>;

    fn try_from(value: &HalfEdgeMesh) -> std::result::Result<Self, Self::Error> {
        let vertices = value.clone_vertices()?;
        let faces = value.clone_faces()?;
        Ok(Mesh3::new(vertices, faces, false))
    }
}

impl TryFrom<&Mesh3> for HalfEdgeMesh {
    type Error = Box<dyn Error>;

    /// Build a half-edge mesh from a `Mesh3`, failing on the first element the structure will not
    /// accept.
    ///
    /// This is the strict path, appropriate when the input is supposed to be clean. Use
    /// [`HalfEdgeMesh::from_mesh_repaired`] on measured data, which is routinely non-manifold.
    ///
    /// Note that this is a test of what can be inserted, not of manifoldness: `alum` refuses an
    /// edge shared by three faces and refuses a fold, but accepts a bowtie whose fans are both
    /// open.
    fn try_from(value: &Mesh3) -> std::result::Result<Self, Self::Error> {
        let mut result = HalfEdgeMesh::new();
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
        let he = HalfEdgeMesh::new();
        assert_eq!(he.vertex_count(), 0);
        assert_eq!(he.face_count(), 0);
        assert_eq!(he.edge_count(), 0);
        assert!(he.clone_vertices().unwrap().is_empty());
        assert!(he.clone_faces().unwrap().is_empty());
    }

    #[test]
    fn a_clean_mesh_round_trips() {
        let mesh = crate::tests::engine_blade();
        let he = HalfEdgeMesh::try_from(&mesh).unwrap();

        assert_eq!(he.vertex_count(), mesh.points().len());
        assert_eq!(he.face_count(), mesh.faces().len());

        let back = Mesh3::try_from(&he).unwrap();
        assert_eq!(back.points().len(), mesh.points().len());
        assert_eq!(back.faces().len(), mesh.faces().len());
    }

    /// The escape hatch has to actually reach `alum`'s own traits, since that is the whole reason
    /// it exists. This is the path the decimater will take.
    #[test]
    fn the_escape_hatch_reaches_alum_decimation() {
        use alum::decimate::HasDecimation;
        use alum::decimate::quadric::{QuadricDecimater, QuadricType};

        let mesh = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh::try_from(&mesh).unwrap();
        let before = he.face_count();

        let mut decimater = QuadricDecimater::new(he.as_alum(), QuadricType::Triangle).unwrap();
        let collapses = he
            .as_alum_mut()
            .decimate_to_face_count(&mut decimater, before / 2)
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
        let he = HalfEdgeMesh::from_mesh_repaired(&mesh, &RepairOpts::default())
            .unwrap()
            .0;
        let faces = he.face_count();

        let rewrapped = HalfEdgeMesh::from_alum(he.into_alum());
        assert_eq!(rewrapped.face_count(), faces);
    }
}
