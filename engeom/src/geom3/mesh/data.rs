//! This module contains `MeshData3`, a plain container for triangle mesh data and its associated
//! per-element attributes.  Mesh information is stored in the common point (vertex) and face
//! (triangle) buffer representation.
//!
//! Unlike `Mesh3`, which is built around a `parry3d` `TriMesh` and pays for a BVH build on every
//! construction, `MeshData3` holds nothing but the point and face buffers and whatever attributes
//! are attached to them. It is the type to reach for when actively building or modifying mesh data
//! in place, or when working directly with serialization/deserialization. Convert it into the
//! accelerated `Mesh3` or the half-edge representation to work with it in those contexts.
//!
//! There are convenience construction methods on `Mesh3` which transparently pass-through to this
//! struct for ergonomic reasons.

mod editing;
mod filtering;
mod operations;
mod primitives;
mod subsets;

pub use crate::geom3::attributes3::Attr3;

use crate::geom3::Mesh3;
use crate::geom3::attributes3::{FaceAttrSet3, PointAttrSet3};
use crate::geom3::mesh::MeshView3;
use crate::geom3::mesh::algorithms;
use crate::geom3::mesh::algorithms::{
    OffsetOpts, compute_face_offset_points, compute_normal_displaced_points,
};
use crate::io::{load_g3d_mesh_data, read_tc_mesh_file, write_tc_mesh_file};
use crate::{Point3, Result, UnitVec3};
use std::fmt;
use std::path::Path;

#[cfg(feature = "ply")]
use crate::io::{PlyWriteOpts, load_ply_mesh_data, write_ply_mesh};
#[cfg(feature = "stl")]
use crate::io::{StlWriteOpts, load_stl_mesh_data, write_stl_mesh};

/// A container for the raw data of a triangle mesh: a buffer of points, a buffer of faces indexing
/// into it, and the per-element attributes attached to either domain.
///
/// This type performs no spatial acceleration of any kind. It is cheap to construct, cheap to edit,
/// and is the form in which mesh data is read from and written to files. Spatial queries require
/// converting it into an accelerated representation, which is where the cost of building a bounding
/// volume hierarchy is paid.
///
/// # Relationship with Other Representations
///
/// The `MeshData3` type serves as the serialization gateway for mesh data and should be the
/// default choice for the following cases:
///
/// - You don't need to do any spatial queries on the mesh at all
/// - You have to do a bunch of edits to the contents of the mesh buffers and don't want to
///   regenerate an entire acceleration structure every operation
/// - You are doing something custom with serialization or deserialization
///
/// If you need to perform spatial queries or do the types of edits requiring a half-edge
/// representation, you should use the `Mesh3` and `HalfEdgeMesh3` types, respectively. The
/// `MeshData3` type has consuming `TryFrom<T>` and `TryInto<T>` implementations for these types.
///
/// # Invariants
///
/// A `MeshData3` is never allowed to exist in an incoherent state. Every face index refers to a
/// point which exists, and every attribute array has a length matching the element count of the
/// domain it belongs to. These are checked on construction and maintained by every method which
/// modifies the mesh, which is why the point and face buffers are private and why attributes are
/// set through methods that supply the counts on the caller's behalf.
///
/// An empty mesh is legal.
#[derive(Clone)]
pub struct MeshData3 {
    points: Vec<Point3>,
    faces: Vec<[u32; 3]>,
    point_attrs: PointAttrSet3,
    face_attrs: FaceAttrSet3,
}

// ===============================================================================================
// Construction
// ===============================================================================================

impl MeshData3 {
    /// Create a new mesh from a buffer of points and a buffer of faces indexing into it, with no
    /// attributes attached.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    /// * `faces`: triangles given as triples of indices into `points`
    ///
    /// returns: `Result<MeshData3>`, failing if any face refers to a point which does not exist
    pub fn new(points: Vec<Point3>, faces: Vec<[u32; 3]>) -> Result<Self> {
        Self::new_with_attrs(points, faces, PointAttrSet3::empty(), FaceAttrSet3::empty())
    }

    /// Create a new mesh from a buffer of points, a buffer of faces indexing into it, and the
    /// per-element attributes of both domains.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    /// * `faces`: triangles given as triples of indices into `points`
    /// * `point_attrs`: the per-point attributes to attach, whose arrays must match the point count
    /// * `face_attrs`: the per-face attributes to attach, whose arrays must match the face count
    ///
    /// returns: `Result<MeshData3>`, failing if any face refers to a point which does not exist or
    /// if any attribute array is the wrong length
    pub fn new_with_attrs(
        points: Vec<Point3>,
        faces: Vec<[u32; 3]>,
        point_attrs: PointAttrSet3,
        face_attrs: FaceAttrSet3,
    ) -> Result<Self> {
        check_face_indices(&faces, points.len())?;
        point_attrs.validate(points.len())?;
        face_attrs.validate(faces.len())?;

        Ok(Self {
            points,
            faces,
            point_attrs,
            face_attrs,
        })
    }

    /// Create an empty mesh, with no points, no faces, and no attributes.
    pub fn empty() -> Self {
        Self {
            points: Vec::new(),
            faces: Vec::new(),
            point_attrs: PointAttrSet3::empty(),
            face_attrs: FaceAttrSet3::empty(),
        }
    }

    /// Load a triangle mesh from a PLY file, preserving every property the file carries.
    ///
    /// The file must have a `vertex` element with `x`, `y`, and `z` properties.
    ///
    /// A `face` element is optional but will result in mesh data with points and no faces, which
    /// is OK for the `MeshData3` type but won't be for `Mesh3`. If you attempted to open a PLY
    /// file that was a save of a point cloud, this is the result you'll get, and you will need
    /// to verify the number of faces independently.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the PLY file
    ///
    /// returns: `Result<MeshData3>`
    #[cfg(feature = "ply")]
    pub fn load_ply(path: &Path) -> Result<Self> {
        load_ply_mesh_data(path)
    }

    /// Write this mesh data to a PLY file, preserving every attribute it carries.
    ///
    /// If you use the default options, the data will be saved in binary format using double
    /// floating point precision.  If you wish to alter this, construct the `PlyWriteOpts`
    /// directly.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `opts`: encoding and header options
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    #[cfg(feature = "ply")]
    pub fn save_ply(&self, path: &Path, opts: &PlyWriteOpts) -> Result<()> {
        write_ply_mesh(path, self, opts)
    }

    /// Load a triangle mesh from an STL file, in either the ascii or binary encoding.
    ///
    /// STL is a triangle soup with no point identity, so the points are recovered by welding on
    /// exact coordinate equality. See `load_stl_mesh_data` for what that does and does not
    /// recover, and for the precision the format costs you.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the STL file
    ///
    /// returns: `Result<MeshData3>`
    #[cfg(feature = "stl")]
    pub fn load_stl(path: &Path) -> Result<Self> {
        load_stl_mesh_data(path)
    }

    /// Write this mesh to an STL file, which carries geometry and nothing else.
    ///
    /// The default options write binary with no attribute loss permitted, so a mesh carrying any
    /// attributes at all is refused rather than silently stripped. Set `allow_attribute_loss` to
    /// accept the loss.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `opts`: encoding, header, and attribute loss options
    ///
    /// returns: `Result<()>`
    #[cfg(feature = "stl")]
    pub fn save_stl(&self, path: &Path, opts: &StlWriteOpts) -> Result<()> {
        write_stl_mesh(path, self, opts)
    }

    /// Load a triangle mesh from a GOM `.g3d` file, the format written by GOM's Atos scanner and
    /// GOM/Zeiss Inspect.
    ///
    /// A `.g3d` file can in principle carry several kinds of view (point clouds, scan sections,
    /// colored meshes), but only the triangle mesh view is read here; any other view present in
    /// the file is skipped. See [`load_g3d_mesh_data`] for the full format support notes.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the `.g3d` file
    ///
    /// returns: `Result<MeshData3>`
    pub fn load_g3d(path: &Path) -> Result<Self> {
        load_g3d_mesh_data(path)
    }

    /// Load a triangle mesh from a tolerance-compressed `.tcmesh` file.
    ///
    /// The recovered positions are guaranteed to be within the tolerance which was given at write
    /// time, and the connectivity is exact. The vertices are not in the order they were in before
    /// the file was written; see [`MeshData3::save_tcmesh`].
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the `.tcmesh` file
    ///
    /// returns: `Result<MeshData3>`
    pub fn load_tcmesh(path: &Path) -> Result<Self> {
        read_tc_mesh_file(path)
    }

    /// Write this mesh to a tolerance-compressed `.tcmesh` file, which carries geometry and
    /// nothing else.
    ///
    /// Vertex positions are quantized to the narrowest bit width per axis which keeps every one of
    /// them within `tol` of where it started, while the connectivity is stored exactly. A smaller
    /// tolerance costs more bytes per vertex.
    ///
    /// Writing **renumbers the vertices**, because reordering them is where most of the format's
    /// advantage comes from. A mesh read back describes the same surface but not with the same
    /// indices, so per-vertex data kept outside the file cannot assume it still lines up. See the
    /// [`crate::io::tol_compress::mesh`] module for how to compute the same ordering up front.
    ///
    /// Unlike the other geometry-only formats here there is no attribute loss option: a mesh
    /// carrying any attribute at all is refused, with the error naming what would have been lost.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `tol`: the largest acceptable round-trip position error for any vertex, in the same units
    ///   as the coordinates
    ///
    /// returns: `Result<()>`, failing if the mesh carries any attributes
    pub fn save_tcmesh(&self, path: &Path, tol: f64) -> Result<()> {
        write_tc_mesh_file(path, self, tol)
    }
}

// ===============================================================================================
// Serialization support
// ===============================================================================================

impl MeshData3 {
    /// Verify that the caller has accepted the loss of this mesh's attributes when using a format
    /// that cannot represent them. This delegates the check to
    /// [`MeshView3::check_attribute_loss`].
    ///
    /// # Arguments
    ///
    /// * `format`: the name of the target format, used in the error message
    /// * `allow_loss`: whether the caller has accepted the loss, which comes from the
    ///   `allow_attribute_loss` field of the format's options struct
    ///
    /// returns: `Result<()>`
    pub fn check_attribute_loss(&self, format: &str, allow_loss: bool) -> Result<()> {
        MeshView3::from(self).check_attribute_loss(format, allow_loss)
    }

    /// Borrow this mesh as the read-only view accepted by the file writers.
    pub fn view(&self) -> MeshView3<'_> {
        MeshView3::from(self)
    }
}

// ===============================================================================================
// Core access
// ===============================================================================================

impl MeshData3 {
    /// Get a reference to the points of the mesh.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Get a reference to the face indices of the mesh.
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// Get the number of points in the mesh.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Get the number of faces in the mesh.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Returns true if the mesh has no points and no faces.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty() && self.faces.is_empty()
    }

    /// Get a reference to the per-point attributes attached to this mesh.
    pub fn point_attrs(&self) -> &PointAttrSet3 {
        &self.point_attrs
    }

    /// Get a reference to the per-face attributes attached to this mesh.
    pub fn face_attrs(&self) -> &FaceAttrSet3 {
        &self.face_attrs
    }

    /// Returns true if the mesh carries any attributes in either domain.
    pub fn has_attrs(&self) -> bool {
        !self.point_attrs.is_empty() || !self.face_attrs.is_empty()
    }
}

// ===============================================================================================
// Derived computations
// ===============================================================================================

impl MeshData3 {
    /// Compute a normal for every point by averaging the normals of the faces which touch it,
    /// weighting each by the interior angle of that face at that point.
    ///
    /// This always computes from the faces. It does not consult any normals already stored on the
    /// mesh, so that a caller which wants measured normals in preference to computed ones can make
    /// that choice explicitly rather than having it made silently here.
    ///
    /// See [`algorithms::normals`] for why the faces are weighted by angle rather than by area,
    /// and for what this normal is and is not suitable for.
    ///
    /// returns: `Result<Vec<UnitVec3>>`, failing if any point has no well-defined normal
    pub fn compute_point_normals(&self) -> Result<Vec<UnitVec3>> {
        algorithms::compute_point_normals(&self.points, &self.faces)
    }

    /// Compute the unit normal of every face, from the winding of its three points.
    ///
    /// returns: `Result<Vec<UnitVec3>>`, one normal per face, failing on a face which has no
    /// well-defined normal because its points are coincident or collinear
    pub fn compute_face_normals(&self) -> Result<Vec<UnitVec3>> {
        algorithms::compute_face_normals(&self.points, &self.faces)
    }

    /// Compute the area of every face.
    ///
    /// returns: `Result<Vec<f64>>`, one area per face. A degenerate face has an area of zero rather
    /// than being an error.
    pub fn compute_face_areas(&self) -> Result<Vec<f64>> {
        algorithms::compute_face_areas(&self.points, &self.faces)
    }

    /// Compute the centroid of every face, which is the mean of its three points.
    ///
    /// returns: `Result<Vec<Point3>>`, one centroid per face
    pub fn compute_face_centers(&self) -> Result<Vec<Point3>> {
        algorithms::compute_face_centers(&self.points, &self.faces)
    }

    /// Offset this mesh's surface by `distance` in place, leaving its faces and attributes alone.
    ///
    /// Only the point positions change, so the topology and every per-element attribute survive
    /// unchanged. Note that any stored `point_normals` are **not** recomputed: an offset can change
    /// them where the surface is curved, so a caller relying on them should recompute them after.
    ///
    /// # Arguments
    ///
    /// * `distance`: how far to move the surface, in the mesh's length units
    /// * `opts`: how to handle unconstrained directions and near-degenerate spikes
    ///
    /// returns: `Result<()>`, leaving the mesh untouched on failure
    pub fn offset_faces_in_place(&mut self, distance: f64, opts: &OffsetOpts) -> Result<()> {
        self.points = compute_face_offset_points(&self.points, &self.faces, distance, opts)?;
        Ok(())
    }

    pub fn offset_faces_copy(&self, distance: f64, opts: &OffsetOpts) -> Result<Self> {
        let mut clone = self.clone();
        clone.offset_faces_in_place(distance, opts)?;
        Ok(clone)
    }

    pub fn offset_points_copy(&self, distance: f64) -> Result<Self> {
        let mut clone = self.clone();
        clone.offset_points_in_place(distance)?;
        Ok(clone)
    }

    /// Offset this mesh's points by `distance` along their surface normals, editing the mesh in
    /// place. The faces and the attributes are left alone, so the topology and per-element
    /// attributes are unchanged.
    ///
    /// If the `MeshData3` does not have point normals stored in its attributes, temporary normals
    /// will be computed for the operation and discarded afterwards.
    ///
    /// # Arguments
    ///
    /// * `distance`: the distance to move, in the mesh's length units
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    pub fn offset_points_in_place(&mut self, distance: f64) -> Result<()> {
        if let Some(normals) = self.point_attrs.normals() {
            self.points = compute_normal_displaced_points(&self.points, normals, distance)?;
        } else {
            let local_normals = self.compute_point_normals()?;
            self.points = compute_normal_displaced_points(&self.points, &local_normals, distance)?;
        }

        Ok(())
    }
}

// ===============================================================================================
// Attribute access
// ===============================================================================================

impl MeshData3 {
    /// Get the per-point unit normals, if present.
    pub fn point_normals(&self) -> Option<&[UnitVec3]> {
        self.point_attrs.normals()
    }

    /// Get the per-point RGB colors, if present.
    pub fn point_colors(&self) -> Option<&[[u8; 3]]> {
        self.point_attrs.colors()
    }

    /// Get the per-point standard deviations, if present. These are 1-sigma values in the mesh's
    /// own length units.
    pub fn point_stdev(&self) -> Option<&[f64]> {
        self.point_attrs.stdev()
    }

    /// Get the per-face RGB colors, if present.
    pub fn face_colors(&self) -> Option<&[[u8; 3]]> {
        self.face_attrs.colors()
    }

    /// Get the per-face labels, if present.
    pub fn face_labels(&self) -> Option<&[u32]> {
        self.face_attrs.labels()
    }

    /// Get the open-map per-point attribute stored under the given name, if present.
    pub fn point_attr(&self, name: &str) -> Option<&Attr3> {
        self.point_attrs.attr(name)
    }

    /// Get the open-map per-face attribute stored under the given name, if present.
    pub fn face_attr(&self, name: &str) -> Option<&Attr3> {
        self.face_attrs.attr(name)
    }
}

// ===============================================================================================
// Attribute mutation
// ===============================================================================================

impl MeshData3 {
    /// Set or clear the per-point unit normals.
    ///
    /// # Arguments
    ///
    /// * `values`: the normals to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_normals(&mut self, values: Option<Vec<UnitVec3>>) -> Result<()> {
        self.point_attrs.set_normals(values, self.points.len())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.point_attrs.set_colors(values, self.points.len())
    }

    /// Set or clear the per-point standard deviations, which must be 1-sigma values in the mesh's
    /// own length units.
    ///
    /// # Arguments
    ///
    /// * `values`: the standard deviations to store, or `None` to clear them. Must match the point
    ///   count, and must be finite and non-negative.
    ///
    /// returns: `Result<()>`
    pub fn set_point_stdev(&mut self, values: Option<Vec<f64>>) -> Result<()> {
        self.point_attrs.set_stdev(values, self.points.len())
    }

    /// Set or clear the per-face RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.face_attrs.set_colors(values, self.faces.len())
    }

    /// Set or clear the per-face labels.
    ///
    /// # Arguments
    ///
    /// * `values`: the labels to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_labels(&mut self, values: Option<Vec<u32>>) -> Result<()> {
        self.face_attrs.set_labels(values, self.faces.len())
    }

    /// Insert an open-map per-point attribute under the given name, replacing any attribute already
    /// stored there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be a reserved name
    /// * `attr`: the attribute array to store, whose length must match the point count
    ///
    /// returns: `Result<()>`
    pub fn insert_point_attr(&mut self, name: &str, attr: Attr3) -> Result<()> {
        self.point_attrs.insert_attr(name, attr, self.points.len())
    }

    /// Insert an open-map per-face attribute under the given name, replacing any attribute already
    /// stored there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be a reserved name
    /// * `attr`: the attribute array to store, whose length must match the face count
    ///
    /// returns: `Result<()>`
    pub fn insert_face_attr(&mut self, name: &str, attr: Attr3) -> Result<()> {
        self.face_attrs.insert_attr(name, attr, self.faces.len())
    }

    /// Remove and return the open-map per-point attribute stored under the given name.
    pub fn remove_point_attr(&mut self, name: &str) -> Option<Attr3> {
        self.point_attrs.remove_attr(name)
    }

    /// Remove and return the open-map per-face attribute stored under the given name.
    pub fn remove_face_attr(&mut self, name: &str) -> Option<Attr3> {
        self.face_attrs.remove_attr(name)
    }

    /// Replace the entire set of per-point attributes without changing the face domain.
    ///
    /// This allows a point-domain set built elsewhere—such as by a reader that handles the two
    /// domains separately or while converting a point cloud into a mesh—to be attached without
    /// calling each typed setter individually.
    ///
    /// # Arguments
    ///
    /// * `attrs`: the per-point attributes to attach, whose arrays must match the point count
    ///
    /// returns: `Result<()>`, leaving the existing attributes untouched on failure
    pub fn set_point_attrs(&mut self, attrs: PointAttrSet3) -> Result<()> {
        attrs.validate(self.points.len())?;
        self.point_attrs = attrs;
        Ok(())
    }

    /// Replace the entire set of per-face attributes without changing the point domain.
    ///
    /// # Arguments
    ///
    /// * `attrs`: the per-face attributes to attach, whose arrays must match the face count
    ///
    /// returns: `Result<()>`, leaving the existing attributes untouched on failure
    pub fn set_face_attrs(&mut self, attrs: FaceAttrSet3) -> Result<()> {
        attrs.validate(self.faces.len())?;
        self.face_attrs = attrs;
        Ok(())
    }

    /// Remove and return all per-point attributes, leaving that domain empty on the mesh.
    pub fn take_point_attrs(&mut self) -> PointAttrSet3 {
        std::mem::take(&mut self.point_attrs)
    }

    /// Remove and return all per-face attributes, leaving that domain empty on the mesh.
    pub fn take_face_attrs(&mut self) -> FaceAttrSet3 {
        std::mem::take(&mut self.face_attrs)
    }

    /// Consume the mesh and return ownership of its four components: the point buffer, the face
    /// buffer, the per-point attributes, and the per-face attributes.
    ///
    /// This is the counterpart to `new_with_attrs` and allows the data to be handed to another
    /// representation without copying the buffers. Once the mesh is decomposed, nothing enforces
    /// the invariants between its components, so callers are responsible for preserving consistency
    /// if they reassemble them.
    ///
    /// returns: `(Vec<Point3>, Vec<[u32; 3]>, PointAttrSet3, FaceAttrSet3)`
    pub fn into_parts(self) -> (Vec<Point3>, Vec<[u32; 3]>, PointAttrSet3, FaceAttrSet3) {
        (self.points, self.faces, self.point_attrs, self.face_attrs)
    }

    /// Consume the mesh data and return an accelerated `Mesh3` built over the same buffers.
    ///
    /// This is the consuming conversion into the accelerated type: the point and face buffers are
    /// moved rather than copied, while the bounding volume hierarchy must still be built. The
    /// opposite conversion is `Mesh3::into_data`, which cannot avoid a copy because `parry3d`
    /// provides no way to move the buffers back out of a `TriMesh`.
    ///
    /// # Arguments
    ///
    /// * `is_solid`: whether distance queries should treat points inside the mesh as having zero
    ///   distance. The plain container has no such concept, so the caller must supply it.
    ///
    /// returns: `Result<Mesh3>`, failing if the mesh has no faces, since there is nothing to build
    /// an acceleration structure over
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{MeshData3, Point3};
    ///
    /// let data = MeshData3::create_box(2.0, 2.0, 2.0);
    /// let counts = (data.point_count(), data.face_count());
    ///
    /// let mesh = data.into_mesh(true).unwrap();
    /// assert!(mesh.is_solid());
    /// assert_eq!((mesh.point_count(), mesh.face_count()), counts);
    /// ```
    pub fn into_mesh(self, is_solid: bool) -> Result<Mesh3> {
        Mesh3::from_data(self, is_solid)
    }
}

// ===============================================================================================
// Helpers
// ===============================================================================================

impl fmt::Debug for MeshData3 {
    /// Summarize the mesh rather than dumping its buffers, which are routinely large enough to make
    /// a derived `Debug` implementation useless.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MeshData3")
            .field("points", &self.points.len())
            .field("faces", &self.faces.len())
            .field("point_attrs", &self.point_attrs)
            .field("face_attrs", &self.face_attrs)
            .finish()
    }
}

/// Verify that every face refers to a point which exists.
pub(crate) fn check_face_indices(faces: &[[u32; 3]], n_points: usize) -> Result<()> {
    for (i, face) in faces.iter().enumerate() {
        for index in face {
            if *index as usize >= n_points {
                return Err(format!(
                    "Face {i} refers to point {index}, but the mesh has only {n_points} points"
                )
                .into());
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use approx::assert_relative_eq;

    fn unit_square() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3]];
        (points, faces)
    }

    fn square_mesh() -> MeshData3 {
        let (points, faces) = unit_square();
        MeshData3::new(points, faces).unwrap()
    }

    /// The container method reaches the same writer and reader the `io` functions do. The format's
    /// own round-trip guarantees are covered in `io::tol_compress::mesh`, so what is checked here
    /// is the mapping: what went out came back, describing the same surface.
    #[test]
    fn a_tcmesh_round_trips_through_the_container_methods() -> Result<()> {
        let mesh = square_mesh();
        let tol = 1e-5;
        let path = std::env::temp_dir().join("engeom_mesh_data_tcmesh_round_trip.tcmesh");

        mesh.save_tcmesh(&path, tol)?;
        let recovered = MeshData3::load_tcmesh(&path)?;

        assert_eq!(recovered.point_count(), mesh.point_count());
        assert_eq!(recovered.face_count(), mesh.face_count());

        // Writing renumbers the vertices, so the assertion is over the set of positions rather than
        // over the buffer order.
        for point in recovered.points() {
            assert!(
                mesh.points().iter().any(|p| (p - point).norm() <= tol),
                "{point:?} is not within {tol} of any original point"
            );
        }

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// The format stores geometry only and refuses anything else outright, and that refusal has to
    /// survive the trip through the container method rather than being swallowed by it.
    #[test]
    fn saving_a_tcmesh_refuses_a_mesh_carrying_attributes() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_point_stdev(Some(vec![0.0, 0.1, 0.2, 0.3]))?;

        let path = std::env::temp_dir().join("engeom_mesh_data_tcmesh_refused.tcmesh");
        let err = mesh.save_tcmesh(&path, 1e-5).unwrap_err().to_string();

        assert!(
            err.contains("stdev"),
            "the error should name what would be lost: {err}"
        );

        Ok(())
    }

    #[test]
    fn new_accepts_a_valid_mesh() {
        let mesh = square_mesh();

        assert_eq!(mesh.point_count(), 4);
        assert_eq!(mesh.face_count(), 2);
        assert!(!mesh.is_empty());
        assert!(!mesh.has_attrs());
    }

    #[test]
    fn new_rejects_an_out_of_range_face_index() {
        let (points, _) = unit_square();
        assert!(MeshData3::new(points, vec![[0, 1, 4]]).is_err());
    }

    #[test]
    fn new_rejects_any_face_against_an_empty_point_buffer() {
        assert!(MeshData3::new(Vec::new(), vec![[0, 0, 0]]).is_err());
    }

    #[test]
    fn an_empty_mesh_is_legal() -> Result<()> {
        let mesh = MeshData3::new(Vec::new(), Vec::new())?;
        assert!(mesh.is_empty());
        assert_eq!(mesh.point_count(), 0);
        assert_eq!(mesh.face_count(), 0);

        assert!(MeshData3::empty().is_empty());

        Ok(())
    }

    #[test]
    fn points_without_faces_are_legal() -> Result<()> {
        let (points, _) = unit_square();
        let mesh = MeshData3::new(points, Vec::new())?;

        assert_eq!(mesh.point_count(), 4);
        assert_eq!(mesh.face_count(), 0);
        assert!(!mesh.is_empty());

        Ok(())
    }

    #[test]
    fn new_with_attrs_rejects_a_mismatched_attribute() {
        let (points, faces) = unit_square();
        let mut attrs = PointAttrSet3::empty();
        attrs.set_colors(Some(vec![[0, 0, 0]]), 1).unwrap();

        // The attribute set is internally consistent for a one-point mesh, but not for this one.
        assert!(MeshData3::new_with_attrs(points, faces, attrs, FaceAttrSet3::empty()).is_err());

        // Check the face domain in the same way.
        let (points, faces) = unit_square();
        let mut attrs = FaceAttrSet3::empty();
        attrs.set_labels(Some(vec![1, 2, 3]), 3).unwrap();
        assert!(MeshData3::new_with_attrs(points, faces, PointAttrSet3::empty(), attrs).is_err());
    }

    #[test]
    fn attribute_setters_supply_the_counts() -> Result<()> {
        let mut mesh = square_mesh();

        mesh.set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]))?;
        mesh.set_face_labels(Some(vec![7, 8]))?;
        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![0.5; 4]))?;
        mesh.insert_face_attr("material_index", Attr3::Label(vec![1, 2]))?;

        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);
        assert_eq!(mesh.face_labels().unwrap(), &[7, 8]);
        assert_eq!(mesh.point_attr("confidence").unwrap().len(), 4);
        assert_eq!(mesh.face_attr("material_index").unwrap().len(), 2);

        // The wrong length is rejected without the caller having to know the count.
        assert!(mesh.set_point_stdev(Some(vec![0.1, 0.2])).is_err());
        assert!(mesh.set_face_labels(Some(vec![1, 2, 3])).is_err());

        Ok(())
    }

    #[test]
    fn set_attrs_validates_against_the_current_counts() -> Result<()> {
        let mut mesh = square_mesh();

        let mut good = PointAttrSet3::empty();
        good.set_normals(
            Some(vec![UnitVec3::new_normalize(Vector3::z()); 4]),
            mesh.point_count(),
        )?;
        mesh.set_point_attrs(good)?;
        assert!(mesh.point_normals().is_some());

        let mut bad = FaceAttrSet3::empty();
        bad.set_labels(Some(vec![1, 2, 3]), 3)?;
        assert!(mesh.set_face_attrs(bad).is_err());

        let mut bad = PointAttrSet3::empty();
        bad.set_stdev(Some(vec![0.1]), 1)?;
        assert!(mesh.set_point_attrs(bad).is_err());

        // The rejected sets must have left the existing attributes in place.
        assert!(mesh.point_normals().is_some());

        Ok(())
    }

    /// A mesh carrying only a face attribute is not bare, which ensures that checking only the
    /// point domain cannot incorrectly classify it as bare.
    #[test]
    fn a_face_only_attribute_counts_as_having_attrs() -> Result<()> {
        let mut mesh = square_mesh();
        assert!(!mesh.has_attrs());

        mesh.set_face_labels(Some(vec![7, 8]))?;
        assert!(mesh.has_attrs());
        assert!(mesh.point_attrs().is_empty());
        assert!(mesh.check_attribute_loss("a test", false).is_err());

        Ok(())
    }

    #[test]
    fn take_attrs_leaves_the_mesh_bare() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_face_labels(Some(vec![7, 8]))?;

        let taken = mesh.take_face_attrs();

        assert_eq!(taken.labels().unwrap(), &[7, 8]);
        assert!(!mesh.has_attrs());
        assert!(mesh.face_labels().is_none());

        Ok(())
    }

    #[test]
    fn debug_summarizes_instead_of_dumping_buffers() {
        let text = format!("{:?}", square_mesh());
        assert!(text.contains("points: 4"), "{text}");
        assert!(text.contains("faces: 2"), "{text}");
    }

    /// The offset moves the points and leaves everything else alone, which is what makes it safe
    /// to apply to a mesh carrying measured data.
    #[test]
    fn offsetting_moves_the_points_and_keeps_the_attributes() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]))?;
        mesh.set_face_labels(Some(vec![7, 9]))?;

        let d = 2.0;
        let moved = mesh.offset_faces_copy(d, &OffsetOpts::default())?;

        // A flat square offsets straight along its normal. Compared approximately, since the
        // answer comes out of an SVD rather than a closed form.
        for (before, after) in mesh.points().iter().zip(moved.points().iter()) {
            assert_relative_eq!(
                after.coords,
                before.coords + Vector3::z() * d,
                epsilon = 1.0e-12
            );
        }

        assert_eq!(moved.faces(), mesh.faces());
        assert_eq!(moved.point_stdev(), Some([0.1, 0.2, 0.3, 0.4].as_slice()));
        assert_eq!(moved.face_labels(), Some([7, 9].as_slice()));

        Ok(())
    }

    /// The two consuming conversions are inverses: the buffers and attributes moved into `Mesh3`
    /// are the same ones returned from it.
    #[test]
    fn the_consuming_conversions_round_trip() -> Result<()> {
        let mut data = square_mesh();
        let computed = data.compute_point_normals()?;
        data.set_point_normals(Some(computed))?;

        let points = data.points().to_vec();
        let faces = data.faces().to_vec();
        let normals = data.point_normals().map(|n| n.to_vec());
        assert!(normals.is_some(), "the fixture must carry an attribute");

        let mesh = data.into_mesh(true)?;
        assert!(mesh.is_solid());

        let back = mesh.into_data();
        assert_eq!(back.points(), points.as_slice());
        assert_eq!(back.faces(), faces.as_slice());
        assert_eq!(back.point_normals().map(|n| n.to_vec()), normals);

        Ok(())
    }

    /// The plain container may hold points without faces, but an acceleration structure has nothing
    /// to build over in that case, so the conversion must report an error.
    #[test]
    fn a_mesh_without_faces_cannot_become_an_accelerated_mesh() {
        let data = MeshData3::new(vec![Point3::origin()], vec![]).unwrap();
        assert!(data.into_mesh(false).is_err());
    }

    /// A failed offset must leave the mesh exactly as it was, not half moved.
    #[test]
    fn a_failed_offset_leaves_the_mesh_untouched() {
        // An orphan point has nothing to offset against, so the whole operation fails.
        let (mut points, faces) = unit_square();
        points.push(Point3::new(9.0, 9.0, 9.0));
        let mut mesh = MeshData3::new(points, faces).unwrap();
        let before = mesh.points().to_vec();

        assert!(
            mesh.offset_faces_in_place(1.0, &OffsetOpts::default())
                .is_err()
        );
        assert_eq!(mesh.points(), before.as_slice());
    }
}
