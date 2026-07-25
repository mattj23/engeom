//! This module contains `MeshData3`, a plain container for triangle mesh data and its associated
//! per-element attributes.  Mesh information is stored in the common point (vertex) and face
//! (triangle) buffer representation.
//!
//! Unlike `Mesh`, which is built around a `parry3d` `TriMesh` and pays for a BVH build on every
//! construction, `MeshData3` holds nothing but the point and face buffers and whatever attributes
//! are attached to them. It is the type to reach for when actively building or modifying mesh data
//! in place, or when working directly with serialization/deserialization. Convert it into the
//! accelerated `Mesh` or the half-edge representation to work with it in those contexts.
//!
//! There are convenience construction methods on `Mesh` which transparently pass-through to this
//! struct for ergonomic reasons.

mod attribute_set;
mod attributes;
mod editing;
mod operations;
mod subsets;

pub use attribute_set::MeshAttrSet3;
pub use attributes::MeshAttr3;

use crate::geom3::mesh::algorithms;
use crate::{Point3, Result, UnitVec3};
use std::fmt;

#[cfg(feature = "ply")]
use crate::io::load_ply_mesh_data;
#[cfg(feature = "ply")]
use std::path::Path;

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
/// representation, you should use the `Mesh` and `HalfEdgeMesh` types, respectively. The
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
    attrs: MeshAttrSet3,
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
        Self::new_with_attrs(points, faces, MeshAttrSet3::empty())
    }

    /// Create a new mesh from a buffer of points, a buffer of faces indexing into it, and a set of
    /// per-element attributes.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    /// * `faces`: triangles given as triples of indices into `points`
    /// * `attrs`: the attributes to attach, whose arrays must match the point and face counts
    ///
    /// returns: `Result<MeshData3>`, failing if any face refers to a point which does not exist or
    /// if any attribute array is the wrong length
    pub fn new_with_attrs(
        points: Vec<Point3>,
        faces: Vec<[u32; 3]>,
        attrs: MeshAttrSet3,
    ) -> Result<Self> {
        check_face_indices(&faces, points.len())?;
        attrs.validate(points.len(), faces.len())?;

        Ok(Self {
            points,
            faces,
            attrs,
        })
    }

    /// Create an empty mesh, with no points, no faces, and no attributes.
    pub fn empty() -> Self {
        Self {
            points: Vec::new(),
            faces: Vec::new(),
            attrs: MeshAttrSet3::empty(),
        }
    }

    /// Load a triangle mesh from a PLY file, preserving every property the file carries.
    ///
    /// The file must have a `vertex` element with `x`, `y`, and `z` properties.
    ///
    /// A `face` element is optional but will result in mesh data with points and no faces, which
    /// is OK for the `MeshData3` type but won't be for `Mesh`. If you attempted to open a PLY
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

    /// Get a reference to the full set of per-element attributes attached to this mesh.
    pub fn attrs(&self) -> &MeshAttrSet3 {
        &self.attrs
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
}

// ===============================================================================================
// Attribute access
// ===============================================================================================

impl MeshData3 {
    /// Get the per-point unit normals, if present.
    pub fn point_normals(&self) -> Option<&[UnitVec3]> {
        self.attrs.point_normals()
    }

    /// Get the per-point RGB colors, if present.
    pub fn point_colors(&self) -> Option<&[[u8; 3]]> {
        self.attrs.point_colors()
    }

    /// Get the per-point standard deviations, if present. These are 1-sigma values in the mesh's
    /// own length units.
    pub fn point_stdev(&self) -> Option<&[f64]> {
        self.attrs.point_stdev()
    }

    /// Get the per-face RGB colors, if present.
    pub fn face_colors(&self) -> Option<&[[u8; 3]]> {
        self.attrs.face_colors()
    }

    /// Get the per-face labels, if present.
    pub fn face_labels(&self) -> Option<&[u32]> {
        self.attrs.face_labels()
    }

    /// Get the open-map per-point attribute stored under the given name, if present.
    pub fn point_attr(&self, name: &str) -> Option<&MeshAttr3> {
        self.attrs.point_attr(name)
    }

    /// Get the open-map per-face attribute stored under the given name, if present.
    pub fn face_attr(&self, name: &str) -> Option<&MeshAttr3> {
        self.attrs.face_attr(name)
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
        self.attrs.set_point_normals(values, self.points.len())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.attrs.set_point_colors(values, self.points.len())
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
        self.attrs.set_point_stdev(values, self.points.len())
    }

    /// Set or clear the per-face RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.attrs.set_face_colors(values, self.faces.len())
    }

    /// Set or clear the per-face labels.
    ///
    /// # Arguments
    ///
    /// * `values`: the labels to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_labels(&mut self, values: Option<Vec<u32>>) -> Result<()> {
        self.attrs.set_face_labels(values, self.faces.len())
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
    pub fn insert_point_attr(&mut self, name: &str, attr: MeshAttr3) -> Result<()> {
        self.attrs.insert_point_attr(name, attr, self.points.len())
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
    pub fn insert_face_attr(&mut self, name: &str, attr: MeshAttr3) -> Result<()> {
        self.attrs.insert_face_attr(name, attr, self.faces.len())
    }

    /// Remove and return the open-map per-point attribute stored under the given name.
    pub fn remove_point_attr(&mut self, name: &str) -> Option<MeshAttr3> {
        self.attrs.remove_point_attr(name)
    }

    /// Remove and return the open-map per-face attribute stored under the given name.
    pub fn remove_face_attr(&mut self, name: &str) -> Option<MeshAttr3> {
        self.attrs.remove_face_attr(name)
    }

    /// Replace the entire set of per-element attributes.
    ///
    /// # Arguments
    ///
    /// * `attrs`: the attribute set to attach, whose arrays must match the point and face counts
    ///
    /// returns: `Result<()>`, leaving the existing attributes untouched on failure
    pub fn set_attrs(&mut self, attrs: MeshAttrSet3) -> Result<()> {
        attrs.validate(self.points.len(), self.faces.len())?;
        self.attrs = attrs;
        Ok(())
    }

    /// Remove and return the entire set of per-element attributes, leaving the mesh with none.
    pub fn take_attrs(&mut self) -> MeshAttrSet3 {
        std::mem::take(&mut self.attrs)
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
            .field("attrs", &self.attrs)
            .finish()
    }
}

/// Verify that every face refers to a point which exists.
fn check_face_indices(faces: &[[u32; 3]], n_points: usize) -> Result<()> {
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

    #[test]
    fn new_accepts_a_valid_mesh() {
        let mesh = square_mesh();

        assert_eq!(mesh.point_count(), 4);
        assert_eq!(mesh.face_count(), 2);
        assert!(!mesh.is_empty());
        assert!(mesh.attrs().is_empty());
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
        let mut attrs = MeshAttrSet3::empty();
        attrs.set_point_colors(Some(vec![[0, 0, 0]]), 1).unwrap();

        // The attribute set is internally consistent for a one-point mesh, but not for this one.
        assert!(MeshData3::new_with_attrs(points, faces, attrs).is_err());
    }

    #[test]
    fn attribute_setters_supply_the_counts() -> Result<()> {
        let mut mesh = square_mesh();

        mesh.set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]))?;
        mesh.set_face_labels(Some(vec![7, 8]))?;
        mesh.insert_point_attr("confidence", MeshAttr3::Scalar(vec![0.5; 4]))?;
        mesh.insert_face_attr("material_index", MeshAttr3::Label(vec![1, 2]))?;

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

        let mut good = MeshAttrSet3::empty();
        good.set_point_normals(
            Some(vec![UnitVec3::new_normalize(Vector3::z()); 4]),
            mesh.point_count(),
        )?;
        mesh.set_attrs(good)?;
        assert!(mesh.point_normals().is_some());

        let mut bad = MeshAttrSet3::empty();
        bad.set_face_labels(Some(vec![1, 2, 3]), 3)?;
        assert!(mesh.set_attrs(bad).is_err());

        // The rejected set must have left the existing attributes in place.
        assert!(mesh.point_normals().is_some());

        Ok(())
    }

    #[test]
    fn take_attrs_leaves_the_mesh_bare() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_face_labels(Some(vec![7, 8]))?;

        let taken = mesh.take_attrs();

        assert_eq!(taken.face_labels().unwrap(), &[7, 8]);
        assert!(mesh.attrs().is_empty());
        assert!(mesh.face_labels().is_none());

        Ok(())
    }

    #[test]
    fn debug_summarizes_instead_of_dumping_buffers() {
        let text = format!("{:?}", square_mesh());
        assert!(text.contains("points: 4"), "{text}");
        assert!(text.contains("faces: 2"), "{text}");
    }
}
