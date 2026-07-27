//! This module contains an abstraction for a mesh of triangles, represented by vertices and their
//! indices into the vertex list.  This abstraction is built around the `TriMesh` type from the
//! `parry3d` crate.

pub mod algorithms;
mod collisions;
mod conformal;
pub mod data;
mod edges;
pub mod filtering;
pub mod half_edge;
mod measurement;
mod nav_structure;
mod outline;
mod queries;
pub mod sampling;
mod section;
mod uv_mapping;

use crate::common::{IndexMask, PCoords};
use crate::geom3::IsoExtensions3;
use crate::io::{deflate_bytes, u_bytes_to_mesh};
use crate::na::SVector;
use crate::{Iso3, Point2, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
pub use collisions::MeshCollisionSet;
pub use data::{Attr3, MeshAttrSet3, MeshData3};
pub use edges::MeshEdges;
pub use half_edge::HalfEdgeMesh;
pub use nav_structure::MeshNav;
use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::shape::{TriMesh, TriMeshFlags};
use parry3d_f64::{shape, transformation};
pub use uv_mapping::UvMapping;

#[cfg(feature = "ply")]
use crate::io::PlyWriteOpts;
#[cfg(feature = "stl")]
use crate::io::StlWriteOpts;
#[cfg(any(feature = "ply", feature = "stl"))]
use std::path::Path;

/// A struct that represents a point on the surface of a mesh, including the index of the face
/// on which it lies, its barycentric coordinates, and the point/normal representation in space.
/// This representation has no link back to the original mesh, so the face index and barycentric
/// coordinates will be invalid if (1) the mesh is modified, or (2) if you attempt to use them on
/// a different mesh.
#[derive(Debug, Clone, Copy)]
pub struct MeshSurfPoint {
    /// The index of the face on which this point lies.
    pub face_index: u32,

    /// The barycentric coordinates of the point on the face.
    pub bc: [f64; 3],

    /// The surface point (point + normal) corresponding to this barycentric coordinate.
    pub sp: SurfacePoint3,
}

impl MeshSurfPoint {
    /// Create a new `MeshSurfPoint` from the given face index, barycentric coordinates, and
    /// surface point.
    pub fn new(face_index: u32, bc: [f64; 3], sp: SurfacePoint3) -> Self {
        Self { face_index, bc, sp }
    }

    /// Get the point in space corresponding to this surface point.
    pub fn point(&self) -> Point3 {
        self.sp.point
    }

    /// Get the normal at this surface point.
    pub fn normal(&self) -> UnitVec3 {
        self.sp.normal
    }

    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            face_index: self.face_index,
            bc: self.bc,
            sp: iso * self.sp,
        }
    }
}

impl Default for MeshSurfPoint {
    fn default() -> Self {
        Self {
            face_index: 0,
            bc: [0.0, 0.0, 0.0],
            sp: SurfacePoint3::default(),
        }
    }
}

impl PCoords<3> for MeshSurfPoint {
    fn coords(&self) -> SVector<f64, 3> {
        self.sp.point.coords
    }
}

/// This is a triangle mesh optimized for collision detection and geometric queries. It is built on
/// top of the `parry3d` library's `TriMesh` type, which provides efficient storage and querying of
/// triangle meshes. This mesh has some basic functionality for interrogating its structure, and
/// some very basic functionality for editing.  However, it is not a structure optimized for
/// editing or modification.
///
/// # Relationship with `MeshData3`
///
/// This is the accelerated half of a complementary pair. `MeshData3` holds the same point and face
/// buffers with no spatial acceleration at all, and is the type to reach for when building or
/// editing mesh data, or when working with serialization. Converting between the two is done with
/// `from_data`/`to_data`, or the equivalent `TryFrom`/`From` implementations.
///
/// Both types carry the same `MeshAttrSet3` of per-element attributes, so a mesh loaded with
/// measured normals, colors, or uncertainties keeps them across the conversion in either direction.
///
/// # Invariants
///
/// The attribute arrays are validated against the point and face counts, exactly as they are on
/// `MeshData3`. That makes the point and face buffers inside the underlying `TriMesh` load bearing:
/// anything which renumbers a point or drops a face has to update the attributes to match, or
/// refuse to run. This is why the `merge_duplicates` and `delete_degenerate` options on
/// `new_with_options` are only reachable from a constructor which attaches no attributes.
#[derive(Clone)]
pub struct Mesh3 {
    shape: TriMesh,
    attrs: MeshAttrSet3,
    is_solid: bool,
    uv: Option<UvMapping>,
}

// ===============================================================================================
// Core access
// ===============================================================================================

impl Mesh3 {
    /// Get a reference to the AABB of the underlying mesh in the local coordinate system.
    pub fn aabb(&self) -> Aabb {
        self.shape.local_aabb()
    }

    /// Gets a reference to the underlying `TriMesh` object to provide direct access to
    /// the `parry3d` API.
    pub fn tri_mesh(&self) -> &TriMesh {
        &self.shape
    }

    /// Return a flag indicating whether the mesh is considered "solid" or not for the purposes of
    /// distance queries. If a mesh is "solid", then distance queries for points on the inside of
    /// the mesh will return a zero distance.
    pub fn is_solid(&self) -> bool {
        self.is_solid
    }

    /// Get a reference to the points of the mesh.
    pub fn points(&self) -> &[Point3] {
        self.shape.vertices()
    }

    /// Get a reference to the face indices of the mesh.
    pub fn faces(&self) -> &[[u32; 3]] {
        self.shape.indices()
    }

    /// Get the number of points in the mesh.
    pub fn point_count(&self) -> usize {
        self.shape.vertices().len()
    }

    /// Get the number of faces in the mesh.
    pub fn face_count(&self) -> usize {
        self.shape.indices().len()
    }
}

// ===============================================================================================
// Attribute access
// ===============================================================================================

impl Mesh3 {
    /// Get a reference to the full set of per-element attributes attached to this mesh.
    pub fn attrs(&self) -> &MeshAttrSet3 {
        &self.attrs
    }

    /// Get the per-point unit normals, if present.
    ///
    /// These are whatever was measured or stored, not a computed quantity. Use
    /// `compute_point_normals` to derive normals from the faces.
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
    pub fn point_attr(&self, name: &str) -> Option<&Attr3> {
        self.attrs.point_attr(name)
    }

    /// Get the open-map per-face attribute stored under the given name, if present.
    pub fn face_attr(&self, name: &str) -> Option<&Attr3> {
        self.attrs.face_attr(name)
    }
}

// ===============================================================================================
// Attribute mutation
// ===============================================================================================

impl Mesh3 {
    /// Set or clear the per-point unit normals.
    ///
    /// # Arguments
    ///
    /// * `values`: the normals to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_normals(&mut self, values: Option<Vec<UnitVec3>>) -> Result<()> {
        self.attrs.set_point_normals(values, self.point_count())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.attrs.set_point_colors(values, self.point_count())
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
        self.attrs.set_point_stdev(values, self.point_count())
    }

    /// Set or clear the per-face RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.attrs.set_face_colors(values, self.face_count())
    }

    /// Set or clear the per-face labels.
    ///
    /// # Arguments
    ///
    /// * `values`: the labels to store, or `None` to clear them. Must match the face count.
    ///
    /// returns: `Result<()>`
    pub fn set_face_labels(&mut self, values: Option<Vec<u32>>) -> Result<()> {
        self.attrs.set_face_labels(values, self.face_count())
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
        self.attrs
            .insert_point_attr(name, attr, self.shape.vertices().len())
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
        self.attrs
            .insert_face_attr(name, attr, self.shape.indices().len())
    }

    /// Remove and return the open-map per-point attribute stored under the given name.
    pub fn remove_point_attr(&mut self, name: &str) -> Option<Attr3> {
        self.attrs.remove_point_attr(name)
    }

    /// Remove and return the open-map per-face attribute stored under the given name.
    pub fn remove_face_attr(&mut self, name: &str) -> Option<Attr3> {
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
        attrs.validate(self.point_count(), self.face_count())?;
        self.attrs = attrs;
        Ok(())
    }

    /// Remove and return the entire set of per-element attributes, leaving the mesh with none.
    pub fn take_attrs(&mut self) -> MeshAttrSet3 {
        std::mem::take(&mut self.attrs)
    }
}

// ===============================================================================================
// General creation methods
// ===============================================================================================
impl Mesh3 {
    /// Create a new mesh from a list of vertices and a list of triangles.  Additional options can
    /// be set to merge duplicate vertices and delete degenerate triangles.
    ///
    /// The resulting mesh carries no attributes. The cleanup options renumber points and drop
    /// faces, so there is no correct way to carry an attribute array through them, which is why
    /// they live here rather than on the `MeshData3` conversion path.
    ///
    /// # Arguments
    ///
    /// * `vertices`:
    /// * `triangles`:
    /// * `is_solid`:
    /// * `merge_duplicates`:
    /// * `delete_degenerate`:
    /// * `uv`:
    ///
    /// returns: Result<Mesh3, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new_with_options(
        vertices: Vec<Point3>,
        triangles: Vec<[u32; 3]>,
        is_solid: bool,
        merge_duplicates: bool,
        delete_degenerate: bool,
        uv: Option<Vec<Point2>>,
    ) -> Result<Self> {
        let mut flags = TriMeshFlags::empty();
        if merge_duplicates {
            flags |= TriMeshFlags::MERGE_DUPLICATE_VERTICES;
            flags |= TriMeshFlags::DELETE_DUPLICATE_TRIANGLES;
        }
        if delete_degenerate {
            flags |= TriMeshFlags::DELETE_BAD_TOPOLOGY_TRIANGLES;
            flags |= TriMeshFlags::DELETE_DEGENERATE_TRIANGLES;
        }

        let uv_mapping = if let Some(uv) = uv {
            Some(UvMapping::new(uv, triangles.clone())?)
        } else {
            None
        };

        let shape = TriMesh::with_flags(vertices, triangles, flags)?;
        Ok(Self {
            shape,
            attrs: MeshAttrSet3::empty(),
            is_solid,
            uv: uv_mapping,
        })
    }

    pub fn new(vertices: Vec<Point3>, triangles: Vec<[u32; 3]>, is_solid: bool) -> Self {
        let shape = TriMesh::new(vertices, triangles).expect("Failed to create TriMesh");
        Self {
            shape,
            attrs: MeshAttrSet3::empty(),
            is_solid,
            uv: None,
        }
    }
    pub fn new_take_trimesh(shape: TriMesh, is_solid: bool) -> Self {
        Self {
            shape,
            attrs: MeshAttrSet3::empty(),
            is_solid,
            uv: None,
        }
    }
}

// ===============================================================================================
// Serialization
// ===============================================================================================
//
// Every one of these goes through `MeshData3`, which is where the format support lives. Loading
// costs a bounding volume hierarchy build on top of the read, and saving costs a copy of the point
// and face buffers, for the reason given on `to_data`. Work directly with `MeshData3` when you have
// no need for the acceleration structure.

impl Mesh3 {
    /// Load a triangle mesh from a PLY file, preserving every property the file carries.
    ///
    /// The file must have a `vertex` element with `x`, `y`, and `z` properties, and a `face`
    /// element with at least one face. A PLY point cloud is refused, since there is nothing for an
    /// acceleration structure to be built over; load those with `PointCloudData3::load_ply`.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the PLY file
    /// * `is_solid`: whether distance queries should treat points inside the mesh as being at zero
    ///   distance
    ///
    /// returns: `Result<Mesh3>`
    #[cfg(feature = "ply")]
    pub fn load_ply(path: &Path, is_solid: bool) -> Result<Self> {
        Self::from_data(MeshData3::load_ply(path)?, is_solid)
    }

    /// Write this mesh to a PLY file, preserving every attribute it carries.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `opts`: encoding and header options
    ///
    /// returns: `Result<()>`
    #[cfg(feature = "ply")]
    pub fn save_ply(&self, path: &Path, opts: &PlyWriteOpts) -> Result<()> {
        self.to_data().save_ply(path, opts)
    }

    /// Load a triangle mesh from an STL file, in either the ascii or binary encoding.
    ///
    /// STL is a triangle soup with no point identity, so the points are recovered by welding on
    /// exact coordinate equality. See `load_stl_mesh_data` for what that does and does not recover,
    /// and for the precision the format costs you.
    ///
    /// The cleanup options renumber points and drop faces, which is why they are available here and
    /// not on the PLY path: an STL carries no attributes for the renumbering to invalidate.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the STL file
    /// * `is_solid`: whether distance queries should treat points inside the mesh as being at zero
    ///   distance
    /// * `merge_duplicates`: additionally merge points which compare equal as `f64`, which after
    ///   the exact weld the reader already performs means only merging `0.0` with `-0.0`
    /// * `delete_degenerate`: drop triangles with zero area or bad topology
    ///
    /// returns: `Result<Mesh3>`
    #[cfg(feature = "stl")]
    pub fn load_stl(
        path: &Path,
        is_solid: bool,
        merge_duplicates: bool,
        delete_degenerate: bool,
    ) -> Result<Self> {
        let (points, faces, _) = MeshData3::load_stl(path)?.into_parts();
        Self::new_with_options(
            points,
            faces,
            is_solid,
            merge_duplicates,
            delete_degenerate,
            None,
        )
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
        self.to_data().save_stl(path, opts)
    }
}

// ===============================================================================================
// Conversion to and from MeshData3
// ===============================================================================================

impl Mesh3 {
    /// Build an accelerated mesh from plain mesh data, taking ownership of its buffers and its
    /// attributes.
    ///
    /// The point and face buffers are moved into the underlying `TriMesh` unchanged, so every
    /// index keeps the meaning it had and the attribute arrays stay valid. The cost paid here is
    /// the bounding volume hierarchy build.
    ///
    /// # Arguments
    ///
    /// * `data`: the mesh data to consume
    /// * `is_solid`: whether distance queries should treat points inside the mesh as being at zero
    ///   distance
    ///
    /// returns: `Result<Mesh3>`, failing if the mesh has no faces, since there is nothing to build
    /// an acceleration structure over
    pub fn from_data(data: MeshData3, is_solid: bool) -> Result<Self> {
        let (points, faces, attrs) = data.into_parts();

        if faces.is_empty() {
            return Err(
                "Cannot build a Mesh3 from mesh data with no faces. A MeshData3 is allowed \
                        to hold points without faces, but there is nothing for an acceleration \
                        structure to be built over."
                    .into(),
            );
        }

        let shape = TriMesh::new(points, faces)?;
        Ok(Self {
            shape,
            attrs,
            is_solid,
            uv: None,
        })
    }

    /// Copy this mesh's buffers and attributes out into plain mesh data.
    ///
    /// Unlike `from_data`, this **clones** the point and face buffers. `parry3d` gives no way to
    /// move them back out of a `TriMesh`, so the copy is unavoidable. Use `into_data` if you are
    /// done with the accelerated mesh, which at least drops the acceleration structure.
    ///
    /// returns: `MeshData3`
    pub fn to_data(&self) -> MeshData3 {
        MeshData3::new_with_attrs(
            self.points().to_vec(),
            self.faces().to_vec(),
            self.attrs.clone(),
        )
        .expect(
            "A Mesh3's attributes are validated against its point and face counts on every \
             operation which can change them, so rebuilding a MeshData3 from them cannot fail",
        )
    }

    /// Consume this mesh and return its buffers and attributes as plain mesh data.
    ///
    /// This still copies the point and face buffers, for the reason given on `to_data`, but it
    /// discards the acceleration structure rather than leaving it alive alongside the copy.
    ///
    /// returns: `MeshData3`
    pub fn into_data(self) -> MeshData3 {
        self.to_data()
    }
}

impl TryFrom<MeshData3> for Mesh3 {
    type Error = Box<dyn std::error::Error>;

    /// Build an accelerated mesh which is **not** solid. Use `Mesh3::from_data` to choose.
    fn try_from(value: MeshData3) -> Result<Self> {
        Mesh3::from_data(value, false)
    }
}

impl From<Mesh3> for MeshData3 {
    fn from(value: Mesh3) -> Self {
        value.into_data()
    }
}

// ===============================================================================================
// Mutation/Transformation
// ===============================================================================================
impl Mesh3 {
    /// Transform the mesh in place by applying the given transformation to all vertices.
    ///
    /// Any stored point normals and `Vector` attributes are rotated with the geometry. Because
    /// those hold directions rather than positions, only the rotation component has any effect.
    pub fn transform_by(&mut self, transform: &Iso3) {
        self.shape.transform_vertices(transform);
        self.attrs.transform_in_place(transform);
    }

    /// Returns a new mesh with all vertices transformed by the given isometry, leaving the
    /// original unchanged.
    pub fn new_transformed_by(&self, transform: &Iso3) -> Self {
        let mut result = self.clone();
        result.transform_by(transform);
        result
    }

    /// Create a new mesh by scaling all points uniformly about the origin.
    ///
    /// Any stored point standard deviations are scaled with the geometry, since they are lengths in
    /// the mesh's own units. Normals are directions and a uniform scale does not move them.
    ///
    /// A **negative** factor is a mirror, which is orientation-reversing: the face normals implied
    /// by the winding keep pointing the same way in space while the solid turns inside out around
    /// them. The winding is therefore reversed and any stored point normals negated, matching what
    /// `MeshData3::scale_copy` does.
    ///
    /// # Arguments
    ///
    /// * `scale`: the factor to scale by, which must be finite and non-zero
    ///
    /// returns: `Result<Mesh3>`
    pub fn scale_copy(&self, scale: f64) -> Result<Self> {
        if !scale.is_finite() || scale == 0.0 {
            return Err(format!(
                "A scale factor must be finite and non-zero, but {scale} was given. Scaling by zero \
                 would collapse the mesh to a point irrecoverably."
            )
            .into());
        }

        let new_shape = self
            .shape
            .clone()
            .scaled(&Vector3::new(scale, scale, scale));

        let mut result = Mesh3::new_take_trimesh(new_shape, self.is_solid);
        result.attrs = self.attrs.clone();
        result.attrs.scale_in_place(scale);

        if scale < 0.0 {
            result.flip_normals();
        }

        Ok(result)
    }

    /// Create a new mesh by offsetting each vertex along its smoothed vertex normal.
    ///
    /// The offset is applied as `vertex + offset * normal`, where the normal is the
    /// normalized per-vertex normal computed from adjacent face normals.
    ///
    /// Positive offsets expand the mesh outward; negative offsets shrink it inward.
    /// The original mesh is not modified.
    ///
    /// # Arguments
    ///
    /// * `offset`: The distance to offset each vertex along its normal.
    ///
    /// returns: Mesh3
    pub fn new_offset_vertices(&self, offset: f64) -> Self {
        // These are already normalized
        let normals = self.get_vertex_normals();

        let updated = self
            .points()
            .iter()
            .zip(normals.iter())
            .map(|(v, n)| v + offset * n)
            .collect();

        Self::new(updated, self.faces().to_vec(), self.is_solid)
    }

    /// Reverse the winding order of every face, turning the surface inside out.
    ///
    /// Any stored point normals are negated to match, since the direction the surface faces has
    /// changed.
    pub fn flip_normals(&mut self) {
        self.shape.reverse();
        self.attrs.flip_in_place();
    }
}

// ===============================================================================================
// Unsorted
// ===============================================================================================

impl Mesh3 {
    pub fn calc_edges(&self) -> Result<MeshEdges<'_>> {
        MeshEdges::new(self)
    }

    /// Return a convex hull of the points in the mesh.
    ///
    /// The hull is new topology: its points are a subset of this mesh's but its faces are not, and
    /// `parry3d` gives no mapping back from one to the other. There is therefore no correct value
    /// to carry any attribute forward with, so a mesh which has any is refused unless the caller
    /// says the loss is acceptable.
    ///
    /// # Arguments
    ///
    /// * `allow_attribute_loss`: accept the loss of every attribute this mesh carries
    ///
    /// returns: `Result<Mesh3>`
    pub fn convex_hull(&self, allow_attribute_loss: bool) -> Result<Self> {
        self.check_attribute_loss("a convex hull", allow_attribute_loss)?;

        let (vertices, faces) = transformation::convex_hull(self.shape.vertices());
        Ok(Self::new(vertices, faces, true))
    }

    /// Append another mesh onto the end of this one.
    ///
    /// The other mesh's points are added after this one's and its faces are re-indexed to match. No
    /// points are welded and no faces are merged, so a shared boundary between the two meshes stays
    /// as two coincident sets of points.
    ///
    /// Attributes are all-or-nothing: a typed field or an open-map key present on one side and
    /// absent on the other is an error, because there is no correct value to pad the missing side
    /// with. The whole append is validated before anything is modified, so a failure leaves this
    /// mesh untouched.
    ///
    /// # Arguments
    ///
    /// * `other`: the mesh to append
    ///
    /// returns: `Result<()>`
    pub fn append(&mut self, other: &Mesh3) -> Result<()> {
        // For now, both meshes must have an empty UV mapping
        if self.uv.is_some() || other.uv.is_some() {
            return Err("Cannot append meshes with UV mappings".into());
        }

        // The attributes are merged into a copy first, because that is the step which can fail.
        // Only once it has succeeded is either the geometry or this mesh's attribute set touched.
        let mut merged = self.attrs.clone();
        merged.extend_from(&other.attrs)?;

        self.shape.append(&other.shape);
        self.attrs = merged;

        Ok(())
    }

    /// Verify that the caller has accepted the loss of this mesh's attributes, for an operation
    /// which cannot carry them.
    ///
    /// See `MeshData3::check_attribute_loss` for why this is an error rather than a warning.
    ///
    /// # Arguments
    ///
    /// * `operation`: what is about to discard them, used in the error message
    /// * `allow_loss`: whether the caller has accepted the loss
    ///
    /// returns: `Result<()>`
    pub fn check_attribute_loss(&self, operation: &str, allow_loss: bool) -> Result<()> {
        if allow_loss || self.attrs.is_empty() {
            return Ok(());
        }

        let mut lost = self.attrs.point_attr_labels();
        lost.extend(self.attrs.face_attr_labels());

        Err(format!(
            "Taking {} would discard the attributes on this mesh ({}), because it produces \
             topology with no index mapping back to the original. Pass `allow_attribute_loss` to \
             accept this.",
            operation,
            lost.join(", ")
        )
        .into())
    }

    pub fn uv(&self) -> Option<&UvMapping> {
        self.uv.as_ref()
    }

    pub fn uv_to_3d(&self, uv: &Point2) -> Option<MeshSurfPoint> {
        let (i, bc) = self.uv()?.triangle(uv)?;
        self.at_barycentric(i, bc).ok()
    }

    pub fn project_to_uv(&self, p: &impl PCoords<3>) -> Option<Point2> {
        let uv_map = self.uv()?;
        let mp = self.surf_closest_to(p);
        Some(uv_map.point(mp.face_index, mp.bc))
    }

    pub fn uv_with_tol(
        &self,
        point: &Point3,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(Point2, f64)> {
        if let Some(uv_map) = self.uv() {
            let point = if let Some(transform) = transform {
                transform * point
            } else {
                *point
            };

            if let Some((prj, id, loc)) = self.project_with_tol(&point, max_dist, max_angle, None) {
                let triangle = self.shape.triangle(id);
                if let Some(normal) = triangle.normal() {
                    let uv = uv_map.point(id, loc.barycentric_coordinates().unwrap());
                    // Now find the depth
                    let sp = SurfacePoint3::new(prj.point, normal);
                    Some((uv, sp.scalar_projection(&point)))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
    /// Create a new `MeshNav` structure for this mesh. This structure is used to efficiently
    /// navigate the mesh through edges and faces.  It is recommended to use this if you will be
    /// performing multiple structural queries on the mesh, so that the structure does not need to
    /// be recomputed each time.
    pub fn nav(&self) -> MeshNav<'_> {
        MeshNav::new(self)
    }

    /// Calculates the patches in the mesh. If you are going to be doing multiple queries of the
    /// structure of the mesh, either use the half-edge representation, or generate a `MeshNav`
    /// through the `nav()` method to avoid having to recompute the mesh structure each time.
    ///
    /// # Arguments
    ///
    /// * `mask`:
    ///
    /// returns: Result<Vec<IndexMask, Global>, Box<dyn Error, Global>>
    pub fn get_patches(&self, mask: Option<&IndexMask>) -> Result<Vec<IndexMask>> {
        let nav = self.nav();
        nav.patches(mask)
    }

    /// Gets the boundary points of each patch in the mesh.  This function will return a list of
    /// lists of points, where each list of points is the boundary of a patch.  Note that this
    /// function will not work on non-manifold meshes.
    ///
    /// returns: Result<Vec<Vec<usize, Global>, Global>>
    pub fn get_patch_boundary_points(&self) -> Result<Vec<Vec<Point3>>> {
        let edges = MeshEdges::new(self)?;

        let mut b_loops = Vec::new();
        for b_loop in edges.boundary_loops.iter() {
            b_loops.push(
                b_loop
                    .iter()
                    .map(|vi| self.points()[*vi as usize])
                    .collect(),
            );
        }

        Ok(b_loops)
    }

    pub fn get_face_normals(&self) -> Result<Vec<UnitVec3>> {
        let mut result = Vec::new();
        for t in self.shape.triangles() {
            if let Some(n) = t.normal() {
                result.push(n);
            } else {
                return Err("Failed to get normal".into());
            }
        }

        Ok(result)
    }

    /// Calculates and returns the areas of all triangular faces in the shape.
    ///
    /// This function iterates over all the triangles in the shape, computes the area
    /// of each triangle using the `area` method, and collects the results into a vector.
    pub fn get_face_areas(&self) -> Vec<f64> {
        self.shape.triangles().map(|t| t.area()).collect::<Vec<_>>()
    }

    /// Calculates and returns the centroid of all triangular faces in the mesh.
    pub fn get_face_centers(&self) -> Vec<Point3> {
        self.shape
            .triangles()
            .map(|t| t.center())
            .collect::<Vec<_>>()
    }

    /// Compute smooth per-vertex normals by averaging the normals of all adjacent triangles
    /// weighted by triangle area. At the end of the computation, the normals are normalized to
    /// have unit length.
    ///
    /// Be aware that vertices that are not referenced by any valid triangle keep the zero vector.
    ///
    /// Also, be aware that this may not produce the results you expect for meshes with large flat
    /// surfaces represented by multiple triangles. For example, on a cube mesh, not all corner
    /// vertices will point along the diagonals, since each vertex will have some faces where it
    /// touches two triangles and may have some faces where it touches only one triangle, making
    /// the weights uneven.
    pub fn get_vertex_normals(&self) -> Vec<Vector3> {
        let mut sums: Vec<Vector3> = vec![Vector3::new(0.0, 0.0, 0.0); self.shape.vertices().len()];
        let mut weights = vec![0.0; self.shape.vertices().len()];

        for (face_i, tri) in self.shape.triangles().enumerate() {
            let indices = self.shape.indices()[face_i];
            let a = tri.area();
            if let Some(n) = tri.normal() {
                for i in indices {
                    sums[i as usize] += n.into_inner() * a;
                    weights[i as usize] += a;
                }
            }
        }

        // Normalize the normals
        for i in 0..sums.len() {
            if weights[i] > 0.0 {
                let v = sums[i] / weights[i];
                sums[i] = v.normalize();
            }
        }

        sums
    }
}

// ===============================================================================================
// Shape creation methods
// ===============================================================================================

impl Mesh3 {
    pub fn create_cone(half_height: f64, radius: f64, steps: usize) -> Self {
        let cone = shape::Cone::new(half_height, radius);
        let (vertices, faces) = cone.to_trimesh(steps as u32);

        Self::new(vertices, faces, true)
    }

    pub fn create_capsule(
        p0: &Point3,
        p1: &Point3,
        radius: f64,
        n_theta: usize,
        n_phi: usize,
    ) -> Self {
        let capsule = shape::Capsule::new(*p0, *p1, radius);
        let (vertices, faces) = capsule.to_trimesh(n_theta as u32, n_phi as u32);

        Self::new(vertices, faces, true)
    }

    /// Create a spherical mesh centered at the origin. The `n_theta` and `n_phi` parameters control
    /// the tessellation density.
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius of the sphere.
    /// * `n_theta` - Number of subdivisions around the polar direction.
    /// * `n_phi` - Number of subdivisions around the azimuthal direction.
    ///
    /// returns: Mesh3
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Mesh3;
    /// use approx::assert_relative_eq;
    ///
    /// let n_t = 14;
    /// let n_p = 15;
    /// let sphere = Mesh3::create_sphere(1.0, n_t, n_p);
    ///
    /// assert_eq!(sphere.points().len(), n_t * (n_p - 1) + 2);
    ///
    /// // Verify that the vertices are on the surface of the unit sphere.
    /// for vertex in sphere.points() {
    ///     let dist_from_origin = vertex.coords.norm();
    ///     assert_relative_eq!(dist_from_origin, 1.0)
    /// }
    /// ```
    pub fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let sphere = shape::Ball::new(radius);
        let (vertices, faces) = sphere.to_trimesh(n_theta as u32, n_phi as u32);

        Self::new(vertices, faces, true)
    }

    /// Create a box mesh with the given dimensions, centered at the origin.
    ///
    /// # Arguments
    ///
    /// * `length`: the dimension of the box in the x direction
    /// * `width`: the dimension of the box in the y direction
    /// * `height`: the dimension of the box in the z direction
    /// * `is_solid`: whether the box is solid or hollow, used for some specific distance queries
    ///   in the underlying parry library
    ///
    /// returns: Mesh3
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Mesh3;
    /// use approx::assert_relative_eq;
    /// let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
    /// assert_relative_eq!(mesh.aabb().maxs.x, 1.0);
    /// assert_relative_eq!(mesh.aabb().maxs.y, 2.0);
    /// assert_relative_eq!(mesh.aabb().maxs.z, 3.0);
    /// assert_relative_eq!(mesh.aabb().mins.x, -1.0);
    /// assert_relative_eq!(mesh.aabb().mins.y, -2.0);
    /// assert_relative_eq!(mesh.aabb().mins.z, -3.0);
    /// ```
    pub fn create_box(length: f64, width: f64, height: f64, is_solid: bool) -> Self {
        let bx = shape::Cuboid::new(Vector3::new(length / 2.0, width / 2.0, height / 2.0));
        let (vertices, triangles) = bx.to_trimesh();
        Self::new(vertices, triangles, is_solid)
    }

    /// Create a cylindrical mesh centered at the origin and aligned with the local `y` axis.
    /// The `radius` controls the cylinder radius, `height` its full height, and `steps`
    /// controls the tessellation density around the circumference.
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius of the cylinder.
    /// * `height` - Full height of the cylinder (along the y-axis).
    /// * `steps` - Number of subdivisions around the cylinder axis.
    ///
    /// returns: Mesh3
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Mesh3, Point3};
    /// use approx::assert_relative_eq;
    ///
    /// let cyl = Mesh3::create_cylinder(1.0, 4.0, 16);
    ///
    /// assert_relative_eq!(cyl.aabb().mins.z, -1.0);
    /// assert_relative_eq!(cyl.aabb().maxs.z,  1.0);
    /// assert_relative_eq!(cyl.aabb().mins.x, -1.0);
    /// assert_relative_eq!(cyl.aabb().maxs.x,  1.0);
    /// assert_relative_eq!(cyl.aabb().mins.y, -2.0);
    /// assert_relative_eq!(cyl.aabb().maxs.y,  2.0);
    ///
    /// for vertex in cyl.points() {
    ///     let proj = Point3::new(vertex.x, 0.0, vertex.z);
    ///     assert_relative_eq!(proj.coords.norm(), 1.0);
    /// }
    /// ```
    pub fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        let cyl = shape::Cylinder::new(height / 2.0, radius);
        let (vertices, faces) = cyl.to_trimesh(steps as u32);

        Self::new(vertices, faces, true)
    }

    pub fn create_rect_beam_between(
        p0: &Point3,
        p1: &Point3,
        width: f64,
        height: f64,
        up: &Vector3,
    ) -> Result<Self> {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;
        let box_geom = shape::Cuboid::new(Vector3::new(width / 2.0, height / 2.0, v.norm() / 2.0));

        // I think this is OK?
        let transform = Iso3::from_basis_zy(&v, up, Some(pc))?;

        let (vertices, faces) = box_geom.to_trimesh();
        let mut mesh = Self::new(vertices, faces, true);
        mesh.transform_by(&transform);
        Ok(mesh)
    }

    pub fn create_cylinder_between(p0: &Point3, p1: &Point3, radius: f64, steps: usize) -> Self {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;
        let cyl = shape::Cylinder::new(v.norm() / 2.0, radius);

        // I think this is OK?
        let transform = Iso3::from_basis_yz(&v, &Vector3::z(), Some(pc))
            .unwrap_or(Iso3::from_basis_yx(&v, &Vector3::x(), Some(pc)).unwrap());

        let (vertices, faces) = cyl.to_trimesh(steps as u32);
        let mut mesh = Self::new(vertices, faces, true);
        mesh.transform_by(&transform);
        mesh
    }

    /// Create a flat, filled circle mesh lying in the XY plane, centered at the origin, with the
    /// normal pointing along +Z. The mesh is a triangle fan from the center to `segments` evenly
    /// spaced perimeter vertices.
    ///
    /// # Arguments
    ///
    /// * `radius` - Radius of the circle.
    /// * `segments` - Number of perimeter vertices (and triangles). Must be at least 3.
    ///
    /// returns: Mesh3
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Mesh3;
    ///
    /// let circle = Mesh3::create_circle(1.0, 32);
    /// assert_eq!(circle.points().len(), 33); // center + 32 perimeter
    /// assert_eq!(circle.faces().len(), 32);
    /// ```
    pub fn create_circle(radius: f64, segments: usize) -> Self {
        use std::f64::consts::TAU;
        let mut vertices = Vec::with_capacity(segments + 1);
        let mut faces = Vec::with_capacity(segments);

        vertices.push(Point3::origin());
        for i in 0..segments {
            let angle = TAU * (i as f64) / (segments as f64);
            vertices.push(Point3::new(radius * angle.cos(), radius * angle.sin(), 0.0));
        }

        for i in 0..segments {
            let a = (i + 1) as u32;
            let b = ((i + 1) % segments + 1) as u32;
            faces.push([0u32, a, b]);
        }

        Self::new(vertices, faces, false)
    }

    /// Load a Stanford bunny mesh embedded in the binary with 453 vertices and 948 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res3.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res4() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_4.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 1889 vertices and 3851 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res3.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res3() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_3.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 8171 vertices and 16301 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same
    /// as the corresponding `bun_zipper_res2.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum vertex deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res2() -> Self {
        let bytes = include_bytes!("../../tests/data/stanford_bun_2.umesh.gz");
        u_bytes_to_mesh(&deflate_bytes(bytes).unwrap()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Plane3;
    use crate::tests::stanford_bun_4;
    use approx::assert_relative_eq;
    use parry3d_f64::query::SplitResult;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn vertex_normals_match_vertex_count_and_are_normalized() {
        let mesh = stanford_bun_4();
        let normals = mesh.get_vertex_normals();

        assert_eq!(normals.len(), mesh.points().len());

        for normal in normals {
            assert_relative_eq!(normal.norm(), 1.0, epsilon = 1.0e-12);
        }
    }

    #[test]
    fn scale_copy_scales_spherical_radius() -> Result<()> {
        let radius = 1.0;
        let scale = 2.5;
        let mesh = Mesh3::create_sphere(radius, 100, 100);
        let scaled = mesh.scale_copy(scale)?;

        assert_eq!(mesh.points().len(), scaled.points().len());

        for vertex in scaled.points() {
            assert_relative_eq!(vertex.coords.norm(), radius * scale, epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn new_offset_vertices_preserves_spherical_radius() {
        let radius = 1.0;
        let offset = 0.1;
        let mesh = Mesh3::create_sphere(radius, 100, 100);
        let offset_mesh = mesh.new_offset_vertices(offset);

        assert_eq!(mesh.points().len(), offset_mesh.points().len());

        for vertex in offset_mesh.points() {
            assert_relative_eq!(vertex.coords.norm(), radius + offset, epsilon = 1.0e-5);
        }
    }

    // ===========================================================================================
    // Conversion to and from MeshData3
    // ===========================================================================================

    /// A single triangle in the xy plane, wound so that its normal is +z, carrying one attribute of
    /// every kind an operation might have to touch.
    fn attributed_data() -> MeshData3 {
        let mut data = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )
        .unwrap();

        data.set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 3]))
            .unwrap();
        data.set_point_stdev(Some(vec![0.1, 0.2, 0.3])).unwrap();
        data.set_face_labels(Some(vec![7])).unwrap();
        data.insert_point_attr("confidence", Attr3::Scalar(vec![0.5, 0.6, 0.7]))
            .unwrap();
        data.insert_point_attr("principal_dir", Attr3::Vector(vec![Vector3::x(); 3]))
            .unwrap();

        data
    }

    #[test]
    fn round_trip_through_mesh_data_preserves_everything() -> Result<()> {
        let before = attributed_data();
        let mesh = Mesh3::from_data(before.clone(), true)?;

        assert!(mesh.is_solid());
        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);

        let after = mesh.to_data();

        assert_eq!(after.points(), before.points());
        assert_eq!(after.faces(), before.faces());
        assert_eq!(after.attrs(), before.attrs());

        Ok(())
    }

    /// The buffers have to come across index-identical, or every attribute array would be pointing
    /// at the wrong elements.
    #[test]
    fn from_data_does_not_renumber_the_buffers() -> Result<()> {
        let data = stanford_bun_4().to_data();
        let mesh = Mesh3::from_data(data.clone(), false)?;

        assert_eq!(mesh.points(), data.points());
        assert_eq!(mesh.faces(), data.faces());

        Ok(())
    }

    /// `MeshData3` allows points without faces, but there is nothing to accelerate.
    #[test]
    fn from_data_rejects_a_mesh_with_no_faces() -> Result<()> {
        let data = MeshData3::new(
            vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
            Vec::new(),
        )?;
        assert!(Mesh3::from_data(data, false).is_err());

        Ok(())
    }

    #[test]
    fn conversion_traits_agree_with_the_named_methods() -> Result<()> {
        let mesh = Mesh3::try_from(attributed_data())?;

        // The trait impl has no way to be told, so it picks the non-solid reading.
        assert!(!mesh.is_solid());

        let data: MeshData3 = mesh.into();
        assert_eq!(data.attrs(), attributed_data().attrs());

        Ok(())
    }

    // ===========================================================================================
    // Attributes under in-place mutation
    // ===========================================================================================

    #[test]
    fn transform_rotates_the_direction_attributes() -> Result<()> {
        let mut mesh = Mesh3::from_data(attributed_data(), false)?;

        // A quarter turn about +z, which maps +x onto +y, plus a translation.
        let iso = Iso3::new(Vector3::new(10.0, 0.0, 0.0), Vector3::z() * FRAC_PI_2);
        mesh.transform_by(&iso);

        assert_relative_eq!(
            mesh.points()[1],
            Point3::new(10.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // +z lies on the axis of rotation and does not move.
        assert_relative_eq!(
            mesh.point_normals().unwrap()[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            mesh.point_attr("principal_dir")
                .unwrap()
                .as_vector()
                .unwrap()[0],
            Vector3::y(),
            epsilon = 1.0e-12
        );

        // Scalars are untouched.
        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3]);

        Ok(())
    }

    #[test]
    fn scaling_scales_the_standard_deviations() -> Result<()> {
        let mesh = Mesh3::from_data(attributed_data(), false)?;
        let scaled = mesh.scale_copy(25.4)?;

        for (actual, expected) in scaled
            .point_stdev()
            .unwrap()
            .iter()
            .zip([0.1, 0.2, 0.3].iter())
        {
            assert_relative_eq!(*actual, expected * 25.4, epsilon = 1.0e-12);
        }

        // Normals are directions and a uniform scale does not move them.
        assert_relative_eq!(
            scaled.point_normals().unwrap()[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn flipping_negates_the_stored_normals() -> Result<()> {
        let mut mesh = Mesh3::from_data(attributed_data(), false)?;
        mesh.flip_normals();

        assert_relative_eq!(
            mesh.point_normals().unwrap()[0].into_inner(),
            -Vector3::z(),
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn append_unions_the_attributes() -> Result<()> {
        let mut mesh = Mesh3::from_data(attributed_data(), false)?;
        let other = Mesh3::from_data(attributed_data(), false)?;

        mesh.append(&other)?;

        assert_eq!(mesh.point_count(), 6);
        assert_eq!(mesh.face_count(), 2);
        assert_eq!(mesh.faces()[1], [3, 4, 5]);
        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.1, 0.2, 0.3]);
        assert_eq!(mesh.face_labels().unwrap(), &[7, 7]);
        assert_eq!(mesh.point_attr("confidence").unwrap().len(), 6);

        mesh.attrs()
            .validate(mesh.point_count(), mesh.face_count())?;

        Ok(())
    }

    /// An attribute on one side and not the other has no correct padding value, and a failure has
    /// to leave the target completely untouched rather than half appended.
    #[test]
    fn append_rejects_a_mismatch_without_modifying_the_target() -> Result<()> {
        let mut mesh = Mesh3::from_data(attributed_data(), false)?;

        let mut bare_data = attributed_data();
        bare_data.set_point_stdev(None)?;
        let other = Mesh3::from_data(bare_data, false)?;

        assert!(mesh.append(&other).is_err());

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);
        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3]);

        Ok(())
    }

    #[test]
    fn appending_two_bare_meshes_works() -> Result<()> {
        let mut mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let other = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let n = mesh.point_count();

        mesh.append(&other)?;
        assert_eq!(mesh.point_count(), n * 2);

        Ok(())
    }

    // ===========================================================================================
    // Attributes on derived meshes
    // ===========================================================================================

    /// A box with a distinct label on every face and a distinct standard deviation on every point,
    /// so a subset can be checked against exactly which elements survived.
    fn labeled_box() -> Result<Mesh3> {
        let mut mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let stdev = (0..mesh.point_count()).map(|i| i as f64).collect();
        let labels = (0..mesh.face_count() as u32).collect();
        mesh.set_point_stdev(Some(stdev))?;
        mesh.set_face_labels(Some(labels))?;
        Ok(mesh)
    }

    #[test]
    fn create_from_mask_carries_the_attributes() -> Result<()> {
        let mesh = labeled_box()?;
        let mask = IndexMask::try_from_indices(&[0, 1], mesh.face_count())?;

        let sub = mesh.create_from_mask(&mask)?;

        assert_eq!(sub.face_count(), 2);
        assert_eq!(sub.face_labels().unwrap(), &[0, 1]);

        // The attribute arrays have to match the counts of the mesh they ended up on.
        assert_eq!(sub.point_stdev().unwrap().len(), sub.point_count());
        sub.attrs().validate(sub.point_count(), sub.face_count())?;

        // The surviving standard deviations are exactly those of the points the faces reference.
        let kept: Vec<f64> = mesh
            .unique_point_mask(&mask)?
            .iter_true()
            .map(|i| mesh.point_stdev().unwrap()[i])
            .collect();
        assert_eq!(sub.point_stdev().unwrap(), kept.as_slice());

        Ok(())
    }

    #[test]
    fn create_from_indices_carries_the_attributes() -> Result<()> {
        let mesh = labeled_box()?;
        let sub = mesh.create_from_indices(&[2, 3])?;

        assert_eq!(sub.face_count(), 2);
        assert_eq!(sub.face_labels().unwrap(), &[2, 3]);
        sub.attrs().validate(sub.point_count(), sub.face_count())?;

        Ok(())
    }

    /// Routing through a mask makes the selection a set, so order does not matter and a repeat
    /// selects its face once rather than duplicating it.
    #[test]
    fn create_from_indices_normalizes_the_selection() -> Result<()> {
        let mesh = labeled_box()?;

        let ordered = mesh.create_from_indices(&[1, 3])?;
        let reversed = mesh.create_from_indices(&[3, 1])?;
        let repeated = mesh.create_from_indices(&[3, 1, 3])?;

        assert_eq!(ordered.faces(), reversed.faces());
        assert_eq!(ordered.faces(), repeated.faces());
        assert_eq!(ordered.face_labels().unwrap(), &[1, 3]);
        assert_eq!(repeated.face_count(), 2);

        Ok(())
    }

    #[test]
    fn create_from_indices_rejects_an_out_of_range_index() -> Result<()> {
        let mesh = labeled_box()?;
        assert!(mesh.create_from_indices(&[0, mesh.face_count()]).is_err());
        Ok(())
    }

    /// The split re-triangulates and introduces new points, with no mapping back, so it has to
    /// refuse rather than silently drop what it cannot carry.
    #[test]
    fn split_refuses_an_attributed_mesh_unless_the_loss_is_accepted() -> Result<()> {
        let mesh = labeled_box()?;
        let plane = Plane3::yz();

        assert!(mesh.split(&plane, false).is_err());

        // The message has to name what would be lost.
        let message = match mesh.split(&plane, false) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("expected the split to be refused"),
        };
        assert!(message.contains("point_stdev"), "{message}");
        assert!(message.contains("face_labels"), "{message}");

        // Accepting the loss lets it through, and the halves come back bare.
        match mesh.split(&plane, true)? {
            SplitResult::Pair(a, b) => {
                assert!(a.attrs().is_empty());
                assert!(b.attrs().is_empty());
            }
            _ => panic!("expected the plane to cut the box in two"),
        }

        // A bare mesh needs no flag at all.
        let bare = Mesh3::create_box(1.0, 1.0, 1.0, false);
        assert!(bare.split(&plane, false).is_ok());

        Ok(())
    }

    #[test]
    fn convex_hull_refuses_an_attributed_mesh_unless_the_loss_is_accepted() -> Result<()> {
        let mesh = labeled_box()?;

        assert!(mesh.convex_hull(false).is_err());
        assert!(mesh.convex_hull(true)?.attrs().is_empty());

        let bare = Mesh3::create_box(1.0, 1.0, 1.0, false);
        assert!(bare.convex_hull(false).is_ok());

        Ok(())
    }

    // ===========================================================================================
    // Scaling
    // ===========================================================================================

    /// The two containers have to agree on what a negative factor means, or the same mesh scaled
    /// through each would come back with opposite orientations.
    #[test]
    fn a_negative_scale_reverses_the_winding_like_mesh_data_does() -> Result<()> {
        let mesh = Mesh3::from_data(attributed_data(), false)?;

        let through_mesh = mesh.scale_copy(-1.0)?;
        let through_data = mesh.to_data().scale_copy(-1.0)?;

        assert_eq!(through_mesh.faces(), through_data.faces());
        assert_eq!(
            through_mesh.point_normals().unwrap(),
            through_data.point_normals().unwrap()
        );

        // The single triangle was wound +z, so mirroring must leave it wound -z.
        let p = through_mesh.points();
        let f = through_mesh.faces()[0];
        let normal = (p[f[1] as usize] - p[f[0] as usize])
            .cross(&(p[f[2] as usize] - p[f[0] as usize]))
            .normalize();
        assert_relative_eq!(normal, -Vector3::z(), epsilon = 1.0e-12);

        Ok(())
    }

    #[test]
    fn scale_copy_rejects_zero_and_non_finite_factors() -> Result<()> {
        let mesh = Mesh3::from_data(attributed_data(), false)?;

        assert!(mesh.scale_copy(0.0).is_err());
        assert!(mesh.scale_copy(f64::NAN).is_err());
        assert!(mesh.scale_copy(f64::INFINITY).is_err());

        Ok(())
    }

    // ===========================================================================================
    // Serialization
    // ===========================================================================================

    /// A file on disk, removed when the test finishes with it.
    #[cfg(any(feature = "ply", feature = "stl"))]
    struct TempFile {
        path: std::path::PathBuf,
    }

    #[cfg(any(feature = "ply", feature = "stl"))]
    impl TempFile {
        fn new(name: &str, ext: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "engeom-mesh3-{}-{}.{}",
                name,
                std::process::id(),
                ext
            ));
            Self { path }
        }

        fn path(&self) -> &Path {
            &self.path
        }
    }

    #[cfg(any(feature = "ply", feature = "stl"))]
    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.path);
        }
    }

    /// The whole point of routing through `MeshData3` is that the attributes survive, which the
    /// deleted shims did not do.
    #[cfg(feature = "ply")]
    #[test]
    fn a_ply_round_trip_keeps_the_attributes_and_the_solid_flag() -> Result<()> {
        let before = Mesh3::from_data(attributed_data(), false)?;
        let file = TempFile::new("round-trip", "ply");

        before.save_ply(file.path(), &PlyWriteOpts::default())?;
        let after = Mesh3::load_ply(file.path(), true)?;

        assert!(after.is_solid());
        assert_eq!(after.points(), before.points());
        assert_eq!(after.faces(), before.faces());
        assert_eq!(after.point_stdev(), before.point_stdev());
        assert_eq!(after.face_labels(), before.face_labels());
        assert_eq!(
            after.point_attr("confidence"),
            before.point_attr("confidence")
        );

        Ok(())
    }

    /// A PLY point cloud has no faces, so there is nothing to accelerate.
    #[cfg(feature = "ply")]
    #[test]
    fn loading_a_point_cloud_ply_as_a_mesh_is_refused() -> Result<()> {
        use crate::geom3::PointCloudData3;

        let cloud = PointCloudData3::new(vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)]);
        let file = TempFile::new("is-a-cloud", "ply");
        cloud.save_ply(file.path(), &PlyWriteOpts::default())?;

        assert!(Mesh3::load_ply(file.path(), false).is_err());

        Ok(())
    }

    #[cfg(feature = "stl")]
    #[test]
    fn an_stl_round_trip_carries_the_geometry() -> Result<()> {
        let before = Mesh3::create_box(1.0, 2.0, 3.0, false);
        let file = TempFile::new("round-trip", "stl");

        before.save_stl(file.path(), &StlWriteOpts::default())?;
        let after = Mesh3::load_stl(file.path(), false, true, false)?;

        assert_eq!(after.face_count(), before.face_count());
        assert_relative_eq!(after.aabb().maxs, before.aabb().maxs, epsilon = 1.0e-6);

        Ok(())
    }

    /// STL carries no attributes, so writing an attributed mesh has to be refused rather than
    /// silently stripping them.
    #[cfg(feature = "stl")]
    #[test]
    fn saving_stl_refuses_to_drop_attributes_silently() -> Result<()> {
        let mesh = Mesh3::from_data(attributed_data(), false)?;
        let file = TempFile::new("lossy", "stl");

        assert!(
            mesh.save_stl(file.path(), &StlWriteOpts::default())
                .is_err()
        );

        let mut opts = StlWriteOpts::default();
        opts.allow_attribute_loss = true;
        mesh.save_stl(file.path(), &opts)?;

        Ok(())
    }

    // ===========================================================================================
    // Attribute setters
    // ===========================================================================================

    #[test]
    fn attribute_setters_supply_the_counts() -> Result<()> {
        let mut mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let n_points = mesh.point_count();
        let n_faces = mesh.face_count();

        mesh.set_point_stdev(Some(vec![0.5; n_points]))?;
        mesh.set_face_labels(Some(vec![3; n_faces]))?;

        assert_eq!(mesh.point_stdev().unwrap().len(), n_points);
        assert_eq!(mesh.face_labels().unwrap().len(), n_faces);

        // The wrong length is rejected without the caller having to know the count.
        assert!(mesh.set_point_stdev(Some(vec![0.5; n_points + 1])).is_err());
        assert!(mesh.set_face_labels(Some(vec![3; n_faces - 1])).is_err());

        // The rejected calls must have left the existing attributes in place.
        assert_eq!(mesh.point_stdev().unwrap().len(), n_points);

        Ok(())
    }

    #[test]
    fn set_attrs_validates_against_the_current_counts() -> Result<()> {
        let mut mesh = Mesh3::from_data(attributed_data(), false)?;

        let mut bad = MeshAttrSet3::empty();
        bad.set_face_labels(Some(vec![1, 2, 3]), 3)?;
        assert!(mesh.set_attrs(bad).is_err());

        // The rejected set must have left the existing attributes in place.
        assert_eq!(mesh.face_labels().unwrap(), &[7]);

        let taken = mesh.take_attrs();
        assert_eq!(taken.face_labels().unwrap(), &[7]);
        assert!(mesh.attrs().is_empty());

        Ok(())
    }
}
