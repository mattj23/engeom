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
pub mod patches;
mod queries;
pub mod sampling;
mod section;
mod uv_mapping;

use crate::common::{IndexMask, PCoords};
use crate::na::SVector;
use crate::{Iso3, Point2, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
pub use collisions::MeshCollisionSet;
pub use data::{Attr3, MeshAttrSet3, MeshData3};
pub use edges::MeshEdges;
pub use half_edge::HalfEdgeMesh;
pub use nav_structure::MeshNav;
use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::shape::{TriMesh, TriMeshFlags};
use parry3d_f64::transformation;
pub use patches::{PatchFilter, PatchLabels, PatchStats};
pub use uv_mapping::UvMapping;

#[cfg(feature = "ply")]
use crate::io::PlyWriteOpts;
#[cfg(feature = "stl")]
use crate::io::StlWriteOpts;
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
    pub fn from_trimesh(shape: TriMesh, is_solid: bool) -> Self {
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

    /// Load a triangle mesh from a GOM `.g3d` file, the format written by GOM's Atos scanner and
    /// GOM/Zeiss Inspect. See [`MeshData3::load_g3d`] for what this does and does not support.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the `.g3d` file
    /// * `is_solid`: whether distance queries should treat points inside the mesh as being at
    ///   zero distance
    ///
    /// returns: `Result<Mesh3>`
    pub fn load_g3d(path: &Path, is_solid: bool) -> Result<Self> {
        Self::from_data(MeshData3::load_g3d(path)?, is_solid)
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
    pub fn transform_in_place(&mut self, transform: &Iso3) {
        self.shape.transform_vertices(transform);
        self.attrs.transform_in_place(transform);
    }

    /// Returns a new mesh with all vertices transformed by the given isometry, leaving the
    /// original unchanged.
    pub fn transform_copy(&self, transform: &Iso3) -> Self {
        let mut result = self.clone();
        result.transform_in_place(transform);
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

        let mut result = Mesh3::from_trimesh(new_shape, self.is_solid);
        result.attrs = self.attrs.clone();
        result.attrs.scale_in_place(scale);

        if scale < 0.0 {
            result.flip_normals_in_place();
        }

        Ok(result)
    }

    /// Create a new mesh by offsetting each point along its smoothed point normal.
    ///
    /// The offset is applied as `point + offset * normal`, where the normal is the per-point normal
    /// computed by `compute_point_normals`.
    ///
    /// Positive offsets expand the mesh outward; negative offsets shrink it inward.
    /// The original mesh is not modified.
    ///
    /// # Arguments
    ///
    /// * `offset`: The distance to offset each point along its normal.
    ///
    /// returns: `Result<Mesh3>`, failing if any point has no well-defined normal
    pub fn offset_points_copy(&self, offset: f64) -> Result<Self> {
        let normals = self.compute_point_normals()?;

        let updated = self
            .points()
            .iter()
            .zip(normals.iter())
            .map(|(v, n)| v + offset * n.into_inner())
            .collect();

        Ok(Self::new(updated, self.faces().to_vec(), self.is_solid))
    }

    /// Reverse the winding order of every face, turning the surface inside out.
    ///
    /// Any stored point normals are negated to match, since the direction the surface faces has
    /// changed.
    pub fn flip_normals_in_place(&mut self) {
        self.shape.reverse();
        self.attrs.flip_in_place();
    }
}

// ===============================================================================================
// Unsorted
// ===============================================================================================

impl Mesh3 {
    pub fn compute_edges(&self) -> Result<MeshEdges<'_>> {
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
    pub fn append_in_place(&mut self, other: &Mesh3) -> Result<()> {
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
        let mp = self.surface_closest_to(p);
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
    pub fn compute_nav(&self) -> MeshNav<'_> {
        MeshNav::new(self)
    }

    /// Label every face with the connected patch it belongs to.
    ///
    /// The labeling costs one `u32` per face regardless of how many patches the mesh splits into.
    /// Use `PatchLabels::mask_where` or `mask_of` to turn the result into face masks, one at a
    /// time, rather than materializing all of them at once. See `MeshNav::patch_labels` for the
    /// connectivity rule.
    ///
    /// This builds a `MeshNav` on every call. Build one yourself with `compute_nav()` if you are
    /// going to query the structure of the mesh more than once.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask restricting which faces take part, as if the mesh had been
    ///   pruned to it.
    ///
    /// returns: `Result<PatchLabels>`
    pub fn compute_patch_labels(&self, mask: Option<&IndexMask>) -> Result<PatchLabels> {
        let nav = self.compute_nav();
        nav.patch_labels(mask)
    }

    /// Build a face mask selecting the connected patches which pass a filter.
    ///
    /// Use this rather than `remove_small_patches` when you want to see what would be discarded,
    /// or to combine the selection with other criteria before extracting.
    ///
    /// # Arguments
    ///
    /// * `filter`: which patches are worth keeping
    ///
    /// returns: `Result<IndexMask>` over the mesh's faces
    pub fn patch_mask(&self, filter: &PatchFilter) -> Result<IndexMask> {
        let nav = self.compute_nav();
        nav.patch_mask(filter, None)
    }

    /// Discard the connected patches which fail a filter, returning what is left.
    ///
    /// This is the tool for cleaning flying patches and speckle out of scan data. Every attribute
    /// survives, because dropping whole faces keeps an index mapping back to the original.
    ///
    /// If no patch fails the filter the mesh is returned unchanged, rather than rebuilt, so a
    /// filter which finds nothing to do is cheap and preserves the UV mapping.
    ///
    /// # Arguments
    ///
    /// * `filter`: which patches are worth keeping
    ///
    /// returns: `Result<Mesh3>`, failing if the filter would discard every face, since a mesh
    /// needs at least one face
    pub fn remove_small_patches(&self, filter: &PatchFilter) -> Result<Self> {
        let mask = self.patch_mask(filter)?;
        let kept = mask.count_true();

        if kept == 0 {
            return Err(
                "Every patch was discarded by the filter, which would leave an empty mesh".into(),
            );
        }

        if kept == self.faces().len() {
            return Ok(self.clone());
        }

        self.extract_subset_faces(&mask)
    }

    /// Gets the boundary points of each patch in the mesh.  This function will return a list of
    /// lists of points, where each list of points is the boundary of a patch.  Note that this
    /// function will not work on non-manifold meshes.
    ///
    /// returns: Result<Vec<Vec<usize, Global>, Global>>
    pub fn compute_patch_boundary_points(&self) -> Result<Vec<Vec<Point3>>> {
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

    /// Compute the unit normal of every face, from the winding of its three points.
    ///
    /// returns: `Result<Vec<UnitVec3>>`, one normal per face, failing on a face which has no
    /// well-defined normal because its points are coincident or collinear
    pub fn compute_face_normals(&self) -> Result<Vec<UnitVec3>> {
        algorithms::compute_face_normals(self.points(), self.faces())
    }

    /// Compute the area of every face.
    ///
    /// returns: `Result<Vec<f64>>`, one area per face. A degenerate face has an area of zero rather
    /// than being an error.
    pub fn compute_face_areas(&self) -> Result<Vec<f64>> {
        algorithms::compute_face_areas(self.points(), self.faces())
    }

    /// Compute the centroid of every face, which is the mean of its three points.
    ///
    /// returns: `Result<Vec<Point3>>`, one centroid per face
    pub fn compute_face_centers(&self) -> Result<Vec<Point3>> {
        algorithms::compute_face_centers(self.points(), self.faces())
    }

    /// Compute a smoothed normal for every point by averaging the normals of the faces which touch
    /// it, weighting each by the interior angle of that face at that point.
    ///
    /// See `algorithms::compute_point_normals` for why the weighting is by angle and not by area,
    /// and for what this quantity is and is not good for.
    ///
    /// returns: `Result<Vec<UnitVec3>>`, one normal per point, failing if any point has no
    /// well-defined normal, which happens when it belongs to no face, belongs only to degenerate
    /// faces, or sits where the weighted contributions cancel exactly
    pub fn compute_point_normals(&self) -> Result<Vec<UnitVec3>> {
        algorithms::compute_point_normals(self.points(), self.faces())
    }
}

// ===============================================================================================
// Shape creation methods
// ===============================================================================================
//
// Every one of these is a pass-through to the `MeshData3` constructor of the same name, adding the
// `is_solid` flag. The tessellations live there because they produce points and faces and nothing
// else; `is_solid` is a query behavior and belongs only to the accelerated type.

impl Mesh3 {
    /// Create a box mesh with the given dimensions, centered at the origin.
    ///
    /// See `MeshData3::create_box`.
    ///
    /// # Arguments
    ///
    /// * `length`: the dimension of the box in the x direction
    /// * `width`: the dimension of the box in the y direction
    /// * `height`: the dimension of the box in the z direction
    /// * `is_solid`: whether distance queries should treat points inside the mesh as being at zero
    ///   distance
    ///
    /// returns: `Mesh3`
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
        Self::from_primitive(MeshData3::create_box(length, width, height), is_solid)
    }

    /// Create a spherical mesh centered at the origin.
    ///
    /// See `MeshData3::create_sphere`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the sphere
    /// * `n_theta`: number of subdivisions around the polar direction
    /// * `n_phi`: number of subdivisions around the azimuthal direction
    ///
    /// returns: `Mesh3`
    pub fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        Self::from_primitive(MeshData3::create_sphere(radius, n_theta, n_phi), true)
    }

    /// Create a cylindrical mesh centered at the origin and aligned with the local `y` axis.
    ///
    /// See `MeshData3::create_cylinder`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the cylinder
    /// * `height`: full height of the cylinder, along the y axis
    /// * `steps`: number of subdivisions around the circumference
    ///
    /// returns: `Mesh3`
    pub fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        Self::from_primitive(MeshData3::create_cylinder(radius, height, steps), true)
    }

    /// Create a conical mesh centered at the origin and aligned with the local `y` axis, with its
    /// apex at `+height/2` and its base at `-height/2`.
    ///
    /// See `MeshData3::create_cone`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the base of the cone
    /// * `height`: full height of the cone, along the y axis
    /// * `steps`: number of subdivisions around the circumference
    ///
    /// returns: `Mesh3`
    pub fn create_cone(radius: f64, height: f64, steps: usize) -> Self {
        Self::from_primitive(MeshData3::create_cone(radius, height, steps), true)
    }

    /// Create a flat, filled circle mesh lying in the XY plane, centered at the origin, with the
    /// normal pointing along +Z.
    ///
    /// See `MeshData3::create_circle`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the circle
    /// * `segments`: number of perimeter points, and of triangles. Must be at least 3.
    ///
    /// returns: `Mesh3`
    pub fn create_circle(radius: f64, segments: usize) -> Self {
        Self::from_primitive(MeshData3::create_circle(radius, segments), false)
    }

    /// Create a capsule mesh spanning the segment between two points.
    ///
    /// See `MeshData3::create_capsule`.
    pub fn create_capsule(
        p0: &Point3,
        p1: &Point3,
        radius: f64,
        n_theta: usize,
        n_phi: usize,
    ) -> Self {
        Self::from_primitive(
            MeshData3::create_capsule(p0, p1, radius, n_theta, n_phi),
            true,
        )
    }

    /// Create a cylindrical mesh spanning the segment between two points.
    ///
    /// See `MeshData3::create_cylinder_between`.
    pub fn create_cylinder_between(p0: &Point3, p1: &Point3, radius: f64, steps: usize) -> Self {
        Self::from_primitive(
            MeshData3::create_cylinder_between(p0, p1, radius, steps),
            true,
        )
    }

    /// Create a rectangular beam spanning the segment between two points.
    ///
    /// See `MeshData3::create_rect_beam_between`.
    ///
    /// returns: `Result<Mesh3>`, failing if `up` is parallel to the segment
    pub fn create_rect_beam_between(
        p0: &Point3,
        p1: &Point3,
        width: f64,
        height: f64,
        up: &Vector3,
    ) -> Result<Self> {
        let data = MeshData3::create_rect_beam_between(p0, p1, width, height, up)?;
        Ok(Self::from_primitive(data, true))
    }

    /// Load a Stanford bunny mesh embedded in the binary with 453 points and 948 faces.
    ///
    /// See `MeshData3::stanford_bunny_res4`.
    pub fn stanford_bunny_res4() -> Self {
        Self::from_primitive(MeshData3::stanford_bunny_res4(), false)
    }

    /// Load a Stanford bunny mesh embedded in the binary with 1889 points and 3851 faces.
    ///
    /// See `MeshData3::stanford_bunny_res3`.
    pub fn stanford_bunny_res3() -> Self {
        Self::from_primitive(MeshData3::stanford_bunny_res3(), false)
    }

    /// Load a Stanford bunny mesh embedded in the binary with 8171 points and 16301 faces.
    ///
    /// See `MeshData3::stanford_bunny_res2`.
    pub fn stanford_bunny_res2() -> Self {
        Self::from_primitive(MeshData3::stanford_bunny_res2(), false)
    }

    /// Accelerate a mesh which came from one of the primitive constructors.
    ///
    /// Those always produce at least one face, which is the only way `from_data` can fail, so the
    /// failure is a bug in the tessellation rather than anything a caller can cause.
    fn from_primitive(data: MeshData3, is_solid: bool) -> Self {
        Self::from_data(data, is_solid)
            .expect("A primitive tessellation always has at least one face")
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
    fn point_normals_match_point_count() -> Result<()> {
        let mesh = Mesh3::create_sphere(1.0, 40, 40);
        let normals = mesh.compute_point_normals()?;

        assert_eq!(normals.len(), mesh.points().len());

        Ok(())
    }

    /// The stanford bunny asset carries a defect: faces 207 and 209 are the same triangle wound in
    /// opposite directions, and point 128 belongs to nothing else. Its two incident normals cancel,
    /// so it genuinely has no normal, and the computation is expected to say so rather than
    /// normalizing the floating point residue into a direction, which is what it used to do.
    ///
    /// The indices moved when the fixture was re-encoded from micro-mesh to tcmesh, which renumbers
    /// vertices. The defect itself is unchanged: the same pair of opposed triangles, the same lone
    /// point, verified rather than assumed when the numbers were updated.
    #[test]
    fn point_normals_report_the_bunny_asset_defect() {
        let mesh = stanford_bun_4();

        let err = mesh.compute_point_normals().unwrap_err();

        assert!(
            err.to_string().starts_with("Point 128 has no well-defined"),
            "unexpected error: {err}"
        );
    }

    /// The outline uses point normals only to nudge itself off the surface, so it tolerates the
    /// defect above rather than failing on a mesh it can perfectly well draw.
    #[test]
    fn visual_outline_survives_a_point_with_no_normal() -> Result<()> {
        let mesh = stanford_bun_4();

        let outline =
            mesh.compute_visual_outline(UnitVec3::new_normalize(Vector3::z()), 1.0, None)?;

        assert!(!outline.is_empty());
        Ok(())
    }

    /// The accelerated type and the plain container are two representations of the same buffers, so
    /// the normals they compute have to agree exactly rather than merely being close.
    #[test]
    fn point_normals_agree_with_mesh_data() -> Result<()> {
        let mesh = Mesh3::create_sphere(1.0, 40, 40);
        let from_mesh = mesh.compute_point_normals()?;
        let from_data = mesh.to_data().compute_point_normals()?;

        assert_eq!(from_mesh.len(), from_data.len());
        for (a, b) in from_mesh.iter().zip(from_data.iter()) {
            assert_relative_eq!(a.into_inner(), b.into_inner(), epsilon = 1.0e-15);
        }

        Ok(())
    }

    /// Angle weighting is invariant to how a flat face is triangulated, so every corner normal of a
    /// box lands exactly on the diagonal. Area weighting, which this method used to do, does not:
    /// a corner touching one triangle on one face and two on another gets pulled off the diagonal.
    #[test]
    fn point_normals_of_a_box_point_along_the_diagonals() -> Result<()> {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, true);
        let normals = mesh.compute_point_normals()?;
        let diagonal = 1.0 / 3.0_f64.sqrt();

        for (point, normal) in mesh.points().iter().zip(normals.iter()) {
            let expected = Vector3::new(
                point.x.signum() * diagonal,
                point.y.signum() * diagonal,
                point.z.signum() * diagonal,
            );
            assert_relative_eq!(normal.into_inner(), &expected, epsilon = 1.0e-12);
        }

        Ok(())
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
    fn offset_points_copy_preserves_spherical_radius() -> Result<()> {
        let radius = 1.0;
        let offset = 0.1;
        let mesh = Mesh3::create_sphere(radius, 100, 100);
        let offset_mesh = mesh.offset_points_copy(offset)?;

        assert_eq!(mesh.points().len(), offset_mesh.points().len());

        for vertex in offset_mesh.points() {
            assert_relative_eq!(vertex.coords.norm(), radius + offset, epsilon = 1.0e-5);
        }

        Ok(())
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
        mesh.transform_in_place(&iso);

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
        mesh.flip_normals_in_place();

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

        mesh.append_in_place(&other)?;

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

        assert!(mesh.append_in_place(&other).is_err());

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

        mesh.append_in_place(&other)?;
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
    fn extract_subset_faces_carries_the_attributes() -> Result<()> {
        let mesh = labeled_box()?;
        let mask = IndexMask::try_from_indices(&[0, 1], mesh.face_count())?;

        let sub = mesh.extract_subset_faces(&mask)?;

        assert_eq!(sub.face_count(), 2);
        assert_eq!(sub.face_labels().unwrap(), &[0, 1]);

        // The attribute arrays have to match the counts of the mesh they ended up on.
        assert_eq!(sub.point_stdev().unwrap().len(), sub.point_count());
        sub.attrs().validate(sub.point_count(), sub.face_count())?;

        // The surviving standard deviations are exactly those of the points the faces reference.
        let kept: Vec<f64> = mesh
            .compute_unique_point_mask(&mask)?
            .iter_true()
            .map(|i| mesh.point_stdev().unwrap()[i])
            .collect();
        assert_eq!(sub.point_stdev().unwrap(), kept.as_slice());

        Ok(())
    }

    #[test]
    fn extract_subset_faces_from_indices_carries_the_attributes() -> Result<()> {
        let mesh = labeled_box()?;
        let sub = mesh.extract_subset_faces_from_indices(&[2, 3])?;

        assert_eq!(sub.face_count(), 2);
        assert_eq!(sub.face_labels().unwrap(), &[2, 3]);
        sub.attrs().validate(sub.point_count(), sub.face_count())?;

        Ok(())
    }

    /// Routing through a mask makes the selection a set, so order does not matter and a repeat
    /// selects its face once rather than duplicating it.
    #[test]
    fn extract_subset_faces_from_indices_normalizes_the_selection() -> Result<()> {
        let mesh = labeled_box()?;

        let ordered = mesh.extract_subset_faces_from_indices(&[1, 3])?;
        let reversed = mesh.extract_subset_faces_from_indices(&[3, 1])?;
        let repeated = mesh.extract_subset_faces_from_indices(&[3, 1, 3])?;

        assert_eq!(ordered.faces(), reversed.faces());
        assert_eq!(ordered.faces(), repeated.faces());
        assert_eq!(ordered.face_labels().unwrap(), &[1, 3]);
        assert_eq!(repeated.face_count(), 2);

        Ok(())
    }

    #[test]
    fn extract_subset_faces_from_indices_rejects_an_out_of_range_index() -> Result<()> {
        let mesh = labeled_box()?;
        assert!(
            mesh.extract_subset_faces_from_indices(&[0, mesh.face_count()])
                .is_err()
        );
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

        let opts = StlWriteOpts {
            allow_attribute_loss: true,
            ..Default::default()
        };
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
