//! This module contains `PointCloud3`, a plain container for a buffer of points and the
//! per-point attributes attached to them.
//!
//! It is the point-cloud counterpart to `MeshData3`, and the two are deliberately shaped the same
//! way: private buffers, a validated attribute set, and the same construction, subsetting, and
//! serialization vocabulary. The attribute half is literally the same code, `PointAttrSet3`, which
//! `MeshAttrSet3` also composes.
//!
//! There is no spatial acceleration here of any kind. Call `compute_index` for a `CloudIndex3`
//! when you need nearest-neighbor queries.

use crate::common::IndexMask;
use crate::common::poisson_disk::sample_poisson_disk_all;
use crate::geom3::Aabb3;
use crate::geom3::attributes3::{Attr3, PointAttrSet3};
use crate::geom3::point_cloud::CloudIndex3;
use crate::{Iso3, KdTree3, Point3, Result, SurfacePoint3, UnitVec3};
use std::fmt;

#[cfg(feature = "ply")]
use crate::io::{PlyWriteOpts, load_ply_points, write_ply_points};
#[cfg(feature = "ply")]
use std::path::Path;

/// A container for the raw data of a point cloud: a buffer of points and the per-point attributes
/// attached to them.
///
/// This type performs no spatial acceleration of any kind. It is cheap to construct, cheap to edit,
/// and is the form in which point data is read from and written to files.
///
/// # Relationship with Other Representations
///
/// `PointCloud3` is the serialization gateway for point data, and is the default choice when:
///
/// - You don't need any spatial queries at all
/// - You are editing the contents of the buffers
/// - You are doing something custom with serialization or deserialization
/// - You are carrying open-map attributes alongside the typed ones
///
/// For nearest-neighbor queries, overlap checks, and Poisson sampling, call `compute_index` to get a
/// `CloudIndex3` borrowing this cloud. The index holds the cloud immutably for its lifetime, so
/// nothing can move the points out from under its tree.
///
/// # Invariants
///
/// A `PointCloud3` is never allowed to exist in an incoherent state. Every attribute array has
/// a length matching the point count. This is checked on construction and maintained by every
/// method which modifies the cloud, which is why the point buffer is private and why attributes are
/// set through methods that supply the count on the caller's behalf.
///
/// An empty cloud is legal.
#[derive(Clone, PartialEq)]
pub struct PointCloud3 {
    points: Vec<Point3>,
    attrs: PointAttrSet3,
}

// ===============================================================================================
// Construction
// ===============================================================================================

impl PointCloud3 {
    /// Create a new point cloud from a buffer of points, with no attributes attached.
    ///
    /// Unlike `MeshData3::new` this cannot fail, because a bare point buffer has no internal
    /// consistency to violate.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    ///
    /// returns: `PointCloud3`
    pub fn new(points: Vec<Point3>) -> Self {
        Self {
            points,
            attrs: PointAttrSet3::empty(),
        }
    }

    /// Create a new point cloud from a buffer of points and a set of per-point attributes.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    /// * `attrs`: the attributes to attach, whose arrays must match the point count
    ///
    /// returns: `Result<PointCloud3>`, failing if any attribute array is the wrong length
    pub fn new_with_attrs(points: Vec<Point3>, attrs: PointAttrSet3) -> Result<Self> {
        attrs.validate(points.len())?;
        Ok(Self { points, attrs })
    }

    /// Create an empty point cloud, with no points and no attributes.
    pub fn empty() -> Self {
        Self::new(Vec::new())
    }

    /// Build a cloud from surface points, keeping their normals.
    ///
    /// returns: `PointCloud3`
    pub fn from_surface_points(points: &[SurfacePoint3]) -> Self {
        let normals: Vec<UnitVec3> = points.iter().map(|p| p.normal).collect();
        let mut attrs = PointAttrSet3::empty();
        attrs
            .set_normals(Some(normals), points.len())
            .expect("one normal per point by construction");

        Self {
            points: points.iter().map(|p| p.point).collect(),
            attrs,
        }
    }

    /// Build a k-d tree over this cloud and return it as a borrowed index.
    ///
    /// The index holds this cloud immutably for its lifetime, which is what makes it impossible for
    /// the tree to go stale. Building costs roughly what `0.4 * N` nearest-neighbor queries do (see
    /// `engeom/tests/kd_tree_backends.md`), so an index is cheap next to any real use of one and
    /// wasted if never queried, which is why it is requested rather than maintained automatically.
    ///
    /// returns: `Result<CloudIndex3>`, failing only if the tree could not be built
    pub fn compute_index(&self) -> Result<CloudIndex3<'_>> {
        CloudIndex3::try_new(self)
    }

    /// Pair this cloud with a k-d tree the caller has already built over it, skipping the build.
    ///
    /// # This is unchecked, and getting it wrong is silent
    ///
    /// **Nothing verifies that `tree` was built from this cloud.** If it was not, every query
    /// returns indices into a point set that is not this one. Those indices are in range whenever
    /// the other cloud was at least as large, so there is no panic and no error: you get plausible,
    /// confidently wrong answers, and a downstream alignment or overlap check quietly uses them.
    /// A tree built over a *smaller* cloud is worse only in that it eventually panics on a
    /// subscript, which at least tells you something is broken.
    ///
    /// The obligation is that `tree` was built from `self.points()` **as they are now**. A tree is
    /// invalidated by anything which adds, removes, reorders or moves a point:
    /// `transform_in_place`, `scale_in_place`, `append_in_place`, `points_mut`, and any subset
    /// extraction, which produces a different cloud entirely. Attribute changes do not invalidate
    /// it, since the tree indexes positions only.
    ///
    /// # Use [`PointCloud3::compute_index`] instead
    ///
    /// The checked path builds the tree itself, so it cannot be mispaired, and the borrow checker
    /// then prevents the cloud from being mutated while the index lives. It is also nearly free
    /// relative to using an index at all: a build costs about what `0.4 * N` nearest-neighbor
    /// queries do (`engeom/tests/kd_tree_backends.md`), and any real consumer does more than that.
    ///
    /// This exists for one situation: a long-lived wrapper which caches a tree across many calls and
    /// drops it on every mutation, where the cache and the invalidation sit together in one place
    /// and can be read side by side. The Python bindings are that wrapper. If you are reaching for
    /// this anywhere else, you almost certainly want `compute_index`.
    ///
    /// This is the same hazard as the `tree_uuid` token that the `CloudIndex3` design replaced, kept
    /// deliberately and narrowly rather than by accident. That token at least checked at runtime;
    /// this does not check at all.
    pub fn index_with_tree_unchecked<'a>(&'a self, tree: &'a KdTree3) -> CloudIndex3<'a> {
        CloudIndex3::with_tree_unchecked(self, tree)
    }

    /// Poisson disk sample the cloud, returning a mask selecting a subset of its points no two of
    /// which are closer together than `radius`.
    ///
    /// This lives on the cloud rather than on [`CloudIndex3`] because it does not use an index. The
    /// sampler voxel-downsamples first and builds its own tree over the survivors, which is a much
    /// smaller tree than one over the whole cloud, so handing it a prebuilt full-cloud index would
    /// waste the build rather than save it.
    pub fn sample_poisson_disk(&self, radius: f64) -> IndexMask {
        sample_poisson_disk_all(self.points(), radius)
    }

    /// Poisson disk sample the cloud and extract the result as a new cloud, carrying attributes.
    ///
    /// Equivalent to [`PointCloud3::extract_subset_points`] with the mask from
    /// [`PointCloud3::sample_poisson_disk`].
    pub fn extract_poisson_sample(&self, radius: f64) -> Result<Self> {
        let mask = self.sample_poisson_disk(radius);
        self.extract_subset_points(&mask)
    }

    /// Consume the cloud and return ownership of its two components: the point buffer and the
    /// attribute set.
    ///
    /// This is the counterpart to `new_with_attrs` and exists so that handing this data to another
    /// representation does not require copying the buffer. Once the cloud is decomposed nothing
    /// enforces the invariant between the two pieces any more, so a caller putting them back
    /// together is responsible for keeping them consistent.
    ///
    /// returns: `(Vec<Point3>, PointAttrSet3)`
    pub fn into_parts(self) -> (Vec<Point3>, PointAttrSet3) {
        (self.points, self.attrs)
    }
}

// ===============================================================================================
// Serialization
// ===============================================================================================

impl PointCloud3 {
    /// Load a point cloud from a PLY file, preserving every property the file carries.
    ///
    /// The file must have a `vertex` element with `x`, `y`, and `z` properties. A file which
    /// declares faces is **refused**, because it is a mesh and loading it here would silently throw
    /// the connectivity away. Use `MeshData3::load_ply` for those.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the PLY file
    ///
    /// returns: `Result<PointCloud3>`
    #[cfg(feature = "ply")]
    pub fn load_ply(path: &Path) -> Result<Self> {
        let (points, attrs, has_faces) = load_ply_points(path)?;

        if has_faces {
            return Err(format!(
                "The PLY file at {} declares faces, so it is a mesh rather than a point cloud. \
                 Loading it here would discard the connectivity; use MeshData3::load_ply instead.",
                path.display()
            )
            .into());
        }

        Ok(Self { points, attrs })
    }

    /// Write this point cloud to a PLY file, preserving every attribute it carries.
    ///
    /// If you use the default options, the data will be saved in binary format using double
    /// floating point precision. If you wish to alter this, construct the `PlyWriteOpts` directly.
    ///
    /// The `vertex` element written here is identical to the one `MeshData3::save_ply` writes for
    /// the same points; the file simply has no `face` element.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `opts`: encoding and header options
    ///
    /// returns: `Result<()>`
    #[cfg(feature = "ply")]
    pub fn save_ply(&self, path: &Path, opts: &PlyWriteOpts) -> Result<()> {
        write_ply_points(path, &self.points, &self.attrs, opts)
    }

    /// Verify that the caller has accepted the loss of this cloud's attributes, for a format which
    /// cannot represent them.
    ///
    /// A writer for a geometry-only format calls this before doing any work. If the cloud carries
    /// no attributes there is nothing to lose and this always succeeds, so the flag only ever
    /// matters when data would actually die. See `MeshData3::check_attribute_loss` for why this is
    /// an error rather than a warning.
    ///
    /// # Arguments
    ///
    /// * `format`: the name of the target format, used in the error message
    /// * `allow_loss`: whether the caller has accepted the loss
    ///
    /// returns: `Result<()>`
    pub fn check_attribute_loss(&self, format: &str, allow_loss: bool) -> Result<()> {
        if allow_loss || self.attrs.is_empty() {
            return Ok(());
        }

        Err(format!(
            "Writing to {} would discard the attributes on this point cloud ({}), because the \
             format cannot represent them. Set `allow_attribute_loss` to accept this.",
            format,
            self.attrs.attr_labels().join(", ")
        )
        .into())
    }
}

// ===============================================================================================
// Core access
// ===============================================================================================

impl PointCloud3 {
    /// Get a reference to the points of the cloud.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Get a mutable reference to the points of the cloud.
    ///
    /// The length cannot change through this, so the attribute arrays stay valid. Note that moving
    /// a point does **not** update any stored normals, which describe the surface the point was
    /// measured on.
    pub fn points_mut(&mut self) -> &mut [Point3] {
        &mut self.points
    }

    /// Get the number of points in the cloud.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// Returns true if the cloud has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Get a reference to the full set of per-point attributes attached to this cloud.
    pub fn attrs(&self) -> &PointAttrSet3 {
        &self.attrs
    }

    /// Compute the axis-aligned bounding box of the points.
    ///
    /// This is a full pass over the buffer with no caching, so hold the result if you need it more
    /// than once.
    pub fn compute_aabb(&self) -> Aabb3 {
        Aabb3::from_points(self.points.clone())
    }
}

// ===============================================================================================
// Attribute access
// ===============================================================================================

impl PointCloud3 {
    /// Get the per-point unit normals, if present.
    pub fn point_normals(&self) -> Option<&[UnitVec3]> {
        self.attrs.normals()
    }

    /// Get the per-point RGB colors, if present.
    pub fn point_colors(&self) -> Option<&[[u8; 3]]> {
        self.attrs.colors()
    }

    /// Get the per-point standard deviations, if present. These are 1-sigma values in the cloud's
    /// own length units.
    pub fn point_stdev(&self) -> Option<&[f64]> {
        self.attrs.stdev()
    }

    /// Get the open-map per-point attribute stored under the given name, if present.
    pub fn point_attr(&self, name: &str) -> Option<&Attr3> {
        self.attrs.attr(name)
    }
}

// ===============================================================================================
// Attribute mutation
// ===============================================================================================

impl PointCloud3 {
    /// Set or clear the per-point unit normals.
    ///
    /// # Arguments
    ///
    /// * `values`: the normals to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_normals(&mut self, values: Option<Vec<UnitVec3>>) -> Result<()> {
        self.attrs.set_normals(values, self.points.len())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them. Must match the point count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_colors(&mut self, values: Option<Vec<[u8; 3]>>) -> Result<()> {
        self.attrs.set_colors(values, self.points.len())
    }

    /// Set or clear the per-point standard deviations, which must be 1-sigma values in the cloud's
    /// own length units.
    ///
    /// # Arguments
    ///
    /// * `values`: the standard deviations to store, or `None` to clear them. Must match the point
    ///   count, and must be finite and non-negative.
    ///
    /// returns: `Result<()>`
    pub fn set_point_stdev(&mut self, values: Option<Vec<f64>>) -> Result<()> {
        self.attrs.set_stdev(values, self.points.len())
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
        self.attrs.insert_attr(name, attr, self.points.len())
    }

    /// Remove and return the open-map per-point attribute stored under the given name.
    pub fn remove_point_attr(&mut self, name: &str) -> Option<Attr3> {
        self.attrs.remove_attr(name)
    }

    /// Replace the entire set of per-point attributes.
    ///
    /// # Arguments
    ///
    /// * `attrs`: the attribute set to attach, whose arrays must match the point count
    ///
    /// returns: `Result<()>`, leaving the existing attributes untouched on failure
    pub fn set_attrs(&mut self, attrs: PointAttrSet3) -> Result<()> {
        attrs.validate(self.points.len())?;
        self.attrs = attrs;
        Ok(())
    }

    /// Remove and return the entire set of per-point attributes, leaving the cloud with none.
    pub fn take_attrs(&mut self) -> PointAttrSet3 {
        std::mem::take(&mut self.attrs)
    }
}

// ===============================================================================================
// Whole-cloud operations
// ===============================================================================================

impl PointCloud3 {
    /// Transform the cloud in place by a rigid isometry.
    ///
    /// The points are moved, and any stored normals and `Vector` attributes are rotated with them.
    /// Because those hold directions rather than positions, only the rotation component of the
    /// isometry affects them.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        for p in self.points.iter_mut() {
            *p = iso * *p;
        }

        self.attrs.transform_in_place(iso);
    }

    /// Return a copy of the cloud transformed by a rigid isometry, leaving this one unchanged.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    ///
    /// returns: `PointCloud3`
    pub fn transform_copy(&self, iso: &Iso3) -> Self {
        let mut result = self.clone();
        result.transform_in_place(iso);
        result
    }

    /// Scale the cloud in place about the origin by a uniform factor.
    ///
    /// Any stored standard deviations are scaled with the geometry, since they are lengths in the
    /// cloud's own units and would otherwise silently come to mean something else. Normals are
    /// unaffected by the magnitude of a uniform scale, being directions.
    ///
    /// A **negative** factor mirrors the cloud. Unlike a mesh there is no winding to reverse, but
    /// the stored normals are negated, because a reflection turns the surface they describe inside
    /// out.
    ///
    /// # Arguments
    ///
    /// * `scale`: the factor to scale by, which must be finite and non-zero
    ///
    /// returns: `Result<()>`
    pub fn scale_in_place(&mut self, scale: f64) -> Result<()> {
        if !scale.is_finite() || scale == 0.0 {
            return Err(format!(
                "A scale factor must be finite and non-zero, but {scale} was given. Scaling by \
                 zero would collapse the cloud to a point irrecoverably."
            )
            .into());
        }

        for p in self.points.iter_mut() {
            *p = Point3::from(p.coords * scale);
        }

        self.attrs.scale_in_place(scale);

        if scale < 0.0 {
            self.attrs.flip_in_place();
        }

        Ok(())
    }

    /// Return a copy of the cloud scaled about the origin by a uniform factor, leaving this one
    /// unchanged.
    ///
    /// # Arguments
    ///
    /// * `scale`: the factor to scale by, which must be finite and non-zero
    ///
    /// returns: `Result<PointCloud3>`
    pub fn scale_copy(&self, scale: f64) -> Result<Self> {
        let mut result = self.clone();
        result.scale_in_place(scale)?;
        Ok(result)
    }

    /// Append another cloud onto the end of this one.
    ///
    /// Attributes are all-or-nothing: a typed field or an open-map key present on one side and
    /// absent on the other is an error, because there is no correct value to pad the missing side
    /// with. The whole append is validated before anything is modified, so a failure leaves this
    /// cloud untouched.
    ///
    /// # Arguments
    ///
    /// * `other`: the cloud to append
    ///
    /// returns: `Result<()>`
    pub fn append_in_place(&mut self, other: &PointCloud3) -> Result<()> {
        // Attributes are checked and merged first, because that is the step which can fail. Once it
        // succeeds, extending the point buffer cannot.
        self.attrs.extend_from(&other.attrs)?;
        self.points.extend_from_slice(&other.points);
        Ok(())
    }
}

// ===============================================================================================
// Subsetting
// ===============================================================================================

impl PointCloud3 {
    /// Create a new cloud containing only the points marked `true` in the given mask, preserving
    /// their original order. Every attribute is carried through.
    ///
    /// # Arguments
    ///
    /// * `mask`: a mask whose length must match the point count
    ///
    /// returns: `Result<PointCloud3>`
    pub fn extract_subset_points(&self, mask: &IndexMask) -> Result<Self> {
        if mask.len() != self.points.len() {
            return Err(format!(
                "The mask has {} entries, but the cloud has {} points",
                mask.len(),
                self.points.len()
            )
            .into());
        }

        Ok(Self {
            points: mask.clone_indices_of(&self.points)?,
            attrs: self.attrs.subset(mask)?,
        })
    }

    /// Create a new cloud containing the points at the given indices, in the order the indices are
    /// given. Indices may repeat. Every attribute is carried through.
    ///
    /// # Arguments
    ///
    /// * `indices`: the points to take, each of which must be less than the point count
    ///
    /// returns: `Result<PointCloud3>`
    pub fn extract_subset_indices(&self, indices: &[usize]) -> Result<Self> {
        let n = self.points.len();
        if let Some(bad) = indices.iter().find(|&&i| i >= n) {
            return Err(format!("Index {bad} is out of bounds for a cloud with {n} points").into());
        }

        Ok(Self {
            points: indices.iter().map(|&i| self.points[i]).collect(),
            attrs: self.attrs.subset_indices(indices)?,
        })
    }
}

// ===============================================================================================
// Helpers
// ===============================================================================================

impl fmt::Debug for PointCloud3 {
    /// Summarize the cloud rather than dumping its buffer, which is routinely large enough to make
    /// a derived `Debug` implementation useless.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PointCloud3")
            .field("points", &self.points.len())
            .field("attrs", &self.attrs)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::common::kd_tree::KdTreeSearch;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    /// A four point cloud carrying one attribute of every kind an operation might have to touch.
    fn loaded_cloud() -> PointCloud3 {
        let mut cloud = PointCloud3::new(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ]);

        cloud
            .set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 4]))
            .unwrap();
        cloud
            .set_point_colors(Some(vec![[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]]))
            .unwrap();
        cloud
            .set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]))
            .unwrap();
        cloud
            .insert_point_attr("confidence", Attr3::Scalar(vec![0.5, 0.6, 0.7, 0.8]))
            .unwrap();
        cloud
            .insert_point_attr("principal_dir", Attr3::Vector(vec![Vector3::x(); 4]))
            .unwrap();

        cloud
    }

    #[test]
    fn a_bare_cloud_has_no_attributes() {
        let cloud = PointCloud3::new(vec![Point3::origin()]);

        assert_eq!(cloud.point_count(), 1);
        assert!(!cloud.is_empty());
        assert!(cloud.attrs().is_empty());
        assert!(cloud.point_normals().is_none());
    }

    #[test]
    fn an_empty_cloud_is_legal() {
        let cloud = PointCloud3::empty();
        assert!(cloud.is_empty());
        assert_eq!(cloud.point_count(), 0);
    }

    #[test]
    fn new_with_attrs_rejects_a_mismatched_attribute() {
        let mut attrs = PointAttrSet3::empty();
        attrs.set_stdev(Some(vec![0.1]), 1).unwrap();

        // The attribute set is internally consistent for a one-point cloud, but not for this one.
        assert!(PointCloud3::new_with_attrs(vec![Point3::origin(); 4], attrs).is_err());
    }

    #[test]
    fn attribute_setters_supply_the_count() -> Result<()> {
        let mut cloud = loaded_cloud();

        assert_eq!(cloud.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);
        assert_eq!(cloud.point_attr("confidence").unwrap().len(), 4);

        // The wrong length is rejected without the caller having to know the count.
        assert!(cloud.set_point_stdev(Some(vec![0.1, 0.2])).is_err());
        assert_eq!(cloud.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);

        Ok(())
    }

    #[test]
    fn compute_aabb_bounds_the_points() {
        let aabb = loaded_cloud().compute_aabb();

        assert_relative_eq!(aabb.mins, Point3::origin(), epsilon = 1.0e-12);
        assert_relative_eq!(aabb.maxs, Point3::new(1.0, 1.0, 1.0), epsilon = 1.0e-12);
    }

    #[test]
    fn debug_summarizes_instead_of_dumping_the_buffer() {
        let text = format!("{:?}", loaded_cloud());
        assert!(text.contains("points: 4"), "{text}");
    }

    // ===========================================================================================
    // Whole-cloud operations
    // ===========================================================================================

    #[test]
    fn transform_moves_points_and_rotates_directions() {
        let mut cloud = loaded_cloud();

        // A quarter turn about +z, which maps +x onto +y, plus a translation.
        let iso = Iso3::new(Vector3::new(10.0, 0.0, 0.0), Vector3::z() * FRAC_PI_2);
        cloud.transform_in_place(&iso);

        assert_relative_eq!(
            cloud.points()[1],
            Point3::new(10.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // +z is on the axis and does not move.
        assert_relative_eq!(
            cloud.point_normals().unwrap()[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            cloud
                .point_attr("principal_dir")
                .unwrap()
                .as_vector()
                .unwrap()[0],
            Vector3::y(),
            epsilon = 1.0e-12
        );

        // Scalars are untouched.
        assert_eq!(cloud.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);
    }

    #[test]
    fn transform_copy_leaves_the_original_alone() {
        let cloud = loaded_cloud();
        let moved = cloud.transform_copy(&Iso3::translation(5.0, 0.0, 0.0));

        assert_relative_eq!(cloud.points()[0], Point3::origin(), epsilon = 1.0e-12);
        assert_relative_eq!(
            moved.points()[0],
            Point3::new(5.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );
    }

    #[test]
    fn scale_scales_points_and_standard_deviations() -> Result<()> {
        let cloud = loaded_cloud().scale_copy(25.4)?;

        assert_relative_eq!(
            cloud.points()[1],
            Point3::new(25.4, 0.0, 0.0),
            epsilon = 1.0e-12
        );

        for (actual, expected) in cloud
            .point_stdev()
            .unwrap()
            .iter()
            .zip([0.1, 0.2, 0.3, 0.4].iter())
        {
            assert_relative_eq!(*actual, expected * 25.4, epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn a_negative_scale_mirrors_and_flips_the_normals() -> Result<()> {
        let mirrored = loaded_cloud().scale_copy(-1.0)?;

        assert_relative_eq!(
            mirrored.points()[1],
            Point3::new(-1.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            mirrored.point_normals().unwrap()[0].into_inner(),
            -Vector3::z(),
            epsilon = 1.0e-12
        );

        // Standard deviations stay positive despite the negative factor.
        for s in mirrored.point_stdev().unwrap() {
            assert!(*s > 0.0, "a standard deviation must not go negative");
        }

        Ok(())
    }

    #[test]
    fn scale_rejects_zero_and_non_finite_factors() {
        let mut cloud = loaded_cloud();

        assert!(cloud.scale_in_place(0.0).is_err());
        assert!(cloud.scale_in_place(f64::NAN).is_err());
        assert!(cloud.scale_in_place(f64::INFINITY).is_err());

        // The rejected calls must have left the cloud alone.
        assert_relative_eq!(
            cloud.points()[1],
            Point3::new(1.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );
    }

    #[test]
    fn append_unions_the_attributes() -> Result<()> {
        let mut cloud = loaded_cloud();
        cloud.append_in_place(&loaded_cloud())?;

        assert_eq!(cloud.point_count(), 8);
        assert_eq!(
            cloud.point_stdev().unwrap(),
            &[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4]
        );
        assert_eq!(cloud.point_attr("confidence").unwrap().len(), 8);
        cloud.attrs().validate(cloud.point_count())?;

        Ok(())
    }

    #[test]
    fn append_rejects_mismatched_attributes_without_modifying_the_target() -> Result<()> {
        let mut cloud = loaded_cloud();
        let mut other = loaded_cloud();
        other.set_point_stdev(None)?;

        assert!(cloud.append_in_place(&other).is_err());

        assert_eq!(cloud.point_count(), 4);
        assert_eq!(cloud.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);

        Ok(())
    }

    // ===========================================================================================
    // Subsetting
    // ===========================================================================================

    #[test]
    fn subset_by_mask_carries_every_attribute() -> Result<()> {
        let cloud = loaded_cloud();
        let mask = IndexMask::try_from_indices(&[1, 3], 4)?;

        let sub = cloud.extract_subset_points(&mask)?;

        assert_eq!(sub.point_count(), 2);
        assert_relative_eq!(sub.points()[0], Point3::new(1.0, 0.0, 0.0), epsilon = 0.0);
        assert_eq!(sub.point_stdev().unwrap(), &[0.2, 0.4]);
        assert_eq!(sub.point_colors().unwrap(), &[[1, 1, 1], [3, 3, 3]]);
        assert_eq!(
            sub.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.6, 0.8]
        );
        sub.attrs().validate(sub.point_count())?;

        Ok(())
    }

    #[test]
    fn subset_by_mask_rejects_a_length_mismatch() {
        let cloud = loaded_cloud();
        assert!(
            cloud
                .extract_subset_points(&IndexMask::new(3, true))
                .is_err()
        );
    }

    #[test]
    fn subset_by_indices_preserves_order_and_allows_repeats() -> Result<()> {
        let cloud = loaded_cloud();
        let sub = cloud.extract_subset_indices(&[3, 0, 3])?;

        assert_eq!(sub.point_count(), 3);
        assert_relative_eq!(sub.points()[0], Point3::new(0.0, 0.0, 1.0), epsilon = 0.0);
        assert_relative_eq!(sub.points()[1], Point3::origin(), epsilon = 0.0);
        assert_eq!(sub.point_stdev().unwrap(), &[0.4, 0.1, 0.4]);
        sub.attrs().validate(sub.point_count())?;

        Ok(())
    }

    #[test]
    fn subset_by_indices_rejects_out_of_bounds() {
        assert!(loaded_cloud().extract_subset_indices(&[0, 4]).is_err());
    }

    // ===========================================================================================
    // Conversion
    // ===========================================================================================

    /// Indexing a cloud no longer copies it into a second type, so nothing can be lost on the way
    /// to a spatial query. This is what replaced the old lossy `to_cloud`/`from_cloud` round trip,
    /// which could not carry the open-map attributes.
    #[test]
    fn indexing_a_cloud_does_not_disturb_its_attributes() -> Result<()> {
        let cloud = loaded_cloud();
        let index = cloud.compute_index()?;

        assert_eq!(index.points(), cloud.points());
        assert_eq!(index.len(), cloud.point_count());

        // Every point is its own nearest neighbor, which is the cheapest check that the tree was
        // built over the buffer the index reports.
        for (i, p) in cloud.points().iter().enumerate() {
            let (found, d) = index.nearest_one(p);
            assert_eq!(found, i);
            assert!(d < 1e-12);
        }

        // The cloud is untouched, open-map attributes included.
        assert!(cloud.point_attr("confidence").is_some());

        Ok(())
    }

    #[test]
    fn check_attribute_loss_only_complains_when_there_is_something_to_lose() -> Result<()> {
        let bare = PointCloud3::new(vec![Point3::origin()]);
        bare.check_attribute_loss("XYZ", false)?;

        let loaded = loaded_cloud();
        assert!(loaded.check_attribute_loss("XYZ", false).is_err());
        loaded.check_attribute_loss("XYZ", true)?;

        Ok(())
    }

    // ===========================================================================================
    // PLY serialization
    // ===========================================================================================

    /// A file on disk, removed when the test finishes with it. `load_ply` and `save_ply` take paths
    /// rather than readers, so these tests cannot stay in memory.
    #[cfg(feature = "ply")]
    struct TempPly {
        path: std::path::PathBuf,
    }

    #[cfg(feature = "ply")]
    impl TempPly {
        fn new(name: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "engeom-cloud-{}-{}.ply",
                name,
                std::process::id()
            ));
            Self { path }
        }

        fn path(&self) -> &Path {
            &self.path
        }
    }

    #[cfg(feature = "ply")]
    impl Drop for TempPly {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.path);
        }
    }

    #[cfg(feature = "ply")]
    #[test]
    fn a_ply_round_trip_preserves_every_attribute() -> Result<()> {
        let before = loaded_cloud();
        let file = TempPly::new("round-trip");

        before.save_ply(file.path(), &PlyWriteOpts::default())?;
        let after = PointCloud3::load_ply(file.path())?;

        assert_eq!(after.point_count(), before.point_count());
        for (p, q) in before.points().iter().zip(after.points()) {
            assert_relative_eq!(p, q, epsilon = 0.0);
        }

        assert_eq!(after.point_colors(), before.point_colors());
        assert_eq!(after.point_stdev(), before.point_stdev());
        assert_eq!(
            after.point_attr("confidence"),
            before.point_attr("confidence")
        );
        assert_eq!(
            after.point_attr("principal_dir"),
            before.point_attr("principal_dir")
        );

        for (a, b) in after
            .point_normals()
            .unwrap()
            .iter()
            .zip(before.point_normals().unwrap())
        {
            assert_relative_eq!(a.into_inner(), b.into_inner(), epsilon = 1.0e-15);
        }

        Ok(())
    }

    /// Loading a mesh here would throw its connectivity away, so it has to be refused rather than
    /// quietly succeeding.
    #[cfg(feature = "ply")]
    #[test]
    fn load_ply_refuses_a_file_which_is_really_a_mesh() -> Result<()> {
        use crate::geom3::MeshData3;

        let mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        let file = TempPly::new("is-a-mesh");
        mesh.save_ply(file.path(), &PlyWriteOpts::default())?;

        let result = PointCloud3::load_ply(file.path());
        assert!(result.is_err());

        // The message has to say what to do instead, not just that it failed.
        let message = result.unwrap_err().to_string();
        assert!(message.contains("MeshData3::load_ply"), "{message}");

        Ok(())
    }

    #[cfg(feature = "ply")]
    #[test]
    fn a_bare_cloud_round_trips_without_inventing_attributes() -> Result<()> {
        let before = PointCloud3::new(vec![
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(-4.0, 5.0, -6.0),
        ]);

        let file = TempPly::new("bare");
        before.save_ply(file.path(), &PlyWriteOpts::default())?;
        let after = PointCloud3::load_ply(file.path())?;

        assert_eq!(after.point_count(), 2);
        assert!(after.attrs().is_empty());

        Ok(())
    }
}
