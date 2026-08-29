//! This module contains `PointCloud3`, a plain container for a buffer of points and the
//! per-point attributes attached to them.
//!
//! It is the point-cloud counterpart to `MeshData3`, and the two are deliberately shaped the same
//! way: private buffers, a validated attribute set, and the same construction, subsetting, and
//! serialization vocabulary. The attribute half is literally the same code, `PointAttrSet3`, which
//! the mesh types also hold for their point domain.
//!
//! There is no spatial acceleration here of any kind. Call `compute_index` for a `CloudIndex3`
//! when you need nearest-neighbor queries.

use crate::common::IndexMask;
use crate::common::poisson_disk::sample_poisson_disk_all;
use crate::common::{VoxelGroups, compute_voxel_groups};
use crate::geom3::Aabb3;
use crate::geom3::attributes3::{Attr3, PointAttrSet3};
use crate::geom3::point_cloud::CloudIndex3;
use crate::{Iso3, KdTree3, Point2, Point3, Result, SurfacePoint3, UnitVec3, Vector3};
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
    /// # This is unchecked!
    ///
    /// <div class="warning">
    ///
    /// **Nothing verifies that `tree` was built from this cloud.** If it was not, every query
    /// returns indices into a point set that is not this one. Those indices are in range whenever
    /// the other cloud was at least as large, so there is no panic and no error: you get plausible,
    /// confidently wrong answers, and a downstream alignment or overlap check happily uses them.
    /// A tree built over a *smaller* cloud is worse only in that it eventually panics on a
    /// subscript, which at least tells you something is broken.
    ///
    /// The obligation is that `tree` was built from `self.points()` **as they are now**. A tree is
    /// invalidated by anything which adds, removes, reorders or moves a point:
    /// `transform_in_place`, `scale_in_place`, `append_in_place`, `points_mut`, and any subset
    /// extraction, which produces a different cloud entirely. Attribute changes do not invalidate
    /// it, since the tree indexes positions only.
    ///
    /// </div>
    ///
    /// # Use [`PointCloud3::compute_index`] instead
    ///
    /// This constructor was created for the Python bindings.  You _probably_ don't mean to be using
    /// it, and if you do, you either know what you're doing and plan on doing it carefully, or
    /// you're going to find out.
    ///
    /// The checked path builds the tree itself, so it cannot be mispaired, and the borrow checker
    /// then prevents the cloud from being mutated while the index lives. It is also nearly free
    /// relative to using an index at all: a build costs about what `0.4 * N` nearest-neighbor
    /// queries do (`engeom/tests/kd_tree_backends.md`), and any real consumer does more than that.
    ///
    /// This exists for use in a long-lived wrapper that caches the tree across many calls and
    /// manages the the complexity of being sure to drop it on every mutation.  Aka, the Python
    /// bindings.  If you're not doing something similar, you almost certainly want `compute_index`
    /// instead of this.
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

    /// Get the per-point flat coordinates, if present: each point's position in a flattened 2D
    /// chart of the surface, expressed in the cloud's own length units. See
    /// [`set_point_flat`](Self::set_point_flat).
    pub fn point_flat(&self) -> Option<&[Point2]> {
        self.attrs.flat()
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

    /// Set or clear the per-point flat coordinates: the position of each point in a flattened 2D
    /// chart of the surface, such as the output of boundary first flattening. These are not
    /// texture coordinates; they express the cloud's own length units in a plane. They scale with
    /// the geometry under a uniform scale while remaining fixed under a rigid transform.
    ///
    /// A cloud projected onto a flattened reference surface can store each point's flat position
    /// here alongside its depth in an open scalar attribute.
    ///
    /// # Arguments
    ///
    /// * `values`: the flat coordinates to store, or `None` to clear them. Must match the point
    ///   count.
    ///
    /// returns: `Result<()>`
    pub fn set_point_flat(&mut self, values: Option<Vec<Point2>>) -> Result<()> {
        self.attrs.set_flat(values, self.points.len())
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
// Voxel reduction
// ===============================================================================================

/// The name under which [`PointCloud3::reduce_by_voxel`] records how well each voxel's normals
/// agreed, as a `Scalar` in `[0, 1]`.
pub const VOXEL_COHERENCE_ATTR: &str = "voxel_coherence";

/// The name under which [`PointCloud3::reduce_by_voxel`] records how many original points went into
/// each output point, as a `Label`.
pub const VOXEL_COUNT_ATTR: &str = "voxel_count";

impl PointCloud3 {
    /// Reduce the cloud onto a coarser grid, replacing the points in each voxel with a single
    /// averaged point.
    ///
    /// This is not the same operation as [`PointCloud3::sample_poisson_disk`] or
    /// `compute_first_per_voxel_mask`, both of which *select* original points. This one _creates_
    /// new ones, so the output positions are not a subset of the input and the noise on them is
    /// lower, since averaging `n` independent samples of a surface reduces the noise by `√n`.
    ///
    /// That makes it a good tool for building a smooth, coarse representation, and a pretty
    /// terrible one to take measurements off of.
    ///
    /// Note that output points are voxel centroids, not voxel centers, so two adjacent outputs can
    /// be closer together than `voxel_size`. There is no minimum spacing guarantee; use Poisson disk
    /// sampling if that is what you need.
    ///
    /// # The output sits inside the curvature
    ///
    /// A centroid of points spread over a curved patch does not lie on that patch, it lies inside
    /// it. For a patch of width `v` on a surface of local radius `R` the offset is roughly
    /// `v² / (16 * R)`, inward along the normal, and it does not shrink with more input points
    /// because it is not noise. Averaging removes noise and adds this in its place.
    ///
    /// The v² is what makes it easy to underestimate. On a 10 mm radius, a 2 mm grid gives about
    /// 25 um, which is negligible next to most scan noise; a 8 mm grid on the same surface gives
    /// about 400 um, which not negligible at all. Measured on a sphere of radius 10 the
    /// estimate is good to three figures at the coarse end and conservative by about 30% at the fine
    /// end.
    ///
    /// This matters most if the result becomes an alignment target, since the offset is
    /// systematic rather than random and therefore does not average out of the fit. Where the points
    /// must be on the surface, sample the surface directly: `Mesh3::sample_voxel_surface` does the
    /// equivalent reduction for a mesh by projecting onto the geometry instead of averaging.
    ///
    /// # How each attribute is combined
    ///
    /// - **normals**: summed and renormalized. Where the sum is degenerate, meaning the voxel's
    ///   normals cancel out, the lowest-indexed point's normal is kept so the output stays a unit
    ///   vector, and the coherence below reports zero.
    /// - **stdev**: `sqrt(sum(sigma_i^2)) / n`, the standard deviation of the mean of `n`
    ///   independent measurements. For equal inputs this is the familiar `sigma / sqrt(n)`.
    /// - **colors**: per-channel mean, rounded.
    /// - **open-map `Scalar`, `Vector`, `Color`**: mean.
    /// - **open-map `Label`**: the most common value, ties broken toward the lower value so the
    ///   result does not depend on iteration order.
    ///
    /// # Two attributes are added
    ///
    /// [`VOXEL_COHERENCE_ATTR`] is the length of the mean normal *before* renormalizing, in
    /// `[0, 1]`. It is near 1 where a voxel's normals all agree and drops toward 0 where the voxel
    /// straddles an edge or a thin wall, so a voxel whose averaged position is meaningless says so
    /// about itself. It is only written when the cloud has normals.
    ///
    /// [`VOXEL_COUNT_ATTR`] is how many input points each output point came from, which is what the
    /// `sqrt(n)` noise argument above depends on, and is useful as a weight.
    ///
    /// Both names would collide if the input already carries them, so an input which does is
    /// rejected rather than silently overwritten.
    ///
    /// # Arguments
    ///
    /// * `voxel_size`: the edge length of the grid cells, which must be finite and positive
    ///
    /// returns: `Result<PointCloud3>`
    pub fn reduce_by_voxel(&self, voxel_size: f64) -> Result<Self> {
        for name in [VOXEL_COHERENCE_ATTR, VOXEL_COUNT_ATTR] {
            if self.attrs.attr(name).is_some() {
                return Err(format!(
                    "The cloud already has an attribute named '{name}', which the reduction would \
                     overwrite. Remove it first if it is no longer wanted."
                )
                .into());
            }
        }

        let groups = compute_voxel_groups(&self.points, voxel_size)?;
        let n_out = groups.len();

        let mut points = Vec::with_capacity(n_out);
        let mut counts = Vec::with_capacity(n_out);

        for g in groups.iter() {
            let mut sum = Vector3::zeros();
            for &i in g {
                sum += self.points[i as usize].coords;
            }
            points.push(Point3::from(sum / g.len() as f64));
            counts.push(g.len() as u32);
        }

        let mut attrs = PointAttrSet3::empty();

        if let Some(normals) = self.attrs.normals() {
            let mut reduced = Vec::with_capacity(n_out);
            let mut coherence = Vec::with_capacity(n_out);

            for g in groups.iter() {
                let mut sum = Vector3::zeros();
                for &i in g {
                    sum += normals[i as usize].into_inner();
                }

                let mean = sum / g.len() as f64;
                let length = mean.norm();

                // The length of the mean of unit vectors is 1 when they all agree and 0 when they
                // cancel, which is the coherence signal. Renormalizing throws it away, so it is
                // captured first.
                coherence.push(length.min(1.0));
                reduced.push(if length > 1e-12 {
                    UnitVec3::new_unchecked(mean / length)
                } else {
                    normals[g[0] as usize]
                });
            }

            attrs.set_normals(Some(reduced), n_out)?;
            attrs.insert_attr(VOXEL_COHERENCE_ATTR, Attr3::Scalar(coherence), n_out)?;
        }

        if let Some(stdev) = self.attrs.stdev() {
            let reduced = groups
                .iter()
                .map(|g| {
                    let sum_sq: f64 = g
                        .iter()
                        .map(|&i| stdev[i as usize] * stdev[i as usize])
                        .sum();
                    sum_sq.sqrt() / g.len() as f64
                })
                .collect();
            attrs.set_stdev(Some(reduced), n_out)?;
        }

        if let Some(colors) = self.attrs.colors() {
            let reduced = groups
                .iter()
                .map(|g| mean_color(g.iter().map(|&i| colors[i as usize])))
                .collect();
            attrs.set_colors(Some(reduced), n_out)?;
        }

        for name in self
            .attrs
            .attr_names()
            .map(str::to_string)
            .collect::<Vec<_>>()
        {
            let attr = self
                .attrs
                .attr(&name)
                .expect("the name came from this attribute set");
            attrs.insert_attr(&name, reduce_attr(attr, &groups), n_out)?;
        }

        attrs.insert_attr(VOXEL_COUNT_ATTR, Attr3::Label(counts), n_out)?;

        Ok(Self { points, attrs })
    }
}

/// The per-channel mean of a set of colors, rounded rather than truncated.
fn mean_color(colors: impl Iterator<Item = [u8; 3]>) -> [u8; 3] {
    let mut sums = [0u32; 3];
    let mut n = 0u32;
    for c in colors {
        for (s, v) in sums.iter_mut().zip(c.iter()) {
            *s += u32::from(*v);
        }
        n += 1;
    }

    let mut out = [0u8; 3];
    for (o, s) in out.iter_mut().zip(sums.iter()) {
        *o = ((*s + n / 2) / n) as u8;
    }
    out
}

/// Combine an open-map attribute across each voxel, by the rules in `reduce_by_voxel`.
fn reduce_attr(attr: &Attr3, groups: &VoxelGroups) -> Attr3 {
    match attr {
        Attr3::Scalar(values) => Attr3::Scalar(
            groups
                .iter()
                .map(|g| g.iter().map(|&i| values[i as usize]).sum::<f64>() / g.len() as f64)
                .collect(),
        ),
        Attr3::Vector(values) => Attr3::Vector(
            groups
                .iter()
                .map(|g| {
                    let sum: Vector3 = g.iter().map(|&i| values[i as usize]).sum();
                    sum / g.len() as f64
                })
                .collect(),
        ),
        Attr3::Color(values) => Attr3::Color(
            groups
                .iter()
                .map(|g| mean_color(g.iter().map(|&i| values[i as usize])))
                .collect(),
        ),
        Attr3::Label(values) => Attr3::Label(
            groups
                .iter()
                .map(|g| modal_label(g.iter().map(|&i| values[i as usize])))
                .collect(),
        ),
    }
}

/// The most common label in a group, ties broken toward the lower value.
///
/// A label names a category, so averaging it would invent a value that means nothing. The mode is
/// the only summary that keeps the result inside the original set of categories.
fn modal_label(labels: impl Iterator<Item = u32>) -> u32 {
    let mut sorted: Vec<u32> = labels.collect();
    sorted.sort_unstable();

    let mut best = sorted[0];
    let mut best_count = 0usize;
    let mut current = sorted[0];
    let mut count = 0usize;

    for v in sorted {
        if v == current {
            count += 1;
        } else {
            current = v;
            count = 1;
        }

        // Strictly greater keeps the first, and the scan is ascending, so ties go to the lower.
        if count > best_count {
            best = current;
            best_count = count;
        }
    }

    best
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
    // Voxel reduction
    // ===========================================================================================

    /// Four points in one voxel collapse to their centroid, and the count records how many.
    #[test]
    fn reduce_by_voxel_averages_positions() -> Result<()> {
        let cloud = PointCloud3::new(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            Point3::new(0.0, 0.4, 0.0),
            Point3::new(0.2, 0.4, 0.0),
            // A second voxel, well clear of the first.
            Point3::new(5.5, 0.0, 0.0),
        ]);

        let out = cloud.reduce_by_voxel(1.0)?;

        assert_eq!(out.point_count(), 2);
        assert_relative_eq!(out.points()[0], Point3::new(0.1, 0.2, 0.0), epsilon = 1e-12);
        assert_relative_eq!(out.points()[1], Point3::new(5.5, 0.0, 0.0), epsilon = 1e-12);

        let counts = out
            .point_attr(VOXEL_COUNT_ATTR)
            .and_then(|a| a.as_label())
            .expect("counts were written");
        assert_eq!(counts, &[4, 1]);

        Ok(())
    }

    /// Coherence is the whole point of the reduction reporting on itself: near 1 where a voxel's
    /// normals agree, near 0 where they cancel.
    #[test]
    fn reduce_by_voxel_reports_normal_coherence() -> Result<()> {
        let mut cloud = PointCloud3::new(vec![
            // A flat voxel: both normals agree.
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            // An edge voxel: normals at right angles.
            Point3::new(5.1, 0.0, 0.0),
            Point3::new(5.2, 0.0, 0.0),
            // A fold voxel: normals directly opposed.
            Point3::new(9.1, 0.0, 0.0),
            Point3::new(9.2, 0.0, 0.0),
        ]);
        cloud.set_point_normals(Some(vec![
            UnitVec3::new_normalize(Vector3::z()),
            UnitVec3::new_normalize(Vector3::z()),
            UnitVec3::new_normalize(Vector3::z()),
            UnitVec3::new_normalize(Vector3::x()),
            UnitVec3::new_normalize(Vector3::z()),
            UnitVec3::new_normalize(-Vector3::z()),
        ]))?;

        let out = cloud.reduce_by_voxel(1.0)?;
        assert_eq!(out.point_count(), 3);

        let coherence = out
            .point_attr(VOXEL_COHERENCE_ATTR)
            .and_then(|a| a.as_scalar())
            .expect("coherence was written");

        assert_relative_eq!(coherence[0], 1.0, epsilon = 1e-12);
        assert_relative_eq!(coherence[1], 0.5_f64.sqrt(), epsilon = 1e-12);
        assert_relative_eq!(coherence[2], 0.0, epsilon = 1e-12);

        // Even where the normals cancel the output must still be a unit vector.
        let normals = out.point_normals().expect("normals were written");
        for n in normals {
            assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-12);
        }

        Ok(())
    }

    /// Averaging `n` independent measurements divides the standard deviation by `sqrt(n)`.
    #[test]
    fn reduce_by_voxel_propagates_stdev() -> Result<()> {
        let mut cloud = PointCloud3::new(vec![
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            Point3::new(0.3, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.0),
        ]);
        cloud.set_point_stdev(Some(vec![0.01; 4]))?;

        let out = cloud.reduce_by_voxel(1.0)?;

        let stdev = out.point_stdev().expect("stdev was written");
        assert_eq!(stdev.len(), 1);
        assert_relative_eq!(stdev[0], 0.01 / 2.0, epsilon = 1e-12);

        Ok(())
    }

    /// A label names a category, so it takes the mode rather than a mean which would invent a
    /// value outside the original set.
    #[test]
    fn reduce_by_voxel_takes_the_modal_label() -> Result<()> {
        let mut cloud = PointCloud3::new(vec![
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            Point3::new(0.3, 0.0, 0.0),
            // A voxel where two labels tie, which must resolve to the lower.
            Point3::new(5.1, 0.0, 0.0),
            Point3::new(5.2, 0.0, 0.0),
        ]);
        cloud.insert_point_attr("pass", Attr3::Label(vec![7, 7, 3, 9, 4]))?;

        let out = cloud.reduce_by_voxel(1.0)?;

        let labels = out
            .point_attr("pass")
            .and_then(|a| a.as_label())
            .expect("the label attribute came through");
        assert_eq!(labels, &[7, 4]);

        Ok(())
    }

    #[test]
    fn reduce_by_voxel_averages_open_map_scalars_and_colors() -> Result<()> {
        let mut cloud =
            PointCloud3::new(vec![Point3::new(0.1, 0.0, 0.0), Point3::new(0.2, 0.0, 0.0)]);
        cloud.insert_point_attr("deviation", Attr3::Scalar(vec![1.0, 3.0]))?;
        cloud.set_point_colors(Some(vec![[0, 10, 255], [2, 11, 255]]))?;

        let out = cloud.reduce_by_voxel(1.0)?;

        let dev = out
            .point_attr("deviation")
            .and_then(|a| a.as_scalar())
            .expect("scalar came through");
        assert_relative_eq!(dev[0], 2.0, epsilon = 1e-12);

        // 10 and 11 average to 10.5, which rounds up rather than truncating down.
        assert_eq!(
            out.point_colors().expect("colors came through"),
            &[[1, 11, 255]]
        );

        Ok(())
    }

    #[test]
    fn reduce_by_voxel_refuses_to_overwrite_its_own_attributes() -> Result<()> {
        let mut cloud = PointCloud3::new(vec![Point3::origin()]);
        cloud.insert_point_attr(VOXEL_COUNT_ATTR, Attr3::Label(vec![1]))?;

        assert!(cloud.reduce_by_voxel(1.0).is_err());

        Ok(())
    }

    #[test]
    fn reduce_by_voxel_rejects_a_nonsense_size() {
        let cloud = PointCloud3::new(vec![Point3::origin()]);
        assert!(cloud.reduce_by_voxel(0.0).is_err());
        assert!(cloud.reduce_by_voxel(f64::NAN).is_err());
    }

    /// A reduction of a sphere sample must still lie on the sphere, within the error a centroid of
    /// points inside one voxel can introduce.
    /// The mesh route wanted no new mesh code: a dense surface sample already carries normals, so
    /// mesh -> sample -> cloud -> reduce works with one reducer rather than a second one written
    /// for meshes. This is the check that the route actually composes.
    #[test]
    fn a_mesh_reduces_through_a_sample_without_mesh_specific_code() -> Result<()> {
        let mesh = crate::Mesh3::create_sphere(50.0, 0.077)?;
        let cloud = mesh.sample_dense(1.0, None);
        assert!(
            cloud.point_normals().is_some(),
            "the sample carried normals"
        );

        let out = cloud.reduce_by_voxel(5.0)?;

        assert!(out.point_count() < cloud.point_count() / 4);
        assert!(out.point_normals().is_some());

        // A sphere is smooth, so every voxel's normals should agree well.
        let coherence = out
            .point_attr(VOXEL_COHERENCE_ATTR)
            .and_then(|a| a.as_scalar())
            .expect("coherence was written");
        let worst = coherence.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            worst > 0.9,
            "worst coherence on a smooth sphere was {worst}"
        );

        Ok(())
    }

    /// The same reduction over a box, where the coherence signal has something to find: voxels
    /// straddling an edge must score lower than voxels in the middle of a face.
    #[test]
    fn reduce_by_voxel_coherence_finds_the_edges_of_a_box() -> Result<()> {
        let mesh = crate::Mesh3::create_box(40.0, 40.0, 40.0, false);
        let cloud = mesh.sample_dense(0.5, None);

        let out = cloud.reduce_by_voxel(4.0)?;
        let coherence = out
            .point_attr(VOXEL_COHERENCE_ATTR)
            .and_then(|a| a.as_scalar())
            .expect("coherence was written");

        let worst = coherence.iter().copied().fold(f64::INFINITY, f64::min);
        let best = coherence.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        assert!(
            best > 0.99,
            "a flat face should be almost perfectly coherent"
        );
        assert!(
            worst < 0.9,
            "an edge voxel should be measurably less coherent, worst was {worst}"
        );

        Ok(())
    }

    #[test]
    fn reduce_by_voxel_keeps_a_sphere_on_its_surface() -> Result<()> {
        let radius = 50.0;
        let n = 20_000;
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        let points: Vec<Point3> = (0..n)
            .map(|i| {
                let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0;
                let r = (1.0 - y * y).max(0.0).sqrt();
                let theta = golden * i as f64;
                Point3::new(
                    radius * theta.cos() * r,
                    radius * y,
                    radius * theta.sin() * r,
                )
            })
            .collect();

        let voxel = 5.0;
        let out = PointCloud3::new(points).reduce_by_voxel(voxel)?;

        assert!(out.point_count() < n / 4, "the reduction did not reduce");

        // A centroid of points inside one voxel can fall at most half a voxel diagonal off the
        // surface, and in practice far less.
        let limit = voxel * 3.0_f64.sqrt() / 2.0;
        for p in out.points() {
            assert!(
                (p.coords.norm() - radius).abs() < limit,
                "reduced point {p:?} left the sphere by more than {limit}"
            );
        }

        Ok(())
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
