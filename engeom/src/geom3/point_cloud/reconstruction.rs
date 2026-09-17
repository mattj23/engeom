//! Building a triangle mesh from a single oriented point cloud.
//!
//! This module provides the main entry point to the reconstruction machinery. The individual steps
//! are available elsewhere in the crate, but assembling them requires several constraints that are
//! easy to overlook. The field's support must be wide enough to retain all eight corners of a cell,
//! a triangle soup from marching cubes requires repair before it becomes a mesh, and an incorrectly
//! oriented normal can create a hole.
//!
//! # What it does, in order
//!
//! 1. Obtain normals from the cloud or estimate them, then orient them. See [`NormalSource3`].
//! 2. If requested, remove points whose normals were fitted from an unreliable neighborhood.
//! 3. Build an implicit field over the remaining points and sample it onto a narrow band of voxels.
//! 4. Extract the zero level set as a triangle soup.
//! 5. Repair the soup into a manifold mesh, then optionally remove small disconnected patches and
//!    smooth the remaining geometry. These operations work directly on point and face buffers.
//!    Repair omits the two edge-topology passes because extraction in step 4 already prevents the
//!    defects that they detect. See
//!    [`RepairOpts::assuming_oriented_edges`].
//!
//! # Why it hands back buffers rather than a `Mesh3`
//!
//! The result is [`MeshData3`]. This pipeline does not query a bounding volume hierarchy: repair
//! works on raw buffers, patch filtering walks the face list, and smoothing uses the half-edge
//! structure. Each operation uses topology and arithmetic over points and faces.
//!
//! `Mesh3` builds the hierarchy eagerly during construction. For more than half a million
//! triangles, this takes longer than all other reconstruction stages combined. Returning `Mesh3`
//! would build an unused hierarchy even when the caller only writes the result to a file.
//!
//! Call `Mesh3::from_data(result, false)` when the caller needs closest-point queries, ray casts, or
//! deviation measurement. Otherwise, retain the unaccelerated result.
//!
//! # What it does not do
//!
//! The result is not closed, and the API makes no claim that it is. A scan sees one side of an
//! object, so reconstruction stops where the data stops. Producing a closed result from partial
//! data requires generating the missing surface, which this module deliberately does not do. The
//! caller must also select `is_solid` during conversion because this process cannot determine
//! whether the result encloses a volume.
//!
//! This module also reconstructs one cloud at a time. If several scans of the same part are merged
//! into one cloud first, disagreements between scans produce a double surface because the field
//! averages the conflicting data. Resolving such conflicts is a separate operation. This separation
//! is why [`super::implicit_field::compute_sdf_grid`] accepts a field rather than a cloud.
//!
//! # What to expect from it
//!
//! Measured against the primitives it was built from, sampled at a quarter unit and reconstructed
//! at a voxel of the same size:
//!
//! | Shape | Reference to result | Result to reference | Pieces |
//! |---|---|---|---|
//! | Sphere | 0.04 voxels | 0.04 voxels | 1 |
//! | Cylinder | 0.32 voxels | 0.45 voxels | 28 |
//! | Box | 0.43 voxels | 0.93 voxels | 9 |
//!
//! This method performs well on the sphere. The cylinder and box demonstrate its behavior at sharp
//! edges: it rounds a crease across the field width, moving the surface and potentially separating
//! it where the rounded surface pulls away from both faces. This behavior results from averaging
//! tangent planes. A field that discounts disagreeing neighbors can preserve these edges. With the
//! current field, parts dominated by creases reconstruct less accurately than mostly smooth parts;
//! use the piece count to detect this problem.

use super::implicit_field::{ImlsField3, compute_sdf_grid};
use super::normal_estimation::NormalEstimates;
use super::orientation::{NormalOrientation3, OrientationReport3};
use super::{CloudIndex3, PointCloud3};
use crate::common::IndexMask;
use crate::common::kd_tree::{KdTree, KdTreeSearch};
use crate::geom3::half_edge3::HalfEdgeMesh3;
use crate::geom3::half_edge3::{RepairOpts, RepairReport, repair_buffers};
use crate::geom3::mesh::PatchFilter;
use crate::raster3::extract_isosurface;
use crate::{MeshData3, Point3, Result, UnitVec3};

/// Where the normals of a reconstruction come from.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum NormalSource3 {
    /// Use the cloud's existing normals. The normals must be present and oriented.
    Existing,

    /// Estimate normals by fitting a plane to the neighbors within `radius`, then orient them with
    /// `orientation`.
    Estimate {
        radius: f64,
        orientation: NormalOrientation3,
    },
}

/// Which implicit field a reconstruction builds.
///
/// The enum currently has one variant. Its structure permits adding a field that preserves sharp
/// edges without changing the options API.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FieldKind3 {
    /// The weighted average of nearby tangent planes. See [`super::implicit_field`] for its
    /// behavior, including the bias it has on curved surfaces.
    #[default]
    Imls,
}

/// How a reconstruction should be carried out.
///
/// Voxel size is the only value without a suitable default, so [`ReconstructOpts3::new`] requires
/// it. Use the `with_` methods to adjust the remaining options.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub struct ReconstructOpts3 {
    /// The spacing of the voxel grid, in the units of the cloud. This sets the resolution of the
    /// result. It should be at least the spacing of the points. See
    /// [`CloudIndex3::estimate_point_spacing`].
    pub voxel_size: f64,

    /// How far the field reaches, in voxels. Below about `1.75` a cell holding a piece of surface
    /// can lose a corner and be dropped, creating holes in the result. Values below `1.75` are
    /// rejected. Larger values increase smoothing and computation cost.
    pub band_width: f64,

    /// Where the normals come from.
    pub normals: NormalSource3,

    /// Which field to build.
    pub field: FieldKind3,

    /// Points whose normal confidence falls below this are left out of the field entirely.
    ///
    /// A neighborhood with fewer than three points receives an arbitrary normal with zero
    /// confidence. A neighborhood that straddles an edge receives a normal between the adjacent
    /// surface normals. Either case can add an incorrectly signed contribution to the field. A
    /// threshold of zero keeps all points.
    pub min_normal_confidence: f64,

    /// Which repair passes to run over the raw triangle soup. `None` skips repair, which leaves a
    /// mesh that may not be manifold.
    ///
    /// The default is [`RepairOpts::assuming_oriented_edges`] rather than
    /// [`RepairOpts::default`] because extraction already prevents the defects that the two omitted
    /// passes detect. Across spheres, a cylinder, a box, and the Stanford bunny, neither pass found
    /// a defect. Together, they account for 83 percent of repair time, while repair accounts for 80
    /// percent of the total reconstruction time. Pass `Some(RepairOpts::default())` to run them.
    pub repair: Option<RepairOpts>,

    /// Which connected patches to keep. `None`, the default, keeps all of them, including the
    /// slivers that collect along the edge of the band.
    ///
    /// # Keeping only the largest patch is not a cleanup
    ///
    /// `PatchFilter::keep_largest` assumes that the surface is one piece and that every other piece
    /// is debris. This assumption does not hold for all shapes with sharp edges. The field rounds a
    /// crease and can pull far enough from both faces that the surfaces do not join. A single part
    /// can therefore produce several valid pieces without any debris.
    ///
    /// A cylinder sampled at a quarter unit reconstructs into 28 pieces this way. Measured against
    /// the original, the whole result sits within 0.32 voxels; keeping only the largest piece
    /// removes a cap and increases that distance to 8.86 voxels, a factor of twenty-seven. The
    /// filter behaves as requested, but it is unsuitable for this shape.
    ///
    /// Prefer a filter that defines debris explicitly, such as a minimum face count or minimum area
    /// fraction, and compare the patch count with the expected part geometry. Use `keep_largest`
    /// only when the shape is one smooth piece.
    pub patch_filter: Option<PatchFilter>,

    /// How many smoothing passes to run at the end. Zero leaves the surface as extracted.
    pub smooth_iterations: usize,
}

impl ReconstructOpts3 {
    /// Default options at a given voxel size.
    ///
    /// By default, the cloud must already contain oriented normals. Reconstruction repairs the
    /// extracted mesh, keeps all patches, and does not smooth the result. The repair omits the two
    /// passes the extraction makes unnecessary; see [`ReconstructOpts3::repair`].
    pub fn new(voxel_size: f64) -> Self {
        Self {
            voxel_size,
            band_width: 2.0,
            normals: NormalSource3::Existing,
            field: FieldKind3::Imls,
            min_normal_confidence: 0.0,
            repair: Some(RepairOpts::assuming_oriented_edges()),
            patch_filter: None,
            smooth_iterations: 0,
        }
    }

    pub fn with_band_width(mut self, band_width: f64) -> Self {
        self.band_width = band_width;
        self
    }

    pub fn with_normals(mut self, normals: NormalSource3) -> Self {
        self.normals = normals;
        self
    }

    pub fn with_field(mut self, field: FieldKind3) -> Self {
        self.field = field;
        self
    }

    pub fn with_min_normal_confidence(mut self, min_normal_confidence: f64) -> Self {
        self.min_normal_confidence = min_normal_confidence;
        self
    }

    pub fn with_repair(mut self, repair: Option<RepairOpts>) -> Self {
        self.repair = repair;
        self
    }

    pub fn with_patch_filter(mut self, patch_filter: Option<PatchFilter>) -> Self {
        self.patch_filter = patch_filter;
        self
    }

    pub fn with_smooth_iterations(mut self, smooth_iterations: usize) -> Self {
        self.smooth_iterations = smooth_iterations;
        self
    }

    /// The distance the field reaches, in the units of the cloud.
    pub fn band(&self) -> f64 {
        self.band_width * self.voxel_size
    }

    fn validate(&self) -> Result<()> {
        if !self.voxel_size.is_finite() || self.voxel_size <= 0.0 {
            return Err(format!(
                "Voxel size must be finite and positive, got {}",
                self.voxel_size
            )
            .into());
        }

        if !self.band_width.is_finite() || self.band_width < 1.75 {
            return Err(format!(
                "Band width must be at least 1.75 voxels, got {}. Below that a cell holding a \
                 piece of surface can lose a corner and be left out, which puts holes in the \
                 result.",
                self.band_width
            )
            .into());
        }

        if !self.min_normal_confidence.is_finite() || self.min_normal_confidence < 0.0 {
            return Err(format!(
                "Minimum normal confidence must be finite and non-negative, got {}",
                self.min_normal_confidence
            )
            .into());
        }

        if let NormalSource3::Estimate { radius, .. } = &self.normals
            && (!radius.is_finite() || *radius <= 0.0)
        {
            return Err(format!(
                "Normal estimation radius must be finite and positive, got {radius}"
            )
            .into());
        }

        Ok(())
    }
}

/// What a reconstruction did at each stage.
///
/// The report helps diagnose reconstruction problems. A large `cells_skipped_unknown` value next
/// to a small `cells_visited` value indicates that the band was too thin. An orientation
/// `components` value greater than one indicates that normal orientation required separate choices
/// in multiple regions. A large `points_dropped` value indicates that the confidence threshold
/// removed more of the cloud than intended.
#[derive(Clone, Debug, Default)]
pub struct ReconstructReport3 {
    /// Points that contributed to the field after confidence filtering.
    pub points_used: usize,

    /// Points left out for falling below the confidence threshold.
    pub points_dropped: usize,

    /// Blocks of voxels the band covered.
    pub active_blocks: usize,

    /// Voxels the field wrote a value into.
    pub known_voxels: usize,

    /// Cells triangulated after all eight corners were evaluated.
    pub cells_visited: usize,

    /// Cells omitted because they touch an unsupported voxel at the edge of the band.
    pub cells_skipped_unknown: usize,

    /// What the normal orientation did, or `None` when the cloud arrived with its own normals.
    pub orientation: Option<OrientationReport3>,

    /// Vertices and faces produced by extraction before repair.
    pub raw_vertices: usize,
    pub raw_faces: usize,

    /// What the repair changed, or `None` when repair was skipped.
    pub repair: Option<RepairReport>,

    /// Faces in the returned mesh data.
    pub faces: usize,
}

impl CloudIndex3<'_> {
    /// The median distance from a point to its nearest neighbor.
    ///
    /// Use this value to select the voxel size. A grid much finer than the point spacing provides no
    /// additional resolution and consumes memory in proportion to the cube of the size ratio. A
    /// much coarser grid discards available scan detail. The median prevents a few distant outliers
    /// from shifting the estimate, as they would shift the mean.
    ///
    /// Returns zero for a cloud of fewer than two points.
    pub fn estimate_point_spacing(&self) -> f64 {
        if self.points().len() < 2 {
            return 0.0;
        }

        let mut distances: Vec<f64> = self
            .points()
            .iter()
            .map(|p| {
                // The nearest point is the query point itself, so request two points.
                let found = self.tree().nearest(p, 2);
                found.last().map(|(_, d)| *d).unwrap_or(0.0)
            })
            .collect();

        distances.sort_unstable_by(f64::total_cmp);
        distances[distances.len() / 2]
    }

    /// Reconstruct mesh data from this cloud.
    ///
    /// See the module documentation for the stages and for what the result does and does not
    /// promise.
    ///
    /// # Arguments
    ///
    /// * `opts`: how to carry out the reconstruction
    ///
    /// returns: `Result<(MeshData3, ReconstructReport3)>`, failing if the options are invalid, if
    /// the cloud has no normals when [`NormalSource3::Existing`] is selected, or if the field and
    /// grid produce no surface
    pub fn reconstruct_surface(
        &self,
        opts: &ReconstructOpts3,
    ) -> Result<(MeshData3, ReconstructReport3)> {
        opts.validate()?;

        if self.points().is_empty() {
            return Err("Cannot reconstruct a surface from an empty point cloud".into());
        }

        let mut report = ReconstructReport3::default();

        // ---- normals -------------------------------------------------------------------------
        let (normals, confidence) = match &opts.normals {
            NormalSource3::Existing => {
                let existing = self.cloud().point_normals().ok_or(
                    "The cloud has no normals. Either set them, or ask for them to be estimated \
                     with NormalSource3::Estimate.",
                )?;
                (existing.to_vec(), None)
            }

            NormalSource3::Estimate {
                radius,
                orientation,
            } => {
                let (estimates, orientation_report) =
                    self.estimate_normals_oriented(*radius, orientation)?;
                report.orientation = Some(orientation_report);
                let NormalEstimates {
                    normals,
                    confidence,
                } = estimates;
                (normals, Some(confidence))
            }
        };

        // ---- confidence filtering ------------------------------------------------------------
        let selected = select_confident(confidence.as_deref(), opts.min_normal_confidence);
        let (points, normals, tree) = match &selected {
            None => (self.points().to_vec(), normals, None),
            Some(mask) => {
                let kept: Vec<Point3> = mask.iter_true().map(|i| self.points()[i]).collect();
                let kept_normals: Vec<UnitVec3> = mask.iter_true().map(|i| normals[i]).collect();

                if kept.is_empty() {
                    return Err(format!(
                        "A minimum normal confidence of {} left no points at all",
                        opts.min_normal_confidence
                    )
                    .into());
                }

                let tree = KdTree::try_new(&kept)?;
                (kept, kept_normals, Some(tree))
            }
        };

        report.points_used = points.len();
        report.points_dropped = self.points().len() - points.len();

        // ---- field and band ------------------------------------------------------------------
        let band = opts.band();
        let sigma = opts.voxel_size;

        let grid = match &tree {
            Some(tree) => {
                let field = build_field(opts.field, &points, &normals, tree, band, sigma)?;
                compute_sdf_grid(&field, &points, opts.voxel_size, Point3::origin())?
            }
            None => {
                let field = build_field(opts.field, &points, &normals, self.tree(), band, sigma)?;
                compute_sdf_grid(&field, &points, opts.voxel_size, Point3::origin())?
            }
        };

        report.active_blocks = grid.block_count();
        report.known_voxels = grid.known_count();

        // ---- extraction ----------------------------------------------------------------------
        let (raw, stats) = extract_isosurface(&grid)?;
        report.cells_visited = stats.cells_visited;
        report.cells_skipped_unknown = stats.cells_skipped_unknown;
        report.raw_vertices = stats.vertices;
        report.raw_faces = stats.faces;

        if raw.faces().is_empty() {
            return Err(format!(
                "The reconstruction produced no surface. The field wrote {} voxels across {} \
                 blocks and {} cells were skipped for having an unevaluated corner. A voxel size \
                 well below the spacing of the points is the usual cause.",
                report.known_voxels, report.active_blocks, report.cells_skipped_unknown
            )
            .into());
        }

        // ---- repair --------------------------------------------------------------------------
        let mesh = match &opts.repair {
            None => MeshData3::new(raw.points().to_vec(), raw.faces().to_vec())?,
            Some(repair_opts) => {
                let repaired = repair_buffers(raw.points(), raw.faces(), repair_opts)?;
                report.repair = Some(repaired.report);
                MeshData3::new(repaired.points, repaired.faces)?
            }
        };

        // ---- patches and smoothing -----------------------------------------------------------
        let mesh = match &opts.patch_filter {
            None => mesh,
            Some(filter) => mesh.remove_small_patches(filter)?,
        };

        let mesh = if opts.smooth_iterations == 0 {
            mesh
        } else {
            let mut half_edge = HalfEdgeMesh3::try_from(mesh.view())?;
            for _ in 0..opts.smooth_iterations {
                half_edge.neighborhood_smooth()?;
            }
            half_edge.to_mesh_data()?
        };

        report.faces = mesh.faces().len();
        Ok((mesh, report))
    }
}

impl PointCloud3 {
    /// Reconstruct mesh data from this cloud, indexing it first.
    ///
    /// This method is equivalent to [`CloudIndex3::reconstruct_surface`] when the caller does not
    /// need the index afterward. For several operations on the same cloud, build the index once
    /// with [`PointCloud3::compute_index`] and use the indexed method.
    pub fn reconstruct_surface(
        &self,
        opts: &ReconstructOpts3,
    ) -> Result<(MeshData3, ReconstructReport3)> {
        self.compute_index()?.reconstruct_surface(opts)
    }
}

/// Select the points that meet the confidence threshold, or return `None` when all points meet it.
///
/// Returning `None` instead of a full mask lets the caller reuse its existing tree.
fn select_confident(confidence: Option<&[f64]>, minimum: f64) -> Option<IndexMask> {
    if minimum <= 0.0 {
        return None;
    }

    let confidence = confidence?;
    if confidence.iter().all(|c| *c >= minimum) {
        return None;
    }

    let mut mask = IndexMask::new(confidence.len(), false);
    for (i, c) in confidence.iter().enumerate() {
        if *c >= minimum {
            mask.set(i, true);
        }
    }
    Some(mask)
}

fn build_field<'a, T: KdTreeSearch<3> + Sync>(
    kind: FieldKind3,
    points: &'a [Point3],
    normals: &'a [UnitVec3],
    tree: &'a T,
    band: f64,
    sigma: f64,
) -> Result<ImlsField3<'a, T>> {
    match kind {
        FieldKind3::Imls => ImlsField3::try_new(points, normals, tree, band, sigma),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Mesh3;

    /// Reconstruct a primitive from a Poisson sample of its surface and return the deviation in
    /// voxels, the mesh, and the reconstruction report.
    fn round_trip(
        source: &Mesh3,
        spacing: f64,
        opts: &ReconstructOpts3,
    ) -> (f64, f64, MeshData3, ReconstructReport3) {
        let cloud = source
            .sample_poisson(spacing, None)
            .expect("sampling failed");
        let (mesh, report) = cloud
            .reconstruct_surface(opts)
            .expect("reconstruction failed");

        // Deviation measurement requires a bounding volume hierarchy, so convert at the point where
        // a caller would need it. Reconstruction returns buffers and leaves this choice to the
        // caller.
        let accelerated = Mesh3::from_data(mesh.clone(), false).expect("mesh build failed");
        let deviation = source
            .measure_surface_deviation(&accelerated, Some(spacing))
            .expect("deviation measurement failed");

        (
            deviation.reference_to_test / opts.voxel_size,
            deviation.test_to_reference / opts.voxel_size,
            mesh,
            report,
        )
    }

    fn patch_count(mesh: &MeshData3) -> usize {
        mesh.compute_patch_labels(None)
            .expect("labels failed")
            .count()
    }

    /// Assert a mesh is closed and consistently wound.
    fn assert_closed(mesh: &MeshData3, context: &str) {
        let mut directed = std::collections::HashMap::new();
        for face in mesh.faces() {
            for k in 0..3 {
                *directed
                    .entry((face[k], face[(k + 1) % 3]))
                    .or_insert(0usize) += 1;
            }
        }
        for (&(a, b), &n) in directed.iter() {
            assert_eq!(n, 1, "{context}: directed edge {a}->{b} used {n} times");
            assert_eq!(
                directed.get(&(b, a)).copied().unwrap_or(0),
                1,
                "{context}: edge {a}-{b} has no opposing face"
            );
        }
    }

    // ===========================================================================================
    // Options
    // ===========================================================================================

    #[test]
    fn options_reject_a_voxel_size_that_cannot_work() {
        let cloud = Mesh3::create_sphere(2.0, 0.01)
            .expect("sphere")
            .sample_poisson(0.3, None)
            .expect("sampling");

        for size in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                cloud
                    .reconstruct_surface(&ReconstructOpts3::new(size))
                    .is_err(),
                "a voxel size of {size} should have been refused"
            );
        }
    }

    /// A band thinner than the diagonal of a cell cannot contain a crossing with all eight corners
    /// evaluated. Rejecting such a band prevents a mesh with holes.
    #[test]
    fn options_reject_a_band_too_thin_to_hold_a_crossing() {
        let cloud = Mesh3::create_sphere(2.0, 0.01)
            .expect("sphere")
            .sample_poisson(0.3, None)
            .expect("sampling");

        for band_width in [0.0, 1.0, 1.7] {
            let opts = ReconstructOpts3::new(0.3).with_band_width(band_width);
            assert!(
                cloud.reconstruct_surface(&opts).is_err(),
                "a band width of {band_width} should have been refused"
            );
        }

        let opts = ReconstructOpts3::new(0.3).with_band_width(1.75);
        assert!(cloud.reconstruct_surface(&opts).is_ok());
    }

    #[test]
    fn a_cloud_without_normals_cannot_be_reconstructed_from_its_own() {
        let mut cloud = Mesh3::create_sphere(2.0, 0.01)
            .expect("sphere")
            .sample_poisson(0.3, None)
            .expect("sampling");
        cloud
            .set_point_normals(None)
            .expect("clearing normals failed");

        let message = match cloud.reconstruct_surface(&ReconstructOpts3::new(0.3)) {
            Ok(_) => panic!("a cloud with no normals should be refused"),
            Err(e) => format!("{e}"),
        };
        assert!(
            message.contains("no normals"),
            "the message should say what is missing, got: {message}"
        );
    }

    #[test]
    fn an_empty_cloud_cannot_be_reconstructed() {
        let cloud = PointCloud3::empty();
        assert!(
            cloud
                .reconstruct_surface(&ReconstructOpts3::new(0.3))
                .is_err()
        );
    }

    #[test]
    fn a_negative_confidence_threshold_is_refused() {
        let cloud = Mesh3::create_sphere(2.0, 0.01)
            .expect("sphere")
            .sample_poisson(0.3, None)
            .expect("sampling");

        let opts = ReconstructOpts3::new(0.3).with_min_normal_confidence(-0.1);
        assert!(cloud.reconstruct_surface(&opts).is_err());
    }

    #[test]
    fn an_estimation_radius_that_cannot_work_is_refused() {
        let cloud = Mesh3::create_sphere(2.0, 0.01)
            .expect("sphere")
            .sample_poisson(0.3, None)
            .expect("sampling");

        let opts = ReconstructOpts3::new(0.3).with_normals(NormalSource3::Estimate {
            radius: 0.0,
            orientation: NormalOrientation3::Propagate { k: 12 },
        });
        assert!(cloud.reconstruct_surface(&opts).is_err());
    }

    // ===========================================================================================
    // Point spacing
    // ===========================================================================================

    /// A Poisson disk sample keeps its points at least one disk radius apart, so the median nearest
    /// neighbor distance must be at or slightly above that radius.
    #[test]
    fn point_spacing_recovers_the_sampling_radius() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");

        for radius in [0.1f64, 0.25, 0.5] {
            let cloud = sphere.sample_poisson(radius, None).expect("sampling");
            let measured = cloud
                .compute_index()
                .expect("index")
                .estimate_point_spacing();

            assert!(
                measured >= radius * 0.98,
                "a disk radius of {radius} came back as a spacing of {measured}"
            );
            assert!(
                measured < radius * 1.5,
                "a disk radius of {radius} came back as a spacing of {measured}, too loose to be \
                 useful for sizing a voxel"
            );
        }
    }

    #[test]
    fn point_spacing_of_a_tiny_cloud_is_zero() {
        let empty = PointCloud3::empty();
        assert_eq!(
            empty
                .compute_index()
                .expect("index")
                .estimate_point_spacing(),
            0.0
        );

        let single = PointCloud3::new(vec![Point3::origin()]);
        assert_eq!(
            single
                .compute_index()
                .expect("index")
                .estimate_point_spacing(),
            0.0
        );
    }

    // ===========================================================================================
    // Known answers
    // ===========================================================================================

    /// A sphere is smooth and closed, with no features for the field to round. The measured
    /// deviations are 0.042 voxels from the reference to the result and 0.041 voxels in the reverse
    /// direction, with one connected piece.
    #[test]
    fn a_sphere_reconstructs_to_within_a_tenth_of_a_voxel() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let opts = ReconstructOpts3::new(spacing);

        let (forward, back, mesh, report) = round_trip(&sphere, spacing, &opts);

        assert!(forward < 0.1, "reference to test was {forward:.4} voxels");
        assert!(back < 0.1, "test to reference was {back:.4} voxels");
        assert_eq!(patch_count(&mesh), 1);
        assert_closed(&mesh, "sphere");

        assert_eq!(report.points_dropped, 0);
        assert!(report.orientation.is_none(), "no orientation was asked for");
        assert!(report.repair.is_some(), "repair runs by default");
        assert_eq!(report.faces, mesh.faces().len());
    }

    /// A cylinder has two sharp circular creases where its caps meet its side, and the field rounds
    /// them. Measured at 0.32 voxels out and 0.45 back.
    ///
    /// The cylinder also produces 28 pieces at this resolution because rounding creates breaks.
    /// Along a crease, the averaged field pulls away from both faces, and the surfaces can fail to
    /// join. See [`ReconstructOpts3::patch_filter`] for why keeping only the largest piece is
    /// unsuitable in this case.
    #[test]
    fn a_cylinder_reconstructs_to_within_half_a_voxel() {
        let cylinder = Mesh3::create_cylinder(3.0, 8.0, 0.005).expect("cylinder");
        let spacing = 0.25;
        let opts = ReconstructOpts3::new(spacing);

        let (forward, back, mesh, _) = round_trip(&cylinder, spacing, &opts);

        assert!(forward < 0.6, "reference to test was {forward:.4} voxels");
        assert!(back < 0.6, "test to reference was {back:.4} voxels");
        assert!(mesh.faces().len() > 10_000);
    }

    /// A box is the most difficult of the three shapes because it consists entirely of sharp edges
    /// and corners. The measured deviations are 0.43 voxels from the reference to the result and
    /// 0.93 voxels in the reverse direction. The reverse deviation is larger because each rounded
    /// edge bows away from the box's flat faces.
    #[test]
    fn a_box_reconstructs_to_within_about_a_voxel() {
        let boxy = Mesh3::create_box(6.0, 5.0, 4.0, false);
        let spacing = 0.25;
        let opts = ReconstructOpts3::new(spacing);

        let (forward, back, mesh, _) = round_trip(&boxy, spacing, &opts);

        assert!(forward < 0.7, "reference to test was {forward:.4} voxels");
        assert!(back < 1.2, "test to reference was {back:.4} voxels");
        assert!(mesh.faces().len() > 5_000);
    }

    /// Reconstruct the same sphere after removing its normals and recovering them through estimation
    /// and propagation. The reconstructed surface must match the result obtained from the sphere's
    /// supplied normals. The measured deviations are 0.044 and 0.043 voxels, compared with 0.042
    /// and 0.041 voxels for the supplied normals.
    #[test]
    fn a_sphere_reconstructs_the_same_from_estimated_normals() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let mut cloud = sphere.sample_poisson(spacing, None).expect("sampling");
        cloud
            .set_point_normals(None)
            .expect("clearing normals failed");

        let opts = ReconstructOpts3::new(spacing).with_normals(NormalSource3::Estimate {
            radius: spacing * 2.5,
            orientation: NormalOrientation3::Propagate { k: 12 },
        });

        let (mesh, report) = cloud
            .reconstruct_surface(&opts)
            .expect("reconstruction failed");
        let accelerated = Mesh3::from_data(mesh.clone(), false).expect("mesh build failed");
        let deviation = sphere
            .measure_surface_deviation(&accelerated, Some(spacing))
            .expect("deviation failed");

        assert!(deviation.reference_to_test / spacing < 0.1);
        assert!(deviation.test_to_reference / spacing < 0.1);
        assert_eq!(patch_count(&mesh), 1);
        assert_closed(&mesh, "estimated sphere");

        let orientation = report.orientation.expect("orientation should be reported");
        assert_eq!(orientation.components, Some(1));
    }

    /// A viewpoint on the far side of a sphere orients only the near half toward itself. This test
    /// verifies the execution path: the call must succeed and report no components.
    #[test]
    fn a_viewpoint_drives_the_estimation_path() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let mut cloud = sphere.sample_poisson(spacing, None).expect("sampling");
        cloud
            .set_point_normals(None)
            .expect("clearing normals failed");

        let opts = ReconstructOpts3::new(spacing).with_normals(NormalSource3::Estimate {
            radius: spacing * 2.5,
            orientation: NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1000.0)),
        });

        let (mesh, report) = cloud
            .reconstruct_surface(&opts)
            .expect("reconstruction failed");
        assert!(mesh.faces().len() > 1000);
        let orientation = report.orientation.expect("orientation should be reported");
        assert_eq!(orientation.components, None);
        assert!(
            orientation.flipped > 0,
            "a scrambled axis set should need flipping"
        );
    }

    /// Reversing every normal reverses the field and produces an inward-facing surface. Downstream
    /// processing does not correct the orientation, so the input normals must be oriented correctly.
    #[test]
    fn reversed_normals_give_an_inward_facing_surface() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let mut cloud = sphere.sample_poisson(spacing, None).expect("sampling");

        let reversed: Vec<UnitVec3> = cloud
            .point_normals()
            .expect("normals")
            .iter()
            .map(|n| -*n)
            .collect();
        cloud
            .set_point_normals(Some(reversed))
            .expect("setting normals failed");

        let (mesh, _) = cloud
            .reconstruct_surface(&ReconstructOpts3::new(spacing))
            .expect("reconstruction failed");

        assert_closed(&mesh, "reversed sphere");

        let mut inward = 0;
        for face in mesh.faces() {
            let p: Vec<Point3> = face.iter().map(|&i| mesh.points()[i as usize]).collect();
            let normal = (p[1] - p[0]).cross(&(p[2] - p[0]));
            let centroid = (p[0].coords + p[1].coords + p[2].coords) / 3.0;
            if normal.dot(&centroid) < 0.0 {
                inward += 1;
            }
        }

        assert_eq!(
            inward,
            mesh.faces().len(),
            "every face should face the middle of the sphere"
        );
    }

    // ===========================================================================================
    // Filtering and post-processing
    // ===========================================================================================

    /// Points whose neighborhoods are too sparse to define a plane receive zero normal confidence.
    /// A positive threshold must remove these points from the field.
    #[test]
    fn a_confidence_threshold_drops_points_with_nothing_around_them() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let cloud = sphere.sample_poisson(spacing, None).expect("sampling");

        // Place several points far from the sphere and each other. Each point is alone within the
        // estimation radius and cannot support a plane fit.
        let mut points = cloud.points().to_vec();
        let strays = 5;
        for i in 0..strays {
            points.push(Point3::new(40.0 + i as f64 * 10.0, 0.0, 0.0));
        }
        let widened = PointCloud3::new(points);

        let opts = ReconstructOpts3::new(spacing)
            .with_normals(NormalSource3::Estimate {
                radius: spacing * 2.5,
                orientation: NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1000.0)),
            })
            .with_min_normal_confidence(0.01);

        let (_, report) = widened
            .reconstruct_surface(&opts)
            .expect("reconstruction failed");

        assert_eq!(
            report.points_dropped, strays,
            "the isolated points should have been the ones dropped"
        );
        assert_eq!(report.points_used, cloud.points().len());
    }

    #[test]
    fn a_confidence_threshold_which_takes_everything_is_an_error() {
        let cloud = Mesh3::create_sphere(5.0, 0.005)
            .expect("sphere")
            .sample_poisson(0.25, None)
            .expect("sampling");

        let opts = ReconstructOpts3::new(0.25)
            .with_normals(NormalSource3::Estimate {
                radius: 0.6,
                orientation: NormalOrientation3::Viewpoint(Point3::new(0.0, 0.0, 1000.0)),
            })
            .with_min_normal_confidence(2.0);

        assert!(cloud.reconstruct_surface(&opts).is_err());
    }

    #[test]
    fn a_patch_filter_leaves_one_piece() {
        let bunny = crate::tests::stanford_bun_3();
        let spacing = 0.003;
        let cloud = bunny.sample_dense(spacing, None).expect("sampling");

        let plain = ReconstructOpts3::new(spacing);
        let (before, _) = cloud
            .reconstruct_surface(&plain)
            .expect("reconstruction failed");
        assert!(
            patch_count(&before) > 1,
            "the bunny should come back in several pieces at this resolution"
        );

        let filtered = plain
            .clone()
            .with_patch_filter(Some(PatchFilter::keep_largest()));
        let (after, _) = cloud
            .reconstruct_surface(&filtered)
            .expect("reconstruction failed");
        assert_eq!(patch_count(&after), 1);
        assert!(after.faces().len() < before.faces().len());
    }

    #[test]
    fn smoothing_runs_and_keeps_the_surface_where_it_was() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;

        let plain = ReconstructOpts3::new(spacing);
        let smoothed = plain.clone().with_smooth_iterations(2);

        let (_, _, rough_mesh, _) = round_trip(&sphere, spacing, &plain);
        let (forward, back, smooth_mesh, _) = round_trip(&sphere, spacing, &smoothed);

        assert!(!smooth_mesh.faces().is_empty());
        assert!(
            forward < 0.2 && back < 0.2,
            "smoothing moved the surface to {forward:.4} and {back:.4} voxels"
        );

        // Smoothing moves geometry without adding or removing it.
        assert_eq!(smooth_mesh.faces().len(), rough_mesh.faces().len());
    }

    /// Skipping repair must preserve the raw extracted mesh.
    #[test]
    fn repair_can_be_skipped() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let spacing = 0.25;
        let opts = ReconstructOpts3::new(spacing).with_repair(None);

        let (_, _, mesh, report) = round_trip(&sphere, spacing, &opts);

        assert!(report.repair.is_none());
        assert_eq!(mesh.faces().len(), report.raw_faces);
    }

    // ===========================================================================================
    // Verify the extraction guarantee required by the default repair preset
    // ===========================================================================================

    #[test]
    fn the_default_repair_omits_the_edge_topology_passes() {
        let opts = ReconstructOpts3::new(1.0);
        assert_eq!(opts.repair, Some(RepairOpts::assuming_oriented_edges()));
    }

    /// The default repair skips two edge-topology passes because `extract_isosurface` cannot produce
    /// the defects that they detect. Run the full repair on each fixture's raw extraction and verify
    /// that these passes find no defects, providing evidence for the default.
    ///
    /// If this test fails, restore [`RepairOpts::default`] as the reconstruction default. The cases
    /// deliberately include sharp-edged and smooth shapes because sharp edges place greater demands
    /// on extraction.
    #[test]
    fn the_skipped_repair_passes_find_nothing_to_fix() {
        let sphere = Mesh3::create_sphere(5.0, 0.005).expect("sphere");
        let cylinder = Mesh3::create_cylinder(3.0, 8.0, 0.005).expect("cylinder");
        let boxy = Mesh3::create_box(6.0, 5.0, 4.0, false);
        let bunny = crate::tests::stanford_bun_3();

        let cases: Vec<(&str, PointCloud3, f64)> = vec![
            (
                "sphere/0.25",
                sphere.sample_poisson(0.25, None).expect("sampling"),
                0.25,
            ),
            (
                "sphere/0.12",
                sphere.sample_poisson(0.12, None).expect("sampling"),
                0.12,
            ),
            (
                "cylinder",
                cylinder.sample_poisson(0.25, None).expect("sampling"),
                0.25,
            ),
            (
                "box",
                boxy.sample_poisson(0.25, None).expect("sampling"),
                0.25,
            ),
            (
                "bunny",
                bunny.sample_dense(0.003, None).expect("sampling"),
                0.003,
            ),
        ];

        for (name, cloud, voxel_size) in cases.iter() {
            // `repair: None` returns the extraction unchanged, which is the input that repair would
            // otherwise receive.
            let raw_opts = ReconstructOpts3::new(*voxel_size).with_repair(None);
            let (raw, _) = cloud
                .reconstruct_surface(&raw_opts)
                .expect("reconstruction failed");

            let full = repair_buffers(raw.points(), raw.faces(), &RepairOpts::default())
                .expect("full repair failed");

            assert_eq!(
                full.report.nonmanifold_edges, 0,
                "{name}: an edge was shared by more than two faces"
            );
            assert_eq!(
                full.report.faces_dropped_at_nonmanifold, 0,
                "{name}: faces had to be dropped to make the edges manifold"
            );
            assert_eq!(
                full.report.faces_reoriented, 0,
                "{name}: faces had to be flipped to agree with their neighbors"
            );
            assert_eq!(
                full.report.faces_dropped_for_orientation, 0,
                "{name}: faces had to be dropped to orient the surface"
            );
            assert_eq!(
                full.report.nonorientable_components, 0,
                "{name}: a component could not be oriented"
            );

            // Because the omitted passes find no defects, the less expensive preset must produce
            // the same mesh.
            let cheap = repair_buffers(
                raw.points(),
                raw.faces(),
                &RepairOpts::assuming_oriented_edges(),
            )
            .expect("cheap repair failed");

            assert_eq!(
                cheap.points.len(),
                full.points.len(),
                "{name}: the two repairs disagree on the point count"
            );
            assert_eq!(
                cheap.faces, full.faces,
                "{name}: the two repairs produced different faces"
            );
        }
    }

    /// Verify that the retained passes detect defects. Sharp-edged shapes exercise these passes.
    #[test]
    fn the_kept_repair_passes_have_work_to_do() {
        let boxy = Mesh3::create_box(6.0, 5.0, 4.0, false);
        let cloud = boxy.sample_poisson(0.25, None).expect("sampling");

        let (mesh, report) = cloud
            .reconstruct_surface(&ReconstructOpts3::new(0.25))
            .expect("reconstruction failed");

        let repair = report.repair.expect("repair should have run");

        // A field value of exactly zero at a grid corner places two vertices at the same position,
        // producing zero-area faces. This fixture originally contained 2,518 such faces.
        assert!(
            repair.degenerate_removed > 100,
            "only {} degenerate faces were found, so this fixture no longer exercises the pass",
            repair.degenerate_removed
        );

        // Marching cubes does not guarantee vertex topology, and a box has the ambiguous faces that
        // produce bowtie vertices. This fixture originally contained 90 of them.
        assert!(
            repair.bowtie_vertices_split > 10,
            "only {} bowtie vertices were split, so this fixture no longer exercises the pass",
            repair.bowtie_vertices_split
        );

        assert!(!mesh.faces().is_empty());
    }

    // ===========================================================================================
    // A real scanned shape
    // ===========================================================================================

    /// The bunny includes thin features, creases, and a flat sawn base that the primitives lack.
    ///
    /// Measured at 0.91 voxels out and 1.76 back at this resolution, in 8 pieces. The outward
    /// figure is the larger one because the reconstruction bridges between the thin parts of the
    /// ears, where the two sides of the surface are closer together than the band is wide.
    #[test]
    fn the_stanford_bunny_reconstructs_within_a_couple_of_voxels() {
        let bunny = crate::tests::stanford_bun_3();
        let spacing = 0.003;
        let cloud = bunny.sample_dense(spacing, None).expect("sampling");

        let opts = ReconstructOpts3::new(spacing);
        let (mesh, report) = cloud
            .reconstruct_surface(&opts)
            .expect("reconstruction failed");

        let accelerated = Mesh3::from_data(mesh.clone(), false).expect("mesh build failed");
        let deviation = bunny
            .measure_surface_deviation(&accelerated, Some(spacing))
            .expect("deviation failed");

        let forward = deviation.reference_to_test / spacing;
        let back = deviation.test_to_reference / spacing;

        assert!(forward < 1.3, "reference to test was {forward:.4} voxels");
        assert!(back < 2.2, "test to reference was {back:.4} voxels");
        assert!(mesh.faces().len() > 20_000);

        assert_eq!(report.points_used, cloud.points().len());
        assert!(report.cells_visited > 0);
        assert!(
            report.cells_skipped_unknown > 0,
            "an open scan should leave a fringe of skipped cells"
        );
    }
}
