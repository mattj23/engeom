//! Quality benchmark for the multi-mesh alignment, run against private laser scan data.
//!
//! The data is a set of samples, each one a part scanned from 15 to 20 directions by a laser
//! profile sensor, together with a `reference.ply` produced by scanning the same physical part on
//! a Zeiss ATOS 5 structured light system. The ATOS mesh is lower resolution and captures less
//! detail, but its bulk shape is the best available source of truth, so it is what a candidate
//! alignment gets measured against.
//!
//! This file establishes that the data resolves and loads and that the stored transforms place
//! their scans where they claim to, builds and validates the frozen measurement set that every
//! candidate alignment is scored against, and then holds the scorer itself along with the
//! regression gates for the baseline candidates and the shipping solver.

#![cfg(all(feature = "private_tests", feature = "ply"))]

mod common;

use crate::common::PathPair;
use engeom::common::DistMode;
use engeom::common::kd_tree::KdTreeSearch;
use engeom::geom3::MultiOutcome3;
use engeom::geom3::XyzQuat;
use engeom::geom3::align3::jacobian::point_surf_jacobian;
use engeom::geom3::align3::{
    AlignInfo3, AlignMesh, AlignOptions3, AlignOrigin3, AlignParams3, AlignStorage3, AlignValues3,
    Dof6, GAPParams, MeshAlignPoint, MultiOptions3, generate_alignment_points,
    multi_mesh_adjustment, multi_mesh_adjustment_with_points, points_to_surface3,
};
use engeom::rayon::prelude::*;
use engeom::{Iso3, KdTree3, Mesh3, Point3, Result, SurfacePoint3, UnitVec3};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// The folder, under the private test data root, holding the multi-mesh alignment samples.
const TEST_DATA_FOLDER: &str = "private-multi-align";

/// The mesh of the CAD model the part was made to. Used as the alignment target in the per-scan
/// step that precedes the simultaneous alignment.
const CAD_FILE: &str = "cad.ply";

/// The independently scanned reference mesh, already best-fit to `CAD_FILE` in Zeiss Inspect.
const REFERENCE_FILE: &str = "reference.ply";

/// None of these meshes are closed volumes. The scans are open surface patches by construction,
/// and the CAD mesh is treated the same way so that every distance query in this benchmark carries
/// the same sign convention.
const IS_SOLID: bool = false;

// ================================================================================================
// Data layout
// ================================================================================================

/// The three files belonging to one laser scan, keyed by the five digit prefix they share.
struct ScanFiles {
    id: String,

    /// The downsampled scan mesh, in the original laser scanner coordinate system.
    mesh: PathBuf,

    /// A static transform associated with the robot pose that held the part, bringing the scan
    /// close enough to the CAD model for a local alignment to converge.
    pre_align: PathBuf,

    /// The result of the multi-mesh alignment as it stood around engeom v0.2.16, collapsed into a
    /// single transform from the scanner frame to the CAD frame.
    cad_align: PathBuf,
}

impl ScanFiles {
    fn mesh(&self) -> Result<Mesh3> {
        Mesh3::load_ply(&self.mesh, IS_SOLID)
    }

    fn pre_align(&self) -> Result<Iso3> {
        load_transform(&self.pre_align)
    }

    fn cad_align(&self) -> Result<Iso3> {
        load_transform(&self.cad_align)
    }
}

/// One sample folder: a part, its reference scan, and the individual laser scans of it.
struct SampleCase {
    name: String,
    dir: PathPair,
    scans: Vec<ScanFiles>,
}

impl SampleCase {
    fn cad(&self) -> Result<Mesh3> {
        Mesh3::load_ply(&self.dir.data().join(CAD_FILE), IS_SOLID)
    }

    fn reference(&self) -> Result<Mesh3> {
        Mesh3::load_ply(&self.dir.data().join(REFERENCE_FILE), IS_SOLID)
    }
}

/// Reads one of the transform files, which are written in the labeled quaternion form that
/// [`XyzQuat`] describes rather than in `nalgebra`'s native serialization.
fn load_transform(path: &Path) -> Result<Iso3> {
    let record: XyzQuat = serde_json::from_reader(std::fs::File::open(path)?)?;
    if !record.is_valid() {
        return Err(format!("{} does not hold a usable rotation", path.display()).into());
    }
    Ok(Iso3::from(&record))
}

/// Discovers the sample folders and the scans within them.
///
/// Scans are found by their file extension rather than from a manifest, with the two whole-part
/// meshes excluded by name. A `.ply` without both of its transform files is an error rather than a
/// skip: a partially populated sample would silently shrink the benchmark.
fn find_cases() -> Result<Vec<SampleCase>> {
    let root = common::find_private_test_data()?.new_joined(TEST_DATA_FOLDER)?;

    let mut sample_dirs = std::fs::read_dir(root.data())?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().is_dir())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    sample_dirs.sort();

    let mut cases = Vec::new();
    for name in sample_dirs {
        let dir = root.new_joined(&name)?;
        let data = dir.data();

        let mut ids = std::fs::read_dir(&data)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().is_some_and(|x| x == "ply"))
            .filter_map(|p| p.file_stem().map(|s| s.to_string_lossy().into_owned()))
            .filter(|stem| stem != "cad" && stem != "reference")
            .collect::<Vec<_>>();
        ids.sort();

        let mut scans = Vec::new();
        for id in ids {
            let scan = ScanFiles {
                mesh: data.join(format!("{id}.ply")),
                pre_align: data.join(format!("{id}.pre-align.json")),
                cad_align: data.join(format!("{id}.cad-align.json")),
                id,
            };

            for path in [&scan.pre_align, &scan.cad_align] {
                if !path.exists() {
                    return Err(format!(
                        "Sample {} scan {} is missing {}",
                        name,
                        scan.id,
                        path.display()
                    )
                    .into());
                }
            }

            scans.push(scan);
        }

        cases.push(SampleCase { name, dir, scans });
    }

    Ok(cases)
}

// ================================================================================================
// Tests
// ================================================================================================

/// Confirms the private data resolves through `find_private_test_data` and that every mesh in it
/// loads.
///
/// This is deliberately about the data rather than the algorithm. It is the tripwire for the two
/// ways the benchmark can be wrong before it computes anything: the data root not being found at
/// all, and a `.ply` that engeom cannot read.
#[test]
fn private_data_loads() -> Result<()> {
    let cases = find_cases()?;
    assert!(
        !cases.is_empty(),
        "No sample folders found under {TEST_DATA_FOLDER}"
    );

    for case in cases.iter() {
        let cad = case.cad()?;
        let reference = case.reference()?;

        assert!(
            !cad.faces().is_empty(),
            "Sample {} has an empty CAD mesh",
            case.name
        );
        assert!(
            !reference.faces().is_empty(),
            "Sample {} has an empty reference mesh",
            case.name
        );

        assert!(
            case.scans.len() >= 2,
            "Sample {} has {} scans, which is too few to align simultaneously",
            case.name,
            case.scans.len()
        );

        for scan in case.scans.iter() {
            let mesh = scan.mesh()?;
            assert!(
                !mesh.faces().is_empty(),
                "Sample {} scan {} loaded with no faces",
                case.name,
                scan.id
            );
        }
    }

    Ok(())
}

/// Takes an evenly strided subset of a mesh's vertices, transformed into a common frame.
///
/// Striding the vertex buffer is enough for a coarse positional check and costs nothing; the
/// benchmark proper uses Poisson sampling so that the measurement points are spread by distance
/// rather than by however the mesh happens to be ordered.
fn strided_points(mesh: &Mesh3, iso: &Iso3, count: usize) -> Vec<Point3> {
    let stride = (mesh.points().len() / count).max(1);
    mesh.points()
        .iter()
        .step_by(stride)
        .map(|p| iso * p)
        .collect()
}

/// Returns the median of the absolute deviations of `points` from `target`.
fn median_abs_deviation(target: &Mesh3, points: &[Point3]) -> f64 {
    let mut devs = target
        .measure_deviations(points, DistMode::ToPoint)
        .into_iter()
        .map(f64::abs)
        .collect::<Vec<_>>();

    devs.sort_by(f64::total_cmp);
    devs[devs.len() / 2]
}

/// Confirms that the stored transforms actually place their scans where they claim to.
///
/// This is the check that the [`XyzQuat`] field mapping agrees with whatever wrote these files.
/// Reading the quaternion components in the wrong order still yields a unit quaternion and a
/// perfectly valid rotation, so no amount of inspecting the numbers can catch it. Applying the
/// transform can: a mis-ordered rotation throws the scan tens of millimetres off a part whose
/// whole bounding box is about 127 x 104 x 56 mm.
///
/// The thresholds are deliberately loose. They are here to separate "correctly placed" from
/// "somewhere else entirely", not to measure alignment quality.
#[test]
fn stored_transforms_place_their_scans() -> Result<()> {
    for case in find_cases()?.iter() {
        let cad = case.cad()?;

        for scan in case.scans.iter() {
            let mesh = scan.mesh()?;

            let pre = median_abs_deviation(&cad, &strided_points(&mesh, &scan.pre_align()?, 2000));
            let post = median_abs_deviation(&cad, &strided_points(&mesh, &scan.cad_align()?, 2000));

            // The pre-alignment is a static robot pose, only ever meant to land close enough for a
            // local alignment to converge.
            assert!(
                pre < 5.0,
                "Sample {} scan {}: pre-alignment leaves a median deviation of {:.3} mm from the \
                 CAD model, which is too far to be a prealignment",
                case.name,
                scan.id,
                pre
            );

            // The cad-alignment is the finished v0.2.16 answer, so it should sit on the part.
            assert!(
                post < 0.5,
                "Sample {} scan {}: cad-alignment leaves a median deviation of {:.3} mm from the \
                 CAD model, which is too far to be an alignment result",
                case.name,
                scan.id,
                post
            );

            // The finished alignment must be an improvement on the pose it started from.
            assert!(
                post < pre,
                "Sample {} scan {}: cad-alignment ({:.3} mm) is no better than the pre-alignment \
                 ({:.3} mm)",
                case.name,
                scan.id,
                post,
                pre
            );
        }
    }

    Ok(())
}

// ================================================================================================
// The frozen benchmark
// ================================================================================================

/// The spacing, in millimetres, of the frozen measurement points across each scan's surface.
///
/// The measurement points exist to sample the residual field densely enough for stable order
/// statistics, not to resolve fine detail, so this is deliberately coarser than the scan meshes.
const SAMPLE_SPACING: f64 = 1.0;

/// The largest distance, in millimetres, at which a measurement point is considered to have landed
/// on a surface another scan also saw.
///
/// The reference configuration this is computed in is a finished alignment, so a genuinely covered
/// point sits far closer than this; the gate is here to reject points over a hole or beyond the
/// edge of the other scan's coverage.
const COVERAGE_MAX_DISTANCE: f64 = 1.0;

/// The largest angle, in radians, between a measurement point's normal and its match's normal for
/// the two to be considered the same piece of surface. This rejects matches through a thin wall,
/// where the geometry is close but facing the wrong way.
const COVERAGE_MAX_NORMAL_ANGLE: f64 = PI_6;

/// Thirty degrees, the normal-agreement gate used both for coverage and in the solver trials.
const PI_6: f64 = std::f64::consts::PI / 6.0;

/// How far, in millimetres, a match must sit from the boundary of the scan it matched against.
///
/// This is the gate that matters most. A measurement point just past the edge of another scan's
/// coverage still projects onto that scan's border and yields a small, entirely meaningless
/// deviation. Without this, the benchmark would quietly reward alignments that slide scans off one
/// another.
const COVERAGE_EDGE_MARGIN: f64 = 1.0;

/// The fewest measurement points a scan may have on the reference mesh and still be fitted to it.
///
/// Six parameters need far fewer than this in principle; the bar is set high because a scan with
/// only a handful of points on the reference has essentially no overlap with it, and a fit from
/// that is not a floor worth quoting.
const MIN_FIT_POINTS: usize = 100;

/// The name of the frozen artifact within each sample folder.
const BENCHMARK_FILE: &str = "benchmark.bin";

const BENCHMARK_MAGIC: &[u8; 4] = b"EGMB";

/// The artifact layout version.
///
/// A mismatch is refused outright rather than migrated. The artifact is cheap to regenerate and
/// is only ever consumed by this file, so silently reading an older layout would buy nothing and
/// risk misinterpreting the bytes as something they are not.
const BENCHMARK_VERSION: u32 = 2;

/// One frozen measurement point, in the coordinate frame of the scan it was sampled from.
#[derive(Clone, Debug)]
struct BenchPoint {
    point: Point3,
    normal: UnitVec3,

    /// The scans this point should be measured against, frozen from the reference configuration.
    ///
    /// Empty is legitimate: a point in a region only one scan saw contributes nothing to mutual
    /// consistency, but may still count toward reference fidelity.
    partners: Vec<u32>,

    /// Whether this point lands on a part of `reference.ply` that the ATOS scanner actually
    /// captured, frozen from the same reference configuration and by the same rules as `partners`.
    ///
    /// The reference mesh has its own holes and borders, and a point projecting past one of them
    /// lands on the border and reports a small, meaningless deviation, exactly as it would against
    /// another scan.
    on_reference: bool,
}

/// The frozen measurement points belonging to one scan.
#[derive(Clone, Debug)]
struct ScanBench {
    id: String,

    /// The vertex and face counts of the mesh these points were sampled from. The points are
    /// stored as raw coordinates, so they do not depend on face ordering, but they do assume the
    /// mesh has not changed underneath them. These are the tripwire for that.
    vertex_count: u32,
    face_count: u32,

    /// The transform placing this scan on `reference.ply` by an individual six degree of freedom
    /// best fit, unconstrained by any of the other scans.
    ///
    /// This is the floor: the closest this scan can possibly come to the reference, and therefore
    /// the best any multi-mesh alignment could achieve on it. A joint alignment must additionally
    /// satisfy mutual consistency, which this fit is free to ignore.
    ///
    /// It is a fixed property of the data rather than of any candidate, which is why it is frozen
    /// here rather than recomputed. Storing the transform rather than a summary statistic is what
    /// lets the floor be re-scored later as a candidate in its own right.
    floor: Iso3,

    points: Vec<BenchPoint>,
}

/// The frozen measurement set for one sample.
///
/// Everything here is chosen once and never recomputed. That is the whole point: a candidate
/// alignment which was allowed to choose its own measurement points, or to decide for itself which
/// of them are in overlap, could improve its score by discarding the points it does badly on.
#[derive(Clone, Debug)]
struct Benchmark {
    spacing: f64,
    max_distance: f64,
    max_normal_angle: f64,
    edge_margin: f64,
    scans: Vec<ScanBench>,
}

impl Benchmark {
    fn point_count(&self) -> usize {
        self.scans.iter().map(|s| s.points.len()).sum()
    }

    fn pair_count(&self) -> usize {
        self.scans
            .iter()
            .flat_map(|s| s.points.iter())
            .map(|p| p.partners.len())
            .sum()
    }
}

// ================================================================================================
// Binary format
// ================================================================================================
//
// Layout, all values little-endian, following the style of `engeom::io::write_mesh_binary`:
//
// - 4 bytes magic: `b"EGMB"`
// - 1 x `u32`: layout version
// - 4 x `f64`: spacing, max distance, max normal angle, edge margin
// - 1 x `u32`: scan count
// - per scan:
//   - 1 x `u32`: byte length of the id, then that many bytes of UTF-8
//   - 2 x `u32`: source mesh vertex count, face count
//   - 7 x `f64`: the floor transform, in the field order of `XyzQuat`
//   - 1 x `u32`: measurement point count
//   - per point: 3 x `f32` position, 3 x `f32` normal, 1 x `u8` reference coverage flag,
//     1 x `u32` partner count, that many `u32`
//
// Positions are stored as `f32`. At the ~100 mm scale of these parts that is a quantization of
// well under a hundredth of a micron, which is four orders of magnitude below the deviations being
// measured, and it is identical for every candidate scored against the artifact.

fn write_benchmark<W: Write>(writer: &mut W, bench: &Benchmark) -> Result<()> {
    writer.write_all(BENCHMARK_MAGIC)?;
    writer.write_all(&BENCHMARK_VERSION.to_le_bytes())?;
    for v in [
        bench.spacing,
        bench.max_distance,
        bench.max_normal_angle,
        bench.edge_margin,
    ] {
        writer.write_all(&v.to_le_bytes())?;
    }

    writer.write_all(&(bench.scans.len() as u32).to_le_bytes())?;
    for scan in bench.scans.iter() {
        writer.write_all(&(scan.id.len() as u32).to_le_bytes())?;
        writer.write_all(scan.id.as_bytes())?;
        writer.write_all(&scan.vertex_count.to_le_bytes())?;
        writer.write_all(&scan.face_count.to_le_bytes())?;

        for v in XyzQuat::from(&scan.floor).to_array() {
            writer.write_all(&v.to_le_bytes())?;
        }

        writer.write_all(&(scan.points.len() as u32).to_le_bytes())?;
        for p in scan.points.iter() {
            for v in [p.point.x, p.point.y, p.point.z] {
                writer.write_all(&(v as f32).to_le_bytes())?;
            }
            for v in [p.normal.x, p.normal.y, p.normal.z] {
                writer.write_all(&(v as f32).to_le_bytes())?;
            }
            writer.write_all(&[u8::from(p.on_reference)])?;
            writer.write_all(&(p.partners.len() as u32).to_le_bytes())?;
            for j in p.partners.iter() {
                writer.write_all(&j.to_le_bytes())?;
            }
        }
    }

    Ok(())
}

fn read_benchmark<R: Read>(reader: &mut R) -> Result<Benchmark> {
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes)?;
    let mut r = ByteReader::new(&bytes);

    if &r.read_bytes::<4>()? != BENCHMARK_MAGIC {
        return Err("Not a benchmark file: invalid magic bytes".into());
    }

    let version = r.read_u32()?;
    if version != BENCHMARK_VERSION {
        return Err(format!(
            "Benchmark file is layout version {version}, but this build expects \
             {BENCHMARK_VERSION}; run the regenerate_benchmark_artifacts producer"
        )
        .into());
    }

    let spacing = r.read_f64()?;
    let max_distance = r.read_f64()?;
    let max_normal_angle = r.read_f64()?;
    let edge_margin = r.read_f64()?;

    let scan_count = r.read_u32()? as usize;
    let mut scans = Vec::with_capacity(scan_count);
    for _ in 0..scan_count {
        let id_len = r.read_u32()? as usize;
        let id = String::from_utf8(r.read_slice(id_len)?.to_vec())?;
        let vertex_count = r.read_u32()?;
        let face_count = r.read_u32()?;

        let mut floor = [0.0_f64; 7];
        for v in floor.iter_mut() {
            *v = r.read_f64()?;
        }
        let floor = Iso3::from(&XyzQuat::from(floor));

        let point_count = r.read_u32()? as usize;
        let mut points = Vec::with_capacity(point_count);
        for _ in 0..point_count {
            let point = Point3::new(
                f64::from(r.read_f32()?),
                f64::from(r.read_f32()?),
                f64::from(r.read_f32()?),
            );
            let normal = UnitVec3::new_normalize(engeom::Vector3::new(
                f64::from(r.read_f32()?),
                f64::from(r.read_f32()?),
                f64::from(r.read_f32()?),
            ));

            let on_reference = r.read_bytes::<1>()?[0] != 0;

            let partner_count = r.read_u32()? as usize;
            let mut partners = Vec::with_capacity(partner_count);
            for _ in 0..partner_count {
                partners.push(r.read_u32()?);
            }

            points.push(BenchPoint {
                point,
                normal,
                partners,
                on_reference,
            });
        }

        scans.push(ScanBench {
            id,
            vertex_count,
            face_count,
            floor,
            points,
        });
    }

    Ok(Benchmark {
        spacing,
        max_distance,
        max_normal_angle,
        edge_margin,
        scans,
    })
}

/// A bounds-checked cursor over the artifact bytes.
///
/// The library's own binary readers index directly and panic on a short file. This one returns an
/// error instead, because the artifact is loaded from a directory a person maintains by hand and a
/// truncated file should say so.
struct ByteReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn read_slice(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.offset + n > self.bytes.len() {
            return Err("Benchmark file ended unexpectedly".into());
        }
        let value = &self.bytes[self.offset..self.offset + n];
        self.offset += n;
        Ok(value)
    }

    fn read_bytes<const N: usize>(&mut self) -> Result<[u8; N]> {
        Ok(self.read_slice(N)?.try_into()?)
    }

    fn read_f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.read_bytes::<4>()?))
    }

    fn read_f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.read_bytes::<8>()?))
    }

    fn read_u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.read_bytes::<4>()?))
    }
}

// ================================================================================================
// The producer
// ================================================================================================

/// A kd-tree over a mesh's boundary vertices, or `None` if the mesh has no boundary at all.
///
/// Testing proximity to boundary vertices rather than to the boundary edges themselves is an
/// approximation, but these meshes triangulate far finer than [`COVERAGE_EDGE_MARGIN`], so a point
/// near an edge is always near one of that edge's endpoints.
fn boundary_tree(mesh: &Mesh3) -> Option<KdTree3> {
    let nav = mesh.compute_nav();
    let verts = nav
        .boundary_vertices(None)
        .into_iter()
        .map(|i| mesh.points()[i as usize])
        .collect::<Vec<_>>();

    if verts.is_empty() {
        None
    } else {
        KdTree3::try_new(&verts).ok()
    }
}

/// Decides whether a query point, expressed in `target`'s own frame, lands on a piece of surface
/// that `target` genuinely captured.
///
/// All three gates matter, but for different reasons. The distance gate rejects a point floating
/// over a hole; the normal gate rejects a match through the far side of a thin wall, where the
/// geometry is close but facing the wrong way; the boundary gate rejects a point just past the
/// edge of the target's coverage, which would otherwise project onto the border and report a
/// meaningless deviation of nearly zero.
fn is_covered(
    target: &Mesh3,
    boundary: Option<&KdTree3>,
    query: &Point3,
    normal: &engeom::Vector3,
) -> bool {
    let match_ = target.surface_closest_to(query);

    if (query - match_.point()).norm() > COVERAGE_MAX_DISTANCE {
        return false;
    }
    if normal.angle(&match_.normal()) > COVERAGE_MAX_NORMAL_ANGLE {
        return false;
    }

    match boundary {
        Some(tree) => tree.nearest_one(&match_.point()).1 >= COVERAGE_EDGE_MARGIN,
        None => true,
    }
}

/// Best fits one scan to the reference mesh on its own, with all six degrees of freedom free.
///
/// `start` is the pose the scan is already in, and the returned transform replaces it. The fit
/// runs on the scan's own frozen measurement points, restricted to those the reference actually
/// covers, so that the floor is measured on the same population as everything else.
///
/// Returns `None` when too few points land on the reference to determine six parameters.
fn fit_to_reference(reference: &Mesh3, points: &[BenchPoint], start: &Iso3) -> Option<Iso3> {
    let fit_points = points
        .iter()
        .filter(|p| p.on_reference)
        .map(|p| start * p.point)
        .collect::<Vec<_>>();

    if fit_points.len() < MIN_FIT_POINTS {
        return None;
    }

    // Rotating about the centroid rather than the world origin keeps the rotation and translation
    // parameters on comparable scales, which matters here because the part sits well off the
    // origin of the CAD frame.
    let centroid = Point3::from(
        fit_points.iter().map(|p| p.coords).sum::<engeom::Vector3>() / fit_points.len() as f64,
    );

    let params = AlignParams3::from_center(centroid, None);
    let outcome =
        points_to_surface3(&fit_points, reference, params, &AlignOptions3::default()).ok()?;

    Some(outcome.alignment().full_transform() * start)
}

/// Builds the frozen measurement set for one sample.
///
/// Coverage is decided in the reference configuration, which is the v0.2.16 alignment, and then
/// never revisited. Note that this deliberately does not use `generate_alignment_points`: that
/// function's filtering depends on the transform it is given, which is exactly the
/// alignment-dependence the frozen mask exists to eliminate.
fn build_benchmark(case: &SampleCase) -> Result<Benchmark> {
    let meshes = case
        .scans
        .iter()
        .map(|s| s.mesh())
        .collect::<Result<Vec<_>>>()?;
    let transforms = case
        .scans
        .iter()
        .map(|s| s.cad_align())
        .collect::<Result<Vec<_>>>()?;

    let boundaries = meshes.par_iter().map(boundary_tree).collect::<Vec<_>>();

    let reference = case.reference()?;
    let reference_boundary = boundary_tree(&reference);

    let samples = meshes
        .par_iter()
        .map(|mesh| mesh.sample_poisson(SAMPLE_SPACING, None))
        .collect::<Vec<_>>();

    let mut scans = Vec::with_capacity(meshes.len());
    for (i, mesh) in meshes.iter().enumerate() {
        // The transform taking a point from scan `i`'s frame into scan `j`'s, in the reference
        // configuration. Computed once per pair rather than once per point.
        let relative = (0..meshes.len())
            .map(|j| (j != i).then(|| transforms[j].inv_mul(&transforms[i])))
            .collect::<Vec<_>>();

        let points = samples[i]
            .par_iter()
            .map(|mp| {
                let partners = relative
                    .iter()
                    .enumerate()
                    .filter_map(|(j, rel)| {
                        let rel = rel.as_ref()?;
                        let covered = is_covered(
                            &meshes[j],
                            boundaries[j].as_ref(),
                            &(rel * mp.point()),
                            &(rel.rotation * mp.normal().into_inner()),
                        );
                        covered.then_some(j as u32)
                    })
                    .collect::<Vec<_>>();

                // The reference mesh is already in the CAD frame, which is where the v0.2.16
                // transforms put the scans, so no relative transform is needed beyond this scan's
                // own.
                let on_reference = is_covered(
                    &reference,
                    reference_boundary.as_ref(),
                    &(transforms[i] * mp.point()),
                    &(transforms[i].rotation * mp.normal().into_inner()),
                );

                BenchPoint {
                    point: mp.point(),
                    normal: mp.normal(),
                    partners,
                    on_reference,
                }
            })
            .collect::<Vec<_>>();

        let floor = fit_to_reference(&reference, &points, &transforms[i]).ok_or_else(|| {
            format!(
                "Scan {} has too few points on the reference mesh to fit against it",
                case.scans[i].id
            )
        })?;

        scans.push(ScanBench {
            id: case.scans[i].id.clone(),
            vertex_count: mesh.points().len() as u32,
            face_count: mesh.faces().len() as u32,
            floor,
            points,
        });
    }

    Ok(Benchmark {
        spacing: SAMPLE_SPACING,
        max_distance: COVERAGE_MAX_DISTANCE,
        max_normal_angle: COVERAGE_MAX_NORMAL_ANGLE,
        edge_margin: COVERAGE_EDGE_MARGIN,
        scans,
    })
}

/// Regenerates the frozen benchmark artifacts and writes them into the sample data folders.
///
/// Ignored by default. This is a producer, not a test: it is run deliberately when the measurement
/// set needs to be rebuilt, and its output is then shipped with the data it describes so that every
/// later run measures in the same places.
///
/// ```text
/// cargo test -r --test multi_mesh_quality_private --features private_tests,ply \
///     -- --ignored --nocapture regenerate_benchmark_artifacts
/// ```
#[test]
#[ignore]
fn regenerate_benchmark_artifacts() -> Result<()> {
    for case in find_cases()?.iter() {
        let bench = build_benchmark(case)?;
        let path = case.dir.data().join(BENCHMARK_FILE);

        let mut file = std::io::BufWriter::new(std::fs::File::create(&path)?);
        write_benchmark(&mut file, &bench)?;
        file.into_inner()?.sync_all()?;

        let bytes = std::fs::metadata(&path)?.len();
        let covered = bench
            .scans
            .iter()
            .flat_map(|s| s.points.iter())
            .filter(|p| !p.partners.is_empty())
            .count();
        let on_ref = bench
            .scans
            .iter()
            .flat_map(|s| s.points.iter())
            .filter(|p| p.on_reference)
            .count();

        println!(
            "{}: {} scans, {} points ({:.1}% paired, {:.1}% on reference), {} pairs, {:.1} MB",
            case.name,
            bench.scans.len(),
            bench.point_count(),
            100.0 * covered as f64 / bench.point_count() as f64,
            100.0 * on_ref as f64 / bench.point_count() as f64,
            bench.pair_count(),
            bytes as f64 / 1.0e6,
        );

        // The per-scan floor: how close each scan gets to the ATOS mesh when fitted to it on its
        // own, with all six degrees of freedom free.
        //
        // Deliberately not printed beside the scan's deviation under its `cad-align.json`
        // transform. That comparison is not meaningful here, because the floor is a fit to
        // `reference.ply` and `cad-align.json` is not: it is a fit to `cad.ply`, and
        // `reference.ply` is separately a fit to `cad.ply`, so composing them leaves a rigid
        // offset of the whole group. Putting the two numbers side by side would charge the
        // alignment for that offset. The scorer removes it, via `gauge_fit`, for every
        // candidate before measuring against the reference.
        let reference = case.reference()?;
        for scan in bench.scans.iter() {
            let pts = scan
                .points
                .iter()
                .filter(|p| p.on_reference)
                .map(|p| scan.floor * p.point)
                .collect::<Vec<_>>();

            println!(
                "    {}: floor {:.4} mm over {} points",
                scan.id,
                median_abs_deviation(&reference, &pts),
                pts.len(),
            );
        }
    }

    Ok(())
}

/// Round trips a synthetic benchmark through the binary format.
///
/// Synthetic rather than generated so that the format is checked without a ten minute producer
/// run, and so that awkward cases the real data may not contain are covered: a point with no
/// partners at all, and a scan with no points.
#[test]
fn benchmark_format_round_trips() -> Result<()> {
    let bench = Benchmark {
        spacing: 1.25,
        max_distance: 0.75,
        max_normal_angle: 0.5,
        edge_margin: 2.0,
        scans: vec![
            ScanBench {
                id: "10001".to_string(),
                vertex_count: 1234,
                face_count: 2345,
                floor: Iso3::translation(1.5, -2.5, 0.75)
                    * Iso3::rotation(engeom::Vector3::new(0.1, -0.2, 0.3)),
                points: vec![
                    BenchPoint {
                        point: Point3::new(1.0, 2.0, 3.0),
                        normal: UnitVec3::new_normalize(engeom::Vector3::new(1.0, 1.0, 0.0)),
                        partners: vec![1, 2],
                        on_reference: true,
                    },
                    BenchPoint {
                        point: Point3::new(-4.5, 0.25, 9.75),
                        normal: UnitVec3::new_normalize(engeom::Vector3::new(0.0, 0.0, -1.0)),
                        partners: vec![],
                        on_reference: false,
                    },
                ],
            },
            ScanBench {
                id: "10002".to_string(),
                vertex_count: 9,
                face_count: 8,
                floor: Iso3::identity(),
                points: vec![],
            },
        ],
    };

    let mut bytes = Vec::new();
    write_benchmark(&mut bytes, &bench)?;
    let back = read_benchmark(&mut bytes.as_slice())?;

    assert_eq!(back.spacing, bench.spacing);
    assert_eq!(back.max_distance, bench.max_distance);
    assert_eq!(back.max_normal_angle, bench.max_normal_angle);
    assert_eq!(back.edge_margin, bench.edge_margin);
    assert_eq!(back.scans.len(), bench.scans.len());

    for (a, b) in bench.scans.iter().zip(back.scans.iter()) {
        assert_eq!(a.id, b.id);
        assert_eq!(a.vertex_count, b.vertex_count);
        assert_eq!(a.face_count, b.face_count);
        assert_eq!(a.points.len(), b.points.len());
        assert!((a.floor.to_matrix() - b.floor.to_matrix()).amax() < 1e-12);

        for (p, q) in a.points.iter().zip(b.points.iter()) {
            // f32 storage, so the tolerance is the quantization rather than a solver epsilon.
            assert!((p.point - q.point).norm() < 1e-5);
            assert!((p.normal.into_inner() - q.normal.into_inner()).norm() < 1e-6);
            assert_eq!(p.partners, q.partners);
            assert_eq!(p.on_reference, q.on_reference);
        }
    }

    Ok(())
}

#[test]
fn a_truncated_benchmark_file_is_rejected() -> Result<()> {
    let bench = Benchmark {
        spacing: 1.0,
        max_distance: 1.0,
        max_normal_angle: 0.5,
        edge_margin: 1.0,
        scans: vec![ScanBench {
            id: "10001".to_string(),
            vertex_count: 1,
            face_count: 1,
            floor: Iso3::identity(),
            points: vec![BenchPoint {
                point: Point3::new(1.0, 2.0, 3.0),
                normal: UnitVec3::new_normalize(engeom::Vector3::new(0.0, 0.0, 1.0)),
                partners: vec![7],
                on_reference: true,
            }],
        }],
    };

    let mut bytes = Vec::new();
    write_benchmark(&mut bytes, &bench)?;

    bytes.truncate(bytes.len() - 4);
    assert!(read_benchmark(&mut bytes.as_slice()).is_err());

    assert!(read_benchmark(&mut b"NOPE0000".as_slice()).is_err());

    Ok(())
}

/// Loads the frozen artifact for every sample and checks it still describes the data it sits
/// beside.
///
/// This is the consumer side of the freeze. The artifact is generated once and then trusted
/// indefinitely, so the ways it can quietly go stale all need to be loud:
///
/// - a source mesh re-exported underneath it, caught by the stored vertex and face counts
/// - a scan added to or removed from the folder, caught by comparing the id lists
/// - the gates or the spacing changed in this file without the artifact being regenerated, caught
///   by comparing the stored parameters against the constants
#[test]
fn benchmark_artifacts_match_their_meshes() -> Result<()> {
    for case in find_cases()?.iter() {
        let path = case.dir.data().join(BENCHMARK_FILE);
        assert!(
            path.exists(),
            "Sample {} has no {}; run the regenerate_benchmark_artifacts producer",
            case.name,
            BENCHMARK_FILE
        );

        let bench = read_benchmark(&mut std::fs::File::open(&path)?)?;

        assert_eq!(
            bench.spacing, SAMPLE_SPACING,
            "Sample {} was generated at a different sample spacing",
            case.name
        );
        assert_eq!(
            bench.max_distance, COVERAGE_MAX_DISTANCE,
            "Sample {} was generated with a different coverage distance gate",
            case.name
        );
        assert_eq!(
            bench.max_normal_angle, COVERAGE_MAX_NORMAL_ANGLE,
            "Sample {} was generated with a different coverage normal gate",
            case.name
        );
        assert_eq!(
            bench.edge_margin, COVERAGE_EDGE_MARGIN,
            "Sample {} was generated with a different coverage edge margin",
            case.name
        );

        let ids = case.scans.iter().map(|s| &s.id).collect::<Vec<_>>();
        let frozen = bench.scans.iter().map(|s| &s.id).collect::<Vec<_>>();
        assert_eq!(
            ids, frozen,
            "Sample {} holds a different set of scans than its benchmark was built from",
            case.name
        );

        for (scan, files) in bench.scans.iter().zip(case.scans.iter()) {
            let mesh = files.mesh()?;

            assert_eq!(
                scan.vertex_count as usize,
                mesh.points().len(),
                "Sample {} scan {} has changed vertex count since its benchmark was built",
                case.name,
                scan.id
            );
            assert_eq!(
                scan.face_count as usize,
                mesh.faces().len(),
                "Sample {} scan {} has changed face count since its benchmark was built",
                case.name,
                scan.id
            );

            assert!(
                !scan.points.is_empty(),
                "Sample {} scan {} has no measurement points",
                case.name,
                scan.id
            );

            for p in scan.points.iter() {
                for &j in p.partners.iter() {
                    assert!(
                        (j as usize) < bench.scans.len(),
                        "Sample {} scan {} references scan {} of {}",
                        case.name,
                        scan.id,
                        j,
                        bench.scans.len()
                    );
                }
            }

            // The measurement points were sampled from this mesh's surface, so they must still sit
            // on it. Only the f32 storage separates them from exactly zero.
            let coords = scan.points.iter().map(|p| p.point).collect::<Vec<_>>();
            let worst = mesh
                .measure_deviations(&coords, DistMode::ToPoint)
                .into_iter()
                .fold(0.0_f64, |acc, d| acc.max(d.abs()));

            assert!(
                worst < 1e-3,
                "Sample {} scan {}: a frozen point sits {:.2e} mm off the mesh it was sampled from",
                case.name,
                scan.id,
                worst
            );
        }
    }

    Ok(())
}

/// Checks that the frozen floor transforms are what they claim to be.
///
/// The floor is an individual six degree of freedom fit of one scan to the reference, started from
/// the v0.2.16 pose. Three things have to hold for it to be usable as a bar:
///
/// - it must be at least as close to the reference as the pose it started from, since it optimized
///   exactly that and had the freedom to stay put
/// - it must be a refinement rather than a relocation, or the fit has slid the scan onto the wrong
///   part of the surface
/// - it must be small in absolute terms, since a floor of hundreds of microns would mean the two
///   scanners disagree so badly that nothing measured against the reference means very much
#[test]
fn floor_transforms_refine_their_starting_pose() -> Result<()> {
    for case in find_cases()?.iter() {
        let bench = read_benchmark(&mut std::fs::File::open(
            case.dir.data().join(BENCHMARK_FILE),
        )?)?;
        let reference = case.reference()?;

        for (scan, files) in bench.scans.iter().zip(case.scans.iter()) {
            let start = files.cad_align()?;

            let on_ref = scan
                .points
                .iter()
                .filter(|p| p.on_reference)
                .map(|p| p.point)
                .collect::<Vec<_>>();

            assert!(
                on_ref.len() >= MIN_FIT_POINTS,
                "Sample {} scan {} has only {} points on the reference",
                case.name,
                scan.id,
                on_ref.len()
            );

            let moved = |iso: &Iso3| on_ref.iter().map(|p| iso * p).collect::<Vec<_>>();
            let floor_dev = median_abs_deviation(&reference, &moved(&scan.floor));
            let start_dev = median_abs_deviation(&reference, &moved(&start));

            assert!(
                floor_dev <= start_dev,
                "Sample {} scan {}: the individual fit ({:.4} mm) is worse than the pose it \
                 started from ({:.4} mm), which it was free to keep",
                case.name,
                scan.id,
                floor_dev,
                start_dev
            );

            assert!(
                floor_dev < 0.1,
                "Sample {} scan {}: the floor is {:.4} mm, far enough that the laser and ATOS \
                 measurements of this part disagree more than any alignment could fix",
                case.name,
                scan.id,
                floor_dev
            );

            // Measured as point motion rather than as a component of the transform, since
            // `full_transform` is taken about a local origin and its translation part alone says
            // little about how far the geometry actually travelled.
            let shift = moved(&scan.floor)
                .iter()
                .zip(moved(&start).iter())
                .fold(0.0_f64, |acc, (a, b)| acc.max((a - b).norm()));

            assert!(
                shift < 1.0,
                "Sample {} scan {}: the individual fit moved the scan by up to {:.3} mm, which is \
                 a relocation rather than a refinement of an already good pose",
                case.name,
                scan.id,
                shift
            );
        }
    }

    Ok(())
}

// ================================================================================================
// The scorer
// ================================================================================================

/// The absolute deviation, in millimetres, beyond which a measurement point is counted as an
/// outlier rather than as a measurement.
///
/// Fixed and geometric, never derived from a solver's weights. Some points simply cannot be
/// reconciled by any alignment, and the count of them is itself a number worth comparing: a
/// candidate that pushes many more points past this bar has done something bad to a region even if
/// its median looks fine.
const OUTLIER_GATE: f64 = 0.1;

/// A summarized distribution of signed deviations.
///
/// Both a robust centre and a tail, deliberately. The median alone would hide a candidate that
/// wrecks five percent of the part; the RMS alone would be dominated by the irreducible tail that
/// every candidate carries.
#[derive(Clone, Copy, Debug, Default)]
struct Stats {
    count: usize,
    rms: f64,
    median_abs: f64,
    p90: f64,
    p95: f64,
    p99: f64,
    max_abs: f64,
    beyond_gate: usize,
}

impl Stats {
    fn from_deviations(devs: &[f64]) -> Self {
        if devs.is_empty() {
            return Self::default();
        }

        let mut abs = devs.iter().map(|d| d.abs()).collect::<Vec<_>>();
        abs.sort_by(f64::total_cmp);

        let at = |q: f64| abs[(((abs.len() - 1) as f64) * q).round() as usize];

        Self {
            count: devs.len(),
            rms: (devs.iter().map(|d| d * d).sum::<f64>() / devs.len() as f64).sqrt(),
            median_abs: at(0.50),
            p90: at(0.90),
            p95: at(0.95),
            p99: at(0.99),
            max_abs: abs[abs.len() - 1],
            beyond_gate: abs.iter().filter(|d| **d > OUTLIER_GATE).count(),
        }
    }

    fn row(&self) -> String {
        format!(
            "n={:<7} rms={:.4} med={:.4} p90={:.4} p95={:.4} p99={:.4} max={:.4} out={}",
            self.count,
            self.rms,
            self.median_abs,
            self.p90,
            self.p95,
            self.p99,
            self.max_abs,
            self.beyond_gate
        )
    }
}

/// What one candidate alignment scored on one sample.
#[derive(Clone, Debug)]
struct Score {
    /// Scan-to-scan disagreement, over the frozen partner mask. The sensitive metric: no reference
    /// term enters it, so it responds directly to registration error.
    consistency: Stats,

    /// Deviation from `reference.ply` after the gauge fit. The guard rail against a solution that
    /// is locally consistent but globally warped.
    reference: Stats,

    /// The same, broken out per scan, which is where a single bad scan becomes visible.
    per_scan_reference: Vec<Stats>,

    /// The rigid motion applied to the whole cloud before measuring against the reference.
    gauge: Iso3,

    /// Raw deviations in a canonical order, so two scores over the same benchmark can be
    /// differenced point by point.
    consistency_devs: Vec<f64>,
    reference_devs: Vec<f64>,
}

impl Score {
    /// The per-point difference in reference deviation against another candidate, scored over the
    /// same benchmark.
    ///
    /// This is what cancels the common-mode disagreement between the laser scans and the ATOS mesh:
    /// where two candidates place a point in nearly the same spot, the shared measurement
    /// difference subtracts out and what is left is the difference the alignments made.
    fn paired_reference_difference(&self, other: &Score) -> Result<Stats> {
        if self.reference_devs.len() != other.reference_devs.len() {
            return Err("Scores were computed over different benchmarks".into());
        }

        Ok(Stats::from_deviations(
            &self
                .reference_devs
                .iter()
                .zip(other.reference_devs.iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>(),
        ))
    }
}

/// Scores one candidate alignment of a sample.
///
/// The signature is the point of this function: it takes a set of per-scan transforms and nothing
/// else. It does not know, and must not learn, what produced them, so the identical code grades the
/// v0.2.16 transforms from disk, the frozen floor, and whatever the solver just returned.
///
/// Two independent choices keep the comparison honest. The measurement points and the coverage mask
/// come from the frozen benchmark, never from the candidate, so a candidate cannot improve its
/// score by discarding the points it does badly on. And every candidate is measured with `ToPoint`
/// deviations, which penalize a point for hanging off the edge of a surface rather than sliding it
/// along the tangent.
///
/// # A coupling worth stating
///
/// The gauge fit uses `points_to_surface3`, so this function is not entirely free of the alignment
/// code it exists to measure. That is deliberate and narrow: it removes six global degrees of
/// freedom that the multi-mesh alignment genuinely does not determine, over a hundred thousand
/// points, from a starting pose already within a tenth of a millimetre. What it must not do is
/// depend on the *multi-mesh* solver, and it does not. The self-test that a global rigid motion
/// leaves both metrics unchanged is what keeps this honest.
fn score(
    bench: &Benchmark,
    scans: &[Mesh3],
    reference: &Mesh3,
    transforms: &[Iso3],
) -> Result<Score> {
    if transforms.len() != bench.scans.len() || scans.len() != bench.scans.len() {
        return Err("The candidate does not have one transform and mesh per frozen scan".into());
    }

    // ------------------------------------------------------------------------------------------
    // Mutual consistency
    // ------------------------------------------------------------------------------------------

    // The work list is built first, in scan-major then point then partner order, so that every
    // candidate emits its deviations in the same positions and two scores can be differenced.
    let mut work = Vec::new();
    for (i, scan) in bench.scans.iter().enumerate() {
        for (k, p) in scan.points.iter().enumerate() {
            for &j in p.partners.iter() {
                work.push((i, k, j as usize));
            }
        }
    }

    // Grouped by ordered pair so the relative transform is formed once per pair rather than once
    // per point, and so the closest point queries batch against a single target mesh.
    let mut groups: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for (idx, &(i, _, j)) in work.iter().enumerate() {
        groups.entry((i, j)).or_default().push(idx);
    }

    let scattered = groups
        .par_iter()
        .flat_map(|(&(i, j), idxs)| {
            let rel = transforms[j].inv_mul(&transforms[i]);
            let pts = idxs
                .iter()
                .map(|&idx| rel * bench.scans[i].points[work[idx].1].point)
                .collect::<Vec<_>>();

            idxs.iter()
                .copied()
                .zip(scans[j].measure_deviations(&pts, DistMode::ToPoint))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut consistency_devs = vec![0.0; work.len()];
    for (idx, d) in scattered {
        consistency_devs[idx] = d;
    }

    // ------------------------------------------------------------------------------------------
    // Reference fidelity
    // ------------------------------------------------------------------------------------------

    let mut cloud = Vec::new();
    let mut scan_spans = Vec::with_capacity(bench.scans.len());
    for (i, scan) in bench.scans.iter().enumerate() {
        let start = cloud.len();
        cloud.extend(
            scan.points
                .iter()
                .filter(|p| p.on_reference)
                .map(|p| transforms[i] * p.point),
        );
        scan_spans.push(start..cloud.len());
    }

    let gauge = gauge_fit(reference, &cloud)?;
    let placed = cloud.par_iter().map(|p| gauge * p).collect::<Vec<_>>();
    let reference_devs = reference.measure_deviations(&placed, DistMode::ToPoint);

    Ok(Score {
        consistency: Stats::from_deviations(&consistency_devs),
        reference: Stats::from_deviations(&reference_devs),
        per_scan_reference: scan_spans
            .into_iter()
            .map(|s| Stats::from_deviations(&reference_devs[s]))
            .collect(),
        gauge,
        consistency_devs,
        reference_devs,
    })
}

/// Finds the single rigid motion placing a candidate's whole point cloud onto the reference mesh.
///
/// This removes the six global degrees of freedom a multi-mesh alignment leaves undetermined. It is
/// not optional and it is not a formality: `reference.ply` is a fit to `cad.ply` and the v0.2.16
/// transforms are separately a fit to `cad.ply`, so composing them leaves a rigid offset of the
/// whole group which was measured at between fifty and a hundred and forty microns. Charging an
/// alignment for that would be charging it for something it neither caused nor can remove.
fn gauge_fit(reference: &Mesh3, cloud: &[Point3]) -> Result<Iso3> {
    let centroid =
        Point3::from(cloud.iter().map(|p| p.coords).sum::<engeom::Vector3>() / cloud.len() as f64);

    let params = AlignParams3::from_center(centroid, None);
    let outcome = points_to_surface3(cloud, reference, params, &AlignOptions3::default())?;

    Ok(*outcome.alignment().full_transform())
}

/// Everything one sample needs to be scored, loaded once.
struct LoadedCase {
    name: String,
    bench: Benchmark,
    scans: Vec<Mesh3>,
    reference: Mesh3,
    cad_align: Vec<Iso3>,
}

impl LoadedCase {
    fn load(case: &SampleCase) -> Result<Self> {
        Ok(Self {
            name: case.name.clone(),
            bench: read_benchmark(&mut std::fs::File::open(
                case.dir.data().join(BENCHMARK_FILE),
            )?)?,
            scans: case
                .scans
                .iter()
                .map(|s| s.mesh())
                .collect::<Result<Vec<_>>>()?,
            reference: case.reference()?,
            cad_align: case
                .scans
                .iter()
                .map(|s| s.cad_align())
                .collect::<Result<Vec<_>>>()?,
        })
    }

    fn floor(&self) -> Vec<Iso3> {
        self.bench.scans.iter().map(|s| s.floor).collect()
    }

    fn score(&self, transforms: &[Iso3]) -> Result<Score> {
        score(&self.bench, &self.scans, &self.reference, transforms)
    }
}

/// Scoring the same transforms twice must give an identically zero paired difference.
///
/// This is the weaker of the two self-tests but it catches a real class of bug: any hidden
/// dependence on iteration order, on hash map traversal, or on a solver that does not land in
/// exactly the same place twice would show up here as a nonzero difference between a candidate and
/// itself.
#[test]
fn scoring_is_repeatable() -> Result<()> {
    let case = LoadedCase::load(&find_cases()?[0])?;

    let a = case.score(&case.cad_align)?;
    let b = case.score(&case.cad_align)?;

    let diff = a.paired_reference_difference(&b)?;
    assert_eq!(
        diff.max_abs, 0.0,
        "Scoring the same transforms twice differed by up to {:.3e} mm",
        diff.max_abs
    );

    assert_eq!(a.consistency_devs, b.consistency_devs);
    assert_eq!(a.reference.count, b.reference.count);

    Ok(())
}

/// Moving every scan together by the same rigid motion must not change either metric.
///
/// This is the self-test that matters. A multi-mesh alignment does not determine the pose of the
/// group as a whole, so a scorer which let that pose leak into its numbers would be grading
/// candidates on something none of them control.
///
/// The two metrics resist this for different reasons, and both need checking. Mutual consistency is
/// invariant by construction, since a common motion cancels out of every relative transform, so a
/// failure there means the frames are being composed wrongly. Reference fidelity is invariant only
/// because the gauge fit absorbs the motion, so a failure there means the gauge fit did not
/// converge or was not given enough freedom.
#[test]
fn a_global_rigid_motion_does_not_change_the_score() -> Result<()> {
    let case = LoadedCase::load(&find_cases()?[0])?;

    let plain = case.score(&case.cad_align)?;

    // Large enough that a scorer which failed to remove it could not possibly be hiding the
    // difference in rounding: two millimetres and about six degrees.
    let motion = Iso3::translation(2.0, -1.5, 0.75)
        * Iso3::rotation(engeom::Vector3::new(0.05, -0.08, 0.03));
    let moved_transforms = case
        .cad_align
        .iter()
        .map(|t| motion * t)
        .collect::<Vec<_>>();
    let moved = case.score(&moved_transforms)?;

    // Algebraically this is exact: `(M T_j)^-1 (M T_i)` reduces to `T_j^-1 T_i`. In floating point
    // it is not, because the composition and inversion are carried out on the moved transforms, so
    // what is checked is that nothing beyond that rounding survives.
    let worst = plain
        .consistency_devs
        .iter()
        .zip(moved.consistency_devs.iter())
        .fold(0.0_f64, |acc, (a, b)| acc.max((a - b).abs()));

    assert!(
        worst < 1e-9,
        "A global rigid motion changed the mutual consistency deviations by up to {worst:.3e} mm, \
         which is far more than the rounding of composing and inverting the moved transforms"
    );

    let diff = moved.paired_reference_difference(&plain)?;
    assert!(
        diff.max_abs < 1e-6,
        "The gauge fit left {:.3e} mm of a global rigid motion in the reference deviations \
         (median {:.3e} mm)",
        diff.max_abs,
        diff.median_abs
    );

    Ok(())
}

// ================================================================================================
// Candidates and reporting
// ================================================================================================

/// One alignment to be graded: a label and one transform per scan.
///
/// Everything the benchmark compares reduces to this. Where the transforms came from, whether from
/// a file on disk or from a solver that just ran, is deliberately not represented.
struct Candidate {
    name: &'static str,
    transforms: Vec<Iso3>,
}

/// A candidate that has been scored.
struct Graded {
    name: &'static str,
    score: Score,
}

/// The registration error implied by a candidate's RMS once the floor is removed in quadrature.
///
/// The floor and the candidate share the laser-to-ATOS measurement difference, and independent
/// error sources add in quadrature, so subtracting that way is what isolates the part of the
/// candidate's deviation that the alignment is actually responsible for. Returns `None` when the
/// candidate is at or below the floor, where the subtraction has nothing left to describe.
fn registration_component(candidate: &Stats, floor: &Stats) -> Option<f64> {
    let excess = candidate.rms * candidate.rms - floor.rms * floor.rms;
    (excess > 0.0).then(|| excess.sqrt())
}

fn print_table(name: &str, design: u64, graded: &[Graded], floor_index: usize) {
    println!("\n=== {name}  [design {design:04x}] ===");
    println!(
        "  {:<10} | {:>7} {:>7} {:>6} | {:>7} {:>7} {:>7} {:>6} | {:>7}",
        "candidate", "cons.med", "p95", "out", "ref.med", "p95", "rms", "out", "reg."
    );

    let floor = &graded[floor_index].score;
    for g in graded.iter() {
        let reg = match registration_component(&g.score.reference, &floor.reference) {
            Some(v) => format!("{v:.4}"),
            None => "  --   ".to_string(),
        };

        println!(
            "  {:<10} | {:>7.4} {:>7.4} {:>6} | {:>7.4} {:>7.4} {:>7.4} {:>6} | {:>7}",
            g.name,
            g.score.consistency.median_abs,
            g.score.consistency.p95,
            g.score.consistency.beyond_gate,
            g.score.reference.median_abs,
            g.score.reference.p95,
            g.score.reference.rms,
            g.score.reference.beyond_gate,
            reg,
        );
    }

    // The full distributions below the summary. The tail percentiles are where a candidate that
    // looks fine at the median gives itself away, so they are printed rather than folded into the
    // table's fixed width.
    for g in graded.iter() {
        let t = g.score.gauge.translation.vector;
        let (i, worst) = g
            .score
            .per_scan_reference
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.median_abs.total_cmp(&b.1.median_abs))
            .expect("a scored candidate always has at least one scan");

        println!(
            "  {:<10} consistency  {}",
            g.name,
            g.score.consistency.row()
        );
        println!("  {:<10} reference    {}", g.name, g.score.reference.row());
        println!(
            "  {:<10} gauge |t| {:.4} mm, {:.4} deg | worst scan #{} {}",
            g.name,
            t.norm(),
            g.score.gauge.rotation.angle().to_degrees(),
            i,
            worst.row(),
        );
    }
}

/// Identifies which part design a sample belongs to, by the bytes of its CAD mesh.
///
/// The four samples are two designs with two instances each. Grouping by this lets the paired
/// samples be checked against each other: they are the same part measured by the same instruments,
/// so a large split between them points at the data or the alignment rather than at the part.
fn design_key(case: &SampleCase) -> Result<u64> {
    use std::hash::{Hash, Hasher};

    let bytes = std::fs::read(case.dir.data().join(CAD_FILE))?;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    bytes.hash(&mut hasher);
    Ok(hasher.finish())
}

/// Scores the two candidates that need no solver at all, and checks the properties that must hold
/// of them.
///
/// These two set the scale everything else is read against. The floor is what the data itself
/// permits, and v0.2.16 is what the previous generation of this library achieved. Neither can
/// regress, because neither is computed by code under development: the floor is frozen in the
/// artifact and v0.2.16 is a file on disk. What this test guards is the *harness*.
#[test]
fn baseline_candidates_score_within_their_bounds() -> Result<()> {
    let cases = find_cases()?;
    let mut floors_by_design: std::collections::HashMap<u64, Vec<(String, f64)>> =
        std::collections::HashMap::new();

    for case in cases.iter() {
        let loaded = LoadedCase::load(case)?;
        let design = design_key(case)?;

        let candidates = [
            Candidate {
                name: "floor",
                transforms: loaded.floor(),
            },
            Candidate {
                name: "v0.2.16",
                transforms: loaded.cad_align.clone(),
            },
        ];

        let graded = candidates
            .iter()
            .map(|c| {
                Ok(Graded {
                    name: c.name,
                    score: loaded.score(&c.transforms)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        print_table(&loaded.name, design, &graded, 0);

        let floor = &graded[0].score;
        let v0216 = &graded[1].score;

        // The floor is an individual per-scan fit to the reference, so nothing constrained to move
        // the scans as one body can beat it on that metric. A candidate that did would mean the
        // harness is measuring the floor and the candidate differently.
        assert!(
            v0216.reference.median_abs >= floor.reference.median_abs,
            "{}: v0.2.16 scored {:.4} mm against the reference, below the individual-fit floor at \
             {:.4} mm, which is not achievable and indicates a harness fault rather than a good \
             alignment",
            loaded.name,
            v0216.reference.median_abs,
            floor.reference.median_abs
        );

        for g in graded.iter() {
            assert!(
                g.score.consistency.count > 0 && g.score.reference.count > 0,
                "{}: candidate {} was scored on no points at all",
                loaded.name,
                g.name
            );

            // The gauge removes the offset between the CAD frame and the reference. It was
            // measured at well under a tenth of a millimetre, so anything approaching a whole
            // millimetre means it is absorbing a real misplacement instead.
            let t = g.score.gauge.translation.vector.norm();
            let r = g.score.gauge.rotation.angle().to_degrees();
            assert!(
                t < 1.0 && r < 1.0,
                "{}: candidate {} needed a gauge of {:.4} mm and {:.4} deg to sit on the \
                 reference, which is a relocation rather than the removal of a frame offset",
                loaded.name,
                g.name,
                t,
                r
            );
        }

        floors_by_design
            .entry(design)
            .or_default()
            .push((loaded.name.clone(), floor.reference.median_abs));
    }

    // Two samples of the same part design, measured on the same two instruments, must agree about
    // how far apart those instruments are. This is the check that the floor is a property of the
    // measurement rather than an artifact of one sample.
    for (design, mut floors) in floors_by_design {
        if floors.len() < 2 {
            continue;
        }
        floors.sort_by(|a, b| a.1.total_cmp(&b.1));

        let (lo_name, lo) = &floors[0];
        let (hi_name, hi) = &floors[floors.len() - 1];
        assert!(
            hi - lo < 0.005,
            "design {design:04x}: the floor is {lo:.4} mm on {lo_name} but {hi:.4} mm on \
             {hi_name}; two instances of the same part should not disagree this much about the \
             difference between the two scanners"
        );
    }

    Ok(())
}

// ================================================================================================
// Solver candidates
// ================================================================================================

/// How much worse than v0.2.16 the solver's mutual consistency may be before it counts as a
/// regression. It was eight to ten percent *better* on all four samples, so this has real margin.
const CONSISTENCY_SLACK: f64 = 1.05;

/// How much worse than v0.2.16 the solver may be against the reference mesh.
///
/// Deliberately loose. Tightening mutual agreement between deformed scans costs reference fidelity,
/// and the worst observed ratio was about 1.25 on a cold start. This is set to catch a collapse,
/// which is what the missing distance gate produced, rather than to police the trade.
const REFERENCE_CEILING: f64 = 1.5;

/// How far a warm start may move its measurement points from the v0.2.16 seed, in millimetres.
/// Observed between 0.098 and 0.164; the ungated defect reached 0.47.
const WARM_DRIFT_BAR: f64 = 0.4;

/// The same for a cold start, which begins further away and so legitimately travels further.
/// Observed between 0.23 and 0.56.
const COLD_DRIFT_BAR: f64 = 1.0;

/// How far apart the warm and cold starts may land on mutual consistency, in millimetres. They
/// agreed to within 0.0001 on every sample.
const WARM_COLD_AGREEMENT: f64 = 0.002;

/// The Poisson spacing, in millimetres, used when the multi-mesh alignment samples its own
/// correspondences.
///
/// Coarser than the benchmark's measurement spacing. The correspondences only have to determine
/// six parameters per body, and the sampling runs over every pair of meshes twice, so this is the
/// dominant cost of the whole run.
const MULTI_MESH_SPACING: f64 = 2.0;

/// The correspondence distance gate handed to the multi-mesh alignment.
///
/// Required rather than optional, and for a reason this benchmark established: without it the
/// initial unweighted solve is dragged most of a millimetre off a correct answer by a handful of
/// correspondences sampled across regions the two scans never both saw. A millimetre is far wider
/// than any real disagreement between these scans, which sit within about twenty microns of each
/// other, and a tighter gate measured no better.
const MULTI_MESH_MAX_DISTANCE: f64 = 1.0;

/// Runs the simultaneous alignment from a given set of starting poses.
fn run_multi_mesh(scans: &[Mesh3], start: &[Iso3]) -> Result<(Vec<Iso3>, MultiOutcome3)> {
    run_multi_mesh_with(scans, start, &MultiOptions3::new(MULTI_MESH_MAX_DISTANCE))
}

fn run_multi_mesh_with(
    scans: &[Mesh3],
    start: &[Iso3],
    opts: &MultiOptions3,
) -> Result<(Vec<Iso3>, MultiOutcome3)> {
    let meshes = scans
        .iter()
        .zip(start.iter())
        .map(|(m, t)| AlignMesh::new(m, None, Some(t), None))
        .collect::<Vec<_>>();

    let outcome = multi_mesh_adjustment(&meshes, opts, &GAPParams::defaults(MULTI_MESH_SPACING))?;

    let transforms = outcome
        .alignments()
        .iter()
        .map(|a| *a.full_transform())
        .collect();

    Ok((transforms, outcome))
}

/// Individually aligns every scan to the CAD model, starting from its robot-pose prealignment.
///
/// This is the first half of the chain that produced the v0.2.16 numbers: the prealignment only
/// gets a scan near the CAD, and a local fit takes it the rest of the way. It runs on the frozen
/// measurement points, all of them rather than only those the ATOS reference covers, since the CAD
/// model spans the whole part.
fn align_each_to_cad(loaded: &LoadedCase, cad: &Mesh3, pre: &[Iso3]) -> Result<Vec<Iso3>> {
    loaded
        .bench
        .scans
        .iter()
        .zip(pre.iter())
        .map(|(scan, start)| {
            let pts = scan
                .points
                .iter()
                .map(|p| start * p.point)
                .collect::<Vec<_>>();

            let g = gauge_fit(cad, &pts)?;
            Ok(g * start)
        })
        .collect()
}

/// Summarizes how a bundle adjustment terminated.
fn outcome_summary(outcome: &MultiOutcome3) -> String {
    format!(
        "{:?} over {} solve(s), halt {:?}",
        outcome.quality(),
        outcome.solves().len(),
        outcome.halt()
    )
}

/// The largest distance any measurement point moved between two sets of transforms.
fn max_point_shift(bench: &Benchmark, a: &[Iso3], b: &[Iso3]) -> f64 {
    bench
        .scans
        .iter()
        .zip(a.iter().zip(b.iter()))
        .flat_map(|(scan, (ta, tb))| {
            scan.points
                .iter()
                .map(move |p| (ta * p.point - tb * p.point).norm())
        })
        .fold(0.0_f64, f64::max)
}

/// Scores the two candidates that require the solver against the two that do not, and holds them
/// to the behaviour established when this benchmark was built.
///
/// The warm start seeds the solver with the v0.2.16 answer and asks whether it holds or improves on
/// it. That is the sharper bug signal of the two, because a degradation cannot be blamed on a bad
/// starting pose and does not depend on any of the subtleties in how the metrics are built. The
/// cold start reproduces the chain that produced the v0.2.16 numbers, prealignment to CAD fit to
/// simultaneous alignment, so it is the like-for-like comparison, and it exercises
/// `points_to_surface3` as well.
///
/// # What is asserted, and what deliberately is not
///
/// The gates below are regression bars, not quality targets. They are set from measurements taken
/// on all four samples with generous margin, so that a change which alters the alignment
/// meaningfully trips one while ordinary variation does not.
///
/// Note what is *not* asserted. The solver is not required to beat v0.2.16 against the reference
/// mesh, because it does not and should not be expected to: it optimizes mutual agreement between
/// the scans, and those scans carry a per-scan deformation of the measurement volume, so pulling
/// them into tighter agreement necessarily bends the group away from an independent measurement of
/// the truth. Requiring both at once would be requiring the impossible. What is asserted is that
/// the reference fidelity does not collapse, which is what the original gate defect looked like.
///
/// Per-scan transforms are not compared either. `reference_priority` chooses the static mesh from
/// the data, and a different choice shifts every transform by a global rigid motion even when the
/// alignment is identical.
///
/// This test runs two bundle adjustments and four scorings per sample, about three minutes each and
/// twelve overall. It is not `#[ignore]`d because a regression gate that does not run is not a gate.
#[test]
fn solver_candidates_meet_their_gates() -> Result<()> {
    for case in find_cases()?.iter() {
        let loaded = LoadedCase::load(case)?;
        let design = design_key(case)?;

        let t = std::time::Instant::now();
        let (warm, warm_outcome) = run_multi_mesh(&loaded.scans, &loaded.cad_align)?;
        let warm_time = t.elapsed();

        let cad = case.cad()?;
        let pre = case
            .scans
            .iter()
            .map(|s| s.pre_align())
            .collect::<Result<Vec<_>>>()?;
        let t = std::time::Instant::now();
        let seeded = align_each_to_cad(&loaded, &cad, &pre)?;
        let seed_time = t.elapsed();

        let t = std::time::Instant::now();
        let (cold, cold_outcome) = run_multi_mesh(&loaded.scans, &seeded)?;
        let cold_time = t.elapsed();

        let mut scoring = std::time::Duration::ZERO;
        let graded = [
            ("floor", loaded.floor()),
            ("v0.2.16", loaded.cad_align.clone()),
            ("warm", warm.clone()),
            ("cold", cold.clone()),
        ]
        .into_iter()
        .map(|(name, transforms)| {
            let t = std::time::Instant::now();
            let score = loaded.score(&transforms)?;
            scoring += t.elapsed();
            Ok(Graded { name, score })
        })
        .collect::<Result<Vec<_>>>()?;

        print_table(&loaded.name, design, &graded, 0);

        println!("  warm: {}", outcome_summary(&warm_outcome));
        println!("  cold: {}", outcome_summary(&cold_outcome));
        println!(
            "  warm moved points up to {:.4} mm from the v0.2.16 seed",
            max_point_shift(&loaded.bench, &warm, &loaded.cad_align)
        );
        println!(
            "  cold moved points up to {:.4} mm from its CAD-aligned seed",
            max_point_shift(&loaded.bench, &cold, &seeded)
        );
        println!(
            "  timing: warm bundle {:.0}s, per-scan CAD fits {:.0}s, cold bundle {:.0}s, \
             scoring {:.0}s over 4 candidates",
            warm_time.as_secs_f64(),
            seed_time.as_secs_f64(),
            cold_time.as_secs_f64(),
            scoring.as_secs_f64(),
        );

        let floor = &graded[0].score;
        let v0216 = &graded[1].score;

        for (label, outcome, score, seed, moved, drift_bar) in [
            (
                "warm",
                &warm_outcome,
                &graded[2].score,
                &loaded.cad_align,
                &warm,
                WARM_DRIFT_BAR,
            ),
            (
                "cold",
                &cold_outcome,
                &graded[3].score,
                &seeded,
                &cold,
                COLD_DRIFT_BAR,
            ),
        ] {
            let where_ = format!("{} {}", loaded.name, label);

            // A bundle adjustment that merely exhausts its budget is a normal outcome and keeps its
            // alignments, but a halt means refinement gave up partway and the reason should be
            // understood rather than absorbed.
            assert!(
                outcome.halt().is_none(),
                "{where_}: refinement halted with {:?}",
                outcome.halt()
            );

            // The solver's own objective. It beat v0.2.16 on every sample by eight percent or more,
            // so a solver that no longer does has lost something real.
            assert!(
                score.consistency.median_abs <= v0216.consistency.median_abs * CONSISTENCY_SLACK,
                "{where_}: mutual consistency is {:.4} mm against v0.2.16 at {:.4} mm, so the                  simultaneous alignment is no longer earning its keep",
                score.consistency.median_abs,
                v0216.consistency.median_abs
            );

            // Not a demand that it beat v0.2.16 against the reference, only that it stays in the
            // same territory. The gate defect showed up here as a factor of two on the median and
            // an order of magnitude on the tail.
            assert!(
                score.reference.median_abs <= v0216.reference.median_abs * REFERENCE_CEILING,
                "{where_}: reference fidelity is {:.4} mm against v0.2.16 at {:.4} mm",
                score.reference.median_abs,
                v0216.reference.median_abs
            );
            assert!(
                score.reference.p95 <= v0216.reference.p95 * REFERENCE_CEILING,
                "{where_}: the reference tail reaches {:.4} mm at p95 against v0.2.16 at {:.4} mm",
                score.reference.p95,
                v0216.reference.p95
            );

            // Nothing constrained to move whole scans can beat a per-scan individual fit to the
            // reference. A candidate that does means the harness is measuring the two differently.
            assert!(
                score.reference.median_abs >= floor.reference.median_abs,
                "{where_}: scored {:.4} mm against the reference, below the individual-fit floor                  at {:.4} mm, which is not achievable",
                score.reference.median_abs,
                floor.reference.median_abs
            );

            // The most direct expression of the failure that this benchmark was built to find: a
            // solve that walks away from a good starting pose while reporting convergence.
            let drift = max_point_shift(&loaded.bench, moved, seed);
            assert!(
                drift < drift_bar,
                "{where_}: the solve moved points up to {drift:.4} mm from its seed, past the                  {drift_bar:.2} mm bar"
            );
        }

        // Two very different starting configurations must reach the same answer. This is the
        // property that distinguishes a solver which converges from one which merely stops.
        let split =
            (graded[2].score.consistency.median_abs - graded[3].score.consistency.median_abs).abs();
        assert!(
            split < WARM_COLD_AGREEMENT,
            "{}: the warm and cold starts disagree by {split:.4} mm on mutual consistency, so the              result depends on where the solve began",
            loaded.name
        );
    }

    Ok(())
}

/// Replicates the solver's own `sigma_max` estimator: a median absolute deviation scaled to a
/// standard deviation.
///
/// These scans carry no per-point uncertainty, so the solver's `inv_sigma` is one throughout and
/// its normalized residuals are the geometric ones an `Align3` reports. That makes this the
/// same number the solver would have computed for itself.
fn estimated_sigma_max(outcome: &MultiOutcome3) -> f64 {
    const MAD_TO_SIGMA: f64 = 1.4826;

    let residuals = outcome
        .alignments()
        .iter()
        .flat_map(|a| a.residuals().iter().copied())
        .collect::<Vec<_>>();

    let median = |v: &[f64]| {
        let mut s = v.to_vec();
        s.sort_by(f64::total_cmp);
        s[s.len() / 2]
    };

    let centre = median(&residuals);
    MAD_TO_SIGMA
        * median(
            &residuals
                .iter()
                .map(|r| (r - centre).abs())
                .collect::<Vec<_>>(),
        )
}

/// Diagnostic: does a warm start collapse its own robust weighting?
///
/// The warm start seeds the solver with an already converged answer, so its initial residuals are
/// small and tightly distributed. `sigma_max` is estimated from those residuals by a median
/// absolute deviation, so a tight distribution yields a small `sigma_max`, and MAGSAC then assigns
/// near-zero weight to everything outside a narrow band. If that is what is happening, the solve is
/// effectively running on a small and badly conditioned subset of its correspondences while still
/// reporting convergence.
///
/// The cold start seeds from per-scan CAD fits which are mutually inconsistent, so its residuals
/// are broader and its `sigma_max` larger. The hypothesis predicts three things: the warm seed's
/// estimated `sigma_max` is much smaller than the cold seed's, disabling refinement removes the
/// failure, and supplying a generous explicit `sigma_max` removes it too.
///
/// Note that the solver's own `Underdetermined` guard cannot catch this. It counts correspondences
/// whose weight is merely greater than zero, which stays large even when the weight distribution
/// has collapsed onto a narrow core.
#[test]
#[ignore]
fn diagnose_warm_start_sigma_collapse() -> Result<()> {
    let cases = find_cases()?;
    let case = cases
        .iter()
        .find(|c| c.name == "sample_00")
        .ok_or("sample_00 not found")?;
    let loaded = LoadedCase::load(case)?;

    let cad = case.cad()?;
    let pre = case
        .scans
        .iter()
        .map(|s| s.pre_align())
        .collect::<Result<Vec<_>>>()?;
    let cold_seed = align_each_to_cad(&loaded, &cad, &pre)?;

    let no_refine = MultiOptions3 {
        refinement_steps: 0,
        ..MultiOptions3::new(MULTI_MESH_MAX_DISTANCE)
    };

    // The two seeds, unrefined, purely to compare the noise scale each one would hand to MAGSAC.
    let (_, warm_raw) = run_multi_mesh_with(&loaded.scans, &loaded.cad_align, &no_refine)?;
    let (_, cold_raw) = run_multi_mesh_with(&loaded.scans, &cold_seed, &no_refine)?;
    println!(
        "\n  estimated sigma_max from the warm seed: {:.5} mm",
        estimated_sigma_max(&warm_raw)
    );
    println!(
        "  estimated sigma_max from the cold seed: {:.5} mm",
        estimated_sigma_max(&cold_raw)
    );

    let trials: Vec<(String, MultiOptions3)> = vec![
        (
            "default, auto sigma".to_string(),
            MultiOptions3::new(MULTI_MESH_MAX_DISTANCE),
        ),
        ("no refinement".to_string(), no_refine),
        (
            "explicit sigma_max 0.05".to_string(),
            MultiOptions3 {
                sigma_max: Some(0.05),
                ..MultiOptions3::new(MULTI_MESH_MAX_DISTANCE)
            },
        ),
        (
            "explicit sigma_max 0.15".to_string(),
            MultiOptions3 {
                sigma_max: Some(0.15),
                ..MultiOptions3::new(MULTI_MESH_MAX_DISTANCE)
            },
        ),
    ];

    println!(
        "\n  {:<24} | {:>8} {:>8} | {:>8} {:>8} {:>7} | {:>8}",
        "warm start variant", "cons.med", "out", "ref.med", "p90", "out", "shift"
    );
    for (name, opts) in trials {
        let (transforms, outcome) = run_multi_mesh_with(&loaded.scans, &loaded.cad_align, &opts)?;
        let score = loaded.score(&transforms)?;

        println!(
            "  {:<24} | {:>8.4} {:>8} | {:>8.4} {:>8.4} {:>7} | {:>8.4}   {}",
            name,
            score.consistency.median_abs,
            score.consistency.beyond_gate,
            score.reference.median_abs,
            score.reference.p90,
            score.reference.beyond_gate,
            max_point_shift(&loaded.bench, &transforms, &loaded.cad_align),
            outcome_summary(&outcome),
        );
    }

    Ok(())
}

/// Diagnostic: does a geometric gate stop the unweighted solve from drifting?
///
/// The previous experiment localized the failure to the initial unweighted solve, which moves a
/// known-good configuration most of a millimetre before robust weighting ever runs. Least squares
/// has no defense against a gross outlier, and `MultiOptions3` leaves both geometric gates off
/// by default, so with seventeen partially overlapping scans any correspondence sampled across a
/// region the other scan never saw pulls on the solve at full strength.
///
/// Refinement is disabled throughout so that only the initial solve is under test. If a gate fixes
/// it here, the defect is in the defaults rather than in the algorithm.
#[test]
#[ignore]
fn diagnose_unweighted_solve_gates() -> Result<()> {
    let cases = find_cases()?;
    let case = cases
        .iter()
        .find(|c| c.name == "sample_00")
        .ok_or("sample_00 not found")?;
    let loaded = LoadedCase::load(case)?;

    // A gate wider than the whole part stands in for the gate being off, which the options type
    // no longer allows to be expressed directly.
    const UNGATED: f64 = 1000.0;

    let gated = |d: f64, a: Option<f64>| MultiOptions3 {
        max_normal_angle: a,
        refinement_steps: 0,
        ..MultiOptions3::new(d)
    };

    let trials: Vec<(String, MultiOptions3)> = vec![
        ("no gates".to_string(), gated(UNGATED, None)),
        ("max_distance 1.0".to_string(), gated(1.0, None)),
        ("max_distance 0.5".to_string(), gated(0.5, None)),
        ("max_distance 0.2".to_string(), gated(0.2, None)),
        ("normal 30 deg only".to_string(), gated(UNGATED, Some(PI_6))),
        (
            "distance 0.5 + normal 30".to_string(),
            gated(0.5, Some(PI_6)),
        ),
    ];

    println!(
        "\n  unweighted solve, warm seed\n  {:<26} | {:>8} {:>8} | {:>8} {:>8} {:>7} | {:>8}",
        "gate", "cons.med", "out", "ref.med", "p90", "out", "shift"
    );
    for (name, opts) in trials {
        let (transforms, outcome) = run_multi_mesh_with(&loaded.scans, &loaded.cad_align, &opts)?;
        let score = loaded.score(&transforms)?;

        println!(
            "  {:<26} | {:>8.4} {:>8} | {:>8.4} {:>8.4} {:>7} | {:>8.4}   {}",
            name,
            score.consistency.median_abs,
            score.consistency.beyond_gate,
            score.reference.median_abs,
            score.reference.p90,
            score.reference.beyond_gate,
            max_point_shift(&loaded.bench, &transforms, &loaded.cad_align),
            outcome_summary(&outcome),
        );
    }

    Ok(())
}

// ================================================================================================
// Correspondence pruning
// ================================================================================================

/// Builds the full correspondence set for a multi-mesh run, one direction per unordered pair.
///
/// This mirrors what `multi_mesh_adjustment` does internally, reproduced here so the same
/// correspondences can be handed to the solver both whole and pruned. The static mesh and the
/// pairing direction are fixed rather than chosen by `reference_priority`, which keeps the pruned
/// and unpruned runs comparing like with like.
fn build_correspondences(
    scans: &[Mesh3],
    start: &[Iso3],
    sample: &GAPParams,
) -> Vec<MeshAlignPoint> {
    let mut pairs = Vec::new();
    for i in 0..scans.len() {
        for j in (i + 1)..scans.len() {
            pairs.push((j, i));
        }
    }

    pairs
        .par_iter()
        .flat_map(|&(mesh_i, ref_i)| {
            let t = start[ref_i].inv_mul(&start[mesh_i]);
            generate_alignment_points(&scans[mesh_i], &scans[ref_i], &t, sample)
                .into_iter()
                .map(|mp| MeshAlignPoint::new(mesh_i, mp, ref_i, 1.0))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// The jacobian row and world-frame correspondence for one alignment point at a given set of poses.
fn forward_row(
    point: &MeshAlignPoint,
    scans: &[Mesh3],
    start: &[Iso3],
    values: &AlignValues3,
) -> AlignStorage3 {
    let t_test = start[point.mesh_i];
    let t_ref = start[point.ref_i];

    let p_world = t_test * point.mp.point();
    let query = t_ref.inverse_transform_point(&p_world);
    let m = scans[point.ref_i].surface_closest_to(&query);
    let c_world = SurfacePoint3::new(t_ref * m.point(), t_ref.rotation * m.normal());

    point_surf_jacobian(&p_world, &c_world, values)
}

/// Prunes the correspondence set by greedy D-optimal selection, body by body.
///
/// Each correspondence belongs to exactly one test mesh, so grouping by that mesh partitions the
/// set and the selections cannot conflict. Within a body the choice is made on the *forward*
/// jacobian rows, the ones describing how that body's own six parameters respond.
///
/// That is an approximation worth naming. A correspondence also constrains the mesh it was matched
/// against, through the reverse jacobian block, and this ignores that when deciding what to drop.
/// The full treatment would be D-optimality over the whole `6 * (n - 1)` parameter space, which
/// `AlignInfo3` is not built for. Whether the approximation costs anything is exactly what
/// scoring the pruned result is meant to reveal.
fn prune_d_optimal(
    points: &[MeshAlignPoint],
    scans: &[Mesh3],
    start: &[Iso3],
    keep_per_body: usize,
) -> Result<Vec<MeshAlignPoint>> {
    let mut by_body = vec![Vec::new(); scans.len()];
    for (k, p) in points.iter().enumerate() {
        by_body[p.mesh_i].push(k);
    }

    let mut kept = Vec::new();
    for (i, idxs) in by_body.iter().enumerate() {
        if idxs.is_empty() {
            continue;
        }

        // The same parameterization the bundle uses for this body: rotation about the mesh's own
        // centre, with the seed pose carried in the working offset so the parameters start at zero.
        let c = scans[i].aabb().center();
        let local = Iso3::translation(c.x, c.y, c.z);
        let params = AlignParams3::new(
            AlignOrigin3::Local(local),
            Some(start[i] * local),
            Some(Dof6::all()),
        );
        let values = params.compute_values();

        let rows = idxs
            .par_iter()
            .map(|&k| forward_row(&points[k], scans, start, &values))
            .collect::<Vec<_>>();

        let weights = vec![1.0; rows.len()];
        let info = AlignInfo3::from_rows(rows, weights, Dof6::all())?;

        for s in info.select_d_optimal(keep_per_body, None) {
            kept.push(points[idxs[s]].clone());
        }
    }

    Ok(kept)
}

/// Diagnostic: does D-optimal pruning preserve the alignment while cutting the correspondence
/// count?
///
/// This is the payoff the stability tooling was rebuilt for. The original motivation was that
/// registering fifteen to twenty scans simultaneously is expensive, and that scoring points
/// independently and keeping the best double-counts redundancy: two points on the same flat patch
/// each look valuable, but the second adds almost nothing once the first is in. Greedy D-optimal
/// selection de-duplicates automatically because a candidate is measured against what has already
/// been chosen.
///
/// Cost matters here for a specific reason. Each Levenberg-Marquardt iteration factors the whole
/// jacobian with a dense pivoted QR inside the solver crate, which is serial and scales with the
/// row count, so it is the correspondence count rather than the thread count that governs runtime.
#[test]
#[ignore]
fn diagnose_pruned_alignment() -> Result<()> {
    let cases = find_cases()?;
    let case = cases
        .iter()
        .find(|c| c.name == "sample_00")
        .ok_or("sample_00 not found")?;
    let loaded = LoadedCase::load(case)?;

    let opts = MultiOptions3::new(MULTI_MESH_MAX_DISTANCE);
    let sample = GAPParams::defaults(MULTI_MESH_SPACING);
    let all = build_correspondences(&loaded.scans, &loaded.cad_align, &sample);
    println!("\n  {} correspondences before pruning", all.len());

    println!(
        "  {:<20} {:>9} {:>7} | {:>8} {:>6} | {:>8} {:>8} {:>6}",
        "kept per body", "total", "secs", "cons.med", "out", "ref.med", "p95", "out"
    );

    for keep in [usize::MAX, 4000, 2000, 1000, 500, 200] {
        let points = if keep == usize::MAX {
            all.clone()
        } else {
            prune_d_optimal(&all, &loaded.scans, &loaded.cad_align, keep)?
        };

        let t = std::time::Instant::now();
        let outcome = multi_mesh_adjustment_with_points(
            &loaded
                .scans
                .iter()
                .zip(loaded.cad_align.iter())
                .map(|(m, t)| AlignMesh::new(m, None, Some(t), None))
                .collect::<Vec<_>>(),
            points.clone(),
            0,
            &opts,
        )?;
        let secs = t.elapsed().as_secs_f64();

        let transforms = outcome
            .alignments()
            .iter()
            .map(|a| *a.full_transform())
            .collect::<Vec<_>>();
        let score = loaded.score(&transforms)?;

        println!(
            "  {:<20} {:>9} {:>7.0} | {:>8.4} {:>6} | {:>8.4} {:>8.4} {:>6}",
            if keep == usize::MAX {
                "all".to_string()
            } else {
                keep.to_string()
            },
            points.len(),
            secs,
            score.consistency.median_abs,
            score.consistency.beyond_gate,
            score.reference.median_abs,
            score.reference.p95,
            score.reference.beyond_gate,
        );
    }

    Ok(())
}
