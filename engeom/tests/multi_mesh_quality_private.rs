//! Quality benchmark for the multi-mesh alignment, run against private laser scan data.
//!
//! The data is a set of samples, each one a part scanned from 15 to 20 directions by a laser
//! profile sensor, together with a `reference.ply` produced by scanning the same physical part on
//! a Zeiss ATOS 5 structured light system. The ATOS mesh is lower resolution and captures less
//! detail, but its bulk shape is the best available source of truth, so it is what a candidate
//! alignment gets measured against.
//!
//! This file currently establishes that the data resolves and loads, that the stored transforms
//! place their scans where they claim to, and it builds and validates the frozen measurement set
//! that every candidate alignment will be scored against. The scorer itself follows.

#![cfg(all(feature = "private_tests", feature = "ply"))]

mod common;

use crate::common::PathPair;
use engeom::common::DistMode;
use engeom::common::kd_tree::KdTreeSearch;
use engeom::geom3::XyzQuat;
use engeom::geom3::align3::{AlignOptions3, AlignParams3, points_to_surface3};
use engeom::rayon::prelude::*;
use engeom::{Iso3, KdTree3, Mesh3, Point3, Result, UnitVec3};
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
const COVERAGE_MAX_NORMAL_ANGLE: f64 = std::f64::consts::PI / 6.0;

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
        // alignment for that offset. See `diagnose_group_fit_to_reference`, which removes it
        // first, and the scorer, which does so for every candidate.
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

/// Diagnostic: how much of the v0.2.16 deviation against the reference is a single global rigid
/// offset rather than registration error.
///
/// `reference.ply` was best fit to `cad.ply` in Zeiss Inspect, and the `cad-align.json` transforms
/// are a fit of the scan group to `cad.ply`. Both are aligned to the CAD, but that does not make
/// them aligned to each other: composing them carries the residual of *both* fits, and the leftover
/// is a rigid motion of the whole group which no multi-mesh alignment could remove or is
/// accountable for. Measuring it is the only way to know whether it matters.
#[test]
#[ignore]
fn diagnose_group_fit_to_reference() -> Result<()> {
    for case in find_cases()?.iter() {
        let bench = read_benchmark(&mut std::fs::File::open(
            case.dir.data().join(BENCHMARK_FILE),
        )?)?;
        let reference = case.reference()?;

        let mut raw = Vec::new();
        let mut floored = Vec::new();
        for (scan, files) in bench.scans.iter().zip(case.scans.iter()) {
            let start = files.cad_align()?;
            for p in scan.points.iter().filter(|p| p.on_reference) {
                raw.push(start * p.point);
                floored.push(scan.floor * p.point);
            }
        }

        let centroid =
            Point3::from(raw.iter().map(|p| p.coords).sum::<engeom::Vector3>() / raw.len() as f64);
        let params = AlignParams3::from_center(centroid, None);
        let outcome = points_to_surface3(&raw, &reference, params, &AlignOptions3::default())?;
        let correction = outcome.alignment().full_transform();

        let fitted = raw.iter().map(|p| correction * p).collect::<Vec<_>>();
        let shift = raw
            .iter()
            .zip(fitted.iter())
            .fold(0.0_f64, |acc, (a, b)| acc.max((a - b).norm()));

        println!(
            "{}: v0.2.16 unfitted {:.4} mm -> group-fitted {:.4} mm | floor {:.4} mm | \
             global correction moves points up to {:.4} mm | {:?}",
            case.name,
            median_abs_deviation(&reference, &raw),
            median_abs_deviation(&reference, &fitted),
            median_abs_deviation(&reference, &floored),
            shift,
            outcome.quality(),
        );
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
