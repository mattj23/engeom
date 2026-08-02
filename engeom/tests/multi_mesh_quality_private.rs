//! Quality benchmark for the multi-mesh alignment, run against private laser scan data.
//!
//! The data is a set of samples, each one a part scanned from 15 to 20 directions by a laser
//! profile sensor, together with a `reference.ply` produced by scanning the same physical part on
//! a Zeiss ATOS 5 structured light system. The ATOS mesh is lower resolution and captures less
//! detail, but its bulk shape is the best available source of truth, so it is what a candidate
//! alignment gets measured against.
//!
//! This file currently establishes that the data resolves, that every mesh in it loads, and that
//! the stored transforms place their scans where they claim to. The measurement machinery follows.

#![cfg(all(feature = "private_tests", feature = "ply"))]

mod common;

use crate::common::PathPair;
use engeom::common::DistMode;
use engeom::common::kd_tree::KdTreeSearch;
use engeom::geom3::XyzQuat;
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

/// The name of the frozen artifact within each sample folder.
const BENCHMARK_FILE: &str = "benchmark.bin";

const BENCHMARK_MAGIC: &[u8; 4] = b"EGMB";

/// One frozen measurement point, in the coordinate frame of the scan it was sampled from.
#[derive(Clone, Debug)]
struct BenchPoint {
    point: Point3,
    normal: UnitVec3,

    /// The scans this point should be measured against, frozen from the reference configuration.
    ///
    /// Empty is legitimate: a point in a region only one scan saw contributes nothing to mutual
    /// consistency, but still counts toward reference fidelity.
    partners: Vec<u32>,
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
// - 4 x `f64`: spacing, max distance, max normal angle, edge margin
// - 1 x `u32`: scan count
// - per scan:
//   - 1 x `u32`: byte length of the id, then that many bytes of UTF-8
//   - 2 x `u32`: source mesh vertex count, face count
//   - 1 x `u32`: measurement point count
//   - per point: 3 x `f32` position, 3 x `f32` normal, 1 x `u32` partner count, that many `u32`
//
// Positions are stored as `f32`. At the ~100 mm scale of these parts that is a quantization of
// well under a hundredth of a micron, which is four orders of magnitude below the deviations being
// measured, and it is identical for every candidate scored against the artifact.

fn write_benchmark<W: Write>(writer: &mut W, bench: &Benchmark) -> Result<()> {
    writer.write_all(BENCHMARK_MAGIC)?;
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

        writer.write_all(&(scan.points.len() as u32).to_le_bytes())?;
        for p in scan.points.iter() {
            for v in [p.point.x, p.point.y, p.point.z] {
                writer.write_all(&(v as f32).to_le_bytes())?;
            }
            for v in [p.normal.x, p.normal.y, p.normal.z] {
                writer.write_all(&(v as f32).to_le_bytes())?;
            }
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

            let partner_count = r.read_u32()? as usize;
            let mut partners = Vec::with_capacity(partner_count);
            for _ in 0..partner_count {
                partners.push(r.read_u32()?);
            }

            points.push(BenchPoint {
                point,
                normal,
                partners,
            });
        }

        scans.push(ScanBench {
            id,
            vertex_count,
            face_count,
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

    // A kd-tree over each scan's boundary vertices, used to keep matches away from the edge of the
    // scan they matched against. Testing proximity to boundary vertices rather than to the
    // boundary edges themselves is an approximation, but these meshes trianglate far finer than
    // the margin, so a point near an edge is always near one of its endpoints.
    let boundaries = meshes
        .par_iter()
        .map(|mesh| {
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
        })
        .collect::<Vec<_>>();

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
                        let q = rel * mp.point();
                        let n = rel.rotation * mp.normal().into_inner();

                        let match_ = meshes[j].surface_closest_to(&q);
                        let dist = (q - match_.point()).norm();
                        if dist > COVERAGE_MAX_DISTANCE {
                            return None;
                        }
                        if n.angle(&match_.normal()) > COVERAGE_MAX_NORMAL_ANGLE {
                            return None;
                        }

                        let tree = boundaries[j].as_ref()?;
                        if tree.nearest_one(&match_.point()).1 < COVERAGE_EDGE_MARGIN {
                            return None;
                        }

                        Some(j as u32)
                    })
                    .collect::<Vec<_>>();

                BenchPoint {
                    point: mp.point(),
                    normal: mp.normal(),
                    partners,
                }
            })
            .collect::<Vec<_>>();

        scans.push(ScanBench {
            id: case.scans[i].id.clone(),
            vertex_count: mesh.points().len() as u32,
            face_count: mesh.faces().len() as u32,
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

        println!(
            "{}: {} scans, {} points ({:.1}% covered), {} measured pairs, {:.1} MB",
            case.name,
            bench.scans.len(),
            bench.point_count(),
            100.0 * covered as f64 / bench.point_count() as f64,
            bench.pair_count(),
            bytes as f64 / 1.0e6,
        );
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
                points: vec![
                    BenchPoint {
                        point: Point3::new(1.0, 2.0, 3.0),
                        normal: UnitVec3::new_normalize(engeom::Vector3::new(1.0, 1.0, 0.0)),
                        partners: vec![1, 2],
                    },
                    BenchPoint {
                        point: Point3::new(-4.5, 0.25, 9.75),
                        normal: UnitVec3::new_normalize(engeom::Vector3::new(0.0, 0.0, -1.0)),
                        partners: vec![],
                    },
                ],
            },
            ScanBench {
                id: "10002".to_string(),
                vertex_count: 9,
                face_count: 8,
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

        for (p, q) in a.points.iter().zip(b.points.iter()) {
            // f32 storage, so the tolerance is the quantization rather than a solver epsilon.
            assert!((p.point - q.point).norm() < 1e-5);
            assert!((p.normal.into_inner() - q.normal.into_inner()).norm() < 1e-6);
            assert_eq!(p.partners, q.partners);
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
            points: vec![BenchPoint {
                point: Point3::new(1.0, 2.0, 3.0),
                normal: UnitVec3::new_normalize(engeom::Vector3::new(0.0, 0.0, 1.0)),
                partners: vec![7],
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
