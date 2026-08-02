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
use engeom::geom3::XyzQuat;
use engeom::{Iso3, Mesh3, Point3, Result};
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
