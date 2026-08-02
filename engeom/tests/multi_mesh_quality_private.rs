//! Quality benchmark for the multi-mesh alignment, run against private laser scan data.
//!
//! The data is a set of samples, each one a part scanned from 15 to 20 directions by a laser
//! profile sensor, together with a `reference.ply` produced by scanning the same physical part on
//! a Zeiss ATOS 5 structured light system. The ATOS mesh is lower resolution and captures less
//! detail, but its bulk shape is the best available source of truth, so it is what a candidate
//! alignment gets measured against.
//!
//! This file currently only establishes that the data resolves and loads. The measurement
//! machinery follows.

#![cfg(all(feature = "private_tests", feature = "ply"))]

mod common;

use crate::common::PathPair;
use engeom::{Mesh3, Result};
use std::path::PathBuf;

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
            let mesh = Mesh3::load_ply(&scan.mesh, IS_SOLID)?;
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
