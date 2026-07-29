#![cfg(feature = "private_tests")]

mod common;
use crate::common::PathPair;
use engeom::io::{
    DiffTanModel, Lptf3DsParams, Lptf3Load, load_lptf3_mesh_data, write_tc_mesh_file,
};
use engeom::{Iso3, Mesh3, Result};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

const TEST_DATA_FOLDER: &str = "multi-align-lptf3";

#[test]
fn multi_align_lptf3_private() -> Result<()> {
    let test_dir = get_test_dir()?;
    let cases = get_cases(&test_dir.data())?;

    for case in cases {
        let manifest_path = test_dir.data().join(case.manifest);
        let manifest: Manifest = serde_json::from_reader(std::fs::File::open(&manifest_path)?)?;
        let case_dir = test_dir.new_joined(
            manifest_path
                .parent()
                .unwrap()
                .file_name()
                .unwrap()
                .to_str()
                .unwrap(),
        )?;

        run_test_case(&manifest, &case_dir)?;
    }

    Ok(())
}

fn run_test_case(manifest: &Manifest, dir: &PathPair) -> Result<()> {
    // Load the reference mesh
    let mesh = Mesh3::load_stl(
        &dir.data().join(&manifest.reference_file),
        false,
        true,
        false,
    )?;

    // Save the reference mesh to the result directory
    write_tc_mesh_file(
        &dir.result().join("reference.tcmesh"),
        &mesh.to_data(),
        1e-5,
        false,
    )?;

    for item in manifest.items.iter() {
        let params = Lptf3DsParams::new(8, 1.0, 1.0, 0.010 * 25.4);
        let load = Lptf3Load::SmoothSample(params);
        let model = DiffTanModel::new(230.0, 80.0, 0.00037716);

        let file_path = dir.data().join(&item.file_name);
        let mesh_data = load_lptf3_mesh_data(&file_path, load, Some(&model))?;

        // The uncertainty model was supplied, so the standard deviations must have landed on the
        // mesh with one value per point.
        let stdev = mesh_data
            .point_stdev()
            .ok_or("Loaded mesh is missing its point standard deviations")?;
        assert_eq!(stdev.len(), mesh_data.point_count());

        let (points, faces, _) = mesh_data.into_parts();
        let mesh = Mesh3::new(points, faces, false);

        let output_path = dir.result().join(&item.file_name).with_extension("tcmesh");
        write_tc_mesh_file(&output_path, &mesh.to_data(), 1e-5, false)?;
        break;
    }

    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ManifestItem {
    file_name: String,
    iso: Iso3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Manifest {
    reference_file: String,
    items: Vec<ManifestItem>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TestCase {
    name: String,
    manifest: String,
}

fn get_cases(test_dir: &Path) -> Result<Vec<TestCase>> {
    let target = test_dir.join("cases.json");
    serde_json::from_reader(std::fs::File::open(target)?).map_err(|e| e.into())
}

fn get_test_dir() -> Result<PathPair> {
    let parent_dir = common::find_private_test_data()?;
    parent_dir.new_joined(TEST_DATA_FOLDER)
}
