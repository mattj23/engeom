// #![cfg(feature = "private_tests")]

mod common;
use crate::common::PathPair;
use engeom::{Curve2, Point2, Result};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use engeom::io::write_tc_curve2_file;

const TEST_DATA_FOLDER: &str = "airfoil-section-measure";

#[test]
fn airfoil_section_measure_private() -> Result<()> {
    let test_dir = get_test_dir()?;
    let cases = get_cases(&test_dir.data())?;

    for case in cases {
        run_test_case(&case, &test_dir)?;
    }

    Ok(())
}

fn run_test_case(case: &TestCase, dir: &PathPair) -> Result<()> {
    for (i, section) in case.items.iter().enumerate() {
        let output_root = format!("{}-sec{:03}", case.name, i);
        let curve = section.curve()?;
        write_tc_curve2_file(&dir.result().join(format!("{}.curve2", output_root)), &curve, 1e-6)?;
    }

    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ThicknessGage {
    x: f64,
    t: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SectionItem {
    xs: Vec<f64>,
    ys: Vec<f64>,
    chord: f64,
    chord_angle: f64,
    t_max: f64,
    t_gages: Vec<ThicknessGage>,
    le_r: Option<f64>,
    te_r: Option<f64>,
    le_x: Option<f64>,
    le_y: Option<f64>,
    te_x: Option<f64>,
    te_y: Option<f64>,
}

impl SectionItem {
    fn curve(&self) -> Result<Curve2> {
        let points = self.xs.iter().zip(self.ys.iter()).map(|(x, y)| Point2::new(*x, *y)).collect::<Vec<_>>();
        Curve2::from_points(&points, 1e-6, false)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TestCase {
    name: String,
    items: Vec<SectionItem>,
}

fn get_cases(test_dir: &Path) -> Result<Vec<TestCase>> {
    let mut test_cases = Vec::new();
    for entry in std::fs::read_dir(test_dir)? {
        let entry = entry?;
        let path = entry.path();
        // Check if the file ends with ".json":
        if path.extension().and_then(|s| s.to_str()) == Some("json") {
            let name = path.file_stem().unwrap().to_str().unwrap().to_string();
            let items: Vec<SectionItem> = serde_json::from_reader(std::fs::File::open(path)?)?;
            test_cases.push(TestCase { name, items })
        }
    }
    Ok(test_cases)
}

fn get_test_dir() -> Result<PathPair> {
    let parent_dir = common::find_private_test_data()?;
    parent_dir.new_joined(TEST_DATA_FOLDER)
}
