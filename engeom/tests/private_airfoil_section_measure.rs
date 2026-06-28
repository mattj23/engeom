#![cfg(feature = "private_tests")]

mod common;
use crate::common::PathPair;
use approx::assert_relative_eq;
use engeom::airfoil2::{AfEdgeSearch, AfGeometry, OrientFwdAft, OrientUpperLower};
use engeom::common::dist;
use engeom::io::write_tc_curve2_file;
use engeom::{Curve2, Point2, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;

const TEST_DATA_FOLDER: &str = "private-airfoil-section-measure";

#[test]
fn private_airfoil_section_measure() -> Result<()> {
    let test_dir = get_test_dir()?;
    let cases = get_cases(&test_dir.data())?;

    for case in cases {
        run_test_case(&case, &test_dir)?;
    }

    Ok(())
}

fn run_test_case(case: &TestCase, dir: &PathPair) -> Result<()> {
    // TODO: This needs some thought about how to unify the measurements and tolerances

    // Delete the subfolder if it exists
    let sub_folder = dir.result().join(&case.name);
    if sub_folder.exists() {
        std::fs::remove_dir_all(&sub_folder)?;
    }
    std::fs::create_dir_all(&sub_folder)?;

    for (i, section) in case.items.iter().enumerate() {
        // Observability output
        let output_root = format!("sec-{:03}", i);
        let curve = section.curve()?;
        write_tc_curve2_file(
            &sub_folder.join(format!("{}.tccurve2", output_root)),
            &curve,
            1e-10,
        )?;

        // Airfoil nominal analysis
        let nominal = AfGeometry::try_from_geometric_analysis(
            &curve,
            1e-3,
            OrientFwdAft::TmaxFwd,
            OrientUpperLower::Curvature,
            AfEdgeSearch::Auto,
            AfEdgeSearch::Auto,
        )?;

        // Chord dimension
        let chord_distance = dist(&nominal.camber.at_front(), &nominal.camber.at_back());
        assert_relative_eq!(chord_distance, section.chord, epsilon = 5e-2);
    }

    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Output {}

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
        let points = self
            .xs
            .iter()
            .zip(self.ys.iter())
            .map(|(x, y)| Point2::new(*x, *y))
            .collect::<Vec<_>>();
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
