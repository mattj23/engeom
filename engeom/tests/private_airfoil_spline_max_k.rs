// #![cfg(feature = "private_tests")]

mod common;
use crate::common::PathPair;
use engeom::airfoil2::camber::extract_inscribed_circles;
use engeom::airfoil2::edges::fit_spline_max_k;
use engeom::airfoil2::{AfEdgeGeometry, OrientFwdAft, OrientUpperLower, SectionInput};
use engeom::geom2::CubicSpline2;
use engeom::io::write_tc_curve2_file;
use engeom::{Curve2, Point2, Result, Vector2};
use serde::{Deserialize, Serialize};
use std::path::Path;

const TEST_DATA_FOLDER: &str = "private-airfoil-spline-max-k";

#[test]
fn private_airfoil_spline_max_k() -> Result<()> {
    let test_dir = get_test_dir()?;
    let cases = get_cases(&test_dir.data())?;

    for case in cases {
        run_test_case(&case, &test_dir)?;
    }

    Ok(())
}

fn run_test_case(case: &TestCase, dir: &PathPair) -> Result<()> {
    // TODO: probably needs to have a subfolder and to clear it before starting

    for (i, section) in case.items.iter().enumerate() {
        // Observability output
        let output_root = format!("{}-sec{:03}", case.name, i);
        let curve = section.curve()?;
        write_tc_curve2_file(
            &dir.result().join(format!("{}.curve2", output_root)),
            &curve,
            1e-10,
        )?;

        let input = SectionInput::new(&curve, 1e-3);
        let unoriented = extract_inscribed_circles(&input)?;

        // Do a forced orientation
        let fwd_oriented = OrientFwdAft::Fwd(-Vector2::x()).apply(unoriented)?;
        let oriented = OrientUpperLower::Upper(Vector2::y()).apply(fwd_oriented)?;

        // Fit the spline max k
        let (fit_result, spline) = fit_spline_max_k(&input, oriented, true)?;
        let (center, radius) = match fit_result.edge.geometry {
            AfEdgeGeometry::SplineMaxK(c, r) => (c, r as f64),
            _ => panic!("Unexpected edge geometry"),
        };

        let output = Output {
            spline,
            center,
            radius,
            expected_r: section.le_r.unwrap(),
            expected_x: section.le_x.unwrap(),
        };

        let output_file = &dir.result().join(format!("{}.json", output_root));
        serde_json::to_writer_pretty(std::fs::File::create(output_file)?, &output)?;
        todo!()
    }

    Ok(())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Output {
    spline: CubicSpline2,
    center: Point2,
    radius: f64,
    expected_r: f64,
    expected_x: f64,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
struct SectionItem {
    xs: Vec<f64>,
    ys: Vec<f64>,
    le_r: Option<f64>,
    le_x: Option<f64>,
    le_y: Option<f64>,
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
