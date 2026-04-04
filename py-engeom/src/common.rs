use engeom::common::DistMode;
use pyo3::prelude::*;

#[pyclass(eq, eq_int, from_py_object, module = "engeom.common")]
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum AngleDir {
    /// Clockwise rotation (negative angular direction).
    Cw = 0,
    /// Counter-clockwise rotation (positive angular direction).
    Ccw = 1,
}

impl From<AngleDir> for engeom::common::AngleDir {
    fn from(val: AngleDir) -> Self {
        match val {
            AngleDir::Cw => engeom::common::AngleDir::Cw,
            AngleDir::Ccw => engeom::common::AngleDir::Ccw,
        }
    }
}

impl From<engeom::common::AngleDir> for AngleDir {
    fn from(val: engeom::common::AngleDir) -> Self {
        match val {
            engeom::common::AngleDir::Cw => AngleDir::Cw,
            engeom::common::AngleDir::Ccw => AngleDir::Ccw,
        }
    }
}

/// Returns the positive angle (in radians) needed to rotate `radians0` to `radians1` in the
/// given direction. The result is always in the range [0, 2π].
#[pyfunction]
pub fn angle_in_direction(radians0: f64, radians1: f64, angle_dir: AngleDir) -> f64 {
    engeom::common::angle_in_direction(radians0, radians1, angle_dir.into())
}

/// Returns the signed shortest angular distance from `radians0` to `radians1`. A positive result
/// means the shortest path is counter-clockwise; negative means clockwise. Magnitude is in [0, π].
#[pyfunction]
pub fn shortest_angle_between(radians0: f64, radians1: f64) -> f64 {
    engeom::common::shortest_angle_between(radians0, radians1)
}

/// Re-expresses an angle in the range (-π, π]. Equivalent angles are preserved; -π maps to π.
#[pyfunction]
pub fn angle_signed_pi(radians: f64) -> f64 {
    engeom::common::angle_signed_pi(radians)
}

/// Re-expresses an angle in the range [0, 2π].
#[pyfunction]
pub fn angle_to_2pi(radians: f64) -> f64 {
    engeom::common::angle_to_2pi(radians)
}

/// Returns the signed complement of an angle with respect to a full rotation (2π). A positive
/// input returns a negative complement, and vice versa. Result is in (-2π, 2π].
#[pyfunction]
pub fn signed_compliment_2pi(radians: f64) -> f64 {
    engeom::common::signed_compliment_2pi(radians)
}

#[pyclass(eq, eq_int, from_py_object, module = "engeom.common")]
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum SelectOp {
    Add = 0,
    Remove = 1,
    Keep = 2,
}

impl From<SelectOp> for engeom::common::SelectOp {
    fn from(val: SelectOp) -> Self {
        match val {
            SelectOp::Add => engeom::common::SelectOp::Add,
            SelectOp::Remove => engeom::common::SelectOp::Remove,
            SelectOp::Keep => engeom::common::SelectOp::KeepOnly,
        }
    }
}

#[pyclass(eq, eq_int, from_py_object, module = "engeom.common")]
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum DeviationMode {
    Point,
    Plane,
}

impl From<DeviationMode> for DistMode {
    fn from(val: DeviationMode) -> Self {
        match val {
            DeviationMode::Point => DistMode::ToPoint,
            DeviationMode::Plane => DistMode::ToPlane,
        }
    }
}

#[pyclass(from_py_object)]
#[derive(Copy, Clone, Debug)]
pub enum Resample {
    Count(usize),
    Spacing(f64),
    MaxSpacing(f64),
}

#[pymethods]
impl Resample {
    fn __repr__(&self) -> String {
        match self {
            Resample::Count(count) => format!("Resample.Count({})", count),
            Resample::Spacing(spacing) => format!("Resample.Spacing({})", spacing),
            Resample::MaxSpacing(max_spacing) => {
                format!("Resample.MaxSpacing({})", max_spacing)
            }
        }
    }
}

impl From<Resample> for engeom::common::Resample {
    fn from(val: Resample) -> Self {
        match val {
            Resample::Count(count) => engeom::common::Resample::ByCount(count),
            Resample::Spacing(spacing) => engeom::common::Resample::BySpacing(spacing),
            Resample::MaxSpacing(max_spacing) => {
                engeom::common::Resample::ByMaxSpacing(max_spacing)
            }
        }
    }
}
