use engeom::IntervalOps;
use engeom::common::DistMode;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Parse an angle-direction string token (`"cw"` or `"ccw"`) into the core `AngleDir` enum.
pub fn angle_dir_from_str(s: &str) -> PyResult<engeom::common::AngleDir> {
    match s {
        "cw" => Ok(engeom::common::AngleDir::Cw),
        "ccw" => Ok(engeom::common::AngleDir::Ccw),
        _ => Err(PyValueError::new_err(format!(
            "Invalid angle direction '{s}', expected 'cw' or 'ccw'"
        ))),
    }
}

/// Render a core `AngleDir` enum as its string token (`"cw"` or `"ccw"`).
pub fn angle_dir_to_str(val: engeom::common::AngleDir) -> &'static str {
    match val {
        engeom::common::AngleDir::Cw => "cw",
        engeom::common::AngleDir::Ccw => "ccw",
    }
}

/// Returns the positive angle (in radians) needed to rotate `radians0` to `radians1` in the
/// given direction. The result is always in the range [0, 2π].
#[pyfunction]
pub fn angle_in_direction(radians0: f64, radians1: f64, angle_dir: &str) -> PyResult<f64> {
    Ok(engeom::common::angle_in_direction(
        radians0,
        radians1,
        angle_dir_from_str(angle_dir)?,
    ))
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

/// Parse a face-selection operation token (`"add"`, `"remove"`, or `"keep"`) into the core
/// `SelectOp` enum.
pub fn select_op_from_str(s: &str) -> PyResult<engeom::SelectOp> {
    match s {
        "add" => Ok(engeom::SelectOp::Add),
        "remove" => Ok(engeom::SelectOp::Remove),
        "keep" => Ok(engeom::SelectOp::KeepOnly),
        _ => Err(PyValueError::new_err(format!(
            "Invalid select operation '{s}', expected 'add', 'remove', or 'keep'"
        ))),
    }
}

/// Parse a deviation-mode token (`"point"` or `"plane"`) into the core `DistMode` enum.
pub fn deviation_mode_from_str(s: &str) -> PyResult<DistMode> {
    match s {
        "point" => Ok(DistMode::ToPoint),
        "plane" => Ok(DistMode::ToPlane),
        _ => Err(PyValueError::new_err(format!(
            "Invalid deviation mode '{s}', expected 'point' or 'plane'"
        ))),
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

impl From<Resample> for engeom::Resample {
    fn from(val: Resample) -> Self {
        match val {
            Resample::Count(count) => engeom::Resample::ByCount(count),
            Resample::Spacing(spacing) => engeom::Resample::BySpacing(spacing),
            Resample::MaxSpacing(max_spacing) => engeom::Resample::ByMaxSpacing(max_spacing),
        }
    }
}

/// A continuous range of angles, specified by a starting angle and a positive (counter-clockwise)
/// included length.
///
/// This is a read-only view returned by other geometric queries (such as
/// `Circle2.intersection_interval` and `Arc2.angle_interval`); it is not directly constructible
/// from Python. Only the core query surface is exposed here; the full interval algebra (overlap,
/// intersection, union, etc.) is not currently bound.
#[pyclass(from_py_object, module = "engeom.common")]
#[derive(Copy, Clone, Debug)]
pub struct AngleInterval {
    inner: engeom::AngleInterval,
}

impl AngleInterval {
    pub fn get_inner(&self) -> &engeom::AngleInterval {
        &self.inner
    }

    pub fn from_inner(inner: engeom::AngleInterval) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AngleInterval {
    /// The starting angle of the interval, in radians, in the range [0, 2π].
    #[getter]
    fn min(&self) -> f64 {
        self.inner.min()
    }

    /// The ending angle of the interval, in radians. If less than `min`, the interval wraps
    /// beyond 2π.
    #[getter]
    fn max(&self) -> f64 {
        self.inner.max()
    }

    /// The angular extent (included angle) of the interval, in radians.
    #[getter]
    fn extent(&self) -> f64 {
        self.inner.extent()
    }

    /// The angle at the midpoint of the interval, in radians.
    #[getter]
    fn center(&self) -> f64 {
        self.inner.center()
    }

    /// Whether the interval contains zero angular extent.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Returns whether the given angle (in radians) falls within the interval.
    fn contains_value(&self, x: f64) -> bool {
        self.inner.contains_value(x)
    }

    fn __repr__(&self) -> String {
        format!(
            "AngleInterval(min={}, max={})",
            self.inner.min(),
            self.inner.max()
        )
    }
}

/// Parse a vector dot-product handling token (`"as_is"`, `"abs"`, or `"clamp_pos"`) into the core
/// `VecDot` enum.
pub fn vec_dot_from_str(s: &str) -> PyResult<engeom::VecDot> {
    match s {
        "as_is" => Ok(engeom::VecDot::AsIs),
        "abs" => Ok(engeom::VecDot::Abs),
        "clamp_pos" => Ok(engeom::VecDot::ClampPos),
        _ => Err(PyValueError::new_err(format!(
            "Invalid vector dot mode '{s}', expected 'as_is', 'abs', or 'clamp_pos'"
        ))),
    }
}
