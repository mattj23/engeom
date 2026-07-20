use engeom::common::DistMode;
use engeom::IntervalOps;
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

impl From<SelectOp> for engeom::SelectOp {
    fn from(val: SelectOp) -> Self {
        match val {
            SelectOp::Add => engeom::SelectOp::Add,
            SelectOp::Remove => engeom::SelectOp::Remove,
            SelectOp::Keep => engeom::SelectOp::KeepOnly,
        }
    }
}

#[pyclass(eq, eq_int, from_py_object, module = "engeom.common")]
#[derive(PartialEq, Copy, Clone, Debug, Default)]
pub enum DeviationMode {
    #[default]
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

#[pyclass(eq, eq_int, from_py_object, module = "engeom")]
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum VecDot {
    /// Use the raw dot product value as-is (can be negative).
    AsIs = 0,
    /// Use the absolute value of the dot product.
    Abs = 1,
    /// Clamp the dot product to zero from below (ignore anti-parallel normals).
    ClampPos = 2,
}

impl From<VecDot> for engeom::VecDot {
    fn from(val: VecDot) -> Self {
        match val {
            VecDot::AsIs => engeom::VecDot::AsIs,
            VecDot::Abs => engeom::VecDot::Abs,
            VecDot::ClampPos => engeom::VecDot::ClampPos,
        }
    }
}
