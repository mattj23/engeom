use engeom::IntervalOps;
use engeom::common::DistMode;
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::{PyIndexError, PyValueError};
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
    #[getter]
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

/// A fixed-length mask of boolean values over a contiguous range of indices, used throughout the
/// library to select subsets of an entity's elements (a mesh's faces, a point cloud's points, and
/// so on).
///
/// A mask has a length fixed at construction, and every index from `0` to `len - 1` is either
/// `True` (selected) or `False` (not selected). It is the counterpart to a plain list of indices:
/// the mask form is what makes set algebra between two selections cheap, which is why the library
/// prefers it.
///
/// The set operations are exposed as Python operators rather than named methods: `~` inverts,
/// `|` unions, `&` intersects, and `-` takes the difference. The augmented forms (`|=`, `&=`, `-=`)
/// modify the mask in place instead of allocating a new one. All of the binary operations require
/// both masks to be the same length.
#[pyclass(from_py_object, module = "engeom.common")]
#[derive(Clone, Debug)]
pub struct IndexMask {
    inner: engeom::common::IndexMask,
}

impl IndexMask {
    pub fn get_inner(&self) -> &engeom::common::IndexMask {
        &self.inner
    }

    pub fn from_inner(inner: engeom::common::IndexMask) -> Self {
        Self { inner }
    }

    /// Turn a Python index (which may be negative, counting from the end) into a mask index,
    /// raising `IndexError` if it falls outside the mask rather than letting the core panic.
    fn resolve_index(&self, index: isize) -> PyResult<usize> {
        let len = self.inner.len() as isize;
        let resolved = if index < 0 { index + len } else { index };

        if resolved < 0 || resolved >= len {
            Err(PyIndexError::new_err(format!(
                "Index {index} is out of range for a mask of length {len}"
            )))
        } else {
            Ok(resolved as usize)
        }
    }
}

#[pymethods]
impl IndexMask {
    /// Create a mask of the given length, with every index set to the same initial value.
    ///
    /// :param length: the number of indices covered by the mask.
    /// :param value: the initial value for every index.
    #[new]
    #[pyo3(signature = (length, value = false))]
    fn new(length: usize, value: bool) -> Self {
        Self::from_inner(engeom::common::IndexMask::new(length, value))
    }

    /// Create a mask of the given length in which only the listed indices are `True`.
    ///
    /// :param indices: the indices to set to `True`.
    /// :param length: the total number of indices covered by the mask.
    /// :raises ValueError: if any index is out of bounds for the given length.
    #[staticmethod]
    fn from_indices(indices: Vec<usize>, length: usize) -> PyResult<Self> {
        let inner = engeom::common::IndexMask::try_from_indices(&indices, length)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Create a mask from a 1D numpy array of booleans, one per index.
    ///
    /// :param values: a 1D array with dtype `bool`. The mask takes its length from the array.
    #[staticmethod]
    fn from_bool_array(values: PyReadonlyArray1<'_, bool>) -> PyResult<Self> {
        let view = values.as_array();
        let mut inner = engeom::common::IndexMask::new(view.len(), false);
        for (i, value) in view.iter().enumerate() {
            if *value {
                inner.set(i, true);
            }
        }
        Ok(Self::from_inner(inner))
    }

    /// Return the indices which are `True`, in ascending order, as a 1D numpy array of `uint64`.
    fn to_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        Array1::from_iter(self.inner.iter_true().map(|i| i as u64)).into_pyarray(py)
    }

    /// Return the mask as a 1D numpy array of booleans with one entry per index.
    fn to_bool_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<bool>> {
        Array1::from_iter((0..self.inner.len()).map(|i| self.inner.get(i))).into_pyarray(py)
    }

    /// The number of indices which are `True`.
    #[getter]
    fn count_true(&self) -> usize {
        self.inner.count_true()
    }

    /// Return whether at least one index is `True`.
    ///
    /// Note that this is about the mask's contents, not its length: `len()` is the size of the
    /// range the mask covers, and a mask of a thousand indices which are all `False` still has a
    /// length of a thousand.
    fn any(&self) -> bool {
        !self.inner.is_empty()
    }

    /// Return whether every index is `True`. A zero length mask returns `True`.
    fn all(&self) -> bool {
        self.inner.count_true() == self.inner.len()
    }

    /// Return an independent copy of this mask. Because a mask is mutable, assigning one to
    /// another name aliases it rather than duplicating it, and this is the way to get a version
    /// which can be modified without disturbing the original.
    ///
    /// Named to match `Mesh3.cloned` and `MeshData3.cloned`.
    fn cloned(&self) -> Self {
        self.clone()
    }

    /// Set every index in the mask to the same value, in place.
    ///
    /// :param value: the value to set every index to.
    fn fill(&mut self, value: bool) {
        self.inner.fill(value);
    }

    /// Set the lowest `True` index to `False` and return it, or return `None` if no index is
    /// `True`. This lets a mask be consumed as a work queue.
    fn pop_true(&mut self) -> Option<usize> {
        self.inner.pop_true()
    }

    /// Return a new mask which is the inverse of this one, leaving this one unchanged.
    fn __invert__(&self) -> Self {
        Self::from_inner(self.inner.not())
    }

    /// Return a new mask holding the indices which are `True` in either mask (set union).
    fn __or__(&self, other: &Self) -> PyResult<Self> {
        let inner = self
            .inner
            .or(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Return a new mask holding the indices which are `True` in both masks (set intersection).
    fn __and__(&self, other: &Self) -> PyResult<Self> {
        let inner = self
            .inner
            .and(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Return a new mask holding the indices which are `True` in this mask but not the other
    /// (set difference).
    fn __sub__(&self, other: &Self) -> PyResult<Self> {
        let inner = self
            .inner
            .and_not(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn __ior__(&mut self, other: &Self) -> PyResult<()> {
        self.inner
            .or_mut(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __iand__(&mut self, other: &Self) -> PyResult<()> {
        self.inner
            .and_mut(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __isub__(&mut self, other: &Self) -> PyResult<()> {
        self.inner
            .and_not_mut(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<bool> {
        Ok(self.inner.get(self.resolve_index(index)?))
    }

    fn __setitem__(&mut self, index: isize, value: bool) -> PyResult<()> {
        let index = self.resolve_index(index)?;
        self.inner.set(index, value);
        Ok(())
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!(
            "IndexMask(length={}, true={})",
            self.inner.len(),
            self.inner.count_true()
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
