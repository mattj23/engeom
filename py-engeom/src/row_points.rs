//! Bindings for the row-organized point format (`.tcrpf3`).
//!
//! This conversion layer wraps `engeom::io::tol_compress::row_points`. Rows cross the boundary as a
//! list of `(n, 3)` float arrays with a separate `uint32` ordinal array, avoiding one Python object
//! per point.

use crate::conversions::{array_to_points3, points_to_array};
use engeom::io::{PointRow3, RowAxis, RowPointsScan3 as Inner};
use engeom::tol_compress::{Metadata, Value};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyBytes, PyDict, PyFloat, PyInt, PyString};
use std::collections::HashMap;
use std::path::PathBuf;

/// Parse the `along` keyword that identifies the world axis varying along a row.
fn axis_from_str(s: &str) -> PyResult<RowAxis> {
    match s {
        "x" => Ok(RowAxis::X),
        "y" => Ok(RowAxis::Y),
        _ => Err(PyValueError::new_err(format!(
            "Invalid along '{s}', expected 'x' or 'y'"
        ))),
    }
}

fn axis_name(axis: RowAxis) -> &'static str {
    match axis {
        RowAxis::X => "x",
        RowAxis::Y => "y",
    }
}

/// Convert a Python value into a metadata value.
///
/// Check `bool` before `int` because Python defines `bool` as an `int` subclass. Otherwise, a flag
/// stored as integer 1 would return as `1` instead of `True`.
fn value_from_py(value: &Bound<'_, PyAny>) -> PyResult<Value> {
    if let Ok(v) = value.cast::<PyBool>() {
        Ok(Value::Bool(v.is_true()))
    } else if let Ok(v) = value.cast::<PyInt>() {
        Ok(Value::I64(v.extract()?))
    } else if let Ok(v) = value.cast::<PyFloat>() {
        Ok(Value::F64(v.extract()?))
    } else if let Ok(v) = value.cast::<PyString>() {
        Ok(Value::Text(v.extract()?))
    } else if let Ok(v) = value.cast::<PyBytes>() {
        Ok(Value::Bytes(v.as_bytes().to_vec()))
    } else {
        Err(PyValueError::new_err(
            "Metadata values must be bool, int, float, str, or bytes",
        ))
    }
}

fn value_to_py<'py>(py: Python<'py>, value: &Value) -> PyResult<Bound<'py, PyAny>> {
    match value {
        Value::Bool(v) => Ok(v.into_pyobject(py)?.to_owned().into_any()),
        Value::I64(v) => Ok(v.into_pyobject(py)?.into_any()),
        Value::F64(v) => Ok(v.into_pyobject(py)?.into_any()),
        Value::Text(v) => Ok(v.into_pyobject(py)?.into_any()),
        Value::Bytes(v) => Ok(PyBytes::new(py, v).into_any()),
    }
}

/// A scan whose points are grouped into rows produced by a rasterizing sensor.
///
/// This is the in-memory form of a `.tcrpf3` file. Use it to read and write the row structure. To
/// process geometry, load the file as `PointCloud3` or `Mesh3`, which applies thinning and, for a
/// mesh, row-strip triangulation.
#[pyclass(module = "engeom.geom3")]
pub struct RowPointsScan3 {
    inner: Inner,
}

impl RowPointsScan3 {
    pub fn from_inner(inner: Inner) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &Inner {
        &self.inner
    }
}

#[pymethods]
impl RowPointsScan3 {
    /// Build a scan from its rows.
    ///
    /// Each row is an `(n, 3)` point array with an ordinal that records its position in the sensor
    /// sweep. Ordinals must strictly increase. A skipped row leaves an ordinal gap instead of
    /// causing later rows to be renumbered.
    #[staticmethod]
    #[pyo3(signature = (rows, ordinals, row_pitch, *, columns=None, col_pitch=None, along="x", name=None, metadata=None))]
    #[allow(clippy::too_many_arguments)]
    fn from_rows(
        rows: Vec<PyReadonlyArray2<f64>>,
        ordinals: PyReadonlyArray1<u32>,
        row_pitch: f64,
        columns: Option<Vec<PyReadonlyArray1<u32>>>,
        col_pitch: Option<f64>,
        along: &str,
        name: Option<String>,
        metadata: Option<HashMap<String, Py<PyAny>>>,
    ) -> PyResult<Self> {
        let ordinals = ordinals.as_slice()?;
        if ordinals.len() != rows.len() {
            return Err(PyValueError::new_err(format!(
                "Got {} rows but {} ordinals, which must match",
                rows.len(),
                ordinals.len()
            )));
        }

        if let Some(c) = &columns
            && c.len() != rows.len()
        {
            return Err(PyValueError::new_err(format!(
                "Got {} rows but {} column arrays, which must match",
                rows.len(),
                c.len()
            )));
        }

        let mut out = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            let points = array_to_points3(&row.as_array())?;
            let mut entry = PointRow3::new(ordinals[i], points);

            if let Some(c) = &columns {
                let cols = c[i].as_slice()?.to_vec();
                if cols.len() != entry.points.len() {
                    return Err(PyValueError::new_err(format!(
                        "Row {i} has {} points but {} column indices",
                        entry.points.len(),
                        cols.len()
                    )));
                }
                entry = entry.with_columns(cols);
            }

            out.push(entry);
        }

        let mut inner = Inner::new(out, row_pitch).with_along_axis(axis_from_str(along)?);
        inner.col_pitch = col_pitch;
        inner.name = name;

        if let Some(map) = metadata {
            let converted = Python::attach(|py| {
                let mut m = Metadata::new();
                for (k, v) in map.iter() {
                    m.insert(k.clone(), value_from_py(v.bind(py))?);
                }
                Ok::<_, PyErr>(m)
            })?;
            inner.metadata = converted;
        }

        Ok(Self { inner })
    }

    /// Read a scan from a `.tcrpf3` file, without thinning or meshing it.
    #[staticmethod]
    fn read(path: PathBuf) -> PyResult<Self> {
        let inner = engeom::io::read_tc_row_points_file(&path)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Write the scan to a `.tcrpf3` file.
    ///
    /// `tol` is the maximum round-trip position error for any point, in the coordinate units. The
    /// format uses the narrowest storage width that guarantees this tolerance, so a looser `tol`
    /// produces a smaller file. When studying source-data quantization, set `tol` well below that
    /// quantization to make storage error negligible.
    fn write(&self, path: PathBuf, tol: f64) -> PyResult<()> {
        engeom::io::write_tc_row_points_file(&path, &self.inner, tol)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    /// The number of rows in the scan.
    #[getter]
    fn row_count(&self) -> usize {
        self.inner.rows.len()
    }

    /// The total number of points across every row.
    #[getter]
    fn point_count(&self) -> usize {
        self.inner.point_count()
    }

    /// The nominal distance between consecutive row ordinals.
    #[getter]
    fn row_pitch(&self) -> f64 {
        self.inner.row_pitch
    }

    /// The nominal spacing between points along a row, if the scan records one.
    #[getter]
    fn col_pitch(&self) -> Option<f64> {
        self.inner.col_pitch
    }

    /// Which world axis varies along a row, `"x"` or `"y"`.
    #[getter]
    fn along(&self) -> &'static str {
        axis_name(self.inner.along)
    }

    /// The scan's name, if it has one.
    #[getter]
    fn name(&self) -> Option<&str> {
        self.inner.name.as_deref()
    }

    /// Each row's position in the sensor's sweep, as a `(rows,)` array of dtype uint32.
    #[getter]
    fn ordinals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.inner
            .rows
            .iter()
            .map(|r| r.ordinal)
            .collect::<Vec<_>>()
            .into_pyarray(py)
    }

    /// Additional scan data recorded by the writer, returned as a dictionary.
    ///
    /// Row pitch, column pitch, and the along axis are exposed through attributes with those names
    /// and do not appear in this dictionary.
    #[getter]
    fn metadata<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new(py);
        for (k, v) in self.inner.metadata.iter() {
            out.set_item(k, value_to_py(py, v)?)?;
        }
        Ok(out)
    }

    /// The points of one row, as an `(n, 3)` array of dtype float64.
    fn row<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let row = self.inner.rows.get(index).ok_or_else(|| {
            PyValueError::new_err(format!(
                "Row {index} is out of range for a scan of {} rows",
                self.inner.rows.len()
            ))
        })?;
        Ok(points_to_array(&row.points).into_pyarray(py))
    }

    /// The column indices of one row, as an `(n,)` array of dtype uint32, or `None` if the scan
    /// has no column indices.
    fn row_columns<'py>(
        &self,
        py: Python<'py>,
        index: usize,
    ) -> PyResult<Option<Bound<'py, PyArray1<u32>>>> {
        let row = self.inner.rows.get(index).ok_or_else(|| {
            PyValueError::new_err(format!(
                "Row {index} is out of range for a scan of {} rows",
                self.inner.rows.len()
            ))
        })?;
        Ok(row.columns.as_ref().map(|c| c.clone().into_pyarray(py)))
    }

    fn __repr__(&self) -> String {
        format!(
            "RowPointsScan3({} rows, {} points, row_pitch={}, along={})",
            self.inner.rows.len(),
            self.inner.point_count(),
            self.inner.row_pitch,
            axis_name(self.inner.along)
        )
    }
}
