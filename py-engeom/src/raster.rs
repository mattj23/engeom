use numpy::ndarray::{Array1, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArrayDyn, PyReadonlyArrayDyn, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyResult, Python, pyclass,
    pyfunction,
};
use std::collections::HashSet;

#[pyfunction]
pub fn clusters_from_sparse<'py>(
    indices: PyReadonlyArrayDyn<'py, i32>,
) -> PyResult<Vec<Bound<'py, PyArrayDyn<i32>>>> {
    let shape = indices.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of indices"));
    }

    let mut working = HashSet::new();
    for r in indices.as_array().rows().into_iter() {
        let idx = (r[0], r[1], r[2]);
        working.insert(idx);
    }

    let mut results = engeom::raster3::clusters_from_sparse(working);
    let mut combined = Vec::new();

    for result in results.drain(..) {
        let mut array = ArrayD::zeros(vec![result.len(), 3]);
        for (i, idx) in result.iter().enumerate() {
            array[[i, 0]] = idx.0;
            array[[i, 1]] = idx.1;
            array[[i, 2]] = idx.2;
        }
        combined.push(array.into_pyarray(indices.py()));
    }

    Ok(combined)
}
