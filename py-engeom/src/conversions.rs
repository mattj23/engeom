//! This module has conversion helpers for numpy arrays and other engeom types

use engeom::na::{Point, SVector};
use engeom::{Point2, Point3, Vector2, Vector3};
use numpy::ndarray::{Array1, Array2, ArrayView2};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::{Bound, PyResult, Python};

/// Build an engeom `DVector` (the dynamically sized column vector used for fitting parameters)
/// from a 1D numpy array.
pub fn dvec_from_array(arr: &PyReadonlyArray1<f64>) -> PyResult<engeom::DVector> {
    let s = arr
        .as_slice()
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(engeom::DVector::from_column_slice(s))
}

/// Convert an engeom `DVector` into a 1D numpy array.
pub fn dvec_to_array<'py>(py: Python<'py>, v: engeom::DVector) -> Bound<'py, PyArray1<f64>> {
    Array1::from_iter(v.iter().copied()).into_pyarray(py)
}

pub fn write_v<const D: usize>(a: &mut Array2<f64>, i: usize, v: &SVector<f64, D>) {
    for j in 0..D {
        a[[i, j]] = v[j];
    }
}

pub fn points_to_array<const D: usize>(points: &[Point<f64, D>]) -> Array2<f64> {
    let mut array = Array2::zeros((points.len(), D));
    for (i, point) in points.iter().enumerate() {
        write_v(&mut array, i, &point.coords);
    }
    array
}

pub fn vectors_to_array<const D: usize>(vectors: &[SVector<f64, D>]) -> Array2<f64> {
    let mut array = Array2::zeros((vectors.len(), D));
    for (i, vector) in vectors.iter().enumerate() {
        write_v(&mut array, i, vector);
    }
    array
}

pub fn array_to_points3(array: &ArrayView2<'_, f64>) -> PyResult<Vec<Point3>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of points"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| Point3::new(row[0], row[1], row[2]))
        .collect())
}
pub fn array_to_vectors3(array: &ArrayView2<'_, f64>) -> PyResult<Vec<Vector3>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of vectors"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| Vector3::new(row[0], row[1], row[2]))
        .collect())
}

pub fn array_to_points2(array: &ArrayView2<'_, f64>) -> PyResult<Vec<Point2>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 2 {
        return Err(PyValueError::new_err("Expected Nx2 array of points"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| Point2::new(row[0], row[1]))
        .collect())
}

pub fn array_to_vectors2(array: &ArrayView2<'_, f64>) -> PyResult<Vec<Vector2>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 2 {
        return Err(PyValueError::new_err("Expected Nx2 array of vectors"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| Vector2::new(row[0], row[1]))
        .collect())
}

pub fn faces_to_array(faces: &[[u32; 3]]) -> Array2<u32> {
    let mut array = Array2::zeros((faces.len(), 3));
    for (i, face) in faces.iter().enumerate() {
        array[[i, 0]] = face[0];
        array[[i, 1]] = face[1];
        array[[i, 2]] = face[2];
    }
    array
}

pub fn array_to_faces(array: &ArrayView2<'_, u32>) -> PyResult<Vec<[u32; 3]>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of faces"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}
