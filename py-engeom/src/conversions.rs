//! This module has conversion helpers for numpy arrays and other engeom types

use engeom::na::{Point, SVector};
use engeom::{Point2, Point3, SurfacePoint3, Vector2, Vector3};
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
/// Build a vector of oriented `SurfacePoint3` samples from a pair of Nx3 numpy arrays of point
/// positions and their corresponding normals. The two arrays must have the same number of rows;
/// each normal is normalized when its surface point is constructed.
pub fn array_to_surface_points3(
    points: &ArrayView2<'_, f64>,
    normals: &ArrayView2<'_, f64>,
) -> PyResult<Vec<SurfacePoint3>> {
    let ps = points.shape();
    let ns = normals.shape();
    if ps.len() != 2 || ps[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of points"));
    }
    if ns.len() != 2 || ns[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of normals"));
    }
    if ps[0] != ns[0] {
        return Err(PyValueError::new_err(
            "points and normals must have the same number of rows",
        ));
    }

    Ok(points
        .rows()
        .into_iter()
        .zip(normals.rows())
        .map(|(p, n)| {
            SurfacePoint3::new_normalize(
                Point3::new(p[0], p[1], p[2]),
                Vector3::new(n[0], n[1], n[2]),
            )
        })
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

// ================================================================================================
// Per-element attribute arrays
// ================================================================================================

/// Convert 8-bit RGB triples into an Nx3 numpy array.
pub fn colors_to_array(colors: &[[u8; 3]]) -> Array2<u8> {
    let mut array = Array2::zeros((colors.len(), 3));
    for (i, color) in colors.iter().enumerate() {
        for j in 0..3 {
            array[[i, j]] = color[j];
        }
    }
    array
}

/// Convert an Nx3 numpy array into 8-bit RGB triples.
pub fn array_to_colors(array: &ArrayView2<'_, u8>) -> PyResult<Vec<[u8; 3]>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of RGB colors"));
    }

    Ok(array
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

/// Convert an Nx3 numpy array into unit vectors, normalizing each row.
///
/// A row of zero length has no direction and is rejected rather than being silently turned into an
/// arbitrary one.
pub fn array_to_unit_vectors3(array: &ArrayView2<'_, f64>) -> PyResult<Vec<engeom::UnitVec3>> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("Expected Nx3 array of directions"));
    }

    array
        .rows()
        .into_iter()
        .enumerate()
        .map(|(i, row)| {
            let v = Vector3::new(row[0], row[1], row[2]);
            engeom::UnitVec3::try_new(v, 1.0e-12).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "Row {i} has zero length, which cannot be a direction"
                ))
            })
        })
        .collect()
}

/// Convert unit vectors into an Nx3 numpy array.
pub fn unit_vectors_to_array(vectors: &[engeom::UnitVec3]) -> Array2<f64> {
    let inner = vectors.iter().map(|v| v.into_inner()).collect::<Vec<_>>();
    vectors_to_array(&inner)
}

/// Convert a slice of scalars into a 1D numpy array.
pub fn scalars_to_array(values: &[f64]) -> Array1<f64> {
    Array1::from_iter(values.iter().copied())
}

/// Convert a slice of unsigned labels into a 1D numpy array.
pub fn labels_to_array(values: &[u32]) -> Array1<u32> {
    Array1::from_iter(values.iter().copied())
}

/// Copy a 1D numpy array into a `Vec`, rejecting a non-contiguous view rather than silently
/// reading the wrong elements.
pub fn array_to_vec<T: numpy::Element + Copy>(array: &PyReadonlyArray1<'_, T>) -> PyResult<Vec<T>> {
    Ok(array
        .as_slice()
        .map_err(|e| PyValueError::new_err(e.to_string()))?
        .to_vec())
}
