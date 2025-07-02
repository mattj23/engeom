use crate::conversions::{array_to_points2, array_to_points3};
use crate::geom2::{Iso2, Vector2};
use crate::geom3::{Iso3, Vector3};
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArrayDyn};
use pyo3::prelude::*;

#[pyclass]
pub struct SvdBasis2 {
    inner: engeom::SvdBasis2,
}

impl SvdBasis2 {
    pub fn get_inner(&self) -> &engeom::SvdBasis2 {
        &self.inner
    }
}

#[pymethods]
impl SvdBasis2 {
    #[new]
    #[pyo3(signature=(points, weights = None))]
    pub fn new<'py>(
        points: PyReadonlyArrayDyn<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points2(&points.as_array())?;

        let basis = match weights {
            Some(weights) => engeom::SvdBasis2::from_points(
                &points,
                Some(weights.as_array().as_slice().unwrap()),
            ),
            None => engeom::SvdBasis2::from_points(&points, None),
        };

        Ok(Self { inner: basis })
    }

    fn rank(&self, tol: f64) -> usize {
        self.inner.rank(tol)
    }

    fn largest(&self) -> Vector2 {
        let largest = self.inner.largest();
        Vector2::from_inner(engeom::Vector2::new(largest[0], largest[1]))
    }

    fn smallest(&self) -> Vector2 {
        let smallest = self.inner.smallest();
        Vector2::from_inner(engeom::Vector2::new(smallest[0], smallest[1]))
    }

    fn basis_variances<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut result = Array1::zeros(2);
        let variances = self.inner.basis_variances();
        result[0] = variances[0];
        result[1] = variances[1];
        result.into_pyarray(py)
    }

    fn basis_stdevs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut result = Array1::zeros(2);
        let stdevs = self.inner.basis_stdevs();
        result[0] = stdevs[0];
        result[1] = stdevs[1];
        result.into_pyarray(py)
    }

    fn to_iso2(&self) -> Iso2 {
        Iso2::from_inner((&self.inner).into())
    }
}

#[pyclass]
pub struct SvdBasis3 {
    inner: engeom::SvdBasis3,
}

impl SvdBasis3 {
    pub fn get_inner(&self) -> &engeom::SvdBasis3 {
        &self.inner
    }
}

#[pymethods]
impl SvdBasis3 {
    #[new]
    #[pyo3(signature=(points, weights = None))]
    pub fn new<'py>(
        points: PyReadonlyArrayDyn<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;

        // TODO: Is there some way to pass it back as a reference?
        let basis = match weights {
            Some(weights) => engeom::SvdBasis3::from_points(
                &points,
                Some(weights.as_array().as_slice().unwrap()),
            ),
            None => engeom::SvdBasis3::from_points(&points, None),
        };

        Ok(Self { inner: basis })
    }

    fn rank(&self, tol: f64) -> usize {
        self.inner.rank(tol)
    }

    fn largest(&self) -> Vector3 {
        let largest = self.inner.largest();
        Vector3::from_inner(engeom::Vector3::new(largest[0], largest[1], largest[2]))
    }

    fn smallest(&self) -> Vector3 {
        let smallest = self.inner.smallest();
        Vector3::from_inner(engeom::Vector3::new(smallest[0], smallest[1], smallest[2]))
    }

    fn basis_variances<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut result = Array1::zeros(3);
        let variances = self.inner.basis_variances();

        result[0] = variances[0];
        result[1] = variances[1];
        result[2] = variances[2];

        result.into_pyarray(py)
    }

    fn basis_stdevs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut result = Array1::zeros(3);
        let stdevs = self.inner.basis_stdevs();

        result[0] = stdevs[0];
        result[1] = stdevs[1];
        result[2] = stdevs[2];

        result.into_pyarray(py)
    }

    fn to_iso3(&self) -> Iso3 {
        Iso3::from_inner((&self.inner).into())
    }
}
