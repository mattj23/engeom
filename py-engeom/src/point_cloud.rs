use crate::bounding::Aabb3;
use crate::common::{DeviationMode, SelectOp};
use crate::conversions::{
    array_to_faces, array_to_points3, array2_to_points3, faces_to_array, points_to_array3,
    points_to_array3_2, vectors_to_array3, vectors_to_array3_2,
};
use crate::geom3::{Curve3, Iso3, Plane3, Point3, SurfacePoint3, Vector3};
use crate::mesh::Mesh;
use crate::metrology::Distance3;
use engeom::PointCloudFeatures;
use engeom::common::points::dist;
use engeom::common::{Selection, SplitResult};
use numpy::ndarray::{Array1, Array2, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayDyn, PyReadonlyArray2, PyReadonlyArrayDyn};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

#[pyclass]
pub struct PointCloud {
    inner: engeom::PointCloud,
    points: Option<Py<PyArray2<f64>>>,
    normals: Option<Py<PyArray2<f64>>>,
    colors: Option<Py<PyArray2<u8>>>,
}

impl PointCloud {
    fn clear_cached(&mut self) {
        self.points = None;
        self.normals = None;
        self.colors = None;
    }

    pub fn get_inner(&self) -> &engeom::PointCloud {
        &self.inner
    }

    pub fn from_inner(inner: engeom::PointCloud) -> Self {
        Self {
            inner,
            points: None,
            normals: None,
            colors: None,
        }
    }
}

impl Clone for PointCloud {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl PointCloud {
    #[new]
    fn new<'py>(points: PyReadonlyArray2<'py, f64>) -> PyResult<Self> {
        let cloud_points = array2_to_points3(&points.as_array())?;
        let cloud = engeom::PointCloud::try_new(cloud_points, None, None)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(cloud))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    fn append(&mut self, other: &PointCloud) -> PyResult<()> {
        self.clear_cached();
        let clone = other.inner.clone();
        self.inner
            .merge(clone)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let array = points_to_array3_2(self.inner.points());
            self.points = Some(array.into_pyarray(py).unbind());
        }
        self.points.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn colors<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<u8>>> {
        if let Some(colors) = self.inner.colors() {
            if self.colors.is_none() {
                let flat_colors = colors.iter().flatten().map(|u| *u).collect::<Vec<_>>();
                let array = Array2::from_shape_vec((self.inner.points().len(), 3), flat_colors)
                    .expect("Failed to create color array");
                self.colors = Some(array.into_pyarray(py).unbind());
            }

            Some(self.colors.as_ref().unwrap().bind(py))
        } else {
            None
        }
    }

    #[getter]
    fn normals<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<f64>>> {
        if let Some(normals) = self.inner.normals() {
            if self.normals.is_none() {
                let n = normals.iter().map(|v| v.into_inner()).collect::<Vec<_>>();
                let array = vectors_to_array3_2(&n);
                self.normals = Some(array.into_pyarray(py).unbind());
            }

            Some(self.normals.as_ref().unwrap().bind(py))
        } else {
            None
        }
    }
}
