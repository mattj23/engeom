use crate::geom2::{Point2, Vector2};
use crate::geom3::{Point3, Vector3};
use numpy::{PyReadonlyArray2, PyReadonlyArrayDyn};
use pyo3::exceptions::PyValueError;
use pyo3::{PyResult, pyclass, pymethods};
// ================================================================================================
// Aabb2
// ================================================================================================

#[pyclass]
#[derive(Clone, Debug)]
pub struct Aabb2 {
    inner: engeom::geom2::Aabb2,
}

impl Aabb2 {
    pub fn get_inner(&self) -> &engeom::geom2::Aabb2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom2::Aabb2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Aabb2 {
    #[new]
    fn new(x_min: f64, y_min: f64, x_max: f64, y_max: f64) -> Self {
        // TODO: check min < max?
        Self {
            inner: engeom::geom2::Aabb2::new(
                engeom::Point2::new(x_min, y_min),
                engeom::Point2::new(x_max, y_max),
            ),
        }
    }

    #[staticmethod]
    fn from_points<'py>(points: PyReadonlyArrayDyn<'py, f64>) -> PyResult<Self> {
        let view = points.as_array();
        if view.shape().len() != 2 || view.shape()[1] != 2 || view.shape()[0] == 0 {
            return Err(PyValueError::new_err("Expected Nx2 array of points"));
        }

        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for i in 0..view.shape()[0] {
            let x = view[[i, 0]];
            let y = view[[i, 1]];
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
        }

        Ok(Self {
            inner: engeom::geom2::Aabb2::new(
                engeom::Point2::new(min_x, min_y),
                engeom::Point2::new(max_x, max_y),
            ),
        })
    }

    #[staticmethod]
    #[pyo3(signature=(x, y, w, h=None))]
    fn at_point(x: f64, y: f64, w: f64, h: Option<f64>) -> Self {
        let h = h.unwrap_or(w) / 2.0;
        let w = w / 2.0;
        let p = engeom::Point2::new(x, y);
        let v = engeom::Vector2::new(w, h);

        Self {
            inner: engeom::geom2::Aabb2::new(p - v, p + v),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Aabb2({}, {}, {}, {})",
            self.inner.mins.x, self.inner.mins.y, self.inner.maxs.x, self.inner.maxs.y,
        )
    }

    #[getter]
    fn min(&self) -> Point2 {
        Point2::from_inner(self.inner.mins.clone())
    }

    #[getter]
    fn max(&self) -> Point2 {
        Point2::from_inner(self.inner.maxs.clone())
    }

    #[getter]
    fn center(&self) -> Point2 {
        Point2::from_inner(self.inner.center())
    }

    #[getter]
    fn extent(&self) -> Vector2 {
        Vector2::from_inner(self.inner.extents())
    }

    fn expand(&self, d: f64) -> Self {
        use parry2d_f64::bounding_volume::BoundingVolume;
        Aabb2::from_inner(self.inner.loosened(d))
    }

    fn shrink(&self, d: f64) -> Self {
        use parry2d_f64::bounding_volume::BoundingVolume;
        Aabb2::from_inner(self.inner.tightened(d))
    }

    fn merged(&self, other: &Self) -> Self {
        use parry2d_f64::bounding_volume::BoundingVolume;
        let merged = self.get_inner().merged(other.get_inner());
        Aabb2::from_inner(merged)
    }

    fn contains_point(&self, point: &Point2) -> bool {
        self.inner.contains_local_point(point.get_inner())
    }

    fn indices_contained<'py>(&self, points: PyReadonlyArray2<'py, f64>) -> PyResult<Vec<usize>> {
        let view = points.as_array();
        if view.shape().len() != 2 || view.shape()[1] != 2 {
            return Err(PyValueError::new_err("Expected Nx2 array of points"));
        }

        let mut indices = Vec::new();
        for (i, point) in view.outer_iter().enumerate() {
            if self.contains_point(&Point2::from_inner(engeom::Point2::new(point[0], point[1]))) {
                indices.push(i);
            }
        }
        Ok(indices)
    }
}

// ================================================================================================
// Aabb3
// ================================================================================================

#[pyclass]
#[derive(Clone, Debug)]
pub struct Aabb3 {
    inner: engeom::geom3::Aabb3,
}

impl Aabb3 {
    pub fn get_inner(&self) -> &engeom::geom3::Aabb3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::Aabb3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Aabb3 {
    #[new]
    fn new(x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64) -> Self {
        // TODO: check min < max?
        Self {
            inner: engeom::geom3::Aabb3::new(
                engeom::Point3::new(x_min, y_min, z_min),
                engeom::Point3::new(x_max, y_max, z_max),
            ),
        }
    }

    #[staticmethod]
    fn from_points<'py>(points: PyReadonlyArrayDyn<'py, f64>) -> PyResult<Self> {
        let view = points.as_array();
        if view.shape().len() != 2 || view.shape()[1] != 3 || view.shape()[0] == 0 {
            return Err(PyValueError::new_err("Expected Nx3 array of points"));
        }

        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut min_z = f64::INFINITY;
        let mut max_z = f64::NEG_INFINITY;

        for i in 0..view.shape()[0] {
            let x = view[[i, 0]];
            let y = view[[i, 1]];
            let z = view[[i, 2]];
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
            min_z = min_z.min(z);
            max_z = max_z.max(z);
        }

        Ok(Self {
            inner: engeom::geom3::Aabb3::new(
                engeom::Point3::new(min_x, min_y, min_z),
                engeom::Point3::new(max_x, max_y, max_z),
            ),
        })
    }

    #[staticmethod]
    #[pyo3(signature=(x, y, z, w, h=None, d=None))]
    fn at_point(x: f64, y: f64, z: f64, w: f64, h: Option<f64>, d: Option<f64>) -> Self {
        let h = h.unwrap_or(w) / 2.0;
        let d = d.unwrap_or(w) / 2.0;
        let w = w / 2.0;
        let p = engeom::Point3::new(x, y, z);
        let v = engeom::Vector3::new(w, h, d);

        Self {
            inner: engeom::geom3::Aabb3::new(p - v, p + v),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Aabb3({}, {}, {}, {}, {}, {})",
            self.inner.mins.x,
            self.inner.mins.y,
            self.inner.mins.z,
            self.inner.maxs.x,
            self.inner.maxs.y,
            self.inner.maxs.z
        )
    }

    #[getter]
    fn min(&self) -> Point3 {
        Point3::from_inner(self.inner.mins.clone())
    }

    #[getter]
    fn max(&self) -> Point3 {
        Point3::from_inner(self.inner.maxs.clone())
    }

    #[getter]
    fn center(&self) -> Point3 {
        Point3::from_inner(self.inner.center())
    }

    #[getter]
    fn extent(&self) -> Vector3 {
        Vector3::from_inner(self.inner.extents())
    }

    fn expand(&self, d: f64) -> Self {
        use parry3d_f64::bounding_volume::BoundingVolume;
        Aabb3::from_inner(self.inner.loosened(d))
    }

    fn shrink(&self, d: f64) -> Self {
        use parry3d_f64::bounding_volume::BoundingVolume;
        Aabb3::from_inner(self.inner.tightened(d))
    }

    fn merged(&self, other: &Self) -> Self {
        use parry3d_f64::bounding_volume::BoundingVolume;
        let merged = self.get_inner().merged(other.get_inner());
        Aabb3::from_inner(merged)
    }

    fn contains_point(&self, point: &Point3) -> bool {
        self.inner.contains_local_point(point.get_inner())
    }

    fn indices_contained<'py>(&self, points: PyReadonlyArray2<'py, f64>) -> PyResult<Vec<usize>> {
        let view = points.as_array();
        if view.shape().len() != 2 || view.shape()[1] != 3 {
            return Err(PyValueError::new_err("Expected Nx3 array of points"));
        }

        let mut indices = Vec::new();
        for (i, point) in view.outer_iter().enumerate() {
            if self.contains_point(&Point3::from_inner(engeom::Point3::new(
                point[0], point[1], point[2],
            ))) {
                indices.push(i);
            }
        }
        Ok(indices)
    }
}
