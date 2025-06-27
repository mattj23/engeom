use std::time::Instant;
use engeom::sensor::SimulatedPointSensor;
use crate::mesh::Mesh;
use crate::geom3::{Iso3, Point3, Vector3};
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use numpy::ndarray::ArrayD;
use pyo3::exceptions::PyValueError;
use pyo3::{pyclass, pymethods, Bound, PyResult, Python};
use crate::conversions::points_to_array3;

#[pyclass]
#[derive(Clone)]
pub struct LaserLine {
    pub inner: engeom::sensor::LaserLine,
}

impl LaserLine {
    pub fn get_inner(&self) -> &engeom::sensor::LaserLine {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensor::LaserLine) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl LaserLine {
    #[new]
    #[pyo3(signature = (ray_origin, detect_origin, line_start, line_end, min_range, max_range, rays, angle_limit = None))]
    fn new(
        ray_origin: Point3,
        detect_origin: Point3,
        line_start: Point3,
        line_end: Point3,
        min_range: f64,
        max_range: f64,
        rays: usize,
        angle_limit: Option<f64>,
    ) -> PyResult<Self> {
        let inner = engeom::sensor::LaserLine::new(
            ray_origin.get_inner().clone(),
            detect_origin.get_inner().clone(),
            line_start.get_inner().clone(),
            line_end.get_inner().clone(),
            min_range,
            max_range,
            rays,
            angle_limit,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self { inner })
    }

    fn get_points<'py>(
        &self,
        py: Python<'py>,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let result = self
            .inner
            .get_points(target.get_inner(), obstruction.map(|o| o.get_inner()), iso.get_inner());
        Ok(points_to_array3(&result).into_pyarray(py))
    }
}


#[pyclass]
#[derive(Clone)]
pub struct PanningLaserLine {
    pub inner: engeom::sensor::PanningLaserLine,
}

impl PanningLaserLine {
    pub fn get_inner(&self) -> &engeom::sensor::PanningLaserLine {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensor::PanningLaserLine) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PanningLaserLine {
    #[new]
    fn new(
        laser_line: LaserLine,
        pan_vector: Vector3,
        steps: usize,
    ) -> PyResult<Self> {
        let inner = engeom::sensor::PanningLaserLine::new(
            laser_line.get_inner().clone(),
            pan_vector.get_inner().clone(),
            steps,
        );

        Ok(Self { inner })
    }

    fn get_points<'py>(
        &self,
        py: Python<'py>,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        // let start = Instant::now();
        let result = self
            .inner
            .get_points(target.get_inner(), obstruction.map(|o| o.get_inner()), iso.get_inner());
        Ok(points_to_array3(&result).into_pyarray(py))
    }
}