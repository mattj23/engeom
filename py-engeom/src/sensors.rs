use crate::conversions::points_to_array3;
use crate::geom3::{Iso3, Point3, Vector3};
use crate::mesh::Mesh;
use crate::point_cloud::PointCloud;
use engeom::PointCloudFeatures;
use engeom::sensors::SimulatedPointSensor;
use numpy::ndarray::ArrayD;
use numpy::{IntoPyArray, PyArrayDyn, PyReadonlyArrayDyn};
use pyo3::exceptions::PyValueError;
use pyo3::{Bound, PyResult, Python, pyclass, pymethods};
use std::time::Instant;

#[pyclass]
#[derive(Clone)]
pub struct LaserLine {
    pub inner: engeom::sensors::LaserLine,
}

impl LaserLine {
    pub fn get_inner(&self) -> &engeom::sensors::LaserLine {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::LaserLine) -> Self {
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
        let inner = engeom::sensors::LaserLine::new(
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
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> PyResult<PointCloud> {
        let (cloud, _) = self.inner.get_points(
            target.get_inner(),
            obstruction.map(|o| o.get_inner()),
            iso.get_inner(),
        );
        Ok(PointCloud::from_inner(cloud))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PanningLaserLine {
    pub inner: engeom::sensors::PanningLaserLine,
}

impl PanningLaserLine {
    pub fn get_inner(&self) -> &engeom::sensors::PanningLaserLine {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::PanningLaserLine) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PanningLaserLine {
    #[new]
    fn new(laser_line: LaserLine, pan_vector: Vector3, steps: usize) -> PyResult<Self> {
        let inner = engeom::sensors::PanningLaserLine::new(
            laser_line.get_inner().clone(),
            pan_vector.get_inner().clone(),
            steps,
        );

        Ok(Self { inner })
    }

    fn get_points<'py>(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> PyResult<PointCloud> {
        let (cloud, _) = self.inner.get_points(
            target.get_inner(),
            obstruction.map(|o| o.get_inner()),
            iso.get_inner(),
        );
        Ok(PointCloud::from_inner(cloud))
    }
}
