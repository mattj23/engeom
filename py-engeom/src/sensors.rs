use crate::geom3::Iso3;
use crate::mesh::Mesh;
use crate::point_cloud::PointCloud;
use engeom::sensors::SimulatedPointSensor;
use pyo3::{PyResult, pyclass, pymethods};

#[pyclass]
#[derive(Clone)]
pub struct LaserProfile {
    pub inner: engeom::sensors::LaserProfile,
}

impl LaserProfile {
    pub fn get_inner(&self) -> &engeom::sensors::LaserProfile {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::LaserProfile) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl LaserProfile {
    #[new]
    #[pyo3(signature = (emitter_z, detector_y, detector_z, volume_width, volume_z_min, volume_z_max, resolution, angle_limit = None))]
    fn new(
        emitter_z: f64,
        detector_y: f64,
        detector_z: f64,
        volume_width: f64,
        volume_z_min: f64,
        volume_z_max: f64,
        resolution: usize,
        angle_limit: Option<f64>,
    ) -> Self {
        let inner = engeom::sensors::LaserProfile::new(
            emitter_z,
            detector_y,
            detector_z,
            volume_width,
            volume_z_min,
            volume_z_max,
            resolution,
            angle_limit,
        );

        Self { inner }
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
pub struct PanningLaserProfile {
    pub inner: engeom::sensors::PanningLaserProfile,
}

impl PanningLaserProfile {
    pub fn get_inner(&self) -> &engeom::sensors::PanningLaserProfile {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::PanningLaserProfile) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PanningLaserProfile {
    #[new]
    fn new(laser_line: LaserProfile, y_step: f64, steps: usize) -> PyResult<Self> {
        let inner = engeom::sensors::PanningLaserProfile::new(
            laser_line.get_inner().clone(),
            y_step,
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
