use pyo3::prelude::*;
use engeom::raster2::Point2I;
use crate::mesh::Mesh;

#[pyclass]
pub struct ScalarRaster {
    inner: engeom::raster2::ScalarRaster,
}

impl ScalarRaster {
    pub fn from_inner(inner: engeom::raster2::ScalarRaster) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::raster2::ScalarRaster {
        &self.inner
    }
}

#[pymethods]
impl ScalarRaster {

    #[getter]
    fn px_size(&self) -> f64 {
        self.inner.px_size
    }

    #[getter]
    fn min_z(&self) -> f64 {
        self.inner.min_z
    }

    #[getter]
    fn max_z(&self) -> f64 {
        self.inner.max_z
    }

    fn f_at(&self, x: i32, y: i32) -> f64 {
        self.inner.f_at(Point2I::new(x, y))
    }

    #[staticmethod]
    fn from_serialized_bytes(bytes: &[u8]) -> PyResult<Self> {
        let raster = engeom::raster2::ScalarRaster::from_serialized_bytes(bytes)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(raster))
    }
    
    fn build_depth_mesh(&self) -> PyResult<Mesh> {
        let inner = self.inner.build_depth_mesh()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Mesh::from_inner(inner))
    }

}