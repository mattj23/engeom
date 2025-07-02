use crate::common::DeviationMode;
use crate::conversions::{array_to_points3, array2_to_points3};
use crate::geom3::Iso3;
use crate::mesh::Mesh;
use crate::point_cloud::PointCloud;
use engeom::PointCloudFeatures;
use numpy::{PyReadonlyArray2, PyReadonlyArrayDyn};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub fn points_to_mesh(
    points: PyReadonlyArrayDyn<'_, f64>,
    mesh: &Mesh,
    initial: &Iso3,
    mode: DeviationMode,
) -> PyResult<Iso3> {
    let points = array_to_points3(&points.as_array())?;

    let result = engeom::geom3::align3::points_to_mesh(
        &points,
        mesh.get_inner(),
        initial.get_inner(),
        mode.into(),
    );

    match result {
        Ok(align) => Ok(Iso3::from_inner(*align.transform())),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

#[pyfunction]
pub fn points_to_cloud(
    points: PyReadonlyArray2<'_, f64>,
    cloud: &mut PointCloud,
    search_radius: f64,
    initial: &Iso3,
) -> PyResult<Iso3> {
    let points = array2_to_points3(&points.as_array())?;

    // Must be a better way to do this
    let with_tree = cloud.with_tree()?;
    let result = engeom::geom3::align3::points_to_cloud(
        &points,
        &with_tree,
        search_radius,
        initial.get_inner(),
    );

    match result {
        Ok(align) => Ok(Iso3::from_inner(*align.transform())),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
