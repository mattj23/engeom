use crate::common::DeviationMode;
use crate::conversions::array_to_points3;
use crate::geom3::{Iso3, Point3};
use crate::mesh::Mesh;
use crate::point_cloud::PointCloud;
use engeom::geom3::align3::GAPParams;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(from_py_object, module="engeom.align")]
#[derive(Clone, Copy, Debug)]
pub struct Dof6 {
    #[pyo3(get, set)]
    pub tx: bool,
    #[pyo3(get, set)]
    pub ty: bool,
    #[pyo3(get, set)]
    pub tz: bool,
    #[pyo3(get, set)]
    pub rx: bool,
    #[pyo3(get, set)]
    pub ry: bool,
    #[pyo3(get, set)]
    pub rz: bool,
}

#[pymethods]
impl Dof6 {
    #[new]
    #[pyo3(signature = (tx=true, ty=true, tz=true, rx=true, ry=true, rz=true))]
    pub fn new(tx: bool, ty: bool, tz: bool, rx: bool, ry: bool, rz: bool) -> Self {
        Self { tx, ty, tz, rx, ry, rz }
    }

    #[staticmethod]
    pub fn all() -> Self {
        Self { tx: true, ty: true, tz: true, rx: true, ry: true, rz: true }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Dof6(tx={}, ty={}, tz={}, rx={}, ry={}, rz={})",
            self.tx, self.ty, self.tz, self.rx, self.ry, self.rz
        )
    }
}

impl From<Dof6> for engeom::geom3::align3::Dof6 {
    fn from(val: Dof6) -> Self {
        engeom::geom3::align3::Dof6::new(val.tx, val.ty, val.tz, val.rx, val.ry, val.rz)
    }
}


#[pyfunction]
#[pyo3(signature = (points, mesh, working=None, mode=None, center=None, dof=None))]
pub fn points_to_mesh(
    points: PyReadonlyArray2<'_, f64>,
    mesh: &Mesh,
    working: Option<Iso3>,
    mode: Option<DeviationMode>,
    center: Option<Point3>,
    dof: Option<Dof6>,
) -> PyResult<Iso3> {
    let points = array_to_points3(&points.as_array())?;

    let working = if let Some(w) = working {
        w.get_inner().clone()
    } else {
        engeom::Iso3::identity()
    };

    let dof = if let Some(d) = dof {
        d.into()
    } else {
        engeom::geom3::align3::Dof6::all()
    };

    let result = engeom::geom3::align3::points_to_mesh(
        &points,
        mesh.get_inner(),
        working,
        mode.unwrap_or_default().into(),
        center.map(|c| *c.get_inner()),
        dof,
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
    let points = array_to_points3(&points.as_array())?;

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

#[allow(clippy::too_many_arguments)]
#[pyfunction]
pub fn mesh_to_mesh_iterative(
    moving: &Mesh,
    reference: &Mesh,
    initial: &Iso3,
    mode: DeviationMode,
    max_iter: usize,
    sample_spacing: f64,
    max_neighbor_angle: f64,       // PI / 3.0
    out_of_plane_ratio: f64,       // 1 /20.0
    centroid_ratio: f64,           // 1.0
    filter_distances: Option<f64>, // Some(3.0)
) -> PyResult<Iso3> {
    let params = GAPParams::new(
        sample_spacing,
        max_neighbor_angle,
        out_of_plane_ratio,
        centroid_ratio,
        filter_distances,
    );
    let result = engeom::geom3::align3::mesh_to_mesh_iterative(
        moving.get_inner(),
        reference.get_inner(),
        initial.get_inner(),
        mode.into(),
        max_iter,
        &params,
    );

    match result {
        Ok(align) => Ok(Iso3::from_inner(*align.transform())),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
