use crate::common::DeviationMode;
use crate::conversions::array_to_points3;
use crate::geom3::{Iso3, Point3};
use crate::mesh::Mesh;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::{PyResult, pyclass, pyfunction, pymethods};

#[pyclass(from_py_object, module = "engeom.align3")]
#[derive(Clone, Debug)]
pub struct AlignParams3 {
    inner: engeom::geom3::align3::AlignParams3,
}

impl AlignParams3 {
    pub fn from_inner(inner: engeom::geom3::align3::AlignParams3) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::geom3::align3::AlignParams3 {
        &self.inner
    }
}

#[pyclass(from_py_object, module = "engeom.align3")]
#[derive(Clone, Debug)]
pub struct Alignment3 {
    inner: engeom::geom3::Alignment3,
}

impl Alignment3 {
    pub fn from_inner(inner: engeom::geom3::Alignment3) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::geom3::Alignment3 {
        &self.inner
    }
}

#[pyclass(from_py_object, module = "engeom.align3")]
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
        Self {
            tx,
            ty,
            tz,
            rx,
            ry,
            rz,
        }
    }

    #[staticmethod]
    pub fn all() -> Self {
        Self {
            tx: true,
            ty: true,
            tz: true,
            rx: true,
            ry: true,
            rz: true,
        }
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
#[pyo3(signature = (points, mesh, params))]
pub fn points_to_mesh(
    points: PyReadonlyArray2<'_, f64>,
    mesh: &Mesh,
    params: AlignParams3,
) -> PyResult<Alignment3> {
    let points = array_to_points3(&points.as_array())?;

    let result = engeom::geom3::align3::points_to_mesh(
        &points,
        mesh.get_inner(),
        params.get_inner().clone(),
    );

    match result {
        Ok(align) => Ok(Alignment3::from_inner(align)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
