use crate::geom3::Iso3;
use crate::mesh::Mesh;
use numpy::ndarray::{Array1, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArrayDyn, PyReadonlyArrayDyn, PyUntypedArrayMethods};
use parry3d_f64::query::{Ray as PRay3, RayCast};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyResult, Python, pyclass,
    pymethods,
};

#[pyclass]
#[derive(Clone)]
pub struct RayBundle3 {
    inner: Vec<PRay3>,
}

impl RayBundle3 {
    pub fn from_inner(inner: Vec<PRay3>) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &Vec<PRay3> {
        &self.inner
    }
}

#[pymethods]
impl RayBundle3 {
    #[new]
    fn new<'py>(array: PyReadonlyArrayDyn<'py, f64>) -> PyResult<Self> {
        let array = array.as_array();
        let shape = array.shape();
        if shape.len() != 2 || shape[1] != 6 {
            return Err(PyValueError::new_err("Expected Nx6 array of points"));
        }

        let mut rays = Vec::with_capacity(shape[0]);

        for row in array.rows().into_iter() {
            let ray = PRay3::new(
                engeom::Point3::new(row[0], row[1], row[2]),
                engeom::Vector3::new(row[3], row[4], row[5]),
            );
            rays.push(ray);
        }

        Ok(Self::from_inner(rays))
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!("<RayBundle3 n={}>", self.inner.len())
    }

    #[pyo3(signature=(mesh, mesh_iso = None, angle = None))]
    fn intersect_mesh<'py>(
        &self,
        py: Python<'py>,
        mesh: &Mesh,
        mesh_iso: Option<&Iso3>,
        angle: Option<f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let mut result = Vec::new();
        let iso = if let Some(mesh_iso) = mesh_iso {
            mesh_iso.get_inner()
        } else {
            &engeom::Iso3::identity()
        };

        if let Some(angle_limit) = angle {
            for ray in self.inner.iter() {
                if let Some(ri) =
                    mesh.get_inner()
                        .tri_mesh()
                        .cast_ray_and_get_normal(iso, ray, f64::MAX, false)
                {
                    let n = ri.normal * -1.0;
                    if ray.dir.angle(&n) < angle_limit {
                        let t = ri.time_of_impact;
                        result.push(ray.point_at(t))
                    }
                }
            }
        } else {
            for ray in self.inner.iter() {
                if let Some(t) = mesh
                    .get_inner()
                    .tri_mesh()
                    .cast_ray(iso, ray, f64::MAX, false)
                {
                    result.push(ray.point_at(t))
                }
            }
        }

        let mut result_array = ArrayD::zeros(vec![result.len(), 3]);
        for (i, point) in result.iter().enumerate() {
            result_array[[i, 0]] = point.x;
            result_array[[i, 1]] = point.y;
            result_array[[i, 2]] = point.z;
        }

        Ok(result_array.into_pyarray(py))
    }
}
