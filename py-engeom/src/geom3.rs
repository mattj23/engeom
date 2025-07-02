use crate::common::Resample;
use crate::conversions::{array_to_points3, array_to_vectors3, points_to_array3};
use engeom::geom3::IsoExtensions3;
use numpy::ndarray::{Array1, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArrayDyn, PyReadonlyArrayDyn, PyUntypedArrayMethods};
use parry3d_f64::na::{Translation3, UnitQuaternion};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyResult, Python, pyclass,
    pymethods,
};

#[derive(FromPyObject)]
enum Vector3OrPoint3 {
    Vector(Vector3),
    Point(Point3),
}

// ================================================================================================
// Vectors
// ================================================================================================

#[pyclass]
#[derive(Clone, Debug)]
pub struct Vector3 {
    inner: engeom::Vector3,
}

impl Vector3 {
    pub fn get_inner(&self) -> &engeom::Vector3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Vector3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Vector3 {
    #[new]
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: engeom::Vector3::new(x, y, z),
        }
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    #[getter]
    fn z(&self) -> f64 {
        self.inner.z
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let o = [self.inner.x, self.inner.y, self.inner.z];
        PyIterator::from_object(&o.into_pyobject(py)?)
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut array = Array1::zeros(3);
        array[0] = self.inner.x;
        array[1] = self.inner.y;
        array[2] = self.inner.z;
        array.into_pyarray(py)
    }

    fn __neg__(&self) -> Self {
        Self { inner: -self.inner }
    }

    fn __mul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self {
            inner: self.inner / other,
        }
    }

    fn __rmul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __add__<'py>(&self, py: Python<'py>, other: Vector3OrPoint3) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Vector3OrPoint3::Vector(other) => {
                Vector3::from_inner(self.inner + other.inner).into_bound_py_any(py)
            }
            Vector3OrPoint3::Point(other) => {
                Point3::from_inner((self.inner + other.inner.coords).into()).into_bound_py_any(py)
            }
        }
    }

    fn __sub__(&self, other: Vector3) -> Self {
        Self::from_inner(self.inner - other.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "Vector3({}, {}, {})",
            self.inner.x, self.inner.y, self.inner.z
        )
    }

    fn dot(&self, other: Vector3) -> f64 {
        self.inner.dot(&other.inner)
    }

    fn cross(&self, other: Vector3) -> Self {
        Self::from_inner(self.inner.cross(&other.inner))
    }

    fn norm(&self) -> f64 {
        self.inner.norm()
    }

    fn normalized(&self) -> Self {
        Self {
            inner: self.inner.normalize(),
        }
    }

    fn angle_to(&self, other: Vector3) -> f64 {
        self.inner.angle(&other.inner)
    }

    fn with_z(&self, z: f64) -> Self {
        Self {
            inner: engeom::Vector3::new(self.inner.x, self.inner.y, z),
        }
    }

    fn with_x(&self, x: f64) -> Self {
        Self {
            inner: engeom::Vector3::new(x, self.inner.y, self.inner.z),
        }
    }

    fn with_y(&self, y: f64) -> Self {
        Self {
            inner: engeom::Vector3::new(self.inner.x, y, self.inner.z),
        }
    }
}

// ================================================================================================
// Points
// ================================================================================================

#[pyclass]
#[derive(Clone, Debug)]
pub struct Point3 {
    inner: engeom::Point3,
}

impl Point3 {
    pub fn get_inner(&self) -> &engeom::Point3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Point3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Point3 {
    #[new]
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: engeom::Point3::new(x, y, z),
        }
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    #[getter]
    fn z(&self) -> f64 {
        self.inner.z
    }

    #[getter]
    fn coords(&self) -> Vector3 {
        Vector3 {
            inner: self.inner.coords,
        }
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let o = [self.inner.x, self.inner.y, self.inner.z];
        PyIterator::from_object(&o.into_pyobject(py)?)
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut array = Array1::zeros(3);
        array[0] = self.inner.x;
        array[1] = self.inner.y;
        array[2] = self.inner.z;
        array.into_pyarray(py)
    }

    fn __add__(&self, other: Vector3) -> Self {
        Self::from_inner(self.inner + other.inner)
    }

    fn __sub__<'py>(&self, py: Python<'py>, other: Vector3OrPoint3) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Vector3OrPoint3::Vector(other) => {
                Point3::from_inner(self.inner - other.inner).into_bound_py_any(py)
            }
            Vector3OrPoint3::Point(other) => {
                Vector3::from_inner(self.inner - other.inner).into_bound_py_any(py)
            }
        }
    }

    fn __neg__(&self) -> Self {
        Self { inner: -self.inner }
    }

    fn __mul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self {
            inner: self.inner / other,
        }
    }

    fn __rmul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Point3({}, {}, {})",
            self.inner.x, self.inner.y, self.inner.z
        )
    }

    #[staticmethod]
    fn mid(a: Point3, b: Point3) -> Self {
        Self::from_inner(engeom::common::points::mid_point(
            a.get_inner(),
            b.get_inner(),
        ))
    }

    fn with_x(&self, x: f64) -> Self {
        Self {
            inner: engeom::Point3::new(x, self.inner.y, self.inner.z),
        }
    }

    fn with_y(&self, y: f64) -> Self {
        Self {
            inner: engeom::Point3::new(self.inner.x, y, self.inner.z),
        }
    }

    fn with_z(&self, z: f64) -> Self {
        Self {
            inner: engeom::Point3::new(self.inner.x, self.inner.y, z),
        }
    }
}

// ================================================================================================
// Surface Point
// ================================================================================================
#[pyclass]
#[derive(Clone, Debug)]
pub struct SurfacePoint3 {
    pub inner: engeom::SurfacePoint3,
}

impl SurfacePoint3 {
    pub fn get_inner(&self) -> &engeom::SurfacePoint3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::SurfacePoint3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl SurfacePoint3 {
    #[new]
    fn new(x: f64, y: f64, z: f64, nx: f64, ny: f64, nz: f64) -> Self {
        Self {
            inner: engeom::SurfacePoint3::new_normalize(
                engeom::Point3::new(x, y, z),
                engeom::Vector3::new(nx, ny, nz),
            ),
        }
    }

    #[getter]
    fn point(&self) -> Point3 {
        Point3::from_inner(self.inner.point)
    }

    #[getter]
    fn normal(&self) -> Vector3 {
        Vector3::from_inner(self.inner.normal.into_inner())
    }

    fn at_distance(&self, distance: f64) -> Point3 {
        Point3::from_inner(self.inner.at_distance(distance))
    }

    fn scalar_projection(&self, other: Point3) -> f64 {
        self.inner.scalar_projection(other.get_inner())
    }

    fn projection(&self, other: Point3) -> Point3 {
        Point3::from_inner(self.inner.projection(other.get_inner()))
    }

    fn reversed(&self) -> Self {
        Self::from_inner(self.inner.reversed())
    }

    fn transformed(&self, iso: Iso3) -> Self {
        Self::from_inner(self.inner.transformed(iso.get_inner()))
    }

    fn __repr__(&self) -> String {
        format!(
            "SurfacePoint3({}, {}, {}, {}, {}, {})",
            self.inner.point.x,
            self.inner.point.y,
            self.inner.point.z,
            self.inner.normal.x,
            self.inner.normal.y,
            self.inner.normal.z
        )
    }

    fn planar_distance(&self, other: Point3) -> f64 {
        self.inner.planar_distance(other.get_inner())
    }

    fn __mul__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint3::new_normalize(
            self.inner.point * other,
            self.inner.normal.into_inner() * other.signum(),
        ))
    }

    fn __rmul__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint3::new_normalize(
            self.inner.point * other,
            self.inner.normal.into_inner() * other.signum(),
        ))
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint3::new_normalize(
            self.inner.point / other,
            self.inner.normal.into_inner() / other.signum(),
        ))
    }

    fn __neg__(&self) -> Self {
        Self::from_inner(engeom::SurfacePoint3::new_normalize(
            -self.inner.point,
            -self.inner.normal.into_inner(),
        ))
    }

    fn get_plane(&self) -> Plane3 {
        Plane3::from_inner(engeom::Plane3::from(&self.inner))
    }

    fn shift(&self, offset: f64) -> Self {
        Self::from_inner(self.inner.shift(offset))
    }
}

// ================================================================================================
// Plane
// ================================================================================================
#[pyclass]
#[derive(Clone, Debug)]
pub struct Plane3 {
    pub inner: engeom::Plane3,
}

impl Plane3 {
    pub fn get_inner(&self) -> &engeom::Plane3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Plane3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Plane3 {
    #[new]
    fn new(a: f64, b: f64, c: f64, d: f64) -> PyResult<Self> {
        let v = engeom::Vector3::new(a, b, c);
        let normal = engeom::UnitVec3::try_new(v, 1.0e-6)
            .ok_or(PyValueError::new_err("Invalid normal vector"))?;

        Ok(Self {
            inner: engeom::Plane3::new(normal, d),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Plane3({}, {}, {}, {})",
            self.inner.normal.x, self.inner.normal.y, self.inner.normal.z, self.inner.d
        )
    }

    fn inverted_normal(&self) -> Self {
        Self::from_inner(self.inner.inverted_normal())
    }

    fn signed_distance_to_point(&self, point: Point3) -> f64 {
        self.inner.signed_distance_to_point(point.get_inner())
    }

    fn project_point(&self, point: Point3) -> Point3 {
        Point3::from_inner(self.inner.project_point(point.get_inner()))
    }

    fn intersection_distance(&self, sp: &SurfacePoint3) -> Option<f64> {
        self.inner.intersection_distance(sp.get_inner())
    }
}

// ================================================================================================
// Curve
// ================================================================================================
#[pyclass]
#[derive(Clone, Debug)]
pub struct CurveStation3 {
    i_point: engeom::Point3,
    i_direction: engeom::Vector3,
    i_index: usize,
    i_fraction: f64,
    i_length_along: f64,
}

impl CurveStation3 {
    pub fn new(
        point: engeom::Point3,
        direction: engeom::Vector3,
        index: usize,
        fraction: f64,
        length_along: f64,
    ) -> Self {
        Self {
            i_point: point,
            i_direction: direction,
            i_index: index,
            i_fraction: fraction,
            i_length_along: length_along,
        }
    }
}

#[pymethods]
impl CurveStation3 {
    #[getter]
    pub fn point(&self) -> Point3 {
        Point3::from_inner(self.i_point)
    }

    #[getter]
    pub fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.i_direction)
    }

    #[getter]
    pub fn direction_point(&self) -> SurfacePoint3 {
        SurfacePoint3::from_inner(engeom::SurfacePoint3::new_normalize(
            self.i_point,
            self.i_direction,
        ))
    }

    #[getter]
    pub fn index(&self) -> usize {
        self.i_index
    }

    #[getter]
    pub fn fraction(&self) -> f64 {
        self.i_fraction
    }

    #[getter]
    pub fn length_along(&self) -> f64 {
        self.i_length_along
    }
}

#[pyclass]
pub struct Curve3 {
    inner: engeom::Curve3,
    points: Option<Py<PyArrayDyn<f64>>>,
}

impl Clone for Curve3 {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            points: None,
        }
    }
}

impl Curve3 {
    pub fn get_inner(&self) -> &engeom::Curve3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Curve3) -> Self {
        Self {
            inner,
            points: None,
        }
    }
}

#[pymethods]
impl Curve3 {
    #[new]
    #[pyo3(signature=(points, tol=1.0e-6))]
    fn new(points: PyReadonlyArrayDyn<'_, f64>, tol: f64) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let inner = engeom::Curve3::from_points(&points, tol)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArrayDyn<f64>> {
        if self.points.is_none() {
            let result = points_to_array3(self.inner.vertices());
            self.points = Some(result.into_pyarray(py).unbind())
        }

        self.points.as_ref().unwrap().bind(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "<Curve3 {} points, {} long>",
            self.inner.vertices().len(),
            self.inner.length()
        )
    }

    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }

    fn length(&self) -> f64 {
        self.inner.length()
    }

    fn at_length(&self, length: f64) -> PyResult<CurveStation3> {
        let station = self
            .inner
            .at_length(length)
            .ok_or(PyValueError::new_err("Invalid length"))?;
        Ok(station.into())
    }

    fn at_fraction(&self, fraction: f64) -> PyResult<CurveStation3> {
        self.at_length(fraction * self.inner.length())
    }

    fn at_closest_to_point(&self, point: Point3) -> PyResult<CurveStation3> {
        let station = self.inner.at_closest_to_point(point.get_inner());
        Ok(station.into())
    }

    fn at_front(&self) -> CurveStation3 {
        self.inner.at_front().into()
    }

    fn at_back(&self) -> CurveStation3 {
        self.inner.at_back().into()
    }

    fn resample(&self, resample: Resample) -> Self {
        Self::from_inner(self.inner.resample(resample.into()))
    }

    fn simplify(&self, tol: f64) -> Self {
        Self::from_inner(self.inner.simplify(tol))
    }

    fn transformed_by(&mut self, iso: Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }
}

impl From<engeom::CurveStation3<'_>> for CurveStation3 {
    fn from(station: engeom::CurveStation3) -> Self {
        Self::new(
            station.point(),
            station.direction().into_inner(),
            station.index(),
            station.fraction(),
            station.length_along(),
        )
    }
}

// ================================================================================================
// Transformations
// ================================================================================================

#[pyclass]
#[derive(Clone, Debug)]
pub struct XyzWpr {
    inner: engeom::geom3::XyzWpr,
}

impl XyzWpr {
    pub fn get_inner(&self) -> &engeom::geom3::XyzWpr {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::XyzWpr) -> Self {
        Self { inner }
    }
}

#[derive(FromPyObject)]
enum Transformable3 {
    Iso(Iso3),
    Vec(Vector3),
    Pnt(Point3),
    Plane(Plane3),
    Sp(SurfacePoint3),
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct Iso3 {
    inner: engeom::Iso3,
}

impl Iso3 {
    pub fn get_inner(&self) -> &engeom::Iso3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Iso3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Iso3 {
    fn __repr__(&self) -> String {
        format!(
            "<Iso3 t=[{}, {}, {}] r=[{}, {}, {}, {}]>",
            self.inner.translation.x,
            self.inner.translation.y,
            self.inner.translation.z,
            self.inner.rotation.i,
            self.inner.rotation.j,
            self.inner.rotation.k,
            self.inner.rotation.w,
        )
    }

    #[getter]
    fn origin(&self) -> Point3 {
        Point3::from_inner(engeom::Point3::new(
            self.inner.translation.x,
            self.inner.translation.y,
            self.inner.translation.z,
        ))
    }

    #[staticmethod]
    fn from_xyzwpr(x: f64, y: f64, z: f64, w: f64, p: f64, r: f64) -> Self {
        Self {
            inner: engeom::geom3::XyzWpr::new(x, y, z, w, p, r).to_isometry(),
        }
    }

    fn to_xyzwpr(&self) -> Vec<f64> {
        let v = engeom::geom3::XyzWpr::from_isometry(&self.inner);
        vec![v.x, v.y, v.z, v.w, v.p, v.r]
    }

    #[new]
    fn new(matrix: PyReadonlyArrayDyn<'_, f64>) -> PyResult<Self> {
        if matrix.shape().len() != 2 || matrix.shape()[0] != 4 || matrix.shape()[1] != 4 {
            return Err(PyValueError::new_err("Expected 4x4 matrix"));
        }

        let mut array = [0.0; 16];
        for (i, value) in matrix.as_array().iter().enumerate() {
            array[i] = *value;
        }

        let inner = engeom::Iso3::try_from_array(&array)
            .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner })
    }

    #[staticmethod]
    fn from_translation(x: f64, y: f64, z: f64) -> Self {
        Self {
            inner: engeom::Iso3::translation(x, y, z),
        }
    }

    #[staticmethod]
    fn from_rotation(angle: f64, a: f64, b: f64, c: f64) -> Self {
        let axis = engeom::UnitVec3::new_normalize(engeom::Vector3::new(a, b, c));
        let rot_vec = axis.into_inner() * angle;

        Self {
            inner: engeom::Iso3::rotation(rot_vec),
        }
    }

    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse(),
        }
    }

    fn __matmul__<'py>(
        &self,
        py: Python<'py>,
        other: Transformable3,
    ) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Transformable3::Iso(other) => {
                Iso3::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable3::Vec(other) => {
                Vector3::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable3::Pnt(other) => {
                Point3::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable3::Plane(other) => {
                Plane3::from_inner(other.inner.transform_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Sp(other) => {
                SurfacePoint3::from_inner(other.inner.transformed(&self.inner))
                    .into_bound_py_any(py)
            }
        }
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        let mut result = ArrayD::zeros(vec![4, 4]);
        let m = self.inner.to_matrix();
        // TODO: In a rush, fix this later
        result[[0, 0]] = m.m11;
        result[[0, 1]] = m.m12;
        result[[0, 2]] = m.m13;
        result[[0, 3]] = m.m14;
        result[[1, 0]] = m.m21;
        result[[1, 1]] = m.m22;
        result[[1, 2]] = m.m23;
        result[[1, 3]] = m.m24;
        result[[2, 0]] = m.m31;
        result[[2, 1]] = m.m32;
        result[[2, 2]] = m.m33;
        result[[2, 3]] = m.m34;
        result[[3, 0]] = m.m41;
        result[[3, 1]] = m.m42;
        result[[3, 2]] = m.m43;
        result[[3, 3]] = m.m44;
        result.into_pyarray(py)
    }

    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: engeom::Iso3::identity(),
        }
    }

    fn flip_around_x(&self) -> Self {
        Self {
            inner: self.inner.flip_around_x(),
        }
    }

    fn flip_around_y(&self) -> Self {
        Self {
            inner: self.inner.flip_around_y(),
        }
    }

    fn flip_around_z(&self) -> Self {
        Self {
            inner: self.inner.flip_around_z(),
        }
    }

    fn translation(&self) -> Iso3 {
        Self {
            inner: engeom::Iso3::from_parts(self.inner.translation, UnitQuaternion::identity()),
        }
    }

    fn rotation(&self) -> Iso3 {
        Self {
            inner: engeom::Iso3::from_parts(Translation3::identity(), self.inner.rotation),
        }
    }

    #[staticmethod]
    #[pyo3(signature=(e0, e1, origin=None))]
    fn from_basis_xy(e0: &Vector3, e1: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_xy(
            e0.get_inner(),
            e1.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    #[pyo3(signature=(e0, e2, origin=None))]
    fn from_basis_xz(e0: &Vector3, e2: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_xz(
            e0.get_inner(),
            e2.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    #[pyo3(signature=(e1, e2, origin=None))]
    fn from_basis_yz(e1: &Vector3, e2: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_yz(
            e1.get_inner(),
            e2.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    #[pyo3(signature=(e1, e0, origin=None))]
    fn from_basis_yx(e1: &Vector3, e0: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_yx(
            e1.get_inner(),
            e0.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    #[pyo3(signature=(e2, e0, origin=None))]
    fn from_basis_zx(e2: &Vector3, e0: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_zx(
            e2.get_inner(),
            e0.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    #[pyo3(signature=(e2, e1, origin=None))]
    fn from_basis_zy(e2: &Vector3, e1: &Vector3, origin: Option<Point3>) -> PyResult<Iso3> {
        let iso = engeom::Iso3::try_from_basis_zy(
            e2.get_inner(),
            e1.get_inner(),
            origin.map(|p| *p.get_inner()),
        )
        .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner: iso })
    }

    #[staticmethod]
    fn from_rx(angle: f64) -> Self {
        Self::from_inner(engeom::Iso3::from_rx(angle))
    }

    #[staticmethod]
    fn from_ry(angle: f64) -> Self {
        Self::from_inner(engeom::Iso3::from_ry(angle))
    }

    #[staticmethod]
    fn from_rz(angle: f64) -> Self {
        Self::from_inner(engeom::Iso3::from_rz(angle))
    }

    fn transform_points<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let points = array_to_points3(&points.as_array())?;
        let mut result = ArrayD::zeros(vec![points.len(), 3]);

        for (i, point) in points.iter().enumerate() {
            let transformed = self.inner * point;
            result[[i, 0]] = transformed.x;
            result[[i, 1]] = transformed.y;
            result[[i, 2]] = transformed.z;
        }

        Ok(result.into_pyarray(py))
    }

    fn transform_vectors<'py>(
        &self,
        py: Python<'py>,
        vectors: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let vectors = array_to_vectors3(&vectors.as_array())?;
        let mut result = ArrayD::zeros(vec![vectors.len(), 3]);

        for (i, vector) in vectors.iter().enumerate() {
            let transformed = self.inner * vector;
            result[[i, 0]] = transformed.x;
            result[[i, 1]] = transformed.y;
            result[[i, 2]] = transformed.z;
        }

        Ok(result.into_pyarray(py))
    }
}
