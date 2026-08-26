use crate::bounding::Aabb3;
use crate::common::Resample;
use crate::conversions::{
    array_to_points3, array_to_surface_points3, array_to_vectors3, dvec_from_array,
    dvec_to_array, points_to_array, vectors_to_array,
};
use crate::geom2::{Point2, SplineProjection, SurfacePoint2, Vector2};
use engeom::common::To2D;
use engeom::geom3::IsoExtensions3;
use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods,
};
use parry3d_f64::na::{Quaternion, Translation3, UnitQuaternion};
use pyo3::exceptions::PyIOError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::PyAnyMethods;
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyRef, PyResult, Python,
    pyclass, pyfunction, pymethods,
};
use std::path::PathBuf;

#[derive(FromPyObject)]
enum Vector3OrPoint3 {
    Vector(Vector3),
    Point(Point3),
}

// ================================================================================================
// Vectors
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
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

    fn __getnewargs__(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    fn __getstate__(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    fn __setstate__(&mut self, state: (f64, f64, f64)) {
        self.inner = engeom::Vector3::new(state.0, state.1, state.2);
    }

    fn __eq__(&self, other: &Vector3) -> bool {
        self.inner == other.inner
    }

    #[staticmethod]
    fn zero() -> Self {
        Self::from_inner(engeom::Vector3::zeros())
    }

    #[staticmethod]
    fn x_axis() -> Self {
        Self::from_inner(engeom::Vector3::x())
    }

    #[staticmethod]
    fn y_axis() -> Self {
        Self::from_inner(engeom::Vector3::y())
    }

    #[staticmethod]
    fn z_axis() -> Self {
        Self::from_inner(engeom::Vector3::z())
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

    fn angle(&self, other: Vector3) -> f64 {
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

    fn to_2d(&self) -> Vector2 {
        Vector2::from_inner(self.inner.to_2d())
    }
}

// ================================================================================================
// Points
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
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

    fn __getnewargs__(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    fn __getstate__(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    fn __setstate__(&mut self, state: (f64, f64, f64)) {
        self.inner = engeom::Point3::new(state.0, state.1, state.2);
    }

    fn __eq__(&self, other: &Point3) -> bool {
        self.inner == other.inner
    }

    #[staticmethod]
    fn origin() -> Self {
        Self::from_inner(engeom::Point3::origin())
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

    fn to_2d(&self) -> Point2 {
        Point2::from_inner(self.inner.to_2d())
    }
}

// ================================================================================================
// Surface Point
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom3")]
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

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.inner.point.x,
            self.inner.point.y,
            self.inner.point.z,
            self.inner.normal.x,
            self.inner.normal.y,
            self.inner.normal.z,
        )
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.inner.point.x,
            self.inner.point.y,
            self.inner.point.z,
            self.inner.normal.x,
            self.inner.normal.y,
            self.inner.normal.z,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64)) {
        let p = engeom::Point3::new(state.0, state.1, state.2);
        let v = engeom::Vector3::new(state.3, state.4, state.5);
        self.inner = engeom::SurfacePoint3::new_normalize(p, v);
    }

    fn __eq__(&self, other: &SurfacePoint3) -> bool {
        self.inner.point == other.inner.point && self.inner.normal == other.inner.normal
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

    fn transformed_by(&self, iso: Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
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
        Plane3::from_inner(engeom::Plane3::from_surface_point(&self.inner))
    }

    fn shifted(&self, offset: f64) -> Self {
        Self::from_inner(self.inner.shifted(offset))
    }

    fn to_2d(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(self.inner.to_2d())
    }

    fn slerp(&self, other: &Self, t: f64) -> Self {
        Self::from_inner(self.inner.slerp(&other.inner, t))
    }
}

// ================================================================================================
// Plane
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom3")]
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

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.normal.x,
            self.inner.normal.y,
            self.inner.normal.z,
            self.inner.d,
        )
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.normal.x,
            self.inner.normal.y,
            self.inner.normal.z,
            self.inner.d,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) -> PyResult<()> {
        let v = engeom::Vector3::new(state.0, state.1, state.2);
        let normal = engeom::UnitVec3::try_new(v, 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;

        self.inner = engeom::Plane3::new(normal, state.3);
        Ok(())
    }

    fn __eq__(&self, other: &Plane3) -> bool {
        self.inner.normal == other.inner.normal && self.inner.d == other.inner.d
    }

    #[staticmethod]
    fn xy() -> Self {
        Self::from_inner(engeom::Plane3::xy())
    }

    #[staticmethod]
    fn xz() -> Self {
        Self::from_inner(engeom::Plane3::xz())
    }

    #[staticmethod]
    fn yz() -> Self {
        Self::from_inner(engeom::Plane3::yz())
    }

    #[staticmethod]
    fn from_point_normal(px: f64, py: f64, pz: f64, nx: f64, ny: f64, nz: f64) -> PyResult<Self> {
        let point = engeom::Point3::new(px, py, pz);
        let normal = engeom::UnitVec3::try_new(engeom::Vector3::new(nx, ny, nz), 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;
        Ok(Self::from_inner(engeom::Plane3::from_point_normal(
            &point, &normal,
        )))
    }

    #[staticmethod]
    fn from_3_points(p1: Point3, p2: Point3, p3: Point3) -> PyResult<Self> {
        engeom::Plane3::from_3_points(p1.get_inner(), p2.get_inner(), p3.get_inner())
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn from_surface_point(surface_point: &SurfacePoint3) -> Self {
        Self::from_inner(engeom::Plane3::from_surface_point(
            surface_point.get_inner(),
        ))
    }

    #[staticmethod]
    #[pyo3(signature=(points, weights=None))]
    fn from_fit<'py>(
        points: PyReadonlyArray2<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let plane = match weights {
            Some(weights) => {
                engeom::Plane3::from_fit(&points, Some(weights.as_array().as_slice().unwrap()))
            }
            None => engeom::Plane3::from_fit(&points, None),
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(plane))
    }

    #[staticmethod]
    #[pyo3(signature=(points, sigma_max, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result = engeom::Plane3::from_consensus(&points, sigma_max, Some(options))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn normal_reversed(&self) -> Self {
        Self::from_inner(self.inner.normal_reversed())
    }

    fn signed_distance_to_point(&self, point: Point3) -> f64 {
        self.inner.signed_distance_to_point(point.get_inner())
    }

    fn distance_to_point(&self, point: Point3) -> f64 {
        self.inner.distance_to_point(point.get_inner())
    }

    fn point_is_positive(&self, point: Point3) -> bool {
        self.inner.point_is_positive(point.get_inner())
    }

    fn project_point(&self, point: Point3) -> Point3 {
        Point3::from_inner(self.inner.project_point(point.get_inner()))
    }

    fn project_vector(&self, v: Vector3) -> Vector3 {
        Vector3::from_inner(self.inner.project_vector(v.get_inner()))
    }

    fn offset_by(&self, shift: f64) -> Self {
        Self::from_inner(self.inner.offset_by(shift))
    }

    fn intersect_distance(&self, sp: &SurfacePoint3) -> Option<f64> {
        self.inner.intersect_distance(sp.get_inner())
    }

    #[getter]
    fn a(&self) -> f64 {
        self.inner.normal.x
    }

    #[getter]
    fn b(&self) -> f64 {
        self.inner.normal.y
    }

    #[getter]
    fn c(&self) -> f64 {
        self.inner.normal.z
    }

    #[getter]
    fn d(&self) -> f64 {
        self.inner.d
    }

    #[getter]
    fn normal(&self) -> Vector3 {
        Vector3::from_inner(self.inner.normal.into_inner())
    }

    fn intersect_plane(&self, other: &Plane3) -> Option<Line3> {
        self.inner
            .intersect_plane(&other.inner)
            .map(Line3::from_inner)
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }
}

// ================================================================================================
// Line3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Line3 {
    inner: engeom::Line3,
}

impl Line3 {
    pub fn get_inner(&self) -> &engeom::Line3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Line3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Line3 {
    #[new]
    fn new(ox: f64, oy: f64, oz: f64, dx: f64, dy: f64, dz: f64) -> Self {
        Self::from_inner(engeom::Line3::new(
            engeom::Point3::new(ox, oy, oz),
            engeom::Vector3::new(dx, dy, dz),
        ))
    }

    #[staticmethod]
    fn from_points(p1: Point3, p2: Point3) -> Self {
        Self::from_inner(engeom::Line3::from_points(p1.get_inner(), p2.get_inner()))
    }

    #[staticmethod]
    fn new_normalize(ox: f64, oy: f64, oz: f64, dx: f64, dy: f64, dz: f64) -> Self {
        Self::from_inner(engeom::Line3::new_normalize(
            engeom::Point3::new(ox, oy, oz),
            engeom::Vector3::new(dx, dy, dz),
        ))
    }

    #[staticmethod]
    #[pyo3(signature=(points, weights=None))]
    fn from_fit<'py>(
        points: PyReadonlyArray2<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let line = match weights {
            Some(weights) => {
                engeom::Line3::from_fit(&points, Some(weights.as_array().as_slice().unwrap()))
            }
            None => engeom::Line3::from_fit(&points, None),
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(line))
    }

    #[staticmethod]
    #[pyo3(signature=(points, sigma_max, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result = engeom::Line3::from_consensus(&points, sigma_max, Some(options))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    fn x_axis() -> Self {
        Self::from_inner(engeom::Line3::x_axis())
    }

    #[staticmethod]
    fn y_axis() -> Self {
        Self::from_inner(engeom::Line3::y_axis())
    }

    #[staticmethod]
    fn z_axis() -> Self {
        Self::from_inner(engeom::Line3::z_axis())
    }

    fn __repr__(&self) -> String {
        let o = self.inner.origin;
        let d = self.inner.direction;
        format!(
            "Line3(origin=({}, {}, {}), direction=({}, {}, {}))",
            o.x, o.y, o.z, d.x, d.y, d.z
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64) {
        let o = self.inner.origin;
        let d = self.inner.direction;
        (o.x, o.y, o.z, d.x, d.y, d.z)
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64) {
        let o = self.inner.origin;
        let d = self.inner.direction;
        (o.x, o.y, o.z, d.x, d.y, d.z)
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64)) {
        self.inner = engeom::Line3::new(
            engeom::Point3::new(state.0, state.1, state.2),
            engeom::Vector3::new(state.3, state.4, state.5),
        );
    }

    fn __eq__(&self, other: &Line3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn origin(&self) -> Point3 {
        Point3::from_inner(self.inner.origin)
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.direction)
    }

    fn at(&self, t: f64) -> Point3 {
        Point3::from_inner(self.inner.at(t))
    }

    fn scalar_project(&self, point: Point3) -> f64 {
        self.inner.scalar_project(point.get_inner())
    }

    fn closest_point(&self, point: Point3) -> Point3 {
        Point3::from_inner(self.inner.closest_point(point.get_inner()))
    }

    fn distance_to(&self, point: Point3) -> f64 {
        self.inner.distance_to(point.get_inner())
    }

    fn intersect_plane(&self, plane: &Plane3) -> Option<f64> {
        self.inner.intersect_plane(&plane.inner)
    }

    fn project_onto_plane(&self, plane: &Plane3) -> Option<Line3> {
        self.inner
            .project_onto_plane(&plane.inner)
            .map(Line3::from_inner)
    }

    fn normalized(&self) -> Line3 {
        Line3::from_inner(self.inner.normalized())
    }

    fn reversed(&self) -> Line3 {
        Line3::from_inner(self.inner.reversed())
    }

    fn shifted_origin(&self, delta_t: f64) -> Line3 {
        Line3::from_inner(self.inner.shifted_origin(delta_t))
    }

    fn slerp(&self, other: &Line3, t: f64) -> Line3 {
        Line3::from_inner(self.inner.slerp(&other.inner, t))
    }

    fn transformed_by(&self, iso: &Iso3) -> Line3 {
        Line3::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    fn intersect_sphere(&self, sphere: &Sphere3) -> Vec<f64> {
        self.inner.intersect_sphere(sphere.get_inner())
    }

    fn intersect_circle(&self, circle: &Circle3) -> Option<f64> {
        self.inner.intersect_circle(circle.get_inner())
    }
}

// ================================================================================================
// Segment3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Segment3 {
    inner: engeom::geom3::Segment3,
}

impl Segment3 {
    pub fn get_inner(&self) -> &engeom::geom3::Segment3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::Segment3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Segment3 {
    #[new]
    fn new(x0: f64, y0: f64, z0: f64, x1: f64, y1: f64, z1: f64) -> PyResult<Self> {
        let p0 = engeom::Point3::new(x0, y0, z0);
        let p1 = engeom::Point3::new(x1, y1, z1);
        Ok(Self {
            inner: engeom::geom3::Segment3::new(&p0, &p1)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        })
    }

    #[staticmethod]
    #[pyo3(signature=(points, weights=None))]
    fn from_fit<'py>(
        points: PyReadonlyArray2<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let result = match weights {
            Some(weights) => engeom::geom3::Segment3::from_fit(
                &points,
                Some(weights.as_array().as_slice().unwrap()),
            ),
            None => engeom::geom3::Segment3::from_fit(&points, None),
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    #[pyo3(signature=(points, sigma_max, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result = engeom::geom3::Segment3::from_consensus(&points, sigma_max, Some(options))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.inner.a.x,
            self.inner.a.y,
            self.inner.a.z,
            self.inner.b.x,
            self.inner.b.y,
            self.inner.b.z,
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.inner.a.x,
            self.inner.a.y,
            self.inner.a.z,
            self.inner.b.x,
            self.inner.b.y,
            self.inner.b.z,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64)) {
        let p0 = engeom::Point3::new(state.0, state.1, state.2);
        let p1 = engeom::Point3::new(state.3, state.4, state.5);
        self.inner =
            engeom::geom3::Segment3::new(&p0, &p1).expect("Invalid segment points in __setstate__");
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.a == other.inner.a && self.inner.b == other.inner.b
    }

    fn __repr__(&self) -> String {
        format!(
            "Segment3({}, {}, {}, {}, {}, {})",
            self.inner.a.x,
            self.inner.a.y,
            self.inner.a.z,
            self.inner.b.x,
            self.inner.b.y,
            self.inner.b.z,
        )
    }

    #[getter]
    fn a(&self) -> Point3 {
        Point3::from_inner(self.inner.a)
    }

    #[getter]
    fn b(&self) -> Point3 {
        Point3::from_inner(self.inner.b)
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.dir())
    }

    #[getter]
    fn length(&self) -> f64 {
        self.inner.length()
    }

    #[getter]
    fn aabb(&self) -> Aabb3 {
        Aabb3::from_inner(self.inner.aabb())
    }

    fn to_line(&self) -> Line3 {
        Line3::from_inner(self.inner.to_line())
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    fn at(&self, t: f64) -> Point3 {
        Point3::from_inner(self.inner.at(t))
    }

    fn reversed(&self) -> Self {
        Self::from_inner(self.inner.reversed())
    }

    fn scalar_projection(&self, other: Point3) -> f64 {
        self.inner.scalar_projection(other.get_inner())
    }

    fn closest_point(&self, other: Point3) -> Point3 {
        Point3::from_inner(self.inner.closest_point(other.get_inner()))
    }

    fn at_t(&self, t: f64) -> Manifold1Pos3 {
        Manifold1Pos3::from_inner(self.inner.at_t(t))
    }

    fn closest_to_point(&self, point: Point3) -> Manifold1Pos3 {
        Manifold1Pos3::from_inner(self.inner.closest_to_point(point.get_inner()))
    }
}

// ================================================================================================
// Sphere3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Sphere3 {
    inner: engeom::Sphere3,
}

impl Sphere3 {
    pub fn get_inner(&self) -> &engeom::Sphere3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Sphere3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Sphere3 {
    #[new]
    fn new(cx: f64, cy: f64, cz: f64, radius: f64) -> Self {
        Self::from_inner(engeom::Sphere3::new(
            &engeom::Point3::new(cx, cy, cz),
            radius,
        ))
    }

    #[staticmethod]
    fn from_4_points(p0: &Point3, p1: &Point3, p2: &Point3, p3: &Point3) -> PyResult<Self> {
        engeom::Sphere3::from_4_points(
            p0.get_inner(),
            p1.get_inner(),
            p2.get_inner(),
            p3.get_inner(),
        )
        .map(Self::from_inner)
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn from_min_enclosing<'py>(points: PyReadonlyArray2<'py, f64>) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        engeom::Sphere3::from_min_enclosing(&points)
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    #[pyo3(signature=(points, weights=None))]
    fn from_fit<'py>(
        points: PyReadonlyArray2<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let result = match weights {
            Some(weights) => {
                engeom::Sphere3::from_fit(&points, Some(weights.as_array().as_slice().unwrap()))
            }
            None => engeom::Sphere3::from_fit(&points, None),
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(points, sigma_max, min_r=None, max_r=None, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result =
            engeom::Sphere3::from_consensus(&points, sigma_max, min_r, max_r, Some(options))
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.center;
        format!(
            "Sphere3(center=({}, {}, {}), radius={})",
            c.x,
            c.y,
            c.z,
            self.inner.r()
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        let c = self.inner.center;
        (c.x, c.y, c.z, self.inner.r())
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        let c = self.inner.center;
        (c.x, c.y, c.z, self.inner.r())
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) {
        self.inner = engeom::Sphere3::new(&engeom::Point3::new(state.0, state.1, state.2), state.3);
    }

    fn __eq__(&self, other: &Sphere3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn center(&self) -> Point3 {
        Point3::from_inner(self.inner.center)
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    fn closest_point(&self, test_point: Point3) -> Option<SurfacePoint3> {
        self.inner
            .closest_point(test_point.get_inner())
            .map(SurfacePoint3::from_inner)
    }

    fn intersect_plane(&self, plane: &Plane3) -> Option<Circle3> {
        self.inner
            .intersect_plane(&plane.inner)
            .map(Circle3::from_inner)
    }

    fn intersect_sphere(&self, other: &Sphere3) -> Option<Circle3> {
        self.inner
            .intersect_sphere(&other.inner)
            .map(Circle3::from_inner)
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    fn intersect_ray(&self, line: &Line3) -> Option<SurfacePoint3> {
        self.inner
            .intersect_ray(line.get_inner())
            .map(SurfacePoint3::from_inner)
    }
}

// ================================================================================================
// Manifold1Pos3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Manifold1Pos3 {
    inner: engeom::geom3::Manifold1Pos3,
}

impl Manifold1Pos3 {
    pub fn from_inner(inner: engeom::geom3::Manifold1Pos3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Manifold1Pos3 {
    fn __repr__(&self) -> String {
        format!(
            "Manifold1Pos3(l={}, point=({}, {}, {}), direction=({}, {}, {}))",
            self.inner.l,
            self.inner.point.x,
            self.inner.point.y,
            self.inner.point.z,
            self.inner.direction.x,
            self.inner.direction.y,
            self.inner.direction.z,
        )
    }

    #[getter]
    fn l(&self) -> f64 {
        self.inner.l
    }

    #[getter]
    fn point(&self) -> Point3 {
        Point3::from_inner(self.inner.point)
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.direction.into_inner())
    }
}

// ================================================================================================
// Circle3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Circle3 {
    inner: engeom::geom3::Circle3,
}

impl Circle3 {
    pub fn get_inner(&self) -> &engeom::geom3::Circle3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::Circle3) -> Self {
        Self { inner }
    }
}

/// Build a `Magsac` consensus configuration from the optional overrides exposed to Python, leaving
/// any unset value at the library default.
fn magsac_options(
    sigma_max: f64,
    max_iterations: Option<usize>,
    refinement_steps: Option<usize>,
    confidence: Option<f64>,
    seed: Option<u64>,
) -> engeom::common::consensus::Magsac {
    let mut options = engeom::common::consensus::Magsac::new(sigma_max);
    options.max_iterations = max_iterations;
    if let Some(steps) = refinement_steps {
        options.refinement_steps = steps;
    }
    if let Some(confidence) = confidence {
        options.confidence = confidence;
    }
    options.seed = seed;
    options
}

#[pymethods]
impl Circle3 {
    #[new]
    fn new(cx: f64, cy: f64, cz: f64, nx: f64, ny: f64, nz: f64, radius: f64) -> PyResult<Self> {
        let center = engeom::Point3::new(cx, cy, cz);
        let normal = engeom::UnitVec3::try_new(engeom::Vector3::new(nx, ny, nz), 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;
        Ok(Circle3::from_inner(engeom::geom3::Circle3::new(
            center, normal, radius,
        )))
    }

    #[staticmethod]
    #[pyo3(signature=(points, weights=None))]
    fn from_fit<'py>(
        points: PyReadonlyArray2<'py, f64>,
        weights: Option<PyReadonlyArray1<'py, f64>>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let result = match weights {
            Some(weights) => engeom::geom3::Circle3::from_fit(
                &points,
                Some(weights.as_array().as_slice().unwrap()),
            ),
            None => engeom::geom3::Circle3::from_fit(&points, None),
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(points, sigma_max, min_r=None, max_r=None, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result =
            engeom::geom3::Circle3::from_consensus(&points, sigma_max, min_r, max_r, Some(options))
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(points, sigma_max, min_r=None, max_r=None, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus_planar<'py>(
        points: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let options = magsac_options(
            sigma_max,
            max_iterations,
            refinement_steps,
            confidence,
            seed,
        );
        let result = engeom::geom3::Circle3::from_consensus_planar(
            &points,
            sigma_max,
            min_r,
            max_r,
            Some(options),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.center;
        let n = self.inner.normal;
        format!(
            "Circle3(center=({}, {}, {}), normal=({}, {}, {}), radius={})",
            c.x,
            c.y,
            c.z,
            n.x,
            n.y,
            n.z,
            self.inner.r()
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        let c = self.inner.center;
        let n = self.inner.normal;
        (c.x, c.y, c.z, n.x, n.y, n.z, self.inner.r())
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        let c = self.inner.center;
        let n = self.inner.normal;
        (c.x, c.y, c.z, n.x, n.y, n.z, self.inner.r())
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64, f64)) -> PyResult<()> {
        let center = engeom::Point3::new(state.0, state.1, state.2);
        let normal =
            engeom::UnitVec3::try_new(engeom::Vector3::new(state.3, state.4, state.5), 1.0e-6)
                .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;
        self.inner = engeom::geom3::Circle3::new(center, normal, state.6);
        Ok(())
    }

    fn __eq__(&self, other: &Circle3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    #[getter]
    fn center(&self) -> Point3 {
        Point3::from_inner(self.inner.center)
    }

    #[getter]
    fn normal(&self) -> Vector3 {
        Vector3::from_inner(self.inner.normal.into_inner())
    }

    #[getter]
    fn plane(&self) -> Plane3 {
        Plane3::from_inner(self.inner.plane())
    }

    fn closest_point(&self, test_point: Point3) -> Option<SurfacePoint3> {
        self.inner
            .closest_point(test_point.get_inner())
            .map(SurfacePoint3::from_inner)
    }

    fn intersect_plane(&self, plane: &Plane3) -> Vec<Point3> {
        self.inner
            .intersect_plane(&plane.inner)
            .into_iter()
            .map(Point3::from_inner)
            .collect()
    }

    fn max_extent_point(&self, dx: f64, dy: f64, dz: f64) -> PyResult<Point3> {
        let dir = engeom::Vector3::new(dx, dy, dz);
        self.inner
            .max_extent_point(&dir)
            .map(Point3::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn flip_normal(&mut self) {
        self.inner = self.inner.normal_reversed();
    }

    fn normal_reversed(&self) -> Self {
        Circle3::from_inner(self.inner.normal_reversed())
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    #[staticmethod]
    fn from_3_points(p0: Point3, p1: Point3, p2: Point3) -> PyResult<Self> {
        engeom::geom3::Circle3::from_3_points(p0.get_inner(), p1.get_inner(), p2.get_inner())
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

// ================================================================================================
// Cylinder3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Cylinder3 {
    inner: engeom::geom3::Cylinder3,
}

impl Cylinder3 {
    pub fn get_inner(&self) -> &engeom::geom3::Cylinder3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::Cylinder3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Cylinder3 {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        cx: f64,
        cy: f64,
        cz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        radius: f64,
        length: f64,
    ) -> PyResult<Self> {
        let center = engeom::Point3::new(cx, cy, cz);
        let direction = engeom::UnitVec3::try_new(engeom::Vector3::new(dx, dy, dz), 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid direction vector"))?;
        Ok(Self::from_inner(engeom::geom3::Cylinder3::new(
            center, direction, radius, length,
        )))
    }

    #[staticmethod]
    fn from_points(p0: Point3, p1: Point3, radius: f64) -> PyResult<Self> {
        engeom::geom3::Cylinder3::from_points(p0.get_inner(), p1.get_inner(), radius)
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(points, normals, sigma_max, min_r=None, max_r=None, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        normals: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        min_r: Option<f64>,
        max_r: Option<f64>,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_surface_points3(&points.as_array(), &normals.as_array())?;
        let options = magsac_options(sigma_max, max_iterations, refinement_steps, confidence, seed);
        let result = engeom::geom3::Cylinder3::from_consensus(
            &points,
            sigma_max,
            min_r,
            max_r,
            Some(options),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.center;
        let d = self.inner.direction;
        format!(
            "Cylinder3(center=({}, {}, {}), direction=({}, {}, {}), radius={}, length={})",
            c.x,
            c.y,
            c.z,
            d.x,
            d.y,
            d.z,
            self.inner.r(),
            self.inner.length
        )
    }

    #[allow(clippy::type_complexity)]
    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let c = self.inner.center;
        let d = self.inner.direction;
        (
            c.x,
            c.y,
            c.z,
            d.x,
            d.y,
            d.z,
            self.inner.r(),
            self.inner.length,
        )
    }

    #[allow(clippy::type_complexity)]
    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let c = self.inner.center;
        let d = self.inner.direction;
        (
            c.x,
            c.y,
            c.z,
            d.x,
            d.y,
            d.z,
            self.inner.r(),
            self.inner.length,
        )
    }

    #[allow(clippy::type_complexity)]
    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64, f64, f64)) -> PyResult<()> {
        let center = engeom::Point3::new(state.0, state.1, state.2);
        let direction =
            engeom::UnitVec3::try_new(engeom::Vector3::new(state.3, state.4, state.5), 1.0e-6)
                .ok_or_else(|| PyValueError::new_err("Invalid direction vector"))?;
        self.inner = engeom::geom3::Cylinder3::new(center, direction, state.6, state.7);
        Ok(())
    }

    fn __eq__(&self, other: &Cylinder3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn center(&self) -> Point3 {
        Point3::from_inner(self.inner.center)
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.direction.into_inner())
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    #[getter]
    fn length(&self) -> f64 {
        self.inner.length
    }

    fn a(&self) -> Point3 {
        Point3::from_inner(self.inner.a())
    }

    fn b(&self) -> Point3 {
        Point3::from_inner(self.inner.b())
    }

    fn axis(&self) -> Line3 {
        Line3::from_inner(self.inner.axis())
    }

    fn start_cap(&self) -> Circle3 {
        Circle3::from_inner(self.inner.start_cap())
    }

    fn end_cap(&self) -> Circle3 {
        Circle3::from_inner(self.inner.end_cap())
    }

    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    fn lateral_area(&self) -> f64 {
        self.inner.lateral_area()
    }

    fn contains_point(&self, test_point: Point3) -> bool {
        self.inner.contains_point(test_point.get_inner())
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    fn reversed(&self) -> Self {
        Self::from_inner(self.inner.reversed())
    }

    fn closest_point(&self, test_point: Point3, infinite: bool) -> Option<SurfacePoint3> {
        self.inner
            .closest_point(test_point.get_inner(), infinite)
            .map(SurfacePoint3::from_inner)
    }
}

// ================================================================================================
// Cone3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct Cone3 {
    inner: engeom::geom3::Cone3,
}

impl Cone3 {
    pub fn get_inner(&self) -> &engeom::geom3::Cone3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::Cone3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Cone3 {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        tx: f64,
        ty: f64,
        tz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        height: f64,
        radius: f64,
    ) -> PyResult<Self> {
        let tip = engeom::Point3::new(tx, ty, tz);
        let direction = engeom::UnitVec3::try_new(engeom::Vector3::new(dx, dy, dz), 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid direction vector"))?;
        Ok(Self::from_inner(engeom::geom3::Cone3::new(
            tip, direction, height, radius,
        )))
    }

    #[staticmethod]
    fn from_points(tip: Point3, base_center: Point3, radius: f64) -> PyResult<Self> {
        engeom::geom3::Cone3::from_points(tip.get_inner(), base_center.get_inner(), radius)
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(points, normals, sigma_max, min_half_angle=None, max_half_angle=None, max_iterations=None, refinement_steps=None, confidence=None, seed=None))]
    fn from_consensus<'py>(
        points: PyReadonlyArray2<'py, f64>,
        normals: PyReadonlyArray2<'py, f64>,
        sigma_max: f64,
        min_half_angle: Option<f64>,
        max_half_angle: Option<f64>,
        max_iterations: Option<usize>,
        refinement_steps: Option<usize>,
        confidence: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        let points = array_to_surface_points3(&points.as_array(), &normals.as_array())?;
        let options = magsac_options(sigma_max, max_iterations, refinement_steps, confidence, seed);
        let result = engeom::geom3::Cone3::from_consensus(
            &points,
            sigma_max,
            min_half_angle,
            max_half_angle,
            Some(options),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    fn __repr__(&self) -> String {
        let t = self.inner.tip;
        let d = self.inner.direction;
        format!(
            "Cone3(tip=({}, {}, {}), direction=({}, {}, {}), height={}, radius={})",
            t.x,
            t.y,
            t.z,
            d.x,
            d.y,
            d.z,
            self.inner.height,
            self.inner.r()
        )
    }

    #[allow(clippy::type_complexity)]
    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let t = self.inner.tip;
        let d = self.inner.direction;
        (
            t.x,
            t.y,
            t.z,
            d.x,
            d.y,
            d.z,
            self.inner.height,
            self.inner.r(),
        )
    }

    #[allow(clippy::type_complexity)]
    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let t = self.inner.tip;
        let d = self.inner.direction;
        (
            t.x,
            t.y,
            t.z,
            d.x,
            d.y,
            d.z,
            self.inner.height,
            self.inner.r(),
        )
    }

    #[allow(clippy::type_complexity)]
    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64, f64, f64)) -> PyResult<()> {
        let tip = engeom::Point3::new(state.0, state.1, state.2);
        let direction =
            engeom::UnitVec3::try_new(engeom::Vector3::new(state.3, state.4, state.5), 1.0e-6)
                .ok_or_else(|| PyValueError::new_err("Invalid direction vector"))?;
        self.inner = engeom::geom3::Cone3::new(tip, direction, state.6, state.7);
        Ok(())
    }

    fn __eq__(&self, other: &Cone3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn tip(&self) -> Point3 {
        Point3::from_inner(self.inner.tip)
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.direction.into_inner())
    }

    #[getter]
    fn height(&self) -> f64 {
        self.inner.height
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    fn base_center(&self) -> Point3 {
        Point3::from_inner(self.inner.base_center())
    }

    fn axis(&self) -> Line3 {
        Line3::from_inner(self.inner.axis())
    }

    fn base(&self) -> Circle3 {
        Circle3::from_inner(self.inner.base())
    }

    fn half_angle(&self) -> f64 {
        self.inner.half_angle()
    }

    fn slant_height(&self) -> f64 {
        self.inner.slant_height()
    }

    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    fn lateral_area(&self) -> f64 {
        self.inner.lateral_area()
    }

    fn contains_point(&self, test_point: Point3) -> bool {
        self.inner.contains_point(test_point.get_inner())
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }

    fn closest_point(&self, test_point: Point3, infinite: bool) -> Option<SurfacePoint3> {
        self.inner
            .closest_point(test_point.get_inner(), infinite)
            .map(SurfacePoint3::from_inner)
    }
}

// ================================================================================================
// Curve
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom3")]
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

#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct Curve3 {
    inner: engeom::Curve3,
    points: Option<Py<PyArray2<f64>>>,
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
    fn new(points: PyReadonlyArray2<'_, f64>, tol: f64) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let inner = engeom::Curve3::from_points(&points, tol)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let result = points_to_array(self.inner.vertices());
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

    fn new_transformed_by(&self, iso: Iso3) -> Self {
        Self::from_inner(self.inner.new_transformed_by(iso.get_inner()))
    }

    fn to_2d(&self) -> PyResult<crate::geom2::Curve2> {
        let c = self
            .inner
            .to_2d()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(crate::geom2::Curve2::from_inner(c))
    }

    #[staticmethod]
    fn load_tccurve3(path: PathBuf) -> PyResult<Self> {
        let curve =
            engeom::Curve3::load_tccurve3(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(curve))
    }

    fn save_tccurve3(&self, path: PathBuf, tol: f64) -> PyResult<()> {
        self.inner
            .save_tccurve3(&path, tol)
            .map_err(|e| PyIOError::new_err(e.to_string()))
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

#[derive(FromPyObject)]
enum Transformable3 {
    Iso(Iso3),
    Vec(Vector3),
    Pnt(Point3),
    Plane(Plane3),
    Sp(SurfacePoint3),
    Line(Line3),
    Sphere(Sphere3),
    Circle(Circle3),
    Seg(Segment3),
}

#[pyclass(from_py_object, module = "engeom.geom3")]
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

    fn __eq__(&self, other: &Iso3) -> PyResult<bool> {
        Ok(self.inner == other.inner)
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
    fn from_quaternion(tx: f64, ty: f64, tz: f64, i: f64, j: f64, k: f64, w: f64) -> Self {
        let translation = Translation3::new(tx, ty, tz);
        let quat = Quaternion::new(w, i, j, k);
        let rotation = UnitQuaternion::from_quaternion(quat);

        Self {
            inner: engeom::Iso3::from_parts(translation, rotation),
        }
    }

    #[staticmethod]
    fn from_xyzwpr(x: f64, y: f64, z: f64, w: f64, p: f64, r: f64) -> Self {
        Self {
            inner: (&engeom::geom3::XyzWpr::new(x, y, z, w, p, r)).into(),
        }
    }

    fn to_quaternion(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        let t = &self.inner.translation;
        let r = &self.inner.rotation.quaternion();
        (t.x, t.y, t.z, r.i, r.j, r.k, r.w)
    }

    fn to_xyzwpr(&self) -> Vec<f64> {
        let v = engeom::geom3::XyzWpr::from(&self.inner);
        vec![v.x, v.y, v.z, v.w, v.p, v.r]
    }

    #[new]
    fn new(matrix: PyReadonlyArray2<'_, f64>) -> PyResult<Self> {
        if matrix.shape().len() != 2 || matrix.shape()[0] != 4 || matrix.shape()[1] != 4 {
            return Err(PyValueError::new_err("Expected 4x4 matrix"));
        }

        let mut array = [0.0; 16];
        for (i, value) in matrix.as_array().iter().enumerate() {
            array[i] = *value;
        }

        let inner = engeom::Iso3::from_array(&array)
            .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner })
    }

    fn __getnewargs__<'py>(&self, py: Python<'py>) -> (Bound<'py, PyArray2<f64>>,) {
        (self.as_numpy(py),)
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        self.to_quaternion()
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64, f64)) {
        let translation = Translation3::new(state.0, state.1, state.2);
        let quat = Quaternion::new(state.6, state.3, state.4, state.5);
        self.inner = engeom::Iso3::from_parts(translation, UnitQuaternion::from_quaternion(quat));
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

    #[staticmethod]
    fn from_rot_axis(axis: &Line3, angle: f64) -> PyResult<Self> {
        let inner = engeom::Iso3::from_rot_axis(axis.get_inner(), angle)
            .map_err(|e| PyValueError::new_err(format!("Error creating Iso3: {}", e)))?;

        Ok(Self { inner })
    }

    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse(),
        }
    }

    #[getter]
    fn x_direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.rotation * engeom::Vector3::x())
    }

    #[getter]
    fn y_direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.rotation * engeom::Vector3::y())
    }

    #[getter]
    fn z_direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.rotation * engeom::Vector3::z())
    }

    #[getter]
    fn x_axis(&self) -> Line3 {
        let origin = self.inner * engeom::Point3::origin();
        let direction = self.inner.rotation * engeom::Vector3::x();
        Line3::from_inner(engeom::Line3::new_normalize(origin, direction))
    }

    #[getter]
    fn y_axis(&self) -> Line3 {
        let origin = self.inner * engeom::Point3::origin();
        let direction = self.inner.rotation * engeom::Vector3::y();
        Line3::from_inner(engeom::Line3::new_normalize(origin, direction))
    }

    #[getter]
    fn z_axis(&self) -> Line3 {
        let origin = self.inner * engeom::Point3::origin();
        let direction = self.inner.rotation * engeom::Vector3::z();
        Line3::from_inner(engeom::Line3::new_normalize(origin, direction))
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
                Plane3::from_inner(other.inner.transformed_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Sp(other) => {
                SurfacePoint3::from_inner(other.inner.transformed_by(&self.inner))
                    .into_bound_py_any(py)
            }
            Transformable3::Line(other) => {
                Line3::from_inner(other.inner.transformed_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Sphere(other) => {
                Sphere3::from_inner(other.inner.transformed_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Circle(other) => {
                Circle3::from_inner(other.inner.transformed_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Seg(other) => {
                Segment3::from_inner(other.inner.transformed_by(&self.inner)).into_bound_py_any(py)
            }
        }
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut result = Array2::zeros((4, 4));
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

    fn flipped_around_x(&self) -> Self {
        Self {
            inner: self.inner.flipped_around_x(),
        }
    }

    fn flipped_around_y(&self) -> Self {
        Self {
            inner: self.inner.flipped_around_y(),
        }
    }

    fn flipped_around_z(&self) -> Self {
        Self {
            inner: self.inner.flipped_around_z(),
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
        let iso = engeom::Iso3::from_basis_xy(
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
        let iso = engeom::Iso3::from_basis_xz(
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
        let iso = engeom::Iso3::from_basis_yz(
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
        let iso = engeom::Iso3::from_basis_yx(
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
        let iso = engeom::Iso3::from_basis_zx(
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
        let iso = engeom::Iso3::from_basis_zy(
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
        points: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let points = array_to_points3(&points.as_array())?;
        let transformed = points
            .iter()
            .map(|point| self.inner * point)
            .collect::<Vec<_>>();
        let result = points_to_array(&transformed);
        Ok(result.into_pyarray(py))
    }

    fn transform_vectors<'py>(
        &self,
        py: Python<'py>,
        vectors: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let vectors = array_to_vectors3(&vectors.as_array())?;
        let transformed: Vec<engeom::Vector3> =
            vectors.iter().map(|vector| self.inner * vector).collect();
        let result = vectors_to_array(&transformed);
        Ok(result.into_pyarray(py))
    }
}

// ================================================================================================
// CubicSpline3
// ================================================================================================

/// Filters a fixed-size `[f64; 2]` root array (where unused slots are `NaN`) down to a `Vec`
/// containing only the finite values, for a cleaner Python-side result than a NaN-padded tuple.
fn finite_roots(roots: [f64; 2]) -> Vec<f64> {
    roots.into_iter().filter(|v| !v.is_nan()).collect()
}

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone, Debug)]
pub struct CubicSpline3 {
    inner: engeom::geom3::CubicSpline3,
}

impl CubicSpline3 {
    pub fn get_inner(&self) -> &engeom::geom3::CubicSpline3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::CubicSpline3) -> Self {
        Self { inner }
    }
}

type CubicSpline3State = (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);

#[pymethods]
impl CubicSpline3 {
    #[new]
    // Four control points as flat coordinates; the signature mirrors the Python API.
    #[allow(clippy::too_many_arguments)]
    fn new(
        x0: f64,
        y0: f64,
        z0: f64,
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        x3: f64,
        y3: f64,
        z3: f64,
    ) -> Self {
        Self {
            inner: engeom::geom3::CubicSpline3::new(
                engeom::Point3::new(x0, y0, z0),
                engeom::Point3::new(x1, y1, z1),
                engeom::Point3::new(x2, y2, z2),
                engeom::Point3::new(x3, y3, z3),
            ),
        }
    }

    fn __getstate__(&self) -> CubicSpline3State {
        (
            self.inner.p0.x,
            self.inner.p0.y,
            self.inner.p0.z,
            self.inner.p1.x,
            self.inner.p1.y,
            self.inner.p1.z,
            self.inner.p2.x,
            self.inner.p2.y,
            self.inner.p2.z,
            self.inner.p3.x,
            self.inner.p3.y,
            self.inner.p3.z,
        )
    }

    fn __getnewargs__(&self) -> CubicSpline3State {
        self.__getstate__()
    }

    fn __setstate__(&mut self, state: CubicSpline3State) {
        self.inner = engeom::geom3::CubicSpline3::new(
            engeom::Point3::new(state.0, state.1, state.2),
            engeom::Point3::new(state.3, state.4, state.5),
            engeom::Point3::new(state.6, state.7, state.8),
            engeom::Point3::new(state.9, state.10, state.11),
        );
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!(
            "CubicSpline3({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
            self.inner.p0.x,
            self.inner.p0.y,
            self.inner.p0.z,
            self.inner.p1.x,
            self.inner.p1.y,
            self.inner.p1.z,
            self.inner.p2.x,
            self.inner.p2.y,
            self.inner.p2.z,
            self.inner.p3.x,
            self.inner.p3.y,
            self.inner.p3.z,
        )
    }

    #[getter]
    fn p0(&self) -> Point3 {
        Point3::from_inner(self.inner.p0)
    }

    #[getter]
    fn p1(&self) -> Point3 {
        Point3::from_inner(self.inner.p1)
    }

    #[getter]
    fn p2(&self) -> Point3 {
        Point3::from_inner(self.inner.p2)
    }

    #[getter]
    fn p3(&self) -> Point3 {
        Point3::from_inner(self.inner.p3)
    }

    fn position(&self, t: f64) -> Point3 {
        Point3::from_inner(self.inner.position(t))
    }

    fn derivative(&self, t: f64) -> Vector3 {
        Vector3::from_inner(self.inner.derivative(t))
    }

    fn tangent(&self, t: f64) -> Vector3 {
        Vector3::from_inner(self.inner.tangent(t).into_inner())
    }

    fn curvature(&self, t: f64) -> f64 {
        self.inner.curvature(t)
    }

    /// Find the point of maximum curvature on the curve over the parameter range `[0, 1]`,
    /// returned as a `(t, curvature)` tuple: the parameter `t` at which the maximum occurs and the
    /// curvature magnitude there. For a fully degenerate curve whose curvature is undefined
    /// everywhere, `t` is `0.0` and the curvature is NaN.
    fn find_max_curvature(&self) -> (f64, f64) {
        let result = self.inner.find_max_curvature();
        (result.t, result.value)
    }

    fn second_derivative(&self, t: f64) -> Vector3 {
        Vector3::from_inner(self.inner.second_derivative(t))
    }

    fn polyline<'py>(&self, py: Python<'py>, tolerance: f64) -> Bound<'py, PyArray2<f64>> {
        let points = self.inner.polyline(tolerance);
        points_to_array(&points).into_pyarray(py)
    }

    /// Build a reusable acceleration structure for closest-point queries against this spline.
    fn query(&self) -> CubicSplineQueries3 {
        CubicSplineQueries3::from_inner(engeom::common::cubic_spline::CubicSplineQueries::from(
            &self.inner,
        ))
    }

    /// Find the closest point on the spline to the given point. This builds a temporary query
    /// structure each call; for repeated queries against the same spline build a `query()` once
    /// and reuse it.
    fn project_point(&self, point: Point3) -> SplineProjection {
        let queries = engeom::common::cubic_spline::CubicSplineQueries::from(&self.inner);
        SplineProjection::from_inner(queries.project_point(point.get_inner()))
    }

    /// Return the position and derivative direction of the curve at parameter `t` in the form of
    /// a parameterized line.
    fn line_at(&self, t: f64) -> Line3 {
        Line3::from_inner(self.inner.line_at(t))
    }

    /// Returns the real roots of the derivative of each component of the curve, as a tuple of
    /// `(x_roots, y_roots, z_roots)` where each is a list of 0, 1, or 2 parameter values.
    fn derivative_roots(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let roots = self.inner.derivative_roots();
        (
            finite_roots(roots[0]),
            finite_roots(roots[1]),
            finite_roots(roots[2]),
        )
    }

    /// Returns the parameter `t` of a cusp if one exists in `[0, 1]`, otherwise `None`.
    fn find_cusp(&self) -> Option<f64> {
        self.inner.find_cusp()
    }

    /// Returns parameter values in `[0, 1]` where the curve's curvature is zero, as a list of 0,
    /// 1, or 2 parameter values.
    fn find_curvature_zeros(&self) -> Vec<f64> {
        finite_roots(self.inner.find_curvature_zeros())
    }

    /// Returns every local maximum of the curvature over `[0, 1]`, each as a `(t, curvature)`
    /// tuple, ordered by ascending `t`.
    fn find_curvature_maxima(&self) -> Vec<(f64, f64)> {
        self.inner
            .find_curvature_maxima()
            .into_iter()
            .map(|v| (v.t, v.value))
            .collect()
    }

    /// Returns the corners `(min, max)` of the tight axis-aligned bounding box of the curve over
    /// the parameter range `[0, 1]`.
    fn compute_bounds(&self) -> (Point3, Point3) {
        let (lo, hi) = self.inner.compute_bounds();
        (Point3::from_inner(lo), Point3::from_inner(hi))
    }

    /// Returns the arc length of the curve over the parameter range `[t0, t1]`.
    fn arc_length_between(&self, t0: f64, t1: f64) -> f64 {
        self.inner.arc_length_between(t0, t1)
    }

    /// Returns the total arc length of the curve over the parameter range `[0, 1]`.
    fn arc_length(&self) -> f64 {
        self.inner.arc_length()
    }

    /// Splits the curve at parameter `t` using de Casteljau's algorithm, returning the left and
    /// right sub-curves.
    fn split(&self, t: f64) -> (Self, Self) {
        let (left, right) = self.inner.split(t);
        (Self::from_inner(left), Self::from_inner(right))
    }

    /// Splits the curve at parameter `t`, returning `None` if `t` is not in `[0, 1]`.
    fn try_split(&self, t: f64) -> Option<(Self, Self)> {
        self.inner
            .try_split(t)
            .map(|(left, right)| (Self::from_inner(left), Self::from_inner(right)))
    }

    /// Returns the axis-aligned bounding box of the curve, computed on demand.
    #[getter]
    fn aabb(&self) -> Aabb3 {
        Aabb3::from_inner(self.inner.aabb())
    }

    fn transformed_by(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
    }
}

// ================================================================================================
// CubicSplineQueries3
// ================================================================================================

/// A prebuilt acceleration structure for running repeated closest-point queries against a
/// `CubicSpline3`. Build it once with `spline.query()` (or the constructor) and reuse it across
/// many queries.
#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct CubicSplineQueries3 {
    inner: engeom::common::cubic_spline::CubicSplineQueries<3>,
}

impl CubicSplineQueries3 {
    pub fn from_inner(inner: engeom::common::cubic_spline::CubicSplineQueries<3>) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl CubicSplineQueries3 {
    #[new]
    fn new(spline: &CubicSpline3) -> Self {
        Self::from_inner(engeom::common::cubic_spline::CubicSplineQueries::from(
            spline.get_inner(),
        ))
    }

    /// Find the closest point on the spline to the given point, returned as a `SplineProjection`.
    fn project_point(&self, point: Point3) -> SplineProjection {
        SplineProjection::from_inner(self.inner.project_point(point.get_inner()))
    }

    /// Project an `Nx3` array of points onto the spline, returning an `Nx2` array whose columns are
    /// the closest-point parameter `t` and the distance to the curve, one row per input point.
    fn project_points<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let pts = array_to_points3(&points.as_array())?;
        let mut out = Array2::zeros((pts.len(), 2));
        for (i, p) in pts.iter().enumerate() {
            let proj = self.inner.project_point(p);
            out[[i, 0]] = proj.t;
            out[[i, 1]] = proj.value;
        }
        Ok(out.into_pyarray(py))
    }
}

// ================================================================================================
// Cubic spline fitting
// ================================================================================================

fn make_spline_builder3(py_builder: Py<PyAny>) -> engeom::common::cubic_spline::SplineBuildFn<3> {
    Box::new(move |params: &engeom::DVector| {
        let arr = Array1::from_iter(params.iter().copied());

        // `Python::attach` is re-entrant in PyO3 0.28 and safe to call here: the LM minimiser runs
        // synchronously inside a #[pyfunction] that already holds the GIL.
        Python::attach(|py| {
            let np_arr = arr.into_pyarray(py);

            let result = py_builder
                .bind(py)
                .call1((np_arr,))
                .map_err(|e| e.to_string())?;

            let spline = result
                .extract::<PyRef<CubicSpline3>>()
                .map_err(|_| "builder must return a CubicSpline3".to_string())?;

            Ok(spline.get_inner().clone())
        })
    })
}

/// Fit a `CubicSpline3` to a set of points using a Levenberg-Marquardt minimization and a
/// user-supplied builder function that maps a parameter vector to a `CubicSpline3`. Returns the
/// optimized parameter vector.
#[pyfunction]
#[pyo3(signature = (points, builder, initial))]
pub fn fit_spline_to_points<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    builder: Py<PyAny>,
    initial: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let pts = array_to_points3(&points.as_array())?;
    let initial_vec = dvec_from_array(&initial)?;
    let bld = make_spline_builder3(builder);

    let result = engeom::common::cubic_spline::fit_spline_to_points(&pts, &bld, initial_vec)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok(dvec_to_array(py, result.params))
}
