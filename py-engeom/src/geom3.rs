use crate::common::Resample;
use crate::conversions::{array_to_points3, array_to_vectors3, points_to_array, vectors_to_array};
use crate::geom2::{Point2, SurfacePoint2, Vector2};
use engeom::common::To2D;
use engeom::geom3::IsoExtensions3;
use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use parry3d_f64::na::{Quaternion, Translation3, UnitQuaternion};
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

    fn new_shifted(&self, offset: f64) -> Self {
        Self::from_inner(self.inner.new_shifted(offset))
    }

    fn to_2d(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(self.inner.to_2d())
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

    fn inverted_normal(&self) -> Self {
        Self::from_inner(self.inner.inverted_normal())
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

    fn shifted(&self, shift: f64) -> Self {
        Self::from_inner(self.inner.shifted(shift))
    }

    fn intersection_distance(&self, sp: &SurfacePoint3) -> Option<f64> {
        self.inner.intersection_distance(sp.get_inner())
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
        Self::from_inner(engeom::Line3::from_points(*p1.get_inner(), *p2.get_inner()))
    }

    fn __repr__(&self) -> String {
        let o = self.inner.origin();
        let d = self.inner.direction();
        format!(
            "Line3(origin=({}, {}, {}), direction=({}, {}, {}))",
            o.x, o.y, o.z, d.x, d.y, d.z
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64, f64) {
        let o = self.inner.origin();
        let d = self.inner.direction();
        (o.x, o.y, o.z, d.x, d.y, d.z)
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64) {
        let o = self.inner.origin();
        let d = self.inner.direction();
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
        Point3::from_inner(self.inner.origin())
    }

    #[getter]
    fn direction(&self) -> Vector3 {
        Vector3::from_inner(self.inner.direction())
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
            engeom::Point3::new(cx, cy, cz),
            radius,
        ))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.center();
        format!(
            "Sphere3(center=({}, {}, {}), radius={})",
            c.x,
            c.y,
            c.z,
            self.inner.r()
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        let c = self.inner.center();
        (c.x, c.y, c.z, self.inner.r())
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        let c = self.inner.center();
        (c.x, c.y, c.z, self.inner.r())
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) {
        self.inner = engeom::Sphere3::new(engeom::Point3::new(state.0, state.1, state.2), state.3);
    }

    fn __eq__(&self, other: &Sphere3) -> bool {
        self.inner == other.inner
    }

    #[getter]
    fn center(&self) -> Point3 {
        Point3::from_inner(self.inner.center())
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
}

// ================================================================================================
// Manifold1Pos3
// ================================================================================================

#[pyclass(module = "engeom.geom3")]
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

#[pymethods]
impl Circle3 {
    #[new]
    fn new(cx: f64, cy: f64, cz: f64, nx: f64, ny: f64, nz: f64, radius: f64) -> PyResult<Self> {
        let center = engeom::Point3::new(cx, cy, cz);
        let normal = engeom::UnitVec3::try_new(engeom::Vector3::new(nx, ny, nz), 1.0e-6)
            .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;
        engeom::geom3::Circle3::from_point_normal(&center, &normal, radius)
            .map(Circle3::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.center();
        let n = self.inner.normal();
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
        let c = self.inner.center();
        let n = self.inner.normal();
        (c.x, c.y, c.z, n.x, n.y, n.z, self.inner.r())
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        let c = self.inner.center();
        let n = self.inner.normal();
        (c.x, c.y, c.z, n.x, n.y, n.z, self.inner.r())
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64, f64, f64)) -> PyResult<()> {
        let center = engeom::Point3::new(state.0, state.1, state.2);
        let normal =
            engeom::UnitVec3::try_new(engeom::Vector3::new(state.3, state.4, state.5), 1.0e-6)
                .ok_or_else(|| PyValueError::new_err("Invalid normal vector"))?;
        self.inner = engeom::geom3::Circle3::from_point_normal(&center, &normal, state.6)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
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
        Point3::from_inner(self.inner.center())
    }

    #[getter]
    fn normal(&self) -> Vector3 {
        Vector3::from_inner(self.inner.normal().into_inner())
    }

    #[getter]
    fn plane(&self) -> Plane3 {
        Plane3::from_inner(self.inner.plane())
    }

    #[getter]
    fn iso(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.iso())
    }

    fn at_angle(&self, angle: f64) -> Manifold1Pos3 {
        Manifold1Pos3::from_inner(self.inner.at_angle(angle))
    }

    fn closest_angle(&self, test_point: Point3) -> f64 {
        self.inner.closest_angle(test_point.get_inner())
    }

    fn closest_position(&self, test_point: Point3) -> Manifold1Pos3 {
        Manifold1Pos3::from_inner(self.inner.closest_position(test_point.get_inner()))
    }

    fn intersect_plane(&self, plane: &Plane3) -> Vec<f64> {
        self.inner.intersect_plane(&plane.inner)
    }

    fn max_extent_angle(&self, dx: f64, dy: f64, dz: f64) -> PyResult<f64> {
        let dir = engeom::Vector3::new(dx, dy, dz);
        self.inner
            .max_extent_angle(&dir)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn set_zero_angle(&mut self, angle: f64) {
        self.inner.set_zero_angle(angle);
    }

    fn at_angles<'py>(
        &self,
        py: Python<'py>,
        angles: PyReadonlyArray1<'_, f64>,
    ) -> Bound<'py, PyArray2<f64>> {
        let angles = angles.as_array();
        let n = angles.len();
        let normal = self.inner.normal().into_inner();
        let mut result = Array2::zeros((n, 9));
        for (i, &angle) in angles.iter().enumerate() {
            let m = self.inner.at_angle(angle);
            result[[i, 0]] = m.point.x;
            result[[i, 1]] = m.point.y;
            result[[i, 2]] = m.point.z;
            result[[i, 3]] = normal.x;
            result[[i, 4]] = normal.y;
            result[[i, 5]] = normal.z;
            result[[i, 6]] = m.direction.x;
            result[[i, 7]] = m.direction.y;
            result[[i, 8]] = m.direction.z;
        }
        result.into_pyarray(py)
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
                Plane3::from_inner(other.inner.new_transformed_by(&self.inner))
                    .into_bound_py_any(py)
            }
            Transformable3::Sp(other) => {
                SurfacePoint3::from_inner(other.inner.transformed(&self.inner))
                    .into_bound_py_any(py)
            }
            Transformable3::Line(other) => {
                Line3::from_inner(other.inner.new_transformed_by(&self.inner)).into_bound_py_any(py)
            }
            Transformable3::Sphere(other) => {
                Sphere3::from_inner(other.inner.new_transformed_by(&self.inner))
                    .into_bound_py_any(py)
            }
            Transformable3::Circle(other) => {
                Circle3::from_inner(other.inner.new_transformed_by(&self.inner))
                    .into_bound_py_any(py)
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
