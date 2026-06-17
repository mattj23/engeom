use crate::bounding::Aabb2;
use crate::common::{AngleDir, Resample};
use crate::conversions::{array_to_points2, array_to_vectors2, points_to_array, vectors_to_array};
use engeom::geom2::HasBounds2;
use engeom::{BestFit, To3D};
use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyResult, Python, pyclass,
    pyfunction, pymethods,
};
use std::path::PathBuf;

#[derive(FromPyObject)]
enum Vector2OrPoint2 {
    Vector(Vector2),
    Point(Point2),
}

/// Returns an isometry representing a 90-degree rotation in the given direction.
#[pyfunction]
pub fn rot90(dir: AngleDir) -> Iso2 {
    Iso2::from_inner(engeom::geom2::rot90(dir.into()))
}

/// Returns an isometry representing a 270-degree rotation in the given direction.
#[pyfunction]
pub fn rot270(dir: AngleDir) -> Iso2 {
    Iso2::from_inner(engeom::geom2::rot270(dir.into()))
}

/// Returns the signed angle from `v1` to `v2` in radians, in the range (-π, π].
/// Positive means the shortest rotation is counter-clockwise; negative means clockwise.
#[pyfunction]
pub fn signed_angle(v1: &Vector2, v2: &Vector2) -> f64 {
    engeom::geom2::signed_angle(v1.get_inner(), v2.get_inner())
}

/// Returns the angle from `v1` to `v2` measured in the given direction, in radians, in [0, 2π].
#[pyfunction]
pub fn directed_angle(v1: &Vector2, v2: &Vector2, direction: AngleDir) -> f64 {
    engeom::geom2::directed_angle(v1.get_inner(), v2.get_inner(), direction.into())
}

// ================================================================================================
// Vectors
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Vector2 {
    inner: engeom::Vector2,
}

impl Vector2 {
    pub fn get_inner(&self) -> &engeom::Vector2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Vector2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Vector2 {
    #[new]
    fn new(x: f64, y: f64) -> Self {
        Self {
            inner: engeom::Vector2::new(x, y),
        }
    }

    fn __getstate__(&self) -> (f64, f64) {
        (self.inner.x, self.inner.y)
    }

    fn __getnewargs__(&self) -> (f64, f64) {
        (self.inner.x, self.inner.y)
    }

    fn __setstate__(&mut self, state: (f64, f64)) {
        self.inner = engeom::Vector2::new(state.0, state.1);
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    #[staticmethod]
    fn zero() -> Self {
        Self::from_inner(engeom::Vector2::zeros())
    }

    #[staticmethod]
    fn x_axis() -> Self {
        Self::from_inner(engeom::Vector2::x())
    }

    #[staticmethod]
    fn y_axis() -> Self {
        Self::from_inner(engeom::Vector2::y())
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let o = [self.inner.x, self.inner.y];
        PyIterator::from_object(&o.into_pyobject(py)?)
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut array = Array1::zeros(2);
        array[0] = self.inner.x;
        array[1] = self.inner.y;
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

    fn __rmul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self {
            inner: self.inner / other,
        }
    }

    fn __add__<'py>(&self, py: Python<'py>, other: Vector2OrPoint2) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Vector2OrPoint2::Vector(other) => {
                let result = self.inner + other.inner;
                Vector2::new(result.x, result.y).into_bound_py_any(py)
            }
            Vector2OrPoint2::Point(other) => {
                let result = self.inner + other.inner.coords;
                Point2::new(result.x, result.y).into_bound_py_any(py)
            }
        }
    }

    fn __sub__(&self, other: Vector2) -> Self {
        Self {
            inner: self.inner - other.inner,
        }
    }

    fn __repr__(&self) -> String {
        format!("Vector2({}, {})", self.inner.x, self.inner.y)
    }

    fn dot(&self, other: Vector2) -> f64 {
        self.inner.dot(&other.inner)
    }

    fn cross(&self, other: Vector2) -> f64 {
        self.inner.cross(&other.inner)[0]
    }

    fn norm(&self) -> f64 {
        self.inner.norm()
    }

    fn normalized(&self) -> Self {
        Self {
            inner: self.inner.normalize(),
        }
    }

    fn angle_to(&self, other: Vector2) -> f64 {
        self.inner.angle(&other.inner)
    }

    fn with_x(&self, x: f64) -> Self {
        Self {
            inner: engeom::Vector2::new(x, self.inner.y),
        }
    }

    fn with_y(&self, y: f64) -> Self {
        Self {
            inner: engeom::Vector2::new(self.inner.x, y),
        }
    }
}

// ================================================================================================
// Points
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Point2 {
    inner: engeom::Point2,
}

impl Point2 {
    pub fn get_inner(&self) -> &engeom::Point2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Point2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Point2 {
    #[new]
    fn new(x: f64, y: f64) -> Self {
        Self {
            inner: engeom::Point2::new(x, y),
        }
    }

    fn __getstate__(&self) -> (f64, f64) {
        (self.inner.x, self.inner.y)
    }

    fn __getnewargs__(&self) -> (f64, f64) {
        (self.inner.x, self.inner.y)
    }

    fn __setstate__(&mut self, state: (f64, f64)) {
        self.inner = engeom::Point2::new(state.0, state.1);
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    #[staticmethod]
    fn origin() -> Self {
        Self::from_inner(engeom::Point2::origin())
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
    fn coords(&self) -> Vector2 {
        Vector2::from_inner(self.inner.coords)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let o = [self.inner.x, self.inner.y];
        PyIterator::from_object(&o.into_pyobject(py)?)
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let mut array = Array1::zeros(2);
        array[0] = self.inner.x;
        array[1] = self.inner.y;
        array.into_pyarray(py)
    }

    fn __add__(&self, other: Vector2) -> Self {
        Self {
            inner: self.inner + other.inner,
        }
    }

    fn __sub__<'py>(&self, py: Python<'py>, other: Vector2OrPoint2) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Vector2OrPoint2::Vector(other) => {
                let result = self.inner - other.inner;
                Point2::new(result.x, result.y).into_bound_py_any(py)
            }
            Vector2OrPoint2::Point(other) => {
                let result = self.inner - other.inner.coords;
                Vector2::new(result.x, result.y).into_bound_py_any(py)
            }
        }
    }

    fn __mul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __rmul__(&self, other: f64) -> Self {
        Self {
            inner: self.inner * other,
        }
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self {
            inner: self.inner / other,
        }
    }

    fn __neg__(&self) -> Self {
        Self { inner: -self.inner }
    }

    fn __repr__(&self) -> String {
        format!("Point2({}, {})", self.inner.x, self.inner.y)
    }

    #[staticmethod]
    fn mid(a: Point2, b: Point2) -> Self {
        Self::from_inner(engeom::common::points::mid_point(
            a.get_inner(),
            b.get_inner(),
        ))
    }

    fn with_x(&self, x: f64) -> Self {
        Self {
            inner: engeom::Point2::new(x, self.inner.y),
        }
    }

    fn with_y(&self, y: f64) -> Self {
        Self {
            inner: engeom::Point2::new(self.inner.x, y),
        }
    }
}

// ================================================================================================
// Surface Point
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct SurfacePoint2 {
    pub inner: engeom::SurfacePoint2,
}

impl SurfacePoint2 {
    pub fn get_inner(&self) -> &engeom::SurfacePoint2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::SurfacePoint2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl SurfacePoint2 {
    #[new]
    fn new(x: f64, y: f64, nx: f64, ny: f64) -> Self {
        Self {
            inner: engeom::SurfacePoint2::new_normalize(
                engeom::Point2::new(x, y),
                engeom::Vector2::new(nx, ny),
            ),
        }
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.point.x,
            self.inner.point.y,
            self.inner.normal.x,
            self.inner.normal.y,
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.point.x,
            self.inner.point.y,
            self.inner.normal.x,
            self.inner.normal.y,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) {
        self.inner = engeom::SurfacePoint2::new_normalize(
            engeom::Point2::new(state.0, state.1),
            engeom::Vector2::new(state.2, state.3),
        );
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.point == other.inner.point && self.inner.normal == other.inner.normal
    }

    #[getter]
    fn point(&self) -> Point2 {
        Point2::from_inner(self.inner.point)
    }

    #[getter]
    fn normal(&self) -> Vector2 {
        Vector2::from_inner(self.inner.normal.into_inner())
    }

    fn at_distance(&self, distance: f64) -> Point2 {
        Point2::from_inner(self.inner.at_distance(distance))
    }

    fn scalar_projection(&self, other: Point2) -> f64 {
        self.inner.scalar_projection(other.get_inner())
    }

    fn projection(&self, other: Point2) -> Point2 {
        Point2::from_inner(self.inner.projection(other.get_inner()))
    }

    fn reversed(&self) -> Self {
        Self::from_inner(self.inner.new_reversed())
    }

    fn transformed(&self, iso: Iso2) -> Self {
        Self::from_inner(self.inner.transformed(iso.get_inner()))
    }

    fn __mul__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint2::new_normalize(
            self.inner.point * other,
            self.inner.normal.into_inner() * other.signum(),
        ))
    }

    fn __rmul__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint2::new_normalize(
            self.inner.point * other,
            self.inner.normal.into_inner() * other.signum(),
        ))
    }

    fn __truediv__(&self, other: f64) -> Self {
        Self::from_inner(engeom::SurfacePoint2::new_normalize(
            self.inner.point / other,
            self.inner.normal.into_inner() / other.signum(),
        ))
    }

    fn __neg__(&self) -> Self {
        Self::from_inner(engeom::SurfacePoint2::new_normalize(
            -self.inner.point,
            -self.inner.normal.into_inner(),
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "SurfacePoint2({}, {}, {}, {})",
            self.inner.point.x, self.inner.point.y, self.inner.normal.x, self.inner.normal.y,
        )
    }

    fn planar_distance(&self, other: Point2) -> f64 {
        self.inner.planar_distance(other.get_inner())
    }

    fn shift_orthogonal(&self, distance: f64) -> Self {
        Self::from_inner(self.inner.shift_orthogonal(distance))
    }

    fn rot_normal(&self, angle: f64) -> Self {
        Self::from_inner(self.inner.rot_normal(angle))
    }

    fn new_shifted(&self, distance: f64) -> Self {
        Self::from_inner(self.inner.new_shifted(distance))
    }

    fn to_line(&self) -> Line2 {
        Line2::from_inner(engeom::geom2::Line2::from(&self.inner))
    }
}

// ================================================================================================
// Circle
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Circle2 {
    inner: engeom::Circle2,
}

impl Circle2 {
    pub fn get_inner(&self) -> &engeom::Circle2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Circle2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Circle2 {
    #[new]
    fn new(x: f64, y: f64, r: f64) -> Self {
        Self {
            inner: engeom::Circle2::new(x, y, r),
        }
    }

    fn __getstate__(&self) -> (f64, f64, f64) {
        (self.inner.center.x, self.inner.center.y, self.inner.r())
    }

    fn __getnewargs__(&self) -> (f64, f64, f64) {
        (self.inner.center.x, self.inner.center.y, self.inner.r())
    }

    fn __setstate__(&mut self, state: (f64, f64, f64)) {
        self.inner = engeom::Circle2::new(state.0, state.1, state.2);
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.center == other.inner.center && self.inner.r() == other.inner.r()
    }

    #[getter]
    fn center(&self) -> Point2 {
        Point2::from_inner(self.inner.center)
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.x()
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.y()
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    fn __repr__(&self) -> String {
        format!(
            "Circle2({}, {}, {})",
            self.inner.x(),
            self.inner.y(),
            self.inner.r()
        )
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(*self.inner.aabb())
    }

    fn point_at_angle(&self, angle: f64) -> Point2 {
        Point2::from_inner(self.inner.point_at_angle(angle))
    }

    fn project_point_to_perimeter(&self, other: &Point2) -> Option<Point2> {
        self.inner
            .project_point_to_perimeter(other.get_inner())
            .map(Point2::from_inner)
    }

    fn angle_of_point(&self, other: &Point2) -> f64 {
        self.inner.angle_of_point(other.get_inner())
    }

    fn intersections_with(&self, other: &Circle2) -> Vec<Point2> {
        self.inner
            .intersections_with(&other.inner)
            .into_iter()
            .map(Point2::from_inner)
            .collect()
    }

    fn distance_to(&self, point: &Point2) -> f64 {
        self.inner.distance_to(point.get_inner())
    }

    fn contains_point(&self, x: f64, y: f64) -> bool {
        self.inner.contains_point(&[x, y])
    }

    fn tangent_points_to(&self, point: &Point2) -> Option<(Point2, Point2)> {
        let pts = self
            .inner
            .tangent_points_to(point.get_inner())
            .map(|(p1, p2)| (Point2::from_inner(p1), Point2::from_inner(p2)))?;
        Some(pts)
    }

    #[staticmethod]
    #[pyo3(signature=(points, guess=None, sigma=None))]
    fn fitting<'py>(
        points: PyReadonlyArray2<'py, f64>,
        guess: Option<Circle2>,
        sigma: Option<f64>,
    ) -> PyResult<Self> {
        let points = array_to_points2(&points.as_array())?;
        let guess = if let Some(c) = guess {
            *c.get_inner()
        } else {
            engeom::Circle2::new(0.0, 0.0, 1.0)
        };

        let mode = if let Some(s) = sigma {
            BestFit::Gaussian(s)
        } else {
            BestFit::All
        };

        let circle = engeom::Circle2::new_fitting_circle(&points, &guess, mode)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(circle))
    }

    #[staticmethod]
    #[pyo3(signature=(points, tol, iterations=None, min_r=None, max_r=None))]
    fn ransac<'py>(
        points: PyReadonlyArray2<'py, f64>,
        tol: f64,
        iterations: Option<usize>,
        min_r: Option<f64>,
        max_r: Option<f64>,
    ) -> PyResult<Self> {
        let points = array_to_points2(&points.as_array())?;
        let result = engeom::Circle2::new_ransac(&points, tol, iterations, min_r, max_r)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    fn from_point(center: &Point2, r: f64) -> Self {
        Self::from_inner(engeom::Circle2::from_point(*center.get_inner(), r))
    }

    #[staticmethod]
    fn from_3_points(p0: &Point2, p1: &Point2, p2: &Point2) -> PyResult<Self> {
        engeom::Circle2::from_3_points(p0.get_inner(), p1.get_inner(), p2.get_inner())
            .map(Self::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn tangent_to_corner(
        corner: &Point2,
        d0: &Vector2,
        d1: &Vector2,
        radius: f64,
    ) -> PyResult<Self> {
        engeom::Circle2::new_tangent_to_corner(
            corner.get_inner(),
            d0.get_inner(),
            d1.get_inner(),
            radius,
        )
        .map(Self::from_inner)
        .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[staticmethod]
    fn tangent_and_point(tangent: &Line2, point: &Point2) -> Self {
        Self::from_inner(engeom::Circle2::new_tangent_and_point(
            tangent.get_inner(),
            point.get_inner(),
        ))
    }

    fn outer_tangents_to(&self, other: &Circle2) -> Option<(Segment2, Segment2)> {
        self.inner
            .outer_tangents_to(other.get_inner())
            .map(|(s0, s1)| (Segment2::from_inner(s0), Segment2::from_inner(s1)))
    }

    fn to_arc(&self) -> Arc2 {
        Arc2::from_inner(self.inner.to_arc())
    }

    fn to_partial_arc(&self, angle0: f64, angle: f64) -> Arc2 {
        Arc2::from_inner(self.inner.to_partial_arc(angle0, angle))
    }

    fn line_direction(&self, line: &Line2) -> AngleDir {
        self.inner.line_direction(line.get_inner()).into()
    }
}

// ================================================================================================
// Line2
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Line2 {
    inner: engeom::geom2::Line2,
}

impl Line2 {
    pub fn get_inner(&self) -> &engeom::geom2::Line2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom2::Line2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Line2 {
    #[new]
    fn new(ox: f64, oy: f64, dx: f64, dy: f64) -> Self {
        Self::from_inner(engeom::geom2::Line2::new(
            engeom::Point2::new(ox, oy),
            engeom::Vector2::new(dx, dy),
        ))
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.origin.x,
            self.inner.origin.y,
            self.inner.direction.x,
            self.inner.direction.y,
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.origin.x,
            self.inner.origin.y,
            self.inner.direction.x,
            self.inner.direction.y,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) {
        self.inner = engeom::geom2::Line2::new(
            engeom::Point2::new(state.0, state.1),
            engeom::Vector2::new(state.2, state.3),
        );
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.origin == other.inner.origin && self.inner.direction == other.inner.direction
    }

    fn __repr__(&self) -> String {
        format!(
            "Line2({}, {}, {}, {})",
            self.inner.origin.x,
            self.inner.origin.y,
            self.inner.direction.x,
            self.inner.direction.y,
        )
    }

    #[staticmethod]
    fn x_axis() -> Self {
        Self::from_inner(engeom::geom2::Line2::x_axis())
    }

    #[staticmethod]
    fn y_axis() -> Self {
        Self::from_inner(engeom::geom2::Line2::y_axis())
    }

    #[staticmethod]
    fn from_points(p1: &Point2, p2: &Point2) -> Self {
        Self::from_inner(engeom::geom2::Line2::from_points(
            p1.get_inner(),
            p2.get_inner(),
        ))
    }

    #[staticmethod]
    fn new_normalize(ox: f64, oy: f64, dx: f64, dy: f64) -> Self {
        Self::from_inner(engeom::geom2::Line2::new_normalize(
            engeom::Point2::new(ox, oy),
            engeom::Vector2::new(dx, dy),
        ))
    }

    #[getter]
    fn origin(&self) -> Point2 {
        Point2::from_inner(self.inner.origin)
    }

    #[getter]
    fn direction(&self) -> Vector2 {
        Vector2::from_inner(self.inner.direction)
    }

    #[getter]
    fn normal(&self) -> Vector2 {
        Vector2::from_inner(self.inner.normal().into_inner())
    }

    fn at(&self, t: f64) -> Point2 {
        Point2::from_inner(self.inner.at(t))
    }

    fn scalar_project(&self, point: &Point2) -> f64 {
        self.inner.scalar_project(point.get_inner())
    }

    fn closest_point(&self, point: &Point2) -> Point2 {
        Point2::from_inner(self.inner.closest_point(point.get_inner()))
    }

    fn distance_to(&self, point: &Point2) -> f64 {
        self.inner.distance_to(point.get_inner())
    }

    fn signed_distance_to(&self, point: &Point2) -> f64 {
        self.inner.signed_distance_to(point.get_inner())
    }

    fn intersect(&self, other: &Line2) -> Option<Point2> {
        self.inner.intersect(&other.inner).map(Point2::from_inner)
    }

    fn normalized(&self) -> Self {
        Self::from_inner(self.inner.normalized())
    }

    fn new_parallel(&self, delta_n: f64) -> Self {
        Self::from_inner(self.inner.new_parallel(delta_n))
    }

    fn new_shifted_along(&self, delta_t: f64) -> Self {
        Self::from_inner(self.inner.new_shifted_along(delta_t))
    }

    fn new_transformed_by(&self, iso: &Iso2) -> Self {
        Self::from_inner(self.inner.new_transformed_by(iso.get_inner()))
    }

    fn to_iso_from_x(&self) -> Iso2 {
        use engeom::geom2::LineOps2;
        Iso2::from_inner(self.inner.to_iso_from_x())
    }

    fn to_iso_from_y(&self) -> Iso2 {
        use engeom::geom2::LineOps2;
        Iso2::from_inner(self.inner.to_iso_from_y())
    }

    fn to_surface_point(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(engeom::SurfacePoint2::from(self.inner))
    }
}

// ================================================================================================
// Segment
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Segment2 {
    inner: engeom::geom2::Segment2,
}

impl Segment2 {
    pub fn get_inner(&self) -> &engeom::geom2::Segment2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom2::Segment2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Segment2 {
    #[new]
    fn new(x0: f64, y0: f64, x1: f64, y1: f64) -> PyResult<Self> {
        let p0 = engeom::Point2::new(x0, y0);
        let p1 = engeom::Point2::new(x1, y1);
        Ok(Self {
            inner: engeom::geom2::Segment2::try_new(&p0, &p1)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        })
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.a.x,
            self.inner.a.y,
            self.inner.b.x,
            self.inner.b.y,
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64) {
        (
            self.inner.a.x,
            self.inner.a.y,
            self.inner.b.x,
            self.inner.b.y,
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64)) {
        let p0 = engeom::Point2::new(state.0, state.1);
        let p1 = engeom::Point2::new(state.2, state.3);
        self.inner = engeom::geom2::Segment2::try_new(&p0, &p1)
            .expect("Invalid segment points in __setstate__");
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.a == other.inner.a && self.inner.b == other.inner.b
    }

    fn __repr__(&self) -> String {
        format!(
            "Segment2({}, {}, {}, {})",
            self.inner.a.x, self.inner.a.y, self.inner.b.x, self.inner.b.y
        )
    }

    #[getter]
    fn a(&self) -> Point2 {
        Point2::from_inner(self.inner.a)
    }

    #[getter]
    fn b(&self) -> Point2 {
        Point2::from_inner(self.inner.b)
    }

    #[getter]
    fn direction(&self) -> Vector2 {
        Vector2::from_inner(self.inner.dir())
    }

    #[getter]
    fn length(&self) -> f64 {
        self.inner.length
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(self.inner.aabb())
    }

    fn to_line(&self) -> Line2 {
        Line2::from_inner(self.inner.to_line())
    }
}

// ================================================================================================
// CubicSpline2
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct CubicSpline2 {
    inner: engeom::geom2::CubicSpline2,
}

impl CubicSpline2 {
    pub fn get_inner(&self) -> &engeom::geom2::CubicSpline2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom2::CubicSpline2) -> Self {
        Self { inner }
    }
}

type CubicSpline2State = (f64, f64, f64, f64, f64, f64, f64, f64);

#[pymethods]
impl CubicSpline2 {
    #[new]
    fn new(x0: f64, y0: f64, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64) -> Self {
        Self {
            inner: engeom::geom2::CubicSpline2::new(
                engeom::Point2::new(x0, y0),
                engeom::Point2::new(x1, y1),
                engeom::Point2::new(x2, y2),
                engeom::Point2::new(x3, y3),
            ),
        }
    }

    fn __getstate__(&self) -> CubicSpline2State {
        (
            self.inner.p0.x,
            self.inner.p0.y,
            self.inner.p1.x,
            self.inner.p1.y,
            self.inner.p2.x,
            self.inner.p2.y,
            self.inner.p3.x,
            self.inner.p3.y,
        )
    }

    fn __getnewargs__(&self) -> CubicSpline2State {
        self.__getstate__()
    }

    fn __setstate__(&mut self, state: CubicSpline2State) {
        self.inner = engeom::geom2::CubicSpline2::new(
            engeom::Point2::new(state.0, state.1),
            engeom::Point2::new(state.2, state.3),
            engeom::Point2::new(state.4, state.5),
            engeom::Point2::new(state.6, state.7),
        );
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!(
            "CubicSpline2({}, {}, {}, {}, {}, {}, {}, {})",
            self.inner.p0.x,
            self.inner.p0.y,
            self.inner.p1.x,
            self.inner.p1.y,
            self.inner.p2.x,
            self.inner.p2.y,
            self.inner.p3.x,
            self.inner.p3.y,
        )
    }

    #[getter]
    fn p0(&self) -> Point2 {
        Point2::from_inner(self.inner.p0)
    }

    #[getter]
    fn p1(&self) -> Point2 {
        Point2::from_inner(self.inner.p1)
    }

    #[getter]
    fn p2(&self) -> Point2 {
        Point2::from_inner(self.inner.p2)
    }

    #[getter]
    fn p3(&self) -> Point2 {
        Point2::from_inner(self.inner.p3)
    }

    fn position(&self, t: f64) -> Point2 {
        Point2::from_inner(self.inner.position(t))
    }

    fn derivative(&self, t: f64) -> Vector2 {
        Vector2::from_inner(self.inner.derivative(t))
    }

    fn tangent(&self, t: f64) -> Vector2 {
        Vector2::from_inner(self.inner.tangent(t).into_inner())
    }

    fn normal(&self, t: f64) -> Vector2 {
        Vector2::from_inner(self.inner.normal(t).into_inner())
    }

    fn curvature(&self, t: f64) -> f64 {
        self.inner.curvature(t)
    }

    fn second_derivative(&self, t: f64) -> Vector2 {
        Vector2::from_inner(self.inner.second_derivative(t))
    }

    fn polyline<'py>(&self, py: Python<'py>, tolerance: f64) -> Bound<'py, PyArray2<f64>> {
        let points = self.inner.polyline(tolerance);
        points_to_array(&points).into_pyarray(py)
    }
}

// ================================================================================================
// Arc
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Arc2 {
    inner: engeom::Arc2,
}

impl Arc2 {
    pub fn get_inner(&self) -> &engeom::Arc2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Arc2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Arc2 {
    fn __repr__(&self) -> String {
        format!(
            "Arc2({}, {}, {}, {}, {})",
            self.inner.center().x,
            self.inner.center().y,
            self.inner.radius(),
            self.inner.angle0(),
            self.inner.angle(),
        )
    }

    fn __getstate__(&self) -> (f64, f64, f64, f64, f64) {
        (
            self.inner.center().x,
            self.inner.center().y,
            self.inner.radius(),
            self.inner.angle0(),
            self.inner.angle(),
        )
    }

    fn __getnewargs__(&self) -> (f64, f64, f64, f64, f64) {
        (
            self.inner.center().x,
            self.inner.center().y,
            self.inner.radius(),
            self.inner.angle0(),
            self.inner.angle(),
        )
    }

    fn __setstate__(&mut self, state: (f64, f64, f64, f64, f64)) {
        self.inner = engeom::Arc2::circle_angles(
            engeom::Point2::new(state.0, state.1),
            state.2,
            state.3,
            state.4,
        );
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner.center() == other.inner.center()
            && self.inner.radius() == other.inner.radius()
            && self.inner.angle0() == other.inner.angle0()
            && self.inner.angle() == other.inner.angle()
    }

    #[new]
    fn new(x: f64, y: f64, r: f64, start_radians: f64, sweep_radians: f64) -> Self {
        Self {
            inner: engeom::Arc2::circle_angles(
                engeom::Point2::new(x, y),
                r,
                start_radians,
                sweep_radians,
            ),
        }
    }

    #[getter]
    fn center(&self) -> Point2 {
        Point2::from_inner(self.inner.center())
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.radius()
    }

    #[getter]
    fn angle0(&self) -> f64 {
        self.inner.angle0()
    }

    #[getter]
    fn angle(&self) -> f64 {
        self.inner.angle()
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(*self.inner.aabb())
    }

    #[getter]
    fn a(&self) -> Point2 {
        Point2::from_inner(self.inner.a())
    }

    #[getter]
    fn b(&self) -> Point2 {
        Point2::from_inner(self.inner.b())
    }

    #[getter]
    fn circle(&self) -> Circle2 {
        Circle2::from_inner(self.inner.circle())
    }

    fn make_points<'py>(&self, py: Python<'py>, tol: f64) -> Bound<'py, PyArray2<f64>> {
        let points = self.inner.make_points(tol);
        points_to_array(&points).into_pyarray(py)
    }
}

// ================================================================================================
// Curve
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct CurveStation2 {
    i_point: engeom::Point2,
    i_direction: engeom::Vector2,
    i_index: usize,
    i_fraction: f64,
    i_length_along: f64,
    i_normal: engeom::Vector2,
}

impl CurveStation2 {
    pub fn new(
        point: engeom::Point2,
        direction: engeom::Vector2,
        index: usize,
        fraction: f64,
        length_along: f64,
        normal: engeom::Vector2,
    ) -> Self {
        Self {
            i_point: point,
            i_direction: direction,
            i_index: index,
            i_fraction: fraction,
            i_length_along: length_along,
            i_normal: normal,
        }
    }
}

#[pymethods]
impl CurveStation2 {
    #[getter]
    pub fn point(&self) -> Point2 {
        Point2::from_inner(self.i_point)
    }

    #[getter]
    pub fn direction(&self) -> Vector2 {
        Vector2::from_inner(self.i_direction)
    }

    #[getter]
    pub fn direction_point(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(engeom::SurfacePoint2::new_normalize(
            self.i_point,
            self.i_direction,
        ))
    }

    #[getter]
    pub fn surface_point(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(engeom::SurfacePoint2::new_normalize(
            self.i_point,
            self.i_normal,
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

    #[getter]
    pub fn normal(&self) -> Vector2 {
        Vector2::from_inner(self.i_normal)
    }
}

impl From<engeom::CurveStation2<'_>> for CurveStation2 {
    fn from(station: engeom::CurveStation2) -> Self {
        Self::new(
            station.point(),
            station.direction().into_inner(),
            station.index(),
            station.fraction(),
            station.length_along(),
            station.normal().into_inner(),
        )
    }
}

#[pyclass(from_py_object, module = "engeom.geom2")]
pub struct Curve2 {
    inner: engeom::Curve2,
    points: Option<Py<PyArray2<f64>>>,
}

impl Curve2 {
    pub fn get_inner(&self) -> &engeom::Curve2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Curve2) -> Self {
        Self {
            inner,
            points: None,
        }
    }
}

impl Clone for Curve2 {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            points: None,
        }
    }
}

#[pymethods]
impl Curve2 {
    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let result = points_to_array(self.inner.points());
            self.points = Some(result.into_pyarray(py).unbind())
        }

        self.points.as_ref().unwrap().bind(py)
    }

    #[new]
    #[pyo3(signature=(points, normals=None, tol=1e-6, force_closed=false, hull_ccw=false))]
    fn new(
        points: PyReadonlyArray2<'_, f64>,
        normals: Option<PyReadonlyArray2<'_, f64>>,
        tol: f64,
        force_closed: bool,
        hull_ccw: bool,
    ) -> PyResult<Self> {
        let points = array_to_points2(&points.as_array())?;

        let curve = if let Some(normal_array) = normals {
            let normals = array_to_vectors2(&normal_array.as_array())?;
            if points.len() != normals.len() {
                return Err(PyValueError::new_err(
                    "Points and normals must have the same length",
                ));
            }

            let surf_points = points
                .iter()
                .zip(normals.iter())
                .map(|(p, n)| engeom::SurfacePoint2::new_normalize(*p, *n))
                .collect::<Vec<_>>();

            engeom::Curve2::from_surf_points(&surf_points, tol, force_closed)
        } else if hull_ccw {
            engeom::Curve2::from_points_ccw(&points, tol, force_closed)
        } else {
            engeom::Curve2::from_points(&points, tol, force_closed)
        }
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(curve))
    }

    fn length(&self) -> f64 {
        self.inner.length()
    }

    fn at_front(&self) -> CurveStation2 {
        self.inner.at_front().into()
    }

    fn at_back(&self) -> CurveStation2 {
        self.inner.at_back().into()
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(self.inner.aabb())
    }

    fn at_length(&self, length: f64) -> PyResult<CurveStation2> {
        self.inner
            .at_length(length)
            .map(|s| s.into())
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn at_fraction(&self, fraction: f64) -> PyResult<CurveStation2> {
        self.inner
            .at_fraction(fraction)
            .map(|s| s.into())
            .ok_or_else(|| PyValueError::new_err("Fraction out of bounds"))
    }

    fn at_closest_to_point(&self, point: Point2) -> CurveStation2 {
        self.inner.at_closest_to_point(point.get_inner()).into()
    }

    #[getter]
    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    fn trim_front(&self, length: f64) -> PyResult<Self> {
        self.inner
            .trim_front(length)
            .map(Self::from_inner)
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn trim_back(&self, length: f64) -> PyResult<Self> {
        self.inner
            .trim_back(length)
            .map(Self::from_inner)
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn between_lengths_by_control(&self, a: f64, b: f64, control: f64) -> PyResult<Self> {
        self.inner
            .between_lengths_by_control(a, b, control)
            .map(Self::from_inner)
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn between_lengths(&self, l0: f64, l1: f64) -> PyResult<Self> {
        self.inner
            .between_lengths(l0, l1)
            .map(Self::from_inner)
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn reversed(&self) -> Self {
        Self::from_inner(self.inner.reversed())
    }

    fn make_hull<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let hull = self
            .inner
            .make_hull()
            .ok_or(PyValueError::new_err("Could not compute convex hull"))?;

        let result = points_to_array(hull.points());
        Ok(result.into_pyarray(py))
    }

    fn max_point_in_direction(&self, direction: Vector2) -> PyResult<(usize, Point2)> {
        let (i, p) = self
            .inner
            .max_point_in_direction(direction.get_inner())
            .ok_or(PyValueError::new_err(
                "Could not compute max point in direction",
            ))?;
        Ok((i, Point2::from_inner(p)))
    }

    fn max_dist_in_direction(&self, surf_point: SurfacePoint2) -> f64 {
        self.inner.max_dist_in_direction(surf_point.get_inner())
    }

    fn simplify(&self, tol: f64) -> Self {
        Self::from_inner(self.inner.simplify(tol))
    }

    fn resample(&self, resample: Resample) -> PyResult<Self> {
        let inner = self
            .inner
            .resample(resample.into())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn new_transformed_by(&self, iso: &Iso2) -> Self {
        Self::from_inner(self.inner.new_transformed_by(iso.get_inner()))
    }

    fn to_3d(&self) -> PyResult<crate::geom3::Curve3> {
        let points = self.inner.points().to_3d();
        let c = engeom::Curve3::from_points(&points, self.inner.tol())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(crate::geom3::Curve3::from_inner(c))
    }

    fn offset_vertices(&self, offset: f64) -> PyResult<Self> {
        let inner = self
            .inner
            .offset_vertices(offset)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn offset_segments(&self, offset: f64) -> PyResult<Self> {
        let inner = self
            .inner
            .offset_segments(offset)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn __add__(&self, other: &Self) -> PyResult<Self> {
        let result = self
            .get_inner()
            .new_appended(other.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }

    #[staticmethod]
    fn load_tccurve2(path: PathBuf) -> PyResult<Self> {
        let curve = engeom::io::read_tc_curve2_file(&path)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(curve))
    }

    fn write_tccurve2(&self, path: PathBuf, tol: f64) -> PyResult<()> {
        engeom::io::write_tc_curve2_file(&path, &self.inner, tol)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "<Curve2 n={}, l={} ({})>",
            self.inner.points().len(),
            self.inner.length(),
            if self.inner.is_closed() {
                "closed"
            } else {
                "open"
            }
        )
    }
}

// ================================================================================================
// Transformations
// ================================================================================================

#[derive(FromPyObject)]
enum Transformable2 {
    Iso(Iso2),
    Vec(Vector2),
    Pnt(Point2),
    Sp(SurfacePoint2),
}

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Iso2 {
    inner: engeom::Iso2,
}

impl Iso2 {
    pub fn get_inner(&self) -> &engeom::Iso2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Iso2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Iso2 {
    #[new]
    fn new(tx: f64, ty: f64, r: f64) -> Self {
        let inner = engeom::Iso2::translation(tx, ty) * engeom::Iso2::rotation(r);
        Self { inner }
    }

    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: engeom::Iso2::identity(),
        }
    }

    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Iso2({}, {}, {})",
            self.inner.translation.x,
            self.inner.translation.y,
            self.inner.rotation.angle()
        )
    }

    fn __matmul__<'py>(
        &self,
        py: Python<'py>,
        other: Transformable2,
    ) -> PyResult<Bound<'py, PyAny>> {
        match other {
            Transformable2::Iso(other) => {
                Iso2::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable2::Vec(other) => {
                Vector2::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable2::Pnt(other) => {
                Point2::from_inner(self.inner * other.inner).into_bound_py_any(py)
            }
            Transformable2::Sp(other) => {
                SurfacePoint2::from_inner(other.inner.transformed(&self.inner))
                    .into_bound_py_any(py)
            }
        }
    }

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut result = Array2::zeros((3, 3));
        let m = self.inner.to_matrix();
        result[[0, 0]] = m.m11;
        result[[0, 1]] = m.m12;
        result[[0, 2]] = m.m13;
        result[[1, 0]] = m.m21;
        result[[1, 1]] = m.m22;
        result[[1, 2]] = m.m23;
        result[[2, 0]] = m.m31;
        result[[2, 1]] = m.m32;
        result[[2, 2]] = m.m33;
        result.into_pyarray(py)
    }

    fn transform_points<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let points = array_to_points2(&points.as_array())?;
        let transformed = points.iter().map(|p| self.inner * p).collect::<Vec<_>>();
        let result = points_to_array(&transformed);
        Ok(result.into_pyarray(py))
    }

    fn transform_vectors<'py>(
        &self,
        py: Python<'py>,
        vectors: PyReadonlyArray2<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let vectors = array_to_vectors2(&vectors.as_array())?;
        let transformed = vectors.iter().map(|v| self.inner * v).collect::<Vec<_>>();
        let result = vectors_to_array(&transformed);
        Ok(result.into_pyarray(py))
    }

    #[getter]
    fn x_direction(&self) -> Vector2 {
        Vector2::from_inner(self.inner.rotation * engeom::Vector2::x())
    }

    #[getter]
    fn y_direction(&self) -> Vector2 {
        Vector2::from_inner(self.inner.rotation * engeom::Vector2::y())
    }

    #[getter]
    fn x_axis(&self) -> SurfacePoint2 {
        let origin = self.inner * engeom::Point2::origin();
        let direction = self.inner.rotation * engeom::Vector2::x();
        SurfacePoint2::from_inner(engeom::SurfacePoint2::new_normalize(origin, direction))
    }

    #[getter]
    fn y_axis(&self) -> SurfacePoint2 {
        let origin = self.inner * engeom::Point2::origin();
        let direction = self.inner.rotation * engeom::Vector2::y();
        SurfacePoint2::from_inner(engeom::SurfacePoint2::new_normalize(origin, direction))
    }
}

// ================================================================================================
// Hull functions
// ================================================================================================

/// Computes the convex hull of a set of 2D points, returning a numpy array of the hull vertex
/// coordinates ordered counter-clockwise.
#[pyfunction]
pub fn convex_hull_2d<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let pts = array_to_points2(&points.as_array())?;
    let indices = engeom::geom2::hull::convex_hull_2d(&pts);
    let hull_points: Vec<engeom::Point2> = indices.iter().map(|&i| pts[i]).collect();
    Ok(points_to_array(&hull_points).into_pyarray(py))
}

/// Computes the convex hull of a set of 2D points, returning a numpy array of indices into
/// the input array ordered counter-clockwise.
#[pyfunction]
pub fn convex_hull_idx<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, PyArray1<usize>>> {
    let pts = array_to_points2(&points.as_array())?;
    let indices = engeom::geom2::hull::convex_hull_2d(&pts);
    Ok(Array1::from_vec(indices).into_pyarray(py))
}
