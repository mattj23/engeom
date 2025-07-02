use crate::bounding::Aabb2;
use crate::common::Resample;
use crate::conversions::{
    array_to_points2, array_to_vectors2, array2_to_points2, points_to_array2,
};
use crate::geom3::Point3;
use engeom::airfoil::OpenEdge;
use engeom::geom2::{HasBounds2, Line2};
use engeom::{BestFit, To3D};
use numpy::ndarray::{Array1, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArrayDyn, PyReadonlyArray2, PyReadonlyArrayDyn};
use pyo3::exceptions::PyValueError;
use pyo3::types::PyIterator;
use pyo3::{
    Bound, FromPyObject, IntoPyObject, IntoPyObjectExt, Py, PyAny, PyResult, Python, pyclass,
    pymethods,
};

#[derive(FromPyObject)]
enum Vector2OrPoint2 {
    Vector(Vector2),
    Point(Point2),
}

// ================================================================================================
// Vectors
// ================================================================================================

#[pyclass]
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

#[pyclass]
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
#[pyclass]
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
        Self::from_inner(self.inner.reversed())
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

    fn shift(&self, distance: f64) -> Self {
        Self::from_inner(self.inner.shift(distance))
    }
}

// ================================================================================================
// Circle
// ================================================================================================
#[pyclass]
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

    #[staticmethod]
    #[pyo3(signature=(points, guess=None, sigma=None))]
    fn fitting<'py>(
        points: PyReadonlyArray2<'py, f64>,
        guess: Option<Circle2>,
        sigma: Option<f64>,
    ) -> PyResult<Self> {
        let points = array2_to_points2(&points.as_array())?;
        let guess = if let Some(c) = guess {
            c.get_inner().clone()
        } else {
            engeom::Circle2::new(0.0, 0.0, 1.0)
        };

        let mode = if let Some(s) = sigma {
            BestFit::Gaussian(s)
        } else {
            BestFit::All
        };

        let circle = engeom::Circle2::fitting_circle(&points, &guess, mode)
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
        let points = array2_to_points2(&points.as_array())?;
        let result = engeom::Circle2::ransac(&points, tol, iterations, min_r, max_r)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
    }
}

// ================================================================================================
// Segment
// ================================================================================================
// TODO: Type hints for this

#[pyclass]
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
            inner: engeom::geom2::Segment2::try_new(p0, p1)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        })
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
}

// ================================================================================================
// Arc
// ================================================================================================
#[pyclass]
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
            self.inner.angle0,
            self.inner.angle
        )
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
    fn x(&self) -> f64 {
        self.inner.center().x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.center().y
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.radius()
    }

    #[getter]
    fn start(&self) -> f64 {
        self.inner.angle0
    }

    #[getter]
    fn sweep(&self) -> f64 {
        self.inner.angle
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(*self.inner.aabb())
    }

    #[getter]
    fn start_point(&self) -> Point2 {
        Point2::from_inner(self.inner.start())
    }

    #[getter]
    fn end_point(&self) -> Point2 {
        Point2::from_inner(self.inner.end())
    }
}

// ================================================================================================
// Curve
// ================================================================================================

#[pyclass]
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

#[pyclass]
pub struct Curve2 {
    inner: engeom::Curve2,
    points: Option<Py<PyArrayDyn<f64>>>,
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
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArrayDyn<f64>> {
        if self.points.is_none() {
            let result = points_to_array2(self.inner.points());
            self.points = Some(result.into_pyarray(py).unbind())
        }

        self.points.as_ref().unwrap().bind(py)
    }

    #[new]
    #[pyo3(signature=(points, normals=None, tol=1e-6, force_closed=false, hull_ccw=false))]
    fn new(
        points: PyReadonlyArrayDyn<'_, f64>,
        normals: Option<PyReadonlyArrayDyn<'_, f64>>,
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
        Aabb2::from_inner(*self.inner.aabb())
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

    fn make_hull<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let hull = self
            .inner
            .make_hull()
            .ok_or(PyValueError::new_err("Could not compute convex hull"))?;

        let mut result = ArrayD::zeros(vec![hull.points().len(), 2]);
        for (i, point) in hull.points().iter().enumerate() {
            result[[i, 0]] = point.x;
            result[[i, 1]] = point.y;
        }
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

    fn transformed_by(&self, iso: &Iso2) -> Self {
        Self::from_inner(self.inner.transformed_by(iso.get_inner()))
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
            .extended(other.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(result))
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

#[pyclass]
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

    fn as_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArrayDyn<f64>> {
        let mut result = ArrayD::zeros(vec![3, 3]);
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
        points: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let points = array_to_points2(&points.as_array())?;
        let transformed = points.iter().map(|p| self.inner * p).collect::<Vec<_>>();
        let mut result = ArrayD::zeros(vec![transformed.len(), 2]);
        for (i, point) in transformed.iter().enumerate() {
            result[[i, 0]] = point.x;
            result[[i, 1]] = point.y;
        }
        Ok(result.into_pyarray(py))
    }

    fn transform_vectors<'py>(
        &self,
        py: Python<'py>,
        vectors: PyReadonlyArrayDyn<'py, f64>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let vectors = array_to_vectors2(&vectors.as_array())?;
        let transformed = vectors.iter().map(|v| self.inner * v).collect::<Vec<_>>();
        let mut result = ArrayD::zeros(vec![transformed.len(), 2]);
        for (i, vector) in transformed.iter().enumerate() {
            result[[i, 0]] = vector.x;
            result[[i, 1]] = vector.y;
        }
        Ok(result.into_pyarray(py))
    }
}
