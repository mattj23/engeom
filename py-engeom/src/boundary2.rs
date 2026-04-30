use crate::bounding::Aabb2;
use crate::conversions::points_to_array;
use crate::geom2::{Iso2, Line2, Point2, SurfacePoint2, Vector2};
use engeom::geom2::BoundaryEditor;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::{Bound, PyResult, Python, pyclass, pymethods};

// ================================================================================================
// ManifoldPosition2
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.geom2")]
#[derive(Clone, Debug)]
pub struct Manifold1Pos2 {
    inner: engeom::geom2::Manifold1Pos2,
}

impl Manifold1Pos2 {
    pub fn from_inner(inner: engeom::geom2::Manifold1Pos2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Manifold1Pos2 {
    #[getter]
    fn l(&self) -> f64 {
        self.inner.l
    }

    #[getter]
    fn point(&self) -> Point2 {
        Point2::from_inner(self.inner.point)
    }

    #[getter]
    fn direction(&self) -> Vector2 {
        Vector2::from_inner(self.inner.direction.into_inner())
    }

    #[getter]
    fn normal(&self) -> Vector2 {
        Vector2::from_inner(self.inner.normal.into_inner())
    }

    #[getter]
    fn surface_point(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(engeom::SurfacePoint2::new_normalize(
            self.inner.point,
            self.inner.normal.into_inner(),
        ))
    }

    #[getter]
    fn direction_line(&self) -> Line2 {
        Line2::from_inner(engeom::Line2::new(
            self.inner.point,
            self.inner.direction.into_inner(),
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "ManifoldPosition2(l={}, point=({}, {}))",
            self.inner.l, self.inner.point.x, self.inner.point.y
        )
    }
}

// ================================================================================================
// BoundaryData2
// ================================================================================================

#[pyclass(module = "engeom.geom2")]
pub struct BoundaryData2 {
    inner: engeom::geom2::BoundaryData2,
}

impl BoundaryData2 {
    pub fn from_inner(inner: engeom::geom2::BoundaryData2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl BoundaryData2 {
    #[new]
    #[pyo3(signature=(x=None, y=None, closed=false))]
    fn new(x: Option<f64>, y: Option<f64>, closed: bool) -> PyResult<Self> {
        if closed {
            Ok(Self::new_closed())
        } else {
            match (x, y) {
                (Some(x), Some(y)) => Ok(Self::new_open(x, y)),
                _ => Err(PyValueError::new_err(
                    "Must provide x and y for open boundary, or specify closed",
                )),
            }
        }
    }

    #[staticmethod]
    fn new_open(x: f64, y: f64) -> Self {
        Self::from_inner(engeom::geom2::BoundaryData2::new_open_xy(x, y))
    }

    #[staticmethod]
    fn new_closed() -> Self {
        Self::from_inner(engeom::geom2::BoundaryData2::new_closed())
    }

    fn add_seg_xy(&mut self, x: f64, y: f64) -> u32 {
        self.inner.add_seg_xy(x, y)
    }

    fn add_arc_xy(&mut self, cx: f64, cy: f64, ex: f64, ey: f64, clockwise: bool) -> u32 {
        self.inner.add_arc_xy(cx, cy, ex, ey, clockwise)
    }

    fn add_corner_fillets(&mut self, points: Vec<(f64, f64)>, radius: f64) -> PyResult<Vec<u32>> {
        let pts = points
            .iter()
            .map(|(x, y)| engeom::Point2::new(*x, *y))
            .collect::<Vec<_>>();
        self.inner
            .add_corner_fillets(&pts, radius)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn transform_by(&mut self, iso: &Iso2) {
        self.inner.transform_by(iso.get_inner());
    }

    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn to_boundary(&self) -> PyResult<Boundary2> {
        self.inner
            .try_to_boundary()
            .map(Boundary2::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "<BoundaryData2 n={} {}>",
            self.inner.len(),
            if self.inner.is_closed() {
                "closed"
            } else {
                "open"
            }
        )
    }
}

// ================================================================================================
// Boundary2
// ================================================================================================

#[pyclass(module = "engeom.geom2")]
pub struct Boundary2 {
    inner: engeom::geom2::Boundary2,
}

impl Boundary2 {
    pub fn from_inner(inner: engeom::geom2::Boundary2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Boundary2 {
    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    fn length(&self) -> f64 {
        self.inner.length()
    }

    fn at_length(&self, length: f64) -> PyResult<Manifold1Pos2> {
        self.inner
            .at_length(length)
            .map(Manifold1Pos2::from_inner)
            .ok_or_else(|| PyValueError::new_err("Length out of bounds"))
    }

    fn at_start(&self) -> Manifold1Pos2 {
        Manifold1Pos2::from_inner(self.inner.at_start())
    }

    fn at_end(&self) -> Manifold1Pos2 {
        Manifold1Pos2::from_inner(self.inner.at_end())
    }

    fn at_closest_to_point(&self, point: &Point2) -> (u32, Manifold1Pos2) {
        let (id, m) = self.inner.at_closest_to_point(point.get_inner());
        (id, Manifold1Pos2::from_inner(m))
    }

    fn to_points<'py>(&self, py: Python<'py>, tol: f64) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let pts = self
            .inner
            .to_points(tol)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(points_to_array(&pts).into_pyarray(py))
    }

    #[getter]
    fn aabb(&self) -> Aabb2 {
        Aabb2::from_inner(self.inner.aabb())
    }

    fn __repr__(&self) -> String {
        format!(
            "<Boundary2 l={} ({})>",
            self.inner.length(),
            if self.inner.is_closed() {
                "closed"
            } else {
                "open"
            }
        )
    }
}
