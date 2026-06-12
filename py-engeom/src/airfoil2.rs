use crate::geom2::{Curve2, Point2, SurfacePoint2, Vector2};
use engeom::airfoil2::SectionInput;
use engeom::airfoil2::edges::AfEdgeFit as InnerAfEdgeFit;
use engeom::airfoil2::inscribed::Inscribed as InnerInscribed;
use engeom::airfoil2::{AfEdge as InnerAfEdge, AfEdgeGeometry as InnerAfEdgeGeometry};
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// ================================================================================================
// Inscribed
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Debug)]
pub struct Inscribed {
    inner: InnerInscribed,
}

impl Inscribed {
    pub fn from_inner(inner: InnerInscribed) -> Self {
        Self { inner }
    }

    pub fn into_inner(self) -> InnerInscribed {
        self.inner
    }
}

#[pymethods]
impl Inscribed {
    #[getter]
    fn c(&self) -> crate::geom2::Circle2 {
        crate::geom2::Circle2::from_inner(self.inner.c)
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
    fn center(&self) -> Point2 {
        Point2::from_inner(self.inner.center())
    }

    #[getter]
    fn radius(&self) -> f64 {
        self.inner.radius()
    }

    fn camber_point(&self) -> SurfacePoint2 {
        SurfacePoint2::from_inner(self.inner.camber_point())
    }

    fn contact_dir(&self) -> Vector2 {
        Vector2::from_inner(self.inner.contact_dir())
    }

    fn __repr__(&self) -> String {
        format!(
            "Inscribed(c=Circle2({}, {}, {}), p0=Point2({}, {}), p1=Point2({}, {}))",
            self.inner.c.center.x,
            self.inner.c.center.y,
            self.inner.c.r(),
            self.inner.p0.x,
            self.inner.p0.y,
            self.inner.p1.x,
            self.inner.p1.y,
        )
    }
}

// ================================================================================================
// AfEdgeGeometry
// ================================================================================================

/// The geometric description of a fitted airfoil edge.
///
/// The `kind` property identifies the variant. Depending on `kind`, other properties are
/// populated and the rest are `None`:
///
/// - `"open"`: no geometric data.
/// - `"sharp"`: `point` is the apex.
/// - `"square"`: `corner0` and `corner1`.
/// - `"rounded_square"`: `corner0`, `corner1`, and `radius`.
/// - `"full_round"`: `point` is the circle center, `radius` is the circle radius.
/// - `"blended_round"`: `point` is the circle center, `radius` is the circle radius.
#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Debug)]
pub struct AfEdgeGeometry {
    inner: InnerAfEdgeGeometry,
}

impl AfEdgeGeometry {
    pub fn from_inner(inner: InnerAfEdgeGeometry) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AfEdgeGeometry {
    /// One of: `"open"`, `"sharp"`, `"square"`, `"rounded_square"`, `"full_round"`,
    /// `"blended_round"`.
    #[getter]
    fn kind(&self) -> &str {
        match &self.inner {
            InnerAfEdgeGeometry::Open => "open",
            InnerAfEdgeGeometry::Sharp(_) => "sharp",
            InnerAfEdgeGeometry::Square(_, _) => "square",
            InnerAfEdgeGeometry::RoundedSquare(_, _, _) => "rounded_square",
            InnerAfEdgeGeometry::FullRound(_, _) => "full_round",
            InnerAfEdgeGeometry::BlendedRound(_, _) => "blended_round",
        }
    }

    /// The apex (`"sharp"`) or circle center (`"full_round"`, `"blended_round"`). `None` for all
    /// other variants.
    #[getter]
    fn point(&self) -> Option<Point2> {
        match &self.inner {
            InnerAfEdgeGeometry::Sharp(p) => Some(Point2::from_inner(*p)),
            InnerAfEdgeGeometry::FullRound(p, _) | InnerAfEdgeGeometry::BlendedRound(p, _) => {
                Some(Point2::from_inner(*p))
            }
            _ => None,
        }
    }

    /// First corner point. Populated for `"square"` and `"rounded_square"` only.
    #[getter]
    fn corner0(&self) -> Option<Point2> {
        match &self.inner {
            InnerAfEdgeGeometry::Square(p, _) | InnerAfEdgeGeometry::RoundedSquare(p, _, _) => {
                Some(Point2::from_inner(*p))
            }
            _ => None,
        }
    }

    /// Second corner point. Populated for `"square"` and `"rounded_square"` only.
    #[getter]
    fn corner1(&self) -> Option<Point2> {
        match &self.inner {
            InnerAfEdgeGeometry::Square(_, p) | InnerAfEdgeGeometry::RoundedSquare(_, p, _) => {
                Some(Point2::from_inner(*p))
            }
            _ => None,
        }
    }

    /// Circle or fillet radius. Populated for `"rounded_square"`, `"full_round"`, and
    /// `"blended_round"` only.
    #[getter]
    fn radius(&self) -> Option<f64> {
        match &self.inner {
            InnerAfEdgeGeometry::RoundedSquare(_, _, r) => Some(*r as f64),
            InnerAfEdgeGeometry::FullRound(_, r) | InnerAfEdgeGeometry::BlendedRound(_, r) => {
                Some(*r as f64)
            }
            _ => None,
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            InnerAfEdgeGeometry::Open => "AfEdgeGeometry(open)".to_string(),
            InnerAfEdgeGeometry::Sharp(p) => {
                format!("AfEdgeGeometry(sharp, apex=({}, {}))", p.x, p.y)
            }
            InnerAfEdgeGeometry::Square(p0, p1) => {
                format!(
                    "AfEdgeGeometry(square, corner0=({}, {}), corner1=({}, {}))",
                    p0.x, p0.y, p1.x, p1.y
                )
            }
            InnerAfEdgeGeometry::RoundedSquare(p0, p1, r) => {
                format!(
                    "AfEdgeGeometry(rounded_square, corner0=({}, {}), corner1=({}, {}), r={})",
                    p0.x, p0.y, p1.x, p1.y, r
                )
            }
            InnerAfEdgeGeometry::FullRound(p, r) => {
                format!(
                    "AfEdgeGeometry(full_round, center=({}, {}), r={})",
                    p.x, p.y, r
                )
            }
            InnerAfEdgeGeometry::BlendedRound(p, r) => {
                format!(
                    "AfEdgeGeometry(blended_round, center=({}, {}), r={})",
                    p.x, p.y, r
                )
            }
        }
    }
}

// ================================================================================================
// AfEdge
// ================================================================================================

/// A fitted airfoil edge: the canonical edge location point and its geometric description.
#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Debug)]
pub struct AfEdge {
    inner: InnerAfEdge,
}

impl AfEdge {
    pub fn from_inner(inner: InnerAfEdge) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AfEdge {
    /// The canonical edge location (e.g. the apex for sharp, midpoint for square, outermost
    /// camber-axis point for round variants).
    #[getter]
    fn point(&self) -> Point2 {
        Point2::from_inner(self.inner.point)
    }

    /// The geometric description of the edge shape.
    #[getter]
    fn geometry(&self) -> AfEdgeGeometry {
        AfEdgeGeometry::from_inner(self.inner.geometry)
    }

    fn __repr__(&self) -> String {
        format!(
            "AfEdge(point=({}, {}), geometry={})",
            self.inner.point.x,
            self.inner.point.y,
            AfEdgeGeometry::from_inner(self.inner.geometry).__repr__()
        )
    }
}

// ================================================================================================
// AfEdgeFit
// ================================================================================================

/// The result of an airfoil edge fitting operation.
#[pyclass(module = "engeom.airfoil2")]
pub struct AfEdgeFit {
    inner: InnerAfEdgeFit,
}

impl AfEdgeFit {
    pub fn from_inner(inner: InnerAfEdgeFit) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AfEdgeFit {
    /// The fitted edge point and geometry.
    #[getter]
    fn edge(&self) -> AfEdge {
        AfEdge::from_inner(self.inner.edge)
    }

    /// The point-to-boundary residuals from the fitting optimisation, as a 1-D numpy array.
    #[getter]
    fn residuals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_iter(self.inner.residuals.iter().copied()).into_pyarray(py)
    }

    /// The inscribed circle stack, potentially refined by the fitting step.
    #[getter]
    fn circles(&self) -> Vec<Inscribed> {
        self.inner
            .circles
            .iter()
            .map(|c| Inscribed::from_inner(c.clone()))
            .collect()
    }
}

// ================================================================================================
// Edge fitting functions
// ================================================================================================

fn to_inner_circles(circles: Vec<Inscribed>) -> Vec<InnerInscribed> {
    circles.into_iter().map(|c| c.into_inner()).collect()
}

fn fit_err(e: Box<dyn std::error::Error>) -> PyErr {
    PyValueError::new_err(e.to_string())
}

#[pyfunction]
pub fn fit_square_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_square_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn fit_rounded_square_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_rounded_square_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn fit_sharp_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_sharp_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn fit_full_round_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_full_round_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn fit_blended_round_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_blended_round_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn extract_inscribed_circles(section: &Curve2, tol: f64) -> PyResult<Vec<Inscribed>> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::camber::extract_inscribed_circles(&input)
        .map(|v| v.into_iter().map(Inscribed::from_inner).collect())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
