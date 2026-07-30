use crate::geom2::{Curve2, Point2, SurfacePoint2, Vector2};
use engeom::airfoil2::SectionInput;
use engeom::airfoil2::edges::AfEdgeFit as InnerAfEdgeFit;
use engeom::airfoil2::inscribed::Inscribed as InnerInscribed;
use engeom::airfoil2::{
    AfEdge as InnerAfEdge, AfEdgeGeometry as InnerAfEdgeGeometry,
    AfEdgeSearch as InnerAfEdgeSearch, AfGeometry as InnerAfGeometry,
    OrientFwdAft as InnerOrientFwdAft, OrientUpperLower as InnerOrientUpperLower,
};
use numpy::ndarray::{Array1, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArrayDyn};
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
            InnerAfEdgeGeometry::SplineMaxK(_, _) => "spline_max_k",
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
            InnerAfEdgeGeometry::SplineMaxK(p, r) => {
                format!(
                    "AfEdgeGeometry(spline_max_k, center=({}, {}), r={})",
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

    /// The point-to-boundary residuals from the fitting optimization, as a 1-D numpy array.
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
pub fn fit_open_edge(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<AfEdgeFit> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::fit_open_edge(&input, to_inner_circles(circles), at_front)
        .map(AfEdgeFit::from_inner)
        .map_err(fit_err)
}

#[pyfunction]
pub fn is_open_at_end(
    section: &Curve2,
    tol: f64,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> PyResult<bool> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::edges::is_open_at_end(&input, &to_inner_circles(circles), at_front)
        .map_err(fit_err)
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

// ================================================================================================
// OrientFwdAft
// ================================================================================================

/// Selects the method used to identify which end of the camber line is the leading edge.
///
/// Variants:
///
/// - `TmaxFwd()`: the end nearer the largest inscribed circle becomes the leading edge.
/// - `Airflow(x, y)`: airflow direction; the end further upstream becomes the leading edge.
/// - `Fwd(x, y)`: vector pointing toward the leading edge; the end further in this direction
///   becomes the leading edge.
#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Copy, Debug)]
pub enum OrientFwdAft {
    TmaxFwd {},
    Airflow { x: f64, y: f64 },
    Fwd { x: f64, y: f64 },
}

#[pymethods]
impl OrientFwdAft {
    fn __repr__(&self) -> String {
        match self {
            OrientFwdAft::TmaxFwd {} => "OrientFwdAft.TmaxFwd".to_string(),
            OrientFwdAft::Airflow { x, y } => format!("OrientFwdAft.Airflow({}, {})", x, y),
            OrientFwdAft::Fwd { x, y } => format!("OrientFwdAft.Fwd({}, {})", x, y),
        }
    }
}

impl From<OrientFwdAft> for InnerOrientFwdAft {
    fn from(value: OrientFwdAft) -> Self {
        match value {
            OrientFwdAft::TmaxFwd {} => InnerOrientFwdAft::TmaxFwd,
            OrientFwdAft::Airflow { x, y } => {
                InnerOrientFwdAft::Airflow(engeom::Vector2::new(x, y))
            }
            OrientFwdAft::Fwd { x, y } => InnerOrientFwdAft::Fwd(engeom::Vector2::new(x, y)),
        }
    }
}

// ================================================================================================
// OrientUpperLower
// ================================================================================================

/// Selects the method used to identify the upper (suction) and lower (pressure) surfaces.
///
/// Variants:
///
/// - `Curvature()`: detect from camber line curvature; the more concave side becomes lower.
/// - `Upper(x, y)`: the side whose points are more positive along (x, y) becomes upper.
/// - `Lower(x, y)`: the side whose points are more positive along (x, y) becomes lower.
#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Copy, Debug)]
pub enum OrientUpperLower {
    Curvature {},
    Upper { x: f64, y: f64 },
    Lower { x: f64, y: f64 },
}

#[pymethods]
impl OrientUpperLower {
    fn __repr__(&self) -> String {
        match self {
            OrientUpperLower::Curvature {} => "OrientUpperLower.Curvature".to_string(),
            OrientUpperLower::Upper { x, y } => format!("OrientUpperLower.Upper({}, {})", x, y),
            OrientUpperLower::Lower { x, y } => format!("OrientUpperLower.Lower({}, {})", x, y),
        }
    }
}

impl From<OrientUpperLower> for InnerOrientUpperLower {
    fn from(value: OrientUpperLower) -> Self {
        match value {
            OrientUpperLower::Curvature {} => InnerOrientUpperLower::Curvature,
            OrientUpperLower::Upper { x, y } => {
                InnerOrientUpperLower::Upper(engeom::Vector2::new(x, y))
            }
            OrientUpperLower::Lower { x, y } => {
                InnerOrientUpperLower::Lower(engeom::Vector2::new(x, y))
            }
        }
    }
}

// ================================================================================================
// AfEdgeSearch
// ================================================================================================

/// Selects which edge geometry to fit at the leading or trailing edge.
///
/// `Auto` will treat the edge as open if the section is open at that end, and will otherwise try
/// every fittable variant and select the one with the lowest average residual.
#[pyclass(eq, eq_int, from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AfEdgeSearch {
    Auto = 0,
    Open = 1,
    Sharp = 2,
    Square = 3,
    RoundedSquare = 4,
    FullRound = 5,
    BlendedRound = 6,
    SplineMaxK = 7,
}

#[pymethods]
impl AfEdgeSearch {
    fn __repr__(&self) -> String {
        match self {
            AfEdgeSearch::Auto => "AfEdgeSearch.Auto".to_string(),
            AfEdgeSearch::Open => "AfEdgeSearch.Open".to_string(),
            AfEdgeSearch::Sharp => "AfEdgeSearch.Sharp".to_string(),
            AfEdgeSearch::Square => "AfEdgeSearch.Square".to_string(),
            AfEdgeSearch::RoundedSquare => "AfEdgeSearch.RoundedSquare".to_string(),
            AfEdgeSearch::FullRound => "AfEdgeSearch.FullRound".to_string(),
            AfEdgeSearch::BlendedRound => "AfEdgeSearch.BlendedRound".to_string(),
            AfEdgeSearch::SplineMaxK => "AfEdgeSearch.SplineMaxK".to_string(),
        }
    }
}

impl From<AfEdgeSearch> for InnerAfEdgeSearch {
    fn from(value: AfEdgeSearch) -> Self {
        match value {
            AfEdgeSearch::Auto => InnerAfEdgeSearch::Auto,
            AfEdgeSearch::Open => InnerAfEdgeSearch::Open,
            AfEdgeSearch::Sharp => InnerAfEdgeSearch::Sharp,
            AfEdgeSearch::Square => InnerAfEdgeSearch::Square,
            AfEdgeSearch::RoundedSquare => InnerAfEdgeSearch::RoundedSquare,
            AfEdgeSearch::FullRound => InnerAfEdgeSearch::FullRound,
            AfEdgeSearch::BlendedRound => InnerAfEdgeSearch::BlendedRound,
            AfEdgeSearch::SplineMaxK => InnerAfEdgeSearch::SplineMaxK,
        }
    }
}

// ================================================================================================
// AfGeometry
// ================================================================================================

/// The result of a geometric analysis of an airfoil section: the leading and trailing edges,
/// the mean camber line, the upper and lower surfaces, and the inscribed circle stack.
#[pyclass(module = "engeom.airfoil2")]
pub struct AfGeometry {
    inner: InnerAfGeometry,
    camber: Option<Py<Curve2>>,
    upper: Option<Py<Curve2>>,
    lower: Option<Py<Curve2>>,
    circle_array: Option<Py<PyArrayDyn<f64>>>,
}

impl AfGeometry {
    pub fn from_inner(inner: InnerAfGeometry) -> Self {
        Self {
            inner,
            camber: None,
            upper: None,
            lower: None,
            circle_array: None,
        }
    }
}

#[pymethods]
impl AfGeometry {
    /// Run a purely geometric analysis of an airfoil section.
    ///
    /// Extracts the mean camber line, identifies the leading and trailing edges using the
    /// supplied edge search strategies, orients the section forward/aft and upper/lower, and
    /// splits the section into upper and lower surfaces.
    #[staticmethod]
    fn from_geometric_analysis(
        section: &Curve2,
        general_tol: f64,
        fwd_aft: OrientFwdAft,
        upper_lower: OrientUpperLower,
        le_search: AfEdgeSearch,
        te_search: AfEdgeSearch,
    ) -> PyResult<Self> {
        let inner = InnerAfGeometry::try_from_geometric_analysis(
            section.get_inner(),
            general_tol,
            fwd_aft.into(),
            upper_lower.into(),
            le_search.into(),
            te_search.into(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(AfGeometry::from_inner(inner))
    }

    /// The leading edge result.
    #[getter]
    fn leading(&self) -> AfEdge {
        AfEdge::from_inner(self.inner.leading)
    }

    /// The trailing edge result.
    #[getter]
    fn trailing(&self) -> AfEdge {
        AfEdge::from_inner(self.inner.trailing)
    }

    /// The mean camber line, oriented from leading edge to trailing edge.
    #[getter]
    fn camber<'py>(&mut self, py: Python<'py>) -> &Bound<'py, Curve2> {
        if self.camber.is_none() {
            let c = Curve2::from_inner(self.inner.camber.clone());
            self.camber = Some(Py::new(py, c).unwrap());
        }
        self.camber.as_ref().unwrap().bind(py)
    }

    /// The upper (suction) surface curve.
    #[getter]
    fn upper<'py>(&mut self, py: Python<'py>) -> &Bound<'py, Curve2> {
        if self.upper.is_none() {
            let u = Curve2::from_inner(self.inner.upper.clone());
            self.upper = Some(Py::new(py, u).unwrap());
        }
        self.upper.as_ref().unwrap().bind(py)
    }

    /// The lower (pressure) surface curve.
    #[getter]
    fn lower<'py>(&mut self, py: Python<'py>) -> &Bound<'py, Curve2> {
        if self.lower.is_none() {
            let l = Curve2::from_inner(self.inner.lower.clone());
            self.lower = Some(Py::new(py, l).unwrap());
        }
        self.lower.as_ref().unwrap().bind(py)
    }

    /// The inscribed circle stack, oriented leading-to-trailing with `p0` on the lower surface
    /// and `p1` on the upper surface.
    #[getter]
    fn circles(&self) -> Vec<Inscribed> {
        self.inner
            .circles
            .iter()
            .map(|c| Inscribed::from_inner(c.clone()))
            .collect()
    }

    /// The inscribed circle with the largest radius, corresponding to the maximum thickness
    /// location along the camber line.
    fn max_thickness_circle(&self) -> Inscribed {
        Inscribed::from_inner(self.inner.max_thickness_circle().clone())
    }

    /// The inscribed circle stack as a `(N, 3)` numpy array of `(center_x, center_y, radius)`.
    #[getter]
    fn circle_array<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArrayDyn<f64>> {
        if self.circle_array.is_none() {
            let mut result = ArrayD::zeros(vec![self.inner.circles.len(), 3]);
            for (i, c) in self.inner.circles.iter().enumerate() {
                result[[i, 0]] = c.c.center.x;
                result[[i, 1]] = c.c.center.y;
                result[[i, 2]] = c.c.r();
            }
            self.circle_array = Some(result.into_pyarray(py).unbind());
        }
        self.circle_array.as_ref().unwrap().bind(py)
    }
}
