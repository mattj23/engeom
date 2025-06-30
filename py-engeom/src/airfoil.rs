use crate::geom2::{Arc2, Circle2, Curve2, Point2};
use crate::metrology::Distance2;
use numpy::ndarray::ArrayD;
use numpy::{IntoPyArray, PyArrayDyn};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rmp_serde::{from_slice, to_vec_named};

// ================================================================================================
// Orientation Methods
// ================================================================================================
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum FaceOrient {
    Detect {},
    UpperDir { x: f64, y: f64 },
}

#[pymethods]
impl FaceOrient {
    fn __repr__(&self) -> String {
        match self {
            FaceOrient::Detect {} => "FaceOrient.Detect".to_string(),
            FaceOrient::UpperDir { x, y } => format!("FaceOrient.UpperDir({}, {})", x, y),
        }
    }
}

impl From<FaceOrient> for engeom::airfoil::FaceOrient {
    fn from(value: FaceOrient) -> Self {
        match value {
            FaceOrient::Detect {} => engeom::airfoil::FaceOrient::Detect,
            FaceOrient::UpperDir { x, y } => {
                engeom::airfoil::FaceOrient::UpperDir(engeom::Vector2::new(x, y))
            }
        }
    }
}

#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum AfGage {
    OnCamber { d: f64 },
    Radius { r: f64 },
}

#[pymethods]
impl AfGage {
    fn __repr__(&self) -> String {
        match self {
            AfGage::OnCamber { d } => format!("AfGage.OnCamber({})", d),
            AfGage::Radius { r } => format!("AfGage.Radius({})", r),
        }
    }
}

impl From<AfGage> for engeom::airfoil::AfGage {
    fn from(value: AfGage) -> Self {
        match value {
            AfGage::OnCamber { d } => engeom::airfoil::AfGage::OnCamber(d),
            AfGage::Radius { r } => engeom::airfoil::AfGage::Radius(r),
        }
    }
}

#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum MclOrient {
    TmaxFwd {},
    DirFwd { x: f64, y: f64 },
}

#[pymethods]
impl MclOrient {
    fn __repr__(&self) -> String {
        match self {
            MclOrient::TmaxFwd {} => "MclOrient.TmaxFwd".to_string(),
            MclOrient::DirFwd { x, y } => format!("MclOrient.DirFwd({}, {})", x, y),
        }
    }
}

impl From<MclOrient> for Box<dyn engeom::airfoil::CamberOrient> {
    fn from(value: MclOrient) -> Self {
        match value {
            MclOrient::TmaxFwd {} => engeom::airfoil::TMaxFwd::make(),
            MclOrient::DirFwd { x, y } => {
                engeom::airfoil::DirectionFwd::make(engeom::Vector2::new(x, y))
            }
        }
    }
}

// ================================================================================================
// Edge Extraction Methods
// ================================================================================================
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub enum EdgeFind {
    Open {},
    OpenIntersect { max_iter: usize },
    Intersect {},
    RansacRadius { in_tol: f64, n: usize },
}

#[pymethods]
impl EdgeFind {
    fn __repr__(&self) -> String {
        match self {
            EdgeFind::Open {} => "EdgeFind.Open".to_string(),
            EdgeFind::OpenIntersect { max_iter } => format!("EdgeFind.OpenIntersect({})", max_iter),
            EdgeFind::Intersect {} => "EdgeFind.Intersect".to_string(),
            EdgeFind::RansacRadius { in_tol, n } => {
                format!("EdgeFind.RansacRadius({}, {})", in_tol, n)
            }
        }
    }
}

impl From<EdgeFind> for Box<dyn engeom::airfoil::EdgeLocate> {
    fn from(value: EdgeFind) -> Self {
        use engeom::airfoil;

        match value {
            EdgeFind::Open {} => airfoil::OpenEdge::make(),
            EdgeFind::OpenIntersect { max_iter } => airfoil::OpenIntersectGap::make(max_iter),
            EdgeFind::Intersect {} => airfoil::IntersectEdge::make(),
            EdgeFind::RansacRadius { in_tol, n } => airfoil::RansacRadiusEdge::make(in_tol, n),
        }
    }
}

// ================================================================================================
// Inscribed Circle
// ================================================================================================

#[pyclass]
#[derive(Clone)]
pub struct InscribedCircle {
    inner: engeom::airfoil::InscribedCircle,
}

impl InscribedCircle {
    pub fn get_inner(&self) -> &engeom::airfoil::InscribedCircle {
        &self.inner
    }

    pub fn from_inner(inner: engeom::airfoil::InscribedCircle) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl InscribedCircle {
    #[getter]
    fn circle(&self) -> Circle2 {
        Circle2::from_inner(self.inner.circle)
    }

    #[getter]
    fn contact_a(&self) -> Point2 {
        Point2::from_inner(self.inner.contact_pos)
    }

    #[getter]
    fn contact_b(&self) -> Point2 {
        Point2::from_inner(self.inner.contact_neg)
    }
}

// ================================================================================================
// Airfoil geometry result
// ================================================================================================

#[pyclass(eq, eq_int)]
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum EdgeType {
    Open,
    Closed,
}

#[pymethods]
impl EdgeType {
    fn __repr__(&self) -> String {
        match self {
            EdgeType::Open => "EdgeType.Open".to_string(),
            EdgeType::Closed => "EdgeType.Closed".to_string(),
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct EdgeResult {
    inner: engeom::airfoil::AirfoilEdge,
}

impl EdgeResult {
    pub fn get_inner(&self) -> &engeom::airfoil::AirfoilEdge {
        &self.inner
    }

    pub fn from_inner(inner: engeom::airfoil::AirfoilEdge) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl EdgeResult {
    #[getter]
    fn point(&self) -> Point2 {
        Point2::from_inner(self.inner.point)
    }

    #[getter]
    fn geometry<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        use engeom::airfoil::EdgeGeometry;

        match self.inner.geometry {
            EdgeGeometry::Open => EdgeType::Open.into_bound_py_any(py),
            EdgeGeometry::Closed => EdgeType::Closed.into_bound_py_any(py),
            EdgeGeometry::Arc(a) => Arc2::from_inner(a).into_bound_py_any(py),
        }
    }
}

#[pyclass]
pub struct AirfoilGeometry {
    inner: engeom::airfoil::AirfoilGeometry,

    leading: Option<EdgeResult>,
    trailing: Option<EdgeResult>,

    camber: Option<Py<Curve2>>,
    upper: Option<Py<Curve2>>,
    lower: Option<Py<Curve2>>,

    sides_failed: bool,
    circle_array: Option<Py<PyArrayDyn<f64>>>,
}

impl AirfoilGeometry {
    pub fn get_inner(&self) -> &engeom::airfoil::AirfoilGeometry {
        &self.inner
    }

    pub fn from_inner(inner: engeom::airfoil::AirfoilGeometry) -> Self {
        let leading = inner
            .leading_edge
            .as_ref()
            .map(|e| EdgeResult::from_inner(e.clone()));
        let trailing = inner
            .trailing_edge
            .as_ref()
            .map(|e| EdgeResult::from_inner(e.clone()));
        Self {
            inner,
            leading,
            trailing,
            camber: None,
            upper: None,
            lower: None,
            sides_failed: false,
            circle_array: None,
        }
    }

    pub fn build_sides(&mut self, py: Python) {
        if self.sides_failed || (self.upper.is_some() && self.lower.is_some()) {
            return;
        }

        if let (Some(u), Some(l)) = (self.inner.upper.as_ref(), self.inner.lower.as_ref()) {
            self.upper = Some(Py::new(py, Curve2::from_inner(u.clone())).unwrap());
            self.lower = Some(Py::new(py, Curve2::from_inner(l.clone())).unwrap());
        } else {
            self.sides_failed = true;
        }
    }
}

impl Clone for AirfoilGeometry {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl AirfoilGeometry {
    fn to_bytes(&self) -> PyResult<Vec<u8>> {
        let bytes = to_vec_named(self.get_inner())
            .map_err(|e| PyValueError::new_err(format!("Failed to serialize: {}", e)))?;

        Ok(bytes)
    }

    #[staticmethod]
    fn from_bytes(bytes: &[u8]) -> PyResult<Self> {
        let inner = from_slice(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(AirfoilGeometry::from_inner(inner))
    }

    #[getter]
    fn camber<'py>(&mut self, py: Python<'py>) -> &Bound<'py, Curve2> {
        if self.camber.is_none() {
            let camber = Curve2::from_inner(self.inner.camber.clone());
            self.camber = Some(Py::new(py, camber).unwrap());
        }
        self.camber.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn leading(&self) -> Option<EdgeResult> {
        self.leading.clone()
    }

    #[getter]
    fn trailing(&self) -> Option<EdgeResult> {
        self.trailing.clone()
    }

    #[getter]
    fn upper<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, Curve2>> {
        self.build_sides(py);
        self.upper.as_ref().map(|u| u.bind(py))
    }

    #[getter]
    fn lower<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, Curve2>> {
        self.build_sides(py);
        self.lower.as_ref().map(|l| l.bind(py))
    }

    #[getter]
    fn circle_array<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArrayDyn<f64>> {
        if self.circle_array.is_none() {
            let mut result = ArrayD::zeros(vec![self.inner.stations.len(), 3]);
            for (i, c) in self.inner.stations.iter().enumerate() {
                result[[i, 0]] = c.circle.center.x;
                result[[i, 1]] = c.circle.center.y;
                result[[i, 2]] = c.circle.r();
            }
            self.circle_array = Some(result.into_pyarray(py).unbind());
        }
        self.circle_array.as_ref().unwrap().bind(py)
    }

    #[staticmethod]
    fn from_analyze(
        section: &Curve2,
        refine_tol: f64,
        camber_orient: MclOrient,
        leading: EdgeFind,
        trailing: EdgeFind,
        face_orient: FaceOrient,
    ) -> PyResult<Self> {
        let result = engeom::airfoil::AirfoilGeometry::try_analyze(
            section.get_inner(),
            refine_tol,
            camber_orient.into(),
            leading.into(),
            trailing.into(),
            face_orient.into(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(AirfoilGeometry::from_inner(result))
    }

    fn get_tmax_circle(&self) -> Circle2 {
        Circle2::from_inner(self.inner.find_tmax().circle)
    }

    fn get_thickness(&self, gage: AfGage) -> PyResult<Distance2> {
        Ok(Distance2::from_inner(
            self.inner
                .get_thickness(gage.into())
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        ))
    }

    fn get_tmax(&self) -> PyResult<Distance2> {
        Ok(Distance2::from_inner(
            self.inner
                .get_thickness_max()
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        ))
    }
}

// ================================================================================================
// Functions
// ================================================================================================
#[pyfunction]
pub fn compute_inscribed_circles(
    section: Curve2,
    refine_tol: f64,
) -> PyResult<Vec<InscribedCircle>> {
    let sec = section.get_inner();
    let hull = sec
        .make_hull()
        .ok_or(PyValueError::new_err("Failed to make convex hull"))?;

    let circles = engeom::airfoil::extract_camber_line(sec, &hull, Some(refine_tol))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let result = circles
        .into_iter()
        .map(InscribedCircle::from_inner)
        .collect();

    Ok(result)
}
