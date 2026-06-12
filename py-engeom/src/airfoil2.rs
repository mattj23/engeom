use crate::geom2::{Circle2, Point2, SurfacePoint2, Vector2};
use engeom::airfoil2::SectionInput;
use engeom::airfoil2::inscribed::Inscribed as InnerInscribed;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::geom2::Curve2;

#[pyclass(from_py_object, module = "engeom.airfoil2")]
#[derive(Clone, Debug)]
pub struct Inscribed {
    inner: InnerInscribed,
}

impl Inscribed {
    pub fn from_inner(inner: InnerInscribed) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Inscribed {
    #[getter]
    fn c(&self) -> Circle2 {
        Circle2::from_inner(self.inner.c)
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

/// Extract the unambiguous inscribed circles of an airfoil section.
///
/// Circles are returned in a consistent order, but may run from leading-to-trailing or the
/// reverse. The `p0`/`p1` contact points are consistently oriented but their upper/lower
/// assignment is also ambiguous until the section is oriented.
#[pyfunction]
pub fn extract_inscribed_circles(section: &Curve2, tol: f64) -> PyResult<Vec<Inscribed>> {
    let input = SectionInput::new(section.get_inner(), tol);
    engeom::airfoil2::camber::extract_inscribed_circles(&input)
        .map(|v| v.into_iter().map(Inscribed::from_inner).collect())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
