use crate::conversions::array_to_points3;
use crate::geom3::{Iso3, Point3};
use crate::mesh::Mesh;
use engeom::geom3::align3::{AlignOrigin, Dof6 as InnerDof6};
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::{Bound, PyResult, Python, pyclass, pyfunction, pymethods};

// ================================================================================================
// Dof6
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.align3")]
#[derive(Clone, Copy, Debug)]
pub struct Dof6 {
    #[pyo3(get, set)]
    pub tx: bool,
    #[pyo3(get, set)]
    pub ty: bool,
    #[pyo3(get, set)]
    pub tz: bool,
    #[pyo3(get, set)]
    pub rx: bool,
    #[pyo3(get, set)]
    pub ry: bool,
    #[pyo3(get, set)]
    pub rz: bool,
}

#[pymethods]
impl Dof6 {
    #[new]
    #[pyo3(signature = (tx=true, ty=true, tz=true, rx=true, ry=true, rz=true))]
    pub fn new(tx: bool, ty: bool, tz: bool, rx: bool, ry: bool, rz: bool) -> Self {
        Self {
            tx,
            ty,
            tz,
            rx,
            ry,
            rz,
        }
    }

    #[staticmethod]
    pub fn all() -> Self {
        Self {
            tx: true,
            ty: true,
            tz: true,
            rx: true,
            ry: true,
            rz: true,
        }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Dof6(tx={}, ty={}, tz={}, rx={}, ry={}, rz={})",
            self.tx, self.ty, self.tz, self.rx, self.ry, self.rz
        )
    }
}

impl From<Dof6> for InnerDof6 {
    fn from(val: Dof6) -> Self {
        InnerDof6::new(val.tx, val.ty, val.tz, val.rx, val.ry, val.rz)
    }
}

// ================================================================================================
// AlignParams3
// ================================================================================================
#[pyclass(from_py_object, module = "engeom.align3")]
#[derive(Clone, Debug)]
pub struct AlignParams3 {
    inner: engeom::geom3::align3::AlignParams3,
}

impl AlignParams3 {
    pub fn from_inner(inner: engeom::geom3::align3::AlignParams3) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::geom3::align3::AlignParams3 {
        &self.inner
    }
}

#[pymethods]
impl AlignParams3 {
    /// Create an `AlignParams3` describing how a 3D alignment is parameterized.
    ///
    /// The local origin $L$ is selected by supplying at most one of `center` or `local`:
    ///
    /// - `center`: rotations happen about this point, and translations act along the world axes.
    /// - `local`: rotations happen about, and translations act along, the axes of this full
    ///   `Iso3` frame. Use this when you want full control over the rotation center and the
    ///   directions of translation, for example when applying DOF constraints along an arbitrary
    ///   direction.
    /// - neither: the world origin is used. Use this when the test geometry is already close to
    ///   the origin and numerical stability of rotations is not a concern.
    ///
    /// Supplying both `center` and `local` raises a `ValueError`.
    ///
    /// If `offset` is not given, it defaults to the local origin frame, so the test geometry
    /// starts in place and the alignment happens about that origin. Only pass an explicit
    /// `offset` (including the identity) if you specifically need the raw `O * A * L^-1`
    /// behavior where the geometry is displaced by `L^-1` before alignment.
    ///
    /// :param center: Optional `Point3` rotation center. Mutually exclusive with `local`.
    /// :param local: Optional `Iso3` local origin frame. Mutually exclusive with `center`.
    /// :param offset: Optional `Iso3` working offset $O$. Defaults to the local origin frame.
    /// :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
    #[new]
    #[pyo3(signature = (center=None, local=None, offset=None, dof=None))]
    pub fn new(
        center: Option<&Point3>,
        local: Option<&Iso3>,
        offset: Option<&Iso3>,
        dof: Option<Dof6>,
    ) -> PyResult<Self> {
        let (origin, frame) = match (center, local) {
            (Some(_), Some(_)) => {
                return Err(PyValueError::new_err(
                    "Supply at most one of `center` or `local`, not both",
                ));
            }
            (Some(p), None) => {
                let p = *p.get_inner();
                (
                    AlignOrigin::Center(p),
                    engeom::Iso3::translation(p.x, p.y, p.z),
                )
            }
            (None, Some(o)) => {
                let t = *o.get_inner();
                (AlignOrigin::Local(t), t)
            }
            (None, None) => (AlignOrigin::Origin, engeom::Iso3::identity()),
        };

        let offset = offset.map(|o| *o.get_inner()).unwrap_or(frame);

        Ok(Self::from_inner(engeom::geom3::align3::AlignParams3::new(
            origin,
            Some(offset),
            dof.map(Into::into),
        )))
    }

    /// The degrees-of-freedom constraint currently active on this alignment.
    #[getter]
    pub fn dof(&self) -> Dof6 {
        let d = self.inner.dof;
        Dof6 {
            tx: d.tx,
            ty: d.ty,
            tz: d.tz,
            rx: d.rx,
            ry: d.ry,
            rz: d.rz,
        }
    }

    /// The local origin transformation $L$.
    #[getter]
    pub fn local(&self) -> Iso3 {
        Iso3::from_inner(self.inner.local)
    }

    /// The working offset transformation $O$.
    #[getter]
    pub fn offset(&self) -> Iso3 {
        Iso3::from_inner(self.inner.offset)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AlignParams3(dof={:?}, local={:?}, offset={:?})",
            self.inner.dof, self.inner.local, self.inner.offset
        )
    }
}

// ================================================================================================
// Alignment3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.align3")]
#[derive(Clone, Debug)]
pub struct Alignment3 {
    inner: engeom::geom3::Alignment3,
}

impl Alignment3 {
    pub fn from_inner(inner: engeom::geom3::Alignment3) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::geom3::Alignment3 {
        &self.inner
    }
}

#[pymethods]
impl Alignment3 {
    /// The full transformation from the test entity's coordinate system to the target's coordinate
    /// system. This is the composite $O * A * L^{-1}$ and is what you typically apply to the test
    /// geometry after alignment completes.
    #[getter]
    pub fn full(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.full())
    }

    /// The alignment transformation $A$, which is the transformation produced by the six
    /// optimized parameters (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`) about the local origin.
    #[getter]
    pub fn alignment(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.alignment())
    }

    /// The local origin transformation $L$ that was used during alignment.
    #[getter]
    pub fn local_origin(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.local_origin())
    }

    /// The working offset transformation $O$ that was used during alignment.
    #[getter]
    pub fn offset(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.offset())
    }

    /// The per-sample residuals from the alignment, as a 1-D numpy array of `float64` values.
    /// Residuals are signed distances between each sampled point and the target surface after
    /// the alignment transformation is applied.
    pub fn residuals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.inner.residuals().to_vec()).into_pyarray(py)
    }

    /// The mean of the residuals.
    pub fn residual_mean(&self) -> f64 {
        self.inner.residual_mean()
    }

    /// The mean and standard deviation of the residuals as a `(mean, std_dev)` tuple.
    pub fn residual_mean_std_dev(&self) -> (f64, f64) {
        self.inner.residual_mean_std_dev()
    }

    pub fn __repr__(&self) -> String {
        let (mean, std) = self.inner.residual_mean_std_dev();
        format!(
            "Alignment3(residual_mean={:.6}, residual_std_dev={:.6})",
            mean, std
        )
    }
}

// ================================================================================================
// Functions
// ================================================================================================

#[pyfunction]
#[pyo3(signature = (points, mesh, params))]
pub fn points_to_mesh(
    points: PyReadonlyArray2<'_, f64>,
    mesh: &Mesh,
    params: AlignParams3,
) -> PyResult<Alignment3> {
    let points = array_to_points3(&points.as_array())?;

    let result = engeom::geom3::align3::points_to_surface3(
        &points,
        mesh.get_inner(),
        params.get_inner().clone(),
        false,
        true,
    );

    match result {
        Ok(align) => Ok(Alignment3::from_inner(align)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
