use crate::conversions::array_to_points3;
use crate::geom3::{Iso3, Point3};
use crate::mesh::Mesh3;
use crate::point_cloud::PointCloud3;
use engeom::geom3::align3::{AlignOrigin3, Dof6 as InnerDof6};
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
                    AlignOrigin3::Center(p),
                    engeom::Iso3::translation(p.x, p.y, p.z),
                )
            }
            (None, Some(o)) => {
                let t = *o.get_inner();
                (AlignOrigin3::Local(t), t)
            }
            (None, None) => (AlignOrigin3::Origin, engeom::Iso3::identity()),
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
    pub fn full_transform(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.full_transform())
    }

    /// The alignment transformation $A$, which is the motion produced by the six optimized
    /// parameters (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`) expressed in the frame of the local origin.
    ///
    /// This is not the transformation to apply to the test geometry; use `full_transform` for
    /// that. Reading $O * A * L^{-1}$ right to left, $L^{-1}$ puts a point into the local origin's
    /// frame, $A$ moves it while it is there, and $O$ maps it back out, so $A$ is only meaningful
    /// applied to local-frame coordinates.
    #[getter]
    pub fn local_transform(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.local_transform())
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
// AlignOutcome3
// ================================================================================================

/// How a single Levenberg-Marquardt solve ended, classified by whether its result can be used.
///
/// `"converged"` means a convergence criterion was met. `"unconverged"` means the solver ran out
/// of its evaluation budget, so the parameters are the best it found but convergence was never
/// demonstrated; the alignment is still valid geometry. A solve that broke down entirely is never
/// reported here, because its result is discarded rather than returned.
fn quality_str(q: engeom::common::SolveQuality) -> &'static str {
    match q {
        engeom::common::SolveQuality::Converged => "converged",
        engeom::common::SolveQuality::Unconverged => "unconverged",
        engeom::common::SolveQuality::Failed => "failed",
    }
}

/// How a single Levenberg-Marquardt solve terminated, as a stable snake_case string rather than
/// the solver crate's `Debug` formatting.
fn termination_str(t: &engeom::common::TerminationReason) -> String {
    use engeom::common::TerminationReason as T;
    match t {
        T::Converged { ftol, xtol } => match (ftol, xtol) {
            (true, true) => "converged(ftol,xtol)".to_string(),
            (true, false) => "converged(ftol)".to_string(),
            (false, true) => "converged(xtol)".to_string(),
            (false, false) => "converged".to_string(),
        },
        T::ResidualsZero => "residuals_zero".to_string(),
        T::Orthogonal => "orthogonal".to_string(),
        T::LostPatience => "lost_patience".to_string(),
        T::Numerical(s) => format!("numerical({s})"),
        T::User(s) => format!("user({s})"),
        T::NoImprovementPossible(s) => format!("no_improvement_possible({s})"),
        T::NoParameters => "no_parameters".to_string(),
        T::NoResiduals => "no_residuals".to_string(),
        T::WrongDimensions(s) => format!("wrong_dimensions({s})"),
    }
}

/// Why robust refinement stopped before completing every requested round.
fn halt_str(h: &engeom::common::RefinementHalt) -> String {
    match h {
        engeom::common::RefinementHalt::NoNoiseEstimate => "no_noise_estimate".to_string(),
        engeom::common::RefinementHalt::Underdetermined { weighted, params } => {
            format!("underdetermined({weighted} weighted points, {params} parameters)")
        }
        engeom::common::RefinementHalt::SolveFailed(t) => {
            format!("solve_failed({})", termination_str(t))
        }
    }
}

/// The full outcome of a 3-D alignment: the `Alignment3` itself, plus a record of how the solves
/// which produced it terminated.
///
/// This is only ever returned, never accepted as an argument, so unlike the other classes here it
/// does not implement `from_py_object`. It cannot: the termination reasons it carries come from
/// the underlying solver crate and are not `Clone`.
#[pyclass(module = "engeom.align3")]
#[derive(Debug)]
pub struct AlignOutcome3 {
    inner: engeom::geom3::AlignOutcome3,
}

impl AlignOutcome3 {
    pub fn from_inner(inner: engeom::geom3::AlignOutcome3) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AlignOutcome3 {
    /// The alignment which was produced.
    #[getter]
    pub fn alignment(&self) -> Alignment3 {
        Alignment3::from_inner(self.inner.alignment().clone())
    }

    /// The quality of the weakest solve that contributed to the result, as `"converged"` or
    /// `"unconverged"`. See the module documentation for why an unconverged result is still usable.
    #[getter]
    pub fn quality(&self) -> &'static str {
        quality_str(self.inner.quality())
    }

    /// Whether every solve that contributed to the result met a convergence criterion.
    #[getter]
    pub fn converged(&self) -> bool {
        self.inner.converged()
    }

    /// The number of robust refinement rounds which completed and contributed to the result.
    #[getter]
    pub fn refinement_rounds(&self) -> usize {
        self.inner.refinement_rounds()
    }

    /// How each contributing solve terminated, as a list of strings, beginning with the initial
    /// solve and followed by one entry per completed refinement round.
    #[getter]
    pub fn solves(&self) -> Vec<String> {
        self.inner.solves().iter().map(termination_str).collect()
    }

    /// Why robust refinement stopped early, or `None` if it ran every round it was asked to.
    #[getter]
    pub fn halt(&self) -> Option<String> {
        self.inner.halt().map(halt_str)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AlignOutcome3(quality={}, refinement_rounds={}, halt={:?})",
            self.quality(),
            self.refinement_rounds(),
            self.halt()
        )
    }
}

// ================================================================================================
// Functions
// ================================================================================================

/// Align a set of 3-D points to a mesh by repeatedly projecting them onto their closest position
/// on the surface as the solver moves them.
///
/// By default this is a robust alignment: an initial unweighted solve followed by
/// `refinement_steps` rounds of iteratively reweighted least squares using MAGSAC++ weights. Pass
/// `refinement_steps=0` for a plain unweighted least-squares alignment.
///
/// A `ValueError` is raised only when there is no answer at all: the arguments were rejected, or
/// the initial solve broke down. A solve that merely exhausts its evaluation budget returns
/// normally with `quality == "unconverged"` on the outcome.
///
/// :param points: an `(n, 3)` array of the points to align, in their own coordinate system.
/// :param mesh: the stationary `Mesh3` to align to. If it carries per-vertex `point_stdev`, that
///     uncertainty is used automatically, interpolated to each match position.
/// :param params: an `AlignParams3` describing how the alignment is parameterized.
/// :param ignore_off_target: weight points at 0.0 when they do not project onto the surface.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
/// :param point_sigma: optional per-point standard deviations, one per input point. Combines in
///     quadrature with any uncertainty the mesh reports.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
#[pyfunction]
#[pyo3(signature = (points, mesh, params, ignore_off_target=false, refinement_steps=4,
                    sigma_max=None, point_sigma=None, patience=100))]
#[allow(clippy::too_many_arguments)]
pub fn points_to_mesh(
    points: PyReadonlyArray2<'_, f64>,
    mesh: &Mesh3,
    params: AlignParams3,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome3> {
    let points = array_to_points3(&points.as_array())?;

    let opts = engeom::geom3::align3::AlignOptions3 {
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma: point_sigma.as_deref(),
        patience,
    };

    let result = engeom::geom3::align3::points_to_surface3(
        &points,
        mesh.get_inner(),
        params.get_inner().clone(),
        &opts,
    );

    match result {
        Ok(outcome) => Ok(AlignOutcome3::from_inner(outcome)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

/// Align a set of 3-D points to a point cloud, by repeatedly projecting them onto the tangent plane
/// at their nearest neighbour as the solver moves them.
///
/// This is the point-cloud counterpart of `points_to_mesh` and takes the same options, but a cloud
/// is only samples of a surface rather than a surface, which brings two differences worth knowing.
///
/// The match is the query projected onto the tangent plane at the nearest cloud point, not that
/// point itself. Matching to the nearest point would leave a residual floor of roughly half the
/// sample spacing even at a perfect pose, since a test point between samples can get no closer than
/// the gap between them. On an engine blade sampled at 2 mm this is the difference between
/// recovering a pose to 4.7e-6 mm and to 2.1e-2 mm.
///
/// Because of that, the cloud **must carry per-point normals**: a normal supplies both the tangent
/// plane and the sign of the residual, and neither can be recovered from positions alone. Use
/// `PointCloud3.estimate_normals` if the cloud does not already have them.
///
/// A `ValueError` is raised only when there is no answer at all: the arguments were rejected, the
/// cloud has no normals, or the initial solve broke down. A solve that merely exhausts its
/// evaluation budget returns normally with `quality == "unconverged"` on the outcome.
///
/// :param points: an `(n, 3)` array of the points to align, in their own coordinate system.
/// :param cloud: the stationary `PointCloud3` to align to, which must carry normals. If it carries
///     `point_stdev` that uncertainty is used automatically, and if it carries `voxel_coherence`
///     from `reduce_by_voxel` that becomes the per-match weight, so voxels which straddled an edge
///     count for less.
/// :param params: an `AlignParams3` describing how the alignment is parameterized.
/// :param max_extrapolation: how far *laterally* a point may sit from the nearest cloud point and
///     still count as on-surface. This bounds the in-plane distance only, never the distance along
///     the normal, because that component is the misalignment being solved for. It exists to catch
///     points which have wandered past the edge of the cloud, where the tangent plane is an
///     extrapolation into space nothing was measured in. Set it at a small multiple of the sample
///     spacing. Points beyond it are still matched, but are only discarded if `ignore_off_target`
///     is set.
/// :param ignore_off_target: weight points at 0.0 when they fall beyond `max_extrapolation`.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
/// :param point_sigma: optional per-point standard deviations, one per input point. Combines in
///     quadrature with any uncertainty the cloud reports.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
#[pyfunction]
#[pyo3(signature = (points, cloud, params, max_extrapolation, ignore_off_target=false,
                    refinement_steps=4, sigma_max=None, point_sigma=None, patience=100))]
#[allow(clippy::too_many_arguments)]
pub fn points_to_cloud(
    points: PyReadonlyArray2<'_, f64>,
    cloud: &PointCloud3,
    params: AlignParams3,
    max_extrapolation: f64,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome3> {
    let points = array_to_points3(&points.as_array())?;

    let target = engeom::geom3::align3::CloudTarget3::try_new(cloud.get_inner(), max_extrapolation)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let opts = engeom::geom3::align3::AlignOptions3 {
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma: point_sigma.as_deref(),
        patience,
    };

    let result = engeom::geom3::align3::points_to_surface3(
        &points,
        &target,
        params.get_inner().clone(),
        &opts,
    );

    match result {
        Ok(outcome) => Ok(AlignOutcome3::from_inner(outcome)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}
