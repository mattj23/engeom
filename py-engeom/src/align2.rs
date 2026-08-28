use crate::conversions::{array_to_points2, array_to_unit_vectors2};
use crate::geom2::{Curve2, CurveGroup2, Iso2, Point2};
use crate::solve_report::{halt_str, quality_str, termination_str};
use engeom::geom2::align2::{AlignOrigin2, Dof3 as InnerDof3};
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::{Bound, PyResult, Python, pyclass, pyfunction, pymethods};

// ================================================================================================
// Dof3
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.align2")]
#[derive(Clone, Copy, Debug)]
pub struct Dof3 {
    #[pyo3(get, set)]
    pub tx: bool,
    #[pyo3(get, set)]
    pub ty: bool,
    #[pyo3(get, set)]
    pub rz: bool,
}

#[pymethods]
impl Dof3 {
    #[new]
    #[pyo3(signature = (tx=true, ty=true, rz=true))]
    pub fn new(tx: bool, ty: bool, rz: bool) -> Self {
        Self { tx, ty, rz }
    }

    #[staticmethod]
    pub fn all() -> Self {
        Self {
            tx: true,
            ty: true,
            rz: true,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("Dof3(tx={}, ty={}, rz={})", self.tx, self.ty, self.rz)
    }
}

impl From<Dof3> for InnerDof3 {
    fn from(val: Dof3) -> Self {
        InnerDof3::new(val.tx, val.ty, val.rz)
    }
}

// ================================================================================================
// AlignParams2
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.align2")]
#[derive(Clone, Debug)]
pub struct AlignParams2 {
    inner: engeom::geom2::align2::AlignParams2,
}

impl AlignParams2 {
    pub fn get_inner(&self) -> &engeom::geom2::align2::AlignParams2 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom2::align2::AlignParams2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AlignParams2 {
    /// Create an `AlignParams2` describing how a 2D alignment is parameterized.
    ///
    /// The local origin $L$ is selected by supplying at most one of `center` or `local`:
    ///
    /// - `center`: rotation happens about this point, and translations act along the world axes.
    /// - `local`: rotation happens about, and translations act along, the axes of this full `Iso2`
    ///   frame. Use this when you want full control over the rotation center and the directions of
    ///   translation, for example when applying DOF constraints along an arbitrary direction.
    /// - neither: the world origin is used. Use this when the test geometry is already close to
    ///   the origin and numerical stability of the rotation is not a concern.
    ///
    /// Supplying both `center` and `local` raises a `ValueError`.
    ///
    /// If `offset` is not given, it defaults to the local origin frame, so the test geometry
    /// starts in place and the alignment happens about that origin. Only pass an explicit
    /// `offset` (including the identity) if you specifically need the raw `O * A * L^-1` behavior
    /// where the geometry is displaced by `L^-1` before alignment.
    ///
    /// :param center: Optional `Point2` rotation center. Mutually exclusive with `local`.
    /// :param local: Optional `Iso2` local origin frame. Mutually exclusive with `center`.
    /// :param offset: Optional `Iso2` working offset $O$. Defaults to the local origin frame.
    /// :param dof: Optional `Dof3` constraint. If `None`, all three degrees of freedom are active.
    #[new]
    #[pyo3(signature = (center=None, local=None, offset=None, dof=None))]
    pub fn new(
        center: Option<&Point2>,
        local: Option<&Iso2>,
        offset: Option<&Iso2>,
        dof: Option<Dof3>,
    ) -> PyResult<Self> {
        let (origin, frame) = match (center, local) {
            (Some(_), Some(_)) => {
                return Err(PyValueError::new_err(
                    "Supply at most one of `center` or `local`, not both",
                ));
            }
            (Some(p), None) => {
                let p = *p.get_inner();
                (AlignOrigin2::Center(p), engeom::Iso2::translation(p.x, p.y))
            }
            (None, Some(o)) => {
                let t = *o.get_inner();
                (AlignOrigin2::Local(t), t)
            }
            (None, None) => (AlignOrigin2::Origin, engeom::Iso2::identity()),
        };

        let offset = offset.map(|o| *o.get_inner()).unwrap_or(frame);

        Ok(Self::from_inner(engeom::geom2::align2::AlignParams2::new(
            origin,
            Some(offset),
            dof.map(Into::into),
        )))
    }

    /// The degrees-of-freedom constraint currently active on this alignment.
    #[getter]
    pub fn dof(&self) -> Dof3 {
        let d = self.inner.dof;
        Dof3 {
            tx: d.tx,
            ty: d.ty,
            rz: d.rz,
        }
    }

    /// The local origin transformation $L$.
    #[getter]
    pub fn local(&self) -> Iso2 {
        Iso2::from_inner(self.inner.local)
    }

    /// The working offset transformation $O$.
    #[getter]
    pub fn offset(&self) -> Iso2 {
        Iso2::from_inner(self.inner.offset)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "AlignParams2(dof={:?}, local={:?}, offset={:?})",
            self.inner.dof, self.inner.local, self.inner.offset
        )
    }
}

// ================================================================================================
// Align2
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.align2")]
#[derive(Clone, Debug)]
pub struct Align2 {
    inner: engeom::geom2::Align2,
}

impl Align2 {
    pub fn from_inner(inner: engeom::geom2::Align2) -> Self {
        Self { inner }
    }

    pub fn get_inner(&self) -> &engeom::geom2::Align2 {
        &self.inner
    }
}

#[pymethods]
impl Align2 {
    /// The full transformation from the test entity's coordinate system to the target's coordinate
    /// system. This is the composite $O * A * L^{-1}$ and is what you typically apply to the test
    /// geometry after alignment completes.
    #[getter]
    pub fn full_transform(&self) -> Iso2 {
        Iso2::from_inner(*self.inner.full_transform())
    }

    /// The alignment transformation $A$, which is the motion produced by the three optimized
    /// parameters (`tx`, `ty`, `rz`) expressed in the frame of the local origin.
    ///
    /// This is not the transformation to apply to the test geometry; use `full_transform` for
    /// that. Reading $O * A * L^{-1}$ right to left, $L^{-1}$ puts a point into the local origin's
    /// frame, $A$ moves it while it is there, and $O$ maps it back out, so $A$ is only meaningful
    /// applied to local-frame coordinates.
    #[getter]
    pub fn local_transform(&self) -> Iso2 {
        Iso2::from_inner(*self.inner.local_transform())
    }

    /// The local origin transformation $L$ that was used during alignment.
    #[getter]
    pub fn local_origin(&self) -> Iso2 {
        Iso2::from_inner(*self.inner.local_origin())
    }

    /// The working offset transformation $O$ that was used during alignment.
    #[getter]
    pub fn offset(&self) -> Iso2 {
        Iso2::from_inner(*self.inner.offset())
    }

    /// The per-sample residuals from the alignment, as a 1-D numpy array of `float64` values.
    /// Residuals are signed distances between each sampled point and the target after the
    /// alignment transformation is applied.
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
            "Align2(residual_mean={:.6}, residual_std_dev={:.6})",
            mean, std
        )
    }
}

// ================================================================================================
// AlignOutcome2
// ================================================================================================

/// The full outcome of a 2-D alignment: the `Align2` itself, plus a record of how the solves
/// which produced it terminated.
///
/// This is only ever returned, never accepted as an argument, so unlike the other classes here it
/// does not implement `from_py_object`. It cannot: the termination reasons it carries come from
/// the underlying solver crate and are not `Clone`.
#[pyclass(module = "engeom.align2")]
#[derive(Debug)]
pub struct AlignOutcome2 {
    inner: engeom::geom2::AlignOutcome2,
}

impl AlignOutcome2 {
    pub fn from_inner(inner: engeom::geom2::AlignOutcome2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl AlignOutcome2 {
    /// The alignment which was produced.
    #[getter]
    pub fn alignment(&self) -> Align2 {
        Align2::from_inner(self.inner.alignment().clone())
    }

    /// The quality of the weakest solve that contributed to the result, as `"converged"` or
    /// `"unconverged"`. An unconverged result is still valid geometry.
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
            "AlignOutcome2(quality={}, refinement_rounds={}, halt={:?})",
            self.quality(),
            self.refinement_rounds(),
            self.halt()
        )
    }
}

// ================================================================================================
// MultiAlignOutcome2
// ================================================================================================

/// The full outcome of a simultaneous alignment of several bodies: one `Align2` per body,
/// plus a record of how the solves which produced them terminated.
///
/// The solve record is shared rather than per body, because a bundle adjustment is one
/// least-squares problem: the bodies converge or fail to converge together.
#[pyclass(module = "engeom.align2")]
#[derive(Debug)]
pub struct MultiAlignOutcome2 {
    inner: engeom::geom2::MultiAlignOutcome2,
}

impl MultiAlignOutcome2 {
    pub fn from_inner(inner: engeom::geom2::MultiAlignOutcome2) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl MultiAlignOutcome2 {
    /// The alignment produced for each body, in the order the bodies were given.
    #[getter]
    pub fn alignments(&self) -> Vec<Align2> {
        self.inner
            .alignments()
            .iter()
            .map(|a| Align2::from_inner(a.clone()))
            .collect()
    }

    /// The alignment produced for one body.
    pub fn alignment(&self, body: usize) -> PyResult<Align2> {
        if body >= self.inner.len() {
            return Err(PyIndexError::new_err(format!(
                "body index {} is out of range for {} bodies",
                body,
                self.inner.len()
            )));
        }
        Ok(Align2::from_inner(self.inner.alignment(body).clone()))
    }

    /// The number of bodies which took part.
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// The quality of the weakest solve that contributed to the result.
    #[getter]
    pub fn quality(&self) -> &'static str {
        quality_str(self.inner.quality())
    }

    /// Whether every solve that contributed to the result met a convergence criterion.
    #[getter]
    pub fn converged(&self) -> bool {
        self.inner.converged()
    }

    /// The number of robust refinement rounds which completed.
    #[getter]
    pub fn refinement_rounds(&self) -> usize {
        self.inner.refinement_rounds()
    }

    /// How each contributing solve terminated.
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
            "MultiAlignOutcome2(bodies={}, quality={}, refinement_rounds={})",
            self.inner.len(),
            self.quality(),
            self.refinement_rounds()
        )
    }
}

// ================================================================================================
// Functions
// ================================================================================================

/// Build the single-body solver options from the loose keyword arguments the bindings take.
fn single_opts(
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<&[f64]>,
    patience: usize,
) -> engeom::geom2::align2::AlignOptions2<'_> {
    engeom::geom2::align2::AlignOptions2 {
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma,
        patience,
    }
}

/// The shared tail of the single-body alignment functions: convert the points, build the
/// options, run the solve, and map the error.
#[allow(clippy::too_many_arguments)]
fn points_to_target(
    points: PyReadonlyArray2<'_, f64>,
    target: &impl engeom::geom2::align2::SurfaceTarget2,
    params: &AlignParams2,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome2> {
    let points = array_to_points2(&points.as_array())?;
    let opts = single_opts(
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma.as_deref(),
        patience,
    );

    engeom::geom2::align2::points_to_surface2(&points, target, params.get_inner().clone(), &opts)
        .map(AlignOutcome2::from_inner)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Align a set of 2-D points to a target by repeatedly projecting them onto their closest position
/// on it as the solver moves them.
///
/// By default this is a robust alignment: an initial unweighted solve followed by
/// `refinement_steps` rounds of iteratively reweighted least squares using MAGSAC++ weights. Pass
/// `refinement_steps=0` for a plain unweighted least-squares alignment.
///
/// A `ValueError` is raised only when there is no answer at all: the arguments were rejected, or
/// the initial solve broke down. A solve that merely exhausts its evaluation budget returns
/// normally with `quality == "unconverged"` on the outcome.
///
/// :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
/// :param curve: the stationary `Curve2` to align to.
/// :param params: an `AlignParams2` describing how the alignment is parameterized.
/// :param ignore_off_target: weight points at 0.0 when they project past an open end of the curve.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
/// :param point_sigma: optional per-point standard deviations, one per input point.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
#[pyfunction]
#[pyo3(signature = (points, curve, params, ignore_off_target=false, refinement_steps=4,
                    sigma_max=None, point_sigma=None, patience=100))]
#[allow(clippy::too_many_arguments)]
pub fn points_to_curve(
    points: PyReadonlyArray2<'_, f64>,
    curve: &Curve2,
    params: AlignParams2,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome2> {
    points_to_target(
        points,
        curve.get_inner(),
        &params,
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma,
        patience,
    )
}

/// Align a set of 2-D points to a `CurveGroup2`, which is a collection of disjoint curves treated
/// as one rigid entity, such as the loops and open segments a planar section of a mesh produces.
///
/// Behaves exactly as `points_to_curve` otherwise: each point matches whichever member of the
/// group is closest, and whether it landed past an open end is judged against that member.
///
/// :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
/// :param group: the stationary `CurveGroup2` to align to.
/// :param params: an `AlignParams2` describing how the alignment is parameterized.
/// :param ignore_off_target: weight points at 0.0 when they project past an open end of a member.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
/// :param point_sigma: optional per-point standard deviations, one per input point.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
#[pyfunction]
#[pyo3(signature = (points, group, params, ignore_off_target=false, refinement_steps=4,
                    sigma_max=None, point_sigma=None, patience=100))]
#[allow(clippy::too_many_arguments)]
pub fn points_to_group(
    points: PyReadonlyArray2<'_, f64>,
    group: &CurveGroup2,
    params: AlignParams2,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome2> {
    points_to_target(
        points,
        group.get_inner(),
        &params,
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma,
        patience,
    )
}

/// Align a set of 2-D points to an unordered set of measured points which carry normals.
///
/// Each match is the query projected onto the tangent line at its nearest neighbor rather than the
/// neighbor itself, which is what lets the alignment resolve below the spacing of the target
/// samples. See `max_extrapolation` for the bound that keeps that projection from being applied
/// where nothing was measured.
///
/// :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
/// :param target_points: an `(m, 2)` array of the stationary measured positions.
/// :param target_normals: an `(m, 2)` array of normals, one per target point. These supply both
///     the tangent line a match lands on and the sign of the residual.
/// :param params: an `AlignParams2` describing how the alignment is parameterized.
/// :param max_extrapolation: how far *along the surface* a query may sit from its nearest target
///     point and still count as on-surface. This is not a total distance: the normal component is
///     excluded, because that component is the residual the alignment exists to remove. Set it at
///     a small multiple of the target's sample spacing.
/// :param target_sigma: optional per-target-point standard deviations, one per target point.
/// :param ignore_off_target: weight points at 0.0 when they sit beyond `max_extrapolation`.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
/// :param point_sigma: optional per-point standard deviations, one per input point. Combines in
///     quadrature with any uncertainty the target reports.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
#[pyfunction]
#[pyo3(signature = (points, target_points, target_normals, params, max_extrapolation,
                    target_sigma=None, ignore_off_target=false, refinement_steps=4,
                    sigma_max=None, point_sigma=None, patience=100))]
#[allow(clippy::too_many_arguments)]
pub fn points_to_cloud(
    points: PyReadonlyArray2<'_, f64>,
    target_points: PyReadonlyArray2<'_, f64>,
    target_normals: PyReadonlyArray2<'_, f64>,
    params: AlignParams2,
    max_extrapolation: f64,
    target_sigma: Option<Vec<f64>>,
    ignore_off_target: bool,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    point_sigma: Option<Vec<f64>>,
    patience: usize,
) -> PyResult<AlignOutcome2> {
    let target_points = array_to_points2(&target_points.as_array())?;
    let normals = array_to_unit_vectors2(&target_normals.as_array())?;

    let target = engeom::geom2::align2::CloudTarget2::try_new(
        &target_points,
        &normals,
        target_sigma.as_deref(),
        max_extrapolation,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    points_to_target(
        points,
        &target,
        &params,
        ignore_off_target,
        refinement_steps,
        sigma_max,
        point_sigma,
        patience,
    )
}

/// Simultaneously align several `CurveGroup2` bodies to each other in one combined solve.
///
/// This is a bundle adjustment rather than a pose graph optimization: one transformation is
/// carried for each body except one, which is held fixed, and all of them are solved for at once
/// against the raw correspondences. It works best on low-noise data which has already been
/// pre-aligned close to each other with substantial overlap. The static body is chosen
/// automatically as the one the others reference most broadly, and correspondences are sampled
/// between every pair.
///
/// :param groups: the bodies taking part. At least two are required.
/// :param max_distance: a hard cutoff on correspondence distance, in the units of the geometry.
///     This is required rather than optional: bodies overlap only partially, so a point in a
///     region another body never saw has no meaningful match at any distance, and the opening
///     unweighted solve has no defense against one. Choose it from the geometry rather than from
///     the expected residual, and err generous.
/// :param initial: an optional starting pose per body, one `Iso2` each. `None` starts every body
///     at the identity.
/// :param max_normal_angle: an optional cutoff, in radians, on the angle between a test point's
///     normal and its match's. Suppresses matches onto the far side of a thin wall.
/// :param refinement_steps: rounds of robust reweighting after the initial solve.
/// :param sigma_max: the MAGSAC++ upper noise bound, in the units of the geometry. Estimated from
///     the data when `None`.
/// :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
///     count.
/// :param dof: an optional `Dof3` constraint applied to every non-static body.
/// :param sample_spacing: the arc length spacing between correspondence samples on each member.
/// :param max_corner_angle: the largest turn, in radians, tolerated within one sample spacing of a
///     sample. This is what keeps samples off corners, where the closest point on another body
///     jumps rather than slides. Defaults to 60 degrees.
/// :param end_margin: the distance from either end of an open member within which samples are
///     discarded. Defaults to twice `sample_spacing`.
/// :param filter_distances: discard candidate correspondences further than this many standard
///     deviations above the mean candidate distance.
#[pyfunction]
#[pyo3(signature = (groups, max_distance, initial=None, max_normal_angle=None,
                    refinement_steps=4, sigma_max=None, patience=100, dof=None,
                    sample_spacing=1.0, max_corner_angle=None,
                    end_margin=None, filter_distances=Some(3.0)))]
#[allow(clippy::too_many_arguments)]
pub fn multi_curve_adjustment(
    groups: Vec<CurveGroup2>,
    max_distance: f64,
    initial: Option<Vec<Iso2>>,
    max_normal_angle: Option<f64>,
    refinement_steps: usize,
    sigma_max: Option<f64>,
    patience: usize,
    dof: Option<Dof3>,
    sample_spacing: f64,
    max_corner_angle: Option<f64>,
    end_margin: Option<f64>,
    filter_distances: Option<f64>,
) -> PyResult<MultiAlignOutcome2> {
    let bodies: Vec<engeom::CurveGroup2> = groups.iter().map(|g| g.get_inner().clone()).collect();
    let poses: Option<Vec<engeom::Iso2>> =
        initial.map(|v| v.iter().map(|t| *t.get_inner()).collect());

    let opts = engeom::geom2::align2::MultiAlignOptions2 {
        max_distance,
        max_normal_angle,
        refinement_steps,
        sigma_max,
        patience,
        dof: dof.map(Into::into),
    };

    // The defaults come from the core so the binding cannot drift from it. `filter_distances`
    // is the exception: `None` means "disabled" here, so its default stays in the signature.
    let defaults = engeom::geom2::align2::CAPParams::defaults(sample_spacing);
    let sample_opts = engeom::geom2::align2::CAPParams {
        sample_spacing,
        max_corner_angle: max_corner_angle.unwrap_or(defaults.max_corner_angle),
        end_margin: end_margin.unwrap_or(defaults.end_margin),
        filter_distances,
    };

    engeom::geom2::align2::multi_curve_adjustment(&bodies, poses.as_deref(), &opts, &sample_opts)
        .map(MultiAlignOutcome2::from_inner)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
