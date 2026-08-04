//! Bindings for the half-edge mesh and the topology edits driven through it.
//!
//! `MeshEditor` is deliberately not bound. It borrows the mesh it started from, which a `#[pyclass]`
//! cannot express, and everything it wraps other than the attribute-loss check is reachable here.
//! The staleness guard which used to live on it now lives on `HalfEdgeMesh3` itself, so the Python
//! path is protected the same way the Rust one is.

use crate::mesh::Mesh3;
use engeom::geom3::half_edge3::{
    BestEffortOpts as InnerBestEffortOpts, DecimateOpts as InnerDecimateOpts, DecimateTarget,
    ErrorMethod, HalfEdgeMesh3 as InnerHalfEdgeMesh3, Placement, QuadricKind,
    RepairOpts as InnerRepairOpts,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// Token parsing ─────────────────────────────────────────────────────────────

fn placement_from_str(s: &str) -> PyResult<Placement> {
    match s {
        "optimal" => Ok(Placement::Optimal),
        "midpoint" => Ok(Placement::Midpoint),
        "endpoint" => Ok(Placement::Endpoint),
        _ => Err(PyValueError::new_err(format!(
            "Invalid placement '{s}', expected 'optimal', 'midpoint', or 'endpoint'"
        ))),
    }
}

fn placement_to_str(v: Placement) -> &'static str {
    match v {
        Placement::Optimal => "optimal",
        Placement::Midpoint => "midpoint",
        Placement::Endpoint => "endpoint",
    }
}

fn quadric_from_str(s: &str) -> PyResult<QuadricKind> {
    match s {
        "triangle" => Ok(QuadricKind::Triangle),
        "probabilistic" => Ok(QuadricKind::Probabilistic),
        _ => Err(PyValueError::new_err(format!(
            "Invalid quadric kind '{s}', expected 'triangle' or 'probabilistic'"
        ))),
    }
}

fn quadric_to_str(v: QuadricKind) -> &'static str {
    match v {
        QuadricKind::Triangle => "triangle",
        QuadricKind::Probabilistic => "probabilistic",
    }
}

fn error_method_from_str(s: &str) -> PyResult<ErrorMethod> {
    match s {
        "contraction" => Ok(ErrorMethod::Contraction),
        "affine_face_map" => Ok(ErrorMethod::AffineFaceMap),
        "projected_overlay" => Ok(ErrorMethod::ProjectedOverlay),
        _ => Err(PyValueError::new_err(format!(
            "Invalid error method '{s}', expected 'contraction', 'affine_face_map', or \
             'projected_overlay'"
        ))),
    }
}

fn error_method_to_str(v: ErrorMethod) -> &'static str {
    match v {
        ErrorMethod::Contraction => "contraction",
        ErrorMethod::AffineFaceMap => "affine_face_map",
        ErrorMethod::ProjectedOverlay => "projected_overlay",
    }
}

/// Resolve the stopping target from the two optional keyword arguments which select it.
///
/// `DecimateTarget` is presence-distinguishable, so it crosses as `face_count=` / `ratio=` with at
/// most one supplied rather than as a wrapper class.
fn target_from_args(face_count: Option<usize>, ratio: Option<f64>) -> PyResult<DecimateTarget> {
    match (face_count, ratio) {
        (None, None) => Ok(DecimateTarget::ToTolerance),
        (Some(n), None) => Ok(DecimateTarget::FaceCount(n)),
        (None, Some(r)) => Ok(DecimateTarget::Ratio(r)),
        (Some(_), Some(_)) => Err(PyValueError::new_err(
            "Supply at most one of 'face_count' or 'ratio'; decimation stops at whichever target \
             is given, or runs to the tolerance if neither is",
        )),
    }
}

fn target_parts(t: DecimateTarget) -> (Option<usize>, Option<f64>) {
    match t {
        DecimateTarget::ToTolerance => (None, None),
        DecimateTarget::FaceCount(n) => (Some(n), None),
        DecimateTarget::Ratio(r) => (None, Some(r)),
    }
}

// Options ───────────────────────────────────────────────────────────────────

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct RepairOpts {
    inner: InnerRepairOpts,
}

impl RepairOpts {
    pub fn get_inner(&self) -> &InnerRepairOpts {
        &self.inner
    }
}

#[pymethods]
impl RepairOpts {
    #[new]
    #[pyo3(signature = (
        *,
        drop_degenerate = true,
        drop_duplicate_faces = true,
        resolve_nonmanifold_edges = true,
        split_bowtie_vertices = true,
        orient_consistently = true,
        drop_isolated_vertices = true,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        drop_degenerate: bool,
        drop_duplicate_faces: bool,
        resolve_nonmanifold_edges: bool,
        split_bowtie_vertices: bool,
        orient_consistently: bool,
        drop_isolated_vertices: bool,
    ) -> Self {
        Self {
            inner: InnerRepairOpts::default()
                .with_drop_degenerate(drop_degenerate)
                .with_drop_duplicate_faces(drop_duplicate_faces)
                .with_resolve_nonmanifold_edges(resolve_nonmanifold_edges)
                .with_split_bowtie_vertices(split_bowtie_vertices)
                .with_orient_consistently(orient_consistently)
                .with_drop_isolated_vertices(drop_isolated_vertices),
        }
    }

    /// Every pass disabled, as a base for turning on only what you want.
    #[staticmethod]
    fn none() -> Self {
        Self {
            inner: InnerRepairOpts::none(),
        }
    }

    #[getter]
    fn drop_degenerate(&self) -> bool {
        self.inner.drop_degenerate
    }

    #[getter]
    fn drop_duplicate_faces(&self) -> bool {
        self.inner.drop_duplicate_faces
    }

    #[getter]
    fn resolve_nonmanifold_edges(&self) -> bool {
        self.inner.resolve_nonmanifold_edges
    }

    #[getter]
    fn split_bowtie_vertices(&self) -> bool {
        self.inner.split_bowtie_vertices
    }

    #[getter]
    fn orient_consistently(&self) -> bool {
        self.inner.orient_consistently
    }

    #[getter]
    fn drop_isolated_vertices(&self) -> bool {
        self.inner.drop_isolated_vertices
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct DecimateOpts {
    inner: InnerDecimateOpts,
}

impl DecimateOpts {
    pub fn get_inner(&self) -> &InnerDecimateOpts {
        &self.inner
    }
}

#[pymethods]
impl DecimateOpts {
    #[new]
    #[pyo3(signature = (
        tol,
        *,
        face_count = None,
        ratio = None,
        max_normal_deviation = None,
        min_aspect = None,
        lock_boundary = None,
        boundary_tol = None,
        placement = None,
        quadric = None,
        error_method = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        tol: f64,
        face_count: Option<usize>,
        ratio: Option<f64>,
        max_normal_deviation: Option<f64>,
        min_aspect: Option<f64>,
        lock_boundary: Option<bool>,
        boundary_tol: Option<f64>,
        placement: Option<&str>,
        quadric: Option<&str>,
        error_method: Option<&str>,
    ) -> PyResult<Self> {
        let mut inner = InnerDecimateOpts::to_tolerance(tol);
        inner.target = target_from_args(face_count, ratio)?;

        if let Some(v) = max_normal_deviation {
            inner = inner.with_max_normal_deviation(v);
        }
        if let Some(v) = min_aspect {
            inner = inner.with_min_aspect(v);
        }
        if let Some(v) = lock_boundary {
            inner = inner.with_lock_boundary(v);
        }
        if let Some(v) = boundary_tol {
            inner = inner.with_boundary_tol(v);
        }
        if let Some(v) = placement {
            inner = inner.with_placement(placement_from_str(v)?);
        }
        if let Some(v) = quadric {
            inner = inner.with_quadric(quadric_from_str(v)?);
        }
        if let Some(v) = error_method {
            inner = inner.with_error_method(error_method_from_str(v)?);
        }

        Ok(Self { inner })
    }

    #[getter]
    fn tol(&self) -> f64 {
        self.inner.tol
    }

    #[getter]
    fn face_count(&self) -> Option<usize> {
        target_parts(self.inner.target).0
    }

    #[getter]
    fn ratio(&self) -> Option<f64> {
        target_parts(self.inner.target).1
    }

    #[getter]
    fn max_normal_deviation(&self) -> f64 {
        self.inner.max_normal_deviation
    }

    #[getter]
    fn min_aspect(&self) -> f64 {
        self.inner.min_aspect
    }

    #[getter]
    fn lock_boundary(&self) -> bool {
        self.inner.lock_boundary
    }

    #[getter]
    fn boundary_tol(&self) -> Option<f64> {
        self.inner.boundary_tol
    }

    #[getter]
    fn placement(&self) -> &'static str {
        placement_to_str(self.inner.placement)
    }

    #[getter]
    fn quadric(&self) -> &'static str {
        quadric_to_str(self.inner.quadric)
    }

    #[getter]
    fn error_method(&self) -> &'static str {
        error_method_to_str(self.inner.error_method)
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct BestEffortOpts {
    inner: InnerBestEffortOpts,
}

impl BestEffortOpts {
    pub fn get_inner(&self) -> &InnerBestEffortOpts {
        &self.inner
    }
}

#[pymethods]
impl BestEffortOpts {
    #[new]
    #[pyo3(signature = (
        tol,
        *,
        face_count = None,
        ratio = None,
        max_normal_deviation = None,
        min_aspect = None,
        lock_boundary = None,
        placement = None,
        quadric = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        tol: f64,
        face_count: Option<usize>,
        ratio: Option<f64>,
        max_normal_deviation: Option<f64>,
        min_aspect: Option<f64>,
        lock_boundary: Option<bool>,
        placement: Option<&str>,
        quadric: Option<&str>,
    ) -> PyResult<Self> {
        let mut inner = InnerBestEffortOpts::to_tolerance(tol);
        inner.target = target_from_args(face_count, ratio)?;

        if let Some(v) = max_normal_deviation {
            inner = inner.with_max_normal_deviation(v);
        }
        if let Some(v) = min_aspect {
            inner = inner.with_min_aspect(v);
        }
        if let Some(v) = lock_boundary {
            inner = inner.with_lock_boundary(v);
        }
        if let Some(v) = placement {
            inner = inner.with_placement(placement_from_str(v)?);
        }
        if let Some(v) = quadric {
            inner = inner.with_quadric(quadric_from_str(v)?);
        }

        Ok(Self { inner })
    }

    #[getter]
    fn tol(&self) -> f64 {
        self.inner.tol
    }

    #[getter]
    fn face_count(&self) -> Option<usize> {
        target_parts(self.inner.target).0
    }

    #[getter]
    fn ratio(&self) -> Option<f64> {
        target_parts(self.inner.target).1
    }

    #[getter]
    fn max_normal_deviation(&self) -> f64 {
        self.inner.max_normal_deviation
    }

    #[getter]
    fn min_aspect(&self) -> f64 {
        self.inner.min_aspect
    }

    #[getter]
    fn lock_boundary(&self) -> bool {
        self.inner.lock_boundary
    }

    #[getter]
    fn placement(&self) -> &'static str {
        placement_to_str(self.inner.placement)
    }

    #[getter]
    fn quadric(&self) -> &'static str {
        quadric_to_str(self.inner.quadric)
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

// Reports ───────────────────────────────────────────────────────────────────

#[pyclass(skip_from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct RepairReport {
    inner: engeom::geom3::half_edge3::RepairReport,
}

#[pymethods]
impl RepairReport {
    #[getter]
    fn degenerate_removed(&self) -> usize {
        self.inner.degenerate_removed
    }

    #[getter]
    fn duplicate_faces_removed(&self) -> usize {
        self.inner.duplicate_faces_removed
    }

    #[getter]
    fn nonmanifold_edges(&self) -> usize {
        self.inner.nonmanifold_edges
    }

    #[getter]
    fn faces_dropped_at_nonmanifold(&self) -> usize {
        self.inner.faces_dropped_at_nonmanifold
    }

    #[getter]
    fn faces_dropped_for_orientation(&self) -> usize {
        self.inner.faces_dropped_for_orientation
    }

    #[getter]
    fn bowtie_vertices_split(&self) -> usize {
        self.inner.bowtie_vertices_split
    }

    #[getter]
    fn points_added_by_splitting(&self) -> usize {
        self.inner.points_added_by_splitting
    }

    #[getter]
    fn faces_reoriented(&self) -> usize {
        self.inner.faces_reoriented
    }

    #[getter]
    fn nonorientable_components(&self) -> usize {
        self.inner.nonorientable_components
    }

    #[getter]
    fn faces_rejected_by_ingest(&self) -> usize {
        self.inner.faces_rejected_by_ingest
    }

    /// Whether the repair had to change anything at all.
    fn is_clean(&self) -> bool {
        self.inner.is_clean()
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

#[pyclass(skip_from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct DecimateStats {
    inner: engeom::geom3::half_edge3::DecimateStats,
}

#[pymethods]
impl DecimateStats {
    #[getter]
    fn evaluations(&self) -> usize {
        self.inner.evaluations
    }

    #[getter]
    fn acceptance_tests(&self) -> usize {
        self.inner.acceptance_tests
    }

    #[getter]
    fn veto_normal(&self) -> usize {
        self.inner.veto_normal
    }

    #[getter]
    fn veto_aspect(&self) -> usize {
        self.inner.veto_aspect
    }

    #[getter]
    fn veto_degenerate(&self) -> usize {
        self.inner.veto_degenerate
    }

    #[getter]
    fn veto_error(&self) -> usize {
        self.inner.veto_error
    }

    #[getter]
    fn method_not_applicable(&self) -> usize {
        self.inner.method_not_applicable
    }

    /// Every veto, however it fired.
    fn vetoes(&self) -> usize {
        self.inner.vetoes()
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

#[pyclass(skip_from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct DecimateReport {
    inner: engeom::geom3::half_edge3::DecimateReport,
}

#[pymethods]
impl DecimateReport {
    #[getter]
    fn collapses(&self) -> usize {
        self.inner.collapses
    }

    #[getter]
    fn faces_before(&self) -> usize {
        self.inner.faces_before
    }

    #[getter]
    fn faces_after(&self) -> usize {
        self.inner.faces_after
    }

    #[getter]
    fn stats(&self) -> DecimateStats {
        DecimateStats {
            inner: self.inner.stats,
        }
    }

    /// The fraction of the starting face count which survived.
    fn ratio(&self) -> f64 {
        self.inner.ratio()
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

// The mesh ──────────────────────────────────────────────────────────────────

/// Pinned to the thread which created it.
///
/// `alum` holds its mesh properties in `Rc<RefCell<..>>`, so the underlying structure is neither
/// `Send` nor `Sync` and cannot be handed to another thread. `unsendable` is how PyO3 expresses
/// that: touching the object from a different thread raises rather than racing. Ordinary
/// single-threaded use is unaffected, and a caller who needs the result elsewhere should convert
/// with `to_mesh` first, since `Mesh3` has no such restriction.
#[pyclass(unsendable, module = "engeom.geom3")]
pub struct HalfEdgeMesh3 {
    inner: InnerHalfEdgeMesh3,
    repair: engeom::geom3::half_edge3::RepairReport,
}

#[pymethods]
impl HalfEdgeMesh3 {
    #[new]
    #[pyo3(signature = (mesh, *, repair = None))]
    fn new(mesh: &Mesh3, repair: Option<RepairOpts>) -> PyResult<Self> {
        let (inner, report) = match repair {
            Some(opts) => engeom::geom3::half_edge3::HalfEdgeMesh3::from_mesh_repaired(
                mesh.get_inner(),
                opts.get_inner(),
            )
            .map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => (
                InnerHalfEdgeMesh3::try_from(mesh.get_inner())
                    .map_err(|e| PyValueError::new_err(e.to_string()))?,
                engeom::geom3::half_edge3::RepairReport::default(),
            ),
        };

        Ok(Self {
            inner,
            repair: report,
        })
    }

    /// What the repairing constructor had to change, or an empty report if none was requested.
    #[getter]
    fn repair_report(&self) -> RepairReport {
        RepairReport { inner: self.repair }
    }

    #[getter]
    fn vertex_count(&self) -> usize {
        self.inner.vertex_count()
    }

    #[getter]
    fn face_count(&self) -> usize {
        self.inner.face_count()
    }

    #[getter]
    fn edge_count(&self) -> usize {
        self.inner.edge_count()
    }

    /// Whether the decimation error volume still accounts for everything done to this mesh.
    ///
    /// False once `decimate_best_effort` has run, after which `decimate_guaranteed` refuses.
    #[getter]
    fn error_volume_is_current(&self) -> bool {
        self.inner.error_volume_is_current()
    }

    /// Decimate the mesh, keeping both surfaces within the tolerance of each other.
    fn decimate_guaranteed(&mut self, opts: &DecimateOpts) -> PyResult<DecimateReport> {
        self.inner
            .decimate_guaranteed(opts.get_inner())
            .map(|inner| DecimateReport { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Decimate the mesh by estimated deviation rather than by a guaranteed bound.
    fn decimate_best_effort(&mut self, opts: &BestEffortOpts) -> PyResult<DecimateReport> {
        self.inner
            .decimate_best_effort(opts.get_inner())
            .map(|inner| DecimateReport { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Smooth the mesh, moving each vertex onto the average height of its one-ring.
    #[pyo3(signature = (iterations = 1))]
    fn smooth(&mut self, iterations: usize) -> PyResult<()> {
        for _ in 0..iterations {
            self.inner
                .neighborhood_smooth()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        }
        Ok(())
    }

    /// Discard elements marked deleted, compacting the structure.
    fn garbage_collect(&mut self) -> PyResult<()> {
        self.inner
            .garbage_collect()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Build a `Mesh3` from the current state.
    ///
    /// Takes `&mut self` because it compacts the structure first, which is a real mutation and the
    /// reason a caller's stale indices are invalid afterwards.
    #[allow(clippy::wrong_self_convention)]
    #[pyo3(signature = (is_solid = false))]
    fn to_mesh(&mut self, is_solid: bool) -> PyResult<Mesh3> {
        self.inner
            .to_mesh(is_solid)
            .map(Mesh3::from_inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "HalfEdgeMesh3({} vertices, {} faces, {} edges)",
            self.inner.vertex_count(),
            self.inner.face_count(),
            self.inner.edge_count()
        )
    }
}
