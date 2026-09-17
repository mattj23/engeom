use crate::bounding::Aabb3;
use crate::common::IndexMask;
use crate::conversions::{
    array_to_colors, array_to_points2, array_to_points3, array_to_unit_vectors3, array_to_vec,
    array_to_vectors3, colors_to_array, labels_to_array, points_to_array, scalars_to_array,
    unit_vectors_to_array,
};
use crate::geom3::Iso3;
use crate::half_edge3::{RepairOpts, RepairReport};
use crate::mesh::{Mesh3, MeshData3, PatchFilter};
use engeom::PointCloudOverlap;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::OnceLock;

/// Build an `engeom::io::Lptf3Load` from the flattened loader keyword arguments.
///
/// The smoothing parameters (`look_scale`, `weight_scale`, `max_move`) form an all-or-nothing
/// group that selects the gaussian smoothing filter. If none of them are given, the load is a
/// plain decimation controlled by `take_every` (`take_every <= 1` loads every point). Supplying
/// some but not all of the smoothing parameters is an error.
pub fn lptf3_load_from_args(
    take_every: u32,
    look_scale: Option<f64>,
    weight_scale: Option<f64>,
    max_move: Option<f64>,
) -> PyResult<engeom::io::Lptf3Load> {
    match (look_scale, weight_scale, max_move) {
        (None, None, None) => {
            if take_every <= 1 {
                Ok(engeom::io::Lptf3Load::All)
            } else {
                Ok(engeom::io::Lptf3Load::TakeEveryN(take_every))
            }
        }
        (Some(look_scale), Some(weight_scale), Some(max_move)) => {
            let p = engeom::io::Lptf3DsParams::new(take_every, look_scale, weight_scale, max_move);
            Ok(engeom::io::Lptf3Load::SmoothSample(p))
        }
        _ => Err(PyValueError::new_err(
            "The smoothing parameters `look_scale`, `weight_scale`, and `max_move` must all be \
             given together or all left unset",
        )),
    }
}

/// Parse the `invalid_normals` value passed to `PointCloud3.load_pcd`.
fn pcd_invalid_normals_from_str(s: &str) -> PyResult<engeom::io::PcdInvalidNormals> {
    match s {
        "error" => Ok(engeom::io::PcdInvalidNormals::Error),
        "drop_points" => Ok(engeom::io::PcdInvalidNormals::DropPoints),
        "drop_normals" => Ok(engeom::io::PcdInvalidNormals::DropNormals),
        _ => Err(PyValueError::new_err(format!(
            "Invalid invalid_normals '{s}', expected 'error', 'drop_points', or 'drop_normals'"
        ))),
    }
}

/// Build an `engeom::NormalOrientation3` from the flattened orientation keyword arguments.
///
/// The three arguments are mutually exclusive, and exactly one is required. Each argument carries
/// a different payload, so the supplied argument also identifies the orientation method without a
/// tagged value.
fn normal_orientation_from_args(
    viewpoint: Option<[f64; 3]>,
    viewpoints: Option<PyReadonlyArray2<f64>>,
    propagate_k: Option<usize>,
) -> PyResult<engeom::NormalOrientation3> {
    let given =
        viewpoint.is_some() as u8 + viewpoints.is_some() as u8 + propagate_k.is_some() as u8;

    if given != 1 {
        return Err(PyValueError::new_err(
            "Supply exactly one of `viewpoint`, `viewpoints`, or `propagate_k` to say where the \
             direction of the normals should come from",
        ));
    }

    if let Some(v) = viewpoint {
        return Ok(engeom::NormalOrientation3::Viewpoint(engeom::Point3::new(
            v[0], v[1], v[2],
        )));
    }

    if let Some(v) = viewpoints {
        let points = array_to_points3(&v.as_array())?;
        return Ok(engeom::NormalOrientation3::Viewpoints(points));
    }

    let k = propagate_k.expect("one of the three was given");
    if k == 0 {
        return Err(PyValueError::new_err(
            "`propagate_k` must be at least 1, since a point linked to no neighbors has nothing to \
             agree with",
        ));
    }
    Ok(engeom::NormalOrientation3::Propagate { k })
}

/// What a reconstruction did at each stage.
///
/// The report helps diagnose reconstruction problems. A `components` value greater than one
/// indicates that normal orientation required a separate seed guess in multiple regions. A large
/// `cells_skipped_unknown` value next to a small `cells_visited` value indicates that the band was
/// too thin to contain the surface.
#[pyclass(skip_from_py_object, module = "engeom.geom3")]
#[derive(Clone)]
pub struct ReconstructReport3 {
    inner: engeom::ReconstructReport3,
}

#[pymethods]
impl ReconstructReport3 {
    #[getter]
    fn points_used(&self) -> usize {
        self.inner.points_used
    }

    #[getter]
    fn points_dropped(&self) -> usize {
        self.inner.points_dropped
    }

    #[getter]
    fn active_blocks(&self) -> usize {
        self.inner.active_blocks
    }

    #[getter]
    fn known_voxels(&self) -> usize {
        self.inner.known_voxels
    }

    #[getter]
    fn cells_visited(&self) -> usize {
        self.inner.cells_visited
    }

    #[getter]
    fn cells_skipped_unknown(&self) -> usize {
        self.inner.cells_skipped_unknown
    }

    #[getter]
    fn normals_flipped(&self) -> Option<usize> {
        self.inner.orientation.as_ref().map(|o| o.flipped)
    }

    #[getter]
    fn components(&self) -> Option<usize> {
        self.inner.orientation.as_ref().and_then(|o| o.components)
    }

    #[getter]
    fn raw_vertices(&self) -> usize {
        self.inner.raw_vertices
    }

    #[getter]
    fn raw_faces(&self) -> usize {
        self.inner.raw_faces
    }

    #[getter]
    fn repair(&self) -> Option<RepairReport> {
        self.inner.repair.map(RepairReport::from_inner)
    }

    #[getter]
    fn faces(&self) -> usize {
        self.inner.faces
    }

    fn __repr__(&self) -> String {
        format!(
            "<ReconstructReport3 points_used={} faces={} known_voxels={} components={:?}>",
            self.inner.points_used,
            self.inner.faces,
            self.inner.known_voxels,
            self.components()
        )
    }
}

/// The pair of arrays `estimate_normals` hands back: an `(n, 3)` of unit normals and an `(n,)` of
/// per-point confidences.
type NormalEstimateArrays<'py> = (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray1<f64>>);

// ================================================================================================
// PointCloud3
// ================================================================================================

/// The unaccelerated point cloud container, holding the point buffer and its per-point attributes.
///
/// Numpy arrays handed out through the getters are cached and invalidated by any method which can
/// change what they hold, so repeated access does not repeatedly copy the buffer.
#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct PointCloud3 {
    inner: engeom::PointCloud3,

    /// The k-d tree over `inner`, built the first time a spatial query needs one and dropped by
    /// `clear_cached` whenever `inner` changes. See `index` for why this is a bare tree rather than
    /// a `CloudIndex3`, and why holding it is sound.
    ///
    /// `OnceLock` rather than `Option` so that a query can build it through `&self`, and rather than
    /// `OnceCell` because pyo3 requires a pyclass to be `Sync`. That keeps the
    /// query methods on shared borrows, which matters because `overlap_points_by_reciprocity` takes
    /// a second cloud: with `&mut self` a caller passing the same object twice would get an
    /// "already borrowed" error out of pyo3 rather than an answer.
    tree: OnceLock<engeom::KdTree3>,

    points: Option<Py<PyArray2<f64>>>,
    point_normals: Option<Py<PyArray2<f64>>>,
    point_colors: Option<Py<PyArray2<u8>>>,
    point_stdev: Option<Py<PyArray1<f64>>>,
    point_flat: Option<Py<PyArray2<f64>>>,
    voxel_coherence: Option<Py<PyArray1<f64>>>,
    voxel_count: Option<Py<PyArray1<u32>>>,
}

impl PointCloud3 {
    /// Drop everything derived from `inner`.
    ///
    /// **Every `&mut self` method must call this**, which is the discipline that keeps the cached
    /// tree honest. It is deliberately coarse: an attribute setter does not move any points and so
    /// does not really invalidate the tree, but distinguishing the cases would buy nothing at this
    /// layer and would give the rule an exception to forget.
    fn clear_cached(&mut self) {
        self.tree.take();
        self.points = None;
        self.point_normals = None;
        self.point_colors = None;
        self.point_stdev = None;
        self.point_flat = None;
        self.voxel_coherence = None;
        self.voxel_count = None;
    }

    /// The cloud paired with its cached tree, building the tree on first use.
    ///
    /// This is the one place in the codebase that calls `index_with_tree_unchecked`, and it upholds
    /// that method's obligation the only way it can be upheld: the tree is built from `self.inner`
    /// here and nowhere else, and `clear_cached` drops it on every mutation. Those two facts sit a
    /// few lines apart on purpose. If you add a `&mut self` method to this class and do not call
    /// `clear_cached`, queries afterwards will return confidently wrong indices.
    fn index(&self) -> PyResult<engeom::CloudIndex3<'_>> {
        if self.tree.get().is_none() {
            let tree = engeom::KdTree3::try_new(self.inner.points())
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            let _ = self.tree.set(tree);
        }

        let tree = self
            .tree
            .get()
            .expect("the tree was just built or already present");

        Ok(self.inner.index_with_tree_unchecked(tree))
    }

    pub fn get_inner(&self) -> &engeom::PointCloud3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::PointCloud3) -> Self {
        Self {
            inner,
            tree: OnceLock::new(),
            points: None,
            point_normals: None,
            point_colors: None,
            point_stdev: None,
            point_flat: None,
            voxel_coherence: None,
            voxel_count: None,
        }
    }
}

impl Clone for PointCloud3 {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl PointCloud3 {
    #[new]
    fn new<'py>(points: PyReadonlyArray2<'py, f64>) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        Ok(Self::from_inner(engeom::PointCloud3::new(points)))
    }

    /// Create an empty point cloud, with no points and no attributes.
    #[staticmethod]
    fn empty() -> Self {
        Self::from_inner(engeom::PointCloud3::empty())
    }

    // --- Serialization ---------------------------------------------------------------------

    #[staticmethod]
    fn load_ply(path: PathBuf) -> PyResult<Self> {
        let inner =
            engeom::PointCloud3::load_ply(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[staticmethod]
    #[pyo3(signature = (path, *, invalid_normals = "error"))]
    fn load_pcd(path: PathBuf, invalid_normals: &str) -> PyResult<Self> {
        let mut opts = engeom::io::PcdReadOpts::default();
        opts.invalid_normals = pcd_invalid_normals_from_str(invalid_normals)?;
        let inner = engeom::PointCloud3::load_pcd(&path, &opts)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[pyo3(signature = (path, binary = true))]
    fn save_ply(&self, path: PathBuf, binary: bool) -> PyResult<()> {
        let mut opts = engeom::io::PlyWriteOpts::default();
        opts.binary = binary;
        self.inner
            .save_ply(&path, &opts)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    #[staticmethod]
    #[pyo3(signature = (path, take_every=1, look_scale=None, weight_scale=None, max_move=None))]
    fn load_lptf3(
        path: PathBuf,
        take_every: u32,
        look_scale: Option<f64>,
        weight_scale: Option<f64>,
        max_move: Option<f64>,
    ) -> PyResult<Self> {
        let load = lptf3_load_from_args(take_every, look_scale, weight_scale, max_move)?;
        let inner =
            engeom::io::load_lptf3(&path, load).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Read a `.tcrpf3` row-organized point file.
    ///
    /// The load keywords have the same meaning as for `load_lptf3`. By default, every point is
    /// loaded. `take_every=n` thins to approximately one point per `n` row pitches in both
    /// directions. Supplying `look_scale`, `weight_scale`, and `max_move` adds Gaussian smoothing
    /// after thinning, using the discarded full-resolution points.
    #[staticmethod]
    #[pyo3(signature = (path, take_every=1, look_scale=None, weight_scale=None, max_move=None))]
    fn load_tc_row_points(
        path: PathBuf,
        take_every: u32,
        look_scale: Option<f64>,
        weight_scale: Option<f64>,
        max_move: Option<f64>,
    ) -> PyResult<Self> {
        let load = lptf3_load_from_args(take_every, look_scale, weight_scale, max_move)?;
        let inner = engeom::io::load_tc_row_points(&path, load)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    // --- Buffers and attributes --------------------------------------------------------------

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let array = points_to_array(self.inner.points());
            self.points = Some(array.into_pyarray(py).unbind());
        }
        self.points.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn point_normals<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<f64>>> {
        let values = self.inner.point_normals()?;
        if self.point_normals.is_none() {
            let array = unit_vectors_to_array(values);
            self.point_normals = Some(array.into_pyarray(py).unbind());
        }
        Some(self.point_normals.as_ref().unwrap().bind(py))
    }

    #[getter]
    fn point_colors<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<u8>>> {
        let values = self.inner.point_colors()?;
        if self.point_colors.is_none() {
            let array = colors_to_array(values);
            self.point_colors = Some(array.into_pyarray(py).unbind());
        }
        Some(self.point_colors.as_ref().unwrap().bind(py))
    }

    #[getter]
    fn point_stdev<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray1<f64>>> {
        let values = self.inner.point_stdev()?;
        if self.point_stdev.is_none() {
            let array = scalars_to_array(values);
            self.point_stdev = Some(array.into_pyarray(py).unbind());
        }
        Some(self.point_stdev.as_ref().unwrap().bind(py))
    }

    #[getter]
    fn point_flat<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<f64>>> {
        let values = self.inner.point_flat()?;
        if self.point_flat.is_none() {
            let array = points_to_array(values);
            self.point_flat = Some(array.into_pyarray(py).unbind());
        }
        Some(self.point_flat.as_ref().unwrap().bind(py))
    }

    #[pyo3(signature = (values=None))]
    fn set_point_normals<'py>(
        &mut self,
        values: Option<PyReadonlyArray2<'py, f64>>,
    ) -> PyResult<()> {
        let values = values
            .map(|v| array_to_unit_vectors3(&v.as_array()))
            .transpose()?;
        self.clear_cached();
        self.inner
            .set_point_normals(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (values=None))]
    fn set_point_colors<'py>(&mut self, values: Option<PyReadonlyArray2<'py, u8>>) -> PyResult<()> {
        let values = values.map(|v| array_to_colors(&v.as_array())).transpose()?;
        self.clear_cached();
        self.inner
            .set_point_colors(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (values=None))]
    fn set_point_stdev<'py>(&mut self, values: Option<PyReadonlyArray1<'py, f64>>) -> PyResult<()> {
        let values = values.map(|v| array_to_vec(&v)).transpose()?;
        self.clear_cached();
        self.inner
            .set_point_stdev(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (values=None))]
    fn set_point_flat<'py>(&mut self, values: Option<PyReadonlyArray2<'py, f64>>) -> PyResult<()> {
        let values = values
            .map(|v| array_to_points2(&v.as_array()))
            .transpose()?;
        self.clear_cached();
        self.inner
            .set_point_flat(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    // --- Operations ------------------------------------------------------------------------

    fn transform_in_place(&mut self, iso: &Iso3) {
        self.clear_cached();
        self.inner.transform_in_place(iso.get_inner());
    }

    fn transform_copy(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transform_copy(iso.get_inner()))
    }

    fn scale_in_place(&mut self, scale: f64) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .scale_in_place(scale)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn scale_copy(&self, scale: f64) -> PyResult<Self> {
        let inner = self
            .inner
            .scale_copy(scale)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn append_in_place(&mut self, other: &PointCloud3) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .append_in_place(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn extract_subset_indices(&self, indices: Vec<usize>) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_indices(&indices)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Extract the points selected by a mask as a new cloud, carrying attributes across.
    fn extract_subset_points(&self, point_mask: &IndexMask) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_points(point_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Reduce the cloud onto a coarser grid, replacing the points in each voxel with a single
    /// averaged point.
    ///
    /// This creates new points rather than selecting existing ones, so the result is less noisy
    /// than the input but is no longer a set of measurements that were actually taken. Output
    /// points are voxel centroids, not voxel centers, so there is no minimum spacing guarantee;
    /// use `sample_poisson_disk` when you need one.
    ///
    /// Every attribute is combined by whatever rule suits it, and two are added: see
    /// `voxel_coherence` and `voxel_count`.
    fn reduce_by_voxel(&self, voxel_size: f64) -> PyResult<Self> {
        self.inner
            .reduce_by_voxel(voxel_size)
            .map_err(|e| PyValueError::new_err(e.to_string()))
            .map(Self::from_inner)
    }

    /// How well the normals within each voxel agreed, in `[0, 1]`, or `None` on a cloud which is
    /// not the output of `reduce_by_voxel` or whose input had no normals.
    ///
    /// Near 1 where a voxel's normals all pointed the same way, falling toward 0 where the voxel
    /// straddled an edge or a thin wall. A low value means the averaged point in that voxel is a
    /// blend of surfaces that face different directions, so it is a natural weight to apply when
    /// using a reduced cloud for fitting or alignment.
    #[getter]
    fn voxel_coherence<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray1<f64>>> {
        let values = self
            .inner
            .point_attr(engeom::VOXEL_COHERENCE_ATTR)
            .and_then(|a| a.as_scalar())?;

        if self.voxel_coherence.is_none() {
            let array = scalars_to_array(values);
            self.voxel_coherence = Some(array.into_pyarray(py).unbind());
        }
        Some(self.voxel_coherence.as_ref().unwrap().bind(py))
    }

    /// How many input points went into each output point, or `None` on a cloud which is not the
    /// output of `reduce_by_voxel`.
    ///
    /// Averaging `n` independent measurements lowers their noise by `sqrt(n)`, so this is what says
    /// how much a given reduced point gained, and is useful as a weight in its own right.
    #[getter]
    fn voxel_count<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray1<u32>>> {
        let values = self
            .inner
            .point_attr(engeom::VOXEL_COUNT_ATTR)
            .and_then(|a| a.as_label())?;

        if self.voxel_count.is_none() {
            let array = labels_to_array(values);
            self.voxel_count = Some(array.into_pyarray(py).unbind());
        }
        Some(self.voxel_count.as_ref().unwrap().bind(py))
    }

    /// The axis-aligned bounding box of the points.
    fn compute_aabb(&self) -> Aabb3 {
        Aabb3::from_inner(self.inner.compute_aabb())
    }

    #[getter]
    fn point_count(&self) -> usize {
        self.inner.point_count()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    // --- Spatial queries -------------------------------------------------------------------
    //
    // All of these go through `index`, which builds the tree once and caches it until the cloud is
    // mutated. A `CloudIndex3` cannot be held in a field here because it borrows the cloud it
    // indexes, so what is cached is the bare tree and the index is re-paired per call.

    /// Poisson disk sample the cloud, returning the indices of a subset of points no two of which
    /// are closer together than `radius`.
    fn sample_poisson_disk(&self, radius: f64) -> PyResult<Vec<usize>> {
        // No index here: the sampler builds its own small tree over a voxel-downsampled subset.
        Ok(self.inner.sample_poisson_disk(radius).to_indices())
    }

    /// Poisson disk sample the cloud and return the result as a new cloud, carrying attributes.
    fn extract_poisson_sample(&self, radius: f64) -> PyResult<Self> {
        self.inner
            .extract_poisson_sample(radius)
            .map_err(|e| PyValueError::new_err(e.to_string()))
            .map(Self::from_inner)
    }

    /// Find the indices of points in this cloud which overlap another cloud, by checking that the
    /// closest point in each direction agrees to within `max_distance`.
    fn overlap_points_by_reciprocity(
        &self,
        other: &PointCloud3,
        max_distance: f64,
    ) -> PyResult<Vec<usize>> {
        Ok(self
            .index()?
            .overlap_by_reciprocity(&other.index()?, max_distance))
    }

    /// Estimate a normal at every point by fitting a plane to the neighbors within `radius`.
    ///
    /// Returns the normals and a per-point confidence in `[0, 1]` saying how plane-like the
    /// neighborhood was, which is low on edges and corners and in sparse regions. Points with
    /// fewer than three neighbors get `+Z` at zero confidence.
    ///
    /// `must_match` supplies one direction per point which the estimate is flipped to agree with.
    /// A plane fit recovers an axis without a direction and cannot resolve the sign independently.
    /// This argument was originally required. For scan data, the usual value is the vector from
    /// each point toward the sensor.
    /// # Why `radius` carries a default it must not take
    ///
    /// `must_match` was originally the only way to resolve a normal's sign. It preceded `radius`,
    /// and both arguments were positional. It is now one of four alternatives and must be optional,
    /// while Python requires every argument after an argument with a default to also have a default.
    /// Reordering the arguments would break existing positional calls, so `radius` has a default
    /// and the method rejects the call if it is omitted.
    #[pyo3(signature = (must_match=None, radius=None, *, viewpoint=None, viewpoints=None, propagate_k=None))]
    fn estimate_normals<'py>(
        &self,
        py: Python<'py>,
        must_match: Option<PyReadonlyArray2<'py, f64>>,
        radius: Option<f64>,
        viewpoint: Option<[f64; 3]>,
        viewpoints: Option<PyReadonlyArray2<'py, f64>>,
        propagate_k: Option<usize>,
    ) -> PyResult<NormalEstimateArrays<'py>> {
        let radius = radius.ok_or_else(|| {
            PyValueError::new_err(
                "`radius` is required. It sits after `must_match` and carries a default only so \
                 that `must_match` can be left out, so pass it by keyword when you are not \
                 supplying `must_match`.",
            )
        })?;

        let estimates = if let Some(must_match) = must_match {
            if viewpoint.is_some() || viewpoints.is_some() || propagate_k.is_some() {
                return Err(PyValueError::new_err(
                    "`must_match` already says which way the normals face, so it cannot be \
                     combined with `viewpoint`, `viewpoints`, or `propagate_k`",
                ));
            }

            let must_match = array_to_vectors3(&must_match.as_array())?;
            self.index()?
                .estimate_normals(&must_match, radius)
                .map_err(|e| PyValueError::new_err(e.to_string()))?
        } else {
            let orientation = normal_orientation_from_args(viewpoint, viewpoints, propagate_k)?;
            self.index()?
                .estimate_normals_oriented(radius, &orientation)
                .map_err(|e| PyValueError::new_err(e.to_string()))?
                .0
        };

        let normals = unit_vectors_to_array(&estimates.normals).into_pyarray(py);
        let confidence = scalars_to_array(&estimates.confidence).into_pyarray(py);

        Ok((normals, confidence))
    }

    /// Build a triangle mesh from this cloud.
    ///
    /// See the Python stub for the stages, what the result promises, and what it does not.
    #[pyo3(signature = (
        voxel_size,
        *,
        band_width=2.0,
        normal_radius=None,
        viewpoint=None,
        viewpoints=None,
        propagate_k=None,
        min_normal_confidence=0.0,
        repair=None,
        patch_filter=None,
        smooth_iterations=0,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn reconstruct_surface(
        &self,
        voxel_size: f64,
        band_width: f64,
        normal_radius: Option<f64>,
        viewpoint: Option<[f64; 3]>,
        viewpoints: Option<PyReadonlyArray2<f64>>,
        propagate_k: Option<usize>,
        min_normal_confidence: f64,
        repair: Option<RepairOpts>,
        patch_filter: Option<PatchFilter>,
        smooth_iterations: usize,
    ) -> PyResult<(MeshData3, ReconstructReport3)> {
        let normals = match normal_radius {
            None => {
                if viewpoint.is_some() || viewpoints.is_some() || propagate_k.is_some() {
                    return Err(PyValueError::new_err(
                        "`viewpoint`, `viewpoints`, and `propagate_k` only apply when \
                         `normal_radius` is given, since without it the cloud's own normals are \
                         used as they are",
                    ));
                }
                engeom::NormalSource3::Existing
            }
            Some(radius) => engeom::NormalSource3::Estimate {
                radius,
                orientation: normal_orientation_from_args(viewpoint, viewpoints, propagate_k)?,
            },
        };

        // An unset `RepairOpts` runs the default passes. Pass `RepairOpts.none()` to disable every
        // pass, which has the same effect as skipping repair.
        //
        // Match the `ReconstructOpts3` default. Reconstruction omits two passes made unnecessary by
        // its extraction, and callers that supply no options must receive the same preset.
        let repair_opts = repair
            .map(|r| *r.get_inner())
            .unwrap_or_else(engeom::geom3::half_edge3::RepairOpts::assuming_oriented_edges);

        let opts = engeom::ReconstructOpts3::new(voxel_size)
            .with_band_width(band_width)
            .with_normals(normals)
            .with_min_normal_confidence(min_normal_confidence)
            .with_repair(Some(repair_opts))
            .with_patch_filter(patch_filter.map(|f| *f.get_inner()))
            .with_smooth_iterations(smooth_iterations);

        let (mesh, report) = self
            .index()?
            .reconstruct_surface(&opts)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok((
            MeshData3::from_inner(mesh),
            ReconstructReport3 { inner: report },
        ))
    }

    /// The median distance from a point to its nearest neighbor, used to select the voxel size.
    fn estimate_point_spacing(&self) -> PyResult<f64> {
        Ok(self.index()?.estimate_point_spacing())
    }

    /// Find the indices of points in this cloud which overlap a mesh, by checking that the closest
    /// point in each direction agrees to within `max_distance`.
    fn overlap_mesh_by_reciprocity(&self, mesh: &Mesh3, max_distance: f64) -> PyResult<Vec<usize>> {
        Ok(self
            .index()?
            .overlap_by_reciprocity(mesh.get_inner(), max_distance))
    }

    fn __len__(&self) -> usize {
        self.inner.point_count()
    }

    fn __repr__(&self) -> String {
        format!("<PointCloud3 {} points>", self.inner.point_count())
    }
}
