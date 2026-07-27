use crate::conversions::{
    array_to_colors, array_to_points3, array_to_unit_vectors3, array_to_vec, colors_to_array,
    points_to_array, scalars_to_array, unit_vectors_to_array, vectors_to_array,
};
use crate::geom3::Iso3;
use crate::mesh::Mesh;
use engeom::{PointCloudFeatures, PointCloudKdTree, PointCloudOverlap};
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

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

#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct PointCloud {
    inner: engeom::PointCloud,
    points: Option<Py<PyArray2<f64>>>,
    normals: Option<Py<PyArray2<f64>>>,
    colors: Option<Py<PyArray2<u8>>>,
    matched_tree: Option<engeom::common::kd_tree::MatchedTree<3>>,
}

impl PointCloud {
    fn clear_cached(&mut self) {
        self.points = None;
        self.normals = None;
        self.colors = None;
        self.matched_tree = None;
    }

    pub fn get_inner(&self) -> &engeom::PointCloud {
        &self.inner
    }

    pub fn from_inner(inner: engeom::PointCloud) -> Self {
        Self {
            inner,
            points: None,
            normals: None,
            colors: None,
            matched_tree: None,
        }
    }

    pub fn with_tree(&mut self) -> PyResult<PointCloudKdTree<'_>> {
        if self.matched_tree.is_none() {
            let tree = self
                .inner
                .create_matched_tree()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            self.matched_tree = Some(tree);
        }

        let tree = self.matched_tree.as_ref().unwrap();
        PointCloudKdTree::try_new(&self.inner, tree)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

impl Clone for PointCloud {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl PointCloud {
    #[new]
    fn new<'py>(points: PyReadonlyArray2<'py, f64>) -> PyResult<Self> {
        let cloud_points = array_to_points3(&points.as_array())?;
        let cloud = engeom::PointCloud::try_new(cloud_points, None, None, None)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(cloud))
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

    #[staticmethod]
    fn load_bxyz(path: PathBuf) -> PyResult<Self> {
        let inner = engeom::io::load_bxyz(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    fn append(&mut self, other: &PointCloud) -> PyResult<()> {
        self.clear_cached();
        let clone = other.inner.clone();
        self.inner
            .merge(clone)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let array = points_to_array(self.inner.points());
            self.points = Some(array.into_pyarray(py).unbind());
        }
        self.points.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn colors<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<u8>>> {
        if let Some(colors) = self.inner.colors() {
            if self.colors.is_none() {
                let flat_colors = colors.iter().flatten().copied().collect::<Vec<_>>();
                let array = Array2::from_shape_vec((self.inner.points().len(), 3), flat_colors)
                    .expect("Failed to create color array");
                self.colors = Some(array.into_pyarray(py).unbind());
            }

            Some(self.colors.as_ref().unwrap().bind(py))
        } else {
            None
        }
    }

    #[getter]
    fn normals<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<f64>>> {
        if let Some(normals) = self.inner.normals() {
            if self.normals.is_none() {
                let n = normals.iter().map(|v| v.into_inner()).collect::<Vec<_>>();
                let array = vectors_to_array(&n);
                self.normals = Some(array.into_pyarray(py).unbind());
            }

            Some(self.normals.as_ref().unwrap().bind(py))
        } else {
            None
        }
    }

    fn transform_by(&mut self, transform: &Iso3) {
        self.clear_cached();
        self.inner.transform_by(transform.get_inner());
    }

    fn create_from_indices(&self, indices: Vec<usize>) -> PyResult<Self> {
        let inner = self
            .inner
            .create_from_indices(&indices)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(inner))
    }

    fn sample_poisson_disk(&mut self, radius: f64) -> PyResult<Vec<usize>> {
        let with_tree = self
            .with_tree()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let mask = with_tree.sample_poisson_disk(radius);

        Ok(mask.to_indices())
    }

    fn create_from_poisson_sample(&mut self, radius: f64) -> PyResult<Self> {
        let with_tree = self
            .with_tree()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        with_tree
            .create_from_poisson_sample(radius)
            .map_err(|e| PyValueError::new_err(e.to_string()))
            .map(Self::from_inner)
    }

    fn overlap_points_by_reciprocity(
        &mut self,
        other: &mut PointCloud,
        max_distance: f64,
    ) -> PyResult<Vec<usize>> {
        let this_with_tree = self
            .with_tree()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let other_with_tree = other
            .with_tree()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(this_with_tree.overlap_by_reciprocity(&other_with_tree, max_distance))
    }

    fn overlap_mesh_by_reciprocity(
        &mut self,
        mesh: &Mesh,
        max_distance: f64,
    ) -> PyResult<Vec<usize>> {
        let this_with_tree = self
            .with_tree()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(this_with_tree.overlap_by_reciprocity(mesh.get_inner(), max_distance))
    }
}

// ================================================================================================
// PointCloudData3
// ================================================================================================

/// The unaccelerated point cloud container, holding the point buffer and its per-point attributes.
///
/// Numpy arrays handed out through the getters are cached and invalidated by any method which can
/// change what they hold, so repeated access does not repeatedly copy the buffer.
#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct PointCloudData3 {
    inner: engeom::PointCloudData3,
    points: Option<Py<PyArray2<f64>>>,
    point_normals: Option<Py<PyArray2<f64>>>,
    point_colors: Option<Py<PyArray2<u8>>>,
    point_stdev: Option<Py<PyArray1<f64>>>,
}

impl PointCloudData3 {
    fn clear_cached(&mut self) {
        self.points = None;
        self.point_normals = None;
        self.point_colors = None;
        self.point_stdev = None;
    }

    pub fn get_inner(&self) -> &engeom::PointCloudData3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::PointCloudData3) -> Self {
        Self {
            inner,
            points: None,
            point_normals: None,
            point_colors: None,
            point_stdev: None,
        }
    }
}

impl Clone for PointCloudData3 {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl PointCloudData3 {
    #[new]
    fn new<'py>(points: PyReadonlyArray2<'py, f64>) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        Ok(Self::from_inner(engeom::PointCloudData3::new(points)))
    }

    // --- Serialization ---------------------------------------------------------------------

    #[staticmethod]
    fn load_ply(path: PathBuf) -> PyResult<Self> {
        let inner = engeom::PointCloudData3::load_ply(&path)
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
        let inner = engeom::io::load_lptf3_point_data(&path, load)
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

    fn append_in_place(&mut self, other: &PointCloudData3) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .append_in_place(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn create_subset_indices(&self, indices: Vec<usize>) -> PyResult<Self> {
        let inner = self
            .inner
            .create_subset_indices(&indices)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    /// Build the queryable `PointCloud` from this data.
    ///
    /// Only the normals, colors, and standard deviations come across; `PointCloud` has nowhere to
    /// put the open-map attributes.
    fn to_cloud(&self) -> PyResult<PointCloud> {
        let inner = self
            .inner
            .to_cloud()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PointCloud::from_inner(inner))
    }

    /// Copy the buffer and attributes out of a queryable `PointCloud`.
    #[staticmethod]
    fn from_cloud(cloud: &PointCloud) -> PyResult<Self> {
        let inner = engeom::PointCloudData3::from_cloud(cloud.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn __len__(&self) -> usize {
        self.inner.point_count()
    }

    fn __repr__(&self) -> String {
        format!("<PointCloudData3 {} points>", self.inner.point_count())
    }
}
