use crate::bounding::Aabb3;
use crate::common::{IndexMask, deviation_mode_from_str, select_op_from_str};
use crate::conversions::{
    array_to_colors, array_to_faces, array_to_points3, array_to_unit_vectors3, array_to_vec,
    colors_to_array, faces_to_array, labels_to_array, points_to_array, scalars_to_array,
    unit_vectors_to_array,
};
use crate::geom3::{Curve3, Iso3, Plane3, Point3, SurfacePoint3, Vector3};
use crate::metrology::Distance3;
use crate::point_cloud::lptf3_load_from_args;
use engeom::Selection;
use engeom::common::SplitResult;
use engeom::geom3::align3::{GAPParams, generate_alignment_points};
use engeom::io::{deflate_bytes, u_bytes_to_mesh_data};
use numpy::ndarray::{Array1, Array2, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayDyn, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

/// The (n, 6) point-pair array and the (n,) per-edge type codes returned by
/// `Mesh3.compute_visual_outline`.
type VisualOutline<'py> = (Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArray1<u8>>);

#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct Mesh3 {
    inner: engeom::Mesh3,
    points: Option<Py<PyArray2<f64>>>,
    faces: Option<Py<PyArray2<u32>>>,
    face_normals: Option<Py<PyArray2<f64>>>,
    computed_point_normals: Option<Py<PyArray2<f64>>>,
    point_normals: Option<Py<PyArray2<f64>>>,
    point_colors: Option<Py<PyArray2<u8>>>,
    point_stdev: Option<Py<PyArray1<f64>>>,
    face_colors: Option<Py<PyArray2<u8>>>,
    face_labels: Option<Py<PyArray1<u32>>>,
}

impl Mesh3 {
    fn clear_cached(&mut self) {
        self.points = None;
        self.faces = None;
        self.face_normals = None;
        self.computed_point_normals = None;
        self.point_normals = None;
        self.point_colors = None;
        self.point_stdev = None;
        self.face_colors = None;
        self.face_labels = None;
    }

    pub fn get_inner(&self) -> &engeom::Mesh3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Mesh3) -> Self {
        Self {
            inner,
            points: None,
            faces: None,
            face_normals: None,
            computed_point_normals: None,
            point_normals: None,
            point_colors: None,
            point_stdev: None,
            face_colors: None,
            face_labels: None,
        }
    }
}

impl Clone for Mesh3 {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl Mesh3 {
    #[new]
    #[pyo3(signature=(points, faces, merge_duplicates = false, delete_degenerate = false, is_solid = false))]
    fn new<'py>(
        points: PyReadonlyArray2<'py, f64>,
        faces: PyReadonlyArray2<'py, u32>,
        merge_duplicates: bool,
        delete_degenerate: bool,
        is_solid: bool,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let faces = array_to_faces(&faces.as_array())?;
        let mesh = engeom::Mesh3::new_with_options(
            points,
            faces,
            is_solid,
            merge_duplicates,
            delete_degenerate,
            None,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(mesh))
    }

    #[getter]
    fn aabb(&self) -> Aabb3 {
        Aabb3::from_inner(self.inner.aabb())
    }

    #[staticmethod]
    #[pyo3(signature=(path, merge_duplicates = false, delete_degenerate = false, is_solid = false))]
    fn load_stl(
        path: PathBuf,
        merge_duplicates: bool,
        delete_degenerate: bool,
        is_solid: bool,
    ) -> PyResult<Self> {
        let mesh = engeom::Mesh3::load_stl(&path, is_solid, merge_duplicates, delete_degenerate)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    #[pyo3(signature=(path, is_solid = false))]
    fn load_ply(path: PathBuf, is_solid: bool) -> PyResult<Self> {
        let mesh = engeom::Mesh3::load_ply(&path, is_solid)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    #[pyo3(signature=(path, is_solid = false))]
    fn load_g3d(path: PathBuf, is_solid: bool) -> PyResult<Self> {
        let mesh = engeom::Mesh3::load_g3d(&path, is_solid)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(mesh))
    }

    #[pyo3(signature = (path, binary = true))]
    fn save_ply(&self, path: PathBuf, binary: bool) -> PyResult<()> {
        let mut opts = engeom::io::PlyWriteOpts::default();
        opts.binary = binary;
        self.inner
            .save_ply(&path, &opts)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    fn transform_in_place(&mut self, iso: &Iso3) {
        self.inner.transform_in_place(iso.get_inner());
        self.clear_cached()
    }

    fn transform_copy(&self, iso: &Iso3) -> Self {
        Self::from_inner(self.inner.transform_copy(iso.get_inner()))
    }

    fn surface_closest_to(&self, x: f64, y: f64, z: f64) -> SurfacePoint3 {
        let p = engeom::Point3::new(x, y, z);
        SurfacePoint3::from_inner(self.inner.surface_closest_to(&p).sp)
    }

    fn barycentric_closest_to(&self, x: f64, y: f64, z: f64) -> (u32, [f64; 3]) {
        let p = engeom::Point3::new(x, y, z);
        let msp = self.inner.surface_closest_to(&p);
        (msp.face_index, msp.bc)
    }

    fn point_closest_to(&self, x: f64, y: f64, z: f64) -> Point3 {
        let p = engeom::Point3::new(x, y, z);
        Point3::from_inner(self.inner.point_closest_to(&p))
    }

    fn distance_closest_to(&self, x: f64, y: f64, z: f64) -> f64 {
        let p = engeom::Point3::new(x, y, z);
        self.inner.distance_closest_to(&p)
    }

    fn face_closest_to(&self, x: f64, y: f64, z: f64) -> u32 {
        let p = engeom::Point3::new(x, y, z);
        self.inner.face_closest_to(&p)
    }

    fn append_in_place(&mut self, other: &Mesh3) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .append_in_place(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    #[pyo3(signature = (path, binary = true, allow_attribute_loss = false))]
    fn save_stl(&self, path: PathBuf, binary: bool, allow_attribute_loss: bool) -> PyResult<()> {
        let mut opts = engeom::io::StlWriteOpts::default();
        opts.binary = binary;
        opts.allow_attribute_loss = allow_attribute_loss;
        self.inner
            .save_stl(&path, &opts)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    #[staticmethod]
    #[pyo3(signature = (path, is_solid = false))]
    fn load_tcmesh(path: PathBuf, is_solid: bool) -> PyResult<Self> {
        let data =
            engeom::io::read_tc_mesh_file(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        let mesh = engeom::Mesh3::from_data(data, is_solid)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(mesh))
    }

    // `allow_attribute_loss` is gone rather than kept as a parameter that only accepts one value:
    // the tcmesh format refuses an attributed mesh outright now.
    #[pyo3(signature = (path, tol))]
    fn save_tcmesh(&self, path: PathBuf, tol: f64) -> PyResult<()> {
        engeom::io::write_tc_mesh_file(&path, &self.inner.to_data(), tol)
            .map_err(|e| PyIOError::new_err(e.to_string()))
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
    fn is_solid(&self) -> bool {
        self.inner.is_solid()
    }

    #[getter]
    fn point_count(&self) -> usize {
        self.inner.point_count()
    }

    #[getter]
    fn face_count(&self) -> usize {
        self.inner.face_count()
    }

    fn __len__(&self) -> usize {
        self.inner.point_count()
    }

    fn flip_normals_in_place(&mut self) {
        self.clear_cached();
        self.inner.flip_normals_in_place();
    }

    // --- Attributes ------------------------------------------------------------------------

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
    fn face_colors<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<u8>>> {
        let values = self.inner.face_colors()?;
        if self.face_colors.is_none() {
            let array = colors_to_array(values);
            self.face_colors = Some(array.into_pyarray(py).unbind());
        }
        Some(self.face_colors.as_ref().unwrap().bind(py))
    }

    #[getter]
    fn face_labels<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray1<u32>>> {
        let values = self.inner.face_labels()?;
        if self.face_labels.is_none() {
            let array = labels_to_array(values);
            self.face_labels = Some(array.into_pyarray(py).unbind());
        }
        Some(self.face_labels.as_ref().unwrap().bind(py))
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
    fn set_face_colors<'py>(&mut self, values: Option<PyReadonlyArray2<'py, u8>>) -> PyResult<()> {
        let values = values.map(|v| array_to_colors(&v.as_array())).transpose()?;
        self.clear_cached();
        self.inner
            .set_face_colors(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (values=None))]
    fn set_face_labels<'py>(&mut self, values: Option<PyReadonlyArray1<'py, u32>>) -> PyResult<()> {
        let values = values.map(|v| array_to_vec(&v)).transpose()?;
        self.clear_cached();
        self.inner
            .set_face_labels(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn compute_point_normals<'py>(
        &mut self,
        py: Python<'py>,
    ) -> PyResult<&Bound<'py, PyArray2<f64>>> {
        if self.computed_point_normals.is_none() {
            let normals = self
                .inner
                .compute_point_normals()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            let array = unit_vectors_to_array(&normals);
            self.computed_point_normals = Some(array.into_pyarray(py).unbind());
        }

        Ok(self.computed_point_normals.as_ref().unwrap().bind(py))
    }

    fn scale_copy(&self, scale: f64) -> PyResult<Self> {
        let inner = self
            .inner
            .scale_copy(scale)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn offset_points_copy(&self, offset: f64) -> PyResult<Self> {
        let inner = self
            .inner
            .offset_points_copy(offset)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn patch_boundaries(&self) -> PyResult<Vec<Curve3>> {
        // TODO: Is this actually used? What's the difference between it and `boundary_curves`?
        let boundaries = self
            .inner
            .compute_patch_boundary_points()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let mut result = Vec::new();
        for b in boundaries.iter() {
            let c = engeom::Curve3::from_points(b, 1.0e-6)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            result.push(Curve3::from_inner(c))
        }

        Ok(result)
    }

    #[pyo3(signature=(facing, max_edge_length, corner_angle = None))]
    fn compute_visual_outline<'py>(
        &self,
        py: Python<'py>,
        facing: Vector3,
        max_edge_length: f64,
        corner_angle: Option<f64>,
    ) -> PyResult<VisualOutline<'py>> {
        let n = engeom::UnitVec3::new_normalize(*facing.get_inner());
        let outline = self
            .inner
            .compute_visual_outline(n, max_edge_length, corner_angle)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut result = ArrayD::zeros(vec![outline.len(), 6]);
        let mut result_type = Array1::zeros(outline.len());
        for (i, (p0, p1, t)) in outline.iter().enumerate() {
            result[[i, 0]] = p0.x;
            result[[i, 1]] = p0.y;
            result[[i, 2]] = p0.z;
            result[[i, 3]] = p1.x;
            result[[i, 4]] = p1.y;
            result[[i, 5]] = p1.z;

            result_type[i] = *t;
        }
        Ok((result.into_pyarray(py), result_type.into_pyarray(py)))
    }

    /// A method rather than a property, matching `compute_point_normals` and the core name. Both
    /// are derived quantities which can fail, and a property which raises reads badly.
    fn compute_face_normals<'py>(
        &mut self,
        py: Python<'py>,
    ) -> PyResult<&Bound<'py, PyArray2<f64>>> {
        if self.face_normals.is_none() {
            let normals = self
                .inner
                .compute_face_normals()
                .map_err(|e| PyValueError::new_err(e.to_string()))?;

            let array = unit_vectors_to_array(&normals);
            self.face_normals = Some(array.into_pyarray(py).unbind());
        }

        Ok(self.face_normals.as_ref().unwrap().bind(py))
    }

    fn compute_face_areas<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let areas = self
            .inner
            .compute_face_areas()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(scalars_to_array(&areas).into_pyarray(py))
    }

    fn compute_face_centers<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let centers = self
            .inner
            .compute_face_centers()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(points_to_array(&centers).into_pyarray(py))
    }

    fn boundary_curves(&self) -> PyResult<Vec<Curve3>> {
        let edges = self
            .inner
            .compute_edges()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let mut curves = Vec::new();
        for item in edges.boundary_loops.iter() {
            let mut points = Vec::new();
            for &index in item.iter() {
                let point = self.inner.points()[index as usize];
                points.push(point);
            }

            let c = engeom::Curve3::from_points(&points, 1.0e-6)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;

            curves.push(Curve3::from_inner(c));
        }
        Ok(curves)
    }

    #[getter]
    fn faces<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<u32>> {
        if self.faces.is_none() {
            let faces = faces_to_array(self.inner.faces());
            self.faces = Some(faces.into_pyarray(py).unbind());
        }

        self.faces.as_ref().unwrap().bind(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "<Mesh3 {} vertices, {} faces>",
            self.inner.points().len(),
            self.inner.faces().len()
        )
    }

    #[pyo3(signature = (plane, allow_attribute_loss = false))]
    fn split(
        &self,
        plane: &Plane3,
        allow_attribute_loss: bool,
    ) -> PyResult<(Option<Self>, Option<Self>)> {
        let result = self
            .inner
            .split(&plane.inner, allow_attribute_loss)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        match result {
            SplitResult::Pair(mesh1, mesh2) => {
                Ok((Some(Self::from_inner(mesh1)), Some(Self::from_inner(mesh2))))
            }
            SplitResult::Negative => Ok((Some(self.clone()), None)),
            SplitResult::Positive => Ok((None, Some(self.clone()))),
        }
    }

    fn measure_deviations<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        mode: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let mode = deviation_mode_from_str(mode)?;
        let points = array_to_points3(&points.as_array())?;
        let values = self.inner.measure_deviations(&points, mode);
        Ok(scalars_to_array(&values).into_pyarray(py))
    }

    fn measure_point_deviation(
        &self,
        x: f64,
        y: f64,
        z: f64,
        dist_mode: &str,
    ) -> PyResult<Distance3> {
        let point = engeom::Point3::new(x, y, z);
        Ok(Distance3::from_inner(self.inner.measure_point_deviation(
            &point,
            deviation_mode_from_str(dist_mode)?,
        )))
    }

    fn boundary_first_flatten<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let edges = self
            .inner
            .compute_edges()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let values = edges
            .boundary_first_flatten()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let result = points_to_array(&values);
        Ok(result.into_pyarray(py))
    }

    fn sample_poisson<'py>(&self, py: Python<'py>, radius: f64) -> Bound<'py, PyArray2<f64>> {
        let mps = self.inner.sample_poisson(radius, None);
        let mut result = Array2::zeros((mps.len(), 6));
        for (i, mp) in mps.iter().enumerate() {
            result[[i, 0]] = mp.sp.point.x;
            result[[i, 1]] = mp.sp.point.y;
            result[[i, 2]] = mp.sp.point.z;
            result[[i, 3]] = mp.sp.normal.x;
            result[[i, 4]] = mp.sp.normal.y;
            result[[i, 5]] = mp.sp.normal.z;
        }
        result.into_pyarray(py)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_alignment_points<'py>(
        &self,
        py: Python<'py>,
        reference: &Mesh3,
        iso: Iso3,
        max_spacing: f64,
        max_neighbor_angle: f64,       // PI / 3.0
        out_of_plane_ratio: f64,       // 1 /20.0
        centroid_ratio: f64,           // 1.0
        filter_distances: Option<f64>, // Some(3.0)
    ) -> Bound<'py, PyArray2<f64>> {
        let params = GAPParams::new(
            max_spacing,
            max_neighbor_angle,
            out_of_plane_ratio,
            centroid_ratio,
            filter_distances,
        );
        let mps =
            generate_alignment_points(&self.inner, &reference.inner, iso.get_inner(), &params);
        let points = mps.into_iter().map(|mp| mp.sp.point).collect::<Vec<_>>();

        let result = points_to_array(&points);
        result.into_pyarray(py)
    }

    #[pyo3(signature=(plane, tol = None, faces = None))]
    fn section_with_plane(
        &self,
        plane: Plane3,
        tol: Option<f64>,
        faces: Option<&IndexMask>,
    ) -> PyResult<Vec<Curve3>> {
        let results = self
            .inner
            .section_with_plane(plane.get_inner(), tol, faces.map(|m| m.get_inner()))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(results.into_iter().map(Curve3::from_inner).collect())
    }

    /// Start a face filter. `start` is either the token `"none"` or `"all"`, or an `IndexMask` over
    /// the faces, covering the three ways the core `Selection` enum can seed a filter.
    #[pyo3(signature = (start = None))]
    fn face_select<'py>(
        slf: PyRef<Self>,
        py: Python<'py>,
        start: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, FaceFilterHandle>> {
        let selection = match start {
            None => Selection::None,
            Some(value) => {
                if let Ok(mask) = value.extract::<PyRef<'_, IndexMask>>() {
                    Selection::Mask(mask.get_inner().clone())
                } else {
                    match value.extract::<String>()?.as_str() {
                        "none" => Selection::None,
                        "all" => Selection::All,
                        other => {
                            return Err(PyValueError::new_err(format!(
                                "Invalid face selection start '{other}', expected 'none', 'all', \
                                 or an IndexMask"
                            )));
                        }
                    }
                }
            }
        };

        // A mask of the wrong length would otherwise reach the core and be zipped against the face
        // buffer, silently selecting the wrong faces rather than complaining.
        if let Selection::Mask(mask) = &selection
            && mask.len() != slf.inner.face_count()
        {
            return Err(PyValueError::new_err(format!(
                "Mask length {} does not match the mesh's face count {}",
                mask.len(),
                slf.inner.face_count()
            )));
        }

        let indices = slf.inner.face_select(selection).collect_indices();
        FaceFilterHandle {
            mesh: slf.into(),
            indices,
        }
        .into_pyobject(py)
    }

    fn extract_subset_faces_from_indices(&self, indices: Vec<usize>) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_faces_from_indices(&indices)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    // --- Masks -------------------------------------------------------------------------------

    #[pyo3(signature = (value = false))]
    fn face_mask(&self, value: bool) -> IndexMask {
        IndexMask::from_inner(self.inner.face_mask(value))
    }

    #[pyo3(signature = (value = false))]
    fn point_mask(&self, value: bool) -> IndexMask {
        IndexMask::from_inner(self.inner.point_mask(value))
    }

    fn extract_subset_faces(&self, mask: &IndexMask) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_faces(mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn compute_unique_point_mask(&self, face_mask: &IndexMask) -> PyResult<IndexMask> {
        let inner = self
            .inner
            .compute_unique_point_mask(face_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(IndexMask::from_inner(inner))
    }

    #[pyo3(signature = (mask = None))]
    fn compute_patches(&self, mask: Option<&IndexMask>) -> PyResult<Vec<IndexMask>> {
        let patches = self
            .inner
            .compute_patches(mask.map(|m| m.get_inner()))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(patches.into_iter().map(IndexMask::from_inner).collect())
    }

    #[pyo3(signature = (points, max_dist, max_angle, transform = None))]
    fn find_points_in_tol<'py>(
        &self,
        points: PyReadonlyArray2<'py, f64>,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> PyResult<IndexMask> {
        let points = array_to_points3(&points.as_array())?;
        let mask = self.inner.find_points_in_tol(
            &points,
            max_dist,
            max_angle,
            transform.map(|t| t.get_inner()),
        );
        Ok(IndexMask::from_inner(mask))
    }

    fn separate_patches(&self) -> PyResult<Vec<Self>> {
        let patch_groups = self
            .inner
            .compute_patches(None)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut results = Vec::with_capacity(patch_groups.len());
        for mask in patch_groups.iter() {
            results.push(
                self.inner
                    .extract_subset_faces(mask)
                    .map_err(|e| PyValueError::new_err(e.to_string()))?,
            );
        }
        Ok(results.into_iter().map(Self::from_inner).collect())
    }

    #[pyo3(signature = (allow_attribute_loss = false))]
    fn convex_hull(&self, allow_attribute_loss: bool) -> PyResult<Self> {
        let inner = self
            .inner
            .convex_hull(allow_attribute_loss)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[staticmethod]
    #[pyo3(signature = (length, width, height, is_solid = true))]
    fn create_box(length: f64, width: f64, height: f64, is_solid: bool) -> Self {
        let mesh = engeom::Mesh3::create_box(length, width, height, is_solid);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        let mesh = engeom::Mesh3::create_cylinder(radius, height, steps);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let mesh = engeom::Mesh3::create_sphere(radius, n_theta, n_phi);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_cone(radius: f64, height: f64, steps: usize) -> Self {
        let mesh = engeom::Mesh3::create_cone(radius, height, steps);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_circle(radius: f64, segments: usize) -> Self {
        let mesh = engeom::Mesh3::create_circle(radius, segments);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_capsule(p0: Point3, p1: Point3, radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let mesh =
            engeom::Mesh3::create_capsule(p0.get_inner(), p1.get_inner(), radius, n_theta, n_phi);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    #[pyo3(signature=(p0, p1, width, height, up=None))]
    fn create_rect_beam_between(
        p0: Point3,
        p1: Point3,
        width: f64,
        height: f64,
        up: Option<Vector3>,
    ) -> PyResult<Self> {
        let up = up.map_or(engeom::Vector3::z(), |v| *v.get_inner());
        let mesh = engeom::Mesh3::create_rect_beam_between(
            p0.get_inner(),
            p1.get_inner(),
            width,
            height,
            &up,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    fn create_cylinder_between(
        p0: Point3,
        p1: Point3,
        radius: f64,
        steps: usize,
    ) -> PyResult<Self> {
        let mesh =
            engeom::Mesh3::create_cylinder_between(p0.get_inner(), p1.get_inner(), radius, steps);
        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    #[pyo3(signature = (file_path, take_every=1, look_scale=None, weight_scale=None, max_move=None))]
    fn load_lptf3(
        file_path: PathBuf,
        take_every: u32,
        look_scale: Option<f64>,
        weight_scale: Option<f64>,
        max_move: Option<f64>,
    ) -> PyResult<Mesh3> {
        let load = lptf3_load_from_args(take_every, look_scale, weight_scale, max_move)?;
        let mesh_data = engeom::io::load_lptf3_mesh_data(&file_path, load, None)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        let (points, faces, _) = mesh_data.into_parts();
        Ok(Self::from_inner(engeom::Mesh3::new(points, faces, false)))
    }

    #[staticmethod]
    fn load_umesh(file_path: PathBuf) -> PyResult<Mesh3> {
        // These files are always small, so we're going to just pull the whole thing into memory
        let file_bytes =
            std::fs::read(&file_path).map_err(|e| PyIOError::new_err(e.to_string()))?;

        let deflated = if let Ok(b) = deflate_bytes(&file_bytes) {
            b
        } else {
            file_bytes
        };

        let data =
            u_bytes_to_mesh_data(&deflated).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mesh = engeom::Mesh3::from_data(data, false)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    fn stanford_bunny_res4() -> Self {
        let mesh = engeom::Mesh3::stanford_bunny_res4();
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn stanford_bunny_res3() -> Self {
        let mesh = engeom::Mesh3::stanford_bunny_res3();
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn stanford_bunny_res2() -> Self {
        let mesh = engeom::Mesh3::stanford_bunny_res2();
        Self::from_inner(mesh)
    }
}

#[pyclass(module = "engeom.geom3")]
pub struct FaceFilterHandle {
    mesh: Py<Mesh3>,
    indices: Vec<usize>,
}

#[pymethods]
impl FaceFilterHandle {
    fn __repr__(&self) -> String {
        format!("<FaceFilterHandle {} triangles>", self.indices.len())
    }

    #[allow(clippy::too_many_arguments)]
    fn facing<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        x: f64,
        y: f64,
        z: f64,
        angle: f64,
        mode: &str,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let normal = engeom::UnitVec3::new_normalize([x, y, z].into());
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .facing(&normal, angle, op)
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(other, all_points, distance_tol, mode, planar_tol = None, angle_tol = None))]
    fn near_mesh<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        other: PyRef<Mesh3>,
        all_points: bool,
        distance_tol: f64,
        mode: &str,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .near_mesh(
                &other.inner,
                all_points,
                distance_tol,
                planar_tol,
                angle_tol,
                op,
            )
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn above_plane<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        plane: &Plane3,
        all_vertices: bool,
        mode: &str,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .above_plane(plane.get_inner(), all_vertices, op)
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[allow(clippy::too_many_arguments)]
    fn vertices_near_point<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        x: f64,
        y: f64,
        z: f64,
        max_dist: f64,
        all_vertices: bool,
        mode: &str,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let point = engeom::Point3::new(x, y, z);
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .vertices_near_point(&point, max_dist, all_vertices, op)
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (mode, exclude = None))]
    fn expand<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        mode: &str,
        exclude: Option<&IndexMask>,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .expand(exclude.map(|m| m.get_inner()), op)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (n, mode, exclude = None))]
    fn expand_n<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        n: usize,
        mode: &str,
        exclude: Option<&IndexMask>,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .expand_n(n, exclude.map(|m| m.get_inner()), op)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[allow(clippy::too_many_arguments)]
    fn faces_overlap<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        other: PyRef<Mesh3>,
        angle_tol: f64,
        distance_tol: f64,
        mode: &str,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .faces_overlap(&other.inner, angle_tol, distance_tol, op)
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn by_mask<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        mask: &IndexMask,
        mode: &str,
    ) -> PyResult<Bound<'py, Self>> {
        let op = select_op_from_str(mode)?;
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .by_mask(mask.get_inner(), op)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn collect_indices(&self) -> Vec<usize> {
        self.indices.clone()
    }

    /// Named `to_mask` rather than mirroring the core's `take_mask`, and `to_mesh` below rather
    /// than the core's `into_mesh`, because both Rust methods consume the filter and neither of
    /// these does: the handle stays usable afterwards.
    fn to_mask(&self, py: Python<'_>) -> PyResult<IndexMask> {
        let mesh = self.mesh.bind(py).borrow();
        let inner =
            engeom::common::IndexMask::try_from_indices(&self.indices, mesh.inner.face_count())
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(IndexMask::from_inner(inner))
    }

    fn to_mesh(&self, py: Python<'_>) -> PyResult<Mesh3> {
        self.mesh
            .bind(py)
            .borrow()
            .extract_subset_faces_from_indices(self.indices.clone())
    }
}

#[pyclass(module = "engeom.geom3")]
pub struct MeshCollisionSet {
    inner: engeom::geom3::MeshCollisionSet,
}

impl MeshCollisionSet {
    pub fn get_inner(&self) -> &engeom::geom3::MeshCollisionSet {
        &self.inner
    }

    pub fn from_inner(inner: engeom::geom3::MeshCollisionSet) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl MeshCollisionSet {
    #[new]
    fn new() -> Self {
        Self::from_inner(engeom::geom3::MeshCollisionSet::new())
    }

    fn add_stationary(&mut self, mesh: &Mesh3) -> usize {
        let inner = mesh.inner.clone();
        self.inner.add_stationary(inner)
    }

    fn add_moving(&mut self, mesh: &Mesh3) -> usize {
        let inner = mesh.inner.clone();
        self.inner.add_moving(inner)
    }

    fn add_exception(&mut self, id1: usize, id2: usize) {
        self.inner.add_exception(id1, id2);
    }

    fn check_all(
        &self,
        transforms: Vec<(usize, Iso3)>,
        stop_at_first: bool,
    ) -> PyResult<Vec<(usize, usize)>> {
        let transforms = transforms
            .into_iter()
            .map(|(id, iso)| (id, *iso.get_inner()))
            .collect::<Vec<_>>();

        let result = self
            .inner
            .check_all(&transforms, stop_at_first)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(result)
    }
}

// ================================================================================================
// MeshData3
// ================================================================================================

/// The unaccelerated mesh container, holding the point and face buffers and their per-element
/// attributes.
///
/// Numpy arrays handed out through the getters are cached and invalidated by any method which can
/// change what they hold, so repeated access does not repeatedly copy the buffers.
#[pyclass(from_py_object, module = "engeom.geom3")]
pub struct MeshData3 {
    inner: engeom::MeshData3,
    points: Option<Py<PyArray2<f64>>>,
    faces: Option<Py<PyArray2<u32>>>,
    point_normals: Option<Py<PyArray2<f64>>>,
    point_colors: Option<Py<PyArray2<u8>>>,
    point_stdev: Option<Py<PyArray1<f64>>>,
    face_colors: Option<Py<PyArray2<u8>>>,
    face_labels: Option<Py<PyArray1<u32>>>,
}

impl MeshData3 {
    fn clear_cached(&mut self) {
        self.points = None;
        self.faces = None;
        self.point_normals = None;
        self.point_colors = None;
        self.point_stdev = None;
        self.face_colors = None;
        self.face_labels = None;
    }

    pub fn get_inner(&self) -> &engeom::MeshData3 {
        &self.inner
    }

    pub fn from_inner(inner: engeom::MeshData3) -> Self {
        Self {
            inner,
            points: None,
            faces: None,
            point_normals: None,
            point_colors: None,
            point_stdev: None,
            face_colors: None,
            face_labels: None,
        }
    }
}

impl Clone for MeshData3 {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl MeshData3 {
    #[new]
    fn new<'py>(
        points: PyReadonlyArray2<'py, f64>,
        faces: PyReadonlyArray2<'py, u32>,
    ) -> PyResult<Self> {
        let points = array_to_points3(&points.as_array())?;
        let faces = array_to_faces(&faces.as_array())?;
        let inner = engeom::MeshData3::new(points, faces)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    /// Create mesh data with no points, no faces, and no attributes, as a starting point for
    /// building one up with `push_point` and `push_face`.
    #[staticmethod]
    fn empty() -> Self {
        Self::from_inner(engeom::MeshData3::empty())
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[getter]
    fn point_count(&self) -> usize {
        self.inner.point_count()
    }

    #[getter]
    fn face_count(&self) -> usize {
        self.inner.face_count()
    }

    // --- Serialization ---------------------------------------------------------------------

    #[staticmethod]
    fn load_ply(path: PathBuf) -> PyResult<Self> {
        let inner =
            engeom::MeshData3::load_ply(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
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
    fn load_stl(path: PathBuf) -> PyResult<Self> {
        let inner =
            engeom::MeshData3::load_stl(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[staticmethod]
    fn load_g3d(path: PathBuf) -> PyResult<Self> {
        let inner =
            engeom::MeshData3::load_g3d(&path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[pyo3(signature = (path, binary = true, allow_attribute_loss = false))]
    fn save_stl(&self, path: PathBuf, binary: bool, allow_attribute_loss: bool) -> PyResult<()> {
        let mut opts = engeom::io::StlWriteOpts::default();
        opts.binary = binary;
        opts.allow_attribute_loss = allow_attribute_loss;
        self.inner
            .save_stl(&path, &opts)
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
        let inner = engeom::io::load_lptf3_mesh_data(&path, load, None)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    // --- Buffers ---------------------------------------------------------------------------

    #[getter]
    fn points<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.points.is_none() {
            let array = points_to_array(self.inner.points());
            self.points = Some(array.into_pyarray(py).unbind());
        }
        self.points.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn faces<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<u32>> {
        if self.faces.is_none() {
            let array = faces_to_array(self.inner.faces());
            self.faces = Some(array.into_pyarray(py).unbind());
        }
        self.faces.as_ref().unwrap().bind(py)
    }

    // --- Attributes ------------------------------------------------------------------------

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
    fn face_colors<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray2<u8>>> {
        let values = self.inner.face_colors()?;
        if self.face_colors.is_none() {
            let array = colors_to_array(values);
            self.face_colors = Some(array.into_pyarray(py).unbind());
        }
        Some(self.face_colors.as_ref().unwrap().bind(py))
    }

    #[getter]
    fn face_labels<'py>(&mut self, py: Python<'py>) -> Option<&Bound<'py, PyArray1<u32>>> {
        let values = self.inner.face_labels()?;
        if self.face_labels.is_none() {
            let array = labels_to_array(values);
            self.face_labels = Some(array.into_pyarray(py).unbind());
        }
        Some(self.face_labels.as_ref().unwrap().bind(py))
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
    fn set_face_colors<'py>(&mut self, values: Option<PyReadonlyArray2<'py, u8>>) -> PyResult<()> {
        let values = values.map(|v| array_to_colors(&v.as_array())).transpose()?;
        self.clear_cached();
        self.inner
            .set_face_colors(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (values=None))]
    fn set_face_labels<'py>(&mut self, values: Option<PyReadonlyArray1<'py, u32>>) -> PyResult<()> {
        let values = values.map(|v| array_to_vec(&v)).transpose()?;
        self.clear_cached();
        self.inner
            .set_face_labels(values)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    // --- Editing ---------------------------------------------------------------------------

    fn push_point(&mut self, point: Point3) -> PyResult<u32> {
        self.clear_cached();
        self.inner
            .push_point(*point.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn push_face(&mut self, face: [u32; 3]) -> PyResult<u32> {
        self.clear_cached();
        self.inner
            .push_face(face)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn set_point(&mut self, index: u32, point: Point3) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .set_point(index, *point.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn set_face(&mut self, index: u32, face: [u32; 3]) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .set_face(index, face)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn remove_point(&mut self, index: u32) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .remove_point(index)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn remove_face(&mut self, index: u32) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .remove_face(index)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    // --- Computed quantities -----------------------------------------------------------------
    //
    // Unlike their counterparts on `Mesh3`, these are not cached. This is the editable type, and a
    // cached derived quantity is one more thing every mutation has to remember to invalidate.

    fn compute_point_normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let normals = self
            .inner
            .compute_point_normals()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(unit_vectors_to_array(&normals).into_pyarray(py))
    }

    fn compute_face_normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let normals = self
            .inner
            .compute_face_normals()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(unit_vectors_to_array(&normals).into_pyarray(py))
    }

    fn compute_face_areas<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let areas = self
            .inner
            .compute_face_areas()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(scalars_to_array(&areas).into_pyarray(py))
    }

    fn compute_face_centers<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let centers = self
            .inner
            .compute_face_centers()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(points_to_array(&centers).into_pyarray(py))
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

    fn flip_faces_in_place(&mut self) {
        self.clear_cached();
        self.inner.flip_faces_in_place();
    }

    fn append_in_place(&mut self, other: &MeshData3) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .append_in_place(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    // --- Masks -------------------------------------------------------------------------------

    fn extract_subset_faces(&self, face_mask: &IndexMask) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_faces(face_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn extract_subset_points(&self, point_mask: &IndexMask) -> PyResult<Self> {
        let inner = self
            .inner
            .extract_subset_points(point_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    fn compute_unique_point_mask(&self, face_mask: &IndexMask) -> PyResult<IndexMask> {
        let inner = self
            .inner
            .compute_unique_point_mask(face_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(IndexMask::from_inner(inner))
    }

    fn compute_unique_face_mask(&self, point_mask: &IndexMask) -> PyResult<IndexMask> {
        let inner = self
            .inner
            .compute_unique_face_mask(point_mask.get_inner())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(IndexMask::from_inner(inner))
    }

    /// Build the accelerated `Mesh3` from this data, carrying every attribute across.
    #[pyo3(signature = (is_solid = false))]
    fn to_mesh(&self, is_solid: bool) -> PyResult<Mesh3> {
        let inner = engeom::Mesh3::from_data(self.inner.clone(), is_solid)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Mesh3::from_inner(inner))
    }

    /// Copy the buffers and attributes out of an accelerated `Mesh3`.
    #[staticmethod]
    fn from_mesh(mesh: &Mesh3) -> Self {
        Self::from_inner(mesh.get_inner().to_data())
    }

    // --- Primitives --------------------------------------------------------------------------

    #[staticmethod]
    fn create_box(length: f64, width: f64, height: f64) -> Self {
        Self::from_inner(engeom::MeshData3::create_box(length, width, height))
    }

    #[staticmethod]
    fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_sphere(radius, n_theta, n_phi))
    }

    #[staticmethod]
    fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_cylinder(radius, height, steps))
    }

    #[staticmethod]
    fn create_cone(radius: f64, height: f64, steps: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_cone(radius, height, steps))
    }

    #[staticmethod]
    fn create_circle(radius: f64, segments: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_circle(radius, segments))
    }

    #[staticmethod]
    fn create_capsule(p0: Point3, p1: Point3, radius: f64, n_theta: usize, n_phi: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_capsule(
            p0.get_inner(),
            p1.get_inner(),
            radius,
            n_theta,
            n_phi,
        ))
    }

    #[staticmethod]
    fn create_cylinder_between(p0: Point3, p1: Point3, radius: f64, steps: usize) -> Self {
        Self::from_inner(engeom::MeshData3::create_cylinder_between(
            p0.get_inner(),
            p1.get_inner(),
            radius,
            steps,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (p0, p1, width, height, up = None))]
    fn create_rect_beam_between(
        p0: Point3,
        p1: Point3,
        width: f64,
        height: f64,
        up: Option<Vector3>,
    ) -> PyResult<Self> {
        let up = up.map_or(engeom::Vector3::z(), |v| *v.get_inner());
        let inner = engeom::MeshData3::create_rect_beam_between(
            p0.get_inner(),
            p1.get_inner(),
            width,
            height,
            &up,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self::from_inner(inner))
    }

    #[staticmethod]
    fn stanford_bunny_res4() -> Self {
        Self::from_inner(engeom::MeshData3::stanford_bunny_res4())
    }

    #[staticmethod]
    fn stanford_bunny_res3() -> Self {
        Self::from_inner(engeom::MeshData3::stanford_bunny_res3())
    }

    #[staticmethod]
    fn stanford_bunny_res2() -> Self {
        Self::from_inner(engeom::MeshData3::stanford_bunny_res2())
    }

    fn __len__(&self) -> usize {
        self.inner.point_count()
    }

    fn __repr__(&self) -> String {
        format!(
            "<MeshData3 {} points, {} faces>",
            self.inner.point_count(),
            self.inner.face_count()
        )
    }
}
