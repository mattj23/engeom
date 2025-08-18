use crate::bounding::Aabb3;
use crate::common::{DeviationMode, SelectOp};
use crate::conversions::{
    array_to_faces, array_to_points3, faces_to_array, points_to_array, vectors_to_array,
};
use crate::geom3::{Curve3, Iso3, Plane3, Point3, SurfacePoint3, Vector3};
use crate::metrology::Distance3;
use crate::point_cloud::Lptf3Load;
use engeom::common::points::dist;
use engeom::common::{Selection, SplitResult};
use engeom::geom3::align3::{GAPParams, generate_alignment_points};
use numpy::ndarray::{Array1, Array2, ArrayD};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayDyn, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyclass]
pub struct Mesh {
    inner: engeom::Mesh,
    vertices: Option<Py<PyArray2<f64>>>,
    faces: Option<Py<PyArray2<u32>>>,
    face_normals: Option<Py<PyArray2<f64>>>,
    vertex_normals: Option<Py<PyArray2<f64>>>,
}

impl Mesh {
    fn clear_cached(&mut self) {
        self.vertices = None;
        self.faces = None;
        self.face_normals = None;
        self.vertex_normals = None;
    }

    pub fn get_inner(&self) -> &engeom::Mesh {
        &self.inner
    }

    pub fn from_inner(inner: engeom::Mesh) -> Self {
        Self {
            inner,
            vertices: None,
            faces: None,
            face_normals: None,
            vertex_normals: None,
        }
    }
}

impl Clone for Mesh {
    fn clone(&self) -> Self {
        Self::from_inner(self.inner.clone())
    }
}

#[pymethods]
impl Mesh {
    #[new]
    #[pyo3(signature=(vertices, faces, merge_duplicates = false, delete_degenerate = false))]
    fn new<'py>(
        vertices: PyReadonlyArray2<'py, f64>,
        faces: PyReadonlyArray2<'py, u32>,
        merge_duplicates: bool,
        delete_degenerate: bool,
    ) -> PyResult<Self> {
        let vertices = array_to_points3(&vertices.as_array())?;
        let faces = array_to_faces(&faces.as_array())?;
        let mesh = engeom::Mesh::new_with_options(
            vertices,
            faces,
            false,
            merge_duplicates,
            delete_degenerate,
            None,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(mesh))
    }

    #[getter]
    fn aabb(&self) -> Aabb3 {
        Aabb3::from_inner(*self.inner.aabb())
    }

    #[staticmethod]
    #[pyo3(signature=(path, merge_duplicates = false, delete_degenerate = false))]
    fn load_stl(path: PathBuf, merge_duplicates: bool, delete_degenerate: bool) -> PyResult<Self> {
        let mesh = engeom::io::read_mesh_stl(&path, merge_duplicates, delete_degenerate)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        Ok(Self::from_inner(mesh))
    }

    fn transform_by(&mut self, iso: &Iso3) {
        self.inner.transform_by(iso.get_inner());

        self.clear_cached()
    }

    fn surface_closest_to(&self, x: f64, y: f64, z: f64) -> SurfacePoint3 {
        let p = engeom::Point3::new(x, y, z);
        SurfacePoint3::from_inner(self.inner.surf_closest_to(&p).sp)
    }

    fn point_closest_to(&self, x: f64, y: f64, z: f64) -> Point3 {
        let p = engeom::Point3::new(x, y, z);
        Point3::from_inner(self.inner.point_closest_to(&p))
    }

    fn append(&mut self, other: &Mesh) -> PyResult<()> {
        self.clear_cached();
        self.inner
            .append(&other.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn cloned(&self) -> Self {
        self.clone()
    }

    fn write_stl(&self, path: PathBuf) -> PyResult<()> {
        engeom::io::write_mesh_stl(&path, &self.inner)
            .map_err(|e| PyIOError::new_err(e.to_string()))
    }

    #[getter]
    fn vertices<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.vertices.is_none() {
            let array = points_to_array(self.inner.vertices());
            self.vertices = Some(array.into_pyarray(py).unbind());
        }
        self.vertices.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn vertex_normals<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.vertex_normals.is_none() {
            let normals = self.inner.get_vertex_normals();
            let array = vectors_to_array(&normals);
            self.vertex_normals = Some(array.into_pyarray(py).unbind());
        }

        self.vertex_normals.as_ref().unwrap().bind(py)
    }

    fn get_patch_boundaries(&self) -> PyResult<Vec<Curve3>> {
        let boundaries = self
            .inner
            .get_patch_boundary_points()
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
    fn visual_outline<'py>(
        &self,
        py: Python<'py>,
        facing: Vector3,
        max_edge_length: f64,
        corner_angle: Option<f64>,
    ) -> (Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArray1<u8>>) {
        let n = engeom::UnitVec3::new_normalize(*facing.get_inner());
        let outline = self.inner.visual_outline(n, max_edge_length, corner_angle);
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
        (result.into_pyarray(py), result_type.into_pyarray(py))
    }

    #[getter]
    fn face_normals<'py>(&mut self, py: Python<'py>) -> PyResult<&Bound<'py, PyArray2<f64>>> {
        if self.face_normals.is_none() {
            let normals = self
                .inner
                .get_face_normals()
                .map_err(|e| PyValueError::new_err(e.to_string()))?
                .into_iter()
                .map(|n| n.into_inner())
                .collect::<Vec<_>>();

            let array = vectors_to_array(&normals);
            self.face_normals = Some(array.into_pyarray(py).unbind());
        }

        Ok(self.face_normals.as_ref().unwrap().bind(py))
    }

    fn boundary_curves(&self) -> PyResult<Vec<Curve3>> {
        let edges = self
            .inner
            .calc_edges()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let mut curves = Vec::new();
        for item in edges.boundary_loops.iter() {
            let mut points = Vec::new();
            for &index in item.iter() {
                let point = self.inner.vertices()[index as usize];
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
            "<Mesh {} vertices, {} faces>",
            self.inner.vertices().len(),
            self.inner.faces().len()
        )
    }

    fn split(&self, plane: &Plane3) -> PyResult<(Option<Self>, Option<Self>)> {
        match self.inner.split(&plane.inner) {
            SplitResult::Pair(mesh1, mesh2) => {
                Ok((Some(Self::from_inner(mesh1)), Some(Self::from_inner(mesh2))))
            }
            SplitResult::Negative => Ok((Some(self.clone()), None)),
            SplitResult::Positive => Ok((None, Some(self.clone()))),
        }
    }

    fn deviation<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        mode: DeviationMode,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let points = array_to_points3(&points.as_array())?;
        let mut result = Array1::zeros(points.len());

        for (i, point) in points.iter().enumerate() {
            let closest = self.inner.surf_closest_to(point);
            let normal_dev = closest.sp.scalar_projection(point);

            result[i] = match mode {
                // Copy the sign of the normal deviation
                DeviationMode::Point => dist(&closest.sp.point, point) * normal_dev.signum(),
                DeviationMode::Plane => normal_dev,
            }
        }

        Ok(result.into_pyarray(py))
    }

    fn measure_point_deviation(
        &self,
        x: f64,
        y: f64,
        z: f64,
        dist_mode: DeviationMode,
    ) -> Distance3 {
        let point = engeom::Point3::new(x, y, z);
        Distance3::from_inner(self.inner.measure_point_deviation(&point, dist_mode.into()))
    }

    fn boundary_first_flatten<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let edges = self
            .inner
            .calc_edges()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let values = edges
            .boundary_first_flatten()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let result = points_to_array(&values);
        Ok(result.into_pyarray(py))
    }

    fn sample_poisson<'py>(&self, py: Python<'py>, radius: f64) -> Bound<'py, PyArray2<f64>> {
        let mps = self.inner.sample_poisson(radius);
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
        reference: &Mesh,
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

    #[pyo3(signature=(plane, tol = None))]
    fn section(&self, plane: Plane3, tol: Option<f64>) -> PyResult<Vec<Curve3>> {
        let results = self
            .inner
            .section(plane.get_inner(), tol)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(results.into_iter().map(Curve3::from_inner).collect())
    }

    fn face_select_all<'py>(
        slf: PyRef<Self>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, FaceFilterHandle>> {
        let indices = slf.inner.face_select(Selection::All).collect_indices();
        FaceFilterHandle {
            mesh: slf.into(),
            indices,
        }
        .into_pyobject(py)
    }

    fn face_select_none<'py>(
        slf: PyRef<Self>,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, FaceFilterHandle>> {
        let indices = slf.inner.face_select(Selection::None).collect_indices();
        FaceFilterHandle {
            mesh: slf.into(),
            indices,
        }
        .into_pyobject(py)
    }

    fn create_from_indices(&self, indices: Vec<usize>) -> Self {
        Self::from_inner(self.inner.create_from_indices(&indices))
    }

    fn separate_patches(&self) -> PyResult<Vec<Self>> {
        let patch_groups = self.inner.get_patches(None)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut results = Vec::with_capacity(patch_groups.len());
        for mask in patch_groups.iter() {
            results.push(
                self.inner
                    .create_from_mask(mask)
                    .map_err(|e| PyValueError::new_err(e.to_string()))?,
            );
        }
        Ok(results.into_iter().map(Self::from_inner).collect())
    }

    fn convex_hull(&self) -> Self {
        Self::from_inner(self.inner.convex_hull())
    }

    #[staticmethod]
    fn create_box(length: f64, width: f64, height: f64) -> Self {
        let mesh = engeom::Mesh::create_box(length, width, height, true);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        let mesh = engeom::Mesh::create_cylinder(radius, height, steps);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let mesh = engeom::Mesh::create_sphere(radius, n_theta, n_phi);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_cone(radius: f64, height: f64, steps: usize) -> Self {
        let mesh = engeom::Mesh::create_cone(radius, height, steps);
        Self::from_inner(mesh)
    }

    #[staticmethod]
    fn create_capsule(p0: Point3, p1: Point3, radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let mesh =
            engeom::Mesh::create_capsule(p0.get_inner(), p1.get_inner(), radius, n_theta, n_phi);
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
        let mesh = engeom::Mesh::create_rect_beam_between(
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
            engeom::Mesh::create_cylinder_between(p0.get_inner(), p1.get_inner(), radius, steps);
        Ok(Self::from_inner(mesh))
    }

    #[staticmethod]
    fn load_lptf3(file_path: PathBuf, params: Lptf3Load) -> PyResult<Mesh> {
        let half_edge_mesh = engeom::io::load_lptf3_mesh(&file_path, params.into())
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        let mesh = engeom::Mesh::try_from(&half_edge_mesh)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok(Self::from_inner(mesh))
    }
}

#[pyclass]
pub struct FaceFilterHandle {
    mesh: Py<Mesh>,
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
        mode: SelectOp,
    ) -> PyResult<Bound<'py, Self>> {
        let normal = engeom::UnitVec3::new_normalize([x, y, z].into());
        let temp = slf.mesh.bind(py).borrow();
        let i = slf.indices.clone();
        slf.indices = temp
            .inner
            .face_select(Selection::Indices(i))
            .facing(&normal, angle, mode.into())
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(other, all_points, distance_tol, mode, planar_tol = None, angle_tol = None))]
    fn near_mesh<'py>(
        mut slf: PyRefMut<'py, Self>,
        py: Python<'py>,
        other: PyRef<Mesh>,
        all_points: bool,
        distance_tol: f64,
        mode: SelectOp,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
    ) -> PyResult<Bound<'py, Self>> {
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
                mode.into(),
            )
            .collect_indices();
        slf.into_pyobject(py)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn collect(&self) -> Vec<usize> {
        self.indices.clone()
    }

    fn create_mesh(&self, py: Python<'_>) -> Mesh {
        self.mesh
            .bind(py)
            .borrow()
            .create_from_indices(self.indices.clone())
    }
}

#[pyclass]
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

    fn add_stationary(&mut self, mesh: &Mesh) -> usize {
        let inner = mesh.inner.clone();
        self.inner.add_stationary(inner)
    }

    fn add_moving(&mut self, mesh: &Mesh) -> usize {
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
