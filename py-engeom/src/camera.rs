//! Bindings for the camera and optics model in `engeom::sensors::camera`.
//!
//! This layer handles conversion and delegation. The Rust API implements the optics, ray casting,
//! and projection. This layer determines which accessors are properties, how grid-shaped and
//! list-shaped data cross the boundary as NumPy arrays, and where a Rust `Result` becomes a
//! `ValueError`.

use crate::bounding::Aabb2;
use crate::common::IndexMask;
use crate::conversions::{array_to_points3, matrix_to_array, points_to_array};
use crate::geom2::Point2;
use crate::geom3::{Iso3, Point3, Vector3};
use crate::mesh::Mesh3;
use crate::point_cloud::PointCloud3;
use crate::ray_casting::RayBundle3;
use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::{Bound, Py, PyErr, PyResult, Python, pyclass, pyfunction, pymethods};

/// The conventional representative wavelength for visible green light, in millimeters.
pub const WAVELENGTH_550NM_MM: f64 = engeom::sensors::camera::WAVELENGTH_550NM_MM;

/// A projected outline: an `(n, 4)` array of `x0, y0, x1, y1` in pixels, and the per-segment
/// visibility codes.
type ProjectedOutline<'py> = (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray1<u8>>);

fn err(e: Box<dyn std::error::Error>) -> PyErr {
    PyValueError::new_err(e.to_string())
}

// ================================================================================================
// Sensor
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.sensors")]
#[derive(Clone, Copy)]
pub struct Sensor {
    inner: engeom::sensors::camera::Sensor,
}

impl Sensor {
    pub fn get_inner(&self) -> &engeom::sensors::camera::Sensor {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::camera::Sensor) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Sensor {
    #[new]
    fn new(width_px: u32, height_px: u32, pitch: f64) -> PyResult<Self> {
        Ok(Self {
            inner: engeom::sensors::camera::Sensor::new(width_px, height_px, pitch).map_err(err)?,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Sensor(width_px={}, height_px={}, pitch={})",
            self.inner.width_px, self.inner.height_px, self.inner.pitch
        )
    }

    #[getter]
    fn width_px(&self) -> u32 {
        self.inner.width_px
    }

    #[getter]
    fn height_px(&self) -> u32 {
        self.inner.height_px
    }

    #[getter]
    fn pitch(&self) -> f64 {
        self.inner.pitch
    }

    #[getter]
    fn width(&self) -> f64 {
        self.inner.width()
    }

    #[getter]
    fn height(&self) -> f64 {
        self.inner.height()
    }

    #[getter]
    fn diagonal(&self) -> f64 {
        self.inner.diagonal()
    }

    #[getter]
    fn aspect_ratio(&self) -> f64 {
        self.inner.aspect_ratio()
    }

    #[getter]
    fn pixel_count(&self) -> u64 {
        self.inner.pixel_count()
    }

    #[getter]
    fn center_px(&self) -> Point2 {
        Point2::from_inner(self.inner.center_px())
    }

    #[getter]
    fn image_aabb(&self) -> Aabb2 {
        Aabb2::from_inner(self.inner.image_aabb())
    }

    fn rescaled(&self, factor: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.rescaled(factor).map_err(err)?,
        })
    }

    fn contains_pixel(&self, pixel: Point2) -> bool {
        self.inner.contains_pixel(pixel.get_inner())
    }

    fn pixel_to_plane(&self, pixel: Point2) -> Point2 {
        Point2::from_inner(self.inner.pixel_to_plane(pixel.get_inner()))
    }

    fn plane_to_pixel(&self, plane: Point2) -> Point2 {
        Point2::from_inner(self.inner.plane_to_pixel(plane.get_inner()))
    }
}

// ================================================================================================
// ThinLens
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.sensors")]
#[derive(Clone, Copy)]
pub struct ThinLens {
    inner: engeom::sensors::camera::ThinLens,
}

impl ThinLens {
    pub fn get_inner(&self) -> &engeom::sensors::camera::ThinLens {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::camera::ThinLens) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl ThinLens {
    #[new]
    fn new(focal_length: f64, f_number: f64, focus_distance: f64) -> PyResult<Self> {
        Ok(Self {
            inner: engeom::sensors::camera::ThinLens::new(focal_length, f_number, focus_distance)
                .map_err(err)?,
        })
    }

    #[staticmethod]
    fn from_magnification(focal_length: f64, f_number: f64, magnification: f64) -> PyResult<Self> {
        Ok(Self {
            inner: engeom::sensors::camera::ThinLens::from_magnification(
                focal_length,
                f_number,
                magnification,
            )
            .map_err(err)?,
        })
    }

    #[staticmethod]
    fn from_field_width(
        f_number: f64,
        focus_distance: f64,
        sensor_width: f64,
        field_width: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: engeom::sensors::camera::ThinLens::from_field_width(
                f_number,
                focus_distance,
                sensor_width,
                field_width,
            )
            .map_err(err)?,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "ThinLens(focal_length={}, f_number={}, focus_distance={})",
            self.inner.focal_length, self.inner.f_number, self.inner.focus_distance
        )
    }

    #[getter]
    fn focal_length(&self) -> f64 {
        self.inner.focal_length
    }

    #[getter]
    fn f_number(&self) -> f64 {
        self.inner.f_number
    }

    #[getter]
    fn focus_distance(&self) -> f64 {
        self.inner.focus_distance
    }

    #[getter]
    fn image_distance(&self) -> f64 {
        self.inner.image_distance()
    }

    #[getter]
    fn magnification(&self) -> f64 {
        self.inner.magnification()
    }

    #[getter]
    fn effective_f_number(&self) -> f64 {
        self.inner.effective_f_number()
    }

    #[getter]
    fn aperture_diameter(&self) -> f64 {
        self.inner.aperture_diameter()
    }

    #[getter]
    fn is_focused_at_infinity(&self) -> bool {
        self.inner.is_focused_at_infinity()
    }

    fn focused_at(&self, focus_distance: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.focused_at(focus_distance).map_err(err)?,
        })
    }

    fn hyperfocal_distance(&self, coc: f64) -> f64 {
        self.inner.hyperfocal_distance(coc)
    }

    fn airy_disk_diameter(&self, wavelength: f64) -> f64 {
        self.inner.airy_disk_diameter(wavelength)
    }

    fn coc_diameter(&self, z: f64) -> Option<f64> {
        self.inner.coc_diameter(z)
    }

    fn depth_of_field(&self, coc: f64) -> Option<(f64, f64)> {
        self.inner.depth_of_field(coc)
    }

    fn total_depth_of_field(&self, coc: f64) -> Option<f64> {
        self.inner.total_depth_of_field(coc)
    }
}

// ================================================================================================
// PinholeCamera
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.sensors")]
#[derive(Clone, Copy)]
pub struct PinholeCamera {
    inner: engeom::sensors::camera::PinholeCamera,
}

impl PinholeCamera {
    pub fn get_inner(&self) -> &engeom::sensors::camera::PinholeCamera {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::camera::PinholeCamera) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PinholeCamera {
    #[new]
    fn new(fx: f64, fy: f64, cx: f64, cy: f64) -> Self {
        Self {
            inner: engeom::sensors::camera::PinholeCamera::new(fx, fy, cx, cy),
        }
    }

    #[staticmethod]
    fn from_focal_length(focal_length: f64, width: u32, height: u32) -> Self {
        Self {
            inner: engeom::sensors::camera::PinholeCamera::from_focal_length(
                focal_length,
                width,
                height,
            ),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PinholeCamera(fx={}, fy={}, cx={}, cy={})",
            self.inner.fx, self.inner.fy, self.inner.cx, self.inner.cy
        )
    }

    #[getter]
    fn fx(&self) -> f64 {
        self.inner.fx
    }

    #[getter]
    fn fy(&self) -> f64 {
        self.inner.fy
    }

    #[getter]
    fn cx(&self) -> f64 {
        self.inner.cx
    }

    #[getter]
    fn cy(&self) -> f64 {
        self.inner.cy
    }
}

// ================================================================================================
// Camera
// ================================================================================================

#[pyclass(from_py_object, module = "engeom.sensors")]
#[derive(Clone, Copy)]
pub struct Camera {
    inner: engeom::sensors::camera::Camera,
}

impl Camera {
    pub fn get_inner(&self) -> &engeom::sensors::camera::Camera {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::camera::Camera) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl Camera {
    #[new]
    fn new(sensor: Sensor, lens: ThinLens) -> Self {
        Self {
            inner: engeom::sensors::camera::Camera::new(*sensor.get_inner(), *lens.get_inner()),
        }
    }

    #[staticmethod]
    fn from_field_width(
        sensor: Sensor,
        f_number: f64,
        working_distance: f64,
        field_width: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: engeom::sensors::camera::Camera::from_field_width(
                *sensor.get_inner(),
                f_number,
                working_distance,
                field_width,
            )
            .map_err(err)?,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "Camera(sensor={}, lens={})",
            Sensor::from_inner(self.inner.sensor).__repr__(),
            ThinLens::from_inner(self.inner.lens).__repr__()
        )
    }

    #[getter]
    fn sensor(&self) -> Sensor {
        Sensor::from_inner(self.inner.sensor)
    }

    #[getter]
    fn lens(&self) -> ThinLens {
        ThinLens::from_inner(self.inner.lens)
    }

    #[getter]
    fn pinhole(&self) -> PinholeCamera {
        PinholeCamera::from_inner(self.inner.pinhole())
    }

    #[getter]
    fn horizontal_fov(&self) -> f64 {
        self.inner.horizontal_fov()
    }

    #[getter]
    fn vertical_fov(&self) -> f64 {
        self.inner.vertical_fov()
    }

    #[getter]
    fn diagonal_fov(&self) -> f64 {
        self.inner.diagonal_fov()
    }

    #[getter]
    fn footprint(&self) -> (f64, f64) {
        self.inner.footprint()
    }

    #[getter]
    fn pixel_size(&self) -> f64 {
        self.inner.pixel_size()
    }

    fn rescaled(&self, factor: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.rescaled(factor).map_err(err)?,
        })
    }

    fn footprint_at(&self, z: f64) -> Option<(f64, f64)> {
        self.inner.footprint_at(z)
    }

    fn pixel_size_at(&self, z: f64) -> Option<f64> {
        self.inner.pixel_size_at(z)
    }

    fn coc_px_at(&self, z: f64) -> Option<f64> {
        self.inner.coc_px_at(z)
    }

    fn depth_of_field_px(&self, coc_px: f64) -> Option<(f64, f64)> {
        self.inner.depth_of_field_px(coc_px)
    }

    fn total_depth_of_field_px(&self, coc_px: f64) -> Option<f64> {
        self.inner.total_depth_of_field_px(coc_px)
    }

    fn hyperfocal_distance_px(&self, coc_px: f64) -> f64 {
        self.inner.hyperfocal_distance_px(coc_px)
    }

    fn airy_disk_px(&self, wavelength: f64) -> f64 {
        self.inner.airy_disk_px(wavelength)
    }

    /// Project world points into the image. Points behind the camera produce rows of NaN, so each
    /// result row corresponds to the same input row.
    fn project<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'_, f64>,
        iso: &Iso3,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let pts = array_to_points3(&points.as_array())?;
        let mut result = Array2::from_elem((pts.len(), 2), f64::NAN);
        for (i, p) in pts.iter().enumerate() {
            if let Some(px) = self.inner.project(p, iso.get_inner()) {
                result[[i, 0]] = px.x;
                result[[i, 1]] = px.y;
            }
        }
        Ok(result.into_pyarray(py))
    }

    /// Back-project image points into world space rays, as a bundle which can be cast at a mesh.
    fn back_project(&self, pixels: PyReadonlyArray2<'_, f64>, iso: &Iso3) -> PyResult<RayBundle3> {
        let view = pixels.as_array();
        if view.shape()[1] != 2 {
            return Err(PyValueError::new_err("pixels must be an (n, 2) array"));
        }
        let pts: Vec<engeom::Point2> = view
            .rows()
            .into_iter()
            .map(|r| engeom::Point2::new(r[0], r[1]))
            .collect();
        Ok(RayBundle3::from_inner(
            self.inner.back_project_many(&pts, iso.get_inner()),
        ))
    }

    /// Project a world-space polyline into the image, split into runs wherever it passes behind
    /// the near plane. Each run is its own `(m, 2)` array.
    #[pyo3(signature = (points, iso, near = None))]
    fn project_polyline<'py>(
        &self,
        py: Python<'py>,
        points: PyReadonlyArray2<'_, f64>,
        iso: &Iso3,
        near: Option<f64>,
    ) -> PyResult<Vec<Bound<'py, PyArray2<f64>>>> {
        let pts = array_to_points3(&points.as_array())?;
        Ok(self
            .inner
            .project_polyline(&pts, iso.get_inner(), near)
            .iter()
            .map(|run| points_to_array(run).into_pyarray(py))
            .collect())
    }

    /// A line drawing of a mesh in this camera's pixel space, as an `(n, 4)` array of
    /// `x0, y0, x1, y1` and the per-segment visibility codes.
    #[pyo3(signature = (mesh, iso, max_edge_px = None, corner_angle = None))]
    fn project_outline<'py>(
        &self,
        py: Python<'py>,
        mesh: &Mesh3,
        iso: &Iso3,
        max_edge_px: Option<f64>,
        corner_angle: Option<f64>,
    ) -> PyResult<ProjectedOutline<'py>> {
        let outline = self
            .inner
            .project_outline(mesh.get_inner(), iso.get_inner(), max_edge_px, corner_angle)
            .map_err(err)?;

        let mut segments = Array2::zeros((outline.len(), 4));
        let mut kinds = Array1::zeros(outline.len());
        for (i, (a, b, k)) in outline.iter().enumerate() {
            segments[[i, 0]] = a.x;
            segments[[i, 1]] = a.y;
            segments[[i, 2]] = b.x;
            segments[[i, 3]] = b.y;
            kinds[i] = *k;
        }

        Ok((segments.into_pyarray(py), kinds.into_pyarray(py)))
    }

    /// The eight corners of the view volume between two depths, as an `(8, 3)` array. The near
    /// face comes first, each face running from the top left of the image clockwise.
    fn frustum_corners<'py>(
        &self,
        py: Python<'py>,
        iso: &Iso3,
        near: f64,
        far: f64,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let corners = self
            .inner
            .frustum_corners(iso.get_inner(), near, far)
            .map_err(err)?;
        Ok(points_to_array(&corners).into_pyarray(py))
    }
}

// ================================================================================================
// ViewBuffer
// ================================================================================================

/// A per-pixel view of a scene rendered through a camera.
///
/// The per-pixel maps are returned as cached `(height, width)` NumPy arrays. Because this API
/// cannot mutate the buffer, these arrays do not require the cache invalidation used by
/// `PointCloud3`.
///
/// `skip_from_py_object` prevents the buffer from being accepted by value. At full sensor
/// resolution, its pixel array can use hundreds of megabytes, and a by-value argument would clone
/// that array on every call.
#[pyclass(skip_from_py_object, module = "engeom.sensors")]
pub struct ViewBuffer {
    inner: engeom::sensors::camera::ViewBuffer,
    depth: Option<Py<PyArray2<f64>>>,
    coc_px: Option<Py<PyArray2<f64>>>,
    incidence: Option<Py<PyArray2<f64>>>,
}

impl ViewBuffer {
    pub fn get_inner(&self) -> &engeom::sensors::camera::ViewBuffer {
        &self.inner
    }

    pub fn from_inner(inner: engeom::sensors::camera::ViewBuffer) -> Self {
        Self {
            inner,
            depth: None,
            coc_px: None,
            incidence: None,
        }
    }
}

#[pymethods]
impl ViewBuffer {
    /// Render a view by casting a ray through the center of every pixel of the camera.
    #[new]
    #[pyo3(signature = (camera, iso, target, obstruction = None))]
    fn new(camera: Camera, iso: &Iso3, target: &Mesh3, obstruction: Option<&Mesh3>) -> Self {
        Self::from_inner(engeom::sensors::camera::ViewBuffer::render(
            camera.get_inner(),
            iso.get_inner(),
            target.get_inner(),
            obstruction.map(|o| o.get_inner()),
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "ViewBuffer(width={}, height={}, hit_fraction={:.3})",
            self.inner.width(),
            self.inner.height(),
            self.inner.hit_fraction()
        )
    }

    #[getter]
    fn camera(&self) -> Camera {
        Camera::from_inner(*self.inner.camera())
    }

    #[getter]
    fn iso(&self) -> Iso3 {
        Iso3::from_inner(*self.inner.iso())
    }

    #[getter]
    fn width(&self) -> u32 {
        self.inner.width()
    }

    #[getter]
    fn height(&self) -> u32 {
        self.inner.height()
    }

    #[getter]
    fn target_face_count(&self) -> usize {
        self.inner.target_face_count()
    }

    #[getter]
    fn hit_count(&self) -> usize {
        self.inner.hit_count()
    }

    #[getter]
    fn hit_fraction(&self) -> f64 {
        self.inner.hit_fraction()
    }

    #[getter]
    fn depth<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.depth.is_none() {
            let array = matrix_to_array(&self.inner.to_depth_matrix());
            self.depth = Some(array.into_pyarray(py).unbind());
        }
        self.depth.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn coc_px<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.coc_px.is_none() {
            let array = matrix_to_array(&self.inner.to_coc_px_matrix());
            self.coc_px = Some(array.into_pyarray(py).unbind());
        }
        self.coc_px.as_ref().unwrap().bind(py)
    }

    #[getter]
    fn incidence<'py>(&mut self, py: Python<'py>) -> &Bound<'py, PyArray2<f64>> {
        if self.incidence.is_none() {
            let array = matrix_to_array(&self.inner.to_incidence_matrix());
            self.incidence = Some(array.into_pyarray(py).unbind());
        }
        self.incidence.as_ref().unwrap().bind(py)
    }

    /// The depth at a single pixel, NaN where the pixel saw nothing or lies outside the buffer.
    fn depth_at(&self, x: i32, y: i32) -> f64 {
        self.inner.depth_at(engeom::raster2::Point2I::new(x, y))
    }

    fn to_point_cloud(&self) -> PointCloud3 {
        PointCloud3::from_inner(self.inner.to_point_cloud())
    }

    #[getter]
    fn face_mask(&self) -> IndexMask {
        IndexMask::from_inner(self.inner.face_mask())
    }

    /// A mask over the target faces which at least one pixel saw well enough to satisfy every
    /// supplied criterion. A criterion left as `None` is not applied.
    #[pyo3(signature = (max_coc_px = None, max_incidence = None, allow_backface = false))]
    fn face_mask_usable(
        &self,
        max_coc_px: Option<f64>,
        max_incidence: Option<f64>,
        allow_backface: bool,
    ) -> IndexMask {
        IndexMask::from_inner(self.inner.face_mask_usable(
            max_coc_px,
            max_incidence,
            allow_backface,
        ))
    }
}

// ================================================================================================
// Free functions
// ================================================================================================

/// The camera-to-world isometry for a camera at `eye` whose optical axis passes through `target`,
/// oriented so `up` points as nearly as possible to the top of the image.
#[pyfunction]
pub fn look_at(eye: Point3, target: Point3, up: Vector3) -> PyResult<Iso3> {
    Ok(Iso3::from_inner(
        engeom::sensors::camera::look_at(eye.get_inner(), target.get_inner(), up.get_inner())
            .map_err(err)?,
    ))
}
