//! This module contains models of cameras: the optics that form an image and the discrete
//! sensor that records it.
//!
//! # Units and conventions
//!
//! Lengths use the caller's chosen unit, usually millimeters, and every derived length uses the
//! same unit. Wavelengths must also use that unit. For millimeters, use
//! [`WAVELENGTH_550NM_MM`]. Angles are in radians.
//!
//! The camera coordinate system is the standard computer vision convention: +X points right in
//! the image, +Y points down, and +Z points out along the optical axis into the scene. A camera
//! does not store its own pose. Every geometric operation takes a camera-to-world [`Iso3`], so
//! a single [`Camera`] can be used at the many poses of a robot or fixture. The
//! [`look_at`] function builds such a pose from an eye point, a target, and an up direction.
//!
//! The plotting-view convention documented in `CONVENTIONS.md` instead has +Y pointing up and
//! +Z pointing out of the image toward the viewer. The two conventions differ by a half turn
//! about X. In addition, a view isometry is world-to-view, while a camera pose is
//! camera-to-world. Convert a view with `iso = view.inverse() * Iso3::from_rx(PI)`.
//!
//! # The image distance, not the focal length
//!
//! The pixel-unit intrinsics of a [`Camera`] are built from the lens's image distance rather
//! than its focal length. The two are the same only when the lens is focused at infinity; a lens
//! focused close projects as though it were longer than its nameplate focal length, which is a
//! small effect at ordinary working distances and a large one in macro imaging.

mod pinhole;
mod sensor;
mod thin_lens;
mod view_buffer;

pub use pinhole::PinholeCamera;
pub use sensor::Sensor;
pub use thin_lens::{ThinLens, WAVELENGTH_550NM_MM};
pub use view_buffer::{PixelHit, ViewBuffer};

use crate::common::PCoords;
use crate::geom3::IsoExtensions3;
use crate::{Iso3, Mesh3, Point2, Point3, Result, Vector3};
use parry3d_f64::query::Ray;

/// A physical camera: a [`Sensor`] positioned behind a [`ThinLens`] at the image distance that
/// brings the lens's focus distance into focus.
///
/// A `Camera` provides the quantities used to design an imaging system. Without a scene, it can
/// determine how much of the world fits in the frame at a given distance, how large a feature is
/// in pixels, how far the depth of field extends, and whether diffraction or pixel sampling
/// limits the system. With a camera-to-world isometry, it can also determine where a world point
/// lands in the image and which world-space ray passes through a pixel.
///
/// The pose is not stored. Every operation that needs one takes it as an argument, so a single
/// camera can be evaluated at many poses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// The discrete sensor which records the image.
    pub sensor: Sensor,

    /// The lens which forms the image on the sensor.
    pub lens: ThinLens,
}

impl Camera {
    /// Create a new camera from a sensor and a lens. Both are already validated by their own
    /// constructors, so this cannot fail.
    ///
    /// # Arguments
    ///
    /// * `sensor`: the discrete sensor which records the image
    /// * `lens`: the lens which forms the image on the sensor
    ///
    /// returns: Camera
    pub fn new(sensor: Sensor, lens: ThinLens) -> Self {
        Self { sensor, lens }
    }

    /// Create a camera with a lens that fills the sensor width with a specified width of the
    /// object plane at a specified working distance. This is the usual way to select a lens when
    /// the part dimensions and standoff are known and the focal length must be calculated.
    ///
    /// # Arguments
    ///
    /// * `sensor`: the discrete sensor which records the image
    /// * `f_number`: the f-number of the lens at infinity focus
    /// * `working_distance`: the object-side distance to the plane of best focus
    /// * `field_width`: the width of the object plane to be imaged across the sensor width
    ///
    /// returns: Result<Camera>
    pub fn from_field_width(
        sensor: Sensor,
        f_number: f64,
        working_distance: f64,
        field_width: f64,
    ) -> Result<Self> {
        let lens =
            ThinLens::from_field_width(f_number, working_distance, sensor.width(), field_width)?;
        Ok(Self::new(sensor, lens))
    }

    /// Create a copy of this camera with its sensor resolution scaled by a factor. Use a reduced
    /// resolution to render a less expensive preview. Because rescaling preserves the sensor's
    /// physical size, the copy has the same field of view and optics as the original. Only the
    /// sampling resolution changes.
    ///
    /// Take care with quantities measured in pixels. Rescaling changes the pixel size, so pixel
    /// measurements of blur or depth of field apply to the sensor on which they were measured.
    /// On a copy scaled by a factor, `coc_px_at` returns the original value multiplied by that
    /// factor. A copy scaled down therefore reports less blur for the same physical defocus. A
    /// focus criterion transferred unchanged from a preview to a full-resolution render will be
    /// much more permissive than intended. Scale the threshold by the same factor, or define the
    /// criterion with `ThinLens::coc_diameter`, which returns a length independent of sampling.
    ///
    /// # Arguments
    ///
    /// * `factor`: the factor to scale the pixel dimensions by, which must be a finite number
    ///   greater than zero
    ///
    /// returns: Result<Camera>
    pub fn rescaled(&self, factor: f64) -> Result<Self> {
        Ok(Self::new(self.sensor.rescaled(factor)?, self.lens))
    }

    /// The pixel-unit intrinsics that project points into the image. The focal lengths come from
    /// the lens's image distance divided by the pixel pitch,
    /// and the principal point is at the center of the sensor.
    pub fn pinhole(&self) -> PinholeCamera {
        let f_px = self.lens.image_distance() / self.sensor.pitch;
        let c = self.sensor.center_px();
        PinholeCamera::new(f_px, f_px, c.x, c.y)
    }

    /// The ratio between a length on the object plane at a distance `z` and the length it
    /// occupies on the sensor. Every object-space quantity is scaled by this.
    ///
    /// Note: this is the only place where perspective projection enters the object-space
    /// calculations. A future telecentric projection would change this function and `pinhole`,
    /// and nothing else.
    fn object_scale(&self, z: f64) -> f64 {
        z / self.lens.image_distance()
    }

    /// Project a single world-space point into the image.
    ///
    /// Returns `None` if the point is behind the camera. A point can project to a
    /// pixel coordinate outside the bounds of the sensor, which means it is within the optical
    /// field but not recorded.
    ///
    /// # Arguments
    ///
    /// * `world_point`: the point to project, expressed in world coordinates
    /// * `iso`: the camera-to-world transform
    ///
    /// returns: Option<Point2>
    pub fn project(&self, world_point: &impl PCoords<3>, iso: &Iso3) -> Option<Point2> {
        self.pinhole().project(world_point, iso)
    }

    /// Project a slice of world-space points into the image, with points behind the camera
    /// producing `None`.
    ///
    /// # Arguments
    ///
    /// * `world_points`: the points to project, expressed in world coordinates
    /// * `iso`: the camera-to-world transform
    ///
    /// returns: Vec<Option<Point2>>
    pub fn project_many(
        &self,
        world_points: &[impl PCoords<3>],
        iso: &Iso3,
    ) -> Vec<Option<Point2>> {
        self.pinhole().project_many(world_points, iso)
    }

    /// Back-project a point in the image into a ray in world space, originating at the camera
    /// center and pointing into the scene.
    ///
    /// # Arguments
    ///
    /// * `image_point`: pixel coordinates in the image
    /// * `iso`: the camera-to-world transform
    ///
    /// returns: Ray
    pub fn back_project(&self, image_point: &Point2, iso: &Iso3) -> Ray {
        self.pinhole().back_project(image_point, iso)
    }

    /// Back-project a slice of points in the image into rays in world space.
    ///
    /// # Arguments
    ///
    /// * `image_points`: pixel coordinates in the image
    /// * `iso`: the camera-to-world transform
    ///
    /// returns: Vec<Ray>
    pub fn back_project_many(&self, image_points: &[Point2], iso: &Iso3) -> Vec<Ray> {
        self.pinhole().back_project_many(image_points, iso)
    }

    /// Resolve the near-plane distance, using the lens focal length when the supplied value is
    /// absent, nonfinite, or at or below zero. Such a value cannot describe a plane in front of
    /// the camera.
    fn near_or_default(&self, near: Option<f64>) -> f64 {
        near.filter(|n| n.is_finite() && *n > 0.0)
            .unwrap_or(self.lens.focal_length)
    }

    /// Project a point already expressed in camera space, without re-testing which side of the
    /// camera it is on. The caller is responsible for having clipped it first.
    fn project_cam(pinhole: &PinholeCamera, p: &Point3) -> Point2 {
        Point2::new(
            pinhole.fx * p.x / p.z + pinhole.cx,
            pinhole.fy * p.y / p.z + pinhole.cy,
        )
    }

    /// Project a segment whose camera-space transform and intrinsics have already been built, so
    /// that a loop over many segments does not rebuild them per call.
    fn project_segment_with(
        pinhole: &PinholeCamera,
        to_cam: &Iso3,
        a: &impl PCoords<3>,
        b: &impl PCoords<3>,
        near: f64,
    ) -> Option<(Point2, Point2)> {
        let a = to_cam * Point3::from(a.coords());
        let b = to_cam * Point3::from(b.coords());
        let (a, b) = clip_near(&a, &b, near)?;
        Some((
            Self::project_cam(pinhole, &a),
            Self::project_cam(pinhole, &b),
        ))
    }

    /// Project a world-space line segment into the image, clipping away any part that
    /// lies behind the near plane.
    ///
    /// Both endpoints of the returned pair lie in front of the near plane, and they keep the
    /// order of the arguments, so the first corresponds to the `a` end of the original segment.
    /// When the segment straddles the near plane, the endpoint behind it is replaced by the point
    /// where the segment crosses. A projected point can lie outside the bounds of the sensor,
    /// which means it is within the optical field but not recorded; use `Sensor::contains_pixel`
    /// to test that.
    ///
    /// A segment with a NaN coordinate is treated as lying behind the camera.
    ///
    /// # Arguments
    ///
    /// * `a`: the first endpoint, expressed in world coordinates
    /// * `b`: the second endpoint, expressed in world coordinates
    /// * `iso`: the camera-to-world transform
    /// * `near`: the camera space distance to the near plane. If `None`, or not a finite number
    ///   greater than zero, the lens focal length is used, because an object nearer than the
    ///   focal length forms no real image.
    ///
    /// returns: Option<(Point2, Point2)>, or `None` if no part of the segment lies in front of
    /// the near plane
    pub fn project_segment(
        &self,
        a: &impl PCoords<3>,
        b: &impl PCoords<3>,
        iso: &Iso3,
        near: Option<f64>,
    ) -> Option<(Point2, Point2)> {
        Self::project_segment_with(
            &self.pinhole(),
            &iso.inverse(),
            a,
            b,
            self.near_or_default(near),
        )
    }

    /// Project a world space polyline into the image, splitting it wherever it passes behind the
    /// near plane.
    ///
    /// The result is a list of runs of consecutive pixel coordinates. A polyline entirely in
    /// front of the near plane produces one run holding every point; one which dips behind it
    /// produces one run per visible stretch, each terminated at the point where the polyline
    /// crosses the plane. Runs appear in the order they occur along the polyline, and the points
    /// within a run keep the input order. Every run holds at least two points, except when the
    /// input itself was a single visible point.
    ///
    /// The result contains raw points because a polyline seen nearly edge-on can project to fewer
    /// than two distinct pixels, which `Curve2` does not accept. Whether a projected run is closed
    /// is also a property of the input rather than the projection.
    ///
    /// A closed polyline must repeat its first point as its last. If such a polyline is clipped
    /// across that repeated point, its first and last runs are two ends of one visible stretch
    /// and are not joined.
    ///
    /// # Arguments
    ///
    /// * `points`: the polyline vertices in order, expressed in world coordinates
    /// * `iso`: the camera-to-world transform
    /// * `near`: the camera space distance to the near plane. If `None`, or not a finite number
    ///   greater than zero, the lens focal length is used.
    ///
    /// returns: Vec<Vec<Point2>>
    pub fn project_polyline(
        &self,
        points: &[impl PCoords<3>],
        iso: &Iso3,
        near: Option<f64>,
    ) -> Vec<Vec<Point2>> {
        let near = self.near_or_default(near);
        let pinhole = self.pinhole();
        let to_cam = iso.inverse();

        let mut runs = Vec::new();
        let mut run: Vec<Point2> = Vec::new();

        if points.is_empty() {
            return runs;
        }

        let mut prev = to_cam * Point3::from(points[0].coords());
        if points.len() == 1 {
            if prev.z >= near {
                runs.push(vec![Self::project_cam(&pinhole, &prev)]);
            }
            return runs;
        }

        for next in points.iter().skip(1) {
            let curr = to_cam * Point3::from(next.coords());
            let prev_in = prev.z >= near;
            let curr_in = curr.z >= near;

            match (prev_in, curr_in) {
                (true, true) => {
                    if run.is_empty() {
                        run.push(Self::project_cam(&pinhole, &prev));
                    }
                    run.push(Self::project_cam(&pinhole, &curr));
                }
                (true, false) => {
                    if run.is_empty() {
                        run.push(Self::project_cam(&pinhole, &prev));
                    }
                    let crossing = lerp_to_near(&prev, &curr, near);
                    run.push(Self::project_cam(&pinhole, &crossing));
                    runs.push(std::mem::take(&mut run));
                }
                (false, true) => {
                    let crossing = lerp_to_near(&curr, &prev, near);
                    run.push(Self::project_cam(&pinhole, &crossing));
                    run.push(Self::project_cam(&pinhole, &curr));
                }
                (false, false) => {}
            }

            prev = curr;
        }

        if !run.is_empty() {
            runs.push(run);
        }

        runs
    }

    /// Compute a line drawing of a mesh in this camera's pixel space from the given pose.
    ///
    /// This composes `Mesh3::compute_perspective_outline` with the near plane clipping
    /// projection, so the outline is found for this camera's own eye point and any part of it
    /// behind the camera is clipped away rather than projected to nonsense. A segment which
    /// straddles the near plane is truncated at the crossing and keeps its classification.
    ///
    /// A projected point can lie outside the bounds of the sensor, which means it is within the
    /// optical field but not recorded; use `Sensor::contains_pixel` or `Sensor::image_aabb` to
    /// test or clip against the frame.
    ///
    /// The subdivision length is chosen once for the whole mesh, from the pixel size at the depth
    /// of the mesh's bounding box center. A part with a large extent along the optical axis is
    /// therefore subdivided a little too finely at its near end and a little too coarsely at its
    /// far end, in proportion to the ratio of those depths. A caller needing better can drive
    /// `Mesh3::compute_perspective_outline` and `Camera::project_segment` directly.
    ///
    /// The tolerance is in the pixels of the camera it is given to. A preview camera from
    /// `Camera::rescaled` therefore gets segments which are coarser in world terms, in
    /// proportion to its coarser pixels and therefore gets proportionally fewer segments. This
    /// reduces the preview cost while preserving the drawn line's apparent quality. By contrast,
    /// carrying a `coc_px_at` threshold between resolutions silently changes the physical quantity
    /// being tested.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh to outline
    /// * `iso`: the camera-to-world transform
    /// * `max_edge_px`: the greatest length of a returned segment, in pixels. If `None`, a
    ///   default of four pixels is used.
    /// * `corner_angle`: the smallest angle in radians between two adjacent faces for their
    ///   shared edge to be drawn as a corner. If `None`, a default of just under 45 degrees is
    ///   used.
    ///
    /// returns: Result<Vec<(Point2, Point2, u8)>>, pixel space segments each tagged with 0 when
    /// the segment is visible from the camera and 1 when the mesh hides it
    pub fn project_outline(
        &self,
        mesh: &Mesh3,
        iso: &Iso3,
        max_edge_px: Option<f64>,
        corner_angle: Option<f64>,
    ) -> Result<Vec<(Point2, Point2, u8)>> {
        if mesh.faces().is_empty() {
            return Ok(Vec::new());
        }

        let eye = iso.origin();

        // Measured along the optical axis rather than as a straight line distance, because that
        // is what `pixel_size_at` is defined on. A straight line distance would overstate the
        // pixel size for a part off to one side of the frame.
        let depth = (mesh.aabb().center() - eye).dot(&iso.z().into_inner());
        let world_px = self
            .pixel_size_at(depth)
            .ok_or("the center of the mesh is not in front of the camera")?;

        let max_edge_length = max_edge_px.unwrap_or(DEFAULT_OUTLINE_EDGE_PX) * world_px;
        let outline = mesh.compute_perspective_outline(&eye, max_edge_length, corner_angle)?;

        let pinhole = self.pinhole();
        let to_cam = iso.inverse();
        let near = self.near_or_default(None);

        Ok(outline
            .iter()
            .filter_map(|(a, b, kind)| {
                Self::project_segment_with(&pinhole, &to_cam, a, b, near)
                    .map(|(p, q)| (p, q, *kind))
            })
            .collect())
    }

    /// Return the eight world-coordinate corners of the view volume between two distances along
    /// the optical axis.
    ///
    /// The corners are ordered as the near face followed by the far face, each running from the
    /// top left of the image clockwise: top left, top right, bottom right, bottom left. Setting
    /// `near` and `far` to the limits from `depth_of_field_px` gives the volume in which the part
    /// is both framed and in focus.
    ///
    /// # Arguments
    ///
    /// * `iso`: the camera-to-world transform
    /// * `near`: the distance along the optical axis to the near face, which must be greater than
    ///   zero
    /// * `far`: the distance along the optical axis to the far face, which must be greater than
    ///   `near` and finite
    ///
    /// returns: Result<[Point3; 8]>
    pub fn frustum_corners(&self, iso: &Iso3, near: f64, far: f64) -> Result<[Point3; 8]> {
        if !near.is_finite() || near <= 0.0 {
            return Err("frustum near distance must be a finite number greater than zero".into());
        }
        if !far.is_finite() || far <= near {
            return Err("frustum far distance must be finite and greater than the near".into());
        }

        let pinhole = self.pinhole();
        let w = self.sensor.width_px as f64;
        let h = self.sensor.height_px as f64;
        let image_corners = [
            Point2::new(0.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w, h),
            Point2::new(0.0, h),
        ];

        let mut out = [Point3::origin(); 8];
        for (face, z) in [near, far].iter().enumerate() {
            for (i, c) in image_corners.iter().enumerate() {
                // Built at the requested depth along the axis rather than along the ray, so the
                // two faces are flat and perpendicular to the optical axis.
                let x = (c.x - pinhole.cx) / pinhole.fx * z;
                let y = (c.y - pinhole.cy) / pinhole.fy * z;
                out[face * 4 + i] = iso * Point3::new(x, y, *z);
            }
        }

        Ok(out)
    }

    /// The full horizontal angle of view, in radians.
    pub fn horizontal_fov(&self) -> f64 {
        2.0 * (self.sensor.width() / (2.0 * self.lens.image_distance())).atan()
    }

    /// The full vertical angle of view, in radians.
    pub fn vertical_fov(&self) -> f64 {
        2.0 * (self.sensor.height() / (2.0 * self.lens.image_distance())).atan()
    }

    /// The full diagonal angle of view, in radians.
    pub fn diagonal_fov(&self) -> f64 {
        2.0 * (self.sensor.diagonal() / (2.0 * self.lens.image_distance())).atan()
    }

    /// The width and height of the region of a plane at distance `z` which is imaged onto the
    /// full sensor.
    ///
    /// # Arguments
    ///
    /// * `z`: the object-side distance from the lens plane to the plane of interest, which must
    ///   be greater than zero
    ///
    /// returns: Option<(f64, f64)>, or `None` if `z` is not greater than zero
    pub fn footprint_at(&self, z: f64) -> Option<(f64, f64)> {
        if z.is_nan() || z <= 0.0 {
            return None;
        }
        let scale = self.object_scale(z);
        Some((self.sensor.width() * scale, self.sensor.height() * scale))
    }

    /// The width and height of the region of the focus plane which is imaged onto the full
    /// sensor. A camera focused at infinity returns infinities.
    pub fn footprint(&self) -> (f64, f64) {
        // The focus distance is always greater than the focal length and so always positive
        self.footprint_at(self.lens.focus_distance)
            .expect("focus distance is always positive")
    }

    /// The size on a plane at distance `z` of the region imaged onto a single pixel, which is
    /// the sampling resolution of the camera at that distance.
    ///
    /// # Arguments
    ///
    /// * `z`: the object-side distance from the lens plane to the plane of interest, which must
    ///   be greater than zero
    ///
    /// returns: Option<f64>, or `None` if `z` is not greater than zero
    pub fn pixel_size_at(&self, z: f64) -> Option<f64> {
        if z.is_nan() || z <= 0.0 {
            return None;
        }
        Some(self.sensor.pitch * self.object_scale(z))
    }

    /// The size on the focus plane of the region imaged onto a single pixel. A camera focused at
    /// infinity returns infinity.
    pub fn pixel_size(&self) -> f64 {
        self.pixel_size_at(self.lens.focus_distance)
            .expect("focus distance is always positive")
    }

    /// The diameter, in pixels, of the defocus blur of a point at distance `z`.
    ///
    /// # Arguments
    ///
    /// * `z`: the object-side distance from the lens plane to the point, which must be greater
    ///   than zero
    ///
    /// returns: Option<f64>, or `None` if `z` is not greater than zero
    pub fn coc_px_at(&self, z: f64) -> Option<f64> {
        self.lens.coc_diameter(z).map(|coc| coc / self.sensor.pitch)
    }

    /// The near and far limits of the depth of field, for a blur tolerance expressed in pixels.
    /// Stating the tolerance in pixels rather than as a length ties the depth of field to the
    /// sampling of the sensor, which is usually what matters for a machine vision system.
    ///
    /// # Arguments
    ///
    /// * `coc_px`: the largest acceptable defocus blur diameter in pixels, which must be zero or
    ///   greater
    ///
    /// returns: Option<(f64, f64)>, or `None` if `coc_px` is negative or NaN
    pub fn depth_of_field_px(&self, coc_px: f64) -> Option<(f64, f64)> {
        self.lens.depth_of_field(coc_px * self.sensor.pitch)
    }

    /// The distance between the near and far limits of the depth of field, for a blur tolerance
    /// expressed in pixels.
    ///
    /// # Arguments
    ///
    /// * `coc_px`: the largest acceptable defocus blur diameter in pixels, which must be zero or
    ///   greater
    ///
    /// returns: Option<f64>, or `None` if `coc_px` is negative or NaN
    pub fn total_depth_of_field_px(&self, coc_px: f64) -> Option<f64> {
        self.lens.total_depth_of_field(coc_px * self.sensor.pitch)
    }

    /// The hyperfocal distance for a blur tolerance expressed in pixels, which is the focus
    /// distance beyond which the far limit of the depth of field reaches infinity.
    ///
    /// # Arguments
    ///
    /// * `coc_px`: the largest acceptable defocus blur diameter in pixels
    ///
    /// returns: f64
    pub fn hyperfocal_distance_px(&self, coc_px: f64) -> f64 {
        self.lens.hyperfocal_distance(coc_px * self.sensor.pitch)
    }

    /// The diameter of the Airy disk in pixels. This is the diffraction limit measured against
    /// the sampling of the sensor: a value well below one means the system is limited by its
    /// pixels and a finer sensor would resolve more, while a value above one means it is limited
    /// by diffraction and a finer sensor would not.
    ///
    /// # Arguments
    ///
    /// * `wavelength`: the wavelength of the light, in the same length unit as the camera.
    ///   `WAVELENGTH_550NM_MM` is the usual choice for visible light in millimeters.
    ///
    /// returns: f64
    pub fn airy_disk_px(&self, wavelength: f64) -> f64 {
        self.lens.airy_disk_diameter(wavelength) / self.sensor.pitch
    }
}

/// The default subdivision length for `Camera::project_outline`, in pixels.
///
/// This sets how finely a long outline edge is sampled before each piece is classified as visible
/// or hidden, so a transition between the two lands within this many pixels of where it belongs.
/// Four pixels is within the width of a drawn line once antialiasing is counted, and the cost is
/// small: the number of ray casts is the outline length in pixels divided by this value, which
/// for any real part is far below the one cast per pixel that `ViewBuffer::render` performs.
const DEFAULT_OUTLINE_EDGE_PX: f64 = 4.0;

/// Clip a segment expressed in camera space against the near plane at `z = near`, returning the
/// portion of it in front of that plane.
///
/// The returned endpoints keep the order of the inputs. A segment with a NaN coordinate fails
/// both side tests and is dropped, the same as one entirely behind the plane.
fn clip_near(a: &Point3, b: &Point3, near: f64) -> Option<(Point3, Point3)> {
    match (a.z >= near, b.z >= near) {
        (true, true) => Some((*a, *b)),
        (true, false) => Some((*a, lerp_to_near(a, b, near))),
        (false, true) => Some((lerp_to_near(b, a, near), *b)),
        (false, false) => None,
    }
}

/// Interpolate along the segment from `inside` to `outside` to the point where it crosses the
/// near plane.
///
/// The crossing's z is assigned rather than left to the interpolation, so that a point produced
/// here always passes the same `z >= near` test which produced it. Rounding in the interpolation
/// can otherwise leave it a few ulps short.
fn lerp_to_near(inside: &Point3, outside: &Point3, near: f64) -> Point3 {
    let t = (near - inside.z) / (outside.z - inside.z);
    let mut p = inside + (outside - inside) * t;
    p.z = near;
    p
}

/// Create the camera-to-world isometry for a camera positioned at `eye` whose optical axis
/// passes through `target`. The isometry orients the camera so that the `up` direction points
/// as nearly as possible toward the top of the image.
///
/// The result is in the camera convention used by this module, with +X right in the image, +Y
/// down, and +Z along the optical axis into the scene. The `up` argument does not need to be
/// perpendicular to the view direction; only its perpendicular component is used.
///
/// # Arguments
///
/// * `eye`: the position of the camera center in world coordinates
/// * `target`: a point in world coordinates which the optical axis passes through, which must
///   not coincide with `eye`
/// * `up`: the world direction to place toward the top of the image, which must not be parallel
///   to the direction from `eye` to `target`
///
/// returns: Result<Iso3>
pub fn look_at(eye: &Point3, target: &Point3, up: &Vector3) -> Result<Iso3> {
    // The camera's +Y is "down" in the image, which is the component of -up perpendicular to the
    // view direction. Its +X then falls out as right, because a rotation must be right-handed.
    Iso3::from_basis_zy(&(target - eye), &(-up), Some(*eye))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// A 2448 x 2048 sensor on a 3.45 micron pitch behind a 50mm f/2.8 lens focused at 1m, in
    /// millimeters. The image distance works out to 1000/19 of the focal length, so the object
    /// plane is 19 times the size of the sensor.
    fn test_camera() -> Camera {
        Camera::new(
            Sensor::new(2448, 2048, 0.00345).unwrap(),
            ThinLens::new(50.0, 2.8, 1000.0).unwrap(),
        )
    }

    #[test]
    fn intrinsics_use_the_image_distance() {
        let cam = test_camera();
        let pinhole = cam.pinhole();

        let expected = 52.631578947368425 / 0.00345;
        assert_relative_eq!(pinhole.fx, expected, epsilon = 1e-6);
        assert_relative_eq!(pinhole.fy, expected, epsilon = 1e-6);
        assert_relative_eq!(pinhole.cx, 1224.0, epsilon = 1e-12);
        assert_relative_eq!(pinhole.cy, 1024.0, epsilon = 1e-12);

        // Focusing close makes the camera project as though the lens were longer
        assert!(pinhole.fx > 50.0 / 0.00345);
    }

    #[test]
    fn intrinsics_reduce_to_the_focal_length_at_infinity_focus() {
        let cam = Camera::new(
            Sensor::new(2448, 2048, 0.00345).unwrap(),
            ThinLens::new(50.0, 2.8, f64::INFINITY).unwrap(),
        );
        let pinhole = cam.pinhole();
        let reference = PinholeCamera::from_focal_length(50.0 / 0.00345, 2448, 2048);

        assert_relative_eq!(pinhole.fx, reference.fx, epsilon = 1e-9);
        assert_relative_eq!(pinhole.fy, reference.fy, epsilon = 1e-9);
        assert_relative_eq!(pinhole.cx, reference.cx, epsilon = 1e-12);
        assert_relative_eq!(pinhole.cy, reference.cy, epsilon = 1e-12);

        assert!(cam.footprint().0.is_infinite());
        assert!(cam.pixel_size().is_infinite());
    }

    #[test]
    fn angles_of_view() {
        let cam = test_camera();
        let s_i = cam.lens.image_distance();

        assert_relative_eq!(
            cam.horizontal_fov(),
            2.0 * (8.4456f64 / (2.0 * s_i)).atan(),
            epsilon = 1e-12
        );
        assert_relative_eq!(cam.horizontal_fov(), 0.16012336, epsilon = 1e-7);
        assert!(cam.vertical_fov() < cam.horizontal_fov());
        assert!(cam.diagonal_fov() > cam.horizontal_fov());
    }

    #[test]
    fn footprint_and_pixel_size_scale_with_distance() {
        let cam = test_camera();

        // The object plane at 1m is 19 times the size of the sensor
        let (w, h) = cam.footprint();
        assert_relative_eq!(w, 160.4664, epsilon = 1e-9);
        assert_relative_eq!(h, 134.2464, epsilon = 1e-9);
        assert_relative_eq!(cam.pixel_size(), 0.06555, epsilon = 1e-12);

        // The footprint is the pixel size times the pixel count, at any distance
        let (w2, h2) = cam.footprint_at(2000.0).unwrap();
        let px = cam.pixel_size_at(2000.0).unwrap();
        assert_relative_eq!(w2, px * 2448.0, epsilon = 1e-9);
        assert_relative_eq!(h2, px * 2048.0, epsilon = 1e-9);

        // And it is linear in the distance
        assert_relative_eq!(w2, 2.0 * w, epsilon = 1e-9);

        assert!(cam.footprint_at(0.0).is_none());
        assert!(cam.footprint_at(-1.0).is_none());
        assert!(cam.pixel_size_at(f64::NAN).is_none());
    }

    #[test]
    fn footprint_agrees_with_projecting_its_own_corners() {
        // A rectangle of the footprint size, centered on the optical axis at the focus distance,
        // should project onto the corners of the sensor.
        let cam = test_camera();
        let (w, h) = cam.footprint();
        let z = cam.lens.focus_distance;
        let corner = Point3::new(-w / 2.0, -h / 2.0, z);

        let projected = cam.project(&corner, &Iso3::identity()).unwrap();
        assert_relative_eq!(projected.x, 0.0, epsilon = 1e-9);
        assert_relative_eq!(projected.y, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn blur_and_depth_of_field_in_pixels() {
        let cam = test_camera();

        // A blur tolerance of one pixel is a tolerance of one pitch
        let by_px = cam.depth_of_field_px(1.0).unwrap();
        let by_length = cam.lens.depth_of_field(0.00345).unwrap();
        assert_relative_eq!(by_px.0, by_length.0, epsilon = 1e-12);
        assert_relative_eq!(by_px.1, by_length.1, epsilon = 1e-12);
        assert_relative_eq!(
            cam.total_depth_of_field_px(1.0).unwrap(),
            by_px.1 - by_px.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            cam.hyperfocal_distance_px(1.0),
            cam.lens.hyperfocal_distance(0.00345),
            epsilon = 1e-12
        );

        // At the limits of a one pixel depth of field, the blur is one pixel
        assert_relative_eq!(cam.coc_px_at(by_px.0).unwrap(), 1.0, epsilon = 1e-9);
        assert_relative_eq!(cam.coc_px_at(by_px.1).unwrap(), 1.0, epsilon = 1e-9);

        assert_relative_eq!(cam.coc_px_at(1200.0).unwrap(), 45.4033634, epsilon = 1e-6);
        assert!(cam.coc_px_at(-1.0).is_none());
    }

    #[test]
    fn airy_disk_against_the_pixel_pitch() {
        let cam = test_camera();
        // At f/2.8 on a 3.45 micron pitch the Airy disk is a little over one pixel, so the
        // system is close to the crossover between diffraction and sampling.
        assert_relative_eq!(
            cam.airy_disk_px(WAVELENGTH_550NM_MM),
            1.14648,
            epsilon = 1e-5
        );

        // Stopping down makes diffraction dominate
        let stopped = Camera::new(cam.sensor, ThinLens::new(50.0, 16.0, 1000.0).unwrap());
        assert!(stopped.airy_disk_px(WAVELENGTH_550NM_MM) > 6.0);
    }

    #[test]
    fn from_field_width_selects_a_focal_length() {
        let sensor = Sensor::new(2448, 2048, 0.00345).unwrap();
        let cam = Camera::from_field_width(sensor, 2.8, 1000.0, 160.4664).unwrap();
        assert_relative_eq!(cam.lens.focal_length, 50.0, epsilon = 1e-9);

        // The camera it produces sees the field width it was asked for
        assert_relative_eq!(cam.footprint().0, 160.4664, epsilon = 1e-9);
    }

    /// A camera whose pinhole intrinsics are mathematically fx = fy = 1 px with cx = 100, cy = 50, so a
    /// camera space point (x, y, z) lands at (x/z + 100, y/z + 50). Its focal length is 1.0,
    /// which is also its default near plane distance.
    fn unit_camera() -> Camera {
        Camera::new(
            Sensor::new(200, 100, 1.0).unwrap(),
            ThinLens::new(1.0, 2.8, f64::INFINITY).unwrap(),
        )
    }

    #[test]
    fn project_segment_passes_a_fully_visible_segment_through() {
        let cam = unit_camera();
        let (a, b) = cam
            .project_segment(
                &Point3::new(-2.0, 0.0, 2.0),
                &Point3::new(2.0, 0.0, 2.0),
                &Iso3::identity(),
                Some(0.5),
            )
            .unwrap();
        assert_relative_eq!(a.x, 99.0, epsilon = 1e-12);
        assert_relative_eq!(b.x, 101.0, epsilon = 1e-12);
        assert_relative_eq!(a.y, 50.0, epsilon = 1e-12);
    }

    #[test]
    fn project_segment_clips_at_the_near_plane() {
        // From (0,0,4) to (8,0,-4) the crossing of z = 1 is at t = 3/8, camera point (3, 0, 1)
        let cam = unit_camera();
        let (a, b) = cam
            .project_segment(
                &Point3::new(0.0, 0.0, 4.0),
                &Point3::new(8.0, 0.0, -4.0),
                &Iso3::identity(),
                Some(1.0),
            )
            .unwrap();
        assert_relative_eq!(a.x, 100.0, epsilon = 1e-12);
        assert_relative_eq!(b.x, 103.0, epsilon = 1e-12);
    }

    #[test]
    fn project_segment_keeps_the_argument_order_when_the_first_end_is_clipped() {
        let cam = unit_camera();
        let forward = cam
            .project_segment(
                &Point3::new(0.0, 0.0, 4.0),
                &Point3::new(8.0, 0.0, -4.0),
                &Iso3::identity(),
                Some(1.0),
            )
            .unwrap();
        let reversed = cam
            .project_segment(
                &Point3::new(8.0, 0.0, -4.0),
                &Point3::new(0.0, 0.0, 4.0),
                &Iso3::identity(),
                Some(1.0),
            )
            .unwrap();
        assert_relative_eq!(forward.0.x, reversed.1.x, epsilon = 1e-12);
        assert_relative_eq!(forward.1.x, reversed.0.x, epsilon = 1e-12);
    }

    #[test]
    fn project_segment_rejects_a_segment_entirely_behind() {
        let cam = unit_camera();
        assert!(
            cam.project_segment(
                &Point3::new(-1.0, 0.0, -1.0),
                &Point3::new(1.0, 0.0, -1.0),
                &Iso3::identity(),
                Some(0.5),
            )
            .is_none()
        );
    }

    #[test]
    fn the_default_near_plane_is_the_focal_length() {
        // unit_camera has a focal length of 1.0, so a segment from z = 0.5 clips at z = 1
        let cam = unit_camera();
        let (a, _) = cam
            .project_segment(
                &Point3::new(0.5, 0.0, 0.5),
                &Point3::new(0.5, 0.0, 4.0),
                &Iso3::identity(),
                None,
            )
            .unwrap();
        assert_relative_eq!(a.x, 100.5, epsilon = 1e-9);

        // With an explicit nearer plane the whole segment survives, landing at x/z = 1.0
        let (a, _) = cam
            .project_segment(
                &Point3::new(0.5, 0.0, 0.5),
                &Point3::new(0.5, 0.0, 4.0),
                &Iso3::identity(),
                Some(0.25),
            )
            .unwrap();
        assert_relative_eq!(a.x, 101.0, epsilon = 1e-9);
    }

    #[test]
    fn an_unusable_near_value_falls_back_to_the_default() {
        let cam = unit_camera();
        let ends = (Point3::new(0.5, 0.0, 0.5), Point3::new(0.5, 0.0, 4.0));
        let reference = cam
            .project_segment(&ends.0, &ends.1, &Iso3::identity(), None)
            .unwrap();

        for bad in [Some(0.0), Some(-1.0), Some(f64::NAN), Some(f64::INFINITY)] {
            let got = cam
                .project_segment(&ends.0, &ends.1, &Iso3::identity(), bad)
                .unwrap();
            assert_relative_eq!(got.0.x, reference.0.x, epsilon = 1e-12);
        }
    }

    #[test]
    fn a_nan_point_is_treated_as_behind_the_camera() {
        let cam = unit_camera();
        assert!(
            cam.project_segment(
                &Point3::new(f64::NAN, 0.0, 2.0),
                &Point3::new(f64::NAN, 0.0, 2.0),
                &Iso3::identity(),
                Some(0.5),
            )
            .is_none()
        );
    }

    #[test]
    fn project_polyline_returns_one_run_when_nothing_is_clipped() {
        let cam = unit_camera();
        let pts = vec![
            Point3::new(-2.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(2.0, 0.0, 2.0),
            Point3::new(4.0, 0.0, 2.0),
        ];
        let runs = cam.project_polyline(&pts, &Iso3::identity(), Some(0.5));
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].len(), 4);
    }

    #[test]
    fn project_polyline_splits_into_runs_around_the_near_plane() {
        let cam = unit_camera();
        let pts = vec![
            Point3::new(0.0, 0.0, 4.0),
            Point3::new(1.0, 0.0, 4.0),
            Point3::new(2.0, 0.0, -4.0),
            Point3::new(3.0, 0.0, 4.0),
            Point3::new(4.0, 0.0, 4.0),
        ];
        let runs = cam.project_polyline(&pts, &Iso3::identity(), Some(1.0));
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].len(), 3);
        assert_eq!(runs[1].len(), 3);
    }

    #[test]
    fn project_polyline_handles_degenerate_input() {
        let cam = unit_camera();
        let empty: Vec<Point3> = Vec::new();
        assert!(
            cam.project_polyline(&empty, &Iso3::identity(), Some(0.5))
                .is_empty()
        );

        let visible = vec![Point3::new(0.0, 0.0, 2.0)];
        let runs = cam.project_polyline(&visible, &Iso3::identity(), Some(0.5));
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].len(), 1);

        let hidden = vec![Point3::new(0.0, 0.0, -2.0)];
        assert!(
            cam.project_polyline(&hidden, &Iso3::identity(), Some(0.5))
                .is_empty()
        );
    }

    #[test]
    fn project_polyline_agrees_with_project_segment() {
        // The two clip paths must not drift apart
        let cam = unit_camera();
        let iso = look_at(
            &Point3::new(3.0, 4.0, 5.0),
            &Point3::origin(),
            &Vector3::z(),
        )
        .unwrap();
        let pts = vec![
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.5),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, -0.5),
        ];
        let runs = cam.project_polyline(&pts, &iso, None);
        assert_eq!(runs.len(), 1);

        for (i, w) in pts.windows(2).enumerate() {
            let (a, _) = cam.project_segment(&w[0], &w[1], &iso, None).unwrap();
            assert_relative_eq!(runs[0][i].x, a.x, epsilon = 1e-9);
            assert_relative_eq!(runs[0][i].y, a.y, epsilon = 1e-9);
        }
    }

    #[test]
    fn project_polyline_agrees_with_project_on_unclipped_points() {
        // Guards the hoisted intrinsics against the per-call version
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, 0.0, 500.0),
            &Point3::origin(),
            &Vector3::y(),
        )
        .unwrap();
        let pts = vec![
            Point3::new(-10.0, -5.0, 0.0),
            Point3::new(10.0, -5.0, 3.0),
            Point3::new(10.0, 5.0, -3.0),
        ];
        let runs = cam.project_polyline(&pts, &iso, None);
        for (p, q) in pts.iter().zip(runs[0].iter()) {
            let direct = cam.project(p, &iso).unwrap();
            assert_relative_eq!(q.x, direct.x, epsilon = 1e-9);
            assert_relative_eq!(q.y, direct.y, epsilon = 1e-9);
        }
    }

    /// The exact image of a sphere is a circle. For radius r seen from distance d the silhouette
    /// subtends asin(r/d), so its pixel radius is f_px * r / sqrt(d^2 - r^2).
    #[test]
    fn the_projected_outline_of_a_sphere_lies_on_its_image_circle() -> Result<()> {
        let cam = Camera::new(
            Sensor::new(2000, 2000, 0.01).unwrap(),
            ThinLens::new(50.0, 2.8, f64::INFINITY).unwrap(),
        );
        let f_px = cam.pinhole().fx;
        let (r, d) = (10.0, 100.0);
        let mesh = Mesh3::create_sphere(r, 0.02)?;
        let iso = look_at(&Point3::new(0.0, 0.0, d), &Point3::origin(), &Vector3::y())?;

        let expected = f_px * r / (d * d - r * r).sqrt();
        let outline = cam.project_outline(&mesh, &iso, None, None)?;
        let visible: Vec<_> = outline.iter().filter(|(_, _, k)| *k == 0).collect();
        assert!(!visible.is_empty());

        for (a, _, _) in &visible {
            let radius = (a.x - 1000.0).hypot(a.y - 1000.0);
            assert_relative_eq!(radius, expected, epsilon = 4.0);
        }
        Ok(())
    }

    #[test]
    fn a_pixel_tolerance_scales_with_the_camera_it_is_given_to() -> Result<()> {
        // An edge length in pixels means pixels of that camera, so a half resolution preview
        // subdivides to twice the world length and produces about half as many segments. This is
        // the intended behavior: the preview stays cheap and its drawn line looks the same in
        // its own smaller image.
        let cam = Camera::new(
            Sensor::new(800, 800, 0.01).unwrap(),
            ThinLens::new(50.0, 2.8, 300.0).unwrap(),
        );
        let mesh = Mesh3::create_box(20.0, 20.0, 20.0, true);
        let iso = look_at(
            &Point3::new(60.0, 80.0, 100.0),
            &Point3::origin(),
            &Vector3::z(),
        )?;

        let full = cam.project_outline(&mesh, &iso, Some(8.0), None)?;
        let preview = cam
            .rescaled(0.5)?
            .project_outline(&mesh, &iso, Some(8.0), None)?;

        let ratio = preview.len() as f64 / full.len() as f64;
        assert!(
            ratio > 0.4 && ratio < 0.6,
            "expected about half, got {ratio}"
        );

        // Asking the preview for half the pixel length restores the original fineness
        let matched = cam
            .rescaled(0.5)?
            .project_outline(&mesh, &iso, Some(4.0), None)?;
        let matched_ratio = matched.len() as f64 / full.len() as f64;
        assert!(
            matched_ratio > 0.9 && matched_ratio < 1.1,
            "expected parity, got {matched_ratio}"
        );
        Ok(())
    }

    #[test]
    fn a_finer_tolerance_produces_more_segments() -> Result<()> {
        let cam = Camera::new(
            Sensor::new(800, 800, 0.01).unwrap(),
            ThinLens::new(50.0, 2.8, 300.0).unwrap(),
        );
        let mesh = Mesh3::create_box(20.0, 20.0, 20.0, true);
        let iso = look_at(
            &Point3::new(60.0, 80.0, 100.0),
            &Point3::origin(),
            &Vector3::z(),
        )?;

        let coarse = cam.project_outline(&mesh, &iso, Some(16.0), None)?;
        let fine = cam.project_outline(&mesh, &iso, Some(4.0), None)?;
        assert!(fine.len() > coarse.len() * 2);
        Ok(())
    }

    #[test]
    fn a_camera_looking_away_from_the_mesh_is_an_error() -> Result<()> {
        let cam = test_camera();
        let mesh = Mesh3::create_box(10.0, 10.0, 10.0, true);
        let eye = Point3::new(0.0, 0.0, 500.0);
        let iso = look_at(&eye, &Point3::new(0.0, 0.0, 1000.0), &Vector3::y())?;
        assert!(cam.project_outline(&mesh, &iso, None, None).is_err());
        Ok(())
    }

    #[test]
    fn a_non_positive_max_edge_px_is_rejected() -> Result<()> {
        let cam = test_camera();
        let mesh = Mesh3::create_box(10.0, 10.0, 10.0, true);
        let iso = look_at(
            &Point3::new(0.0, 0.0, 500.0),
            &Point3::origin(),
            &Vector3::y(),
        )?;
        assert!(cam.project_outline(&mesh, &iso, Some(0.0), None).is_err());
        assert!(cam.project_outline(&mesh, &iso, Some(-1.0), None).is_err());
        Ok(())
    }

    #[test]
    fn frustum_corners_project_back_to_the_image_corners() -> Result<()> {
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(3.0, 4.0, 5.0),
            &Point3::origin(),
            &Vector3::z(),
        )?;
        let corners = cam.frustum_corners(&iso, 100.0, 400.0)?;

        let w = cam.sensor.width_px as f64;
        let h = cam.sensor.height_px as f64;
        let expected = [
            Point2::new(0.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w, h),
            Point2::new(0.0, h),
        ];

        // Both faces project onto the same four image corners and therefore form a frustum.
        for (i, c) in corners.iter().enumerate() {
            let px = cam.project(c, &iso).unwrap();
            assert_relative_eq!(px.x, expected[i % 4].x, epsilon = 1e-6);
            assert_relative_eq!(px.y, expected[i % 4].y, epsilon = 1e-6);
        }

        // The near face sits at the requested depth along the optical axis
        let to_cam = iso.inverse();
        assert_relative_eq!((to_cam * corners[0]).z, 100.0, epsilon = 1e-9);
        assert_relative_eq!((to_cam * corners[4]).z, 400.0, epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn frustum_corners_reject_an_invalid_range() -> Result<()> {
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, 0.0, 500.0),
            &Point3::origin(),
            &Vector3::y(),
        )?;
        assert!(cam.frustum_corners(&iso, 0.0, 100.0).is_err());
        assert!(cam.frustum_corners(&iso, -1.0, 100.0).is_err());
        assert!(cam.frustum_corners(&iso, 100.0, 100.0).is_err());
        assert!(cam.frustum_corners(&iso, 100.0, 50.0).is_err());
        assert!(cam.frustum_corners(&iso, 100.0, f64::INFINITY).is_err());
        Ok(())
    }

    #[test]
    fn look_at_puts_the_target_at_the_principal_point() {
        let eye = Point3::new(3.0, 4.0, 5.0);
        let target = Point3::origin();
        let iso = look_at(&eye, &target, &Vector3::z()).unwrap();
        let cam = test_camera();
        let pinhole = cam.pinhole();

        let img = cam.project(&target, &iso).unwrap();
        assert_relative_eq!(img.x, pinhole.cx, epsilon = 1e-9);
        assert_relative_eq!(img.y, pinhole.cy, epsilon = 1e-9);

        // The ray through the principal point starts at the eye and looks at the target
        let ray = cam.back_project(&Point2::new(pinhole.cx, pinhole.cy), &iso);
        assert_relative_eq!(ray.origin, eye, epsilon = 1e-9);
        assert_relative_eq!(ray.dir, (target - eye).normalize(), epsilon = 1e-9);
    }

    #[test]
    fn look_at_orients_up_and_right_correctly() {
        let eye = Point3::new(3.0, 4.0, 5.0);
        let target = Point3::origin();
        let up = Vector3::z();
        let iso = look_at(&eye, &target, &up).unwrap();
        let cam = test_camera();
        let pinhole = cam.pinhole();

        // A point above the target lands above the center of the image, which is a smaller v,
        // and stays on the vertical center line
        let above = cam.project(&(target + up * 0.5), &iso).unwrap();
        assert!(above.y < pinhole.cy);
        assert_relative_eq!(above.x, pinhole.cx, epsilon = 1e-9);

        // A point displaced along (forward x up) lands to the right of center
        let right = (target - eye).cross(&up).normalize();
        let to_the_right = cam.project(&(target + right * 0.5), &iso).unwrap();
        assert!(to_the_right.x > pinhole.cx);
        assert_relative_eq!(to_the_right.y, pinhole.cy, epsilon = 1e-9);
    }

    #[test]
    fn look_at_frame_is_a_right_handed_rotation() {
        let iso = look_at(
            &Point3::new(3.0, 4.0, 5.0),
            &Point3::origin(),
            &Vector3::z(),
        )
        .unwrap();
        let m = iso.rotation.to_rotation_matrix();
        assert_relative_eq!(m.matrix().determinant(), 1.0, epsilon = 1e-12);

        // The camera z axis is the view direction
        let expected = (Point3::origin() - Point3::new(3.0, 4.0, 5.0)).normalize();
        assert_relative_eq!(iso.z().into_inner(), expected, epsilon = 1e-9);
    }

    #[test]
    fn look_at_rejects_degenerate_input() {
        let eye = Point3::new(3.0, 4.0, 5.0);
        let target = Point3::origin();
        assert!(look_at(&eye, &eye, &Vector3::z()).is_err());
        assert!(look_at(&eye, &target, &(target - eye)).is_err());
    }

    #[test]
    fn camera_pose_relates_to_a_plotting_view_by_a_half_turn() {
        // A plotting view is world-to-view with +Y up and +Z toward the viewer, so it is the
        // inverse of a camera pose composed with a half turn about X.
        let eye = Point3::new(3.0, 4.0, 5.0);
        let iso = look_at(&eye, &Point3::origin(), &Vector3::z()).unwrap();
        let view = Iso3::from_rx(PI) * iso.inverse();
        let round_trip = view.inverse() * Iso3::from_rx(PI);

        assert_relative_eq!(
            round_trip.translation.vector,
            iso.translation.vector,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            round_trip.z().into_inner(),
            iso.z().into_inner(),
            epsilon = 1e-9
        );

        // In the view convention, up in the world is up in the view
        let up_in_view = view * (Point3::origin() + Vector3::z()) - view * Point3::origin();
        assert!(up_in_view.y > 0.0);
    }
}
