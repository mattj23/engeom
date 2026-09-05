//! This module stores the per-pixel results of viewing scene geometry through a camera, linking
//! the optical model to the scene geometry.

use super::Camera;
use crate::common::IndexMask;
use crate::geom3::IsoExtensions3;
use crate::na::DMatrix;
use crate::raster2::{Point2I, RasterMask};
use crate::{Iso3, Mesh3, Point3, PointCloud3, SurfacePoint3, UnitVec3, Vector3};
use parry3d_f64::query::{Ray, RayCast};
use parry3d_f64::shape::FeatureId;
use rayon::prelude::*;

/// Data recorded by a single pixel. A pixel that saw nothing has a depth of `NaN`; its other
/// fields have no meaning.
#[derive(Debug, Clone, Copy)]
struct PixelData {
    /// Distance from the lens plane along the optical axis, not along the ray.
    depth: f64,

    /// World-space surface normal, oriented back toward the camera.
    normal: Vector3,

    /// Index of the triangle which was hit.
    face: u32,

    /// Whether the triangle was seen from the side its winding faces away from.
    backface: bool,
}

impl PixelData {
    fn miss() -> Self {
        Self {
            depth: f64::NAN,
            normal: Vector3::zeros(),
            face: u32::MAX,
            backface: false,
        }
    }

    fn is_hit(&self) -> bool {
        !self.depth.is_nan()
    }
}

/// A per-pixel view of a scene. The view is rendered by casting a ray through each pixel center
/// and recording the nearest surface that the ray strikes.
///
/// A `ViewBuffer` is a snapshot. It carries the camera and the pose it was rendered from, so
/// every quantity derived from it (defocus blur, incidence angle, world position, surface
/// coverage) is available without holding on to the scene it came from.
///
/// Pixel depth is measured from the lens plane along the optical axis rather than along the ray.
/// The thin-lens model uses this axial distance, so the stored depth can be passed directly to
/// [`Camera::coc_px_at`] and the depth-of-field methods.
///
/// # Cost
///
/// Rendering casts one ray per pixel, distributes image rows across threads, and uses roughly 48
/// bytes of storage per pixel. A large sensor at full resolution is therefore expensive in both
/// time and memory. For images other than the final result, a reduced resolution obtained from
/// [`Camera::rescaled`] is usually appropriate. Rescaling preserves the sensor's physical size,
/// so the preview and full-resolution render have the same field of view. See
/// [`Camera::rescaled`] for the warning about thresholds measured in pixels.
#[derive(Debug, Clone)]
pub struct ViewBuffer {
    camera: Camera,
    iso: Iso3,
    target_face_count: usize,
    pixels: Vec<PixelData>,
}

impl ViewBuffer {
    /// Render a scene by casting a ray through each camera pixel's center and recording the
    /// nearest target surface that the ray strikes.
    ///
    /// The camera optics do not participate in the ray cast. Each ray follows the pinhole
    /// intrinsics, which describe the path of an unaberrated ray through the aperture center.
    /// Defocus is calculated afterward from the depth, without tracing the full aperture. This
    /// calculation permits the renderer to use one ray per pixel.
    ///
    /// # Arguments
    ///
    /// * `camera`: the camera to look through, whose sensor sets the size of the buffer
    /// * `iso`: the camera-to-world transform placing the camera in the scene
    /// * `target`: the mesh whose surface is recorded
    /// * `obstruction`: an optional mesh that can block the target without being recorded. A
    ///   pixel whose view passes through the obstruction is left empty, like a pixel that saw no
    ///   surface.
    ///
    /// returns: ViewBuffer
    pub fn render(
        camera: &Camera,
        iso: &Iso3,
        target: &Mesh3,
        obstruction: Option<&Mesh3>,
    ) -> Self {
        let pinhole = camera.pinhole();
        let width = camera.sensor.width_px as usize;
        let height = camera.sensor.height_px as usize;
        let face_count = target.faces().len();
        let origin = iso.origin();

        // The rays and the meshes are both in world coordinates, so the shapes need no pose of
        // their own.
        let at_origin = Iso3::identity();

        let mut pixels = vec![PixelData::miss(); width * height];

        pixels
            .par_chunks_mut(width)
            .enumerate()
            .for_each(|(y, row)| {
                let dy = (y as f64 + 0.5 - pinhole.cy) / pinhole.fy;

                for (x, pixel) in row.iter_mut().enumerate() {
                    let dx = (x as f64 + 0.5 - pinhole.cx) / pinhole.fx;

                    // The camera-space direction through the pixel center, normalized so that
                    // the time of impact is a distance. Its z component is the reciprocal of the
                    // length, which converts that distance back into a depth along the axis.
                    let axis_fraction = 1.0 / (dx * dx + dy * dy + 1.0).sqrt();
                    let dir_cam =
                        Vector3::new(dx * axis_fraction, dy * axis_fraction, axis_fraction);
                    let ray = Ray::new(origin, iso * dir_cam);

                    let Some(hit) = target.tri_mesh().cast_ray_and_get_normal(
                        &at_origin,
                        &ray,
                        f64::MAX,
                        false,
                    ) else {
                        continue;
                    };

                    // The obstruction only matters nearer than the surface that was found, so
                    // the second cast is bounded by the first.
                    if let Some(ob) = obstruction
                        && ob
                            .tri_mesh()
                            .cast_ray(&at_origin, &ray, hit.time_of_impact, false)
                            .is_some()
                    {
                        continue;
                    }

                    let FeatureId::Face(feature) = hit.feature else {
                        continue;
                    };

                    // A trimesh reports a back face by offsetting the index by the face count,
                    // which is the only record of which side of the surface was seen.
                    let backface = feature as usize >= face_count;
                    let face = if backface {
                        feature - face_count as u32
                    } else {
                        feature
                    };

                    *pixel = PixelData {
                        depth: hit.time_of_impact * axis_fraction,
                        normal: hit.normal.normalize(),
                        face,
                        backface,
                    };
                }
            });

        Self {
            camera: *camera,
            iso: *iso,
            target_face_count: face_count,
            pixels,
        }
    }

    /// The camera this view was rendered through.
    pub fn camera(&self) -> &Camera {
        &self.camera
    }

    /// The camera-to-world transform this view was rendered from.
    pub fn iso(&self) -> &Iso3 {
        &self.iso
    }

    /// The width of the buffer in pixels.
    pub fn width(&self) -> u32 {
        self.camera.sensor.width_px
    }

    /// The height of the buffer in pixels.
    pub fn height(&self) -> u32 {
        self.camera.sensor.height_px
    }

    /// The number of faces in the target mesh which was rendered, which sets the length of the
    /// masks returned by `face_mask` and `face_mask_where`.
    pub fn target_face_count(&self) -> usize {
        self.target_face_count
    }

    /// Convert a pixel coordinate into an index into the buffer, or `None` if it lies outside.
    fn index_of(&self, p: Point2I) -> Option<usize> {
        if p.x < 0 || p.y < 0 || p.x >= self.width() as i32 || p.y >= self.height() as i32 {
            None
        } else {
            Some(p.y as usize * self.width() as usize + p.x as usize)
        }
    }

    /// What a single pixel saw, or `None` if the pixel lies outside the buffer or saw nothing.
    ///
    /// # Arguments
    ///
    /// * `p`: the pixel coordinate
    ///
    /// returns: Option<PixelHit>
    pub fn at(&self, p: Point2I) -> Option<PixelHit<'_>> {
        let index = self.index_of(p)?;
        if self.pixels[index].is_hit() {
            Some(PixelHit {
                buffer: self,
                index,
            })
        } else {
            None
        }
    }

    /// The depth recorded at a single pixel, which is `NaN` if the pixel saw nothing or lies
    /// outside the buffer.
    ///
    /// # Arguments
    ///
    /// * `p`: the pixel coordinate
    ///
    /// returns: f64
    pub fn depth_at(&self, p: Point2I) -> f64 {
        self.index_of(p)
            .map(|i| self.pixels[i].depth)
            .unwrap_or(f64::NAN)
    }

    /// Iterate over every pixel which saw a surface.
    pub fn iter_hits(&self) -> impl Iterator<Item = PixelHit<'_>> {
        self.pixels
            .iter()
            .enumerate()
            .filter(|(_, p)| p.is_hit())
            .map(move |(index, _)| PixelHit {
                buffer: self,
                index,
            })
    }

    /// The number of pixels which saw a surface.
    pub fn hit_count(&self) -> usize {
        self.pixels.iter().filter(|p| p.is_hit()).count()
    }

    /// The fraction of the sensor that saw a surface, from zero to one. This value measures how
    /// much of the frame the part fills and whether the pose wastes image area.
    pub fn hit_fraction(&self) -> f64 {
        self.hit_count() as f64 / self.pixels.len() as f64
    }

    /// Return the target silhouette as a raster mask, with every pixel that saw a surface set to
    /// true.
    ///
    /// This mask connects a rendered view to the `raster2` tools. Their mask and distance
    /// operations can compare the shape seen by the camera with a segmented image.
    pub fn hit_mask(&self) -> RasterMask {
        let width = self.width() as usize;
        let mut mask = RasterMask::empty(self.width(), self.height());
        for (index, pixel) in self.pixels.iter().enumerate() {
            if pixel.is_hit() {
                let p = Point2I::new((index % width) as i32, (index / width) as i32);
                mask.set_point_unchecked(p, true);
            }
        }
        mask
    }

    /// Build a matrix with one value per pixel and `NaN` wherever a pixel saw nothing. Each
    /// matrix row corresponds to an image row, as required by the `raster2` tools. The result can
    /// be passed to `raster2::render_d_matrix` to write an image.
    ///
    /// # Arguments
    ///
    /// * `f`: computes the value for a pixel which saw a surface
    ///
    /// returns: DMatrix<f64>
    pub fn to_matrix(&self, f: impl Fn(&PixelHit) -> f64) -> DMatrix<f64> {
        let mut matrix =
            DMatrix::from_element(self.height() as usize, self.width() as usize, f64::NAN);
        let width = self.width() as usize;

        for hit in self.iter_hits() {
            let row = hit.index / width;
            let col = hit.index % width;
            matrix[(row, col)] = f(&hit);
        }

        matrix
    }

    /// A matrix of the depth at each pixel, with `NaN` wherever the pixel saw nothing.
    pub fn to_depth_matrix(&self) -> DMatrix<f64> {
        self.to_matrix(|h| h.depth())
    }

    /// A matrix of the defocus-blur diameter in pixels at each pixel, with `NaN` wherever the
    /// pixel saw nothing. This matrix identifies the parts of the frame that are in focus.
    pub fn to_coc_px_matrix(&self) -> DMatrix<f64> {
        self.to_matrix(|h| h.coc_px())
    }

    /// A matrix of the incidence angle in radians at each pixel, with `NaN` wherever the pixel
    /// saw nothing. Zero means the camera is looking straight into the surface and a right angle
    /// means it is grazing along it.
    pub fn to_incidence_matrix(&self) -> DMatrix<f64> {
        self.to_matrix(|h| h.incidence())
    }

    /// Collect every pixel that saw a surface into a point cloud of world-space points and
    /// normals. The result contains the same data that a depth sensor at this pose would produce.
    pub fn to_point_cloud(&self) -> PointCloud3 {
        let points: Vec<SurfacePoint3> = self.iter_hits().map(|h| h.surface_point()).collect();
        PointCloud3::from_surface_points(&points)
    }

    /// A mask over the faces of the target mesh marking every face which was seen by at least
    /// one pixel.
    pub fn face_mask(&self) -> IndexMask {
        self.face_mask_where(|_| true)
    }

    /// A mask over the target mesh faces, marking every face which at least one pixel saw well
    /// enough to satisfy all of the supplied criteria.
    ///
    /// This numeric form of `face_mask_where` supports the thresholds on recorded quantities that
    /// are commonly used for shot planning. A criterion set to `None` is not applied. Calling
    /// this method with no thresholds and `allow_backface` set to `true` is therefore equivalent
    /// to calling `face_mask`.
    ///
    /// A pixel whose blur or incidence is not a number fails any threshold on that quantity.
    ///
    /// Blur is measured in the pixels of the camera this view was rendered through, so a
    /// threshold applies to the resolution it was written for. See `Camera::rescaled`.
    ///
    /// # Arguments
    ///
    /// * `max_coc_px`: the largest acceptable defocus blur diameter in pixels; a pixel qualifies
    ///   when its blur is at or below this. If `None`, blur is not considered.
    /// * `max_incidence`: the largest acceptable angle in radians between the surface normal and
    ///   the direction back to the camera; a pixel qualifies when its incidence is at or below
    ///   this. If `None`, incidence is not considered.
    /// * `allow_backface`: whether a pixel which saw the face from the side its winding faces
    ///   away from may qualify. Pass `false` to require that the camera saw the outside of the
    ///   surface.
    ///
    /// returns: IndexMask
    pub fn face_mask_usable(
        &self,
        max_coc_px: Option<f64>,
        max_incidence: Option<f64>,
        allow_backface: bool,
    ) -> IndexMask {
        self.face_mask_where(|h| {
            (allow_backface || !h.is_backface())
                && max_coc_px.is_none_or(|m| h.coc_px() <= m)
                && max_incidence.is_none_or(|m| h.incidence() <= m)
        })
    }

    /// A mask over the target mesh faces. A face is marked when at least one pixel that sees it
    /// satisfies the predicate. For example, the predicate can require pixels to be in focus,
    /// below an incidence-angle limit, or both, to identify the faces usefully captured by a
    /// pose.
    ///
    /// # Arguments
    ///
    /// * `predicate`: decides whether a pixel counts as having usefully seen its face
    ///
    /// returns: IndexMask
    pub fn face_mask_where(&self, predicate: impl Fn(&PixelHit) -> bool) -> IndexMask {
        let mut mask = IndexMask::new(self.target_face_count, false);
        for hit in self.iter_hits() {
            if predicate(&hit) {
                mask.set(hit.face() as usize, true);
            }
        }
        mask
    }
}

/// A borrowed view of the data recorded by one [`ViewBuffer`] pixel. Each quantity is derived on
/// demand from the stored depth and normal and from the camera and pose used to render the
/// buffer.
#[derive(Debug, Clone, Copy)]
pub struct PixelHit<'a> {
    buffer: &'a ViewBuffer,
    index: usize,
}

impl PixelHit<'_> {
    fn data(&self) -> &PixelData {
        &self.buffer.pixels[self.index]
    }

    /// The coordinate of this pixel in the image.
    pub fn pixel(&self) -> Point2I {
        let width = self.buffer.width() as usize;
        Point2I::new((self.index % width) as i32, (self.index / width) as i32)
    }

    /// The distance from the lens plane to the surface, measured along the optical axis.
    pub fn depth(&self) -> f64 {
        self.data().depth
    }

    /// The index of the face of the target mesh which was seen.
    pub fn face(&self) -> u32 {
        self.data().face
    }

    /// Whether the face was seen from the side its winding faces away from, which on a closed
    /// mesh means the camera is inside the part and on an open one means it is behind the
    /// surface.
    pub fn is_backface(&self) -> bool {
        self.data().backface
    }

    /// The surface normal in world space, oriented back toward the camera. On a back face this
    /// is the reverse of the face's own normal.
    pub fn normal(&self) -> UnitVec3 {
        UnitVec3::new_normalize(self.data().normal)
    }

    /// The position of the surface in world space.
    pub fn point(&self) -> Point3 {
        let pinhole = self.buffer.camera.pinhole();
        let p = self.pixel();
        let z = self.depth();
        let x = (p.x as f64 + 0.5 - pinhole.cx) / pinhole.fx * z;
        let y = (p.y as f64 + 0.5 - pinhole.cy) / pinhole.fy * z;
        self.buffer.iso * Point3::new(x, y, z)
    }

    /// The position and normal of the surface in world space.
    pub fn surface_point(&self) -> SurfacePoint3 {
        SurfacePoint3::new_normalize(self.point(), self.data().normal)
    }

    /// The diameter of the defocus blur of this point, in pixels of the camera this view was
    /// rendered through. A view rendered with a camera from [`Camera::rescaled`] reports the
    /// diameter in its own coarser pixels. A threshold on this value therefore applies to the
    /// resolution for which it was defined.
    pub fn coc_px(&self) -> f64 {
        self.buffer
            .camera
            .coc_px_at(self.depth())
            .unwrap_or(f64::NAN)
    }

    /// The angle in radians between the surface normal and the direction back to the camera.
    /// Zero means the camera is looking straight into the surface, and a right angle means it is
    /// grazing along it. Data captured at a high incidence angle is generally less trustworthy,
    /// so this angle helps determine whether a pose captured the surface usefully.
    pub fn incidence(&self) -> f64 {
        let to_camera = self.buffer.iso.origin() - self.point();
        self.data().normal.angle(&to_camera)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sensors::camera::{Sensor, ThinLens, look_at};
    use crate::{Point2, Point3};
    use approx::assert_relative_eq;

    /// A small camera, coarse enough to render quickly in a test but with an odd pixel count so
    /// that a center pixel exists.
    fn test_camera() -> Camera {
        Camera::new(
            Sensor::new(101, 101, 0.01).unwrap(),
            ThinLens::new(50.0, 2.8, 1000.0).unwrap(),
        )
    }

    /// A camera with the same optics but a much larger sensor, giving a wide enough angle of
    /// view that the difference between a depth and a ray length is obvious.
    fn wide_camera() -> Camera {
        Camera::new(
            Sensor::new(101, 101, 0.5).unwrap(),
            ThinLens::new(50.0, 2.8, 1000.0).unwrap(),
        )
    }

    /// A single square in the z = 0 plane spanning [-size, size] in x and y, wound so that its
    /// normals point toward +z.
    fn square(size: f64) -> Mesh3 {
        let vertices = vec![
            Point3::new(-size, -size, 0.0),
            Point3::new(size, -size, 0.0),
            Point3::new(size, size, 0.0),
            Point3::new(-size, size, 0.0),
        ];
        Mesh3::new(vertices, vec![[0, 1, 2], [0, 2, 3]], false)
    }

    /// A camera on the +z axis looking down at the z = 0 plane from `height`, with world +y up
    /// in the image.
    fn looking_down(height: f64) -> Iso3 {
        look_at(
            &Point3::new(0.0, 0.0, height),
            &Point3::origin(),
            &Vector3::y(),
        )
        .unwrap()
    }

    #[test]
    fn depth_is_measured_along_the_optical_axis() {
        // A camera looking straight down at a flat plane should record the same depth at every
        // pixel, even at the corners where the ray is longer than the depth.
        let cam = wide_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1000.0), None);

        assert_eq!(buffer.hit_count(), 101 * 101);
        assert_relative_eq!(buffer.hit_fraction(), 1.0, epsilon = 1e-12);

        let center = buffer.depth_at(Point2I::new(50, 50));
        let corner = buffer.depth_at(Point2I::new(0, 0));
        assert_relative_eq!(center, 500.0, epsilon = 1e-9);
        assert_relative_eq!(corner, 500.0, epsilon = 1e-9);

        // The distance along the ray to the corner is much longer, confirming that ray distance
        // and axial depth are different quantities.
        let corner_point = buffer.at(Point2I::new(0, 0)).unwrap().point();
        let eye = Point3::new(0.0, 0.0, 500.0);
        assert!((corner_point - eye).norm() > 550.0);
    }

    #[test]
    fn misses_are_recorded_as_nan() {
        // A square smaller than the field of view leaves the edges of the frame empty
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1.0), None);

        assert!(buffer.at(Point2I::new(50, 50)).is_some());
        assert!(buffer.depth_at(Point2I::new(0, 0)).is_nan());
        assert!(buffer.at(Point2I::new(0, 0)).is_none());

        assert!(buffer.hit_count() > 0);
        assert!(buffer.hit_count() < 101 * 101);

        // Out of bounds behaves the same as an empty pixel
        assert!(buffer.depth_at(Point2I::new(-1, 0)).is_nan());
        assert!(buffer.at(Point2I::new(101, 0)).is_none());
        assert!(buffer.at(Point2I::new(0, -1)).is_none());
    }

    #[test]
    fn world_points_land_on_the_surface_and_reproject_to_their_own_pixel() {
        let cam = test_camera();
        let iso = looking_down(500.0);
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);

        for p in [
            Point2I::new(50, 50),
            Point2I::new(0, 0),
            Point2I::new(100, 17),
        ] {
            let hit = buffer.at(p).unwrap();

            // The surface is the z = 0 plane
            assert_relative_eq!(hit.point().z, 0.0, epsilon = 1e-9);

            // And projecting it back through the camera returns the center of its own pixel
            let px = cam.project(&hit.point(), &iso).unwrap();
            assert_relative_eq!(px.x, p.x as f64 + 0.5, epsilon = 1e-6);
            assert_relative_eq!(px.y, p.y as f64 + 0.5, epsilon = 1e-6);
        }
    }

    #[test]
    fn normals_face_the_camera_on_a_front_face() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1000.0), None);
        let hit = buffer.at(Point2I::new(50, 50)).unwrap();

        // The square is wound toward +z and the camera is above it, so this is a front face and
        // the normal is the face's own
        assert!(!hit.is_backface());
        assert_relative_eq!(hit.normal().into_inner(), Vector3::z(), epsilon = 1e-9);

        // Looking straight down at it, the incidence angle is zero
        assert_relative_eq!(hit.incidence(), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn normals_are_flipped_toward_the_camera_on_a_back_face() {
        // The same square seen from below is a back face. Its normal is reported reversed so
        // that it still faces the camera, and the flag records that it was reversed.
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(-500.0), &square(1000.0), None);
        let hit = buffer.at(Point2I::new(50, 50)).unwrap();

        assert!(hit.is_backface());
        assert_relative_eq!(hit.normal().into_inner(), -Vector3::z(), epsilon = 1e-9);
        assert_relative_eq!(hit.incidence(), 0.0, epsilon = 1e-9);

        // The face index remains a valid mesh index after removing the offset that a trimesh
        // uses to report a back face.
        assert!((hit.face() as usize) < buffer.target_face_count());
        assert_relative_eq!(hit.depth(), 500.0, epsilon = 1e-9);
    }

    #[test]
    fn hit_mask_marks_the_pixels_that_saw_the_surface() {
        // A square that covers only part of the frame: the mask must agree with the hit count,
        // be true where the square is, and be false in a corner that looks past it.
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(2.0), None);
        let mask = buffer.hit_mask();

        assert_eq!(mask.width(), buffer.width());
        assert_eq!(mask.height(), buffer.height());
        assert_eq!(mask.count_true(), buffer.hit_count());
        assert!(buffer.hit_count() > 0);
        assert!(mask.get_point(Point2I::new(50, 50)));
        assert!(!mask.get_point(Point2I::new(0, 0)));
    }

    #[test]
    fn incidence_grows_as_the_surface_tilts_away() {
        // Looking at the plane from 45 degrees off its normal
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, -500.0, 500.0),
            &Point3::origin(),
            &Vector3::z(),
        )
        .unwrap();
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);
        let hit = buffer.at(Point2I::new(50, 50)).unwrap();

        assert_relative_eq!(hit.incidence(), std::f64::consts::FRAC_PI_4, epsilon = 1e-9);
    }

    #[test]
    fn nearer_surfaces_win() {
        // Two stacked planes: only the nearer one should be recorded
        let cam = test_camera();
        let mut near = square(1000.0);
        near.transform_in_place(&Iso3::from_translation(0.0, 0.0, 100.0));

        let mut both = square(1000.0);
        both.append_in_place(&near).unwrap();

        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &both, None);
        assert_relative_eq!(buffer.depth_at(Point2I::new(50, 50)), 400.0, epsilon = 1e-9);
    }

    #[test]
    fn an_obstruction_blocks_without_being_recorded() {
        let cam = test_camera();
        let iso = looking_down(500.0);
        let target = square(1000.0);

        // A small blocker halfway between the camera and the target
        let mut blocker = square(2.0);
        blocker.transform_in_place(&Iso3::from_translation(0.0, 0.0, 250.0));

        let clear = ViewBuffer::render(&cam, &iso, &target, None);
        let blocked = ViewBuffer::render(&cam, &iso, &target, Some(&blocker));

        // The center is lost to the blocker
        assert!(clear.at(Point2I::new(50, 50)).is_some());
        assert!(blocked.at(Point2I::new(50, 50)).is_none());

        // The corners still see the target at its own depth, not the blocker's
        assert_relative_eq!(blocked.depth_at(Point2I::new(0, 0)), 500.0, epsilon = 1e-9);
        assert!(blocked.hit_count() < clear.hit_count());

        // A blocker behind the target changes nothing
        let mut behind = square(2.0);
        behind.transform_in_place(&Iso3::from_translation(0.0, 0.0, -250.0));
        let unblocked = ViewBuffer::render(&cam, &iso, &target, Some(&behind));
        assert_eq!(unblocked.hit_count(), clear.hit_count());
    }

    #[test]
    fn defocus_blur_follows_the_depth() {
        // A camera focused at 500 looking at a plane at 500 sees no blur; the same camera
        // focused elsewhere sees the blur its own optics predict.
        let cam = Camera::new(
            Sensor::new(21, 21, 0.01).unwrap(),
            ThinLens::new(50.0, 2.8, 500.0).unwrap(),
        );
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1000.0), None);
        let hit = buffer.at(Point2I::new(10, 10)).unwrap();
        assert_relative_eq!(hit.coc_px(), 0.0, epsilon = 1e-9);

        let far = ViewBuffer::render(&cam, &looking_down(600.0), &square(1000.0), None);
        let far_hit = far.at(Point2I::new(10, 10)).unwrap();
        assert_relative_eq!(
            far_hit.coc_px(),
            cam.coc_px_at(600.0).unwrap(),
            epsilon = 1e-12
        );
        assert!(far_hit.coc_px() > 0.0);
    }

    #[test]
    fn matrices_are_laid_out_by_row_and_mark_misses() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1.0), None);

        let depth = buffer.to_depth_matrix();
        assert_eq!(depth.nrows(), 101);
        assert_eq!(depth.ncols(), 101);
        assert_relative_eq!(depth[(50, 50)], 500.0, epsilon = 1e-9);
        assert!(depth[(0, 0)].is_nan());

        // The other matrices agree with the per-pixel values
        let hit = buffer.at(Point2I::new(50, 50)).unwrap();
        assert_relative_eq!(
            buffer.to_coc_px_matrix()[(50, 50)],
            hit.coc_px(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            buffer.to_incidence_matrix()[(50, 50)],
            hit.incidence(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn pixel_coordinates_round_trip_through_the_iterator() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1.0), None);

        assert_eq!(buffer.iter_hits().count(), buffer.hit_count());

        for hit in buffer.iter_hits() {
            let p = hit.pixel();
            assert_relative_eq!(buffer.depth_at(p), hit.depth(), epsilon = 1e-12);
            assert!(p.x >= 0 && p.x < 101 && p.y >= 0 && p.y < 101);
        }
    }

    #[test]
    fn point_cloud_carries_the_hits_with_their_normals() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(500.0), &square(1.0), None);
        let cloud = buffer.to_point_cloud();

        assert_eq!(cloud.point_count(), buffer.hit_count());
        let normals = cloud.point_normals().expect("cloud should carry normals");
        assert_eq!(normals.len(), cloud.point_count());

        for (p, n) in cloud.points().iter().zip(normals.iter()) {
            assert_relative_eq!(p.z, 0.0, epsilon = 1e-9);
            assert_relative_eq!(n.into_inner(), Vector3::z(), epsilon = 1e-9);
        }
    }

    #[test]
    fn face_masks_report_coverage() {
        let cam = test_camera();
        let target = square(1000.0);

        // Rendered from the camera's own focus distance, so the view is genuinely in focus
        let buffer = ViewBuffer::render(&cam, &looking_down(1000.0), &target, None);

        // Both triangles of the square fill the frame, so both are seen
        let mask = buffer.face_mask();
        assert_eq!(mask.len(), 2);
        assert_eq!(mask.count_true(), 2);

        // A predicate no pixel satisfies leaves nothing covered
        assert_eq!(buffer.face_mask_where(|_| false).count_true(), 0);

        // And one that only the sharpest pixels satisfy still covers this in-focus view
        let sharp = buffer.face_mask_where(|h| h.coc_px() < 1.0 && h.incidence() < 0.1);
        assert_eq!(sharp.count_true(), 2);
    }

    #[test]
    fn no_criteria_reproduces_the_plain_face_mask() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(1000.0), &square(1000.0), None);
        let plain = buffer.face_mask();
        let usable = buffer.face_mask_usable(None, None, true);
        assert_eq!(plain.count_true(), usable.count_true());
        assert_eq!(usable.count_true(), 2);
    }

    #[test]
    fn a_blur_threshold_excludes_an_out_of_focus_view() {
        // test_camera is focused at 1000, so rendering the plane there is sharp and rendering it
        // at 500 is not
        let cam = test_camera();
        let sharp = ViewBuffer::render(&cam, &looking_down(1000.0), &square(1000.0), None);
        assert_eq!(
            sharp.face_mask_usable(Some(1.0), None, false).count_true(),
            2
        );

        let blurred = ViewBuffer::render(&cam, &looking_down(500.0), &square(1000.0), None);
        assert!(blurred.at(Point2I::new(50, 50)).unwrap().coc_px() > 1.0);
        assert_eq!(
            blurred
                .face_mask_usable(Some(1.0), None, false)
                .count_true(),
            0
        );
    }

    #[test]
    fn a_zero_blur_threshold_still_marks_a_perfectly_focused_face() {
        // The comparison is inclusive, so asking for zero blur returns the pixels which have
        // none rather than returning nothing at all.
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(1000.0), &square(1000.0), None);
        assert_eq!(
            buffer.face_mask_usable(Some(0.0), None, false).count_true(),
            2
        );
    }

    #[test]
    fn an_incidence_threshold_excludes_a_grazing_view() {
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, -1000.0, 1000.0),
            &Point3::origin(),
            &Vector3::z(),
        )
        .unwrap();
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);

        let limit = std::f64::consts::FRAC_PI_4;
        assert_eq!(
            buffer
                .face_mask_usable(None, Some(limit + 0.01), true)
                .count_true(),
            2
        );
        assert_eq!(
            buffer
                .face_mask_usable(None, Some(limit - 0.01), true)
                .count_true(),
            0
        );
    }

    #[test]
    fn backfaces_are_excluded_unless_allowed() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(-1000.0), &square(1000.0), None);
        assert!(buffer.at(Point2I::new(50, 50)).unwrap().is_backface());
        assert_eq!(buffer.face_mask_usable(None, None, false).count_true(), 0);
        assert_eq!(buffer.face_mask_usable(None, None, true).count_true(), 2);
    }

    #[test]
    fn the_numeric_mask_agrees_with_the_equivalent_predicate() {
        let cam = test_camera();
        let buffer = ViewBuffer::render(&cam, &looking_down(900.0), &square(1000.0), None);

        for (coc, inc, back) in [
            (Some(5.0), None, true),
            (None, Some(0.5), false),
            (Some(20.0), Some(1.0), false),
            (None, None, true),
        ] {
            let numeric = buffer.face_mask_usable(coc, inc, back);
            let predicate = buffer.face_mask_where(|h| {
                (back || !h.is_backface())
                    && coc.is_none_or(|m| h.coc_px() <= m)
                    && inc.is_none_or(|m| h.incidence() <= m)
            });
            assert_eq!(numeric.to_indices(), predicate.to_indices());
        }
    }

    #[test]
    fn an_empty_view_gives_an_empty_usable_mask() {
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, 0.0, 500.0),
            &Point3::new(0.0, 0.0, 1000.0),
            &Vector3::y(),
        )
        .unwrap();
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);
        let mask = buffer.face_mask_usable(None, None, true);
        assert_eq!(mask.len(), 2);
        assert_eq!(mask.count_true(), 0);
    }

    #[test]
    fn coverage_accumulates_across_poses() {
        // Two poses that each see half of a long strip provide complete coverage together, as
        // required for shot planning.
        let cam = test_camera();
        let mut strip = square(20.0);
        strip.transform_in_place(&Iso3::from_translation(60.0, 0.0, 0.0));
        let mut target = square(20.0);
        target.transform_in_place(&Iso3::from_translation(-60.0, 0.0, 0.0));
        target.append_in_place(&strip).unwrap();

        let left = ViewBuffer::render(
            &cam,
            &look_at(
                &Point3::new(-60.0, 0.0, 500.0),
                &Point3::new(-60.0, 0.0, 0.0),
                &Vector3::y(),
            )
            .unwrap(),
            &target,
            None,
        );
        let right = ViewBuffer::render(
            &cam,
            &look_at(
                &Point3::new(60.0, 0.0, 500.0),
                &Point3::new(60.0, 0.0, 0.0),
                &Vector3::y(),
            )
            .unwrap(),
            &target,
            None,
        );

        let left_mask = left.face_mask();
        let right_mask = right.face_mask();
        assert_eq!(left_mask.len(), 4);

        // Each pose provides partial coverage, and their union provides complete coverage.
        assert!(left_mask.count_true() < 4);
        assert!(right_mask.count_true() < 4);
        assert_eq!(left_mask.or(&right_mask).unwrap().count_true(), 4);
    }

    #[test]
    fn a_rescaled_camera_previews_the_same_view() {
        let cam = test_camera();
        let iso = looking_down(500.0);
        let target = square(3.0);

        let full = ViewBuffer::render(&cam, &iso, &target, None);
        let preview = ViewBuffer::render(&cam.rescaled(0.2).unwrap(), &iso, &target, None);

        assert_eq!(preview.width(), 20);
        assert_eq!(preview.height(), 20);

        // The field of view is preserved, so the part fills the same fraction of the frame, to
        // within the coarser sampling of the preview along the edges of the part
        assert_relative_eq!(preview.hit_fraction(), full.hit_fraction(), epsilon = 0.05);

        // And the depth it reports is the same
        assert_relative_eq!(
            preview.depth_at(Point2I::new(10, 10)),
            full.depth_at(Point2I::new(50, 50)),
            epsilon = 1e-9
        );
    }

    #[test]
    fn blur_in_pixels_is_tied_to_the_resolution_it_was_measured_at() {
        // A rescaled camera has larger pixels, so the same physical defocus covers fewer of
        // them. Each reading applies to its own sensor, and the values differ by the scale
        // factor. A pixel threshold therefore cannot be transferred between the sensors
        // unchanged.
        let cam = test_camera();
        let factor = 0.2;
        let preview = cam.rescaled(factor).unwrap();

        let full_blur = cam.coc_px_at(500.0).unwrap();
        let preview_blur = preview.coc_px_at(500.0).unwrap();
        assert_relative_eq!(preview_blur, full_blur * factor, epsilon = 1e-12);
        assert!(preview_blur < full_blur);

        // The physical blur is the same for both, being a property of the optics alone
        assert_relative_eq!(
            cam.lens.coc_diameter(500.0).unwrap(),
            preview.lens.coc_diameter(500.0).unwrap(),
            epsilon = 1e-12
        );

        // And a rendered view reports the blur of the camera it was rendered through
        let target = square(1000.0);
        let iso = looking_down(500.0);
        let full_view = ViewBuffer::render(&cam, &iso, &target, None);
        let preview_view = ViewBuffer::render(&preview, &iso, &target, None);
        assert_relative_eq!(
            full_view.at(Point2I::new(50, 50)).unwrap().coc_px(),
            full_blur,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            preview_view.at(Point2I::new(10, 10)).unwrap().coc_px(),
            preview_blur,
            epsilon = 1e-9
        );
    }

    #[test]
    fn buffer_keeps_the_camera_and_pose_it_was_rendered_from() {
        let cam = test_camera();
        let iso = looking_down(500.0);
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);

        assert_eq!(*buffer.camera(), cam);
        assert_relative_eq!(buffer.iso().origin(), iso.origin(), epsilon = 1e-12);
        assert_eq!(buffer.width(), 101);
        assert_eq!(buffer.height(), 101);
        assert_eq!(buffer.target_face_count(), 2);
    }

    #[test]
    fn an_empty_view_is_empty_everywhere() {
        // Pointed away from the part
        let cam = test_camera();
        let iso = look_at(
            &Point3::new(0.0, 0.0, 500.0),
            &Point3::new(0.0, 0.0, 1000.0),
            &Vector3::y(),
        )
        .unwrap();
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);

        assert_eq!(buffer.hit_count(), 0);
        assert_relative_eq!(buffer.hit_fraction(), 0.0, epsilon = 1e-12);
        assert_eq!(buffer.iter_hits().count(), 0);
        assert_eq!(buffer.to_point_cloud().point_count(), 0);
        assert_eq!(buffer.face_mask().count_true(), 0);
        assert!(buffer.to_depth_matrix().iter().all(|v| v.is_nan()));
    }

    #[test]
    fn rescaling_preserves_the_physical_sensor() {
        let sensor = Sensor::new(2448, 2048, 0.00345).unwrap();
        let half = sensor.rescaled(0.5).unwrap();

        assert_eq!(half.width_px, 1224);
        assert_eq!(half.height_px, 1024);
        assert_relative_eq!(half.width(), sensor.width(), epsilon = 1e-12);
        assert_relative_eq!(half.height(), sensor.height(), epsilon = 1e-12);

        // Which means the camera built on it sees the same field of view
        let lens = ThinLens::new(50.0, 2.8, 1000.0).unwrap();
        assert_relative_eq!(
            Camera::new(half, lens).horizontal_fov(),
            Camera::new(sensor, lens).horizontal_fov(),
            epsilon = 1e-12
        );

        // Rounding never produces a sensor with no pixels
        assert_eq!(sensor.rescaled(1e-9).unwrap().width_px, 1);
        assert!(sensor.rescaled(0.0).is_err());
        assert!(sensor.rescaled(-1.0).is_err());
        assert!(sensor.rescaled(f64::NAN).is_err());
    }

    #[test]
    fn a_point_projects_to_the_pixel_that_saw_it() {
        // The whole loop: pick a world point on the part, project it to a pixel, and confirm the
        // buffer recorded that same point at that pixel.
        let cam = test_camera();
        let iso = looking_down(500.0);
        let buffer = ViewBuffer::render(&cam, &iso, &square(1000.0), None);

        let world = Point3::new(3.0, -2.0, 0.0);
        let px: Point2 = cam.project(&world, &iso).unwrap();
        let hit = buffer.at(Point2I::new(px.x as i32, px.y as i32)).unwrap();

        // Within half a pixel of the point we asked about, since the buffer samples centers
        let px_size = cam.pixel_size_at(500.0).unwrap();
        assert!((hit.point() - world).norm() < px_size);
    }
}
