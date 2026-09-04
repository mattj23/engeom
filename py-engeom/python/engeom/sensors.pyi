from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy
from numpy.typing import NDArray

from .bounding import Aabb2
from .common import IndexMask
from .geom2 import Point2
from .geom3 import Point3, Mesh3, Iso3, Vector3, PointCloud3, RayBundle3

WAVELENGTH_550NM_MM: float
"""The conventional representative wavelength for visible green light, in millimeters."""


def look_at(eye: Point3, target: Point3, up: Vector3) -> Iso3:
    """
    Build the camera-to-world isometry for a camera at ``eye`` whose optical axis passes through
    ``target``. The isometry orients the camera so that ``up`` points as nearly as possible toward
    the top of the image.

    The result is in the camera convention used by this module: +X right in the image, +Y down,
    and +Z along the optical axis into the scene. This is not the plotting view convention, which
    has +Y up and +Z toward the viewer.

    :param eye: the position of the camera center in world coordinates.
    :param target: a point the optical axis passes through, which must not coincide with ``eye``.
    :param up: the world direction to place toward the top of the image, which must not be
        parallel to the direction from ``eye`` to ``target``.
    :return: the camera-to-world isometry.
    :raises ValueError: if ``eye`` and ``target`` coincide, or ``up`` is parallel to the view
        direction.
    """
    ...


class Sensor:
    """
    A digital imaging sensor: a rectangular grid of square pixels with a known physical pitch.

    Pixel coordinates are continuous rather than integer, and the pixel with integer index ``i``
    covers ``[i, i + 1)``, so the center of the image lies at ``(width_px / 2, height_px / 2)``.

    The caller chooses the length unit, which is usually millimeters. All physical dimensions use
    that same unit.
    """

    def __init__(self, width_px: int, height_px: int, pitch: float):
        """
        Create a sensor from its pixel dimensions and physical pixel pitch.
        :param width_px: the number of pixel columns, which must not be zero.
        :param height_px: the number of pixel rows, which must not be zero.
        :param pitch: the physical center-to-center distance between adjacent pixels, which must
            be a finite number greater than zero.
        :raises ValueError: if either dimension is zero or the pitch is not positive and finite.
        """
        ...

    @property
    def width_px(self) -> int:
        """The number of pixel columns."""
        ...

    @property
    def height_px(self) -> int:
        """The number of pixel rows."""
        ...

    @property
    def pitch(self) -> float:
        """The physical center-to-center distance between adjacent pixels."""
        ...

    @property
    def width(self) -> float:
        """The physical width of the active sensor area."""
        ...

    @property
    def height(self) -> float:
        """The physical height of the active sensor area."""
        ...

    @property
    def diagonal(self) -> float:
        """The physical diagonal of the active sensor area."""
        ...

    @property
    def aspect_ratio(self) -> float:
        """The ratio of width to height, the same in pixels or physical units."""
        ...

    @property
    def pixel_count(self) -> int:
        """The total number of pixels in the sensor."""
        ...

    @property
    def center_px(self) -> Point2:
        """The center of the image in continuous pixel coordinates."""
        ...

    @property
    def image_aabb(self) -> Aabb2:
        """
        The region of pixel coordinates covered by the sensor, as a closed box from ``(0, 0)`` to
        ``(width_px, height_px)``. It includes its far edges, while ``contains_pixel`` does not.
        """
        ...

    def rescaled(self, factor: float) -> Sensor:
        """
        Create a copy with scaled pixel dimensions and an inversely scaled pitch. This preserves
        the physical sensor size and the field of view of a camera built with it.
        :param factor: the factor to scale the pixel dimensions by, finite and greater than zero.
        :return: the rescaled sensor.
        :raises ValueError: if the factor is not positive and finite.
        """
        ...

    def contains_pixel(self, pixel: Point2) -> bool:
        """
        Whether a point in continuous pixel coordinates is recorded by a pixel of this sensor.

        Half open on the far edges, so a point at ``x = width_px`` is not contained even though it
        lies on the boundary of ``image_aabb``.
        :param pixel: a point in continuous pixel coordinates.
        :return: True if a pixel of the sensor records the point.
        """
        ...

    def pixel_to_plane(self, pixel: Point2) -> Point2:
        """
        Convert pixel coordinates to a physical position on the sensor plane, measured from the
        center of the sensor with +x right and +y down.
        :param pixel: a point in continuous pixel coordinates.
        :return: the physical position relative to the sensor center.
        """
        ...

    def plane_to_pixel(self, plane: Point2) -> Point2:
        """
        The inverse of ``pixel_to_plane``.
        :param plane: a physical position relative to the sensor center.
        :return: the point in continuous pixel coordinates.
        """
        ...


class LaserProfile:
    def __init__(
            self,
            emitter_z: float,
            detector_y: float,
            detector_z: float,
            volume_width: float,
            volume_z_min: float,
            volume_z_max: float,
            resolution: int,
            angle_limit: float | None = None,
    ):
        """
        Create the base geometry of a laser profile line sensor, which emits a laser line into a
        scene and detects the reflection of that line to triangulate the distance to points on a
        surface.
       
        The general coordinate system is specified in X and Z. The center of the detection volume
        is at the origin, with the laser line ranging from the -X direction to the +X direction.
        The +Z direction points directly up towards the emitter.  The +Y direction is orthogonal to
        laser line and is typically the direction which the sensor will be panned.
       
        The geometry is specified with the following assumptions:
          - The laser line is emitted from a point directly on the +Z axis, with no offset in
            the X or Y direction.
          - The detector is not offset in the X direction, and can be specified with a Y and
            Z offset from the center of the detection volume.
          - The detection volume is trapezoidal, and its flat top and bottom are specified by a
            maximum and minimum Z value.
          - The detection volume's with is specified at Z=0, and is symmetrical around X=0.
       
        # Arguments
       
        :param emitter_z: The Z coordinate of the laser emitter. This is the height from the volume
          center where the laser fans into a triangle.
        :param detector_y: The Y coordinate of the detector's optical center. This is the out-of-plane
          offset from the plane of the laser line.
        :param detector_z: The Z coordinate of the detector's optical center. This is the height from
          the volume center where the detector's optical center is located.
        :param volume_width: The width of the detection volume at Z=0. The volume is assumed to be
          symmetrical around the X axis, ranging from -volume_width/2 to +volume_width/2.
        :param volume_z_min: The minimum Z value of the detection volume. This is the bottom of the
          trapezoidal volume, the farthest distance from the emitter where the sensor will still
          return points.
        :param volume_z_max: The maximum Z value of the detection volume. This is the top of the
          trapezoidal volume, the closest distance to the emitter where the sensor will still
          return points.
        :param resolution: The number of rays to cast across the laser line. This is the number of
          points that will be returned in the point cloud.
        :param angle_limit: An optional angle limit in radians. If specified, the sensor will only
          return a point if the angle between the surface normal at the point and the detector is
          less than this limit.
        """
        ...

    def get_points(self, target: Mesh3, obstruction: Mesh3 | None, iso: Iso3) -> PointCloud3:
        """

        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...


class PanningLaserProfile:
    def __init__(self, laser_line: LaserProfile, y_step: float, steps: int):
        """
        :param laser_line:
        :param y_step:
        :param steps:
        """
        ...

    def get_points(self, target: Mesh3, obstruction: Mesh3 | None, iso: Iso3) -> PointCloud3:
        """
        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...


class ThinLens:
    """
    A thin lens with a circular aperture, focused so that objects at ``focus_distance`` are imaged
    sharply onto the sensor plane.

    All distances are measured from the lens plane along the optical axis. A lens focused at
    infinity is represented by a ``focus_distance`` of ``float('inf')``, and every method handles
    that without special treatment by the caller.
    """

    def __init__(self, focal_length: float, f_number: float, focus_distance: float):
        """
        Create a thin lens.
        :param focal_length: the focal length, finite and greater than zero.
        :param f_number: the f-number at infinity focus, finite and greater than zero.
        :param focus_distance: the object-side distance to the plane of best focus, which must be
            greater than the focal length. ``float('inf')`` means focused at infinity.
        :raises ValueError: if any argument is out of range, in particular a focus distance at or
            inside the focal length, which forms no real image.
        """
        ...

    @staticmethod
    def from_magnification(focal_length: float, f_number: float, magnification: float) -> ThinLens:
        """
        Create a lens focused to achieve a specified magnification, as commonly used to specify a
        macro system.
        :param focal_length: the focal length, finite and greater than zero.
        :param f_number: the f-number at infinity focus, finite and greater than zero.
        :param magnification: the ratio of image size to object size, finite and greater than zero.
        :return: the lens focused to that magnification.
        :raises ValueError: if any argument is out of range.
        """
        ...

    @staticmethod
    def from_field_width(
            f_number: float,
            focus_distance: float,
            sensor_width: float,
            field_width: float,
    ) -> ThinLens:
        """
        Create the lens whose focal length fills a sensor of a given width with a given width of
        the object plane at a given working distance.
        :param f_number: the f-number at infinity focus.
        :param focus_distance: the object-side working distance, finite and greater than zero.
        :param sensor_width: the physical width of the sensor.
        :param field_width: the width of the object plane to image across that sensor width.
        :return: the lens with the solved focal length.
        :raises ValueError: if any argument is not positive and finite.
        """
        ...

    @property
    def focal_length(self) -> float:
        """The focal length of the lens."""
        ...

    @property
    def f_number(self) -> float:
        """The f-number of the lens at infinity focus."""
        ...

    @property
    def focus_distance(self) -> float:
        """The object-side distance to the plane of best focus. May be infinite."""
        ...

    @property
    def image_distance(self) -> float:
        """
        The image-side distance from the lens plane to the sensor plane. A lens focused at
        infinity returns its focal length; any closer focus returns something longer.
        """
        ...

    @property
    def magnification(self) -> float:
        """The ratio of image size to object size at the focus distance. Zero at infinity focus."""
        ...

    @property
    def effective_f_number(self) -> float:
        """
        The f-number at the focus distance, ``N * (1 + m)``. At 1:1 magnification this is twice
        the marked f-number, and it is the value governing both diffraction and exposure.
        """
        ...

    @property
    def aperture_diameter(self) -> float:
        """The diameter of the entrance pupil, the focal length divided by the f-number."""
        ...

    @property
    def is_focused_at_infinity(self) -> bool:
        """Whether the lens is focused at infinity."""
        ...

    def focused_at(self, focus_distance: float) -> ThinLens:
        """
        Create a copy refocused to a different distance, leaving the focal length and f-number
        unchanged.
        :param focus_distance: the new object-side distance to the plane of best focus.
        :return: the refocused lens.
        :raises ValueError: if the distance is at or inside the focal length.
        """
        ...

    def hyperfocal_distance(self, coc: float) -> float:
        """
        The focus distance beyond which the far limit of the depth of field reaches infinity.
        :param coc: the largest acceptable circle of confusion diameter on the sensor plane. Zero
            or less returns infinity.
        :return: the hyperfocal distance.
        """
        ...

    def airy_disk_diameter(self, wavelength: float) -> float:
        """
        The diameter of the Airy disk on the sensor plane, to the first zero, using the effective
        f-number so the result is correct at high magnification.
        :param wavelength: the wavelength of the light, in the same length unit as the lens.
        :return: the Airy disk diameter.
        """
        ...

    def coc_diameter(self, z: float) -> float | None:
        """
        The diameter of the defocus blur on the sensor plane for a point at a given distance.
        :param z: the object-side distance, greater than zero. Infinity is allowed.
        :return: the blur diameter, or None if ``z`` is not greater than zero.
        """
        ...

    def depth_of_field(self, coc: float) -> Tuple[float, float] | None:
        """
        The near and far limits of the depth of field. The far limit is infinite when the lens is
        focused at or beyond its hyperfocal distance.
        :param coc: the largest acceptable circle of confusion diameter, zero or greater.
        :return: the near and far limits, or None if ``coc`` is negative or NaN.
        """
        ...

    def total_depth_of_field(self, coc: float) -> float | None:
        """
        The distance between the near and far limits of the depth of field.
        :param coc: the largest acceptable circle of confusion diameter, zero or greater.
        :return: the total depth of field, or None if ``coc`` is negative or NaN.
        """
        ...


class PinholeCamera:
    """
    The pixel-unit intrinsics of a camera: focal lengths and a principal point, in pixels.

    This is the low-level projection model. Prefer ``Camera`` when physical optics are available;
    it derives these intrinsics from the sensor and lens.
    """

    def __init__(self, fx: float, fy: float, cx: float, cy: float):
        """
        Create intrinsics from focal lengths and a principal point, all in pixels.
        :param fx: the horizontal focal length in pixels.
        :param fy: the vertical focal length in pixels.
        :param cx: the horizontal coordinate of the principal point in pixels.
        :param cy: the vertical coordinate of the principal point in pixels.
        """
        ...

    @staticmethod
    def from_focal_length(focal_length: float, width: int, height: int) -> PinholeCamera:
        """
        Create intrinsics with equal focal lengths and the principal point at the image center.
        :param focal_length: the focal length in pixels, applied to both axes.
        :param width: the image width in pixels.
        :param height: the image height in pixels.
        :return: the intrinsics.
        """
        ...

    @property
    def fx(self) -> float:
        """The horizontal focal length in pixels."""
        ...

    @property
    def fy(self) -> float:
        """The vertical focal length in pixels."""
        ...

    @property
    def cx(self) -> float:
        """The horizontal coordinate of the principal point in pixels."""
        ...

    @property
    def cy(self) -> float:
        """The vertical coordinate of the principal point in pixels."""
        ...


class Camera:
    """
    A physical camera: a ``Sensor`` behind a ``ThinLens``.

    The pose is not stored. Every operation needing one takes a camera-to-world ``Iso3``, so a
    single camera can be evaluated at many poses. Build a pose with ``look_at``.

    The camera convention is +X right in the image, +Y down, and +Z along the optical axis into
    the scene, which is not the plotting view convention.
    """

    def __init__(self, sensor: Sensor, lens: ThinLens):
        """
        Create a camera from a sensor and a lens.
        :param sensor: the discrete sensor which records the image.
        :param lens: the lens which forms the image on the sensor.
        """
        ...

    @staticmethod
    def from_field_width(
            sensor: Sensor,
            f_number: float,
            working_distance: float,
            field_width: float,
    ) -> Camera:
        """
        Create a camera whose lens fills the sensor width with a specified width of the object
        plane at a given working distance.
        :param sensor: the discrete sensor which records the image.
        :param f_number: the f-number of the lens at infinity focus.
        :param working_distance: the object-side distance to the plane of best focus.
        :param field_width: the width of the object plane to image across the sensor width.
        :return: the camera with the solved focal length.
        :raises ValueError: if any argument is not positive and finite.
        """
        ...

    @property
    def sensor(self) -> Sensor:
        """The sensor which records the image."""
        ...

    @property
    def lens(self) -> ThinLens:
        """The lens which forms the image."""
        ...

    @property
    def pinhole(self) -> PinholeCamera:
        """
        The pixel-unit intrinsics. The focal lengths come from the lens's image distance divided
        by the pixel pitch, not from the focal length, so a lens focused close projects as though
        it were longer than its nameplate.
        """
        ...

    @property
    def horizontal_fov(self) -> float:
        """The full horizontal angle of view, in radians."""
        ...

    @property
    def vertical_fov(self) -> float:
        """The full vertical angle of view, in radians."""
        ...

    @property
    def diagonal_fov(self) -> float:
        """The full diagonal angle of view, in radians."""
        ...

    @property
    def footprint(self) -> Tuple[float, float]:
        """
        The width and height of the region of the focus plane imaged onto the full sensor. A
        camera focused at infinity returns infinities.
        """
        ...

    @property
    def pixel_size(self) -> float:
        """
        The size on the focus plane of the region imaged onto a single pixel. A camera focused at
        infinity returns infinity.
        """
        ...

    def rescaled(self, factor: float) -> Camera:
        """
        Create a copy with a scaled sensor resolution for rendering a less expensive preview. The
        physical sensor size is preserved, so the field of view is unchanged.

        Take care with quantities denominated in pixels: ``coc_px_at`` on a copy scaled by a
        factor returns the original value multiplied by that factor. A focus criterion carried
        unchanged from a preview to a full render is therefore much more permissive than intended.
        Scale the threshold too, or define the criterion with ``ThinLens.coc_diameter``, which
        returns a length.
        :param factor: the factor to scale the pixel dimensions by, finite and greater than zero.
        :return: the rescaled camera.
        :raises ValueError: if the factor is not positive and finite.
        """
        ...

    def footprint_at(self, z: float) -> Tuple[float, float] | None:
        """
        The width and height of the region of a plane at distance ``z`` imaged onto the sensor.
        :param z: the object-side distance, greater than zero.
        :return: the width and height, or None if ``z`` is not greater than zero.
        """
        ...

    def pixel_size_at(self, z: float) -> float | None:
        """
        The size on a plane at distance ``z`` of the region imaged onto a single pixel, which is
        the sampling resolution of the camera at that distance.
        :param z: the object-side distance, greater than zero.
        :return: the object-space pixel size, or None if ``z`` is not greater than zero.
        """
        ...

    def coc_px_at(self, z: float) -> float | None:
        """
        The diameter of the defocus blur of a point at distance ``z``, in pixels of this camera.
        :param z: the object-side distance, greater than zero.
        :return: the blur diameter in pixels, or None if ``z`` is not greater than zero.
        """
        ...

    def depth_of_field_px(self, coc_px: float) -> Tuple[float, float] | None:
        """
        The near and far limits of the depth of field, for a blur tolerance in pixels. Stating the
        tolerance in pixels ties the depth of field to the sampling of the sensor.
        :param coc_px: the largest acceptable blur diameter in pixels, zero or greater.
        :return: the near and far limits, or None if ``coc_px`` is negative or NaN.
        """
        ...

    def total_depth_of_field_px(self, coc_px: float) -> float | None:
        """
        The distance between the near and far limits, for a blur tolerance in pixels.
        :param coc_px: the largest acceptable blur diameter in pixels, zero or greater.
        :return: the total depth of field, or None if ``coc_px`` is negative or NaN.
        """
        ...

    def hyperfocal_distance_px(self, coc_px: float) -> float:
        """
        The hyperfocal distance for a blur tolerance in pixels.
        :param coc_px: the largest acceptable blur diameter in pixels.
        :return: the hyperfocal distance.
        """
        ...

    def airy_disk_px(self, wavelength: float) -> float:
        """
        The diameter of the Airy disk in pixels. Well below one means the system is limited by its
        pixels and a finer sensor would resolve more; above one means it is limited by diffraction
        and a finer sensor would not.
        :param wavelength: the wavelength of the light, in the same length unit as the camera.
        :return: the Airy disk diameter in pixels.
        """
        ...

    def project(self, points: NDArray[float], iso: Iso3) -> NDArray[float]:
        """
        Project world points into the image.

        Points behind the camera come back as rows of NaN, so the result lines up row for row with
        the input. A projected point can lie outside the sensor bounds, meaning it is within the
        optical field but not recorded; test with ``Sensor.contains_pixel``.
        :param points: an ``(n, 3)`` array of world coordinates.
        :param iso: the camera-to-world transform.
        :return: an ``(n, 2)`` array of pixel coordinates, NaN where the point was behind.
        """
        ...

    def back_project(self, pixels: NDArray[float], iso: Iso3) -> RayBundle3:
        """
        Back-project image points into world space rays which can be cast at a mesh.
        :param pixels: an ``(n, 2)`` array of pixel coordinates.
        :param iso: the camera-to-world transform.
        :return: a bundle of rays originating at the camera center.
        :raises ValueError: if the array is not ``(n, 2)``.
        """
        ...

    def project_polyline(
            self,
            points: NDArray[float],
            iso: Iso3,
            near: float | None = None,
    ) -> List[NDArray[float]]:
        """
        Project a world-space polyline into the image, split into runs wherever it passes behind
        the near plane.

        The result contains raw point runs because a polyline seen nearly edge-on
        can project to fewer than two distinct pixels.
        :param points: an ``(n, 3)`` array of polyline vertices in order.
        :param iso: the camera-to-world transform.
        :param near: the camera-space distance to the near plane. If None, or not positive and
            finite, the lens focal length is used, because an object nearer than the focal length
            forms no real image.
        :return: a list of ``(m, 2)`` pixel-coordinate arrays, one per visible run.
        """
        ...

    def project_outline(
            self,
            mesh: Mesh3,
            iso: Iso3,
            max_edge_px: float | None = None,
            corner_angle: float | None = None,
    ) -> Tuple[NDArray[float], NDArray[numpy.uint8]]:
        """
        Compute a line drawing of a mesh in this camera's pixel space, with hidden lines marked.

        The subdivision length is in pixels of this camera, so a preview from ``rescaled`` gets
        proportionally fewer and coarser segments, reducing the preview cost.
        :param mesh: the mesh to outline.
        :param iso: the camera-to-world transform.
        :param max_edge_px: the greatest length of a returned segment in pixels. If None, four.
        :param corner_angle: the smallest angle in radians between adjacent faces for their shared
            edge to be drawn as a corner. If None, just under 45 degrees.
        :return: an ``(n, 4)`` array of ``x0, y0, x1, y1`` in pixels, and an ``(n,)`` array of
            visibility codes where 0 is visible and 1 is hidden.
        :raises ValueError: if the mesh center is not in front of the camera, or ``max_edge_px``
            is not positive.
        """
        ...

    def frustum_corners(self, iso: Iso3, near: float, far: float) -> NDArray[float]:
        """
        The eight corners of the view volume between two depths along the optical axis.

        Passing the limits from ``depth_of_field_px`` gives the volume in which a part is both
        framed and in focus.
        :param iso: the camera-to-world transform.
        :param near: the depth of the near face, greater than zero.
        :param far: the depth of the far face, finite and greater than ``near``.
        :return: an ``(8, 3)`` array, the near face first, each running from the top left of the
            image clockwise.
        :raises ValueError: if the near or far distance is out of range.
        """
        ...


class ViewBuffer:
    """
    A per-pixel view of a scene, rendered by casting a ray through each pixel center and recording
    the nearest surface that the ray strikes.

    A ``ViewBuffer`` is a snapshot. It carries the camera and pose it was rendered from, so every
    derived quantity is available without holding on to the scene.

    Rendering a large sensor at full resolution is expensive in time and memory. A reduced
    resolution from ``Camera.rescaled`` is usually appropriate for previews.
    """

    def __init__(
            self,
            camera: Camera,
            iso: Iso3,
            target: Mesh3,
            obstruction: Mesh3 | None = None,
    ):
        """
        Render a view of a scene.
        :param camera: the camera to look through, whose sensor sets the size of the buffer.
        :param iso: the camera-to-world transform placing the camera in the scene.
        :param target: the mesh whose surface is recorded.
        :param obstruction: an optional mesh which can block the view of the target without being
            recorded itself. A pixel whose view passes through it is left empty.
        """
        ...

    @property
    def camera(self) -> Camera:
        """The camera this view was rendered through."""
        ...

    @property
    def iso(self) -> Iso3:
        """The camera-to-world transform this view was rendered from."""
        ...

    @property
    def width(self) -> int:
        """The width of the buffer in pixels."""
        ...

    @property
    def height(self) -> int:
        """The height of the buffer in pixels."""
        ...

    @property
    def target_face_count(self) -> int:
        """The number of faces in the target mesh, the length of the masks this returns."""
        ...

    @property
    def hit_count(self) -> int:
        """The number of pixels which saw a surface."""
        ...

    @property
    def hit_fraction(self) -> float:
        """
        The fraction of the sensor which saw a surface, from zero to one. This is how much of the
        frame the part fills.
        """
        ...

    @property
    def depth(self) -> NDArray[float]:
        """
        A ``(height, width)`` array of the distance from the lens plane along the optical axis,
        NaN where the pixel saw nothing. This is the quantity the thin lens model is written in
        terms of, so it can be passed to ``Camera.coc_px_at`` directly.
        """
        ...

    @property
    def coc_px(self) -> NDArray[float]:
        """
        A ``(height, width)`` array of the defocus-blur diameter in pixels of this view's camera,
        with NaN where the pixel saw nothing. This array identifies which parts of the frame are
        in focus.
        """
        ...

    @property
    def incidence(self) -> NDArray[float]:
        """
        A ``(height, width)`` array of the angle in radians between the surface normal and the
        direction back to the camera, NaN where the pixel saw nothing. Zero means looking straight
        into the surface and a right angle means grazing along it.
        """
        ...

    def depth_at(self, x: int, y: int) -> float:
        """
        The depth recorded at a single pixel.
        :param x: the pixel column.
        :param y: the pixel row.
        :return: the depth, or NaN if the pixel saw nothing or lies outside the buffer.
        """
        ...

    def to_point_cloud(self) -> PointCloud3:
        """
        Collect every pixel that saw a surface into a point cloud of world points and normals.

        Note that this samples the image plane uniformly, not the surface, so density is biased
        toward near and face-on geometry. Use ``face_mask_usable`` for coverage accounting.
        :return: the point cloud.
        """
        ...

    @property
    def face_mask(self) -> IndexMask:
        """
        A mask over the target mesh faces marking every face seen by at least one pixel.
        :return: the mask, of length ``target_face_count``.
        """
        ...

    def face_mask_usable(
            self,
            max_coc_px: float | None = None,
            max_incidence: float | None = None,
            allow_backface: bool = False,
    ) -> IndexMask:
        """
        A mask over the target mesh faces. A face is marked when at least one pixel saw it well
        enough to satisfy every supplied criterion, identifying what the pose usefully captured.

        A criterion left as None is not applied. A pixel whose blur or incidence is NaN fails any
        threshold on that quantity.
        :param max_coc_px: the largest acceptable blur diameter in pixels of this view's camera.
        :param max_incidence: the largest acceptable incidence angle in radians.
        :param allow_backface: whether a pixel which saw the face from behind may qualify.
        :return: the mask, of length ``target_face_count``.
        """
        ...
