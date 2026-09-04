"""
Tests for the camera and optics bindings.

These focus on what the binding layer is responsible for: the numpy shapes, dtypes and
orientations crossing the boundary, the mapping of Rust `Result` and `Option` onto `ValueError`
and `None`, and the property versus method split. The optics and the rendering themselves are
tested on the Rust side.

The depth array is the one place a silent bug could hide: `DMatrix` is column major and numpy is
row major, so a scene that is deliberately not symmetric is used to catch a transpose.
"""

from __future__ import annotations

import math

import numpy
import pytest

from engeom.geom2 import Point2
from engeom.geom3 import Iso3, Mesh3, Point3, Vector3
from engeom.sensors import (
    WAVELENGTH_550NM_MM,
    Camera,
    PinholeCamera,
    Sensor,
    ThinLens,
    ViewBuffer,
    look_at,
)


def machine_vision_sensor() -> Sensor:
    """A common 5MP sensor: 2448 x 2048 on a 3.45 micron pitch, in mm."""
    return Sensor(2448, 2048, 0.00345)


def small_camera() -> Camera:
    """Small enough to render quickly, with an odd pixel count so a center pixel exists."""
    return Camera(Sensor(41, 41, 0.01), ThinLens(50.0, 2.8, 500.0))


def flat_square(size: float) -> Mesh3:
    """A square in the z = 0 plane wound so its normals point toward +z."""
    return Mesh3(
        numpy.array(
            [[-size, -size, 0.0], [size, -size, 0.0], [size, size, 0.0], [-size, size, 0.0]]
        ),
        numpy.array([[0, 1, 2], [0, 2, 3]], dtype=numpy.uint32),
    )


def looking_down(height: float) -> Iso3:
    return look_at(Point3(0.0, 0.0, height), Point3(0.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0))


# ================================================================================================
# Construction and validation
# ================================================================================================


def test_constructors_reject_invalid_arguments():
    with pytest.raises(ValueError):
        Sensor(0, 100, 0.01)
    with pytest.raises(ValueError):
        Sensor(100, 100, 0.0)
    with pytest.raises(ValueError):
        ThinLens(50.0, 2.8, 40.0)
    with pytest.raises(ValueError):
        ThinLens(50.0, 2.8, 50.0)
    with pytest.raises(ValueError):
        ThinLens(-50.0, 2.8, 500.0)
    with pytest.raises(ValueError):
        ThinLens.from_magnification(50.0, 2.8, 0.0)
    with pytest.raises(ValueError):
        Sensor(100, 100, 0.01).rescaled(0.0)


def test_a_lens_focused_at_infinity_is_representable():
    lens = ThinLens(50.0, 2.8, float("inf"))
    assert lens.is_focused_at_infinity
    assert lens.image_distance == pytest.approx(50.0)
    assert lens.magnification == pytest.approx(0.0)
    assert math.isinf(Camera(machine_vision_sensor(), lens).footprint[0])


def test_optional_returns_become_none():
    lens = ThinLens(50.0, 2.8, 500.0)
    assert lens.coc_diameter(-1.0) is None
    assert lens.depth_of_field(-1.0) is None
    cam = small_camera()
    assert cam.footprint_at(0.0) is None
    assert cam.pixel_size_at(-1.0) is None
    assert cam.coc_px_at(0.0) is None


def test_look_at_rejects_degenerate_input():
    eye = Point3(3.0, 4.0, 5.0)
    with pytest.raises(ValueError):
        look_at(eye, eye, Vector3(0.0, 0.0, 1.0))


def test_the_wavelength_constant_is_exposed():
    assert WAVELENGTH_550NM_MM == pytest.approx(0.55e-3)


# ================================================================================================
# Accessors: what is a property and what is a method
# ================================================================================================


def test_state_and_cheap_derived_values_are_properties():
    sensor = machine_vision_sensor()
    assert sensor.width_px == 2448
    assert sensor.height_px == 2048
    assert sensor.pitch == pytest.approx(0.00345)
    assert sensor.width == pytest.approx(8.4456)
    assert sensor.height == pytest.approx(7.0656)
    assert sensor.pixel_count == 2448 * 2048
    assert isinstance(sensor.center_px, Point2)

    lens = ThinLens(50.0, 2.8, 1000.0)
    assert lens.image_distance == pytest.approx(52.631578947368425)
    assert lens.magnification == pytest.approx(50.0 / 950.0)
    assert lens.effective_f_number == pytest.approx(2.8 * (1.0 + 50.0 / 950.0))

    cam = Camera(sensor, lens)
    assert isinstance(cam.sensor, Sensor)
    assert isinstance(cam.lens, ThinLens)
    assert isinstance(cam.pinhole, PinholeCamera)
    assert cam.horizontal_fov == pytest.approx(0.16012336, abs=1e-6)
    assert len(cam.footprint) == 2


def test_accessors_taking_arguments_are_methods():
    cam = Camera(machine_vision_sensor(), ThinLens(50.0, 2.8, 1000.0))
    assert cam.coc_px_at(1200.0) == pytest.approx(45.4033634, abs=1e-5)
    assert cam.pixel_size_at(1000.0) == pytest.approx(0.06555)
    near, far = cam.depth_of_field_px(1.0)
    assert near < 1000.0 < far


def test_from_field_width_solves_a_focal_length():
    cam = Camera.from_field_width(machine_vision_sensor(), 2.8, 1000.0, 160.4664)
    assert cam.lens.focal_length == pytest.approx(50.0, abs=1e-9)
    assert cam.footprint[0] == pytest.approx(160.4664, abs=1e-9)


def test_repr_is_informative():
    assert "Sensor(" in repr(machine_vision_sensor())
    assert "ThinLens(" in repr(ThinLens(50.0, 2.8, 500.0))
    assert "Camera(" in repr(small_camera())


# ================================================================================================
# Array shapes, dtypes and orientation
# ================================================================================================


def test_project_returns_one_row_per_input_point():
    cam = small_camera()
    iso = looking_down(500.0)
    points = numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1000.0]])
    px = cam.project(points, iso)

    assert px.shape == (3, 2)
    assert px.dtype == numpy.float64

    # The third point is behind the camera, so its row is NaN and the rows still line up
    assert numpy.all(numpy.isfinite(px[:2]))
    assert numpy.all(numpy.isnan(px[2]))

    # The target of the view lands on the principal point
    assert px[0, 0] == pytest.approx(cam.pinhole.cx, abs=1e-9)
    assert px[0, 1] == pytest.approx(cam.pinhole.cy, abs=1e-9)


def test_view_buffer_arrays_have_the_documented_shape_and_dtype():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1000.0))

    assert view.width == 41
    assert view.height == 41
    for array in (view.depth, view.coc_px, view.incidence):
        assert array.shape == (41, 41)
        assert array.dtype == numpy.float64


def test_the_depth_array_is_not_transposed():
    """
    `DMatrix` is column major and numpy is row major, so handing the buffer over directly would
    transpose the result. A scene tilted about the x axis makes depth vary down the rows and stay
    constant across the columns, which a transpose would swap.
    """
    cam = small_camera()
    iso = look_at(Point3(0.0, 0.0, 500.0), Point3(0.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0))
    tilted = flat_square(1000.0)
    tilted.transform_in_place(Iso3.from_rotation(0.4, 1.0, 0.0, 0.0))

    depth = ViewBuffer(cam, iso, tilted).depth
    assert numpy.all(numpy.isfinite(depth))

    row_spread = float(numpy.ptp(depth[:, 20]))
    col_spread = float(numpy.ptp(depth[20, :]))
    assert row_spread > 10.0 * max(col_spread, 1e-9), (
        f"depth should vary down rows and not across columns, got {row_spread} vs {col_spread}"
    )


def test_the_depth_array_agrees_with_depth_at():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1000.0))
    depth = view.depth

    # depth_at is indexed by column and row, the array by row and column
    for x, y in [(0, 0), (20, 20), (40, 7), (13, 40)]:
        assert depth[y, x] == pytest.approx(view.depth_at(x, y))

    # Out of bounds is NaN rather than an error
    assert math.isnan(view.depth_at(-1, 0))
    assert math.isnan(view.depth_at(41, 0))


def test_misses_are_nan_in_the_arrays():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1.0))
    depth = view.depth
    assert numpy.any(numpy.isnan(depth))
    assert numpy.any(numpy.isfinite(depth))
    assert view.hit_count == int(numpy.sum(numpy.isfinite(depth)))
    assert view.hit_fraction == pytest.approx(view.hit_count / (41 * 41))


def test_cached_arrays_are_returned_repeatedly():
    view = ViewBuffer(small_camera(), looking_down(500.0), flat_square(1000.0))
    first = view.depth
    second = view.depth
    assert numpy.shares_memory(first, second)


def test_project_polyline_returns_a_run_per_visible_stretch():
    cam = small_camera()
    iso = looking_down(500.0)
    points = numpy.array([[-5.0, 0.0, 0.0], [0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    runs = cam.project_polyline(points, iso)

    assert len(runs) == 1
    assert runs[0].shape == (3, 2)
    assert runs[0].dtype == numpy.float64


def test_project_outline_returns_segments_and_codes():
    cam = small_camera()
    mesh = Mesh3.create_box(20.0, 20.0, 20.0, True)
    iso = look_at(Point3(60.0, 80.0, 100.0), Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0))

    segments, codes = cam.project_outline(mesh, iso)
    assert segments.shape[1] == 4
    assert codes.shape == (segments.shape[0],)
    assert codes.dtype == numpy.uint8
    assert set(numpy.unique(codes)).issubset({0, 1})

    with pytest.raises(ValueError):
        cam.project_outline(mesh, iso, max_edge_px=0.0)


def test_frustum_corners_are_eight_points():
    cam = small_camera()
    iso = looking_down(500.0)
    corners = cam.frustum_corners(iso, 100.0, 400.0)
    assert corners.shape == (8, 3)

    # Both faces project onto the same four image corners
    px = cam.project(corners, iso)
    assert px[:4] == pytest.approx(px[4:], abs=1e-6)

    with pytest.raises(ValueError):
        cam.frustum_corners(iso, 400.0, 100.0)


def test_back_project_produces_a_castable_ray_bundle():
    cam = small_camera()
    iso = looking_down(500.0)
    pixels = numpy.array([[20.5, 20.5], [0.5, 0.5]])
    bundle = cam.back_project(pixels, iso)
    assert len(bundle) == 2

    hits = bundle.intersect_mesh(flat_square(1000.0))
    assert hits.shape == (2, 3)
    assert hits[:, 2] == pytest.approx(numpy.zeros(2), abs=1e-9)

    with pytest.raises(ValueError):
        cam.back_project(numpy.zeros((2, 3)), iso)


# ================================================================================================
# Coverage masks and the point cloud
# ================================================================================================


def test_face_masks_have_the_length_of_the_target():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1000.0))
    assert view.target_face_count == 2

    mask = view.face_mask
    assert len(mask) == 2
    assert mask.count_true == 2


def test_face_mask_usable_applies_only_the_criteria_it_is_given():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1000.0))

    assert view.face_mask_usable(allow_backface=True).count_true == 2
    assert view.face_mask_usable(max_coc_px=1.0).count_true == 2
    assert view.face_mask_usable(max_incidence=0.01).count_true == 2

    # Seen from below the square is a back face
    behind = ViewBuffer(cam, looking_down(-500.0), flat_square(1000.0))
    assert behind.face_mask_usable().count_true == 0
    assert behind.face_mask_usable(allow_backface=True).count_true == 2


def test_to_point_cloud_carries_points_and_normals():
    cam = small_camera()
    view = ViewBuffer(cam, looking_down(500.0), flat_square(1000.0))
    cloud = view.to_point_cloud()

    assert cloud.point_count == view.hit_count
    assert cloud.points.shape == (view.hit_count, 3)
    assert cloud.point_normals.shape == (view.hit_count, 3)
    assert cloud.points[:, 2] == pytest.approx(numpy.zeros(view.hit_count), abs=1e-9)


def test_the_buffer_keeps_the_camera_and_pose_it_was_rendered_from():
    cam = small_camera()
    iso = looking_down(500.0)
    view = ViewBuffer(cam, iso, flat_square(1000.0))
    assert view.camera.sensor.width_px == cam.sensor.width_px
    assert view.iso.origin.z == pytest.approx(500.0)


def test_an_obstruction_blocks_without_being_recorded():
    cam = small_camera()
    iso = looking_down(500.0)
    target = flat_square(1000.0)

    # Sized to cover the middle of the frame but not its corners: at z = 250 the frame is only
    # about 0.92mm from the axis to the edge, so a larger blocker would hide everything.
    blocker = flat_square(0.3)
    blocker.transform_in_place(Iso3.from_translation(0.0, 0.0, 250.0))

    clear = ViewBuffer(cam, iso, target)
    blocked = ViewBuffer(cam, iso, target, blocker)
    assert blocked.hit_count < clear.hit_count
    assert math.isnan(blocked.depth_at(20, 20))
    assert blocked.depth_at(0, 0) == pytest.approx(500.0, abs=1e-9)
