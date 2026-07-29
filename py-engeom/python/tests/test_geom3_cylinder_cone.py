"""
    Tests of the Cylinder3 and Cone3 geometric primitives in the geom3 module.
"""
import math
import pytest
from engeom.geom3 import (
    Vector3,
    Point3,
    SurfacePoint3,
    Line3,
    Circle3,
    Cylinder3,
    Cone3,
    Iso3,
)


# ==============================================================================
# Cylinder3 tests
# ==============================================================================

def test_cylinder3_construction():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    assert isinstance(cyl.center, Point3)
    assert cyl.center.x == pytest.approx(0.0, abs=1e-6)
    assert cyl.center.y == pytest.approx(0.0, abs=1e-6)
    assert cyl.center.z == pytest.approx(0.0, abs=1e-6)
    assert cyl.r == pytest.approx(2.0)
    assert cyl.length == pytest.approx(10.0)


def test_cylinder3_direction_normalized():
    # Supply un-normalized direction (0, 0, 5); should normalize to (0, 0, 1)
    cyl = Cylinder3(0, 0, 0, 0, 0, 5, 2.0, 10.0)
    assert isinstance(cyl.direction, Vector3)
    assert cyl.direction.z == pytest.approx(1.0)


def test_cylinder3_zero_direction_raises():
    with pytest.raises(ValueError):
        Cylinder3(0, 0, 0, 0, 0, 0, 2.0, 10.0)


def test_cylinder3_from_points():
    cyl = Cylinder3.from_points(Point3(0, 0, 0), Point3(0, 0, 4), 1.5)
    assert cyl.r == pytest.approx(1.5)
    assert cyl.length == pytest.approx(4.0)
    assert cyl.a().z == pytest.approx(0.0, abs=1e-6)
    assert cyl.b().z == pytest.approx(4.0)


def test_cylinder3_from_points_coincident_raises():
    with pytest.raises(ValueError):
        Cylinder3.from_points(Point3(1, 1, 1), Point3(1, 1, 1), 1.0)


def test_cylinder3_a_b_endpoints():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    a = cyl.a()
    b = cyl.b()
    assert isinstance(a, Point3)
    assert isinstance(b, Point3)
    assert a.z == pytest.approx(0.0, abs=1e-6)
    assert b.z == pytest.approx(10.0)


def test_cylinder3_axis_returns_line3():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    axis = cyl.axis()
    assert isinstance(axis, Line3)
    assert axis.origin.z == pytest.approx(0.0, abs=1e-6)
    assert axis.direction.z == pytest.approx(1.0)


def test_cylinder3_start_end_caps():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    start = cyl.start_cap()
    end = cyl.end_cap()
    assert isinstance(start, Circle3)
    assert isinstance(end, Circle3)
    assert start.r == pytest.approx(2.0)
    assert end.r == pytest.approx(2.0)
    assert start.normal.z == pytest.approx(-1.0)
    assert end.normal.z == pytest.approx(1.0)


def test_cylinder3_volume_and_lateral_area():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    assert cyl.volume() == pytest.approx(math.pi * 4.0 * 10.0)
    assert cyl.lateral_area() == pytest.approx(2.0 * math.pi * 2.0 * 10.0)


def test_cylinder3_contains_point():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    assert cyl.contains_point(Point3(0, 0, 5))
    assert not cyl.contains_point(Point3(5, 0, 5))
    assert not cyl.contains_point(Point3(0, 0, -1))


def test_cylinder3_transformed_by():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    iso = Iso3.from_translation(1.0, 2.0, 3.0)
    moved = cyl.transformed_by(iso)
    assert isinstance(moved, Cylinder3)
    assert moved.center.x == pytest.approx(1.0)
    assert moved.center.y == pytest.approx(2.0)
    assert moved.center.z == pytest.approx(3.0)
    assert moved.r == pytest.approx(cyl.r)
    assert moved.length == pytest.approx(cyl.length)


def test_cylinder3_reversed():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    reversed_cyl = cyl.reversed()
    assert isinstance(reversed_cyl, Cylinder3)
    assert reversed_cyl.direction.z == pytest.approx(-1.0)
    assert reversed_cyl.a().z == pytest.approx(cyl.b().z)
    assert reversed_cyl.b().z == pytest.approx(cyl.a().z)


def test_cylinder3_closest_point_clamped():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    sp = cyl.closest_point(Point3(5, 0, 5), False)
    assert isinstance(sp, SurfacePoint3)
    assert sp.point.x == pytest.approx(2.0)
    assert sp.point.z == pytest.approx(5.0)
    assert sp.normal.x == pytest.approx(1.0)


def test_cylinder3_closest_point_infinite_extends_past_ends():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    clamped = cyl.closest_point(Point3(3, 0, 15), False)
    infinite = cyl.closest_point(Point3(3, 0, 15), True)
    assert clamped.point.z == pytest.approx(10.0)
    assert infinite.point.z == pytest.approx(15.0)


def test_cylinder3_closest_point_on_axis_returns_none():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    assert cyl.closest_point(Point3(0, 0, 5), False) is None


def test_cylinder3_pickle_roundtrip():
    import pickle
    cyl = Cylinder3(1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 4.0, 5.0)
    unpickled = pickle.loads(pickle.dumps(cyl))
    assert cyl == unpickled


def test_cylinder3_repr():
    cyl = Cylinder3(0, 0, 0, 0, 0, 1, 2.0, 10.0)
    assert "Cylinder3" in repr(cyl)


# ==============================================================================
# Cone3 tests
# ==============================================================================

def test_cone3_construction():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    assert isinstance(cone.tip, Point3)
    assert cone.tip.x == pytest.approx(0.0, abs=1e-6)
    assert cone.height == pytest.approx(10.0)
    assert cone.r == pytest.approx(2.0)


def test_cone3_direction_normalized():
    cone = Cone3(0, 0, 0, 0, 0, 5, 10.0, 2.0)
    assert isinstance(cone.direction, Vector3)
    assert cone.direction.z == pytest.approx(1.0)


def test_cone3_zero_direction_raises():
    with pytest.raises(ValueError):
        Cone3(0, 0, 0, 0, 0, 0, 10.0, 2.0)


def test_cone3_from_points():
    cone = Cone3.from_points(Point3(0, 0, 0), Point3(0, 0, 4), 1.5)
    assert cone.r == pytest.approx(1.5)
    assert cone.height == pytest.approx(4.0)
    assert cone.base_center().z == pytest.approx(4.0)


def test_cone3_from_points_coincident_raises():
    with pytest.raises(ValueError):
        Cone3.from_points(Point3(1, 1, 1), Point3(1, 1, 1), 1.0)


def test_cone3_base_center():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    base_center = cone.base_center()
    assert isinstance(base_center, Point3)
    assert base_center.z == pytest.approx(10.0)


def test_cone3_axis_returns_line3():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    axis = cone.axis()
    assert isinstance(axis, Line3)
    assert axis.direction.z == pytest.approx(1.0)


def test_cone3_base_returns_circle3():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    base = cone.base()
    assert isinstance(base, Circle3)
    assert base.r == pytest.approx(2.0)
    assert base.normal.z == pytest.approx(1.0)


def test_cone3_half_angle():
    # height 1, radius 1 -> 45 degree half angle
    cone = Cone3(0, 0, 0, 0, 0, 1, 1.0, 1.0)
    assert cone.half_angle() == pytest.approx(math.pi / 4.0)


def test_cone3_slant_height():
    # height 3, radius 4 -> slant 5
    cone = Cone3(0, 0, 0, 0, 0, 1, 3.0, 4.0)
    assert cone.slant_height() == pytest.approx(5.0)


def test_cone3_volume_and_lateral_area():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    assert cone.volume() == pytest.approx(math.pi * 4.0 * 10.0 / 3.0)
    assert cone.lateral_area() == pytest.approx(math.pi * 2.0 * cone.slant_height())


def test_cone3_contains_point():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    assert cone.contains_point(Point3(0, 0, 5))
    assert not cone.contains_point(Point3(2, 0, 5))
    assert not cone.contains_point(Point3(0, 0, -1))


def test_cone3_transformed_by():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    iso = Iso3.from_translation(1.0, 2.0, 3.0)
    moved = cone.transformed_by(iso)
    assert isinstance(moved, Cone3)
    assert moved.tip.x == pytest.approx(1.0)
    assert moved.tip.y == pytest.approx(2.0)
    assert moved.tip.z == pytest.approx(3.0)
    assert moved.height == pytest.approx(cone.height)
    assert moved.r == pytest.approx(cone.r)


def test_cone3_closest_point_clamped():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    sp = cone.closest_point(Point3(0.5, 0, 20), False)
    assert isinstance(sp, SurfacePoint3)
    # Far beyond the base, within the cap radius: clamps to the rim/base plane.
    assert sp.point.z == pytest.approx(10.0)


def test_cone3_closest_point_infinite_extends_past_base():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    # A point exactly on the infinite extension of the slant line beyond the base.
    test_point = Point3(4.0, 0.0, 20.0)
    clamped = cone.closest_point(test_point, False)
    infinite = cone.closest_point(test_point, True)
    assert clamped.point.x == pytest.approx(2.0)
    assert clamped.point.z == pytest.approx(10.0)
    assert infinite.point.x == pytest.approx(4.0)
    assert infinite.point.z == pytest.approx(20.0)


def test_cone3_closest_point_on_axis_returns_none():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    assert cone.closest_point(Point3(0, 0, 5), False) is None


def test_cone3_pickle_roundtrip():
    import pickle
    cone = Cone3(1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 5.0, 4.0)
    unpickled = pickle.loads(pickle.dumps(cone))
    assert cone == unpickled


def test_cone3_repr():
    cone = Cone3(0, 0, 0, 0, 0, 1, 10.0, 2.0)
    assert "Cone3" in repr(cone)
