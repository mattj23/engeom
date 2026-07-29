"""
    Tests of geometric primitives (vectors, points, etc.) in geom2 module.
"""
import math

import pytest
import numpy
from engeom.geom2 import Vector2, Point2, SurfacePoint2, Segment2, Arc2, Iso2


def test_vector_mul_scalar():
    v = Vector2(1, 2)
    result = v * 3
    assert abs(result.x - 3) < 1e-6
    assert abs(result.y - 6) < 1e-6


def test_vector_div_scalar():
    v = Vector2(3, 6)
    result = v / 3
    assert abs(result.x - 1) < 1e-6
    assert abs(result.y - 2) < 1e-6


def test_point_mul_scalar():
    p = Point2(1, 2)
    result = p * 3
    assert abs(result.x - 3) < 1e-6
    assert abs(result.y - 6) < 1e-6


def test_point_div_scalar():
    p = Point2(3, 6)
    result = p / 3
    assert abs(result.x - 1) < 1e-6
    assert abs(result.y - 2) < 1e-6


def test_sp_mul_scalar_pos():
    sp = SurfacePoint2(1, 2, 1, 0)
    result = sp * 3
    assert abs(result.point.x - 3) < 1e-6
    assert abs(result.point.y - 6) < 1e-6
    assert abs(result.normal.x - 1) < 1e-6
    assert abs(result.normal.y) < 1e-6


def test_sp_mul_scalar_neg():
    sp = SurfacePoint2(1, 2, 1, 0)
    result = sp * -3
    assert abs(result.point.x + 3) < 1e-6
    assert abs(result.point.y + 6) < 1e-6
    assert abs(result.normal.x + 1) < 1e-6
    assert abs(result.normal.y) < 1e-6


def test_sp_div_scalar():
    sp = SurfacePoint2(3, 6, 1, 0)
    result = sp / 3
    assert abs(result.point.x - 1) < 1e-6
    assert abs(result.point.y - 2) < 1e-6
    assert abs(result.normal.x - 1) < 1e-6
    assert abs(result.normal.y) < 1e-6


def test_sp_rotate():
    from math import pi
    a = SurfacePoint2(0, 0, 0, 1)
    b = a.normal_rotated(pi / 2)

    assert abs(b.point.x) < 1e-6
    assert abs(b.point.y) < 1e-6
    assert abs(b.normal.x + 1) < 1e-6
    assert abs(b.normal.y) < 1e-6


# Test that a vector plus a vector is a vector.
def test_vector_plus_vector():
    v1 = Vector2(1, 2)
    v2 = Vector2(3, 4)
    v3 = v1 + v2
    assert isinstance(v3, Vector2)


def test_vector_plus_point():
    v = Vector2(1, 2)
    p = Point2(3, 4)
    result = v + p
    assert isinstance(result, Point2)


# Test that a point plus a vector is a point.
def test_point_plus_vector():
    p = Point2(1, 2)
    v = Vector2(3, 4)
    result = p + v
    assert isinstance(result, Point2)


# Test that a point minus a point is a vector.
def test_point_minus_point():
    p1 = Point2(1, 2)
    p2 = Point2(3, 4)
    result = p1 - p2
    assert isinstance(result, Vector2)


# Test that a point minus a vector is a point.
def test_point_minus_vector():
    p = Point2(3, 4)
    v = Vector2(1, 2)
    result = p - v
    assert isinstance(result, Point2)


# Test that a vector minus a vector is a vector.
def test_vector_minus_vector():
    v1 = Vector2(3, 4)
    v2 = Vector2(1, 2)
    result = v1 - v2
    assert isinstance(result, Vector2)


# Test that an Iso2 matmul by a vector returns a vector.
def test_iso2_matmul_vector():
    iso = Iso2(1, 2, 0.5)
    v = Vector2(3, 4)
    result = iso @ v
    assert isinstance(result, Vector2)


# Test that an Iso2 matmul by a point returns a point.
def test_iso2_matmul_point():
    iso = Iso2(1, 2, 0.5)
    p = Point2(3, 4)
    result = iso @ p
    assert isinstance(result, Point2)


# Test that an Iso2 matmul by a surface point returns a surface point.
def test_iso2_matmul_surfacepoint():
    iso = Iso2(1, 2, 0.5)
    sp = SurfacePoint2(3, 4, 1, 0)
    result = iso @ sp
    assert isinstance(result, SurfacePoint2)


# Test that an Iso2 matmul by another Iso2 returns an Iso2.
def test_iso2_matmul_iso2():
    iso1 = Iso2(1, 2, 0.5)
    iso2 = Iso2(3, 4, 0.5)
    result = iso1 @ iso2
    assert isinstance(result, Iso2)


def test_segment2_transformed_by():
    s = Segment2(0, 0, 1, 0)
    iso = Iso2(1, 2, 0)
    moved = s.transformed_by(iso)
    assert moved.a.x == pytest.approx(1.0)
    assert moved.a.y == pytest.approx(2.0)
    assert moved.b.x == pytest.approx(2.0)
    assert moved.b.y == pytest.approx(2.0)


# Test that an Iso2 matmul by a Segment2 returns a Segment2.
def test_iso2_matmul_segment2():
    iso = Iso2(1, 2, 0.5)
    s = Segment2(0, 0, 1, 0)
    result = iso @ s
    assert isinstance(result, Segment2)


def test_iso2_matmul_segment2_matches_transformed_by():
    s = Segment2(0, 0, 1, 2)
    iso = Iso2(1, 2, 0.5)
    assert (iso @ s) == s.transformed_by(iso)


def test_arc2_construction():
    a = Arc2(0, 0, 1, 0, math.pi / 2)
    assert a.center.x == pytest.approx(0.0)
    assert a.center.y == pytest.approx(0.0)
    assert a.r == pytest.approx(1.0)
    assert a.angle0 == pytest.approx(0.0)
    assert a.angle == pytest.approx(math.pi / 2)


def test_arc2_to_points_includes_endpoints():
    a = Arc2(0, 0, 1, 0, math.pi / 2)
    points = a.to_points(0.01)
    assert points.shape[1] == 2
    assert points[0, 0] == pytest.approx(a.a.x)
    assert points[0, 1] == pytest.approx(a.a.y)
    assert points[-1, 0] == pytest.approx(a.b.x)
    assert points[-1, 1] == pytest.approx(a.b.y)
