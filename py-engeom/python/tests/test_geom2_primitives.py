"""
    Tests of geometric primitives (vectors, points, etc.) in geom2 module.
"""
import math

import pytest
import numpy
from engeom.geom2 import Vector2, Point2, SurfacePoint2, Segment2, Arc2, Circle2, Iso2


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


def test_iso2_from_translation():
    iso = Iso2.from_translation(1.0, 2.0)
    assert iso.origin.x == pytest.approx(1.0)
    assert iso.origin.y == pytest.approx(2.0)
    assert iso == Iso2(1.0, 2.0, 0.0)


def test_iso2_from_rotation():
    iso = Iso2.from_rotation(math.pi / 2)
    rotated = iso @ Point2(1, 0)
    assert rotated.x == pytest.approx(0.0)
    assert rotated.y == pytest.approx(1.0)
    assert iso == Iso2(0.0, 0.0, math.pi / 2)


# The translation and rotation components should recompose into the original isometry, in the
# same order the Iso2 constructor applies them (rotate first, then translate).
def test_iso2_translation_rotation_decomposition():
    iso = Iso2(3.0, -4.0, 0.7)

    assert iso.translation() == Iso2.from_translation(3.0, -4.0)
    assert iso.rotation() == Iso2.from_rotation(0.7)
    assert iso.translation() @ iso.rotation() == iso


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


def test_arc2_from_consensus_bounds_inlier_sector():
    rng = numpy.random.default_rng(0)
    cx, cy, r = 2.0, -1.0, 1.3

    # Inliers over the counter-clockwise sector [0, pi] (an upper half), small radial noise.
    n = 120
    t = numpy.linspace(0.0, math.pi, n)
    rr = r + rng.normal(0.0, 0.004, n)
    inliers = numpy.column_stack([cx + rr * numpy.cos(t), cy + rr * numpy.sin(t)])

    # A cluster of gross outliers near the empty lower sector.
    m = 40
    ot = rng.uniform(0.0, 2.0 * math.pi, m)
    outliers = numpy.column_stack([cx + 0.4 * numpy.cos(ot), cy - 4.0 + 0.4 * numpy.sin(ot)])

    points = numpy.vstack([inliers, outliers])
    arc = Arc2.from_consensus(points, 0.02, seed=42)

    assert arc.center.x == pytest.approx(cx, abs=5e-3)
    assert arc.center.y == pytest.approx(cy, abs=5e-3)
    assert arc.r == pytest.approx(r, abs=5e-3)

    # The arc spans only the inlier sector [0, pi], not the outlier-adjacent lower half. The start
    # angle is normalized to (-pi, pi] because a noisy endpoint near t=0 can dip just below zero and
    # wrap to just under 2*pi.
    assert arc.angle > 0.0
    start = ((arc.angle0 + math.pi) % (2.0 * math.pi)) - math.pi
    assert start == pytest.approx(0.0, abs=2e-2)
    assert arc.angle == pytest.approx(math.pi, abs=3e-2)


def test_circle2_transformed_by():
    c = Circle2(1.0, 0.0, 2.0)
    iso = Iso2(0.0, 3.0, math.pi / 2)
    moved = c.transformed_by(iso)
    expected = iso @ c.center
    assert moved.center.x == pytest.approx(expected.x)
    assert moved.center.y == pytest.approx(expected.y)
    assert moved.r == pytest.approx(2.0)

    # The matmul operator is an alias for transformed_by
    by_op = iso @ c
    assert by_op.center.x == pytest.approx(moved.center.x)
    assert by_op.center.y == pytest.approx(moved.center.y)
    assert by_op.r == pytest.approx(moved.r)


def test_circle2_resized_by():
    c = Circle2(1.0, 2.0, 3.0)
    bigger = c.resized_by(0.5)
    smaller = c.resized_by(-0.5)
    assert bigger.center.x == pytest.approx(1.0)
    assert bigger.center.y == pytest.approx(2.0)
    assert bigger.r == pytest.approx(3.5)
    assert smaller.r == pytest.approx(2.5)
    assert c.r == pytest.approx(3.0)


def test_circle2_contained_points():
    c = Circle2(1.0, 1.0, 1.0)
    points = numpy.array(
        [
            [1.0, 1.0],  # center, inside
            [5.0, 5.0],  # far outside
            [2.0, 1.0],  # exactly on the boundary, counts as inside
            [1.0, -1.0],  # outside
            [1.5, 1.0],  # inside
        ]
    )

    result = c.contained_points(points)
    expected = numpy.array([[1.0, 1.0], [2.0, 1.0], [1.5, 1.0]])
    assert result.shape == (3, 2)
    assert numpy.allclose(result, expected)


def test_circle2_contained_points_none_inside():
    c = Circle2(0.0, 0.0, 1.0)
    points = numpy.array([[5.0, 5.0], [-3.0, 0.0]])
    result = c.contained_points(points)
    assert result.shape == (0, 2)


def test_circle2_from_min_enclosing():
    points = numpy.array([[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0], [1.0, 1.5]])
    circle = Circle2.from_min_enclosing(points)
    assert abs(circle.x - 1.0) < 1e-12
    assert abs(circle.y - 1.0) < 1e-12
    assert abs(circle.r - math.sqrt(2.0)) < 1e-12


def test_circle2_from_min_enclosing_empty_raises():
    points = numpy.empty((0, 2))
    with pytest.raises(ValueError):
        Circle2.from_min_enclosing(points)
