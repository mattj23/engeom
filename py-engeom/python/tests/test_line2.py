"""Tests for Line2."""
import math
import pickle

import pytest
from engeom.geom2 import Line2, Point2, Vector2, Iso2


EPS = 1e-10


# ── helpers ──────────────────────────────────────────────────────────────────

def approx_eq(a: float, b: float) -> bool:
    return abs(a - b) < EPS


def point_eq(p: Point2, x: float, y: float) -> bool:
    return approx_eq(p.x, x) and approx_eq(p.y, y)


def vec_eq(v: Vector2, x: float, y: float) -> bool:
    return approx_eq(v.x, x) and approx_eq(v.y, y)


# ── construction ─────────────────────────────────────────────────────────────

def test_new_stores_origin_and_direction():
    line = Line2(1.0, 2.0, 3.0, 4.0)
    assert point_eq(line.origin, 1.0, 2.0)
    assert vec_eq(line.direction, 3.0, 4.0)


def test_x_axis():
    line = Line2.x_axis()
    assert point_eq(line.origin, 0.0, 0.0)
    assert vec_eq(line.direction, 1.0, 0.0)


def test_y_axis():
    line = Line2.y_axis()
    assert point_eq(line.origin, 0.0, 0.0)
    assert vec_eq(line.direction, 0.0, 1.0)


def test_from_points_direction_is_difference():
    line = Line2.from_points(Point2(1.0, 0.0), Point2(4.0, 0.0))
    assert vec_eq(line.direction, 3.0, 0.0)
    assert point_eq(line.origin, 1.0, 0.0)


def test_new_normalize_unit_direction():
    line = Line2.new_normalize(0.0, 0.0, 3.0, 0.0)
    assert approx_eq(line.direction.norm(), 1.0)
    assert vec_eq(line.direction, 1.0, 0.0)


def test_normalized_returns_unit_direction():
    line = Line2(0.0, 0.0, 3.0, 4.0).normalized()
    assert approx_eq(line.direction.norm(), 1.0)


# ── normal ───────────────────────────────────────────────────────────────────

def test_normal_is_clockwise_rotation_of_direction():
    # X-axis direction → normal should be (0, -1) (CW 90°)
    line = Line2.new_normalize(0.0, 0.0, 1.0, 0.0)
    assert approx_eq(line.normal.norm(), 1.0)
    assert vec_eq(line.normal, 0.0, -1.0)


def test_normal_y_axis():
    # Y-axis direction → normal should be (1, 0) (CW 90°)
    line = Line2.new_normalize(0.0, 0.0, 0.0, 1.0)
    assert vec_eq(line.normal, 1.0, 0.0)


# ── at ───────────────────────────────────────────────────────────────────────

def test_at_zero_returns_origin():
    line = Line2(1.0, 2.0, 1.0, 0.0)
    assert point_eq(line.at(0.0), 1.0, 2.0)


def test_at_one_returns_origin_plus_direction():
    line = Line2(1.0, 2.0, 0.0, 1.0)
    assert point_eq(line.at(1.0), 1.0, 3.0)
    assert point_eq(line.at(-1.0), 1.0, 1.0)


# ── scalar_project ────────────────────────────────────────────────────────────

def test_scalar_project_on_line():
    line = Line2.x_axis()
    assert approx_eq(line.scalar_project(Point2(5.0, 0.0)), 5.0)


def test_scalar_project_perpendicular_offset_unchanged():
    line = Line2.x_axis()
    assert approx_eq(line.scalar_project(Point2(3.0, 7.0)), 3.0)


def test_scalar_project_normalized_equals_arc_length():
    line = Line2.new_normalize(0.0, 0.0, 1.0, 0.0)
    assert approx_eq(line.scalar_project(Point2(4.0, 0.0)), 4.0)


# ── closest_point ─────────────────────────────────────────────────────────────

def test_closest_point_drops_perpendicular():
    line = Line2.x_axis()
    cp = line.closest_point(Point2(4.0, 3.0))
    assert point_eq(cp, 4.0, 0.0)


def test_closest_point_on_line_is_itself():
    line = Line2.x_axis()
    p = Point2(2.0, 0.0)
    assert point_eq(line.closest_point(p), 2.0, 0.0)


# ── distance_to ───────────────────────────────────────────────────────────────

def test_distance_to_known_value():
    line = Line2.x_axis()
    assert approx_eq(line.distance_to(Point2(0.0, 3.0)), 3.0)


def test_distance_to_is_always_non_negative():
    line = Line2.x_axis()
    assert line.distance_to(Point2(0.0, -5.0)) >= 0.0


# ── signed_distance_to ────────────────────────────────────────────────────────

def test_signed_distance_right_of_travel_positive():
    # X-axis, point below (CW normal → positive to the right / below)
    line = Line2.x_axis()
    assert approx_eq(line.signed_distance_to(Point2(0.0, -3.0)), 3.0)


def test_signed_distance_left_of_travel_negative():
    line = Line2.x_axis()
    assert approx_eq(line.signed_distance_to(Point2(0.0, 3.0)), -3.0)


def test_signed_distance_on_line_is_zero():
    line = Line2.x_axis()
    assert approx_eq(line.signed_distance_to(Point2(5.0, 0.0)), 0.0)


# ── intersect ─────────────────────────────────────────────────────────────────

def test_intersect_perpendicular_lines():
    a = Line2(0.0, 0.0, 1.0, 0.0)
    b = Line2(3.0, 0.0, 0.0, 1.0)
    pt = a.intersect(b)
    assert pt is not None
    assert point_eq(pt, 3.0, 0.0)


def test_intersect_diagonal_lines():
    a = Line2(0.0, 0.0, 1.0, 1.0)
    b = Line2(1.0, 0.0, -1.0, 1.0)
    pt = a.intersect(b)
    assert pt is not None
    assert point_eq(pt, 0.5, 0.5)


def test_intersect_parallel_returns_none():
    a = Line2(0.0, 0.0, 1.0, 0.0)
    b = Line2(0.0, 1.0, 1.0, 0.0)
    assert a.intersect(b) is None


def test_intersect_coincident_returns_none():
    a = Line2(0.0, 0.0, 1.0, 0.0)
    b = Line2(1.0, 0.0, 2.0, 0.0)
    assert a.intersect(b) is None


# ── new_parallel ──────────────────────────────────────────────────────────────

def test_new_parallel_shifts_origin_along_normal():
    # X-axis, CW normal = (0, -1); delta_n=1 shifts origin by (0, -1)
    line = Line2.new_normalize(0.0, 0.0, 1.0, 0.0)
    shifted = line.new_parallel(1.0)
    assert point_eq(shifted.origin, 0.0, -1.0)
    assert vec_eq(shifted.direction, 1.0, 0.0)


def test_new_parallel_preserves_direction():
    line = Line2.new_normalize(2.0, 3.0, 0.0, 1.0)
    shifted = line.new_parallel(2.0)
    assert vec_eq(shifted.direction, 0.0, 1.0)


# ── new_shifted_along ─────────────────────────────────────────────────────────

def test_new_shifted_along_moves_origin():
    line = Line2(0.0, 0.0, 1.0, 0.0)
    shifted = line.new_shifted_along(3.0)
    assert point_eq(shifted.origin, 3.0, 0.0)
    assert vec_eq(shifted.direction, 1.0, 0.0)


def test_new_shifted_along_negative():
    line = Line2(5.0, 0.0, 1.0, 0.0)
    shifted = line.new_shifted_along(-2.0)
    assert point_eq(shifted.origin, 3.0, 0.0)


# ── new_transformed_by ───────────────────────────────────────────────────────

def test_new_transformed_by_translation():
    line = Line2.x_axis()
    iso = Iso2(0.0, 5.0, 0.0)
    t = line.new_transformed_by(iso)
    assert point_eq(t.origin, 0.0, 5.0)
    assert vec_eq(t.direction, 1.0, 0.0)


def test_new_transformed_by_rotation():
    line = Line2.x_axis()
    iso = Iso2(0.0, 0.0, math.pi / 2)
    t = line.new_transformed_by(iso)
    assert approx_eq(t.direction.x, 0.0)
    assert approx_eq(t.direction.y, 1.0)


def test_new_transformed_preserves_points_on_line():
    line = Line2(1.0, 0.0, 1.0, 0.0)
    iso = Iso2(3.0, -1.0, math.pi / 4)
    tline = line.new_transformed_by(iso)
    for t in [-1.0, 0.0, 1.0, 2.5]:
        expected = iso @ line.at(t)
        actual = tline.at(t)
        assert point_eq(actual, expected.x, expected.y)


# ── to_iso_from_x / to_iso_from_y ────────────────────────────────────────────

def test_to_iso_from_x_inverse_maps_point_to_x_axis():
    line = Line2.new_normalize(2.0, 3.0, 0.0, 1.0)
    iso = line.to_iso_from_x()
    # A point on the line at t=1 should map to (1, 0) in local space
    world_pt = line.at(1.0)
    local_pt = iso.inverse() @ world_pt
    assert approx_eq(local_pt.x, 1.0)
    assert approx_eq(local_pt.y, 0.0)


def test_to_iso_from_y_inverse_maps_point_to_y_axis():
    line = Line2.new_normalize(2.0, 3.0, 1.0, 0.0)
    iso = line.to_iso_from_y()
    world_pt = line.at(1.0)
    local_pt = iso.inverse() @ world_pt
    assert approx_eq(local_pt.x, 0.0)
    assert approx_eq(local_pt.y, 1.0)


# ── pickle ────────────────────────────────────────────────────────────────────

def test_pickle_roundtrip():
    line = Line2(1.0, 2.0, 3.0, 4.0)
    assert line == pickle.loads(pickle.dumps(line))


def test_pickle_normalized_roundtrip():
    line = Line2.new_normalize(1.0, 2.0, 1.0, 0.0)
    restored = pickle.loads(pickle.dumps(line))
    assert approx_eq(restored.origin.x, line.origin.x)
    assert approx_eq(restored.direction.x, line.direction.x)


# ── __eq__ ────────────────────────────────────────────────────────────────────

def test_eq_same_values():
    assert Line2(1.0, 2.0, 3.0, 4.0) == Line2(1.0, 2.0, 3.0, 4.0)


def test_eq_different_direction():
    assert Line2(1.0, 2.0, 3.0, 4.0) != Line2(1.0, 2.0, 3.0, 5.0)


def test_eq_different_origin():
    assert Line2(1.0, 2.0, 3.0, 4.0) != Line2(1.0, 3.0, 3.0, 4.0)
