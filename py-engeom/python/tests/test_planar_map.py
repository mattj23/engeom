"""
Tests for `PlanarMap` and the frame choices of `Mesh3.section_with_plane`. The Rust implementation
contains the numerical tests; these tests verify that the bindings expose the same API and that
the frame keyword arguments behave as documented in the stub.
"""

from __future__ import annotations

import numpy
import pytest

from engeom.geom2 import Curve2, Line2, Point2, Segment2, SurfacePoint2, Vector2
from engeom.geom3 import (
    Curve3,
    Line3,
    Mesh3,
    PlanarMap,
    Plane3,
    Point3,
    Segment3,
    SurfacePoint3,
    Vector3,
)


def tilted_plane() -> Plane3:
    return Plane3.from_point_normal(1.0, -2.0, 3.0, 1.0, 2.0, 3.0)


def close(a, b, tol=1e-9) -> bool:
    return numpy.allclose(numpy.array(list(a)), numpy.array(list(b)), atol=tol)


# =================================================================================================
# Construction and accessors
# =================================================================================================


def test_xy_plane_map_is_the_identity():
    m = PlanarMap.from_plane(Plane3(0.0, 0.0, 1.0, 0.0))
    assert close(m.point_to_2d(Point3(1.5, -2.5, 7.0)), Point2(1.5, -2.5))
    assert close(m.point_to_3d(Point2(1.5, -2.5)), Point3(1.5, -2.5, 0.0))


def test_frame_inverse_and_plane_agree():
    plane = tilted_plane()
    m = PlanarMap.from_plane(plane)

    p = Point3(4.0, 5.0, 6.0)
    assert close(m.inverse() @ (m.frame @ p), p)
    assert close(m.plane.normal, plane.normal)
    assert m.plane.d == pytest.approx(plane.d)


def test_oriented_map_places_the_origin_and_x_axis():
    plane = Plane3(0.0, 0.0, 1.0, 0.0)
    m = PlanarMap.from_plane_oriented(plane, Point3(1.0, 1.0, 7.0), Vector3(0.0, 1.0, 3.0))

    assert close(m.point_to_2d(Point3(1.0, 1.0, 0.0)), Point2(0.0, 0.0))
    assert close(m.point_to_2d(Point3(1.0, 3.0, 0.0)), Point2(2.0, 0.0))


def test_oriented_map_rejects_x_along_the_normal():
    with pytest.raises(ValueError):
        PlanarMap.from_plane_oriented(Plane3(0.0, 0.0, 1.0, 0.0), Point3(0, 0, 0), Vector3(0, 0, 1))


def test_svd_map_centers_the_points_with_x_along_the_long_side():
    plane = Plane3(0.0, 0.0, 1.0, 0.0)
    # A 2 by 10 rectangle centered on (3, -1), long side along y, first vertex repeated as last.
    pts = numpy.array(
        [[2, -6, 0], [4, -6, 0], [4, 4, 0], [2, 4, 0], [2, -6, 0]], dtype=float
    )
    # Half the length of each adjacent segment, which is what the section uses.
    weights = numpy.array([1.0, 6.0, 6.0, 6.0, 5.0])

    m = PlanarMap.from_plane_svd(plane, pts, weights)
    flat = m.points_to_2d(pts)

    assert close(m.frame.origin, Point3(3.0, -1.0, 0.0))
    assert flat[:, 0].max() - flat[:, 0].min() == pytest.approx(10.0)
    assert flat[:, 1].max() - flat[:, 1].min() == pytest.approx(2.0)


def test_svd_map_rejects_mismatched_weights():
    pts = numpy.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    with pytest.raises(ValueError):
        PlanarMap.from_plane_svd(Plane3(0.0, 0.0, 1.0, 0.0), pts, numpy.array([1.0]))


# =================================================================================================
# Conversions
# =================================================================================================


def test_points_round_trip_in_bulk():
    m = PlanarMap.from_plane(tilted_plane())
    flat = numpy.array([[0.0, 0.0], [1.0, 2.0], [-3.0, 0.5]])
    lifted = m.points_to_3d(flat)

    assert lifted.shape == (3, 3)
    assert numpy.allclose(m.points_to_2d(lifted), flat, atol=1e-9)


def test_vectors_round_trip_and_the_normal_collapses():
    plane = tilted_plane()
    m = PlanarMap.from_plane(plane)

    v3 = m.vector_to_3d(Vector2(0.3, -0.7))
    assert close(m.vector_to_2d(v3), Vector2(0.3, -0.7))

    u3 = m.unit_vector_to_3d(Vector2(2.0, 1.0))
    assert u3.norm() == pytest.approx(1.0)
    assert close(m.unit_vector_to_2d(u3), Vector2(2.0, 1.0).normalized())

    with pytest.raises(ValueError):
        m.unit_vector_to_2d(plane.normal)


def test_surface_points_lines_and_segments_round_trip():
    m = PlanarMap.from_plane(tilted_plane())

    sp2 = SurfacePoint2(1.0, 2.0, -1.0, 1.0)
    back = m.surface_point_to_2d(m.surface_point_to_3d(sp2))
    assert close(back.point, sp2.point)
    assert close(back.normal, sp2.normal)

    l2 = Line2(1.0, 1.0, 2.0, -1.0)
    back = m.line_to_2d(m.line_to_3d(l2))
    assert close(back.origin, l2.origin)
    assert close(back.direction, l2.direction)

    s2 = Segment2(0.0, 0.0, 3.0, 4.0)
    s3 = m.segment_to_3d(s2)
    assert s3.length == pytest.approx(5.0)
    back = m.segment_to_2d(s3)
    assert close(back.a, s2.a)
    assert close(back.b, s2.b)


def test_directions_along_the_normal_raise():
    plane = tilted_plane()
    m = PlanarMap.from_plane(plane)
    n = plane.normal
    o = m.point_to_3d(Point2(0.0, 0.0))

    with pytest.raises(ValueError):
        m.surface_point_to_2d(SurfacePoint3(o.x, o.y, o.z, n.x, n.y, n.z))
    with pytest.raises(ValueError):
        m.line_to_2d(Line3(o.x, o.y, o.z, n.x, n.y, n.z))
    with pytest.raises(ValueError):
        m.segment_to_2d(Segment3(o.x, o.y, o.z, o.x + n.x, o.y + n.y, o.z + n.z))


def test_curves_round_trip_and_keep_closure():
    m = PlanarMap.from_plane(tilted_plane())
    square = Curve2(numpy.array([[0, 0], [2, 0], [2, 1], [0, 1]], dtype=float), tol=1e-9, force_closed=True)
    assert square.is_closed

    lifted = m.curve_to_3d(square)
    assert isinstance(lifted, Curve3)
    assert lifted.length == pytest.approx(square.length)

    back = m.curve_to_2d(lifted)
    assert back.is_closed
    assert numpy.allclose(back.points, square.points, atol=1e-9)


# =================================================================================================
# Frame arguments of Mesh3.section_with_plane
# =================================================================================================


def test_svd_frame_centers_an_off_center_section():
    mesh = Mesh3.create_box(4.0, 2.0, 6.0)
    # A plane whose auto frame origin is well away from the section.
    plane = Plane3.from_point_normal(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)
    mesh = mesh.transform_copy(engeom_translation(10.0, -5.0, 0.0))

    flat = mesh.section_with_plane(plane, tol=1e-9, frame="svd").to_2d()
    aabb = flat[0].aabb
    assert close(aabb.center, Point2(0.0, 0.0))
    assert aabb.max.x - aabb.min.x == pytest.approx(4.0)
    assert aabb.max.y - aabb.min.y == pytest.approx(2.0)


def test_oriented_frame_places_the_origin_and_x_axis():
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    plane = Plane3(0.0, 0.0, 1.0, 0.0)

    flat = mesh.section_with_plane(
        plane, tol=1e-9, origin=Point3(-1.0, -2.0, 5.0), x=Vector3(0.0, 1.0, 0.5)
    ).to_2d()
    aabb = flat[0].aabb
    assert aabb.min.x == pytest.approx(0.0)
    assert aabb.max.x == pytest.approx(4.0)
    assert aabb.min.y == pytest.approx(-2.0)
    assert aabb.max.y == pytest.approx(0.0)


def test_frame_argument_validation():
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    plane = Plane3(0.0, 0.0, 1.0, 0.0)

    with pytest.raises(ValueError):
        mesh.section_with_plane(plane, origin=Point3(0, 0, 0))
    with pytest.raises(ValueError):
        mesh.section_with_plane(plane, frame="svd", origin=Point3(0, 0, 0), x=Vector3(1, 0, 0))
    with pytest.raises(ValueError):
        mesh.section_with_plane(plane, frame="principal")
    with pytest.raises(ValueError):
        mesh.section_with_plane(plane, origin=Point3(0, 0, 0), x=Vector3(0, 0, 1))


def test_a_2d_result_lifts_back_onto_the_section():
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    plane = Plane3.from_point_normal(0.1, 0.2, 0.3, 0.3, -0.5, 0.8)

    section = mesh.section_with_plane(plane, tol=1e-9)
    flat = section.to_2d()
    back = section.map.curve_group_to_3d(flat)

    assert back.length == pytest.approx(section.curves.length)
    for p in back[0].points:
        assert plane.distance_to_point(Point3(*p)) == pytest.approx(0.0, abs=1e-9)


def engeom_translation(x: float, y: float, z: float):
    from engeom.geom3 import Iso3

    return Iso3.from_translation(x, y, z)
