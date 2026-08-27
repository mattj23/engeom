"""
Tests for the `CurveGroup2` and `CurveGroup3` bindings.

A curve group is one rigid body made of several disjoint polylines, which is what a planar section
of a mesh naturally is. These check the sequence behavior the bindings promise, the projection from
3D into a plane, and the file round trip.
"""

from __future__ import annotations

import numpy
import pytest

from engeom.geom2 import Curve2, CurveGroup2, Iso2, Point2
from engeom.geom3 import Curve3, CurveGroup3, Mesh3, Plane3, Point3


def square2(x0: float = 0.0, y0: float = 0.0) -> Curve2:
    """A closed unit square, counter-clockwise, with a length of 4."""
    points = numpy.array(
        [[x0, y0], [x0 + 1, y0], [x0 + 1, y0 + 1], [x0, y0 + 1], [x0, y0]], dtype=float
    )
    return Curve2(points, tol=1e-8)


def segment2() -> Curve2:
    """An open two-vertex segment one unit long."""
    return Curve2(numpy.array([[0.0, 5.0], [1.0, 5.0]]), tol=1e-8)


def square3(z0: float = 0.0) -> Curve3:
    """A closed unit square in the plane z = z0, its first vertex repeated as its last."""
    points = numpy.array(
        [[0, 0, z0], [1, 0, z0], [1, 1, z0], [0, 1, z0], [0, 0, z0]], dtype=float
    )
    return Curve3(points, tol=1e-8)


# =================================================================================================
# Construction and the sequence protocol
# =================================================================================================


def test_group_rejects_an_empty_collection():
    """A group with no members has no bounding box and no closest point, so it is refused."""
    with pytest.raises(ValueError):
        CurveGroup2([])
    with pytest.raises(ValueError):
        CurveGroup3([])


def test_group_behaves_as_a_sequence():
    group = CurveGroup2([square2(), segment2(), square2(10.0)])

    assert len(group) == 3
    assert len(group.curves) == 3

    # Indexing, including from the end.
    assert group[0].length() == pytest.approx(4.0)
    assert group[1].length() == pytest.approx(1.0)
    assert group[-1].length() == pytest.approx(4.0)
    assert group[-2].length() == pytest.approx(1.0)

    # Iteration, in member order.
    lengths = [c.length() for c in group]
    assert lengths == pytest.approx([4.0, 1.0, 4.0])

    # And unpacking, which is how the plotting helpers take a section.
    unpacked = [*group]
    assert len(unpacked) == 3


def test_group_index_out_of_range_raises():
    group = CurveGroup2([square2()])

    with pytest.raises(IndexError):
        _ = group[1]
    with pytest.raises(IndexError):
        _ = group[-2]


def test_group_reports_whole_body_geometry():
    group = CurveGroup2([square2(), square2(10.0)])

    assert group.length() == pytest.approx(8.0)

    aabb = group.aabb
    assert aabb.min.x == pytest.approx(0.0)
    assert aabb.max.x == pytest.approx(11.0)
    assert aabb.min.y == pytest.approx(0.0)
    assert aabb.max.y == pytest.approx(1.0)


def test_closest_point_reports_its_owning_member():
    group = CurveGroup2([square2(), square2(10.0)])

    member, station = group.at_closest_to_point(Point2(-1.0, 0.5))
    assert member == 0
    assert station.point.x == pytest.approx(0.0)

    member, _ = group.at_closest_to_point(Point2(12.0, 0.5))
    assert member == 1


def test_transforming_a_group_moves_every_member():
    group = CurveGroup2([square2(), segment2()])
    moved = group.new_transformed_by(Iso2(100.0, -50.0, 0.0))

    assert len(moved) == 2
    assert moved.length() == pytest.approx(group.length())
    assert moved.aabb.min.x == pytest.approx(group.aabb.min.x + 100.0)


# =================================================================================================
# Sections and the projection into a plane
# =================================================================================================


def test_mesh_section_returns_a_curve_group():
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    section = mesh.section_with_plane(Plane3(0.0, 0.0, 1.0, 0.0), tol=1e-9)

    assert isinstance(section, CurveGroup3)
    assert len(section) == 1

    # The old list-shaped call sites still read the same way.
    assert section[0].length() == pytest.approx(12.0, abs=1e-5)


def test_a_section_that_misses_the_mesh_raises():
    """A plane which does not intersect is an error rather than an empty group."""
    mesh = Mesh3.create_box(2.0, 2.0, 2.0)

    with pytest.raises(ValueError):
        mesh.section_with_plane(Plane3(0.0, 0.0, 1.0, 50.0))


def test_projecting_onto_the_xy_plane_matches_dropping_z():
    """
    The x-y plane's frame is the identity, so the group projection and the per-curve `to_2d` have
    to agree exactly rather than merely closely.
    """
    group = CurveGroup3([square3(3.0)])
    projected = group.to_2d_in_plane(Plane3(0.0, 0.0, 1.0, 0.0))
    dropped = group[0].to_2d()

    assert len(projected) == 1
    assert numpy.allclose(projected[0].points, dropped.points, atol=1e-15)


def test_a_section_loop_projects_to_a_closed_curve():
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    plane = Plane3(0.0, 0.0, 1.0, 0.0)

    flat = mesh.section_with_plane(plane, tol=1e-9).to_2d_in_plane(plane)

    assert len(flat) == 1
    assert flat[0].is_closed
    assert flat[0].length() == pytest.approx(12.0, abs=1e-5)


def test_a_member_collapsing_under_projection_raises():
    """A curve running along the plane normal has no 2D form, and fails the whole group."""
    along_normal = Curve3(numpy.array([[5.0, 5.0, 0.0], [5.0, 5.0, 9.0]]), tol=1e-8)
    group = CurveGroup3([square3(), along_normal])

    with pytest.raises(ValueError):
        group.to_2d_in_plane(Plane3(0.0, 0.0, 1.0, 0.0))


# =================================================================================================
# Files
# =================================================================================================


def test_group2_round_trips_through_a_file(tmp_path):
    group = CurveGroup2([square2(), segment2(), square2(10.0)])
    path = tmp_path / "group.tccurve2"

    group.save_tccurve2(path, 1e-7)
    back = CurveGroup2.load_tccurve2(path)

    assert len(back) == len(group)
    # Order, closedness and shape all survive.
    assert [c.is_closed for c in back] == [True, False, True]
    for a, b in zip(group, back):
        assert a.length() == pytest.approx(b.length(), abs=1e-6)


def test_group3_round_trips_through_a_file(tmp_path):
    group = CurveGroup3([square3(0.0), square3(5.0)])
    path = tmp_path / "group.tccurve3"

    group.save_tccurve3(path, 1e-7)
    back = CurveGroup3.load_tccurve3(path)

    assert len(back) == 2
    for a, b in zip(group, back):
        assert numpy.allclose(a.points, b.points, atol=1e-6)


def test_a_single_curve_file_loads_as_a_group_of_one(tmp_path):
    """A single-curve file is just a collection of one, so the group reader accepts it."""
    path = tmp_path / "single.tccurve2"
    square2().save_tccurve2(path, 1e-7)

    group = CurveGroup2.load_tccurve2(path)
    assert len(group) == 1
    assert group[0].is_closed


def test_a_group_file_is_refused_by_the_single_curve_reader(tmp_path):
    """
    The reverse is not true and must not silently become true: reading one curve out of a
    multi-curve file would discard the rest with nothing to show for it.
    """
    path = tmp_path / "group.tccurve2"
    CurveGroup2([square2(), segment2()]).save_tccurve2(path, 1e-7)

    with pytest.raises(IOError):
        Curve2.load_tccurve2(path)


def test_a_saved_section_still_projects_to_a_closed_curve(tmp_path):
    """
    A loop is a loop only because its first and last vertices coincide, which is the only record of
    closure a `Curve3` has. If the file did not preserve it, a restored section would silently
    project to open curves.
    """
    mesh = Mesh3.create_box(2.0, 4.0, 6.0)
    plane = Plane3(0.0, 0.0, 1.0, 0.0)
    path = tmp_path / "section.tccurve3"

    mesh.section_with_plane(plane, tol=1e-9).save_tccurve3(path, 1e-9)
    restored = CurveGroup3.load_tccurve3(path)

    assert restored.to_2d_in_plane(plane)[0].is_closed
