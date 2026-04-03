"""
    Tests of basic geometric elements (planes, circles, spheres, lines, etc.) in geom3 module.
"""
import math
import pytest
from engeom.geom3 import Vector3, Point3, SurfacePoint3, Plane3, Line3, Sphere3, Circle3, Iso3


# ==============================================================================
# Plane3 tests
# ==============================================================================

def test_plane3_construction():
    # XY plane: 0x + 0y + 1z + 0 = 0
    plane = Plane3(0, 0, 1, 0)
    assert abs(plane.a) < 1e-6
    assert abs(plane.b) < 1e-6
    assert abs(plane.c - 1.0) < 1e-6
    assert abs(plane.d) < 1e-6


def test_plane3_normal_normalized():
    # Supply un-normalized normal (0, 0, 2) - should be normalized to (0, 0, 1)
    plane = Plane3(0, 0, 2, 0)
    n = plane.normal
    assert abs(n.x) < 1e-6
    assert abs(n.y) < 1e-6
    assert abs(n.z - 1.0) < 1e-6


def test_plane3_static_xy():
    plane = Plane3.xy()
    assert abs(plane.a) < 1e-6
    assert abs(plane.b) < 1e-6
    assert abs(plane.c - 1.0) < 1e-6


def test_plane3_static_xz():
    plane = Plane3.xz()
    assert abs(plane.a) < 1e-6
    assert abs(plane.b - 1.0) < 1e-6
    assert abs(plane.c) < 1e-6


def test_plane3_static_yz():
    plane = Plane3.yz()
    assert abs(plane.a - 1.0) < 1e-6
    assert abs(plane.b) < 1e-6
    assert abs(plane.c) < 1e-6


def test_plane3_normal_property_type():
    plane = Plane3(0, 0, 1, 0)
    assert isinstance(plane.normal, Vector3)


def test_plane3_inverted_normal():
    plane = Plane3(0, 0, 1, 0)
    inv = plane.inverted_normal()
    assert abs(inv.c + 1.0) < 1e-6


def test_plane3_signed_distance_positive_side():
    # XY plane at z=0; point at z=3 is on positive side
    plane = Plane3.xy()
    point = Point3(0, 0, 3)
    assert abs(plane.signed_distance_to_point(point) - 3.0) < 1e-6


def test_plane3_signed_distance_negative_side():
    plane = Plane3.xy()
    point = Point3(0, 0, -5)
    assert abs(plane.signed_distance_to_point(point) + 5.0) < 1e-6


def test_plane3_distance_to_point_unsigned():
    plane = Plane3.xy()
    point = Point3(1, 2, -4)
    assert abs(plane.distance_to_point(point) - 4.0) < 1e-6


def test_plane3_point_is_positive():
    plane = Plane3.xy()
    assert plane.point_is_positive(Point3(0, 0, 1))
    assert plane.point_is_positive(Point3(0, 0, 0))  # on plane → signed dist = 0, >= 0
    assert not plane.point_is_positive(Point3(0, 0, -1))


def test_plane3_project_point():
    # Project (1, 2, 5) onto XY plane → (1, 2, 0)
    plane = Plane3.xy()
    projected = plane.project_point(Point3(1, 2, 5))
    assert isinstance(projected, Point3)
    assert abs(projected.x - 1.0) < 1e-6
    assert abs(projected.y - 2.0) < 1e-6
    assert abs(projected.z) < 1e-6


def test_plane3_shifted():
    # Shift XY plane by +3 → plane at z=3
    plane = Plane3.xy()
    shifted = plane.shifted(3.0)
    assert isinstance(shifted, Plane3)
    # Point at z=3 should now be on the plane (distance ~0)
    assert abs(shifted.distance_to_point(Point3(0, 0, 3))) < 1e-6


def test_plane3_intersection_distance():
    # Surface point at origin pointing along +Z, intersects XY plane at distance 0
    plane = Plane3.xy()
    sp = SurfacePoint3(0, 0, 5, 0, 0, 1)
    dist = plane.intersection_distance(sp)
    assert dist is not None
    assert abs(dist + 5.0) < 1e-6


def test_plane3_intersection_distance_parallel_returns_none():
    # Surface point normal is parallel to XY plane, no intersection
    plane = Plane3.xy()
    sp = SurfacePoint3(0, 0, 1, 1, 0, 0)
    assert plane.intersection_distance(sp) is None


def test_plane3_intersect_plane_returns_line3():
    # XY and XZ planes intersect along the X axis
    xy = Plane3.xy()
    xz = Plane3.xz()
    line = xy.intersect_plane(xz)
    assert line is not None
    assert isinstance(line, Line3)


def test_plane3_intersect_parallel_planes_returns_none():
    # Two parallel planes (both horizontal, different offsets) should not intersect
    p1 = Plane3(0, 0, 1, 0)
    p2 = Plane3(0, 0, 1, -3)
    assert p1.intersect_plane(p2) is None


# ==============================================================================
# Line3 tests
# ==============================================================================

def test_line3_construction():
    line = Line3(1, 2, 3, 4, 5, 6)
    assert isinstance(line.origin, Point3)
    assert isinstance(line.direction, Vector3)
    assert abs(line.origin.x - 1) < 1e-6
    assert abs(line.origin.y - 2) < 1e-6
    assert abs(line.origin.z - 3) < 1e-6
    assert abs(line.direction.x - 4) < 1e-6
    assert abs(line.direction.y - 5) < 1e-6
    assert abs(line.direction.z - 6) < 1e-6


def test_line3_from_points():
    p1 = Point3(0, 0, 0)
    p2 = Point3(1, 0, 0)
    line = Line3.from_points(p1, p2)
    assert isinstance(line, Line3)
    assert abs(line.origin.x) < 1e-6
    assert abs(line.direction.x - 1.0) < 1e-6


def test_line3_at():
    # Line along X axis from origin with direction (1, 0, 0)
    line = Line3(0, 0, 0, 1, 0, 0)
    p = line.at(5.0)
    assert isinstance(p, Point3)
    assert abs(p.x - 5.0) < 1e-6
    assert abs(p.y) < 1e-6
    assert abs(p.z) < 1e-6


def test_line3_scalar_project():
    line = Line3(0, 0, 0, 1, 0, 0)
    # Point (3, 4, 0): closest point on line is (3, 0, 0), t = 3
    t = line.scalar_project(Point3(3, 4, 0))
    assert abs(t - 3.0) < 1e-6


def test_line3_closest_point():
    line = Line3(0, 0, 0, 1, 0, 0)
    pt = line.closest_point(Point3(3, 4, 0))
    assert isinstance(pt, Point3)
    assert abs(pt.x - 3.0) < 1e-6
    assert abs(pt.y) < 1e-6
    assert abs(pt.z) < 1e-6


def test_line3_distance_to():
    # Line along X axis; point at (0, 3, 4) is 5 units away
    line = Line3(0, 0, 0, 1, 0, 0)
    dist = line.distance_to(Point3(0, 3, 4))
    assert abs(dist - 5.0) < 1e-6


def test_line3_intersect_plane():
    # Line along +Z axis intersects XY plane at z=0, t=0 from z=5 → t should be -5
    line = Line3(0, 0, 5, 0, 0, 1)
    plane = Plane3.xy()
    t = line.intersect_plane(plane)
    assert t is not None
    assert abs(t + 5.0) < 1e-6


def test_line3_intersect_plane_parallel_returns_none():
    # Line parallel to XY plane; no intersection
    line = Line3(0, 0, 1, 1, 0, 0)
    plane = Plane3.xy()
    assert line.intersect_plane(plane) is None


def test_line3_project_onto_plane():
    # Line with a Z component projected onto XY plane loses Z from direction
    line = Line3(0, 0, 5, 1, 0, 1)
    plane = Plane3.xy()
    projected = line.project_onto_plane(plane)
    assert projected is not None
    assert isinstance(projected, Line3)
    assert abs(projected.direction.z) < 1e-6


def test_line3_project_onto_plane_perpendicular_returns_none():
    # Line perpendicular to XY plane projects to a point, not a line
    line = Line3(0, 0, 5, 0, 0, 1)
    plane = Plane3.xy()
    assert line.project_onto_plane(plane) is None


def test_line3_normalized():
    line = Line3(0, 0, 0, 3, 0, 0)
    norm = line.normalized()
    assert isinstance(norm, Line3)
    length = math.sqrt(norm.direction.x**2 + norm.direction.y**2 + norm.direction.z**2)
    assert abs(length - 1.0) < 1e-6


# ==============================================================================
# Sphere3 tests
# ==============================================================================

def test_sphere3_construction():
    s = Sphere3(1, 2, 3, 5)
    assert isinstance(s.center, Point3)
    assert abs(s.center.x - 1) < 1e-6
    assert abs(s.center.y - 2) < 1e-6
    assert abs(s.center.z - 3) < 1e-6
    assert abs(s.r - 5) < 1e-6


def test_sphere3_closest_point():
    # Sphere at origin, radius 5; point at (10, 0, 0) → closest point at (5, 0, 0)
    s = Sphere3(0, 0, 0, 5)
    sp = s.closest_point(Point3(10, 0, 0))
    assert sp is not None
    assert isinstance(sp, SurfacePoint3)
    assert abs(sp.point.x - 5.0) < 1e-6
    assert abs(sp.point.y) < 1e-6
    assert abs(sp.point.z) < 1e-6


def test_sphere3_closest_point_normal_outward():
    # Normal at (5, 0, 0) should point outward (+X direction)
    s = Sphere3(0, 0, 0, 5)
    sp = s.closest_point(Point3(10, 0, 0))
    assert sp is not None
    assert abs(sp.normal.x - 1.0) < 1e-6
    assert abs(sp.normal.y) < 1e-6
    assert abs(sp.normal.z) < 1e-6


def test_sphere3_closest_point_at_center_returns_none():
    s = Sphere3(0, 0, 0, 5)
    assert s.closest_point(Point3(0, 0, 0)) is None


def test_sphere3_intersect_plane_returns_circle3():
    # Unit sphere at origin intersected by XY plane → circle at z=0, radius=1
    s = Sphere3(0, 0, 0, 1)
    circle = s.intersect_plane(Plane3.xy())
    assert circle is not None
    assert isinstance(circle, Circle3)
    assert abs(circle.r - 1.0) < 1e-6
    assert abs(circle.center.z) < 1e-6


def test_sphere3_intersect_plane_miss_returns_none():
    # Plane that doesn't touch the sphere
    s = Sphere3(0, 0, 0, 1)
    plane = Plane3(0, 0, 1, -5)  # z = 5, far from unit sphere at origin
    assert s.intersect_plane(plane) is None


def test_sphere3_intersect_sphere_returns_circle3():
    # Two unit spheres offset by 1 along X; they overlap
    s1 = Sphere3(0, 0, 0, 1)
    s2 = Sphere3(1, 0, 0, 1)
    circle = s1.intersect_sphere(s2)
    assert circle is not None
    assert isinstance(circle, Circle3)


def test_sphere3_intersect_sphere_no_overlap_returns_none():
    # Two unit spheres 10 apart; no intersection
    s1 = Sphere3(0, 0, 0, 1)
    s2 = Sphere3(10, 0, 0, 1)
    assert s1.intersect_sphere(s2) is None


# ==============================================================================
# Circle3 tests
# ==============================================================================

def test_circle3_construction():
    # Circle in XY plane, center at origin, radius 3
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    assert abs(c.r - 3.0) < 1e-6
    assert abs(c.center.x) < 1e-6
    assert abs(c.center.y) < 1e-6
    assert abs(c.center.z) < 1e-6


def test_circle3_normal_normalized():
    # Supply un-normalized normal (0, 0, 5); should normalize to (0, 0, 1)
    c = Circle3(0, 0, 0, 0, 0, 5, 3)
    n = c.normal
    assert abs(n.z - 1.0) < 1e-6


def test_circle3_normal_type():
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    assert isinstance(c.normal, Vector3)


def test_circle3_plane_property():
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    assert isinstance(c.plane, Plane3)


def test_circle3_zero_normal_raises():
    with pytest.raises(Exception):
        Circle3(0, 0, 0, 0, 0, 0, 3)


def test_circle3_at_angle_returns_surface_point():
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    sp = c.at_angle(0.0)
    assert isinstance(sp, SurfacePoint3)


def test_circle3_at_angle_point_on_circumference():
    # Circle at origin in XY plane, radius 3; angle=0 → point at (3, 0, 0)
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    sp = c.at_angle(0.0)
    dist = math.sqrt((sp.point.x - 0)**2 + (sp.point.y - 0)**2 + (sp.point.z - 0)**2)
    assert abs(dist - 3.0) < 1e-6


def test_circle3_at_angle_full_revolution():
    # angle=0 and angle=2π should give the same point
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    sp0 = c.at_angle(0.0)
    sp2pi = c.at_angle(2 * math.pi)
    assert abs(sp0.point.x - sp2pi.point.x) < 1e-6
    assert abs(sp0.point.y - sp2pi.point.y) < 1e-6
    assert abs(sp0.point.z - sp2pi.point.z) < 1e-6


def test_circle3_closest_position():
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    sp = c.closest_position(Point3(3, 0, 5))
    assert isinstance(sp, SurfacePoint3)
    dist = math.sqrt(sp.point.x**2 + sp.point.y**2 + sp.point.z**2)
    assert abs(dist - 3.0) < 1e-6


def test_circle3_intersect_plane_two_intersections():
    # Circle in XY plane at origin, radius 3; YZ plane cuts through it at two points
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    angles = c.intersect_plane(Plane3.yz())
    assert len(angles) == 2


def test_circle3_intersect_plane_no_intersection():
    # Circle in XY plane at origin, radius 1; plane at z=5 doesn't intersect
    c = Circle3(0, 0, 0, 0, 0, 1, 1)
    plane = Plane3(0, 0, 1, -5)
    angles = c.intersect_plane(plane)
    assert len(angles) == 0


def test_circle3_intersect_plane_angles_are_floats():
    c = Circle3(0, 0, 0, 0, 0, 1, 3)
    angles = c.intersect_plane(Plane3.yz())
    for a in angles:
        assert isinstance(a, float)


# ==============================================================================
# Iso3 transformation tests
# ==============================================================================

def test_iso3_matmul_plane3_returns_plane3():
    iso = Iso3.identity()
    plane = Plane3.xy()
    result = iso @ plane
    assert isinstance(result, Plane3)


def test_iso3_matmul_line3_returns_line3():
    iso = Iso3.identity()
    line = Line3(0, 0, 0, 1, 0, 0)
    result = iso @ line
    assert isinstance(result, Line3)


def test_iso3_matmul_sphere3_returns_sphere3():
    iso = Iso3.identity()
    sphere = Sphere3(0, 0, 0, 5)
    result = iso @ sphere
    assert isinstance(result, Sphere3)


def test_iso3_matmul_circle3_returns_circle3():
    iso = Iso3.identity()
    circle = Circle3(0, 0, 0, 0, 0, 1, 3)
    result = iso @ circle
    assert isinstance(result, Circle3)


def test_iso3_translation_transforms_plane3():
    # Translate XY plane (z=0) by +5 along Z → plane should now sit at z=5
    iso = Iso3.from_translation(0, 0, 5)
    plane = Plane3.xy()
    result = iso @ plane
    assert isinstance(result, Plane3)
    # Origin (0,0,0) should now be 5 units away (negative side)
    assert abs(result.signed_distance_to_point(Point3(0, 0, 0)) + 5.0) < 1e-6
    # Point at z=5 should lie on the new plane
    assert abs(result.distance_to_point(Point3(0, 0, 5))) < 1e-6


def test_iso3_translation_transforms_line3():
    # Line along X axis translated by (0, 3, 0) → origin moves to (0, 3, 0)
    iso = Iso3.from_translation(0, 3, 0)
    line = Line3(0, 0, 0, 1, 0, 0)
    result = iso @ line
    assert abs(result.origin.x) < 1e-6
    assert abs(result.origin.y - 3.0) < 1e-6
    assert abs(result.origin.z) < 1e-6
    # Direction should be unchanged
    assert abs(result.direction.x - 1.0) < 1e-6
    assert abs(result.direction.y) < 1e-6
    assert abs(result.direction.z) < 1e-6


def test_iso3_translation_transforms_sphere3():
    # Sphere at origin translated by (1, 2, 3) → center moves to (1, 2, 3)
    iso = Iso3.from_translation(1, 2, 3)
    sphere = Sphere3(0, 0, 0, 5)
    result = iso @ sphere
    assert abs(result.center.x - 1.0) < 1e-6
    assert abs(result.center.y - 2.0) < 1e-6
    assert abs(result.center.z - 3.0) < 1e-6
    # Radius is preserved
    assert abs(result.r - 5.0) < 1e-6


def test_iso3_translation_transforms_circle3():
    # Circle at origin translated by (1, 2, 3) → center moves to (1, 2, 3)
    iso = Iso3.from_translation(1, 2, 3)
    circle = Circle3(0, 0, 0, 0, 0, 1, 3)
    result = iso @ circle
    assert pytest.approx(1.0) == result.center.x
    assert pytest.approx(2.0) == result.center.y
    assert pytest.approx(3.0) == result.center.z

    # Radius and normal are preserved
    assert pytest.approx(3.0) == result.r
    assert pytest.approx(1.0) == result.normal.z


def test_iso3_rotation_transforms_plane3():
    # Rotate XY plane (normal +Z) by 90° around X → normal becomes -Y
    iso = Iso3.from_rx(math.pi / 2)
    plane = Plane3.xy()
    result = iso @ plane
    assert pytest.approx(0.0) == result.a
    assert pytest.approx(-1.0) == result.b
    assert pytest.approx(0.0) == result.c


def test_iso3_rotation_transforms_line3_direction():
    # Line along +X rotated 90° around Z → direction becomes +Y
    iso = Iso3.from_rz(math.pi / 2)
    line = Line3(0, 0, 0, 1, 0, 0)
    result = iso @ line
    assert abs(result.direction.x) < 1e-6
    assert abs(result.direction.y - 1.0) < 1e-6
    assert abs(result.direction.z) < 1e-6


def test_iso3_rotation_transforms_sphere3_center():
    # Sphere at (1, 0, 0) rotated 90° around Z → center moves to (0, 1, 0)
    iso = Iso3.from_rz(math.pi / 2)
    sphere = Sphere3(1, 0, 0, 2)
    result = iso @ sphere
    assert abs(result.center.x) < 1e-6
    assert abs(result.center.y - 1.0) < 1e-6
    assert abs(result.center.z) < 1e-6
    assert abs(result.r - 2.0) < 1e-6


def test_iso3_rotation_transforms_circle3_normal():
    # Circle with normal +Z rotated 90° around X → normal becomes -Y
    iso = Iso3.from_rx(math.pi / 2)
    circle = Circle3(0, 0, 0, 0, 0, 1, 3)
    result = iso @ circle
    assert pytest.approx(0.0) == result.normal.x
    assert pytest.approx(-1.0) == result.normal.y
    assert pytest.approx(0.0) == result.normal.z


def test_iso3_identity_preserves_plane3():
    iso = Iso3.identity()
    plane = Plane3(1, 2, 3, -4)
    result = iso @ plane
    assert abs(result.a - plane.a) < 1e-6
    assert abs(result.b - plane.b) < 1e-6
    assert abs(result.c - plane.c) < 1e-6
    assert abs(result.d - plane.d) < 1e-6


def test_iso3_identity_preserves_sphere3():
    iso = Iso3.identity()
    sphere = Sphere3(1, 2, 3, 7)
    result = iso @ sphere
    assert abs(result.center.x - 1.0) < 1e-6
    assert abs(result.center.y - 2.0) < 1e-6
    assert abs(result.center.z - 3.0) < 1e-6
    assert abs(result.r - 7.0) < 1e-6
