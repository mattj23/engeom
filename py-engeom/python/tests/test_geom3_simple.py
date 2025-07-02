"""
    Tests of simple objects in geom3 module.
"""
import pytest
import numpy
from engeom.geom3 import Vector3, Point3, SurfacePoint3, Iso3


def test_unpacking():
    def some_function(a: float, b: float, c: float) -> float:
        return a + b + c

    p = Point3(1, 2, 3)

    x, y, z = p
    assert x == 1
    assert y == 2
    assert z == 3

    coords = list(p)
    assert coords == [1, 2, 3]

    value = some_function(*p)
    assert value == 6


def test_vector_mul_scalar():
    v = Vector3(1, 2, 3)
    result = v * 3
    assert abs(result.x - 3) < 1e-6
    assert abs(result.y - 6) < 1e-6
    assert abs(result.z - 9) < 1e-6


def test_vector_div_scalar():
    v = Vector3(3, 6, 9)
    result = v / 3
    assert abs(result.x - 1) < 1e-6
    assert abs(result.y - 2) < 1e-6
    assert abs(result.z - 3) < 1e-6


def test_point_mul_scalar():
    p = Point3(1, 2, 3)
    result = p * 3
    assert abs(result.x - 3) < 1e-6
    assert abs(result.y - 6) < 1e-6
    assert abs(result.z - 9) < 1e-6


def test_point_div_scalar():
    p = Point3(3, 6, 9)
    result = p / 3
    assert abs(result.x - 1) < 1e-6
    assert abs(result.y - 2) < 1e-6
    assert abs(result.z - 3) < 1e-6


def test_sp_mul_scalar_pos():
    sp = SurfacePoint3(1, 2, 3, 1, 0, 0)
    result = sp * 3
    assert abs(result.point.x - 3) < 1e-6
    assert abs(result.point.y - 6) < 1e-6
    assert abs(result.point.z - 9) < 1e-6

    assert abs(result.normal.x - 1) < 1e-6
    assert abs(result.normal.y) < 1e-6
    assert abs(result.normal.z) < 1e-6


def test_sp_mul_scalar_neg():
    sp = SurfacePoint3(1, 2, 3, 1, 0, 0)
    result = sp * -3
    assert abs(result.point.x + 3) < 1e-6
    assert abs(result.point.y + 6) < 1e-6
    assert abs(result.point.z + 9) < 1e-6

    assert abs(result.normal.x + 1) < 1e-6
    assert abs(result.normal.y) < 1e-6
    assert abs(result.normal.z) < 1e-6


def test_sp_div_scalar():
    sp = SurfacePoint3(3, 6, 9, 1, 0, 0)
    result = sp / 3
    assert abs(result.point.x - 1) < 1e-6
    assert abs(result.point.y - 2) < 1e-6
    assert abs(result.point.z - 3) < 1e-6

    assert abs(result.normal.x - 1) < 1e-6
    assert abs(result.normal.y) < 1e-6
    assert abs(result.normal.z) < 1e-6


# Test that a vector plus a vector is a vector.
def test_vector_plus_vector():
    v1 = Vector3(1, 2, 3)
    v2 = Vector3(4, 5, 6)
    v3 = v1 + v2
    assert isinstance(v3, Vector3)


def test_vector_plus_point():
    v = Vector3(1, 2, 3)
    p = Point3(4, 5, 6)
    result = v + p
    assert isinstance(result, Point3)


# Test that a point plus a vector is a point.
def test_point_plus_vector():
    p = Point3(1, 2, 3)
    v = Vector3(4, 5, 6)
    result = p + v
    assert isinstance(result, Point3)


# Test that a point minus a point is a vector.
def test_point_minus_point():
    p1 = Point3(1, 2, 3)
    p2 = Point3(4, 5, 6)
    result = p1 - p2
    assert isinstance(result, Vector3)


# Test that a point minus a vector is a point.
def test_point_minus_vector():
    p = Point3(4, 5, 6)
    v = Vector3(1, 2, 3)
    result = p - v
    assert isinstance(result, Point3)


# Test that a vector minus a vector is a vector.
def test_vector_minus_vector():
    v1 = Vector3(4, 5, 6)
    v2 = Vector3(1, 2, 3)
    result = v1 - v2
    assert isinstance(result, Vector3)


# Test that an Iso3 matmul by a vector returns a vector.
def test_iso3_matmul_vector():
    iso = Iso3.identity()
    v = Vector3(4, 5, 6)
    result = iso @ v
    assert isinstance(result, Vector3)


# Test that an Iso3 matmul by a point returns a point.
def test_iso3_matmul_point():
    iso = Iso3.identity()
    p = Point3(4, 5, 6)
    result = iso @ p
    assert isinstance(result, Point3)


# Test that an Iso3 matmul by a surface point returns a surface point.
def test_iso3_matmul_surfacepoint():
    iso = Iso3.identity()
    sp = SurfacePoint3(4, 5, 6, 1, 0, 0)
    result = iso @ sp
    assert isinstance(result, SurfacePoint3)


# Test that an Iso3 matmul by another Iso3 returns an Iso3.
def test_iso3_matmul_iso3():
    iso1 = Iso3.identity()
    iso2 = Iso3.identity()
    result = iso1 @ iso2
    assert isinstance(result, Iso3)
