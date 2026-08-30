"""
Tests for the `Curve2` binding.

These cover the enclosed-region properties of a closed curve. The rest of the `Curve2` API is
exercised indirectly by the curve group, alignment, and airfoil tests.
"""

from __future__ import annotations

import math

import numpy
import pytest

from engeom.geom2 import Curve2


def square(x0: float = 0.0, y0: float = 0.0, side: float = 1.0) -> Curve2:
    """A counter-clockwise closed square with its lower-left corner at (x0, y0)."""
    points = numpy.array(
        [[x0, y0], [x0 + side, y0], [x0 + side, y0 + side], [x0, y0 + side]], dtype=float
    )
    return Curve2(points, tol=1e-9, force_closed=True)


def l_shape() -> Curve2:
    """A closed L shape with area 3 and centroid (2.5/3, 2.5/3), rather than the vertex mean (1, 1)."""
    points = numpy.array([[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]], dtype=float)
    return Curve2(points, tol=1e-9, force_closed=True)


# =================================================================================================
# Area
# =================================================================================================


def test_signed_area_carries_the_winding():
    ccw = square()
    assert ccw.signed_area == pytest.approx(1.0)
    assert ccw.reversed().signed_area == pytest.approx(-1.0)


def test_area_ignores_the_winding():
    ccw = square()
    assert ccw.area == pytest.approx(1.0)
    assert ccw.reversed().area == pytest.approx(1.0)


def test_area_of_an_open_curve_raises():
    open_curve = Curve2(numpy.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]))
    assert not open_curve.is_closed
    with pytest.raises(ValueError):
        open_curve.area
    with pytest.raises(ValueError):
        open_curve.signed_area


# =================================================================================================
# Area centroid
# =================================================================================================


def test_area_centroid_is_not_the_vertex_mean():
    c = l_shape()
    assert c.area == pytest.approx(3.0)
    centroid = c.area_centroid
    assert centroid.x == pytest.approx(2.5 / 3.0)
    assert centroid.y == pytest.approx(2.5 / 3.0)


def test_area_centroid_holds_still_under_resampling():
    """Denser sampling moves a circle's vertex mean but must not move its area centroid."""
    center = (3.0, -2.0)
    r = 1.5
    for n in (8, 50, 1000):
        t = numpy.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        points = numpy.column_stack([center[0] + r * numpy.cos(t), center[1] + r * numpy.sin(t)])
        c = Curve2(points, tol=1e-9, force_closed=True)
        centroid = c.area_centroid
        assert centroid.x == pytest.approx(center[0], abs=1e-9)
        assert centroid.y == pytest.approx(center[1], abs=1e-9)


def test_area_centroid_of_an_open_curve_raises():
    open_curve = Curve2(numpy.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]))
    with pytest.raises(ValueError):
        open_curve.area_centroid


def test_area_centroid_of_a_zero_area_loop_raises():
    doubled_back = Curve2(
        numpy.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.0, 0.0]]), tol=1e-9, force_closed=True
    )
    assert doubled_back.is_closed
    assert doubled_back.area == pytest.approx(0.0)
    with pytest.raises(ValueError):
        doubled_back.area_centroid
