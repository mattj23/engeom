"""
    Tests of the airfoil2 measurement bindings.

    These use a synthetic elliptical section rather than a data file so that the expected values
    are known analytically. An ellipse is a symmetric airfoil: its camber line is the major axis,
    so the thickness at a station is just the ellipse's vertical extent there, and the maximum
    thickness is the minor axis.
"""

import math

import numpy
import pytest
from engeom.airfoil2 import AfGeometry, OrientFwdAft, OrientUpperLower
from engeom.geom2 import Curve2

# Semi-major (chordwise) and semi-minor (thickness) axes of the test section
A = 10.0
B = 3.0


@pytest.fixture(scope="module")
def geom() -> AfGeometry:
    # Points around the ellipse, counter-clockwise, with the first point repeated at the end so
    # that the curve closes
    t = numpy.linspace(0.0, 2.0 * math.pi, 721)
    points = numpy.column_stack((A * numpy.cos(t), B * numpy.sin(t)))
    section = Curve2(points, tol=1e-8)

    # The section is symmetric, so the maximum thickness cannot pick out a leading edge and the
    # camber has no curvature to read. Both orientations are given explicitly instead.
    return AfGeometry.from_geometric_analysis(
        section=section,
        general_tol=1e-3,
        fwd_aft=OrientFwdAft.fwd(-1.0, 0.0),
        upper_lower=OrientUpperLower.upper(0.0, 1.0),
        le_search="auto",
        te_search="auto",
    )


def test_af_point_lands_on_the_requested_side(geom: AfGeometry):
    mid = geom.camber.length / 2.0
    upper = geom.af_point("upper", "on_camber", mid)
    lower = geom.af_point("lower", "on_camber", mid)

    assert upper is not None and lower is not None

    # At mid-chord the ellipse is at its full height
    assert upper.point.y == pytest.approx(B, abs=1e-3)
    assert lower.point.y == pytest.approx(-B, abs=1e-3)
    assert upper.point.x == pytest.approx(0.0, abs=1e-3)


def test_af_point_off_the_surface_is_none(geom: AfGeometry):
    assert geom.af_point("upper", "on_camber", geom.camber.length * 1.5) is None


def test_af_point_rejects_bad_tokens(geom: AfGeometry):
    with pytest.raises(ValueError, match="Invalid airfoil side"):
        geom.af_point("sideways", "on_camber", 1.0)

    with pytest.raises(ValueError, match="Invalid position method"):
        geom.af_point("upper", "somewhere", 1.0)


def test_thickness_at_mid_chord(geom: AfGeometry):
    mid = geom.camber.length / 2.0
    d = geom.thickness_at("on_camber", mid)

    assert d is not None
    assert d.value == pytest.approx(2.0 * B, abs=1e-3)

    # `a` is the lower point and `b` the upper one
    assert d.a.y == pytest.approx(-B, abs=1e-3)
    assert d.b.y == pytest.approx(B, abs=1e-3)


def test_thickness_at_matches_its_own_gage_points(geom: AfGeometry):
    value = geom.camber.length * 0.3
    d = geom.thickness_at("on_camber", value)
    upper = geom.af_point("upper", "on_camber", value)
    lower = geom.af_point("lower", "on_camber", value)

    assert d is not None and upper is not None and lower is not None
    expected = math.hypot(upper.point.x - lower.point.x, upper.point.y - lower.point.y)
    assert d.value == pytest.approx(expected, abs=1e-9)


def test_thickness_at_rejects_bad_token(geom: AfGeometry):
    with pytest.raises(ValueError, match="Invalid position method"):
        geom.thickness_at("somewhere", 1.0)


def test_max_thickness_is_the_minor_axis(geom: AfGeometry):
    d = geom.max_thickness()
    assert d.value == pytest.approx(2.0 * B, abs=1e-3)
    assert d.value == pytest.approx(geom.max_thickness_circle().radius * 2.0, abs=1e-3)


def test_max_thickness_bounds_a_camber_sweep(geom: AfGeometry):
    """
    On a symmetric section the camber line is straight, which is the case where the inscribed
    circle and the orthogonal chord definitions of maximum thickness do agree. A sweep should not
    find anything meaningfully thicker.
    """
    best = geom.max_thickness().value
    length = geom.camber.length

    sweep = 0.0
    for i in range(1, 400):
        d = geom.thickness_at("on_camber", length * i / 400.0)
        if d is not None:
            sweep = max(sweep, d.value)

    assert sweep > 0.0
    assert sweep <= best + 1e-3
