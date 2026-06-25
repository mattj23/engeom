"""
    Tests of cubic spline closest-point queries and fitting in the geom2 and geom3 modules.
"""
import numpy
import pytest

from engeom.geom2 import (
    CubicSpline2,
    CubicSplineQueries2,
    SplineProjection,
    Point2,
    fit_spline_to_points as fit_spline_to_points_2,
)
from engeom.geom3 import (
    CubicSpline3,
    CubicSplineQueries3,
    Point3,
    fit_spline_to_points as fit_spline_to_points_3,
)


# ==============================================================================
# 2D queries
# ==============================================================================

def test_spline2_query_project_point():
    spline = CubicSpline2(0, 0, 1, 2, 2, 2, 3, 0)
    queries = spline.query()
    assert isinstance(queries, CubicSplineQueries2)

    proj = queries.project_point(Point2(1.5, 1.0))
    assert isinstance(proj, SplineProjection)
    # The curve is symmetric about x = 1.5, so the closest point is the apex at t = 0.5.
    assert proj.t == pytest.approx(0.5, abs=1e-9)
    assert proj.distance == pytest.approx(spline.position(0.5).y - 1.0, abs=1e-9)


def test_spline2_constructor_matches_method():
    spline = CubicSpline2(0, 0, 1, 2, 2, 2, 3, 0)
    from_method = spline.query().project_point(Point2(0.7, 0.3))
    from_ctor = CubicSplineQueries2(spline).project_point(Point2(0.7, 0.3))
    assert from_ctor.t == pytest.approx(from_method.t)
    assert from_ctor.distance == pytest.approx(from_method.distance)


def test_spline2_project_point_convenience():
    spline = CubicSpline2(0, 0, 1, 2, 2, 2, 3, 0)
    direct = spline.project_point(Point2(1.5, 1.0))
    via_query = spline.query().project_point(Point2(1.5, 1.0))
    assert direct.t == pytest.approx(via_query.t)
    assert direct.distance == pytest.approx(via_query.distance)


def test_spline2_project_points_vectorized():
    spline = CubicSpline2(0, 0, 1, 2, 2, 2, 3, 0)
    queries = spline.query()
    pts = numpy.array([[1.5, 1.0], [0.0, 0.0], [3.0, 0.0]])
    out = queries.project_points(pts)
    assert out.shape == (3, 2)
    # Endpoints project onto the curve ends at zero distance.
    assert out[1, 0] == pytest.approx(0.0, abs=1e-9)
    assert out[1, 1] == pytest.approx(0.0, abs=1e-9)
    assert out[2, 0] == pytest.approx(1.0, abs=1e-9)
    assert out[2, 1] == pytest.approx(0.0, abs=1e-9)


# ==============================================================================
# 3D queries
# ==============================================================================

def test_spline3_query_project_point():
    spline = CubicSpline3(0, 0, 0, 1, 2, 1, 2, 2, -1, 3, 0, 0)
    queries = spline.query()
    assert isinstance(queries, CubicSplineQueries3)

    apex = spline.position(0.5)
    proj = queries.project_point(Point3(apex.x, apex.y, apex.z))
    assert proj.t == pytest.approx(0.5, abs=1e-6)
    assert proj.distance == pytest.approx(0.0, abs=1e-9)


def test_spline3_project_points_vectorized():
    spline = CubicSpline3(0, 0, 0, 1, 2, 1, 2, 2, -1, 3, 0, 0)
    pts = numpy.array([[apex_c for apex_c in (spline.position(t / 10.0).x,
                                              spline.position(t / 10.0).y,
                                              spline.position(t / 10.0).z)]
                       for t in range(11)])
    out = spline.query().project_points(pts)
    assert out.shape == (11, 2)
    # Every sample lies on the curve, so all distances are ~zero.
    assert numpy.max(out[:, 1]) < 1e-9


# ==============================================================================
# Fitting
# ==============================================================================

def _samples_2d(spline, n=21):
    return numpy.array([[spline.position(i / (n - 1)).x, spline.position(i / (n - 1)).y]
                        for i in range(n)])


def _samples_3d(spline, n=41):
    return numpy.array([[spline.position(i / (n - 1)).x,
                         spline.position(i / (n - 1)).y,
                         spline.position(i / (n - 1)).z]
                        for i in range(n)])


def test_fit_spline2_recovers_control_points():
    truth = CubicSpline2(0, 0, 1, 2, 2, 2, 3, 0)
    points = _samples_2d(truth)

    # Endpoints fixed, inner control points free.
    def builder(p):
        return CubicSpline2(0.0, 0.0, p[0], p[1], p[2], p[3], 3.0, 0.0)

    result = fit_spline_to_points_2(points, builder, numpy.array([1.0, 0.0, 2.0, 0.0]))
    assert numpy.allclose(result, [1.0, 2.0, 2.0, 2.0], atol=1e-5)


def test_fit_spline3_recovers_control_points():
    truth = CubicSpline3(0, 0, 0, 1, 2, 1, 2, 2, -1, 3, 0, 0)
    points = _samples_3d(truth)

    def builder(p):
        return CubicSpline3(0, 0, 0, p[0], p[1], p[2], p[3], p[4], p[5], 3, 0, 0)

    result = fit_spline_to_points_3(
        points, builder, numpy.array([1.0, 1.5, 0.5, 2.0, 1.5, -0.5])
    )
    assert numpy.allclose(result, [1, 2, 1, 2, 2, -1], atol=1e-5)


def test_fit_spline2_bad_builder_raises():
    points = numpy.array([[0.0, 0.0], [1.0, 0.0]])

    def builder(_p):
        return "not a spline"

    with pytest.raises(ValueError):
        fit_spline_to_points_2(points, builder, numpy.array([0.0, 0.0, 0.0, 0.0]))
