"""
Tests for the pure-Python plotting helpers in `engeom.plot`.

These focus on the binding-adjacent behavior the helpers own: coercing loose point-like arguments,
the bookkeeping that turns disjoint runs into a single Matplotlib artist, and the contracts the
helpers advertise in their signatures. Rendering correctness is not checked; the goal is that every
public entry point is exercised and that its documented edge cases hold.

Matplotlib is an optional dependency, so the whole module skips when it is absent.
"""

import numpy
import pytest

pytest.importorskip("matplotlib")

import matplotlib

matplotlib.use("Agg")

from engeom.geom2 import Point2, Vector2
from engeom.geom3 import Point3
from engeom.plot.matplotlib import TraceBuilder

TOL = 1e-12


# ---------------------------------------------------------------------------
# TraceBuilder: list alignment
# ---------------------------------------------------------------------------

def test_trace_builder_starts_empty():
    b = TraceBuilder()
    assert b.xs == []
    assert b.ys == []
    assert b.c == []


def test_add_points_keeps_all_three_lists_aligned():
    b = TraceBuilder()
    b.add_points((0.0, 1.0), (2.0, 3.0))
    assert b.xs == [0.0, 2.0]
    assert b.ys == [1.0, 3.0]
    assert b.c == [None, None]


def test_mixing_colored_and_uncolored_points_keeps_lists_aligned():
    """ Regression: add_points used to skip `c`, silently desyncing it from `xs` and `ys`. """
    b = TraceBuilder()
    b.add_points((0.0, 0.0))
    b.add_point_and_color((1.0, 1.0), 5.0)
    b.add_points((2.0, 2.0))
    b.add_blank()
    b.add_point_and_color((3.0, 3.0), 7.0)

    assert len(b.xs) == len(b.ys) == len(b.c)
    assert b.c == [None, 5.0, None, None, 7.0]
    # The color of a point is found at the same index as its coordinates.
    assert b.xs[b.c.index(5.0)] == 1.0
    assert b.xs[b.c.index(7.0)] == 3.0


def test_add_segment_terminates_the_run_with_a_break():
    b = TraceBuilder()
    b.add_segment((0.0, 0.0), (1.0, 1.0))
    b.add_segment((2.0, 2.0), (3.0, 3.0))

    assert b.xs == [0.0, 1.0, None, 2.0, 3.0, None]
    assert b.ys == [0.0, 1.0, None, 2.0, 3.0, None]
    assert b.c == [None] * 6


def test_add_blank_appends_to_all_three_lists():
    b = TraceBuilder()
    b.add_blank()
    assert b.xs == [None]
    assert b.ys == [None]
    assert b.c == [None]


# ---------------------------------------------------------------------------
# TraceBuilder: accepted point types
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("point", [
    Point2(1.0, 2.0),
    Vector2(1.0, 2.0),
    Point3(1.0, 2.0, 9.0),
    (1.0, 2.0),
    (1.0, 2.0, 9.0),
    [1.0, 2.0],
    numpy.array([1.0, 2.0]),
    numpy.array([1.0, 2.0, 9.0]),
])
def test_add_points_accepts_any_point_like_and_ignores_z(point):
    b = TraceBuilder()
    b.add_points(point)
    assert b.xs[0] == pytest.approx(1.0, abs=TOL)
    assert b.ys[0] == pytest.approx(2.0, abs=TOL)


def test_add_point_and_color_accepts_an_engeom_point():
    """ Regression: this used to subscript the point, which `Point2` does not support. """
    b = TraceBuilder()
    b.add_point_and_color(Point2(1.0, 2.0), 3.0)
    assert b.xs == [1.0]
    assert b.ys == [2.0]
    assert b.c == [3.0]


# ---------------------------------------------------------------------------
# TraceBuilder: bounds
# ---------------------------------------------------------------------------

def test_bounds_of_an_empty_builder_is_none():
    """ Regression: this used to raise a bare ValueError out of min() on an empty sequence. """
    assert TraceBuilder().bounds() is None


def test_bounds_of_a_builder_holding_only_breaks_is_none():
    b = TraceBuilder()
    b.add_blank()
    b.add_blank()
    assert b.bounds() is None


def test_bounds_covers_every_point_and_ignores_breaks():
    b = TraceBuilder()
    b.add_segment((0.0, 1.0), (4.0, -2.0))
    b.add_segment((-1.0, 3.0), (2.0, 0.0))

    box = b.bounds()
    assert box.min.x == pytest.approx(-1.0, abs=TOL)
    assert box.max.x == pytest.approx(4.0, abs=TOL)
    assert box.min.y == pytest.approx(-2.0, abs=TOL)
    assert box.max.y == pytest.approx(3.0, abs=TOL)


def test_bounds_of_a_single_point_is_degenerate_but_valid():
    b = TraceBuilder()
    b.add_points((2.0, 5.0))

    box = b.bounds()
    assert box.min.x == pytest.approx(2.0, abs=TOL)
    assert box.max.x == pytest.approx(2.0, abs=TOL)
    assert box.min.y == pytest.approx(5.0, abs=TOL)
    assert box.max.y == pytest.approx(5.0, abs=TOL)


# ---------------------------------------------------------------------------
# TraceBuilder: transforms and output
# ---------------------------------------------------------------------------

def test_invert_y_negates_coordinates_and_preserves_breaks():
    b = TraceBuilder()
    b.add_segment((0.0, 1.0), (2.0, -3.0))
    b.invert_y()

    assert b.ys == [-1.0, 3.0, None]
    assert b.xs == [0.0, 2.0, None]


def test_invert_y_is_its_own_inverse():
    b = TraceBuilder()
    b.add_points((0.0, 1.0), (2.0, -3.0))
    before = list(b.ys)
    b.invert_y()
    b.invert_y()
    assert b.ys == before


def test_xy_returns_the_coordinate_lists_for_unpacking():
    b = TraceBuilder()
    b.add_segment((0.0, 1.0), (2.0, 3.0))

    xs, ys = b.xy
    assert xs is b.xs
    assert ys is b.ys


def test_xy_unpacks_into_a_matplotlib_plot_call():
    from matplotlib.pyplot import figure

    b = TraceBuilder()
    b.add_segment((0.0, 0.0), (1.0, 1.0))
    b.add_segment((2.0, 2.0), (3.0, 3.0))

    ax = figure().subplots()
    lines = ax.plot(*b.xy, color="black")
    assert len(lines) == 1
