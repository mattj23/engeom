"""
Tests for the pure-Python plotting helpers in `engeom.plot`.

These focus on the binding-adjacent behavior the helpers own: coercing loose point-like arguments,
the bookkeeping that turns disjoint runs into a single Matplotlib artist, and the contracts the
helpers advertise in their signatures. Rendering correctness is not checked; the goal is that every
public entry point is exercised and that its documented edge cases hold.

Matplotlib is an optional dependency, so the whole module skips when it is absent.
"""

import warnings
from typing import get_args

import numpy
import pytest

pytest.importorskip("matplotlib")

import matplotlib

matplotlib.use("Agg")

from matplotlib.figure import Figure
from matplotlib.lines import AxLine, Line2D
from matplotlib.text import Annotation

from engeom.airfoil2 import AfGeometry, OrientFwdAft, OrientUpperLower
from engeom.geom2 import (Aabb2, Arc2, BoundaryData2, Circle2, CubicSpline2, Curve2, Line2, Point2,
                          Segment2, SurfacePoint2, Vector2)
from engeom.geom3 import (Aabb3, Circle3, Curve3, Iso3, Line3, Mesh3, Plane3, Point3, PointCloud,
                          Vector3)
from engeom.metrology import Distance2, Distance3
from engeom.plot import LabelPlace
from engeom.plot._coerce import to_point2, to_point3, to_tuple2, to_tuple3
from engeom.plot._common import LABEL_PLACES
from engeom.plot.matplotlib import (GOM_CMAP, AxesHelper, GomColorMap, TraceBuilder, ViewPort3,
                                    deviation_limit, deviation_norm, extend_for, has_extremes)
from engeom.plot.matplotlib._style import element_style, merge_style
from engeom.plot.matplotlib.viewport import _plane_basis

TOL = 1e-12


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

def new_helper() -> AxesHelper:
    """
    Build a helper over a standalone Figure. A bare `Figure` is used rather than `pyplot.figure` so
    that the tests do not accumulate state in the pyplot registry.
    """
    ax = Figure().subplots()
    helper = AxesHelper(ax)
    helper.set_bounds(Aabb2(x_min=-10.0, x_max=10.0, y_min=-10.0, y_max=10.0))
    return helper


def artist_count(ax) -> int:
    """
    Total artists on an axes across every container the helpers draw into. Counted in aggregate so
    the tests assert that something was drawn without pinning which Matplotlib container it lands in.
    """
    return len(ax.lines) + len(ax.patches) + len(ax.texts) + len(ax.collections)


def sample_curve() -> Curve2:
    return Curve2(numpy.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0]]), tol=1e-8)


def sample_boundary():
    data = BoundaryData2(0.0, 0.0)
    data.add_seg_xy(1.0, 0.0)
    data.add_arc_xy(1.0, 0.5, 1.0, 1.0, False)
    data.add_seg_xy(0.0, 1.0)
    return data.to_boundary()


# ---------------------------------------------------------------------------
# Coercion helpers
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
def test_to_tuple2_accepts_any_point_like_and_drops_extra_values(point):
    assert to_tuple2(point) == pytest.approx((1.0, 2.0), abs=TOL)


@pytest.mark.parametrize("point", [
    Point3(1.0, 2.0, 3.0),
    Vector3(1.0, 2.0, 3.0),
    (1.0, 2.0, 3.0),
    (1.0, 2.0, 3.0, 9.0),
    [1.0, 2.0, 3.0],
    numpy.array([1.0, 2.0, 3.0]),
])
def test_to_tuple3_accepts_any_point_like_and_drops_extra_values(point):
    assert to_tuple3(point) == pytest.approx((1.0, 2.0, 3.0), abs=TOL)


def test_to_point3_promotes_a_2d_point_onto_the_z_zero_plane():
    assert to_point3(Point2(1.0, 2.0)) == Point3(1.0, 2.0, 0.0)


def test_to_point3_passes_a_3d_point_through_unchanged():
    p = Point3(1.0, 2.0, 3.0)
    assert to_point3(p) is p


def test_to_point2_truncates_the_z_of_a_3d_point():
    assert to_point2(Point3(1.0, 2.0, 3.0)) == Point2(1.0, 2.0)


def test_to_point2_passes_a_2d_point_through_unchanged():
    p = Point2(1.0, 2.0)
    assert to_point2(p) is p


@pytest.mark.parametrize("short,expected", [([1.0], Point2(1.0, 0.0)), ([], Point2(0.0, 0.0))])
def test_to_point2_pads_a_short_iterable_with_zeros(short, expected):
    assert to_point2(short) == expected


@pytest.mark.parametrize("short,expected", [
    ([1.0], Point3(1.0, 0.0, 0.0)),
    ([1.0, 2.0], Point3(1.0, 2.0, 0.0)),
])
def test_to_point3_pads_a_short_iterable_with_zeros(short, expected):
    assert to_point3(short) == expected


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
    b = TraceBuilder()
    b.add_segment((0.0, 0.0), (1.0, 1.0))
    b.add_segment((2.0, 2.0), (3.0, 3.0))

    ax = Figure().subplots()
    lines = ax.plot(*b.xy, color="black")
    assert len(lines) == 1


# ---------------------------------------------------------------------------
# Color maps
# ---------------------------------------------------------------------------

def test_gom_color_map_has_eight_discrete_colors():
    assert GomColorMap().N == 8


def test_gom_color_map_sets_over_and_under_colors():
    from matplotlib.colors import to_rgba

    cmap = GomColorMap()
    assert tuple(cmap.get_over()) == pytest.approx(to_rgba("darkred"), abs=TOL)
    assert tuple(cmap.get_under()) == pytest.approx(to_rgba("magenta"), abs=TOL)


def test_gom_cmap_is_a_ready_made_instance():
    assert isinstance(GOM_CMAP, GomColorMap)


# ---------------------------------------------------------------------------
# AxesHelper: construction and configuration
# ---------------------------------------------------------------------------

def test_helper_wraps_the_axes_and_exposes_it():
    ax = Figure().subplots()
    assert AxesHelper(ax).ax is ax


def test_helper_enforces_an_equal_aspect_ratio_by_default():
    ax = Figure().subplots()
    AxesHelper(ax)
    assert ax.get_aspect() == 1.0


def test_skip_aspect_leaves_the_aspect_ratio_alone():
    ax = Figure().subplots()
    before = ax.get_aspect()
    AxesHelper(ax, skip_aspect=True)
    assert ax.get_aspect() == before


def test_hide_axes_turns_the_axis_off():
    ax = Figure().subplots()
    AxesHelper(ax, hide_axes=True)
    assert not ax.axison


def test_set_bounds_applies_the_box_to_both_limits():
    helper = new_helper()
    helper.set_bounds(Aabb2(x_min=-1.0, x_max=2.0, y_min=-3.0, y_max=4.0))
    assert helper.ax.get_xlim() == pytest.approx((-1.0, 2.0), abs=TOL)
    assert helper.ax.get_ylim() == pytest.approx((-3.0, 4.0), abs=TOL)


def test_viewport_returns_a_viewport_bound_to_the_helper():
    helper = new_helper()
    view = helper.viewport(Iso3.identity())
    assert isinstance(view, ViewPort3)
    assert view.helper is helper


# ---------------------------------------------------------------------------
# AxesHelper: 2D draw methods
# ---------------------------------------------------------------------------

def test_draw_curve_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_curve(sample_curve(), color="black")
    assert artist_count(helper.ax) > before


def test_fill_curve_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.fill_curve(sample_curve(), alpha=0.2)
    assert artist_count(helper.ax) > before


def test_draw_circle_adds_one_patch_per_circle():
    helper = new_helper()
    helper.draw_circle(Circle2(0.0, 0.0, 1.0), Circle2(2.0, 0.0, 1.0))
    assert len(helper.ax.patches) == 2


def test_draw_circle_accepts_a_raw_xyr_iterable():
    helper = new_helper()
    helper.draw_circle((0.0, 0.0, 1.0))
    assert len(helper.ax.patches) == 1


@pytest.mark.parametrize("bad", [(0.0, 0.0), (0.0, 0.0, 1.0, 9.0), ()])
def test_draw_circle_rejects_an_iterable_that_is_not_three_values(bad):
    """ Regression: a 4-value input silently dropped the extras, and a 2-value one raised a bare
    unpacking error that never mentioned circles. """
    helper = new_helper()
    with pytest.raises(ValueError, match="x, y, r"):
        helper.draw_circle(bad)


# ---------------------------------------------------------------------------
# AxesHelper: zero-valued arguments must not be mistaken for "not supplied"
# ---------------------------------------------------------------------------

def test_draw_boundary_rejects_a_non_positive_tolerance():
    """ Regression: `tol or default` silently replaced an explicit 0.0 with the computed default. """
    helper = new_helper()
    with pytest.raises(ValueError, match="positive"):
        helper.draw_boundary(sample_boundary(), tol=0.0)


def test_draw_distance_honors_a_zero_label_offset():
    """
    Regression: `label_offset or default` treated an explicit 0.0 as "not supplied", so a label
    pinned to the leader end silently jumped to the default standoff.
    """
    helper = new_helper()
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), label_offset=0.0)
    zero_offset = {t.get_text(): t.xy for t in helper.ax.texts}["3.000"]

    helper = new_helper()
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)))
    defaulted = {t.get_text(): t.xy for t in helper.ax.texts}["3.000"]

    assert zero_offset != defaulted


def test_draw_normals_omits_color_when_none_so_the_default_applies():
    helper = new_helper()
    helper.draw_normals(sample_boundary(), count=3, length=0.1, color=None)
    assert artist_count(helper.ax) == 3


def test_draw_circle_honors_the_fill_flag():
    helper = new_helper()
    helper.draw_circle(Circle2(0.0, 0.0, 1.0), fill=True)
    assert helper.ax.patches[0].get_fill()


def test_draw_arc_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_arc(Arc2(0.0, 0.0, 1.0, 0.0, numpy.pi / 2))
    assert artist_count(helper.ax) > before


def test_draw_arc_handles_a_negative_sweep():
    helper = new_helper()
    helper.draw_arc(Arc2(0.0, 0.0, 1.0, numpy.pi, -numpy.pi / 2))
    assert len(helper.ax.patches) == 1


def test_draw_segment_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_segment(Segment2(0.0, 0.0, 1.0, 1.0))
    assert artist_count(helper.ax) > before


def test_draw_boundary_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_boundary(sample_boundary(), color="black")
    assert artist_count(helper.ax) > before


def test_draw_boundary_accepts_an_explicit_tolerance():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_boundary(sample_boundary(), tol=0.001)
    assert artist_count(helper.ax) > before


def test_draw_aabb_adds_one_patch_per_box():
    helper = new_helper()
    helper.draw_aabb(Aabb2(0.0, 0.0, 1.0, 2.0), Aabb2(3.0, 3.0, 4.0, 4.0))
    assert len(helper.ax.patches) == 2


def test_draw_aabb_places_the_rectangle_on_the_box():
    helper = new_helper()
    patch = helper.draw_aabb(Aabb2(-1.0, 2.0, 4.0, 8.0))[0]
    assert patch.get_xy() == pytest.approx((-1.0, 2.0), abs=TOL)
    assert patch.get_width() == pytest.approx(5.0, abs=TOL)
    assert patch.get_height() == pytest.approx(6.0, abs=TOL)


def test_draw_aabb_honors_the_fill_flag():
    helper = new_helper()
    helper.draw_aabb(Aabb2(0.0, 0.0, 1.0, 1.0), fill=True)
    assert helper.ax.patches[0].get_fill()


def test_draw_line_adds_one_artist_per_line():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_line(Line2.x_axis(), Line2.y_axis())
    assert artist_count(helper.ax) == before + 2


def test_draw_line_spans_the_axes_when_no_extent_is_given():
    """ An infinite line is drawn as an AxLine, which re-clips itself as the view changes. """
    helper = new_helper()
    artist = helper.draw_line(Line2(0.0, 0.0, 1.0, 1.0))[0]
    assert isinstance(artist, AxLine)


def test_draw_line_cut_to_an_extent_lands_on_the_parameter_endpoints():
    helper = new_helper()
    line = Line2.new_normalize(1.0, 1.0, 1.0, 0.0)
    artist = helper.draw_line(line, t=(2.0, 5.0))[0]
    assert list(artist.get_xdata()) == pytest.approx([3.0, 6.0], abs=TOL)
    assert list(artist.get_ydata()) == pytest.approx([1.0, 1.0], abs=TOL)


def test_draw_line_cut_to_an_extent_is_a_plain_line_not_an_axline():
    """ A finite piece is real geometry, so it autoscales like anything else drawn. """
    helper = new_helper()
    artist = helper.draw_line(Line2.x_axis(), t=(0.0, 1.0))[0]
    assert isinstance(artist, Line2D) and not isinstance(artist, AxLine)


@pytest.mark.parametrize("bad", [(0.0,), (0.0, 1.0, 2.0), ()])
def test_draw_line_rejects_an_extent_that_is_not_two_values(bad):
    helper = new_helper()
    with pytest.raises(ValueError, match="two values"):
        helper.draw_line(Line2.x_axis(), t=bad)


def test_draw_line_honors_a_zero_valued_extent_endpoint():
    """
    Regression guard on the same falsy-zero shape found elsewhere: `t=(0.0, 5.0)` must start at the
    origin, not be treated as "not supplied".
    """
    helper = new_helper()
    artist = helper.draw_line(Line2.new_normalize(0.0, 0.0, 1.0, 0.0), t=(0.0, 5.0))[0]
    assert list(artist.get_xdata()) == pytest.approx([0.0, 5.0], abs=TOL)


def test_an_infinite_line_does_not_drag_the_autoscaled_view():
    """
    `Axes.axline` registers its two reference points in the data limits, so a construction line
    through a far-off origin would otherwise rescale a plot of distant geometry. An infinite line is
    an annotation and must not influence autoscaling.
    """
    ax = Figure().subplots()
    helper = AxesHelper(ax)
    helper.draw_curve(Curve2(numpy.array([[100.0, 100.0], [200.0, 200.0]]), tol=1.0e-6))
    ax.autoscale_view()
    before = ax.get_xlim()

    helper.draw_line(Line2(0.0, 0.0, 1.0, 1.0))
    ax.autoscale_view()
    assert ax.get_xlim() == pytest.approx(before, abs=TOL)


def test_a_finite_line_piece_does_take_part_in_autoscaling():
    ax = Figure().subplots()
    helper = AxesHelper(ax)
    helper.draw_curve(Curve2(numpy.array([[0.0, 0.0], [1.0, 1.0]]), tol=1.0e-6))
    ax.autoscale_view()
    before = ax.get_xlim()

    helper.draw_line(Line2.x_axis(), t=(0.0, 50.0))
    ax.autoscale_view()
    assert ax.get_xlim()[1] > before[1]


def test_draw_normals_adds_one_arrow_per_sample():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_normals(sample_boundary(), count=5, length=0.1, color="red")
    assert artist_count(helper.ax) == before + 5


@pytest.mark.parametrize("source", [sample_curve(), sample_boundary()])
def test_draw_normals_accepts_both_curves_and_boundaries(source):
    """
    `Curve2` and `Boundary2` are unrelated types, but both are arc-length parameterized and both
    report a surface point with a normal, which is the entire contract this method needs.
    """
    helper = new_helper()
    arrows = helper.draw_normals(source, count=4, length=0.1)
    assert len(arrows) == 4


def test_draw_normals_accepts_a_mix_of_sources_in_one_call():
    helper = new_helper()
    arrows = helper.draw_normals(sample_curve(), sample_boundary(), count=3, length=0.1)
    assert len(arrows) == 6


def test_draw_normals_points_along_the_surface_normal():
    """ The arrow must run from the sampled point to that point offset along its own normal. """
    helper = new_helper()
    source = sample_boundary()
    arrows = helper.draw_normals(source, count=3, length=0.5)

    for arrow, t in zip(arrows, numpy.linspace(0, source.length(), 3)):
        station = source.at_length(t)
        expected_tail = (station.point.x, station.point.y)
        expected_tip = tuple(station.surface_point.at_distance(0.5))
        assert arrow.xyann == pytest.approx(expected_tail, abs=1e-9)
        assert tuple(arrow.xy) == pytest.approx(expected_tip, abs=1e-9)


@pytest.mark.parametrize("bad", ["not geometry", Point2(0.0, 0.0), Circle2(0.0, 0.0, 1.0)])
def test_draw_normals_rejects_a_source_that_is_neither(bad):
    helper = new_helper()
    with pytest.raises(TypeError, match="Curve2 or Boundary2"):
        helper.draw_normals(bad, count=3, length=0.1)


def test_draw_spline_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_spline(CubicSpline2(0.0, 0.0, 1.0, 2.0, 2.0, -1.0, 3.0, 0.0), tol=0.01)
    assert artist_count(helper.ax) > before


def test_draw_point_adds_an_artist():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_point(Point2(0.0, 0.0), (1.0, 1.0), [2.0, 2.0])
    assert artist_count(helper.ax) > before


def test_draw_surface_point_adds_markers_and_one_arrow_each():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_surface_point(SurfacePoint2(0.0, 0.0, 0.0, 1.0), SurfacePoint2(1.0, 0.0, 0.0, 1.0))
    # One `plot` call for the markers, plus one annotation arrow per point.
    assert artist_count(helper.ax) == before + 3


def test_draw_text_adds_and_returns_the_annotation():
    helper = new_helper()
    result = helper.draw_text("hello", Point2(0.0, 0.0))
    assert result is not None
    assert helper.ax.texts[-1].get_text() == "hello"


def test_draw_text_applies_the_shift_to_the_position():
    helper = new_helper()
    plain = helper.draw_text("a", (1.0, 1.0))
    shifted = helper.draw_text("b", (1.0, 1.0), shift=(2.0, 3.0))
    assert shifted.xy == pytest.approx((plain.xy[0] + 2.0, plain.xy[1] + 3.0), abs=TOL)


def test_draw_arrow_adds_and_returns_the_annotation():
    helper = new_helper()
    before = artist_count(helper.ax)
    result = helper.draw_arrow(Point2(0.0, 0.0), Point2(1.0, 1.0))
    assert result is not None
    assert artist_count(helper.ax) > before


def test_draw_labeled_arrow_adds_the_arrow_and_the_label():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_labeled_arrow(Point2(0.0, 0.0), Point2(2.0, 0.0), "label")
    assert artist_count(helper.ax) == before + 2
    assert helper.ax.texts[-1].get_text() == "label"


# ---------------------------------------------------------------------------
# AxesHelper: distance and label placement
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label_place", LABEL_PLACES)
def test_draw_distance_adds_artists_for_every_label_placement(label_place):
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), label_place=label_place)
    assert artist_count(helper.ax) > before


def test_label_place_tokens_match_the_literal_alias():
    """ The valid-token tuple is derived from the alias, so the two can never disagree. """
    assert set(LABEL_PLACES) == set(get_args(LabelPlace))


@pytest.mark.parametrize("bad", ["Outside", "OUTSIDE", "outside_reversed", "", "middle"])
def test_draw_distance_rejects_an_unknown_label_placement(bad):
    """
    Regression: an out-of-range placement used to fall through every branch and surface as an
    UnboundLocalError on `label_coords`, rather than naming the bad argument.
    """
    helper = new_helper()
    with pytest.raises(ValueError) as excinfo:
        helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), label_place=bad)

    message = str(excinfo.value)
    assert "label_place" in message
    # The message must name the alternatives, not just reject the input.
    for token in LABEL_PLACES:
        assert token in message


def test_draw_distance_validates_before_drawing_anything():
    """ A rejected placement must not leave a half-drawn measurement on the axes. """
    helper = new_helper()
    before = artist_count(helper.ax)
    with pytest.raises(ValueError):
        helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), label_place="nope")
    assert artist_count(helper.ax) == before


def test_draw_distance_formats_the_value_with_the_template():
    helper = new_helper()
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), template="{value:.1f} mm")
    assert "3.0 mm" in {t.get_text() for t in helper.ax.texts}


def test_draw_distance_applies_the_value_scale_to_the_label_only():
    helper = new_helper()
    d = Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0))
    helper.draw_distance(d, template="{value:.1f}", scale_value=10.0)
    assert "30.0" in {t.get_text() for t in helper.ax.texts}
    assert d.value == pytest.approx(3.0, abs=TOL)


def test_draw_distance_accepts_a_side_shift():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), side_shift=1.0)
    assert artist_count(helper.ax) > before


def test_draw_distance_of_a_negative_value_still_draws():
    helper = new_helper()
    d = Distance2(Point2(3.0, 0.0), Point2(0.0, 0.0), Vector2(1.0, 0.0))
    before = artist_count(helper.ax)
    helper.draw_distance(d)
    assert artist_count(helper.ax) > before


# ---------------------------------------------------------------------------
# ViewPort3: 3D entities in parallel projection
# ---------------------------------------------------------------------------

def new_viewport() -> ViewPort3:
    return new_helper().viewport(Iso3.identity())


def test_viewport_draw_line_adds_an_artist():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_line(Line3(0.0, 0.0, 0.0, 1.0, 1.0, 0.0), t=2.0)
    assert artist_count(view.helper.ax) > before


def test_viewport_draw_line_accepts_an_explicit_start():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_line(Line3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0), t=2.0, t0=1.0)
    assert artist_count(view.helper.ax) > before


def test_viewport_draw_line_honors_a_zero_start():
    """
    Regression: `t0 or -t` treated an explicit t0=0.0 as "not supplied", so a ray asked to start at
    the line origin was silently drawn from -t instead, doubling its length in the wrong direction.
    """
    line = Line3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0)

    view = new_viewport()
    view.draw_line(line, t=2.0, t0=0.0)
    from_origin = view.helper.ax.lines[-1].get_xdata()

    view = new_viewport()
    view.draw_line(line, t=2.0)
    symmetric = view.helper.ax.lines[-1].get_xdata()

    assert tuple(from_origin) == pytest.approx((0.0, 2.0), abs=TOL)
    assert tuple(symmetric) == pytest.approx((-2.0, 2.0), abs=TOL)


def test_viewport_draw_curve_projects_the_curve_vertices():
    view = new_viewport()
    curve = Curve3(numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 5.0], [3.0, 1.0, 9.0]]), tol=1.0e-6)
    line = view.draw_curve(curve)[0]
    # Under the identity view, projecting is just dropping z.
    assert list(line.get_xdata()) == pytest.approx([0.0, 1.0, 3.0], abs=TOL)
    assert list(line.get_ydata()) == pytest.approx([0.0, 2.0, 1.0], abs=TOL)


def test_viewport_draw_point_cloud_draws_markers_without_connecting_them():
    view = new_viewport()
    line = view.draw_point_cloud(PointCloud(numpy.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])))[0]
    assert line.get_linestyle() == "None"
    assert len(line.get_xdata()) == 2


def test_viewport_draw_aabb_draws_all_twelve_edges():
    """ A box has twelve edges, each contributing two points and a break to the single artist. """
    view = new_viewport()
    line = view.draw_aabb(Aabb3(0.0, 0.0, 0.0, 1.0, 2.0, 3.0))[0]
    xs = list(line.get_xdata())
    assert len(xs) == 36
    assert sum(1 for x in xs if x is None or numpy.isnan(x)) == 12


def test_viewport_draw_plane_places_a_quad_of_the_requested_size():
    view = new_viewport()
    patch = view.draw_plane(Plane3.from_point_normal(0.0, 0.0, 4.0, 0.0, 0.0, 1.0),
                            center=Point3(0.0, 0.0, 0.0), size=2.0)[0]
    corners = patch.get_xy()[:4]
    # Seen face-on, the quad projects to a 2x2 square centered where the anchor projects onto it.
    assert corners[:, 0].min() == pytest.approx(-1.0, abs=TOL)
    assert corners[:, 0].max() == pytest.approx(1.0, abs=TOL)
    assert corners[:, 1].min() == pytest.approx(-1.0, abs=TOL)
    assert corners[:, 1].max() == pytest.approx(1.0, abs=TOL)


def test_viewport_draw_plane_rejects_a_non_positive_size():
    view = new_viewport()
    with pytest.raises(ValueError, match="positive"):
        view.draw_plane(Plane3.from_point_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
                        center=Point3(0.0, 0.0, 0.0), size=0.0)


def test_viewport_draw_circle_seen_face_on_is_a_true_circle():
    view = new_viewport()
    patch = view.draw_circle(Circle3(1.0, 2.0, 0.0, 0.0, 0.0, 1.0, 3.0))[0]
    assert patch.get_center() == pytest.approx((1.0, 2.0), abs=TOL)
    assert patch.get_width() == pytest.approx(6.0, abs=TOL)
    assert patch.get_height() == pytest.approx(6.0, abs=TOL)


def test_viewport_draw_circle_foreshortens_by_the_cosine_of_the_tilt():
    """
    The projected minor axis of a tilted circle is exactly r*cos(tilt), which is the property the
    ellipse construction exists to get right.
    """
    tilt = numpy.deg2rad(60.0)
    view = new_viewport()
    patch = view.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, numpy.sin(tilt), numpy.cos(tilt), 2.0))[0]
    assert patch.get_width() == pytest.approx(4.0, abs=1.0e-9)
    assert patch.get_height() == pytest.approx(4.0 * numpy.cos(tilt), abs=1.0e-9)


def test_viewport_draw_circle_seen_edge_on_collapses_to_a_line():
    view = new_viewport()
    patch = view.draw_circle(Circle3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0))[0]
    assert patch.get_width() == pytest.approx(6.0, abs=TOL)
    assert patch.get_height() == pytest.approx(0.0, abs=TOL)


def test_viewport_draw_distance_labels_the_true_3d_value_not_the_projected_one():
    """
    Projection foreshortens the arrow, but the number beside it has to stay the real measurement. A
    dimension whose label disagreed with the geometry would be silently wrong rather than merely
    ugly.
    """
    view = new_helper().viewport(Iso3.from_rotation(numpy.deg2rad(30.0), 1.0, 0.0, 0.0))
    view.draw_distance(Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 10.0)))
    labels = [t.get_text() for t in view.helper.ax.texts if t.get_text()]
    # The measurement projects to half its length in this view, so an uncorrected label would read
    # "5.000".
    assert labels == ["10.000"]


def test_viewport_draw_distance_still_honors_a_unit_scale_factor():
    view = new_helper().viewport(Iso3.from_rotation(numpy.deg2rad(30.0), 1.0, 0.0, 0.0))
    view.draw_distance(Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 10.0)), scale_value=25.4)
    labels = [t.get_text() for t in view.helper.ax.texts if t.get_text()]
    assert labels == ["254.000"]


def test_viewport_draw_distance_draws_a_genuinely_zero_length_measurement():
    """
    A zero-length measurement projects to zero too, but it is not a foreshortening failure and must
    not be rejected as one. It also must not divide by zero recovering the true value.
    """
    view = new_viewport()
    view.draw_distance(Distance3(Point3(1.0, 1.0, 0.0), Point3(1.0, 1.0, 0.0),
                                 Vector3(1.0, 0.0, 0.0)))
    labels = [t.get_text() for t in view.helper.ax.texts if t.get_text()]
    assert labels == ["0.000"]


def test_viewport_draw_distance_rejects_a_measurement_along_the_view_direction():
    view = new_viewport()
    with pytest.raises(ValueError, match="view direction"):
        view.draw_distance(Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 5.0)))


def test_viewport_draw_distance_passes_a_bad_label_place_through_to_validation():
    view = new_viewport()
    with pytest.raises(ValueError, match="invalid label_place"):
        view.draw_distance(Distance3(Point3(0.0, 0.0, 0.0), Point3(3.0, 0.0, 0.0)),
                           label_place="sideways")


def test_plane_basis_returns_an_orthonormal_pair_for_any_normal():
    for normal in [Vector3(0.0, 0.0, 1.0), Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0),
                   Vector3(1.0, 1.0, 1.0), Vector3(-3.0, 0.2, 0.0)]:
        u, v = _plane_basis(normal)
        n = normal.normalized()
        assert u.norm() == pytest.approx(1.0, abs=1.0e-12)
        assert v.norm() == pytest.approx(1.0, abs=1.0e-12)
        assert u.dot(v) == pytest.approx(0.0, abs=1.0e-12)
        assert u.dot(n) == pytest.approx(0.0, abs=1.0e-12)
        assert v.dot(n) == pytest.approx(0.0, abs=1.0e-12)


def test_viewport_draw_mesh_outline_honors_an_empty_kwargs_dict():
    """ Regression: `visible_kwargs or default` replaced an explicit empty dict with the default. """
    view = new_viewport()
    view.draw_mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0), visible_kwargs={}, no_hidden=True)
    # An empty dict means "use Matplotlib's own defaults", which are not the helper's black.
    assert view.helper.ax.lines[-1].get_color() != "black"


def test_viewport_draw_coordinate_system_adds_the_visible_axes():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_coordinate_system(Iso3.identity(), 1.0)
    # Looking down Z, the X and Y axes are visible and the Z axis is hidden; each visible axis
    # contributes an arrow and a label.
    assert artist_count(view.helper.ax) == before + 4


def test_viewport_draw_labeled_point_adds_an_artist():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_labeled_point(Point3(0.0, 0.0, 0.0), "P")
    assert artist_count(view.helper.ax) > before


def test_viewport_draw_labeled_point_can_add_a_leader_arrow():
    view = new_viewport()
    plain_before = artist_count(view.helper.ax)
    view.draw_labeled_point(Point3(0.0, 0.0, 0.0), "P")
    plain = artist_count(view.helper.ax) - plain_before

    view = new_viewport()
    arrow_before = artist_count(view.helper.ax)
    view.draw_labeled_point(Point3(0.0, 0.0, 0.0), "P", offset_2d=(1.0, 1.0), arrow=True)
    assert artist_count(view.helper.ax) - arrow_before == plain + 1


def test_viewport_draw_dimension_arrow_adds_an_artist():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0")
    assert artist_count(view.helper.ax) > before


def test_viewport_draw_dimension_arrow_adds_leaders_when_shifted():
    view = new_viewport()
    plain_before = artist_count(view.helper.ax)
    view.draw_dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0")
    plain = artist_count(view.helper.ax) - plain_before

    view = new_viewport()
    shifted_before = artist_count(view.helper.ax)
    view.draw_dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0", leader_shift=(0.0, 1.0, 0.0))
    assert artist_count(view.helper.ax) - shifted_before == plain + 2


def test_viewport_draw_mesh_outline_adds_an_artist():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0))
    assert artist_count(view.helper.ax) > before


def test_viewport_draw_mesh_outline_can_suppress_hidden_edges():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.draw_mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0), no_hidden=True)
    # Only the visible trace is plotted, rather than a visible and a hidden trace.
    assert artist_count(view.helper.ax) == before + 1


def test_viewport_find_mesh_edge_point_returns_a_point_on_the_mesh():
    view = new_viewport()
    mesh = Mesh3.create_box(2.0, 2.0, 2.0)
    found = view.find_mesh_edge_point(1.0, 0.0, mesh)
    assert isinstance(found, Point3)
    # The box spans -1 to 1 on each axis, and the query asks for the +X extreme.
    assert found.x == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Varargs and returned artists
# ---------------------------------------------------------------------------

def two_of(factory):
    return [factory(), factory()]


VARARGS_CALLS = {
    "draw_aabb": lambda h: h.draw_aabb(Aabb2(0.0, 0.0, 1.0, 1.0), Aabb2(2.0, 2.0, 3.0, 3.0)),
    "draw_arc": lambda h: h.draw_arc(Arc2(0.0, 0.0, 1.0, 0.0, 1.0), Arc2(3.0, 0.0, 1.0, 0.0, 1.0)),
    "draw_circle": lambda h: h.draw_circle(Circle2(0.0, 0.0, 1.0), Circle2(3.0, 0.0, 1.0)),
    "draw_line": lambda h: h.draw_line(Line2.x_axis(), Line2.y_axis()),
    "draw_line_extent": lambda h: h.draw_line(Line2.x_axis(), Line2.y_axis(), t=(0.0, 1.0)),
    "draw_curve": lambda h: h.draw_curve(sample_curve(), sample_curve()),
    "draw_segment": lambda h: h.draw_segment(Segment2(0.0, 0.0, 1.0, 1.0), Segment2(2.0, 2.0, 3.0, 3.0)),
    "draw_boundary": lambda h: h.draw_boundary(sample_boundary(), sample_boundary()),
    "fill_curve": lambda h: h.fill_curve(sample_curve(), sample_curve()),
    "draw_spline": lambda h: h.draw_spline(
        CubicSpline2(0.0, 0.0, 1.0, 2.0, 2.0, -1.0, 3.0, 0.0),
        CubicSpline2(0.0, 5.0, 1.0, 7.0, 2.0, 4.0, 3.0, 5.0),
        tol=0.01,
    ),
}


@pytest.mark.parametrize("name", sorted(VARARGS_CALLS))
def test_entity_draw_methods_accept_varargs_and_return_one_artist_each(name):
    helper = new_helper()
    result = VARARGS_CALLS[name](helper)
    assert isinstance(result, list)
    assert len(result) == 2


# Methods whose no-entity call still needs an argument, or whose entry above is a second calling
# form of a method already keyed by its own name.
EMPTY_CALLS = {
    "draw_spline": lambda h: h.draw_spline(tol=0.01),
    "draw_line_extent": lambda h: h.draw_line(t=(0.0, 1.0)),
}


@pytest.mark.parametrize("name", sorted(VARARGS_CALLS))
def test_entity_draw_methods_return_an_empty_list_when_given_nothing(name):
    """ Drawing a computed and possibly empty collection should not need a special case. """
    helper = new_helper()
    call = EMPTY_CALLS.get(name, lambda h, n=name: getattr(h, n)())
    assert call(helper) == []


def test_draw_point_returns_a_single_artist_because_it_makes_one():
    helper = new_helper()
    result = helper.draw_point(Point2(0.0, 0.0), Point2(1.0, 1.0))
    assert isinstance(result, Line2D)


def test_draw_point_accepts_a_single_point_array():
    """
    Regression: an (n, 2) array used to be treated as one coordinate, unpacking its first two rows
    as x and y. That silently drew a single wrong point instead of raising, and an array is the form
    the rest of the library hands point sets back in.
    """
    points = numpy.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
    artist = new_helper().draw_point(points)
    assert list(artist.get_xdata()) == pytest.approx([0.0, 2.0, 4.0], abs=TOL)
    assert list(artist.get_ydata()) == pytest.approx([1.0, 3.0, 5.0], abs=TOL)


def test_draw_point_of_an_array_matches_the_varargs_form():
    points = numpy.array([[0.0, 1.0], [2.0, 3.0]])
    as_array = new_helper().draw_point(points)
    as_varargs = new_helper().draw_point(*points)
    assert list(as_array.get_xdata()) == pytest.approx(list(as_varargs.get_xdata()), abs=TOL)
    assert list(as_array.get_ydata()) == pytest.approx(list(as_varargs.get_ydata()), abs=TOL)


def test_draw_point_ignores_the_extra_columns_of_a_wider_array():
    points = numpy.array([[0.0, 1.0, 9.0], [2.0, 3.0, 9.0]])
    artist = new_helper().draw_point(points)
    assert list(artist.get_xdata()) == pytest.approx([0.0, 2.0], abs=TOL)


def test_draw_point_accepts_a_nested_sequence_of_points():
    artist = new_helper().draw_point([[0.0, 1.0], [2.0, 3.0]])
    assert list(artist.get_xdata()) == pytest.approx([0.0, 2.0], abs=TOL)


def test_draw_point_rejects_a_point_array_that_is_too_narrow():
    with pytest.raises(ValueError, match="at least two columns"):
        new_helper().draw_point(numpy.array([[0.0], [1.0]]))


@pytest.mark.parametrize("point", [Point2(1.0, 2.0), Vector2(1.0, 2.0), (1.0, 2.0), [1.0, 2.0],
                                   numpy.array([1.0, 2.0])])
def test_draw_point_still_takes_a_single_loose_coordinate(point):
    """ The array form must not swallow the one-point case, which reports a lower dimension. """
    artist = new_helper().draw_point(point)
    assert list(artist.get_xdata()) == pytest.approx([1.0], abs=TOL)
    assert list(artist.get_ydata()) == pytest.approx([2.0], abs=TOL)


def test_draw_point_of_nothing_still_returns_an_empty_artist():
    helper = new_helper()
    result = helper.draw_point()
    assert isinstance(result, Line2D)
    assert len(result.get_xdata()) == 0


@pytest.mark.parametrize("call", [
    lambda h: h.draw_text("t", Point2(0.0, 0.0)),
    lambda h: h.draw_arrow(Point2(0.0, 0.0), Point2(1.0, 1.0)),
])
def test_coordinate_draw_methods_return_a_single_annotation(call):
    assert isinstance(call(new_helper()), Annotation)


COMPOSITE_CALLS = {
    "draw_distance": lambda h: h.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0))),
    "draw_labeled_arrow": lambda h: h.draw_labeled_arrow(Point2(0.0, 0.0), Point2(2.0, 0.0), "x"),
    "draw_surface_point": lambda h: h.draw_surface_point(SurfacePoint2(0.0, 0.0, 0.0, 1.0)),
    "draw_normals": lambda h: h.draw_normals(sample_boundary(), count=3, length=0.1),
}


@pytest.mark.parametrize("name", sorted(COMPOSITE_CALLS))
def test_composite_draw_methods_return_every_artist_they_added(name):
    helper = new_helper()
    before = artist_count(helper.ax)
    result = COMPOSITE_CALLS[name](helper)
    assert isinstance(result, list)
    assert len(result) == artist_count(helper.ax) - before


ALL_RETURNING_CALLS = {**VARARGS_CALLS, **COMPOSITE_CALLS}


@pytest.mark.parametrize("name", sorted(ALL_RETURNING_CALLS))
def test_returned_artists_are_really_attached_to_the_axes(name):
    """
    A returned artist is only useful if it is the one on the plot, so that restyling it through the
    return value actually changes the rendered figure.
    """
    helper = new_helper()
    on_axes = set()
    for artist in ALL_RETURNING_CALLS[name](helper):
        assert artist.axes is helper.ax
        on_axes.add(id(artist))
    assert len(on_axes) == len(ALL_RETURNING_CALLS[name](new_helper()))


def test_viewport_draw_line_returns_the_line_artist():
    view = new_viewport()
    assert isinstance(view.draw_line(Line3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0)), Line2D)


@pytest.mark.parametrize("call", [
    lambda v: v.draw_coordinate_system(Iso3.identity(), 1.0),
    lambda v: v.draw_labeled_point(Point3(0.0, 0.0, 0.0), "P"),
    lambda v: v.draw_dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0"),
    lambda v: v.draw_mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0)),
])
def test_viewport_composites_return_every_artist_they_added(call):
    view = new_viewport()
    before = artist_count(view.helper.ax)
    result = call(view)
    assert isinstance(result, list)
    assert len(result) == artist_count(view.helper.ax) - before


# ---------------------------------------------------------------------------
# Named styling arguments
# ---------------------------------------------------------------------------

def test_merge_style_drops_the_arguments_left_as_none():
    assert merge_style({}, color="red", linewidth=None, alpha=0.5) == {"color": "red", "alpha": 0.5}


def test_merge_style_leaves_the_open_kwargs_intact():
    assert merge_style({"solid_capstyle": "round"}, color=None) == {"solid_capstyle": "round"}


def test_omitting_color_leaves_the_axes_color_cycle_in_charge():
    """
    The named styling arguments must not be forwarded when left as None. Passing `color=None` to
    Matplotlib is not the same as omitting it: omitting it draws the next color from the cycle.
    """
    helper = new_helper()
    first = helper.draw_curve(sample_curve())[0]
    second = helper.draw_curve(sample_curve())[0]
    assert first.get_color() != second.get_color()


@pytest.mark.parametrize("name,call,check", [
    ("draw_curve", lambda h: h.draw_curve(sample_curve(), color="red")[0],
     lambda a: a.get_color() == "red"),
    ("draw_segment", lambda h: h.draw_segment(Segment2(0.0, 0.0, 1.0, 1.0), linewidth=4.0)[0],
     lambda a: a.get_linewidth() == 4.0),
    ("draw_spline", lambda h: h.draw_spline(
        CubicSpline2(0.0, 0.0, 1.0, 2.0, 2.0, -1.0, 3.0, 0.0), tol=0.01, linestyle="--")[0],
     lambda a: a.get_linestyle() == "--"),
    ("draw_boundary", lambda h: h.draw_boundary(sample_boundary(), alpha=0.25)[0],
     lambda a: a.get_alpha() == 0.25),
    ("draw_circle", lambda h: h.draw_circle(Circle2(0.0, 0.0, 1.0), linewidth=3.0)[0],
     lambda a: a.get_linewidth() == 3.0),
    ("draw_arc", lambda h: h.draw_arc(Arc2(0.0, 0.0, 1.0, 0.0, 1.0), alpha=0.5)[0],
     lambda a: a.get_alpha() == 0.5),
    ("draw_aabb", lambda h: h.draw_aabb(Aabb2(0.0, 0.0, 1.0, 1.0), linestyle="--")[0],
     lambda a: a.get_linestyle() == "--"),
    ("draw_line", lambda h: h.draw_line(Line2.x_axis(), color="purple")[0],
     lambda a: a.get_color() == "purple"),
    ("draw_line_extent", lambda h: h.draw_line(Line2.x_axis(), t=(0.0, 1.0), linewidth=6.0)[0],
     lambda a: a.get_linewidth() == 6.0),
    ("fill_curve", lambda h: h.fill_curve(sample_curve(), alpha=0.3)[0],
     lambda a: a.get_alpha() == 0.3),
    ("draw_point", lambda h: h.draw_point(Point2(0.0, 0.0), color="green"),
     lambda a: a.get_color() == "green"),
    ("draw_text", lambda h: h.draw_text("x", Point2(0.0, 0.0), fontsize=22),
     lambda a: a.get_fontsize() == 22),
])
def test_named_styling_arguments_reach_the_artist(name, call, check):
    assert check(call(new_helper()))


@pytest.mark.parametrize("name,call", [
    ("draw_curve", lambda h: h.draw_curve(sample_curve(), solid_capstyle="round")[0]),
    ("draw_segment", lambda h: h.draw_segment(Segment2(0.0, 0.0, 1.0, 1.0), dash_capstyle="round")[0]),
    ("draw_circle", lambda h: h.draw_circle(Circle2(0.0, 0.0, 1.0), hatch="//")[0]),
    ("draw_aabb", lambda h: h.draw_aabb(Aabb2(0.0, 0.0, 1.0, 1.0), hatch="//")[0]),
    ("draw_line", lambda h: h.draw_line(Line2.x_axis(), dash_capstyle="round")[0]),
    ("draw_text", lambda h: h.draw_text("x", Point2(0.0, 0.0), rotation=45)),
])
def test_uncommon_matplotlib_arguments_still_pass_through(name, call):
    """
    Naming the common styling arguments must not close the set. Anything else Matplotlib accepts
    has to keep working, which is why **kwargs stays open and untyped.
    """
    assert call(new_helper()) is not None


def test_labels_reach_the_legend():
    helper = new_helper()
    helper.draw_curve(sample_curve(), label="upper")
    helper.draw_curve(sample_curve(), label="lower")
    legend_texts = [t.get_text() for t in helper.ax.legend().get_texts()]
    assert legend_texts == ["upper", "lower"]


# ---------------------------------------------------------------------------
# Composite element styling
# ---------------------------------------------------------------------------

def test_element_style_suppresses_on_false():
    assert element_style(False, {"color": "red"}) is None


@pytest.mark.parametrize("value", [None, True])
def test_element_style_accepts_the_defaults(value):
    assert element_style(value, {"color": "red"}) == {"color": "red"}


def test_element_style_merges_over_the_defaults_instead_of_replacing_them():
    """
    Restyling one thing about an element must not silently discard the rest of its appearance,
    which is what replacing the defaults wholesale would do.
    """
    resolved = element_style({"color": "blue"}, {"color": "red", "linestyle": "--"})
    assert resolved == {"color": "blue", "linestyle": "--"}


def test_element_style_does_not_mutate_the_defaults():
    defaults = {"color": "red"}
    element_style({"color": "blue"}, defaults)
    element_style(None, defaults)["color"] = "green"
    assert defaults == {"color": "red"}


def test_an_empty_style_dict_is_not_mistaken_for_a_suppressed_element():
    """ An empty dict is falsy, so a truth test rather than an `is False` check would drop it. """
    assert element_style({}, {"color": "red"}) == {"color": "red"}


# ---------------------------------------------------------------------------
# AxesHelper: airfoils
# ---------------------------------------------------------------------------

def sample_airfoil() -> AfGeometry:
    """
    An elliptical section, matching the synthetic airfoil in test_airfoil2.py, so that the plot
    tests need no data file and the geometry is known: the leading edge sits at (-10, 0), the
    trailing edge at (10, 0), and the chord is 20 units long.
    """
    t = numpy.linspace(0.0, 2.0 * numpy.pi, 721)
    section = Curve2(numpy.column_stack((10.0 * numpy.cos(t), 3.0 * numpy.sin(t))), tol=1e-8)
    return AfGeometry.from_geometric_analysis(
        section=section,
        general_tol=1e-3,
        fwd_aft=OrientFwdAft.fwd(-1.0, 0.0),
        upper_lower=OrientUpperLower.upper(0.0, 1.0),
        le_search="auto",
        te_search="auto",
    )


AIRFOIL_ELEMENTS = ["circles", "max_thickness", "upper", "lower", "camber", "edge_labels"]


def test_draw_airfoil_draws_every_element_by_default():
    result = new_helper().draw_airfoil(sample_airfoil())
    assert sorted(result) == sorted(AIRFOIL_ELEMENTS)


@pytest.mark.parametrize("element", AIRFOIL_ELEMENTS)
def test_draw_airfoil_can_suppress_any_single_element(element):
    helper = new_helper()
    before = artist_count(helper.ax)
    result = helper.draw_airfoil(sample_airfoil(), **{element: False})

    assert element not in result
    assert sorted(result) == sorted(e for e in AIRFOIL_ELEMENTS if e != element)
    assert artist_count(helper.ax) > before


def test_draw_airfoil_draws_one_patch_per_inscribed_circle():
    geom = sample_airfoil()
    result = new_helper().draw_airfoil(geom)
    assert len(result["circles"]) == len(geom.circle_array)


def test_draw_airfoil_merges_a_style_over_the_element_defaults():
    result = new_helper().draw_airfoil(sample_airfoil(), camber={"color": "red"})
    camber = result["camber"][0]
    assert camber.get_color() == "red"
    # The dashed default survives, rather than being reset by the override.
    assert camber.get_linestyle() == "--"


def test_draw_airfoil_places_the_edge_labels_outside_the_section():
    """ Both labels are pushed outward along the chord so neither lands on the geometry. """
    geom = sample_airfoil()
    result = new_helper().draw_airfoil(geom, label_offset=1.0)
    labels = {a.get_text(): a.xy for a in result["edge_labels"] if isinstance(a, Annotation)}

    assert labels["LE"][0] == pytest.approx(geom.leading.point.x - 1.0, abs=1e-6)
    assert labels["TE"][0] == pytest.approx(geom.trailing.point.x + 1.0, abs=1e-6)


def test_draw_airfoil_honors_a_zero_label_offset():
    """ Regression guard on the falsy-zero pattern: 0.0 must pin the labels to their points. """
    geom = sample_airfoil()
    result = new_helper().draw_airfoil(geom, label_offset=0.0)
    labels = {a.get_text(): a.xy for a in result["edge_labels"] if isinstance(a, Annotation)}

    assert labels["LE"][0] == pytest.approx(geom.leading.point.x, abs=1e-6)
    assert labels["TE"][0] == pytest.approx(geom.trailing.point.x, abs=1e-6)


def test_draw_airfoil_defaults_the_label_offset_to_a_fraction_of_the_chord():
    """ A scale-independent default keeps the layout the same on any size of section. """
    geom = sample_airfoil()
    result = new_helper().draw_airfoil(geom)
    labels = {a.get_text(): a.xy for a in result["edge_labels"] if isinstance(a, Annotation)}
    chord = (geom.trailing.point - geom.leading.point).norm()

    assert geom.leading.point.x - labels["LE"][0] == pytest.approx(chord * 0.04, abs=1e-6)


def test_draw_airfoil_returns_artists_that_are_really_on_the_axes():
    helper = new_helper()
    result = helper.draw_airfoil(sample_airfoil())
    for artists in result.values():
        for artist in artists:
            assert artist.axes is helper.ax


# ---------------------------------------------------------------------------
# Deviation scales and colorbars
# ---------------------------------------------------------------------------

def test_gom_color_map_builds_without_a_deprecation_warning():
    """
    Matplotlib 3.11 pends the deprecation of set_under/set_over in favor of constructor arguments,
    which older versions do not accept. Whichever path is taken must be the quiet one.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", PendingDeprecationWarning)
        warnings.simplefilter("error", DeprecationWarning)
        GomColorMap()


def test_deviation_limit_is_the_largest_magnitude_by_default():
    assert deviation_limit([-4.0, 1.0, 2.5]) == pytest.approx(4.0, abs=TOL)


def test_deviation_limit_ignores_non_finite_values():
    assert deviation_limit([1.0, numpy.nan, -2.0, numpy.inf]) == pytest.approx(2.0, abs=TOL)


def test_deviation_limit_can_clip_outliers_with_a_percentile():
    """ One wild point should not be able to flatten the whole scale. """
    values = numpy.concatenate([numpy.linspace(-1.0, 1.0, 99), [500.0]])
    assert deviation_limit(values, percentile=95.0) < 2.0
    assert deviation_limit(values) == pytest.approx(500.0, abs=TOL)


@pytest.mark.parametrize("bad", [-1.0, 100.5])
def test_deviation_limit_rejects_a_percentile_out_of_range(bad):
    with pytest.raises(ValueError, match="between 0 and 100"):
        deviation_limit([1.0, 2.0], percentile=bad)


@pytest.mark.parametrize("values", [[], [numpy.nan, numpy.inf]])
def test_deviation_limit_rejects_data_with_nothing_finite_in_it(values):
    with pytest.raises(ValueError, match="no finite value"):
        deviation_limit(values)


def test_deviation_norm_is_symmetric_about_zero():
    """ Zero has to land dead center, or the color that reads as nominal quietly moves. """
    norm = deviation_norm(0.75)
    assert norm.vmin == pytest.approx(-0.75, abs=TOL)
    assert norm.vmax == pytest.approx(0.75, abs=TOL)
    assert norm(0.0) == pytest.approx(0.5, abs=TOL)


def test_deviation_norm_ignores_the_sign_of_the_limit():
    assert deviation_norm(-2.0).vmin == deviation_norm(2.0).vmin


def test_deviation_norm_rejects_a_zero_limit():
    with pytest.raises(ValueError, match="non-zero"):
        deviation_norm(0.0)


def test_has_extremes_detects_the_gom_out_of_range_colors():
    assert has_extremes(GOM_CMAP) == (True, True)


def test_has_extremes_is_false_for_a_map_that_sets_no_extremes():
    assert has_extremes(matplotlib.pyplot.get_cmap("viridis")) == (False, False)


@pytest.mark.parametrize("cmap,expected", [
    (GOM_CMAP, "both"),
    (matplotlib.pyplot.get_cmap("viridis"), "neither"),
    (matplotlib.pyplot.get_cmap("viridis").with_extremes(over="red"), "max"),
    (matplotlib.pyplot.get_cmap("viridis").with_extremes(under="red"), "min"),
])
def test_extend_for_maps_the_extremes_to_a_colorbar_token(cmap, expected):
    assert extend_for(cmap) == expected


def test_draw_colorbar_from_a_norm_alone():
    helper = new_helper()
    bar = helper.draw_colorbar(norm=deviation_norm(1.0), label="deviation [mm]")
    assert bar.extend == "both"
    assert "deviation" in bar.ax.get_ylabel()


def test_draw_colorbar_takes_the_scale_from_a_mappable():
    helper = new_helper()
    scatter = helper.ax.scatter([0.0, 1.0], [0.0, 1.0], c=[-1.0, 1.0], cmap=GOM_CMAP,
                                norm=deviation_norm(1.0))
    assert helper.draw_colorbar(scatter).extend == "both"


def test_draw_colorbar_extends_only_where_the_color_map_calls_values_out():
    """
    Square ends beside a map that reserves colors for out-of-range data read as though the scale
    covers everything, which is the opposite of what those colors mean.
    """
    helper = new_helper()
    bar = helper.draw_colorbar(norm=deviation_norm(1.0),
                               cmap=matplotlib.pyplot.get_cmap("viridis"))
    assert bar.extend == "neither"


def test_draw_colorbar_honors_an_explicit_extend():
    helper = new_helper()
    assert helper.draw_colorbar(norm=deviation_norm(1.0), extend="max").extend == "max"


def test_draw_colorbar_needs_something_to_draw_a_scale_from():
    with pytest.raises(ValueError, match="mappable or a norm"):
        new_helper().draw_colorbar()


# ---------------------------------------------------------------------------
# AxesHelper: layout
# ---------------------------------------------------------------------------

def test_fill_available_space_equalizes_the_two_axis_scales():
    """
    The point of this over set_aspect("equal") is that the geometry keeps a 1:1 aspect without the
    plot shrinking inside its box, so both axes must end up at the same pixels per data unit.
    """
    figure = Figure(figsize=(12.0, 4.0))
    helper = AxesHelper(figure.subplots())
    helper.set_bounds(Aabb2(0.0, 0.0, 10.0, 10.0))
    figure.canvas.draw()

    helper.fill_available_space()

    x0, x1 = helper.ax.get_xlim()
    y0, y1 = helper.ax.get_ylim()
    box = helper.ax.get_window_extent()
    assert box.width / (x1 - x0) == pytest.approx(box.height / (y1 - y0), rel=1e-9)


def test_fill_available_space_only_ever_widens_the_limits():
    """ Filling space must not crop anything that was already in view. """
    figure = Figure(figsize=(12.0, 4.0))
    helper = AxesHelper(figure.subplots())
    helper.set_bounds(Aabb2(0.0, 0.0, 10.0, 10.0))
    figure.canvas.draw()

    helper.fill_available_space()

    assert helper.ax.get_xlim()[0] <= 0.0 and helper.ax.get_xlim()[1] >= 10.0
    assert helper.ax.get_ylim()[0] <= 0.0 and helper.ax.get_ylim()[1] >= 10.0


# ---------------------------------------------------------------------------
# Coverage drift
# ---------------------------------------------------------------------------

# Every public method on the helper classes, each of which must be exercised above. Adding a method
# without adding a test will fail these, in the spirit of test_stub_drift.py.
EXERCISED_HELPER_MEMBERS = {
    "draw_aabb", "draw_arc", "draw_arrow", "draw_boundary", "draw_circle", "draw_line",
    "draw_curve", "draw_distance", "draw_labeled_arrow", "draw_point", "draw_segment",
    "draw_airfoil", "draw_colorbar", "draw_normals", "draw_spline", "draw_surface_point",
    "draw_text", "fill_available_space", "fill_curve", "set_bounds",
    "viewport",
}

EXERCISED_VIEWPORT_MEMBERS = {
    "draw_aabb", "draw_circle", "draw_coordinate_system", "draw_curve", "draw_dimension_arrow",
    "draw_distance", "draw_labeled_point", "draw_line", "draw_mesh_outline", "draw_plane",
    "draw_point_cloud", "find_mesh_edge_point",
}

EXERCISED_TRACE_MEMBERS = {
    "add_blank", "add_point_and_color", "add_points", "add_segment", "bounds", "invert_y", "xy",
}


def public_members(cls) -> set:
    return {name for name in dir(cls) if not name.startswith("_")}


@pytest.mark.parametrize("cls,exercised", [
    (AxesHelper, EXERCISED_HELPER_MEMBERS),
    (ViewPort3, EXERCISED_VIEWPORT_MEMBERS),
    (TraceBuilder, EXERCISED_TRACE_MEMBERS),
])
def test_every_public_member_is_covered_by_a_test(cls, exercised):
    actual = public_members(cls)
    missing = actual - exercised
    stale = exercised - actual
    assert not missing, f"{cls.__name__} gained {sorted(missing)} with no test covering it"
    assert not stale, f"{cls.__name__} no longer has {sorted(stale)}; update the exercised set"
