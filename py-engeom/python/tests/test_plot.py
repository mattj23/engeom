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

from matplotlib.figure import Figure

from engeom.geom2 import (Aabb2, Arc2, BoundaryData2, Circle2, CubicSpline2, Curve2, Point2,
                          Segment2, SurfacePoint2, Vector2)
from engeom.geom3 import Iso3, Line3, Mesh3, Point3, Vector3
from engeom.metrology import Distance2
from engeom.plot import LabelPlace
from engeom.plot._coerce import to_point2, to_point3, to_tuple2, to_tuple3
from engeom.plot.matplotlib import GOM_CMAP, GomColorMap, AxesHelper, TraceBuilder, ViewPort3

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


def test_draw_boundary_normals_adds_one_arrow_per_sample():
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_boundary_normals(sample_boundary(), 5, 0.1, color="red")
    assert artist_count(helper.ax) == before + 5


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

@pytest.mark.parametrize("label_place", [LabelPlace.Outside, LabelPlace.Inside, LabelPlace.OutsideRev])
def test_draw_distance_adds_artists_for_every_label_placement(label_place):
    helper = new_helper()
    before = artist_count(helper.ax)
    helper.draw_distance(Distance2(Point2(0.0, 0.0), Point2(3.0, 0.0)), label_place=label_place)
    assert artist_count(helper.ax) > before


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


def test_viewport_line_draws():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.line(Line3(0.0, 0.0, 0.0, 1.0, 1.0, 0.0), t=2.0)
    assert artist_count(view.helper.ax) > before


def test_viewport_line_accepts_an_explicit_start():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.line(Line3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0), t=2.0, t0=1.0)
    assert artist_count(view.helper.ax) > before


def test_viewport_coordinate_system_draws_the_visible_axes():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.coordinate_system(Iso3.identity(), 1.0)
    # Looking down Z, the X and Y axes are visible and the Z axis is hidden; each visible axis
    # contributes an arrow and a label.
    assert artist_count(view.helper.ax) == before + 4


def test_viewport_labeled_point_draws():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.labeled_point(Point3(0.0, 0.0, 0.0), "P")
    assert artist_count(view.helper.ax) > before


def test_viewport_labeled_point_can_draw_a_leader_arrow():
    view = new_viewport()
    plain_before = artist_count(view.helper.ax)
    view.labeled_point(Point3(0.0, 0.0, 0.0), "P")
    plain = artist_count(view.helper.ax) - plain_before

    view = new_viewport()
    arrow_before = artist_count(view.helper.ax)
    view.labeled_point(Point3(0.0, 0.0, 0.0), "P", offset_2d=(1.0, 1.0), arrow=True)
    assert artist_count(view.helper.ax) - arrow_before == plain + 1


def test_viewport_dimension_arrow_draws():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0")
    assert artist_count(view.helper.ax) > before


def test_viewport_dimension_arrow_adds_leaders_when_shifted():
    view = new_viewport()
    plain_before = artist_count(view.helper.ax)
    view.dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0")
    plain = artist_count(view.helper.ax) - plain_before

    view = new_viewport()
    shifted_before = artist_count(view.helper.ax)
    view.dimension_arrow(Point3(0.0, 0.0, 0.0), Point3(2.0, 0.0, 0.0), "2.0", leader_shift=(0.0, 1.0, 0.0))
    assert artist_count(view.helper.ax) - shifted_before == plain + 2


def test_viewport_mesh_outline_draws():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0))
    assert artist_count(view.helper.ax) > before


def test_viewport_mesh_outline_can_suppress_hidden_edges():
    view = new_viewport()
    before = artist_count(view.helper.ax)
    view.mesh_outline(Mesh3.create_box(1.0, 1.0, 1.0), no_hidden=True)
    # Only the visible trace is plotted, rather than a visible and a hidden trace.
    assert artist_count(view.helper.ax) == before + 1


def test_viewport_mesh_edge_point_in_dir_returns_a_point_on_the_mesh():
    view = new_viewport()
    mesh = Mesh3.create_box(2.0, 2.0, 2.0)
    found = view.mesh_edge_point_in_dir(1.0, 0.0, mesh)
    assert isinstance(found, Point3)
    # The box spans -1 to 1 on each axis, and the query asks for the +X extreme.
    assert found.x == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Coverage drift
# ---------------------------------------------------------------------------

# Every public method on the helper classes, each of which must be exercised above. Adding a method
# without adding a test will fail these, in the spirit of test_stub_drift.py.
EXERCISED_HELPER_MEMBERS = {
    "draw_arc", "draw_arrow", "draw_boundary", "draw_boundary_normals", "draw_circle",
    "draw_curve", "draw_distance", "draw_labeled_arrow", "draw_point", "draw_segment",
    "draw_spline", "draw_surface_point", "draw_text", "fill_curve", "set_bounds", "viewport",
}

EXERCISED_VIEWPORT_MEMBERS = {
    "coordinate_system", "dimension_arrow", "labeled_point", "line", "mesh_edge_point_in_dir",
    "mesh_outline",
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
