"""
Tests for the PyVista plotting helper in `engeom.plot.pyvista`.

Split out from `test_plot.py` so that the two optional backends skip independently: with matplotlib
installed and PyVista absent, the matplotlib tests still run.

Nothing here renders. Building an actor does not need an OpenGL context, so these run on a headless
machine, but anything that would rasterize (`screenshot`, `show`) does, and is not exercised. As in
`test_plot.py`, the goal is that every public entry point is reached and that its documented edge
cases hold, not that the pixels are right.
"""

import numpy
import pytest

pytest.importorskip("pyvista")

import pyvista

from engeom.geom3 import (Circle3, Curve3, Iso3, Mesh3, Point3, PointCloud3, PointCloud3,
                          Sphere3, Vector3)
from engeom.metrology import Distance3
from engeom.plot._common import LABEL_PLACES
from engeom.plot.pyvista import LineBuilder, PlotterHelper, _cmap_extremes

TOL = 1e-12


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

def new_helper() -> PlotterHelper:
    """
    Build a helper over an off-screen plotter, so that no window is created and the tests do not
    require a display.
    """
    return PlotterHelper(pyvista.Plotter(off_screen=True))


def actor_count(helper: PlotterHelper) -> int:
    return len(helper.pv.renderer.actors)


def a_mesh() -> Mesh3:
    return Mesh3.create_box(2.0, 3.0, 4.0, True)


def a_curve() -> Curve3:
    return Curve3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]), tol=1e-6)


def a_cloud(normals: bool = False, colors: bool = False) -> PointCloud3:
    data = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))
    if normals:
        data.set_point_normals(numpy.array([[0.0, 0.0, 1.0]] * 3))
    if colors:
        data.set_point_colors(numpy.array([[255, 0, 0], [0, 255, 0], [0, 0, 255]],
                                          dtype=numpy.uint8))
    return data


def a_distance() -> Distance3:
    return Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Construction and delegation
# ---------------------------------------------------------------------------

def test_the_helper_exposes_the_plotter_it_wraps():
    """
    The helper wraps rather than subclasses, so the host object stays reachable for anything the
    helper does not cover.
    """
    plotter = pyvista.Plotter(off_screen=True)
    helper = PlotterHelper(plotter)
    assert helper.pv is plotter


def test_with_new_plotter_builds_its_own_plotter():
    helper = PlotterHelper.with_new_plotter(off_screen=True)
    assert isinstance(helper, PlotterHelper)
    assert isinstance(helper.pv, pyvista.Plotter)


def test_show_delegates_to_the_wrapped_plotter():
    """
    `show` cannot be called for real without a rendering context, so this checks only that it
    forwards, which is the whole of its behavior.
    """

    class FakePlotter:
        def __init__(self):
            self.shown = 0

        def show(self):
            self.shown += 1

    fake = FakePlotter()
    PlotterHelper(fake).show()
    assert fake.shown == 1


# ---------------------------------------------------------------------------
# Meshes, points, and curves
# ---------------------------------------------------------------------------

def test_draw_mesh_adds_an_actor():
    helper = new_helper()
    actor = helper.draw_mesh(a_mesh(), color="red")
    assert actor is not None
    assert actor_count(helper) == 1


def test_draw_mesh_carries_the_mesh_triangles_across():
    """
    The helper builds PyVista's padded face array by hand, so check the triangle count survives.
    """
    mesh = a_mesh()
    helper = new_helper()
    actor = helper.draw_mesh(mesh)
    assert actor.mapper.dataset.n_cells == mesh.faces.shape[0]


def test_draw_point_accepts_mixed_coordinate_forms():
    helper = new_helper()
    actor = helper.draw_point((0.0, 0.0, 0.0), Point3(1.0, 1.0, 1.0), [2.0, 2.0, 2.0])
    assert actor.mapper.dataset.n_points == 3


def test_draw_curve_accepts_several_curves_in_one_call():
    helper = new_helper()
    actor = helper.draw_curve(a_curve(), a_curve(), color="w")
    assert actor is not None
    assert actor_count(helper) == 1


def test_draw_sphere_adds_an_actor():
    sphere = Sphere3.from_fit(numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                           [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]))
    helper = new_helper()
    assert helper.draw_sphere(sphere) is not None
    assert actor_count(helper) == 1


# ---------------------------------------------------------------------------
# Point clouds
# ---------------------------------------------------------------------------

def test_draw_point_cloud_adds_one_actor_for_a_plain_cloud():
    helper = new_helper()
    actors = helper.draw_point_cloud(a_cloud())
    assert len(actors) == 1
    assert actor_count(helper) == 1


def test_draw_point_cloud_adds_arrows_only_for_a_positive_size():
    """
    A size of zero means "no arrows", which is the documented default. The check is `> 0.0` rather
    than `>= 0.0` precisely so that the default does not draw degenerate arrows.
    """
    cloud = a_cloud(normals=True)
    assert len(new_helper().draw_point_cloud(cloud)) == 1
    assert len(new_helper().draw_point_cloud(cloud, normal_arrow_size=0.0)) == 1
    assert len(new_helper().draw_point_cloud(cloud, normal_arrow_size=0.5)) == 2


def test_draw_point_cloud_skips_arrows_when_the_cloud_has_no_normals():
    helper = new_helper()
    assert len(helper.draw_point_cloud(a_cloud(), normal_arrow_size=0.5)) == 1


def test_draw_point_cloud_prefers_the_clouds_own_colors():
    """
    When the cloud carries colors and `use_colors` is on, an explicit `color` has to be dropped:
    PyVista rejects being given both.
    """
    helper = new_helper()
    actors = helper.draw_point_cloud(a_cloud(colors=True), color="red")
    assert len(actors) == 1


def test_draw_point_cloud_honors_an_explicit_color_when_colors_are_off():
    helper = new_helper()
    actors = helper.draw_point_cloud(a_cloud(colors=True), use_colors=False, color="red")
    assert len(actors) == 1


# ---------------------------------------------------------------------------
# Circles
# ---------------------------------------------------------------------------

def test_draw_circle_draws_both_an_edge_and_a_face_by_default():
    helper = new_helper()
    actors = helper.draw_circle(Circle3(1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 2.0))
    assert len(actors) == 2
    assert actor_count(helper) == 2


@pytest.mark.parametrize("kwargs,expected", [
    ({}, 2),
    ({"face_color": None}, 1),
    ({"edge_color": None}, 1),
    ({"edge_color": None, "face_color": None}, 0),
])
def test_draw_circle_suppresses_each_part_with_a_none_color(kwargs, expected):
    helper = new_helper()
    assert len(helper.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0), **kwargs)) == expected


def test_draw_circle_puts_its_vertices_on_the_circle():
    """
    A `Circle3` gives a center, a normal, and a radius but no in-plane orientation, so the helper
    picks a basis itself. Whichever basis it picks, every vertex must lie at the radius and in the
    plane; this fails outright if the construction is wrong, which it was when the helper called a
    `Circle3.iso` property that does not exist.
    """
    center = numpy.array([1.0, 2.0, 3.0])
    normal = numpy.array([1.0, 1.0, 1.0]) / numpy.sqrt(3.0)
    circle = Circle3(*center, *normal, 2.0)

    helper = new_helper()
    edge, _ = helper.draw_circle(circle, n=64)
    points = numpy.asarray(edge.mapper.dataset.points)

    radii = numpy.linalg.norm(points - center, axis=1)
    assert radii == pytest.approx(2.0, abs=1e-9)
    assert numpy.abs((points - center) @ normal).max() == pytest.approx(0.0, abs=1e-9)


def test_draw_circle_puts_its_face_where_the_circle_is():
    """
    The face is a separate construction from the edge, built by transforming an origin-centered
    mesh, so it gets its own check that the two land in the same place.
    """
    circle = Circle3(1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 2.0)
    helper = new_helper()
    _, face = helper.draw_circle(circle, n=64)
    points = numpy.asarray(face.mapper.dataset.points)

    assert points.mean(axis=0) == pytest.approx([1.0, 2.0, 3.0], abs=1e-9)
    assert numpy.abs(points[:, 1] - 2.0).max() == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label_place", LABEL_PLACES)
def test_draw_distance_handles_every_label_placement(label_place):
    helper = new_helper()
    actors = helper.draw_distance(a_distance(), label_place=label_place)
    assert len(actors) == 3
    assert actor_count(helper) == 3


def test_draw_distance_rejects_an_unknown_label_placement():
    """
    Validating up front is what keeps an unknown token from falling through the placement branches
    into an `UnboundLocalError` on the label coordinates.
    """
    with pytest.raises(ValueError, match="invalid label_place"):
        new_helper().draw_distance(a_distance(), label_place="sideways")


def test_draw_distance_scales_the_value_in_the_label_only():
    helper = new_helper()
    distance = a_distance()
    helper.draw_distance(distance, template="{value:.2f}", scale_value=25.4)
    assert distance.value == pytest.approx(5.0, abs=TOL)


def test_draw_distance_handles_a_negative_measurement():
    """
    The leader direction flips with the sign of the value, so the negative case is a separate path.
    """
    distance = Distance3(Point3(0.0, 0.0, 5.0), Point3(0.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0))
    assert distance.value < 0.0
    assert len(new_helper().draw_distance(distance)) == 3


# ---------------------------------------------------------------------------
# Coordinate systems, labels, and arrows
# ---------------------------------------------------------------------------

def test_draw_coordinate_system_adds_one_actor_per_axis():
    helper = new_helper()
    assert len(helper.draw_coordinate_system(Iso3.identity(), size=1.0)) == 3


def test_draw_coordinate_system_adds_a_label_when_asked():
    helper = new_helper()
    assert len(helper.draw_coordinate_system(Iso3.identity(), size=1.0, label="A")) == 4


def test_draw_coordinate_system_accepts_a_4x4_array():
    helper = new_helper()
    assert len(helper.draw_coordinate_system(numpy.eye(4))) == 3


def test_draw_coordinate_system_accepts_anything_with_an_as_numpy_method():
    """
    The duck-typed path exists so an alignment result can be handed straight in without unwrapping.
    """

    class HasAsNumpy:
        def as_numpy(self):
            return numpy.eye(4)

    assert len(new_helper().draw_coordinate_system(HasAsNumpy())) == 3


def test_draw_coordinate_system_rejects_a_wrongly_shaped_array():
    with pytest.raises(ValueError, match=r"expected \(4, 4\)"):
        new_helper().draw_coordinate_system(numpy.eye(3))


def test_draw_coordinate_system_rejects_an_unusable_type():
    with pytest.raises(TypeError, match="expected Iso3 or numpy.ndarray"):
        new_helper().draw_coordinate_system("not an isometry")


def test_draw_coordinate_system_puts_the_axes_at_the_isometrys_origin():
    helper = new_helper()
    actors = helper.draw_coordinate_system(Iso3.from_translation(10.0, 0.0, 0.0), size=2.0)
    x_axis = numpy.asarray(actors[0].mapper.dataset.points)
    assert x_axis[0] == pytest.approx([10.0, 0.0, 0.0], abs=1e-9)
    assert x_axis[1] == pytest.approx([12.0, 0.0, 0.0], abs=1e-9)


def test_draw_label_returns_the_label_it_added():
    helper = new_helper()
    label = helper.draw_label((1.0, 2.0, 3.0), "hello")
    assert isinstance(label, pyvista.Label)
    assert actor_count(helper) == 1


def test_draw_arrow_adds_an_actor():
    helper = new_helper()
    assert helper.draw_arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)) is not None
    assert actor_count(helper) == 1


def test_draw_arrow_accepts_an_explicit_color():
    """
    The color used to be passed as a second literal `color=` alongside `**kwargs`, so supplying one
    raised a duplicate-argument `TypeError` rather than recoloring the arrow.
    """
    helper = new_helper()
    assert helper.draw_arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), color="red") is not None


# ---------------------------------------------------------------------------
# Color map extremes
# ---------------------------------------------------------------------------

def test_cmap_extremes_carries_a_matplotlib_maps_over_and_under_across():
    """
    `GOM_CMAP` calls out data beyond the scale with distinct over and under colors, and PyVista has
    its own arguments for that; without this the map would silently clamp instead.
    """
    pytest.importorskip("matplotlib")
    from engeom.plot.matplotlib import GOM_CMAP

    extremes = _cmap_extremes(GOM_CMAP)
    assert set(extremes) == {"above_color", "below_color"}


def test_cmap_extremes_returns_colors_pyvista_can_use():
    """
    `Colormap.get_over` returns a numpy array, and PyVista tests these arguments for truthiness
    internally, which raises on an array. They have to come back as plain tuples.
    """
    pytest.importorskip("matplotlib")
    from engeom.plot.matplotlib import GOM_CMAP

    for value in _cmap_extremes(GOM_CMAP).values():
        assert isinstance(value, tuple)
        assert all(isinstance(c, float) for c in value)
        assert bool(value)  # the operation PyVista performs, which an array cannot answer


def test_cmap_extremes_carries_nothing_across_for_a_map_that_sets_neither():
    """
    Matplotlib has no public way to ask whether an extreme was set, so an unset one reports the
    map's own end color. Those must not be forwarded, or every map would get spurious extremes.
    """
    pytest.importorskip("matplotlib")
    from matplotlib import colormaps

    assert _cmap_extremes(colormaps["viridis"]) == {}


def test_cmap_extremes_ignores_anything_that_is_not_a_colormap():
    assert _cmap_extremes("viridis") == {}
    assert _cmap_extremes(None) == {}


def test_draw_mesh_forwards_the_color_map_extremes():
    """
    The end-to-end version of the above: this is the call that failed outright while the extremes
    were being handed to PyVista as numpy arrays.
    """
    pytest.importorskip("matplotlib")
    from engeom.plot.matplotlib import GOM_CMAP

    mesh = a_mesh()
    helper = new_helper()
    actor = helper.draw_mesh(mesh, scalars=numpy.zeros(mesh.points.shape[0]), cmap=GOM_CMAP)
    assert actor is not None


# ---------------------------------------------------------------------------
# LineBuilder
# ---------------------------------------------------------------------------

def test_line_builder_repeats_interior_vertices_for_disjoint_segments():
    """
    `Plotter.add_lines` takes vertices in pairs, so a run of connected points has to have its
    interior vertices doubled.
    """
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    builder.add((1.0, 0.0, 0.0))
    builder.add((2.0, 0.0, 0.0))
    vertices = builder.build()
    assert vertices.shape == (4, 3)
    assert vertices[1] == pytest.approx(vertices[2], abs=TOL)


def test_line_builder_skip_starts_a_new_run():
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    builder.add((1.0, 0.0, 0.0))
    builder.skip()
    builder.add((5.0, 0.0, 0.0))
    builder.add((6.0, 0.0, 0.0))
    vertices = builder.build()
    assert vertices.shape == (4, 3)
    assert vertices[2] == pytest.approx([5.0, 0.0, 0.0], abs=TOL)


# ---------------------------------------------------------------------------
# Coverage drift
# ---------------------------------------------------------------------------

# Every public method on the helper, each of which must be exercised above. Adding a method without
# adding a test will fail this, in the spirit of test_stub_drift.py.
EXERCISED_HELPER_MEMBERS = {
    "draw_arrow", "draw_circle", "draw_coordinate_system", "draw_curve", "draw_distance",
    "draw_label", "draw_mesh", "draw_point", "draw_point_cloud", "draw_sphere", "show",
    "with_new_plotter",
}


def test_every_public_member_is_covered_by_a_test():
    actual = {name for name in dir(PlotterHelper) if not name.startswith("_")}
    missing = actual - EXERCISED_HELPER_MEMBERS
    stale = EXERCISED_HELPER_MEMBERS - actual
    assert not missing, f"PlotterHelper gained {sorted(missing)} with no test covering it"
    assert not stale, f"PlotterHelper no longer has {sorted(stale)}; update the exercised set"


def test_every_public_member_is_a_draw_method_or_documented_otherwise():
    """
    The naming convention is what makes the surface discoverable by autocomplete: anything that adds
    an actor is prefixed `draw_`, and anything without the prefix configures or queries.
    """
    not_draw = {name for name in dir(PlotterHelper)
                if not name.startswith("_") and not name.startswith("draw_")}
    assert not_draw == {"show", "with_new_plotter"}
