"""
Tests for the PyVista plotting helper in `engeom.plot.pyvista`.

Split out from `test_plot.py` so that the two optional backends skip independently: with matplotlib
installed and PyVista absent, the matplotlib tests still run.

Nothing here renders. Building an actor does not need an OpenGL context, so these run on a headless
machine, but anything that would rasterize (`screenshot`, `show`) does, and is not exercised. As in
`test_plot.py`, the goal is that every public entry point is reached and that its documented edge
cases hold, not that the pixels are right.
"""

import importlib.metadata
import subprocess
import sys
import textwrap
import warnings

import numpy
import pytest

pytest.importorskip("pyvista")

import pyvista

from engeom.geom3 import (Aabb3, Circle3, Cone3, CubicSpline3, Curve3, CurveGroup3, Cylinder3,
                          Iso3, Line3, Mesh3, Plane3, Point3, PointCloud3, Segment3, Sphere3,
                          SurfacePoint3, SvdBasis3, Vector3)
from engeom.metrology import Distance3
from engeom.common import IndexMask
from engeom.plot._common import LABEL_PLACES
import engeom.plot.pyvista as engeom_pyvista
from engeom.plot.pyvista import EngeomPlotter, LineBuilder, PlotterHelper
from engeom.plot.pyvista.convert import FACE_ID, POINT_ID, _lines_polydata, to_mesh3, to_polydata
from engeom.plot.pyvista.extent import clip_line_to_aabb, clip_plane_to_aabb, resolve_extent
from engeom.plot.pyvista.helper import _cmap_extremes

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


def test_draw_curve_gives_each_curve_its_own_actor():
    helper = new_helper()
    actors = helper.draw_curve(a_curve(), a_curve(), color="w")
    assert len(actors) == 2
    assert actor_count(helper) == 2


def test_draw_curve_draws_a_curve_group_as_a_single_actor():
    """
    A `CurveGroup3` is one entity that moves as a rigid body, so it draws as one actor holding a
    cell per curve rather than as several actors.
    """
    helper = new_helper()
    group = CurveGroup3([a_curve(), a_curve(), a_curve()])
    actors = helper.draw_curve(group)
    assert len(actors) == 1
    assert actors[0].mapper.dataset.n_cells == 3


def test_draw_sphere_adds_an_actor():
    sphere = Sphere3.from_fit(numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                                           [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]]))
    helper = new_helper()
    assert len(helper.draw_sphere(sphere)) == 1
    assert actor_count(helper) == 1


def test_draw_sphere_accepts_several_spheres():
    sphere = Sphere3(0.0, 0.0, 0.0, 1.0)
    helper = new_helper()
    assert len(helper.draw_sphere(sphere, sphere, sphere)) == 3
    assert actor_count(helper) == 3


# ---------------------------------------------------------------------------
# Point clouds
# ---------------------------------------------------------------------------

def test_draw_point_cloud_adds_one_actor_for_a_plain_cloud():
    helper = new_helper()
    actors = helper.draw_point_cloud(a_cloud())
    assert len(actors) == 1
    assert actor_count(helper) == 1


@pytest.mark.parametrize("normals,expected", [
    (False, 1),
    (True, 2),
    ({"mag": 0.5}, 2),
])
def test_draw_point_cloud_draws_normals_only_when_asked(normals, expected):
    """
    Normal arrows are a sub-element styled the same way as any other composite element: `False`
    leaves them out, which is the default, and `True` or a dict of `add_arrows` keyword arguments
    draws them.
    """
    cloud = a_cloud(normals=True)
    assert len(new_helper().draw_point_cloud(cloud, normals=normals)) == expected


def test_draw_point_cloud_skips_arrows_when_the_cloud_has_no_normals():
    helper = new_helper()
    assert len(helper.draw_point_cloud(a_cloud(), normals=True)) == 1


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

def test_draw_circle_draws_one_actor_holding_face_and_outline():
    """
    The face and the outline are a single polygon actor, so they cannot z-fight with each other the
    way two coplanar actors did.
    """
    helper = new_helper()
    actors = helper.draw_circle(Circle3(1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 2.0))
    assert len(actors) == 1
    assert actor_count(helper) == 1
    assert actors[0].mapper.dataset.n_cells == 1
    assert actors[0].prop.show_edges


@pytest.mark.parametrize("kwargs,expect_edges", [
    ({}, True),
    ({"show_edges": False}, False),
])
def test_draw_circle_can_leave_the_outline_off(kwargs, expect_edges):
    helper = new_helper()
    actors = helper.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0), **kwargs)
    assert actors[0].prop.show_edges is expect_edges


def test_draw_circle_draws_an_outline_alone_as_a_wireframe():
    helper = new_helper()
    actors = helper.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0), style="wireframe")
    assert actors[0].prop.style == "Wireframe"


def test_draw_circle_accepts_an_opacity():
    """
    Opacity used to reach `Plotter.add_lines`, which has no `opacity` parameter, so asking for a
    translucent circle raised a `TypeError` instead of drawing one.
    """
    helper = new_helper()
    actors = helper.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0), opacity=0.5)
    assert actors[0].prop.opacity == pytest.approx(0.5)


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
    actor, = helper.draw_circle(circle, tol=1.0e-3)
    points = numpy.asarray(actor.mapper.dataset.points)

    radii = numpy.linalg.norm(points - center, axis=1)
    assert radii == pytest.approx(2.0, abs=1e-9)
    assert numpy.abs((points - center) @ normal).max() == pytest.approx(0.0, abs=1e-9)


def test_draw_circle_puts_its_polygon_where_the_circle_is():
    circle = Circle3(1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 2.0)
    helper = new_helper()
    actor, = helper.draw_circle(circle, tol=1.0e-3)
    points = numpy.asarray(actor.mapper.dataset.points)

    assert points.mean(axis=0) == pytest.approx([1.0, 2.0, 3.0], abs=1e-9)
    assert numpy.abs(points[:, 1] - 2.0).max() == pytest.approx(0.0, abs=1e-9)


def test_draw_circle_does_not_repeat_the_closing_vertex():
    """
    A polygon cell closes itself, so a duplicated first and last point would leave a degenerate
    zero-length edge in the outline.
    """
    helper = new_helper()
    actor, = helper.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0), tol=1.0e-3)
    points = numpy.asarray(actor.mapper.dataset.points)
    assert numpy.linalg.norm(points[0] - points[-1]) > 1e-6


# ---------------------------------------------------------------------------
# The remaining bounded primitives
# ---------------------------------------------------------------------------

def test_draw_cylinder_lands_on_its_own_axis():
    """
    A tessellated primitive is built at the origin and moved into place, so the check that matters
    is that it ends up where the entity says it is.
    """
    cylinder = Cylinder3(1.0, 2.0, 3.0, 0.0, 0.0, 1.0, 2.0, 10.0)
    helper = new_helper()
    actor, = helper.draw_cylinder(cylinder)
    points = numpy.asarray(actor.mapper.dataset.points)

    assert points[:, 2].min() == pytest.approx(3.0, abs=1e-6)
    assert points[:, 2].max() == pytest.approx(13.0, abs=1e-6)
    radii = numpy.linalg.norm(points[:, :2] - numpy.array([1.0, 2.0]), axis=1)
    assert radii.max() == pytest.approx(2.0, abs=1e-3)


def test_draw_cone_puts_its_apex_on_the_tip():
    """
    `create_cone` builds a cone centered on the origin with its apex towards +Z, so placing one
    means both rotating onto the entity's axis and shifting off the midpoint of that axis.
    """
    cone = Cone3(1.0, 2.0, 13.0, 0.0, 0.0, -1.0, 10.0, 2.0)
    helper = new_helper()
    actor, = helper.draw_cone(cone)
    points = numpy.asarray(actor.mapper.dataset.points)

    apex = points[numpy.argmax(points[:, 2])]
    assert apex == pytest.approx([1.0, 2.0, 13.0], abs=1e-6)
    assert points[:, 2].min() == pytest.approx(3.0, abs=1e-6)
    base = points[numpy.abs(points[:, 2] - 3.0) < 1e-6]
    assert numpy.linalg.norm(base[:, :2] - numpy.array([1.0, 2.0]), axis=1).max() == pytest.approx(
        2.0, abs=1e-3)


def test_draw_cone_follows_a_tilted_axis():
    cone = Cone3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0, 1.0)
    actor, = new_helper().draw_cone(cone)
    points = numpy.asarray(actor.mapper.dataset.points)
    assert points[numpy.argmin(points[:, 0])] == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)
    assert points[:, 0].max() == pytest.approx(10.0, abs=1e-6)


def test_draw_spline_flattens_to_a_polyline():
    spline = CubicSpline3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0)
    actor, = new_helper().draw_spline(spline)
    assert actor.mapper.dataset.n_cells == 1
    assert actor.mapper.dataset.n_points > 2


def test_draw_spline_takes_its_tolerance_from_the_spline_when_none_is_given():
    """
    A tolerance relative to the spline's own size keeps a small one and a large one equally smooth,
    which a fixed default could not do.
    """
    small = CubicSpline3(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, -1.0, 0.0, 3.0, 0.0, 0.0)
    large = CubicSpline3(0.0, 0.0, 0.0, 100.0, 100.0, 0.0, 200.0, -100.0, 0.0, 300.0, 0.0, 0.0)
    helper = new_helper()
    a, = helper.draw_spline(small)
    b, = helper.draw_spline(large)
    assert a.mapper.dataset.n_points == b.mapper.dataset.n_points


def test_draw_spline_accepts_several_splines():
    spline = CubicSpline3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0)
    assert len(new_helper().draw_spline(spline, spline)) == 2


def test_draw_surface_point_draws_a_point_and_an_arrow():
    helper = new_helper()
    actors = helper.draw_surface_point(SurfacePoint3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
                                       SurfacePoint3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0), length=1.0)
    assert len(actors) == 2
    assert actors[0].mapper.dataset.n_points == 2


def test_draw_surface_point_sizes_its_arrows_to_the_scene():
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 100.0, 100.0, True))
    _, arrows = helper.draw_surface_point(SurfacePoint3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    expected = numpy.linalg.norm([100.0, 100.0, 100.0]) * 0.1
    assert arrows.bounds[5] - arrows.bounds[4] == pytest.approx(expected, rel=0.05)


def test_draw_surface_point_needs_no_special_case_for_nothing():
    assert new_helper().draw_surface_point() == []


def test_draw_basis_scales_each_axis_by_its_own_spread():
    """
    The point of drawing a basis rather than a plain frame is that the axis lengths say how far the
    data actually spreads along each direction.
    """
    # A cloud that is wide in x, narrower in y, and almost flat in z.
    rng = numpy.random.default_rng(0)
    points = rng.normal(size=(400, 3)) * numpy.array([10.0, 3.0, 0.1])
    basis = SvdBasis3(points)

    actors = new_helper().draw_basis(basis)
    assert len(actors) == 3

    lengths = []
    for actor in actors:
        ends = numpy.asarray(actor.mapper.dataset.points)
        lengths.append(float(numpy.linalg.norm(ends[1] - ends[0])))
    assert lengths == sorted(lengths, reverse=True)
    assert lengths[0] == pytest.approx(basis.basis_stdevs[0], rel=1e-6)


# ---------------------------------------------------------------------------
# Extents, and the unbounded entities that need them
# ---------------------------------------------------------------------------

def a_box() -> Aabb3:
    return Aabb3(0.0, 0.0, 0.0, 10.0, 10.0, 10.0)


def test_an_extent_can_come_from_the_scene():
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    box = resolve_extent(None, helper.pv, pad=0.0)
    assert box.extent.as_numpy() == pytest.approx([10.0, 10.0, 10.0], abs=1e-9)


def test_an_extent_from_an_empty_scene_is_refused():
    """
    PyVista reports a two-unit cube for a renderer with nothing in it, which is a placeholder and
    not a measurement. Sizing a plane against it would draw something plausible and wrong.
    """
    with pytest.raises(ValueError, match="nothing in the scene"):
        resolve_extent(None, new_helper().pv)


def test_an_extent_ignores_a_scene_holding_only_flat_annotation():
    """
    A text label is a 2-D actor and contributes no bounds, so a scene holding only labels reports
    the same placeholder as an empty one.
    """
    helper = new_helper()
    helper.draw_label((0.0, 0.0, 0.0), "note")
    with pytest.raises(ValueError, match="nothing in the scene"):
        resolve_extent(None, helper.pv)


@pytest.mark.parametrize("build", [
    lambda: a_box(),
    lambda: Mesh3.create_box(10.0, 10.0, 10.0, True),
    lambda: PointCloud3(numpy.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]])),
    lambda: numpy.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]]),
    lambda: Segment3(0.0, 0.0, 0.0, 10.0, 10.0, 10.0),
])
def test_an_extent_can_come_from_anything_carrying_a_box(build):
    box = resolve_extent(build(), pad=0.0)
    assert box.extent.as_numpy() == pytest.approx([10.0, 10.0, 10.0], abs=1e-9)


def test_an_extent_is_padded_by_a_fraction_of_its_diagonal():
    grown = resolve_extent(a_box(), pad=0.1)
    expected = numpy.linalg.norm([10.0, 10.0, 10.0]) * 0.1
    assert grown.min.x == pytest.approx(-expected, abs=1e-9)
    assert grown.max.x == pytest.approx(10.0 + expected, abs=1e-9)


def test_an_extent_with_no_size_is_refused():
    """ A single point, or a set of coincident ones, gives a box a plane cannot be drawn across. """
    with pytest.raises(ValueError, match="no size"):
        resolve_extent(numpy.zeros((4, 3)))


def test_an_extent_refuses_something_that_carries_no_box():
    with pytest.raises(TypeError):
        resolve_extent("not an entity")


@pytest.mark.parametrize("plane,vertices", [
    # Straight through the middle, meeting four faces.
    (Plane3.from_point_normal(5.0, 5.0, 5.0, 0.0, 0.0, 1.0), 4),
    # Slicing a corner off, meeting three.
    (Plane3.from_point_normal(1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 3),
    # Through the body diagonal, which is the hexagonal section.
    (Plane3.from_point_normal(5.0, 5.0, 5.0, 1.0, 1.0, 1.0), 6),
    # Lying exactly in a face, where an edge-intersection approach degenerates.
    (Plane3.from_point_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0), 4),
])
def test_a_plane_clips_to_the_polygon_where_it_crosses_the_box(plane, vertices):
    polygon = clip_plane_to_aabb(plane, a_box())
    assert polygon.shape == (vertices, 3)
    # Every vertex has to be on the plane and inside the box.
    assert numpy.abs([plane.signed_distance_to_point(Point3(*p)) for p in polygon]).max() < 1e-9
    assert (polygon >= -1e-9).all() and (polygon <= 10.0 + 1e-9).all()


@pytest.mark.parametrize("plane", [
    # Clear of the box entirely.
    Plane3.from_point_normal(0.0, 0.0, 99.0, 0.0, 0.0, 1.0),
    # Touching a single corner, so the intersection has no area to draw.
    Plane3.from_point_normal(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
])
def test_a_plane_that_does_not_cross_the_box_clips_to_nothing(plane):
    assert clip_plane_to_aabb(plane, a_box()) is None


def test_a_line_clips_to_the_span_inside_the_box():
    span = clip_line_to_aabb(Line3(5.0, 5.0, -100.0, 0.0, 0.0, 1.0), a_box())
    assert span == pytest.approx(numpy.array([[5.0, 5.0, 0.0], [5.0, 5.0, 10.0]]), abs=1e-9)


def test_a_line_clips_the_same_way_whatever_its_direction_scale():
    """ A `Line3` stores its direction as given, so clipping cannot assume it is a unit vector. """
    unit = clip_line_to_aabb(Line3(5.0, 5.0, -100.0, 0.0, 0.0, 1.0), a_box())
    scaled = clip_line_to_aabb(Line3(5.0, 5.0, -100.0, 0.0, 0.0, 7.0), a_box())
    assert unit == pytest.approx(scaled, abs=1e-9)


def test_a_line_lying_along_an_edge_still_clips():
    span = clip_line_to_aabb(Line3(0.0, 0.0, -1.0, 0.0, 0.0, 1.0), a_box())
    assert span == pytest.approx(numpy.array([[0.0, 0.0, 0.0], [0.0, 0.0, 10.0]]), abs=1e-9)


def test_a_line_that_misses_the_box_clips_to_nothing():
    assert clip_line_to_aabb(Line3(-5.0, -5.0, 0.0, 0.0, 0.0, 1.0), a_box()) is None


def test_a_line_pointing_nowhere_is_refused():
    with pytest.raises(ValueError, match="no direction"):
        clip_line_to_aabb(Line3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), a_box())


def test_draw_plane_adds_one_polygon_actor_per_plane():
    helper = new_helper()
    planes = [Plane3.from_point_normal(5.0, 5.0, 5.0, 0.0, 0.0, 1.0),
              Plane3.from_point_normal(5.0, 5.0, 5.0, 1.0, 0.0, 0.0)]
    actors = helper.draw_plane(*planes, extent=a_box())
    assert len(actors) == 2
    assert all(a.mapper.dataset.n_cells == 1 for a in actors)
    assert actors[0].prop.opacity == pytest.approx(0.5)


def test_draw_plane_takes_its_extent_from_the_scene_by_default():
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    actor, = helper.draw_plane(Plane3.from_point_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    points = numpy.asarray(actor.mapper.dataset.points)
    # The mesh is centered on the origin, so the padded extent reaches past its faces.
    assert points[:, 0].max() > 5.0


def test_draw_plane_says_which_plane_missed_the_extent():
    """
    Silently skipping would shift every later actor in the returned list out of step with the
    planes that were asked for.
    """
    helper = new_helper()
    good = Plane3.from_point_normal(5.0, 5.0, 5.0, 0.0, 0.0, 1.0)
    bad = Plane3.from_point_normal(0.0, 0.0, 999.0, 0.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="Plane 1"):
        helper.draw_plane(good, bad, extent=a_box())


def test_draw_line_adds_one_actor_per_line():
    helper = new_helper()
    actors = helper.draw_line(Line3(5.0, 5.0, 0.0, 0.0, 0.0, 1.0),
                              Line3(5.0, 0.0, 5.0, 0.0, 1.0, 0.0), extent=a_box())
    assert len(actors) == 2


def test_draw_line_says_which_line_missed_the_extent():
    helper = new_helper()
    with pytest.raises(ValueError, match="Line 0"):
        helper.draw_line(Line3(-99.0, -99.0, 0.0, 0.0, 0.0, 1.0), extent=a_box())


def test_drawing_an_unbounded_entity_into_an_empty_scene_explains_itself():
    helper = new_helper()
    with pytest.raises(ValueError, match="nothing in the scene"):
        helper.draw_plane(Plane3.from_point_normal(0.0, 0.0, 0.0, 0.0, 0.0, 1.0))


def test_drawing_no_unbounded_entities_needs_no_extent():
    """ A computed and possibly empty collection should not have to be special-cased. """
    helper = new_helper()
    assert helper.draw_plane() == []
    assert helper.draw_line() == []


def test_draw_segment_adds_one_actor_per_segment():
    helper = new_helper()
    actors = helper.draw_segment(Segment3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
                                 Segment3(0.0, 0.0, 0.0, 0.0, 1.0, 0.0))
    assert len(actors) == 2
    assert actors[0].mapper.dataset.n_points == 2


def test_draw_aabb_draws_the_twelve_edges():
    helper = new_helper()
    actor, = helper.draw_aabb(a_box())
    assert actor.prop.style == "Wireframe"
    assert actor.bounds == pytest.approx((0.0, 10.0, 0.0, 10.0, 0.0, 10.0), abs=1e-6)


def test_draw_aabb_accepts_anything_carrying_a_box():
    """ Drawing the extent of a mesh should not need the box unwrapped out of it first. """
    helper = new_helper()
    actor, = helper.draw_aabb(Mesh3.create_box(2.0, 4.0, 6.0, True))
    assert actor.bounds == pytest.approx((-1.0, 1.0, -2.0, 2.0, -3.0, 3.0), abs=1e-6)


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


def test_draw_distance_stands_the_label_off_against_the_scene():
    """
    Every dimension on a part gets the same standoff whatever it measures, rather than one
    proportional to its own value, which used to throw a long measurement's label clear off screen.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 60.0, 20.0, True))
    expected = numpy.linalg.norm([100.0, 60.0, 20.0]) * 0.1

    for value in (2.0, 100.0):
        distance = Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, value), Vector3(0.0, 0.0, 1.0))
        _, leaders, _ = helper.draw_distance(distance)
        # The leader that runs out to the label starts at `b` and ends at the label position.
        points = numpy.asarray(leaders.mapper.dataset.points)
        standoff = float(numpy.linalg.norm(points[-1] - numpy.array([0.0, 0.0, value])))
        assert standoff == pytest.approx(expected, rel=1e-6)


def test_an_inside_label_is_sized_by_the_measurement_not_the_scene():
    """
    An inside label belongs between the two anchor points, so a standoff taken from the scene would
    put it outside the very span it is meant to sit within. It is sized against the measurement
    instead, which makes it independent of everything else in the scene, unlike an outside label.
    """
    distance = Distance3(Point3(0.0, 0.0, -1.0), Point3(0.0, 0.0, 1.0), Vector3(0.0, 0.0, 1.0))

    def leader_points(box: float, label_place: str) -> numpy.ndarray:
        helper = new_helper()
        helper.draw_mesh(Mesh3.create_box(box, box, box, True))
        _, leaders, _ = helper.draw_distance(distance, label_place=label_place)
        return numpy.asarray(leaders.mapper.dataset.points)

    assert leader_points(10.0, "inside") == pytest.approx(leader_points(1000.0, "inside"), abs=TOL)
    assert numpy.abs(leader_points(10.0, "outside")
                     - leader_points(1000.0, "outside")).max() > 1.0


def test_drawing_a_dimension_does_not_grow_the_scene_for_the_next_one():
    """
    The leaders and label are annotation rather than geometry, so they stay out of the bounds. If
    they did not, each dimension would enlarge the scene that the following one measures its own
    standoff against, and the standoffs would compound.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 60.0, 20.0, True))
    distance = Distance3(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 5.0), Vector3(0.0, 0.0, 1.0))

    helper.draw_distance(distance)
    settled = helper._scene_diagonal()
    for _ in range(4):
        helper.draw_distance(distance)
    assert helper._scene_diagonal() == pytest.approx(settled, rel=1e-9)


def test_draw_distance_survives_a_zero_length_measurement():
    """
    An inside offset taken from the value would be zero for a zero-length measurement, leaving the
    label on top of the anchors, so it falls back to the scene in that case.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    distance = Distance3(Point3(1.0, 0.0, 0.0), Point3(1.0, 0.0, 0.0), Vector3(0.0, 0.0, 1.0))
    assert len(helper.draw_distance(distance, label_place="inside")) == 3


# ---------------------------------------------------------------------------
# Coordinate systems, labels, and arrows
# ---------------------------------------------------------------------------

def test_draw_coordinate_system_adds_one_actor_per_axis():
    helper = new_helper()
    assert len(helper.draw_coordinate_system(Iso3.identity(), length=1.0)) == 3


def test_draw_coordinate_system_adds_a_label_when_asked():
    helper = new_helper()
    assert len(helper.draw_coordinate_system(Iso3.identity(), length=1.0, text="A")) == 4


def test_draw_coordinate_system_accepts_several_frames():
    helper = new_helper()
    frames = helper.draw_coordinate_system(Iso3.identity(), Iso3.from_translation(1, 0, 0),
                                           length=1.0)
    assert len(frames) == 6


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
    actors = helper.draw_coordinate_system(Iso3.from_translation(10.0, 0.0, 0.0), length=2.0)
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

def test_line_builder_joins_added_points_into_one_polyline():
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    builder.add((1.0, 0.0, 0.0))
    builder.add((2.0, 0.0, 0.0))
    data = builder.build()
    assert data.n_points == 3
    assert data.n_cells == 1


def test_line_builder_skip_starts_a_new_run():
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    builder.add((1.0, 0.0, 0.0))
    builder.skip()
    builder.add((5.0, 0.0, 0.0))
    builder.add((6.0, 0.0, 0.0))
    data = builder.build()
    assert data.n_cells == 2
    assert numpy.asarray(data.points)[2] == pytest.approx([5.0, 0.0, 0.0], abs=TOL)


def test_line_builder_trailing_skip_adds_no_empty_run():
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    builder.add((1.0, 0.0, 0.0))
    builder.skip()
    assert builder.build().n_cells == 1


def test_line_builder_raises_when_there_is_nothing_to_draw():
    builder = LineBuilder()
    builder.add((0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        builder.build()


def test_lines_polydata_packs_several_runs_into_one_dataset():
    a = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    b = numpy.array([[0.0, 5.0, 0.0], [1.0, 5.0, 0.0]])
    data = _lines_polydata(a, b)
    assert data.n_points == 5
    assert data.n_cells == 2
    assert numpy.asarray(data.points)[3] == pytest.approx([0.0, 5.0, 0.0], abs=TOL)


def test_lines_polydata_skips_runs_too_short_to_draw():
    data = _lines_polydata(numpy.zeros((1, 3)), numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    assert data.n_cells == 1


# ---------------------------------------------------------------------------
# The `plotter.engeom` accessor
# ---------------------------------------------------------------------------

# The accessor is attached to `pyvista.BasePlotter` at import, so a test that detaches it has to
# put it back or every later test in the session loses it.
@pytest.fixture
def accessor_restored():
    yield
    engeom_pyvista.register()


def test_the_accessor_is_attached_when_the_module_is_imported():
    assert engeom_pyvista._attached()


def test_the_accessor_gives_a_helper_wrapping_that_plotter():
    """
    The point of the accessor is that a plotter and the helper that draws onto it are one object,
    rather than two similar ones sitting side by side in the same scope.
    """
    plotter = pyvista.Plotter(off_screen=True)
    assert isinstance(plotter.engeom, PlotterHelper)
    assert plotter.engeom.pv is plotter


def test_the_accessor_is_built_once_per_plotter():
    first = pyvista.Plotter(off_screen=True)
    second = pyvista.Plotter(off_screen=True)
    assert first.engeom is first.engeom
    assert first.engeom is not second.engeom
    assert second.engeom.pv is second


def test_the_accessor_is_discoverable_on_the_plotter():
    """ Autocomplete is most of the value, so the name has to show up in `dir`. """
    assert "engeom" in dir(pyvista.Plotter(off_screen=True))


def test_drawing_through_the_accessor_adds_to_the_plotter():
    plotter = pyvista.Plotter(off_screen=True)
    plotter.engeom.draw_mesh(a_mesh())
    assert len(plotter.renderer.actors) == 1


def test_the_accessor_draws_into_whichever_renderer_is_active():
    """
    One helper serves the whole plotter, so in a subplot layout it has to follow the active
    renderer rather than being bound to the one that happened to be current when it was built.
    """
    plotter = pyvista.Plotter(off_screen=True, shape=(1, 2))
    before = [len(r.actors) for r in plotter.renderers]

    plotter.subplot(0, 1)
    plotter.engeom.draw_mesh(a_mesh())
    plotter.subplot(0, 0)
    plotter.engeom.draw_mesh(a_mesh())
    plotter.engeom.draw_mesh(a_mesh())

    after = [len(r.actors) for r in plotter.renderers]
    assert [b - a for a, b in zip(before, after)] == [2, 1]


def test_the_helper_can_still_be_built_directly():
    """ The accessor is an addition, not a replacement. """
    plotter = pyvista.Plotter(off_screen=True)
    assert PlotterHelper(plotter).pv is plotter


def test_the_accessor_can_be_detached_and_put_back(accessor_restored):
    """
    Attaching to a third-party class is worth doing only if it can be undone, so the way out is
    part of the API rather than something to be reverse engineered.
    """
    assert engeom_pyvista.unregister()
    assert not hasattr(pyvista.Plotter(off_screen=True), "engeom")
    assert not engeom_pyvista.unregister(), "detaching twice should report that there was nothing to do"

    assert engeom_pyvista.register()
    assert isinstance(pyvista.Plotter(off_screen=True).engeom, PlotterHelper)


def test_registering_twice_does_not_warn(accessor_restored):
    """
    Re-registering is how the accessor is restored, so it must not be reported as an accidental
    collision with another plugin.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        engeom_pyvista.register()


def test_a_helper_already_built_survives_detaching(accessor_restored):
    plotter = pyvista.Plotter(off_screen=True)
    helper = plotter.engeom
    engeom_pyvista.unregister()
    helper.draw_mesh(a_mesh())
    assert len(plotter.renderer.actors) == 1


def test_the_typed_plotter_declares_the_accessor_without_shadowing_it():
    """
    The subclass exists so that editors and type checkers can see the accessor. The declaration has
    to stay an annotation: assigning to it would shadow the descriptor that does the work.
    """
    assert "engeom" not in vars(EngeomPlotter)
    assert EngeomPlotter.__annotations__["engeom"] is PlotterHelper

    plotter = EngeomPlotter(off_screen=True)
    assert isinstance(plotter.engeom, PlotterHelper)
    assert plotter.engeom.pv is plotter


def test_the_accessor_reaches_every_kind_of_plotter():
    """
    Registering against `BasePlotter` rather than `Plotter` is what puts the accessor on the Qt and
    background plotters too, which are not importable here but do subclass it.
    """
    assert "engeom" in vars(pyvista.BasePlotter)


def test_the_accessor_resolves_without_importing_engeom_first():
    """
    Declared as a `pyvista.plotter_components` entry point, so PyVista imports the helper the first
    time the attribute is looked up. This is what makes `import pyvista` on its own enough.

    Skipped on a source tree whose metadata predates the entry point, since the installed
    `dist-info` is what PyVista reads, not the working tree.
    """
    entries = importlib.metadata.entry_points(group="pyvista.plotter_components")
    if not any(e.name == "engeom" for e in entries):
        pytest.skip("the installed metadata declares no entry point; re-run `maturin develop`")

    script = textwrap.dedent("""
        import sys
        import pyvista
        pyvista.OFF_SCREEN = True
        assert not [m for m in sys.modules if m.startswith("engeom")], "engeom imported too early"
        helper = pyvista.Plotter(off_screen=True).engeom
        assert type(helper).__name__ == "PlotterHelper", type(helper)
        assert "engeom.plot.pyvista" in sys.modules
        print("ok")
    """)
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "ok" in result.stdout


# ---------------------------------------------------------------------------
# Interaction: widgets and pickers
# ---------------------------------------------------------------------------

# Nothing here drives the mouse. What is tested is the wiring either side of it: that the entity
# handed to a callback is the right `engeom` type in the right place, and that a selection made in
# the window can be mapped back onto the mesh it came from.

def test_the_plane_widget_reports_an_engeom_plane():
    """
    PyVista reports a moved plane as a pair of loose tuples, and the widget fires once on creation,
    so whatever the callback draws is present before the first interaction.
    """
    seen = []
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    helper.plane_widget(seen.append, normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 2.0))

    assert len(seen) == 1
    plane = seen[0]
    assert isinstance(plane, Plane3)
    assert numpy.asarray(plane.normal.as_numpy()) == pytest.approx([0.0, 0.0, 1.0], abs=1e-6)
    assert plane.signed_distance_to_point(Point3(0.0, 0.0, 2.0)) == pytest.approx(0.0, abs=1e-6)


def test_the_plane_widget_sections_a_mesh_live():
    """
    The reason the widget exists: section in the callback, draw under a fixed actor name, and each
    move replaces the last section rather than piling them up.
    """
    mesh = Mesh3.create_box(10.0, 10.0, 10.0, True)
    helper = new_helper()
    helper.draw_mesh(mesh, opacity=0.5, name="part")

    def show_section(plane):
        helper.draw_curve(*mesh.section_with_plane(plane).curves, color="black", name="section")

    helper.plane_widget(show_section)
    assert "section" in helper.pv.renderer.actors

    before = actor_count(helper)
    show_section(Plane3.from_point_normal(0.0, 0.0, 1.0, 0.0, 0.0, 1.0))
    assert actor_count(helper) == before


def test_the_plane_widget_needs_something_to_span():
    helper = new_helper()
    with pytest.raises(ValueError, match="nothing in the scene"):
        helper.plane_widget(lambda plane: None)


def test_picking_a_surface_point_reports_the_normal_as_well(monkeypatch):
    """
    A bare coordinate is much less useful than a `SurfacePoint3`: the normal is what makes a picked
    point usable as an alignment seed or a measurement anchor.
    """
    captured = {}
    monkeypatch.setattr(pyvista.Plotter, "enable_surface_point_picking",
                        lambda self, callback=None, **kwargs: captured.update(callback=callback))

    seen = []
    helper = new_helper()
    helper.pick_surface_point(seen.append)

    class FakePicker:
        def GetPickNormal(self):
            return (0.0, 0.0, 1.0)

    captured["callback"]((1.0, 2.0, 3.0), FakePicker())
    assert len(seen) == 1
    assert isinstance(seen[0], SurfacePoint3)
    assert numpy.asarray(seen[0].point.as_numpy()) == pytest.approx([1.0, 2.0, 3.0], abs=1e-9)
    assert numpy.asarray(seen[0].normal.as_numpy()) == pytest.approx([0.0, 0.0, 1.0], abs=1e-9)


def _capture_cell_picking(monkeypatch) -> dict:
    captured = {}
    monkeypatch.setattr(pyvista.Plotter, "enable_cell_picking",
                        lambda self, callback=None, **kwargs: captured.update(callback=callback))
    return captured


def test_picking_faces_maps_the_selection_back_onto_the_mesh(monkeypatch):
    """
    The selection PyVista returns has been through an extraction that renumbers the faces and drops
    VTK's own original-cell-ids, so the stamped `FACE_ID` array is what relates it back to the
    `engeom` mesh.
    """
    captured = _capture_cell_picking(monkeypatch)
    mesh = a_mesh()
    seen = []
    helper = new_helper()
    helper.pick_faces(mesh, seen.append)

    wanted = [2, 3, 7]
    captured["callback"](to_polydata(mesh).extract_cells(wanted))

    assert len(seen) == 1
    mask = seen[0]
    assert isinstance(mask, IndexMask)
    assert len(mask) == mesh.face_count
    assert mask.to_indices().tolist() == wanted


def test_picking_faces_across_several_actors_keeps_only_what_it_can_map(monkeypatch):
    """
    One rubber-band selection can cross several actors, which PyVista returns as a `MultiBlock`.
    Only the blocks carrying face ids came from a mesh drawn through this helper.
    """
    captured = _capture_cell_picking(monkeypatch)
    mesh = a_mesh()
    seen = []
    helper = new_helper()
    helper.pick_faces(mesh, seen.append)

    stranger = pyvista.Sphere()          # never went through `to_polydata`, so carries no ids
    captured["callback"](pyvista.MultiBlock([to_polydata(mesh).extract_cells([4, 5]), stranger]))

    assert seen[0].to_indices().tolist() == [4, 5]


def test_picking_nothing_reports_an_empty_mask(monkeypatch):
    """ A selection that lands on nothing should not need a special case at the call site. """
    captured = _capture_cell_picking(monkeypatch)
    mesh = a_mesh()
    seen = []
    helper = new_helper()
    helper.pick_faces(mesh, seen.append)

    captured["callback"](pyvista.MultiBlock([]))
    assert seen[0].count_true == 0
    assert len(seen[0]) == mesh.face_count


def test_a_picked_mask_goes_straight_back_into_the_mesh(monkeypatch):
    """
    The point of returning an `IndexMask` rather than a list of indices: it is the same currency
    the mesh's own operations take.
    """
    captured = _capture_cell_picking(monkeypatch)
    mesh = a_mesh()
    seen = []
    helper = new_helper()
    helper.pick_faces(mesh, seen.append)
    captured["callback"](to_polydata(mesh).extract_cells([0, 1, 2, 3]))

    mask = seen[0]
    assert mesh.extract_subset_faces(mask).face_count == 4
    assert helper.draw_mesh(mesh, highlight=mask) is not None


def test_closing_the_plotter_releases_the_interaction_it_switched_on():
    """
    The component's close hook runs on every close, and a plotter can be closed more than once, so
    it has to be safe to run again with nothing left to do.
    """
    helper = new_helper()
    helper.draw_mesh(a_mesh())
    helper.plane_widget(lambda plane: None)
    helper.pick_surface_point(lambda point: None)
    assert helper._widgets_added and helper._picking_enabled

    helper.__plotter_close__()
    assert not helper._widgets_added and not helper._picking_enabled
    helper.__plotter_close__()


def test_the_close_hook_runs_when_the_plotter_closes():
    """ The hook is only useful if PyVista's component lifecycle actually calls it. """
    plotter = pyvista.Plotter(off_screen=True)
    plotter.engeom.draw_mesh(a_mesh())
    plotter.engeom.pick_surface_point(lambda point: None)
    assert plotter.engeom._picking_enabled

    plotter.close()
    assert not plotter.engeom._picking_enabled


# ---------------------------------------------------------------------------
# The camera
# ---------------------------------------------------------------------------

def a_posed_plotter() -> PlotterHelper:
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    helper.pv.camera_position = [(30.0, -40.0, 25.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    helper.pv.enable_parallel_projection()
    return helper


def test_a_camera_pose_is_a_rigid_transform():
    """
    VTK does not keep the camera's up vector square to its view direction, so a frame built from
    the raw values is not a rotation and `Iso3` refuses it outright.
    """
    helper = a_posed_plotter()
    helper.pv.camera.up = (1.0, 0.0, 0.5)   # deliberately not square to the view direction

    matrix = helper.camera_pose().as_numpy()
    rotation = matrix[:3, :3]
    assert numpy.linalg.det(rotation) == pytest.approx(1.0, abs=1e-9)
    assert rotation @ rotation.T == pytest.approx(numpy.eye(3), abs=1e-9)


def test_a_camera_pose_describes_what_the_render_window_shows():
    """
    The pose is only worth having if its X and Y are the image's X and Y. A frame with +X right and
    +Y up must have +Z coming towards the viewer, and getting that backwards would leave the pose a
    valid rotation whose projection is the mirror of the render, which no round-trip test would
    notice.
    """
    helper = a_posed_plotter()
    helper.pv.render()
    pose = helper.camera_pose()

    probes = numpy.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0],
                          [-3.0, 2.0, 4.0], [1.0, -2.0, 3.0]])
    renderer = helper.pv.renderer
    screen = []
    for point in probes:
        renderer.SetWorldPoint(*point, 1.0)
        renderer.WorldToDisplay()
        screen.append(renderer.GetDisplayPoint()[:2])
    screen = numpy.array(screen)
    projected = numpy.array([(pose @ Point3(*p)).as_numpy()[:2] for p in probes])

    # Screen position must be the projected position under one positive scale plus an offset. A
    # negative term on either diagonal would be a mirrored axis.
    fitted, *_ = numpy.linalg.lstsq(numpy.hstack([projected, numpy.ones((len(probes), 1))]),
                                    screen, rcond=None)
    scale = fitted[:2].T
    assert numpy.abs(scale - numpy.diag(numpy.diag(scale))).max() < 1e-6
    assert scale[0, 0] > 0.0 and scale[1, 1] > 0.0
    assert scale[0, 0] == pytest.approx(scale[1, 1], rel=1e-9)


def test_a_camera_pose_round_trips_through_view_pose():
    helper = a_posed_plotter()
    pose = helper.camera_pose()
    helper.view_pose(pose)
    assert helper.camera_pose().as_numpy() == pytest.approx(pose.as_numpy(), abs=1e-9)


def test_view_pose_leaves_the_camera_distance_alone():
    """
    A pose carries orientation and what is being looked at, not how far off the camera sits, so
    applying one to a scene of another size must not throw the geometry out of frame.
    """
    helper = a_posed_plotter()
    pose = helper.camera_pose()
    before = numpy.linalg.norm(numpy.array(helper.pv.camera.position)
                               - numpy.array(helper.pv.camera.focal_point))

    helper.pv.camera_position = [(1.0, 1.0, 1.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    helper.view_pose(pose)
    after = numpy.linalg.norm(numpy.array(helper.pv.camera.position)
                              - numpy.array(helper.pv.camera.focal_point))
    assert after == pytest.approx(numpy.sqrt(3.0), abs=1e-9)
    assert after != pytest.approx(before, abs=1e-6)


def test_view_from_puts_the_camera_on_the_given_side():
    """
    The argument is where you are looking from, not the way the camera points, so `(0, 0, 1)` is a
    view from above.
    """
    helper = a_posed_plotter()
    helper.view_from((0.0, 0.0, 1.0))
    position = numpy.array(helper.pv.camera.position)
    focal = numpy.array(helper.pv.camera.focal_point)
    assert position[2] > focal[2]
    assert position[:2] == pytest.approx(focal[:2], abs=1e-6)


def test_view_from_can_switch_the_projection():
    helper = new_helper()
    helper.draw_mesh(a_mesh())
    helper.view_from((1.0, 0.0, 0.0), parallel=True)
    assert helper.pv.camera.parallel_projection
    helper.view_from((1.0, 0.0, 0.0), parallel=False)
    assert not helper.pv.camera.parallel_projection


def test_view_from_takes_an_up_direction():
    helper = a_posed_plotter()
    helper.view_from((0.0, -1.0, 0.0), up=(0.0, 0.0, 1.0))
    assert numpy.array(helper.pv.camera.up) == pytest.approx([0.0, 0.0, 1.0], abs=1e-6)


def test_fit_to_frames_only_what_it_is_given():
    """
    Framing one feature of a large part is the point, so the camera has to end up looking at that
    feature rather than at the whole scene.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 100.0, 100.0, True))
    feature = Aabb3(40.0, 40.0, 40.0, 45.0, 45.0, 45.0)

    helper.fit_to(feature)
    assert numpy.array(helper.pv.camera.focal_point) == pytest.approx([42.5, 42.5, 42.5], abs=1e-6)


def test_fit_to_merges_several_regions():
    helper = new_helper()
    helper.draw_mesh(a_mesh())
    helper.fit_to(Aabb3(0.0, 0.0, 0.0, 1.0, 1.0, 1.0), Aabb3(9.0, 9.0, 9.0, 10.0, 10.0, 10.0))
    assert numpy.array(helper.pv.camera.focal_point) == pytest.approx([5.0, 5.0, 5.0], abs=1e-6)


def test_fit_to_survives_the_next_read_of_the_camera():
    """
    PyVista defers its own framing until something marks the camera as set, so a plain
    `reset_camera(bounds=...)` is quietly undone by the next read of `plotter.camera`, which is
    what `show` and `screenshot` both do. The framing has to still be there when the scene is
    finally drawn, or `fit_to` would appear to work and then not.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 100.0, 100.0, True))
    helper.fit_to(Aabb3(40.0, 40.0, 40.0, 45.0, 45.0, 45.0))

    for _ in range(3):
        assert numpy.array(helper.pv.camera.focal_point) == pytest.approx([42.5, 42.5, 42.5],
                                                                          abs=1e-6)
    helper.pv.render()
    assert numpy.array(helper.pv.camera.focal_point) == pytest.approx([42.5, 42.5, 42.5], abs=1e-6)


def test_fit_to_with_nothing_frames_the_whole_scene():
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(10.0, 10.0, 10.0, True))
    helper.pv.camera_position = [(500.0, 500.0, 500.0), (400.0, 400.0, 400.0), (0.0, 0.0, 1.0)]
    helper.fit_to()
    assert numpy.array(helper.pv.camera.focal_point) == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)


# ---------------------------------------------------------------------------
# Conversions
# ---------------------------------------------------------------------------

def test_to_polydata_carries_the_mesh_triangles_across():
    mesh = a_mesh()
    data = to_polydata(mesh)
    assert data.n_cells == mesh.faces.shape[0]
    assert data.n_points == mesh.points.shape[0]


def test_to_polydata_stamps_element_ids_without_making_them_active_scalars():
    """
    An active scalar array is what `add_mesh` reaches for when given no color of its own, so
    stamping these by assignment would silently color every mesh drawn by its face index.
    """
    data = to_polydata(a_mesh())
    assert FACE_ID in data.cell_data
    assert POINT_ID in data.point_data
    assert data.active_scalars_name is None


def test_stamped_ids_survive_a_filter_that_renumbers_elements():
    """
    The ids are what make a selection made in the render window, or the output of any PyVista
    filter, relatable back to the original `engeom` entity.
    """
    kept = [3, 7, 8]
    extracted = to_polydata(a_mesh()).extract_cells(kept)
    assert numpy.asarray(extracted.cell_data[FACE_ID]).tolist() == kept


@pytest.mark.parametrize("build,cells,points", [
    (lambda: a_curve(), 1, 3),
    (lambda: CurveGroup3([a_curve(), a_curve()]), 2, 6),
    # A cloud and a bare array become vertex cells, one per point, which is what makes the
    # individual points pickable in the render window.
    (lambda: a_cloud(), 3, 3),
    (lambda: Segment3(0.0, 0.0, 0.0, 1.0, 0.0, 0.0), 1, 2),
    (lambda: numpy.zeros((7, 3)), 7, 7),
])
def test_to_polydata_accepts_every_documented_entity(build, cells, points):
    data = to_polydata(build())
    assert data.n_points == points
    assert data.n_cells == cells


def test_to_polydata_rejects_something_it_cannot_convert():
    with pytest.raises(TypeError):
        to_polydata("not an entity")


def test_to_mesh3_round_trips_a_mesh():
    mesh = a_mesh()
    back = to_mesh3(to_polydata(mesh))
    assert back.face_count == mesh.face_count
    assert back.point_count == mesh.point_count
    assert numpy.asarray(back.points) == pytest.approx(numpy.asarray(mesh.points), abs=TOL)


def test_to_mesh3_triangulates_a_dataset_that_is_not_already_triangles():
    """
    A `Mesh3` is a triangle mesh, so anything coming back from PyVista has to be triangulated,
    which is why the stamped ids rather than the face order are the way to relate the two.
    """
    box = pyvista.Box()
    assert box.n_cells == 6
    assert to_mesh3(box).face_count == 12


def test_to_mesh3_extracts_a_surface_from_a_volume():
    assert to_mesh3(pyvista.ImageData(dimensions=(3, 3, 3))).face_count > 0


# ---------------------------------------------------------------------------
# Coloring meshes and clouds
# ---------------------------------------------------------------------------

def test_draw_mesh_colors_by_scalars():
    mesh = a_mesh()
    actor = new_helper().draw_mesh(mesh, scalars=numpy.arange(mesh.face_count, dtype=float))
    assert actor.mapper.scalar_map_mode == "cell"


def test_draw_mesh_colors_by_point_scalars():
    mesh = a_mesh()
    actor = new_helper().draw_mesh(mesh, scalars=numpy.arange(mesh.point_count, dtype=float))
    assert actor.mapper.scalar_map_mode == "point"


def test_draw_mesh_highlights_a_face_selection():
    """
    The highlight is per-face color on the one actor, so it cannot z-fight with the surface it sits
    on the way a second actor drawn over the selection would.
    """
    mesh = a_mesh()
    selection = mesh.face_select().facing(0.0, 0.0, 1.0, 0.1, "add").to_mask()
    assert 0 < selection.count_true < mesh.face_count

    actor = new_helper().draw_mesh(mesh, color="white", highlight=selection,
                                   highlight_color="red")
    colors = numpy.asarray(actor.mapper.dataset.cell_data[actor.mapper.array_name])
    assert colors.shape == (mesh.face_count, 3)

    picked = selection.to_bool_array()
    assert (colors[picked] == numpy.array([255, 0, 0], dtype=numpy.uint8)).all()
    assert (colors[~picked] == numpy.array([255, 255, 255], dtype=numpy.uint8)).all()


def test_draw_mesh_highlights_points_when_the_mask_covers_points():
    """
    A mask knows how many elements it was built against, which is what tells a selection of faces
    apart from one of points without the caller having to say.
    """
    mesh = a_mesh()
    actor = new_helper().draw_mesh(mesh, highlight=mesh.point_mask(True))
    assert actor.mapper.scalar_map_mode == "point"


def test_draw_mesh_highlights_faces_for_a_bare_array_of_indices():
    actor = new_helper().draw_mesh(a_mesh(), highlight=[0, 2, 4])
    assert actor.mapper.scalar_map_mode == "cell"


def test_draw_mesh_rejects_a_highlight_that_fits_neither_faces_nor_points():
    with pytest.raises(ValueError, match="matches neither"):
        new_helper().draw_mesh(a_mesh(), highlight=numpy.zeros(7, dtype=bool))


def test_draw_mesh_rejects_scalars_and_a_highlight_together():
    mesh = a_mesh()
    with pytest.raises(ValueError, match="not both"):
        new_helper().draw_mesh(mesh, scalars=numpy.zeros(mesh.face_count), highlight=[0])


def test_draw_mesh_uses_the_meshs_own_stored_colors():
    mesh = a_mesh()
    colors = numpy.zeros((mesh.point_count, 3), dtype=numpy.uint8)
    colors[:, 1] = 255
    mesh.set_point_colors(colors)
    actor = new_helper().draw_mesh(mesh)
    assert actor.mapper.scalar_visibility


def test_an_explicit_color_wins_over_stored_colors():
    """
    PyVista will not take a color and a set of scalars at once, so one has to give way. An explicit
    argument is the more specific of the two, so the stored colors are the ones that yield.
    """
    mesh = a_mesh()
    mesh.set_point_colors(numpy.zeros((mesh.point_count, 3), dtype=numpy.uint8))
    actor = new_helper().draw_mesh(mesh, color="red")
    assert not actor.mapper.scalar_visibility


def test_draw_mesh_poses_the_actor_without_moving_the_mesh():
    """
    Posing through the renderer is what lets the same mesh be drawn at several candidate poses
    without copying its points, and lets a pose be changed afterwards on the actor.
    """
    mesh = a_mesh()
    before = numpy.asarray(mesh.points).copy()
    actor = new_helper().draw_mesh(mesh, pose=Iso3.from_translation(100.0, 0.0, 0.0))

    assert numpy.asarray(actor.user_matrix)[0, 3] == pytest.approx(100.0, abs=TOL)
    assert actor.bounds[0] > 90.0
    assert numpy.asarray(mesh.points) == pytest.approx(before, abs=TOL)


def test_draw_point_cloud_highlights_a_point_selection():
    cloud = a_cloud()
    actor = new_helper().draw_point_cloud(cloud, color="white", highlight=[1])[0]
    colors = numpy.asarray(actor.mapper.dataset.point_data[actor.mapper.array_name])
    assert (colors[1] == numpy.array([255, 0, 0], dtype=numpy.uint8)).all()
    assert (colors[0] == numpy.array([255, 255, 255], dtype=numpy.uint8)).all()


def test_draw_point_cloud_colors_by_scalars():
    cloud = a_cloud()
    actor = new_helper().draw_point_cloud(cloud, scalars=numpy.arange(3, dtype=float))[0]
    assert actor.mapper.scalar_map_mode == "point"


def test_draw_point_colors_by_scalars():
    actor = new_helper().draw_point(numpy.zeros((4, 3)), scalars=numpy.arange(4, dtype=float))
    assert actor.mapper.scalar_visibility


# ---------------------------------------------------------------------------
# Feature edges
# ---------------------------------------------------------------------------

def test_draw_feature_edges_finds_the_creases_of_a_box():
    """
    A box has twelve sharp creases, which is what makes it read as a box rather than as a blob.
    """
    actors = new_helper().draw_feature_edges(a_mesh())
    assert len(actors) == 1
    assert actors[0].mapper.dataset.n_cells == 12


def test_draw_feature_edges_finds_the_boundary_of_an_open_mesh():
    mesh = Mesh3.create_circle(1.0, 1.0e-3)
    edges = new_helper().draw_feature_edges(mesh, angle=None, non_manifold=False)
    assert edges[0].mapper.dataset.n_cells > 0


def test_draw_feature_edges_accepts_several_meshes():
    assert len(new_helper().draw_feature_edges(a_mesh(), a_mesh())) == 2


# ---------------------------------------------------------------------------
# Shared argument handling
# ---------------------------------------------------------------------------

def test_draw_point_accepts_an_array_of_points():
    """
    An `(n, 3)` array is the form the rest of the library hands point sets back in, so it is told
    apart from a loose coordinate by its shape rather than needing a separate method.
    """
    helper = new_helper()
    actor = helper.draw_point(numpy.zeros((5, 3)))
    assert actor.mapper.dataset.n_points == 5


def test_draw_point_accepts_a_point_cloud():
    helper = new_helper()
    actor = helper.draw_point(a_cloud())
    assert actor.mapper.dataset.n_points == 3


def test_draw_point_mixes_every_accepted_form_into_one_actor():
    helper = new_helper()
    actor = helper.draw_point(Point3(1.0, 1.0, 1.0), numpy.zeros((4, 3)), a_cloud(),
                              (2.0, 2.0, 2.0))
    assert actor.mapper.dataset.n_points == 1 + 4 + 3 + 1


def test_draw_point_rejects_being_given_nothing():
    with pytest.raises(ValueError):
        new_helper().draw_point()


def test_an_actor_name_is_suffixed_when_several_entities_are_drawn():
    """
    PyVista replaces any existing actor of the same name, so reusing one name across a varargs call
    would leave only the last entity on the plotter.
    """
    helper = new_helper()
    helper.draw_curve(a_curve(), a_curve(), a_curve(), name="section")
    assert actor_count(helper) == 3
    assert set(helper.pv.renderer.actors) == {"section_0", "section_1", "section_2"}


def test_an_actor_name_is_left_alone_for_a_single_entity():
    helper = new_helper()
    helper.draw_curve(a_curve(), name="section")
    assert set(helper.pv.renderer.actors) == {"section"}


def test_a_legend_label_is_applied_to_the_first_entity_only():
    """
    One legend entry per call, not one per entity, so that drawing a computed group of curves under
    a single label does not fill the legend with identical rows.
    """
    helper = new_helper()
    helper.draw_curve(a_curve(), a_curve(), a_curve(), label="nominal")
    assert len(helper.pv.renderer._labels) == 1


@pytest.mark.parametrize("draw", [
    lambda h: h.draw_mesh(a_mesh()),
    lambda h: h.draw_curve(a_curve()),
    lambda h: h.draw_point((0.0, 0.0, 0.0)),
    lambda h: h.draw_sphere(Sphere3(0.0, 0.0, 0.0, 1.0)),
    lambda h: h.draw_circle(Circle3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)),
    lambda h: h.draw_point_cloud(a_cloud()),
    lambda h: h.draw_arrow((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
])
def test_every_draw_method_works_with_no_styling_arguments_at_all(draw):
    """
    A styling argument left as `None` must not be forwarded, since PyVista's defaults are not all
    fixed values and passing `None` through is not the same as omitting it.
    """
    assert draw(new_helper()) is not None


def test_a_coordinate_system_sizes_itself_to_the_scene():
    """
    A frame carries no size of its own, so an unspecified length is taken from what has already
    been drawn, which keeps it legible whether the scene is millimeters or meters across.
    """
    helper = new_helper()
    helper.draw_mesh(Mesh3.create_box(100.0, 100.0, 100.0, True))
    axes = helper.draw_coordinate_system(Iso3.identity())

    ends = numpy.asarray(axes[0].mapper.dataset.points)
    drawn = float(numpy.linalg.norm(ends[1] - ends[0]))
    assert drawn == pytest.approx(numpy.linalg.norm([100.0, 100.0, 100.0]) * 0.1, rel=1e-6)


def test_a_coordinate_system_falls_back_to_unit_length_on_an_empty_scene():
    """
    PyVista reports a two-unit cube for an empty renderer, which is a placeholder rather than a
    measurement, so it must not be mistaken for a scene to scale against.
    """
    helper = new_helper()
    axes = helper.draw_coordinate_system(Iso3.identity())
    ends = numpy.asarray(axes[0].mapper.dataset.points)
    assert float(numpy.linalg.norm(ends[1] - ends[0])) == pytest.approx(1.0, abs=TOL)


# ---------------------------------------------------------------------------
# Coverage drift
# ---------------------------------------------------------------------------

# Every public method on the helper, each of which must be exercised above. Adding a method without
# adding a test will fail this, in the spirit of test_stub_drift.py.
EXERCISED_HELPER_MEMBERS = {
    "draw_arrow", "draw_circle", "draw_coordinate_system", "draw_curve", "draw_distance",
    "camera_pose", "fit_to", "pick_faces", "pick_surface_point", "plane_widget", "view_from",
    "view_pose",
    "draw_aabb", "draw_basis", "draw_cone", "draw_cylinder", "draw_feature_edges", "draw_label",
    "draw_line", "draw_mesh", "draw_plane", "draw_point", "draw_point_cloud", "draw_segment",
    "draw_sphere", "draw_spline", "draw_surface_point",
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
    assert not_draw == {"camera_pose", "fit_to", "pick_faces", "pick_surface_point",
                        "plane_widget", "view_from", "view_pose"}
