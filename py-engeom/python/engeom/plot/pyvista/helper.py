"""
The `PlotterHelper` itself, which wraps a PyVista `Plotter` and draws `engeom` entities onto it.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy
import pyvista
from pyvista import ColorLike

from engeom.common import IndexMask
from engeom.geom3 import Circle3, Cone3, CubicSpline3, Curve3, CurveGroup3, Cylinder3, \
    Iso3, Line3, Mesh3, Plane3, PointCloud3, Segment3, Sphere3, SurfacePoint3, SvdBasis3
from engeom.metrology import Distance3

from .._coerce import Coords3, to_tuple3
from .._common import LabelPlace, check_label_place, plane_basis
from .._style import ElementStyle, element_style, merge_style
from .convert import FACE_ID, LineBuilder, _lines_polydata, to_polydata
from .extent import clip_line_to_aabb, clip_plane_to_aabb, resolve_extent

__all__ = ["PlotterHelper"]

# When a draw method's chordal tolerance is left as `None`, it defaults to this fraction of the
# entity's radius. Perceived smoothness in a plot scales with the entity's size on screen, so this
# radius-relative default keeps entities equally smooth without requiring the caller to calculate
# it.
_DEFAULT_REL_TOL = 1.0e-3

# The default width, in pixels, of any line the helper draws. PyVista's own default is 1.0, which
# is hard to pick out against a rendered surface, so the lines drawn here are wider by default.
# Passing `line_width=None` defers to PyVista and gets the thin line back.
_DEFAULT_LINE_WIDTH = 3.0

# How far a dimension's label sits from the midpoint when it is placed inside the measurement, as
# a fraction of the measured value. An inside label belongs between the two anchor points, so it is
# sized against the measurement itself rather than against the scene.
_INSIDE_LABEL_FRACTION = 0.25

# The length of a point cloud's normal arrows, as a fraction of the diagonal of the cloud itself.
# Scaling to the cloud rather than to the whole scene keeps the arrows readable when the cloud is a
# small patch of a much larger part.
_DEFAULT_REL_ARROW = 0.05

# How transparent a plane is by default. A plane is drawn to show where it cuts through
# something, so an opaque one would hide the very thing it was drawn to explain.
_DEFAULT_PLANE_OPACITY = 0.5

# The length of an entity drawn with no size of its own, such as a coordinate frame, as a fraction
# of the diagonal of everything drawn so far. A frame sized to the scene stays legible whether the
# scene is a 2mm feature or a 2m part.
_DEFAULT_REL_LENGTH = 0.1


def _segments_for_tol(radius: float, tol: float) -> int:
    """
    Return the number of segments needed so that a chord inscribed on a full circle of `radius`
    deviates by no more than `tol`, with a minimum of 8 segments. This is a small counterpart to
    the Rust `arc_segments_for_tol` function for polylines built directly in NumPy.
    """
    half_angle = numpy.arcsin(min(1.0, numpy.sqrt(tol / (2.0 * radius))))
    if half_angle <= 0.0:
        raise ValueError(f"Tolerance {tol} is not usable with radius {radius}")
    return max(8, int(numpy.ceil(numpy.pi / (2.0 * half_angle))))


def _bounds_center(box) -> tuple[float, float, float]:
    """The center of a bounding box, as the plain tuple PyVista's widgets expect."""
    center = box.center
    return (center.x, center.y, center.z)


def _selection_length(selection) -> int | None:
    """
    The number of elements a selection was built against, or `None` if it carries no length.

    An `IndexMask` and a boolean array both know how many elements they cover, which is what tells
    a selection over a mesh's faces apart from one over its points. A plain array of indices
    carries no such information.
    """
    if isinstance(selection, IndexMask):
        return len(selection)
    array = numpy.asarray(selection)
    return int(array.size) if array.dtype == numpy.bool_ else None


def _selection_mask(selection, count: int) -> numpy.ndarray:
    """
    Resolve a selection of elements into a boolean mask over `count` elements.

    :param selection: an `IndexMask`, a boolean array, or an array of indices.
    :param count: the number of elements the resulting mask covers.
    :return: a boolean array of length `count`.
    :raises ValueError: if a mask or a boolean array does not cover exactly `count` elements.
    """
    if isinstance(selection, IndexMask):
        selection = selection.to_bool_array()

    array = numpy.asarray(selection)
    if array.dtype == numpy.bool_:
        if array.size != count:
            raise ValueError(f"A selection over {array.size} elements cannot be applied to "
                             f"{count} elements")
        return array

    mask = numpy.zeros(count, dtype=bool)
    mask[array.astype(numpy.int64)] = True
    return mask


def _highlight_scalars(mask: numpy.ndarray, base: ColorLike,
                       highlight: ColorLike) -> numpy.ndarray:
    """
    Build the per-element colors that draw a highlighted subset of an entity.

    Coloring the elements of the one actor is what keeps a highlight from z-fighting with the
    entity underneath it, which a second actor drawn over the selection would, and keeps a scalar
    bar out of it, which encoding the selection as a scalar range would not.

    :param mask: a boolean array selecting the elements to highlight.
    :param base: the color of the elements that are not selected.
    :param highlight: the color of the elements that are.
    :return: an `(n, 3)` array of `uint8` colors, one per element.
    """
    colors = numpy.empty((mask.size, 3), dtype=numpy.uint8)
    colors[:] = pyvista.Color(base).int_rgb
    colors[mask] = pyvista.Color(highlight).int_rgb
    return colors


def _color_scalars(scalars) -> bool:
    """
    Whether an array of scalars holds colors rather than values to be mapped through a color map.

    An `(n, 3)` or `(n, 4)` array of `uint8` is unambiguously RGB or RGBA, which is the form the
    stored colors on a mesh or a point cloud take, and the form a highlight is built in. Floating
    point arrays of the same shape are left alone, since those are as likely to be vectors.
    """
    return (isinstance(scalars, numpy.ndarray) and scalars.ndim == 2
            and scalars.shape[1] in (3, 4) and scalars.dtype == numpy.uint8)


def _per_entity(kwargs: Dict[str, Any], index: int, count: int) -> Dict[str, Any]:
    """
    Adjust the keyword arguments shared by a varargs draw method for one entity of several.

    Two of PyVista's arguments cannot simply be repeated. An actor `name` replaces any existing
    actor of the same name, so reusing one across a call would leave only the last entity drawn; a
    suffix is added to keep them distinct. A legend `label` would produce one identical entry per
    entity, so only the first entity carries it and the group gets a single legend entry.

    :param kwargs: the keyword arguments to adjust, which are not mutated.
    :param index: the index of the entity being drawn.
    :param count: how many entities the call is drawing in total.
    :return: the keyword arguments to draw this entity with.
    """
    if count == 1:
        return kwargs
    working = dict(kwargs)
    if "name" in working:
        working["name"] = f"{working['name']}_{index}"
    if index > 0:
        working.pop("label", None)
    return working


class PlotterHelper:
    """
    A helper that wraps a PyVista `Plotter` and provides methods for drawing `engeom` entities.

    Every method that adds an actor is prefixed `draw_`, matching the Matplotlib helper, so
    `helper.draw_<tab>` lists everything that can be drawn. The wrapped plotter stays public as
    `helper.pv`, for anything the helper does not cover.

    Each `draw_` method names the styling arguments that get used constantly, and forwards an open
    `**kwargs` on to the PyVista call named in its docstring, so anything PyVista accepts can still
    be passed. A named styling argument left as `None` is not forwarded at all, so PyVista's own
    default applies. The one deliberate departure is `line_width`, which defaults to 3.0 rather
    than PyVista's 1.0, because a one-pixel line is hard to pick out in a rendered 3D scene.

    Methods that draw one kind of entity are named in the singular and take that entity as varargs,
    returning one actor per entity in the order given. `draw_mesh` and `draw_point_cloud` are the
    exceptions, taking a single entity, because their coloring arguments describe one entity's
    per-point or per-face data and could not be shared across several.

    !!! example
        ```python
        import pyvista
        plotter = pyvista.Plotter()
        helper = PlotterHelper(plotter)
        ```
    """

    def __init__(self, plotter: pyvista.Plotter):
        """
        Initialize the helper with a PyVista `Plotter`.

        :param plotter: the PyVista `Plotter` to wrap.
        """
        self.pv = plotter
        self._picking_enabled = False
        self._widgets_added = False

    def _scene_diagonal(self) -> float | None:
        """
        The diagonal of everything drawn into the active renderer so far, or `None` if there is
        nothing bounded to measure.

        PyVista reports a two-unit cube centered on the origin for a renderer holding no 3D actors,
        which is a placeholder rather than a measurement, so it is rejected here rather than
        silently scaling later entities against it.
        """
        if not self.pv.renderer.actors:
            return None

        bounds = tuple(float(v) for v in self.pv.bounds)
        if bounds == (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0):
            return None

        spans = [bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]]
        diagonal = float(numpy.linalg.norm(spans))
        return diagonal if diagonal > 0.0 else None

    def _relative_length(self, fraction: float = _DEFAULT_REL_LENGTH) -> float:
        """
        A length scaled to the scene, for entities that carry no size of their own.

        :param fraction: the fraction of the scene diagonal to return.
        :return: that fraction of the scene diagonal, or 1.0 if the scene has nothing to measure.
        """
        diagonal = self._scene_diagonal()
        return 1.0 if diagonal is None else diagonal * fraction

    def plane_widget(self, callback, extent=None, pad: float = 0.1, normal: Coords3 = (0.0, 0.0, 1.0),
                     origin: Coords3 | None = None, **kwargs):
        """
        Add a plane that can be dragged and rotated in the render window, reporting a `Plane3` as
        it moves.

        This is what turns a plot into an instrument. Sectioning a mesh in the callback and drawing
        the result under a fixed actor `name`, so each redraw replaces the last, gives a section
        that follows the plane live:

        ```python
        def show_section(plane):
            plotter.engeom.draw_curve(*mesh.section_with_plane(plane).curves, name="section")

        plotter.engeom.draw_mesh(mesh, opacity=0.5)
        plotter.engeom.plane_widget(show_section)
        ```

        :param callback: called with a `Plane3` whenever the widget moves, and once immediately so
            that whatever it draws is there before the first interaction.
        :param extent: the region the widget spans. `None` uses everything already drawn into the
            active renderer; otherwise an `Aabb3`, an entity carrying a bounding box, or an
            `(n, 3)` array of points.
        :param pad: how far to grow that region, as a fraction of its diagonal.
        :param normal: the plane's starting normal.
        :param origin: the plane's starting origin. If None, the center of the extent is used.
        :param kwargs: any other keyword argument accepted by PyVista's
            `Plotter.add_plane_widget`.
        :return: the VTK widget, for anything PyVista exposes on it directly.
        :raises ValueError: if the extent cannot be resolved.
        """
        box = resolve_extent(extent, self.pv, pad)
        bounds = (box.min.x, box.max.x, box.min.y, box.max.y, box.min.z, box.max.z)

        def on_move(widget_normal, widget_origin):
            callback(Plane3.from_point_normal(*to_tuple3(widget_origin), *to_tuple3(widget_normal)))

        kwargs.setdefault("bounds", bounds)
        widget = self.pv.add_plane_widget(
            on_move,
            normal=to_tuple3(normal),
            origin=_bounds_center(box) if origin is None else to_tuple3(origin),
            **kwargs,
        )
        self._widgets_added = True
        return widget

    def pick_surface_point(self, callback, **kwargs):
        """
        Report a `SurfacePoint3` wherever the user clicks on a surface.

        The point comes back with the surface normal at the point of contact, not just the
        position, which is what makes it usable as an alignment seed or a measurement anchor rather
        than only as a coordinate.

        :param callback: called with a `SurfacePoint3` for each click that lands on geometry.
        :param kwargs: any other keyword argument accepted by PyVista's
            `Plotter.enable_surface_point_picking`.
        """

        def on_pick(point, picker):
            normal = picker.GetPickNormal()
            callback(SurfacePoint3(*to_tuple3(point), *normal))

        kwargs.setdefault("show_message", False)
        self.pv.enable_surface_point_picking(callback=on_pick, use_picker=True, **kwargs)
        self._picking_enabled = True

    def pick_faces(self, mesh: Mesh3, callback, **kwargs):
        """
        Report the faces of a mesh that the user selects in the render window, as an `IndexMask`
        over that mesh's faces.

        The mask indexes the `engeom` mesh, not the PyVista dataset, so it can go straight back
        into `Mesh3.face_select`, `extract_subset_faces`, `section_with_plane`, or `draw_mesh`'s
        own `highlight`. That works because `to_polydata` records each face's original index on
        the dataset, and PyVista's filters carry cell data through the extraction a rubber-band
        selection performs; VTK's own original-cell-ids are discarded by it.

        The mesh has to have been drawn through this helper for the selection to be relatable back
        to it.

        :param mesh: the mesh whose faces are being selected, used for the length of the mask.
        :param callback: called with an `IndexMask` over the mesh's faces for each selection. A
            selection that lands on nothing calls it with an empty mask.
        :param kwargs: any other keyword argument accepted by PyVista's
            `Plotter.enable_cell_picking`.
        """

        def on_pick(picked):
            indices = []
            # Several actors can fall inside one rubber-band selection, and PyVista hands back a
            # `MultiBlock` when they do. Only the blocks carrying face ids came from a mesh drawn
            # through this helper, and only those can be mapped back.
            blocks = picked if isinstance(picked, pyvista.MultiBlock) else [picked]
            for block in blocks:
                if block is not None and FACE_ID in block.cell_data:
                    indices.append(numpy.asarray(block.cell_data[FACE_ID], dtype=numpy.int64))

            found = numpy.concatenate(indices) if indices else numpy.empty(0, dtype=numpy.int64)
            callback(IndexMask.from_indices(numpy.unique(found), mesh.face_count))

        kwargs.setdefault("show_message", False)
        self.pv.enable_cell_picking(callback=on_pick, **kwargs)
        self._picking_enabled = True

    def __plotter_close__(self) -> None:
        """
        Release the interaction the helper switched on, when the plotter closes.

        PyVista calls this on every close, and a plotter can be closed more than once, so it has to
        be safe to run again with nothing left to do.
        """
        if self._picking_enabled:
            self._picking_enabled = False
            self.pv.disable_picking()

        if self._widgets_added:
            self._widgets_added = False
            self.pv.clear_plane_widgets()

    def camera_pose(self) -> Iso3:
        """
        The current view, as an isometry transforming world coordinates into view coordinates.

        In view coordinates +X is to the right, +Y is up, and +Z points out of the image towards
        the viewer. Those three cannot be chosen freely: with +X right and +Y up, a rigid transform
        has to have +Z coming towards the viewer, since the other choice makes the frame
        left-handed and no longer a rotation at all. This is the same convention VTK and OpenGL
        use for a camera, so the pose describes what the render window is actually showing.

        Pass it to `view_pose` to come back to a view later, so that one found by dragging the
        window is reusable rather than something to be recreated by hand.

        It is also the convention `ViewPort3` uses in the Matplotlib backend, so a pose taken here
        can be handed straight to `AxesHelper.viewport` to draw a line diagram from the same angle
        the render window is showing.

        VTK does not keep the camera's up vector square to its view direction, so the up vector is
        orthogonalized against the direction before the frame is built. An `Iso3` is a rigid
        transform and refuses a matrix that is not a rotation, and a camera that has been dragged
        around will not have supplied one.

        The pose carries orientation and the point being looked at. It does not carry how far away
        the camera is or how much it has zoomed, since neither survives the trip through a parallel
        projection; use `fit_to` to frame the scene.

        :return: the world-to-view isometry.
        """
        camera = self.pv.camera
        focal = numpy.array(camera.focal_point, dtype=numpy.float64)
        direction = focal - numpy.array(camera.position, dtype=numpy.float64)
        direction = direction / numpy.linalg.norm(direction)

        up = numpy.array(camera.up, dtype=numpy.float64)
        up = up - direction * float(up @ direction)
        up = up / numpy.linalg.norm(up)

        matrix = numpy.eye(4)
        matrix[0, :3] = numpy.cross(direction, up)
        matrix[1, :3] = up
        matrix[2, :3] = -direction
        matrix[:3, 3] = -matrix[:3, :3] @ focal
        return Iso3(matrix)

    def view_pose(self, pose: Iso3) -> None:
        """
        Point the camera according to a world-to-view isometry, in the convention `camera_pose`
        returns and `ViewPort3` takes.

        Only the orientation and the point being looked at are set. The camera keeps its distance
        and its zoom, so a pose can be applied to a scene of a different size without the geometry
        flying out of frame; follow with `fit_to` to frame it.

        :param pose: the world-to-view isometry to adopt.
        """
        matrix = pose.as_numpy()
        rotation = matrix[:3, :3]

        # The world point that lands on the view origin is what the camera is looking at, and the
        # rows of a world-to-view rotation are the view axes expressed in world coordinates. The
        # third row points out of the image towards the viewer, so the camera looks the other way
        # along it.
        focal = -rotation.T @ matrix[:3, 3]
        direction = -rotation[2]
        up = rotation[1]

        camera = self.pv.camera
        distance = float(numpy.linalg.norm(
            numpy.array(camera.focal_point, dtype=numpy.float64)
            - numpy.array(camera.position, dtype=numpy.float64)))

        camera.focal_point = tuple(focal)
        camera.position = tuple(focal - direction * distance)
        camera.up = tuple(up)

    def view_from(self, direction: Coords3, up: Coords3 | None = None,
                  parallel: bool | None = None) -> None:
        """
        Look at the scene from a given direction.

        :param direction: the direction from the scene towards the camera, so ``(0, 0, 1)`` looks
            down from above and ``(1, 0, 0)`` looks back along the X axis. Any coordinate form the
            helpers accept.
        :param up: which way is up in the resulting image. If None, PyVista chooses one.
        :param parallel: whether to use a parallel projection, which is what makes equal lengths
            measure equally anywhere in the image. If None, whatever the plotter is already using
            is kept.
        """
        if parallel is True:
            self.pv.enable_parallel_projection()
        elif parallel is False:
            self.pv.disable_parallel_projection()

        self.pv.view_vector(to_tuple3(direction),
                            viewup=None if up is None else to_tuple3(up))

    def fit_to(self, *entities, pad: float = 0.0) -> None:
        """
        Frame the camera on particular geometry rather than on the whole scene.

        With nothing given this is `Plotter.reset_camera`, framing everything. Given entities, only
        those are framed, which is how to zoom in on one feature of a large part without hiding the
        rest of it.

        :param entities: what to frame. Each is an `Aabb3`, an entity carrying a bounding box
            (`Mesh3`, `Segment3`, `CurveGroup3`, `CubicSpline3`, `PointCloud3`), or an `(n, 3)`
            array of points.
        :param pad: how far to grow the framed region, as a fraction of its diagonal, to leave a
            margin around it.
        :raises TypeError: if something is given that carries no bounding box.
        """
        # Reading the camera is what settles any framing PyVista still has pending: until something
        # marks it as set, the next read of `plotter.camera` resets it to fit the whole scene, and
        # would throw away the framing done here. Taking the reference first, then framing, then
        # marking it set, is what makes the result stick.
        camera = self.pv.camera

        if not entities:
            self.pv.reset_camera()
        else:
            box = resolve_extent(entities[0], self.pv, pad)
            for entity in entities[1:]:
                box = box.merged(resolve_extent(entity, self.pv, pad))
            self.pv.reset_camera(bounds=(box.min.x, box.max.x, box.min.y, box.max.y,
                                         box.min.z, box.max.z))

        camera.is_set = True

    def draw_mesh(
            self,
            mesh: Mesh3,
            color: ColorLike | None = None,
            opacity: float | None = None,
            show_edges: bool | None = None,
            scalars: numpy.ndarray | str | None = None,
            cmap=None,
            clim: tuple[float, float] | None = None,
            show_scalar_bar: bool | None = None,
            categorical: bool = False,
            highlight=None,
            highlight_color: ColorLike = "red",
            use_colors: bool = True,
            pose: Iso3 | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> pyvista.Actor:
        """
        Add an `engeom` mesh to be plotted.

        Unlike most of the draw methods this takes a single mesh rather than varargs, because the
        arguments describing how to color it refer to that one mesh's points and faces.

        The color of the surface comes from whichever of `color`, `scalars`, `highlight`, or the
        mesh's own stored colors is given, in that order of specificity. `scalars` and `highlight`
        cannot both be given, since each sets the color of every element.

        :param mesh: the mesh to add to the plotter scene.
        :param color: a single color for the whole mesh. If None, the PyVista default applies,
            unless another of the color arguments supplies one.
        :param opacity: the opacity of the mesh, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param show_edges: whether to draw the mesh's triangle edges over its surface. If None, the
            PyVista default applies.
        :param scalars: a value per point or per face to color the mesh by, or the name of an array
            already on the dataset. Which of the two it indexes is inferred from its length. An
            `(n, 3)` or `(n, 4)` array of `uint8` is taken as RGB or RGBA colors rather than values
            to map.
        :param cmap: the color map for `scalars`, as a name or a Matplotlib `Colormap`. A color
            map's out-of-range colors are carried across to PyVista, so a map like `GOM_CMAP` keeps
            calling out data beyond the scale rather than clamping it.
        :param clim: the `(low, high)` limits of the color map. If None, the range of the data is
            used.
        :param show_scalar_bar: whether to draw the scalar bar. If None, the PyVista default
            applies.
        :param categorical: whether the scalars are category labels rather than a continuous
            quantity, which draws one flat color per distinct value. Use this for label arrays such
            as the mesh's face labels or the output of `compute_patch_labels`.
        :param highlight: a selection of faces or points to pick out in `highlight_color`, as an
            `IndexMask`, a boolean array, or an array of indices. A mask or boolean array is taken
            as indexing faces or points according to its length; a bare array of indices is taken
            as indexing faces. The whole mesh stays a single actor, so the highlight cannot
            z-fight with the surface under it.
        :param highlight_color: the color of the highlighted elements.
        :param use_colors: whether to use the mesh's own stored point or face colors when no other
            color argument is given. Point colors are preferred where the mesh carries both.
        :param pose: an isometry to draw the mesh at, applied by the renderer rather than by moving
            the mesh, so the same mesh can be drawn at several poses without copying its points.
            The actor's `user_matrix` can be reassigned afterwards to move it again.
        :param label: the legend label. If None, the mesh is not labeled.
        :param name: the actor name. Drawing again under the same name replaces the earlier actor,
            which is what makes a redraw in a callback or an animation loop update in place.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: the PyVista actor that was added to the plotter.
        :raises ValueError: if both `scalars` and `highlight` are given, or if a highlight covers a
            number of elements matching neither the mesh's faces nor its points.
        """
        if scalars is not None and highlight is not None:
            raise ValueError("A mesh takes either `scalars` or `highlight`, not both, since each "
                             "one sets the color of every element")

        if highlight is not None:
            length = _selection_length(highlight)
            if length is None or length == mesh.face_count:
                count = mesh.face_count
            elif length == mesh.point_count:
                count = mesh.point_count
            else:
                raise ValueError(
                    f"A highlight over {length} elements matches neither the mesh's "
                    f"{mesh.face_count} faces nor its {mesh.point_count} points")
            base = color if color is not None else pyvista.global_theme.color
            scalars = _highlight_scalars(_selection_mask(highlight, count), base, highlight_color)
            color = None
        elif scalars is None and color is None and use_colors:
            scalars = mesh.point_colors if mesh.point_colors is not None else mesh.face_colors

        if cmap is not None:
            kwargs.update(_cmap_extremes(cmap))
        if _color_scalars(scalars):
            kwargs.setdefault("rgb", True)
        if pose is not None:
            kwargs.setdefault("user_matrix", pose.as_numpy())
        if categorical:
            kwargs.setdefault("categories", True)

        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             scalars=scalars, cmap=cmap, clim=clim,
                             show_scalar_bar=show_scalar_bar, label=label, name=name)
        return self.pv.add_mesh(to_polydata(mesh), **kwargs)

    def draw_feature_edges(
            self,
            *meshes: Mesh3,
            angle: float | None = 30.0,
            boundary: bool = True,
            non_manifold: bool = True,
            manifold: bool = False,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            pose: Iso3 | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Draw the edges that describe a mesh's shape: the creases between faces meeting at more than
        `angle`, the open boundaries, and the edges where the surface is not manifold.

        This reads a mesh's form far better than drawing every triangle edge does, and it is the
        quickest way to see where a mesh has holes or non-manifold construction.

        :param meshes: the meshes whose edges to draw.
        :param angle: the angle, in degrees, beyond which the crease between two faces counts as a
            feature. If None, creases are left out and only the other kinds of edge are drawn.
        :param boundary: whether to include edges with only one face, which are the mesh's open
            boundaries.
        :param non_manifold: whether to include edges shared by more than two faces.
        :param manifold: whether to include every remaining interior edge, which draws the whole
            triangulation.
        :param color: the color of the edges. If None, the PyVista default applies.
        :param line_width: the pixel width of the edges. If None, the PyVista default applies.
        :param opacity: the opacity of the edges, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param pose: an isometry to draw the edges at, applied by the renderer rather than by
            moving the mesh.
        :param label: the legend label, applied to the first mesh only.
        :param name: the actor name, suffixed with the mesh's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per mesh, in the order given. A mesh with no edges of the kinds
            asked for still returns an actor, holding nothing.
        """
        if pose is not None:
            kwargs.setdefault("user_matrix", pose.as_numpy())
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             label=label, name=name)

        actors = []
        for i, mesh in enumerate(meshes):
            edges = to_polydata(mesh).extract_feature_edges(
                feature_angle=0.0 if angle is None else angle,
                boundary_edges=boundary,
                non_manifold_edges=non_manifold,
                feature_edges=angle is not None,
                manifold_edges=manifold,
            )
            actors.append(self.pv.add_mesh(edges, **_per_entity(kwargs, i, len(meshes))))
        return actors

    def draw_point_cloud(
            self,
            cloud: PointCloud3,
            use_colors: bool = True,
            normals: ElementStyle = False,
            color: ColorLike | None = None,
            point_size: float | None = None,
            opacity: float | None = None,
            render_points_as_spheres: bool = True,
            scalars: numpy.ndarray | str | None = None,
            cmap=None,
            clim: tuple[float, float] | None = None,
            show_scalar_bar: bool | None = None,
            highlight=None,
            highlight_color: ColorLike = "red",
            pose: Iso3 | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add a `PointCloud3` to the plotter, optionally colored by the cloud's own per-point colors
        and optionally overlaid with arrows showing the per-point normals.

        Unlike most of the draw methods this takes a single cloud rather than varargs, because the
        arguments describing how to color it refer to that one cloud's points.

        The color of the points comes from whichever of `color`, `scalars`, `highlight`, or the
        cloud's own stored colors is given, in that order of specificity.

        :param cloud: the point cloud to plot.
        :param use_colors: whether to use the cloud's own stored per-point colors when no other
            color argument is given.
        :param normals: whether to draw an arrow at each point along its normal, and how. `False`
            leaves them out, `True` draws them with the defaults, and a dict of keyword arguments
            for PyVista's `Plotter.add_arrows` is merged over those defaults. The arrow length is
            `mag`, defaulting to a twentieth of the cloud's own diagonal. Nothing is drawn if the
            cloud has no normals.
        :param color: a single color for every point. If None, the PyVista default applies, unless
            another of the color arguments supplies one.
        :param point_size: the size of the points. If None, the PyVista default applies.
        :param opacity: the opacity of the points, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param render_points_as_spheres: whether to render the points as spheres rather than flat
            squares.
        :param scalars: a value per point to color the cloud by, or the name of an array already on
            the dataset. An `(n, 3)` or `(n, 4)` array of `uint8` is taken as RGB or RGBA colors
            rather than values to map. Useful for the per-point results of a measurement, such as
            `Mesh3.measure_deviations` or an alignment's residuals.
        :param cmap: the color map for `scalars`, as a name or a Matplotlib `Colormap`.
        :param clim: the `(low, high)` limits of the color map. If None, the range of the data is
            used.
        :param show_scalar_bar: whether to draw the scalar bar. If None, the PyVista default
            applies.
        :param highlight: a selection of points to pick out in `highlight_color`, as an
            `IndexMask`, a boolean array, or an array of indices.
        :param highlight_color: the color of the highlighted points.
        :param pose: an isometry to draw the cloud at, applied by the renderer rather than by
            moving the points, so the same cloud can be drawn at several poses without copying it.
        :param label: the legend label. If None, the cloud is not labeled.
        :param name: the actor name. Drawing again under the same name replaces the earlier actor.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_points`.
        :return: the PyVista actors that were added to the plotter, the normal arrows first if they
            were drawn, then the points.
        :raises ValueError: if both `scalars` and `highlight` are given, or if a highlight does not
            cover the cloud's points.
        """
        if scalars is not None and highlight is not None:
            raise ValueError("A point cloud takes either `scalars` or `highlight`, not both, "
                             "since each one sets the color of every point")

        actors = []

        if normals is not False and cloud.point_normals is not None:
            diagonal = float(numpy.linalg.norm(cloud.compute_aabb().extent.as_numpy()))
            normal_style = element_style(normals, {
                "mag": diagonal * _DEFAULT_REL_ARROW if diagonal > 0.0 else 1.0,
                "color": color if color is not None else "gray",
                "reset_camera": False,
            })
            actors.append(self.pv.add_arrows(cloud.points, cloud.point_normals, **normal_style))

        if highlight is not None:
            base = color if color is not None else pyvista.global_theme.color
            # `PointCloud3.point_count` is a method rather than a property, unlike its `Mesh3`
            # counterpart, so the length is taken through `len` to keep the two the same here.
            scalars = _highlight_scalars(_selection_mask(highlight, len(cloud)), base,
                                         highlight_color)
            color = None
        elif scalars is None and color is None and use_colors:
            scalars = cloud.point_colors

        if cmap is not None:
            kwargs.update(_cmap_extremes(cmap))
        if _color_scalars(scalars):
            kwargs.setdefault("rgb", True)
        if pose is not None:
            kwargs.setdefault("user_matrix", pose.as_numpy())

        kwargs = merge_style(kwargs, color=color, point_size=point_size, opacity=opacity,
                             render_points_as_spheres=render_points_as_spheres, scalars=scalars,
                             cmap=cmap, clim=clim, show_scalar_bar=show_scalar_bar, label=label,
                             name=name)
        actors.append(self.pv.add_points(cloud.points, **kwargs))
        return actors

    def draw_point(
            self,
            *points: Coords3 | PointCloud3 | numpy.ndarray,
            color: ColorLike | None = None,
            point_size: float | None = None,
            render_points_as_spheres: bool = True,
            opacity: float | None = None,
            scalars: numpy.ndarray | None = None,
            cmap=None,
            clim: tuple[float, float] | None = None,
            show_scalar_bar: bool | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> pyvista.Actor:
        """
        Add one or more discrete points to be plotted.

        Every point given goes into a single actor, so this returns one actor rather than one per
        argument. Alongside loose coordinates, an `(n, 3)` array or a `PointCloud3` may be given,
        since those are the forms the rest of the library hands point sets back in; the three are
        told apart by their shape.

        :param points: the points to add, each a single coordinate in any form the helpers accept,
            an `(n, 3)` array of coordinates, or a `PointCloud3`.
        :param color: the color of the points. If None, the PyVista default applies.
        :param point_size: the size of the points. If None, the PyVista default applies.
        :param render_points_as_spheres: whether to render the points as spheres rather than flat
            squares.
        :param opacity: the opacity of the points, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param scalars: a value per point to color by, in the order the points were given. An
            `(n, 3)` or `(n, 4)` array of `uint8` is taken as RGB or RGBA colors instead.
        :param cmap: the color map for `scalars`, as a name or a Matplotlib `Colormap`.
        :param clim: the `(low, high)` limits of the color map. If None, the range of the data is
            used.
        :param show_scalar_bar: whether to draw the scalar bar. If None, the PyVista default
            applies.
        :param label: the legend label. If None, the points are not labeled.
        :param name: the actor name. Drawing again under the same name replaces the earlier actor.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_points`.
        :return: the single PyVista actor holding every point given.
        :raises ValueError: if no points were given.
        """
        collected = []
        for item in points:
            if isinstance(item, PointCloud3):
                collected.append(numpy.asarray(item.points, dtype=numpy.float64))
            elif numpy.ndim(item) == 2:
                collected.append(numpy.asarray(item, dtype=numpy.float64).reshape(-1, 3))
            else:
                collected.append(numpy.array([to_tuple3(item)], dtype=numpy.float64))

        if not collected:
            raise ValueError("No points were given, so there is nothing to draw")

        if cmap is not None:
            kwargs.update(_cmap_extremes(cmap))
        if _color_scalars(scalars):
            kwargs.setdefault("rgb", True)

        kwargs = merge_style(kwargs, color=color, point_size=point_size, opacity=opacity,
                             render_points_as_spheres=render_points_as_spheres, scalars=scalars,
                             cmap=cmap, clim=clim, show_scalar_bar=show_scalar_bar, label=label,
                             name=name)
        return self.pv.add_points(numpy.concatenate(collected), **kwargs)

    def draw_curve(
            self,
            *curves: Curve3 | CurveGroup3,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more curves to be plotted, each as its own actor.

        A `CurveGroup3` draws as a single actor holding all of its curves, since the group is one
        entity that moves as a rigid body, which is also the way to draw a large number of curves
        without an actor apiece.

        :param curves: the curves to add, each either a `Curve3` or a `CurveGroup3`.
        :param color: the color of the curves. If None, the PyVista default applies.
        :param line_width: the pixel width of the curves. If None, the PyVista default applies.
        :param opacity: the opacity of the curves, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param label: the legend label, applied to the first curve only, so that a group of curves
            drawn in one call produces a single legend entry.
        :param name: the actor name, suffixed with the curve's index when more than one is drawn.
            Drawing again under the same name replaces the earlier actor.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per curve, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             label=label, name=name)

        actors = []
        for i, curve in enumerate(curves):
            runs = [c.points for c in curve.curves] if isinstance(curve, CurveGroup3) \
                else [curve.points]
            actors.append(self.pv.add_mesh(_lines_polydata(*runs),
                                           **_per_entity(kwargs, i, len(curves))))
        return actors

    def draw_sphere(
            self,
            *spheres: Sphere3,
            tol: float | None = None,
            color: ColorLike | None = None,
            opacity: float | None = None,
            show_edges: bool | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more `Sphere3` to the plotter, each tessellated into its own actor.

        :param spheres: the spheres to plot.
        :param tol: the maximum deviation of the tessellation from the true sphere. Smaller values
            produce a smoother result. If None, a thousandth of each sphere's radius is used.
        :param color: the color of the spheres. If None, the PyVista default applies.
        :param opacity: the opacity of the spheres, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param show_edges: whether to draw the tessellation's edges over the surface. If None, the
            PyVista default applies.
        :param label: the legend label, applied to the first sphere only.
        :param name: the actor name, suffixed with the sphere's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per sphere, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             label=label, name=name)

        actors = []
        for i, sphere in enumerate(spheres):
            mesh = Mesh3.create_sphere(sphere.r, sphere.r * _DEFAULT_REL_TOL if tol is None else tol)
            mesh.transform_in_place(Iso3.from_translation(*sphere.center))
            actors.append(self.pv.add_mesh(to_polydata(mesh),
                                           **_per_entity(kwargs, i, len(spheres))))
        return actors

    def draw_circle(
            self,
            *circles: Circle3,
            tol: float | None = None,
            color: ColorLike | None = None,
            opacity: float | None = None,
            show_edges: bool = True,
            edge_color: ColorLike | None = None,
            edge_opacity: float | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            style: str | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more `Circle3` to the plotter, each as a single filled polygon with an outline.

        The face and the outline are one actor rather than two, so they cannot be drawn at
        different depths and cannot z-fight with each other. Pass ``style="wireframe"`` for an
        outline with no face, or ``show_edges=False`` for a face with no outline.

        :param circles: the circles to plot.
        :param tol: the maximum chordal deviation of the polygon from the true circle. If None, a
            thousandth of each circle's radius is used.
        :param color: the color of the face. If None, the PyVista default applies.
        :param opacity: the opacity of the face, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param show_edges: whether to draw the circle's outline over its face.
        :param edge_color: the color of the outline. If None, the PyVista default applies.
        :param edge_opacity: the opacity of the outline, from 0.0 to 1.0. If None, the PyVista
            default applies.
        :param line_width: the pixel width of the outline. If None, the PyVista default applies.
        :param style: the PyVista draw style: ``"surface"``, ``"wireframe"``, or ``"points"``. If
            None, the PyVista default applies. Use ``"wireframe"`` to draw the outline alone.
        :param label: the legend label, applied to the first circle only.
        :param name: the actor name, suffixed with the circle's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per circle, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             edge_color=edge_color, edge_opacity=edge_opacity,
                             line_width=line_width, style=style, label=label, name=name)

        actors = []
        for i, circle in enumerate(circles):
            # A `Circle3` carries only a center, a normal, and a radius, so an in-plane orientation
            # has to be chosen before it can be turned into vertices. Which one does not matter, as
            # the result is rotationally symmetric about the normal either way.
            n = _segments_for_tol(circle.r, circle.r * _DEFAULT_REL_TOL if tol is None else tol)
            u, v = plane_basis(circle.normal)
            center = numpy.array(to_tuple3(circle.center), dtype=numpy.float64)
            u_r = numpy.array(to_tuple3(u), dtype=numpy.float64) * circle.r
            v_r = numpy.array(to_tuple3(v), dtype=numpy.float64) * circle.r

            theta = numpy.linspace(0.0, 2.0 * numpy.pi, n, endpoint=False)
            points = center + numpy.outer(numpy.cos(theta), u_r) + numpy.outer(numpy.sin(theta), v_r)
            faces = numpy.concatenate([[n], numpy.arange(n)])

            actors.append(self.pv.add_mesh(pyvista.PolyData(points, faces=faces),
                                           **_per_entity(kwargs, i, len(circles))))
        return actors

    def draw_cylinder(
            self,
            *cylinders: Cylinder3,
            tol: float | None = None,
            color: ColorLike | None = None,
            opacity: float | None = None,
            show_edges: bool | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more `Cylinder3` to the plotter, each tessellated into its own actor.

        :param cylinders: the cylinders to draw.
        :param tol: the maximum deviation of the tessellation from the true surface. If None, a
            thousandth of each cylinder's radius is used.
        :param color: the color of the cylinders. If None, the PyVista default applies.
        :param opacity: the opacity, from 0.0 to 1.0. If None, the PyVista default applies.
        :param show_edges: whether to draw the tessellation's edges. If None, the PyVista default
            applies.
        :param label: the legend label, applied to the first cylinder only.
        :param name: the actor name, suffixed with the cylinder's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per cylinder, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             label=label, name=name)
        actors = []
        for i, cylinder in enumerate(cylinders):
            mesh = Mesh3.create_cylinder_between(
                cylinder.a, cylinder.b, cylinder.r,
                cylinder.r * _DEFAULT_REL_TOL if tol is None else tol)
            actors.append(self.pv.add_mesh(to_polydata(mesh),
                                           **_per_entity(kwargs, i, len(cylinders))))
        return actors

    def draw_cone(
            self,
            *cones: Cone3,
            tol: float | None = None,
            color: ColorLike | None = None,
            opacity: float | None = None,
            show_edges: bool | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more `Cone3` to the plotter, each tessellated into its own actor.

        :param cones: the cones to draw.
        :param tol: the maximum deviation of the tessellation from the true surface. If None, a
            thousandth of each cone's radius is used.
        :param color: the color of the cones. If None, the PyVista default applies.
        :param opacity: the opacity, from 0.0 to 1.0. If None, the PyVista default applies.
        :param show_edges: whether to draw the tessellation's edges. If None, the PyVista default
            applies.
        :param label: the legend label, applied to the first cone only.
        :param name: the actor name, suffixed with the cone's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per cone, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             label=label, name=name)
        actors = []
        for i, cone in enumerate(cones):
            mesh = Mesh3.create_cone(cone.r, cone.height,
                                     cone.r * _DEFAULT_REL_TOL if tol is None else tol)

            # `create_cone` builds a cone centered on the origin with its apex towards +Z, so the
            # local +Z axis is the direction from the base back to the tip, and the local origin is
            # the midpoint of the axis rather than either end.
            # The basis has to be taken of the axis itself and not of the cone's direction: the
            # two point opposite ways, and a frame built from the wrong one is left-handed, which
            # is not a rotation and will not convert into an `Iso3`.
            towards_tip = -cone.direction
            u, v = plane_basis(towards_tip)
            axis = numpy.array(to_tuple3(towards_tip), dtype=numpy.float64)
            axis = axis / numpy.linalg.norm(axis)
            tip = numpy.array(to_tuple3(cone.tip), dtype=numpy.float64)

            matrix = numpy.eye(4)
            matrix[:3, 0] = to_tuple3(u)
            matrix[:3, 1] = to_tuple3(v)
            matrix[:3, 2] = axis
            matrix[:3, 3] = tip - axis * (cone.height / 2.0)
            mesh.transform_in_place(Iso3(matrix))

            actors.append(self.pv.add_mesh(to_polydata(mesh),
                                           **_per_entity(kwargs, i, len(cones))))
        return actors

    def draw_spline(
            self,
            *splines: CubicSpline3,
            tol: float | None = None,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more cubic splines to the plotter, each flattened into a polyline.

        :param splines: the splines to draw.
        :param tol: the maximum deviation of the polyline from the true curve. If None, a
            thousandth of the diagonal of each spline's own bounding box is used, which keeps a
            spline equally smooth whatever its size.
        :param color: the color of the splines. If None, the PyVista default applies.
        :param line_width: the pixel width. If None, the PyVista default applies.
        :param opacity: the opacity, from 0.0 to 1.0. If None, the PyVista default applies.
        :param label: the legend label, applied to the first spline only.
        :param name: the actor name, suffixed with the spline's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per spline, in the order given.
        :raises ValueError: if a spline has no extent for a default tolerance to be derived from.
        """
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             label=label, name=name)
        actors = []
        for i, spline in enumerate(splines):
            if tol is None:
                diagonal = float(numpy.linalg.norm(spline.aabb.extent.as_numpy()))
                if diagonal <= 0.0:
                    raise ValueError(
                        f"Spline {i} has no extent, so no tolerance can be derived from it. Give "
                        "`tol` explicitly."
                    )
                spline_tol = diagonal * _DEFAULT_REL_TOL
            else:
                spline_tol = tol
            actors.append(self.pv.add_mesh(_lines_polydata(spline.polyline(spline_tol)),
                                           **_per_entity(kwargs, i, len(splines))))
        return actors

    def draw_surface_point(
            self,
            *points: SurfacePoint3,
            length: float | None = None,
            color: ColorLike | None = None,
            point_size: float | None = None,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more surface points, each drawn as a point with an arrow along its normal.

        Every point given goes into one actor and every arrow into another, rather than an actor
        apiece, so that a large set of them restyles as one thing.

        :param points: the surface points to draw.
        :param length: the length of the normal arrows. If None, a tenth of the diagonal of
            everything drawn so far is used, or 1.0 if nothing bounded has been drawn yet.
        :param color: the color of the points and arrows. If None, the PyVista default applies.
        :param point_size: the size of the points. If None, the PyVista default applies.
        :param opacity: the opacity, from 0.0 to 1.0. If None, the PyVista default applies.
        :param label: the legend label, applied to the points.
        :param name: the actor name. The arrows take the same name with `_normals` appended.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_points`.
        :return: two PyVista actors, the points and then the arrows, or an empty list if no points
            were given.
        """
        if not points:
            return []

        if length is None:
            length = self._relative_length()

        positions = numpy.array([to_tuple3(p.point) for p in points], dtype=numpy.float64)
        normals = numpy.array([to_tuple3(p.normal) for p in points], dtype=numpy.float64)

        marker_kwargs = merge_style(dict(kwargs), color=color, point_size=point_size,
                                    opacity=opacity, label=label, name=name)
        actors = [self.pv.add_points(positions, render_points_as_spheres=True, **marker_kwargs)]

        arrow_kwargs = merge_style({}, color=color, opacity=opacity,
                                   name=None if name is None else f"{name}_normals")
        actors.append(self.pv.add_arrows(positions, normals, mag=length, reset_camera=False,
                                         **arrow_kwargs))
        return actors

    def draw_basis(
            self,
            *bases: SvdBasis3,
            scale: float = 1.0,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more `SvdBasis3` to the plotter, as three axis lines whose lengths show how far
        the data spreads along each principal direction.

        This is a coordinate frame whose axes are scaled rather than equal: the first axis is the
        direction of greatest variation and the third the least, so a basis fitted to a flat patch
        of points draws as a nearly flat cross rather than as a cube's corner.

        :param bases: the bases to draw.
        :param scale: multiplies each axis length, which is one standard deviation of the data
            along that direction. The default of 1.0 draws one sigma.
        :param line_width: the pixel width of the axis lines. If None, the PyVista default applies.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: three PyVista actors per basis, X then Y then Z, in the order given.
        """
        kwargs = merge_style(kwargs, line_width=line_width)
        actors = []
        for basis in bases:
            iso = basis.to_iso3()
            lengths = numpy.asarray(basis.basis_stdevs, dtype=numpy.float64) * scale
            origin = numpy.array([[0.0, 0.0, 0.0]], dtype=numpy.float64)
            for axis, color in enumerate(("red", "green", "blue")):
                tip = numpy.zeros((1, 3))
                tip[0, axis] = lengths[axis]
                ends = iso.transform_points(numpy.concatenate([origin, tip]))
                actors.append(self.pv.add_mesh(_lines_polydata(ends), color=color, **kwargs))
        return actors

    def draw_segment(
            self,
            *segments: Segment3,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more line segments to the plotter.

        :param segments: the segments to draw.
        :param color: the color of the segments. If None, the PyVista default applies.
        :param line_width: the pixel width of the segments. If None, the PyVista default applies.
        :param opacity: the opacity of the segments, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param label: the legend label, applied to the first segment only.
        :param name: the actor name, suffixed with the segment's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per segment, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             label=label, name=name)
        actors = []
        for i, segment in enumerate(segments):
            points = numpy.array([to_tuple3(segment.a), to_tuple3(segment.b)], dtype=numpy.float64)
            actors.append(self.pv.add_mesh(_lines_polydata(points),
                                           **_per_entity(kwargs, i, len(segments))))
        return actors

    def draw_aabb(
            self,
            *boxes,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            style: str | None = "wireframe",
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more axis-aligned bounding boxes to the plotter, as wireframes by default.

        Anything carrying a bounding box may be given as well as an `Aabb3` itself, so the extent
        of a mesh or a point cloud can be drawn without unwrapping it first.

        :param boxes: the boxes to draw. Each is an `Aabb3`, an entity carrying one (`Mesh3`,
            `Segment3`, `CurveGroup3`, `CubicSpline3`, `PointCloud3`), or an `(n, 3)` array of
            points to take the box of.
        :param color: the color of the boxes. If None, the PyVista default applies.
        :param line_width: the pixel width of the box edges. If None, the PyVista default applies.
        :param opacity: the opacity of the boxes, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param style: the PyVista draw style. The default of ``"wireframe"`` draws the twelve edges;
            pass ``"surface"`` for solid faces, which is usually worth combining with an opacity.
        :param label: the legend label, applied to the first box only.
        :param name: the actor name, suffixed with the box's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per box, in the order given.
        :raises TypeError: if something is given that carries no bounding box.
        """
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             style=style, label=label, name=name)
        actors = []
        for i, value in enumerate(boxes):
            box = resolve_extent(value, pad=0.0)
            data = pyvista.Box(bounds=(box.min.x, box.max.x, box.min.y, box.max.y,
                                       box.min.z, box.max.z))
            actors.append(self.pv.add_mesh(data, **_per_entity(kwargs, i, len(boxes))))
        return actors

    def draw_line(
            self,
            *lines: Line3,
            extent=None,
            pad: float = 0.1,
            color: ColorLike | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more infinite lines to the plotter, each drawn across a finite extent.

        A `Line3` has no ends, so how much of it to show has to come from somewhere. By default
        that is the scene: draw the geometry the line relates to first, and the line is clipped to
        a box around it.

        :param lines: the lines to draw.
        :param extent: the region to draw the lines across. `None` uses everything already drawn
            into the active renderer; otherwise an `Aabb3`, an entity carrying a bounding box, or
            an `(n, 3)` array of points. For a specific span instead, draw a `Segment3` built from
            `Line3.at`.
        :param pad: how far to grow the extent, as a fraction of its diagonal, so a line crossing a
            part runs past it rather than stopping flush with its surface.
        :param color: the color of the lines. If None, the PyVista default applies.
        :param line_width: the pixel width of the lines. If None, the PyVista default applies.
        :param opacity: the opacity of the lines, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param label: the legend label, applied to the first line only.
        :param name: the actor name, suffixed with the line's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per line, in the order given.
        :raises ValueError: if the extent cannot be resolved, or if a line misses it entirely, in
            which case the line's position is named so it can be identified.
        """
        if not lines:
            return []

        box = resolve_extent(extent, self.pv, pad)
        kwargs = merge_style(kwargs, color=color, line_width=line_width, opacity=opacity,
                             label=label, name=name)

        actors = []
        for i, line in enumerate(lines):
            span = clip_line_to_aabb(line, box)
            if span is None:
                raise ValueError(
                    f"Line {i} does not pass through the extent, so there is nothing to draw. "
                    "Give an `extent` that the line crosses."
                )
            actors.append(self.pv.add_mesh(_lines_polydata(span),
                                           **_per_entity(kwargs, i, len(lines))))
        return actors

    def draw_plane(
            self,
            *planes: Plane3,
            extent=None,
            pad: float = 0.1,
            color: ColorLike | None = None,
            opacity: float | None = _DEFAULT_PLANE_OPACITY,
            show_edges: bool = True,
            edge_color: ColorLike | None = None,
            edge_opacity: float | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            style: str | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more infinite planes to the plotter, each drawn as the polygon where it crosses
        a finite extent.

        A `Plane3` carries no origin, no in-plane orientation, and no size, so all three come from
        the extent: the plane is clipped against that box, which is what makes the drawn polygon
        follow the shape of the region rather than being an arbitrary square floating in it.

        Planes default to partly transparent, since an opaque one hides whatever it cuts through,
        which is usually the thing being looked at. Overlapping transparent surfaces are drawn in
        whatever order they were added unless `Plotter.enable_depth_peeling` is switched on.

        :param planes: the planes to draw.
        :param extent: the region to draw the planes across. `None` uses everything already drawn
            into the active renderer; otherwise an `Aabb3`, an entity carrying a bounding box, or
            an `(n, 3)` array of points.
        :param pad: how far to grow the extent, as a fraction of its diagonal, so a plane cutting
            through a part protrudes past it rather than stopping flush with its surface.
        :param color: the color of the plane's face. If None, the PyVista default applies.
        :param opacity: the opacity of the face, from 0.0 to 1.0.
        :param show_edges: whether to outline the polygon where the plane meets the extent.
        :param edge_color: the color of that outline. If None, the PyVista default applies.
        :param edge_opacity: the opacity of the outline. If None, the PyVista default applies.
        :param line_width: the pixel width of the outline. If None, the PyVista default applies.
        :param style: the PyVista draw style. Pass ``"wireframe"`` for the outline alone.
        :param label: the legend label, applied to the first plane only.
        :param name: the actor name, suffixed with the plane's index when more than one is drawn.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: one PyVista actor per plane, in the order given.
        :raises ValueError: if the extent cannot be resolved, or if a plane misses it, which
            includes a plane meeting it at a single corner or edge and so having no area to draw.
        """
        if not planes:
            return []

        box = resolve_extent(extent, self.pv, pad)
        kwargs = merge_style(kwargs, color=color, opacity=opacity, show_edges=show_edges,
                             edge_color=edge_color, edge_opacity=edge_opacity,
                             line_width=line_width, style=style, label=label, name=name)

        actors = []
        for i, plane in enumerate(planes):
            polygon = clip_plane_to_aabb(plane, box)
            if polygon is None:
                raise ValueError(
                    f"Plane {i} does not cross the extent, so there is no polygon to draw. "
                    "Give an `extent` the plane passes through."
                )
            faces = numpy.concatenate([[len(polygon)], numpy.arange(len(polygon))])
            actors.append(self.pv.add_mesh(pyvista.PolyData(polygon, faces=faces),
                                           **_per_entity(kwargs, i, len(planes))))
        return actors

    def draw_distance(
            self,
            distance: Distance3,
            template: str = "{value:.3f}",
            label_place: LabelPlace = "outside",
            label_offset: float | None = None,
            font_size: int = 12,
            font_family: str | None = None,
            scale_value: float = 1.0,
            color: ColorLike | None = "black",
            line_width: float | None = 1.5,
            point_size: float | None = 4.0,
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add a distance entity to the plotter, as a pair of anchor points, leader lines running out
        from them, and a text label carrying the measured value.

        :param distance: the distance entity to add.
        :param template: a format string for the label. The `value` key is replaced with the actual
            value read from the measurement.
        :param label_place: the placement of the label relative to the distance entity's anchor
            points. See `LabelPlace`.
        :param label_offset: how far the label sits from its anchor, which means different things
            depending on `label_place`. If None, the outside placements stand the label off by a
            tenth of the diagonal of everything drawn so far, so that every dimension on a part
            gets the same standoff whatever it measures, while `"inside"` offsets it from the
            midpoint by a quarter of the measured value, so that it stays between the two anchors.
            The leader lines and the label are left out of the scene's bounds, so drawing one
            dimension does not change the standoff of the next.
        :param font_size: the size of the text to use for the label.
        :param font_family: the font family to use for the label. If None, the PyVista default
            applies.
        :param scale_value: a scaling factor applied to the value before it is displayed in the
            label. Use this to convert between units of measurement without modifying the actual
            value or the coordinate system.
        :param color: the color of the anchor points and the leader lines. If None, the PyVista
            default applies.
        :param line_width: the pixel width of the leader lines. If None, the PyVista default
            applies.
        :param point_size: the size of the anchor points. If None, the PyVista default applies.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_point_labels`,
            which draws the label.
        :return: the PyVista actors that were added to the plotter: the anchor points, the leader
            lines, and the label.
        :raises ValueError: if `label_place` is not one of the valid placement tokens.
        """
        label_place = check_label_place(label_place)
        if label_offset is None:
            if label_place == "inside":
                label_offset = abs(distance.value) * _INSIDE_LABEL_FRACTION
                if label_offset == 0.0:
                    label_offset = self._relative_length()
            else:
                label_offset = self._relative_length()

        # The offset_dir is the direction from `a` to `b` projected so that it's parallel to the
        # measurement direction.
        offset_dir = distance.direction if distance.value >= 0 else -distance.direction

        # Rather than arrows, we'll use spheres to indicate the anchor points at the end of the
        # leader lines
        spheres = [distance.a, distance.b]
        builder = LineBuilder()

        if label_place == "inside":
            c = SurfacePoint3(*distance.center.point, *offset_dir)
            label_coords = c.at_distance(label_offset)

            builder.add(distance.a)
            builder.add(distance.a - offset_dir * label_offset * 0.25)
            builder.skip()

            builder.add(distance.b)
            builder.add(distance.b + offset_dir * label_offset * 0.25)

        elif label_place == "outside":
            label_coords = distance.b + offset_dir * label_offset

            builder.add(distance.a)
            builder.add(distance.a - offset_dir * label_offset * 0.25)
            builder.skip()

            builder.add(distance.b)
            builder.add(label_coords)

        else:  # "outside_rev", the only remaining token after validation
            label_coords = distance.a - offset_dir * label_offset

            builder.add(distance.b)
            builder.add(distance.b + offset_dir * label_offset * 0.25)
            builder.skip()

            builder.add(distance.a)
            builder.add(label_coords)

        points = numpy.array([to_tuple3(p) for p in spheres], dtype=numpy.float64)
        marker_style = merge_style({}, color=color, point_size=point_size)
        actors = [self.pv.add_points(points, render_points_as_spheres=True, **marker_style)]

        leader_style = merge_style({}, color=color, line_width=line_width)
        leaders = self.pv.add_mesh(builder.build(), **leader_style)
        actors.append(leaders)

        value = distance.value * scale_value
        kwargs = merge_style(kwargs, font_family=font_family)
        label = self.pv.add_point_labels(
            [to_tuple3(label_coords)],
            [template.format(value=value)],
            show_points=False,
            background_color="white",
            font_size=font_size,
            bold=False,
            **kwargs,
        )
        actors.append(label)

        # The leaders and the label are annotation standing off from the part rather than geometry
        # belonging to it, so they are kept out of the renderer's bounds. Without this a dimension
        # would enlarge the scene that the next one sizes its own standoff against, and a long
        # measurement would pull the camera back off the part it annotates. The anchor points are
        # left contributing, since those sit at positions actually being measured.
        leaders.use_bounds = False
        label.use_bounds = False

        return actors

    def draw_coordinate_system(
            self,
            *isos,
            length: float | None = None,
            line_width: float | None = _DEFAULT_LINE_WIDTH,
            text: str | None = None,
            font_size: int = 12,
            font_family: str | None = "courier",
            **kwargs,
    ) -> list[pyvista.Actor]:
        """
        Add one or more coordinate frames to the plotter. Each appears as three lines, with X in
        red, Y in green, and Z in blue.

        :param isos: the isometries to use as the origin and orientation of each coordinate frame.
            Each may be an `Iso3`, a 4x4 `numpy.ndarray` that validly converts into an `Iso3`, or
            anything with an `as_numpy` method returning a valid 4x4 `numpy.ndarray`.
        :param length: the length of each axis line. If None, a tenth of the diagonal of everything
            drawn so far is used, or 1.0 if nothing bounded has been drawn yet, so draw the geometry
            the frame belongs to before drawing the frame.
        :param line_width: the pixel width of the axis lines. If None, the PyVista default applies.
        :param text: a label to display at the origin of every frame drawn. If None, no label is
            drawn.
        :param font_size: the size of the label text.
        :param font_family: the font family of the label text. If None, the PyVista default applies.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`, applied
            to the axis lines.
        :return: the PyVista actors that were added to the plotter: three axis lines per frame,
            each followed by that frame's label if one was requested.
        :raises ValueError: if an `iso` is a `numpy.ndarray` whose shape is not (4, 4).
        :raises TypeError: if an `iso` is neither an `Iso3` nor something convertible to a 4x4
            `numpy.ndarray`.
        """
        if length is None:
            length = self._relative_length()

        kwargs = merge_style(kwargs, line_width=line_width)

        actors = []
        for iso in isos:
            if not isinstance(iso, Iso3):
                if hasattr(iso, "as_numpy"):
                    iso = iso.as_numpy()

                if isinstance(iso, numpy.ndarray):
                    if iso.shape == (4, 4):
                        iso = Iso3(iso)
                    else:
                        raise ValueError(f"Invalid shape for iso: expected (4, 4), got {iso.shape}")
                else:
                    raise TypeError(
                        f"Invalid type for iso: expected Iso3 or numpy.ndarray, got {type(iso)}")

            points = numpy.array([[0, 0, 0], [length, 0, 0], [0, length, 0], [0, 0, length]],
                                 dtype=numpy.float64)
            points = iso.transform_points(points)

            for axis, color in ((1, "red"), (2, "green"), (3, "blue")):
                actors.append(self.pv.add_mesh(_lines_polydata(points[[0, axis]]), color=color,
                                               **kwargs))

            if text:
                label_style = merge_style({}, font_family=font_family)
                actors.append(self.pv.add_point_labels(
                    [points[0]],
                    [text],
                    show_points=False,
                    background_color="white",
                    font_size=font_size,
                    bold=False,
                    **label_style,
                ))

        return actors

    def draw_label(self, point: Coords3, text: str, **kwargs) -> pyvista.Label:
        """
        Add a text label to the plotter.

        :param point: the position of the label in 3D space.
        :param text: the text to display in the label.
        :param kwargs: any keyword argument accepted by the `pyvista.Label` constructor.
        :return: the `pyvista.Label` that was added to the plotter.
        """
        label = pyvista.Label(text=text, position=to_tuple3(point), **kwargs)
        self.pv.add_actor(label)
        return label

    def draw_arrow(
            self,
            start: Coords3,
            direction: Coords3,
            tip_length: float = 0.25,
            tip_radius: float = 0.1,
            shaft_radius: float = 0.05,
            color: ColorLike | None = None,
            opacity: float | None = None,
            label: str | None = None,
            name: str | None = None,
            **kwargs,
    ) -> pyvista.Actor:
        """
        Add an arrow to the plotter, drawn as a solid mesh rather than a line, so that it reads as
        an arrow from any viewing angle.

        :param start: the position of the arrow's tail in 3D space.
        :param direction: the direction the arrow points. Its magnitude sets the arrow's overall
            length.
        :param tip_length: the length of the arrow's head, as a fraction of the overall length.
        :param tip_radius: the radius of the arrow's head, as a fraction of the overall length.
        :param shaft_radius: the radius of the arrow's shaft, as a fraction of the overall length.
        :param color: the color of the arrow. If None, the PyVista default applies.
        :param opacity: the opacity of the arrow, from 0.0 to 1.0. If None, the PyVista default
            applies.
        :param label: the legend label. If None, the arrow is not labeled.
        :param name: the actor name. Drawing again under the same name replaces the earlier actor.
        :param kwargs: any other keyword argument accepted by PyVista's `Plotter.add_mesh`.
        :return: the PyVista actor that was added to the plotter.
        """
        # PyVista 0.50 makes `start` and `direction` keyword-only, and warns about positional use
        # before then.
        pd = pyvista.Arrow(start=to_tuple3(start), direction=to_tuple3(direction),
                           tip_length=tip_length, tip_radius=tip_radius, shaft_radius=shaft_radius)
        kwargs = merge_style(kwargs, color=color, opacity=opacity, label=label, name=name)
        return self.pv.add_mesh(pd, **kwargs)


def _cmap_extremes(item: Any) -> Dict[str, pyvista.ColorLike]:
    # Carry a Matplotlib color map's out-of-range colors across to PyVista's own equivalents, so a
    # map like `GOM_CMAP` keeps calling out data beyond the scale instead of clamping it.
    #
    # Matplotlib has no public way to ask whether an extreme was explicitly set, so the end colors
    # are compared against the map's own first and last: a map that has set neither returns its
    # ends, and there is nothing to carry across. This deliberately avoids reading the private
    # `_rgba_over`/`_rgba_under` attributes, which are not part of Matplotlib's API.
    #
    # The result is converted to a plain tuple because `Colormap.get_over` returns a numpy array,
    # and PyVista tests these arguments for truthiness internally, which raises on an array.
    working = {}
    try:
        import numpy
        from matplotlib.colors import Colormap
    except ImportError:
        return working

    if isinstance(item, Colormap):
        for key, extreme, end in (("above_color", item.get_over(), item(1.0)),
                                  ("below_color", item.get_under(), item(0.0))):
            if not numpy.array_equal(numpy.asarray(extreme), numpy.asarray(end)):
                working[key] = tuple(float(c) for c in numpy.asarray(extreme).ravel())
    return working
