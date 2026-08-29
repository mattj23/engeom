"""
Helpers for drawing `engeom` entities into a PyVista `Plotter` object.

Requires `pyvista`, which is an optional dependency: `pip install engeom[pyvista]`.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy

try:
    import pyvista
    from pyvista import ColorLike
except ImportError as _e:
    raise ImportError(
        "engeom.plot.pyvista requires the `pyvista` package, which is not installed. "
        "Install it with `pip install engeom[pyvista]`."
    ) from _e

from engeom.geom3 import Mesh3, Curve3, Iso3, SurfacePoint3, PointCloud3, Circle3, Sphere3
from engeom.metrology import Distance3

from ._coerce import Coords3, to_tuple3
from ._common import LabelPlace, check_label_place, plane_basis

__all__ = ["PlotterHelper"]

# When a draw method's chordal tolerance is left as `None`, it defaults to this fraction of the
# entity's radius. Perceived smoothness in a plot scales with the entity's size on screen, so this
# radius-relative default keeps entities equally smooth without requiring the caller to calculate
# it.
_DEFAULT_REL_TOL = 1.0e-3


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


class PlotterHelper:
    """
    A helper class for working with PyVista. It wraps around a PyVista `Plotter` object and provides direct methods
    for plotting some `engeom` entities.

    Every method that adds an actor to the plotter is prefixed `draw_`, matching the Matplotlib
    helper, so `helper.draw_<tab>` lists everything that can be drawn.

    !!! example
        ```python
        import pyvista
        plotter = pyvista.Plotter()
        helper = PlotterHelper(plotter)
        ```
    """

    def __init__(self, plotter: pyvista.Plotter):
        """
        Initialize the helper with a PyVista `Plotter` object.

        :param plotter: The PyVista `Plotter` object to wrap around.
        """
        self.pv = plotter

    @staticmethod
    def with_new_plotter(**kwargs) -> PlotterHelper:
        """
        Create a new PyVista plotter and return a helper for it.

        :param kwargs: Additional keyword arguments to pass to the PyVista `Plotter` constructor.
        :return: A new `PlotterHelper` instance.
        """
        plotter = pyvista.Plotter(**kwargs)
        return PlotterHelper(plotter)

    def show(self):
        """
        Show the plotter.
        """
        self.pv.show()

    def draw_sphere(
            self,
            sphere: Sphere3,
            tol: float | None = None,
            color: ColorLike | None = "green",
            opacity: float | None = None,
    ) -> pyvista.vtkActor:
        """
        Add a `Sphere3` to the plotter

        :param sphere: The sphere to plot.
        :param tol: The maximum deviation of the tessellation from the true sphere. Smaller values produce a smoother
        result. If `None`, a thousandth of the sphere's radius is used.
        :param color: The color of the sphere. Set to `None` to use the default PyVista color.
        :param opacity: Opacity of the sphere (0.0–1.0). If `None`, the default PyVista opacity is used.
        :return: The PyVista `vtkActor` representing the sphere.
        """
        if tol is None:
            tol = sphere.r * _DEFAULT_REL_TOL

        mesh = Mesh3.create_sphere(sphere.r, tol)
        mesh.transform_in_place(Iso3.from_translation(*sphere.center))
        return self.draw_mesh(mesh, color=color, opacity=opacity)

    def draw_circle(
            self,
            circle: Circle3,
            tol: float | None = None,
            edge_color: ColorLike | None = "black",
            face_color: ColorLike | None = "green",
            edge_width: float = 4.0,
            line_opacity: float | None = None,
            face_opacity: float | None = None,
    ) -> list[pyvista.vtkActor]:
        """
        Add a `Circle3` to the plotter, optionally rendering a filled face, an outline edge, or both.

        :param circle: The circle to plot.
        :param tol: The maximum chordal deviation of the segments approximating the circle. If `None`,
            one-thousandth of the circle's radius is used.
        :param edge_color: The color of the outline edge. Set to `None` to suppress the edge entirely.
        :param face_color: The color of the filled face. Set to `None` to suppress the face entirely.
        :param edge_width: The pixel width of the outline edge.
        :param line_opacity: Opacity of the outline edge (0.0–1.0). If `None`, the default PyVista opacity is used.
        :param face_opacity: Opacity of the filled face (0.0–1.0). If `None`, the default PyVista opacity is used.
        :return: A list of the PyVista actors that were added to the plotter, which will be empty if both the edge
        and the face were suppressed.
        """
        # A `Circle3` carries only a center, a normal, and a radius, so an in-plane orientation has to
        # be chosen before it can be turned into vertices. Which one does not matter, as the result is
        # rotationally symmetric about the normal either way.
        if tol is None:
            tol = circle.r * _DEFAULT_REL_TOL

        u, v = plane_basis(circle.normal)
        center = numpy.array(to_tuple3(circle.center), dtype=numpy.float64)
        u_r = numpy.array(to_tuple3(u), dtype=numpy.float64) * circle.r
        v_r = numpy.array(to_tuple3(v), dtype=numpy.float64) * circle.r

        actors = []
        if edge_color is not None:
            n = _segments_for_tol(circle.r, tol)
            theta = numpy.linspace(0.0, 2.0 * numpy.pi, n + 1)
            points = center + numpy.outer(numpy.cos(theta), u_r) + numpy.outer(numpy.sin(theta), v_r)

            kwargs = {"color": edge_color, "width": edge_width}
            if line_opacity is not None:
                kwargs["opacity"] = line_opacity
            actors.append(self.pv.add_lines(points, connected=True, **kwargs))

        if face_color is not None:
            # `create_circle` builds the circle in the XY plane at the origin, so it needs the same
            # basis assembled into a transform to put it where the entity actually is.
            matrix = numpy.eye(4)
            matrix[:3, 0] = to_tuple3(u)
            matrix[:3, 1] = to_tuple3(v)
            matrix[:3, 2] = to_tuple3(circle.normal)
            matrix[:3, 3] = center
            mesh = Mesh3.create_circle(circle.r, tol)
            mesh.transform_in_place(Iso3(matrix))

            kwargs = {"color": face_color}
            if face_opacity is not None:
                kwargs["opacity"] = face_opacity
            actors.append(self.draw_mesh(mesh, **kwargs))
        return actors

    def draw_point(self, *points, color: pyvista.ColorLike = "b", point_size: float = 5.0,
                   render_points_as_spheres: bool = True, **kwargs) -> pyvista.vtkActor:
        """
        Add one or more discrete points to be plotted.
        :param points: The points to add.
        :param color: The color to use for the point(s).
        :param point_size: The size of the point(s).
        :param render_points_as_spheres: Whether to render the points as spheres or not.
        :param kwargs: Additional keyword arguments to pass to the PyVista `Plotter.add_points` method.
        :return: The PyVista actor that was added to the plotter.
        """
        return self.pv.add_points(
            numpy.array([to_tuple3(p) for p in points], dtype=numpy.float64),
            color=color,
            point_size=point_size,
            render_points_as_spheres=render_points_as_spheres,
            **kwargs
        )

    def draw_curve(
            self,
            *curves: Curve3,
            color: pyvista.ColorLike = "w",
            width: float = 3.0,
            label: str | None = None,
            name: str | None = None,
    ) -> pyvista.vtkActor:
        """
        Add one or more curves to be plotted.
        :param curves: The curves to add.
        :param color: The color to use for the curve(s).
        :param width: The line width to use for the curve(s).
        :param label: The label to use for the curve(s) if a legend is present.
        :param name: The name to use for the actor in the scene.
        :return: The PyVista actor that was added to the plotter.
        """
        curve_vertices = []
        for curve in curves:
            b = curve.points[1:-1]
            c = numpy.zeros((len(curve.points) + len(b), 3), dtype=curve.points.dtype)
            c[0::2, :] = curve.points[0:-1]
            c[1:-1:2, :] = b
            c[-1] = curve.points[-1]
            curve_vertices.append(c)

        vertices = numpy.concatenate(curve_vertices, axis=0)
        return self.pv.add_lines(
            vertices,
            color=color,
            width=width,
            label=label,
            name=name,
        )

    def draw_mesh(self, mesh: Mesh3, **kwargs) -> pyvista.vtkActor:
        """
        Add an `engeom` mesh to be plotted. Additional keyword arguments will be passed directly to the PyVista
        `Plotter.add_mesh` method, allowing for customization of the mesh appearance.

        :param mesh: The mesh object to add to the plotter scene
        :return: The PyVista actor that was added to the plotter.
        """
        if "cmap" in kwargs:
            cmap_extremes = _cmap_extremes(kwargs["cmap"])
            kwargs.update(cmap_extremes)

        prefix = numpy.ones((mesh.faces.shape[0], 1), dtype=mesh.faces.dtype)
        faces = numpy.hstack((prefix * 3, mesh.faces))
        data = pyvista.PolyData(mesh.points, faces)
        return self.pv.add_mesh(data, **kwargs)

    def draw_point_cloud(self, cloud: PointCloud3, use_colors: bool = True,
                         normal_arrow_size: float = 0.0, **kwargs) -> list[pyvista.vtkActor]:
        """
        Add a `PointCloud3` to the plotter, optionally colored by the cloud's own per-point colors and
        optionally overlaid with arrows showing the per-point normals.

        :param cloud: The point cloud to plot.
        :param use_colors: If True and the cloud carries per-point colors, they are used as RGBA scalars
        and any `color` keyword argument is discarded.
        :param normal_arrow_size: The length of the normal arrows. Leave at the default of 0.0 to suppress
        the arrows entirely; they are only drawn for a positive size, and only if the cloud has normals.
        :param kwargs: Additional keyword arguments to pass to the PyVista `Plotter.add_points` method.
        :return: A list of the PyVista actors that were added to the plotter.
        """
        actors = []
        if normal_arrow_size > 0.0 and cloud.point_normals is not None:
            arrow_color = kwargs.get("color", "gray")
            arrow_actor = self.pv.add_arrows(cloud.points, cloud.point_normals, mag=normal_arrow_size,
                                             color=arrow_color, reset_camera=False)
            actors.append(arrow_actor)

        if use_colors and cloud.point_colors is not None:
            kwargs.pop("color", None)  # Remove color if it's set, as colors will be used from the cloud
            kwargs["scalars"] = cloud.point_colors
            kwargs["rgba"] = True

        point_actor = self.pv.add_points(cloud.points, **kwargs)
        actors.append(point_actor)

        return actors

    def draw_distance(
            self,
            distance: Distance3,
            template: str = "{value:.3f}",
            label_place: LabelPlace = "outside",
            label_offset: float | None = None,
            font_size: int = 12,
            scale_value: float = 1.0,
            font_family=None,
    ) -> list[pyvista.vtkActor]:
        """
        Add a distance entity to the plotter.
        :param distance: The distance entity to add.
        :param template: A format string to use for the label. The `value` key will be replaced with the actual
        value read from the measurement.
        :param label_place: The placement of the label relative to the distance entity's anchor points.
        :param label_offset: The distance offset to use for the label. Will have different meanings depending on
        the `label_place` parameter.
        :param font_size: The size of the text to use for the label.
        :param scale_value: A scaling factor to apply to the value before displaying it in the label. Use this to
        convert between different units of measurement without having to modify the actual value or the coordinate
        system.
        :param font_family: The font family to use for the label.
        :return: A list of the PyVista actors that were added to the plotter: the anchor spheres, the leader lines,
        and the label.
        :raises ValueError: if `label_place` is not one of the valid placement tokens.
        """
        label_place = check_label_place(label_place)
        label_offset = label_offset if label_offset is not None else max(abs(distance.value), 1.0) * 3

        # The offset_dir is the direction from `a` to `b` projected so that it's parallel to the measurement
        # direction.
        offset_dir = distance.direction if distance.value >= 0 else -distance.direction

        # Rather than arrows, we'll use spheres to indicate the anchor points at the end of the leader lines
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
        actors = [self.pv.add_points(points, color="black", point_size=4,
                                     render_points_as_spheres=True)]

        lines = builder.build()
        actors.append(self.pv.add_lines(lines, color="black", width=1.5))

        value = distance.value * scale_value
        actors.append(self.pv.add_point_labels(
            [to_tuple3(label_coords)],
            [template.format(value=value)],
            show_points=False,
            background_color="white",
            font_family=font_family,
            # justification_vertical="center",
            # justification_horizontal="center",
            font_size=font_size,
            bold=False,
        ))
        return actors

    def draw_coordinate_system(self, iso, size: float = 1.0, line_width=3.0, label: str | None = None,
                               label_size: int = 12) -> list[pyvista.vtkActor]:
        """
        Add a coordinate frame to the plotter.  This will appear as three lines, with X in red, Y in green,
        and Z in blue.  The length of each line is determined by the `size` parameter.
        :param iso: The isometry to use as the origin and orientation of the coordinate frame. May be an `Iso3`, a
        4x4 `numpy.ndarray` that validly converts into an `Iso3`, or anything with an `as_numpy` method that
        returns a valid 4x4 `numpy.ndarray`.
        :param size: The length of each line in the coordinate frame.
        :param line_width: The width of the lines in the coordinate frame.
        :param label: An optional label to display at the origin of the coordinate frame.
        :param label_size: The size of the label text.
        :return: A list of the PyVista actors that were added to the plotter: one per axis, plus the label if one
        was requested.
        :raises ValueError: if `iso` is a `numpy.ndarray` whose shape is not (4, 4).
        :raises TypeError: if `iso` is neither an `Iso3` nor something convertible to a 4x4 `numpy.ndarray`.
        """
        if not isinstance(iso, Iso3):
            if hasattr(iso, "as_numpy"):
                iso = iso.as_numpy()

            if isinstance(iso, numpy.ndarray):
                if iso.shape == (4, 4):
                    iso = Iso3(iso)
                else:
                    raise ValueError(f"Invalid shape for iso: expected (4, 4), got {iso.shape}")
            else:
                raise TypeError(f"Invalid type for iso: expected Iso3 or numpy.ndarray, got {type(iso)}")

        points = numpy.array([[0, 0, 0], [size, 0, 0], [0, size, 0], [0, 0, size]], dtype=numpy.float64)
        points = iso.transform_points(points)

        actors = [self.pv.add_lines(points[[0, 1]], color="red", width=line_width),
                  self.pv.add_lines(points[[0, 2]], color="green", width=line_width),
                  self.pv.add_lines(points[[0, 3]], color="blue", width=line_width)]

        if label:
            actors.append(self.pv.add_point_labels(
                [points[0]],
                [label],
                show_points=False,
                background_color="white",
                font_family="courier",
                font_size=label_size,
                bold=False,
            ))

        return actors

    def draw_label(self, point: Coords3, text: str, **kwargs) -> pyvista.Label:
        """
        Add a text label to the plotter.
        :param point: The position of the label in 3D space.
        :param text: The text to display in the label.
        :param kwargs: Additional keyword arguments to pass to the `pyvista.Label` constructor.
        :return: The `pyvista.Label` that was added to the plotter.
        """
        label = pyvista.Label(text=text, position=to_tuple3(point), **kwargs)
        self.pv.add_actor(label)
        return label

    def draw_arrow(self, start: Coords3, direction: Coords3,
                   tip_length: float = 0.25,
                   tip_radius: float = 0.1,
                   shaft_radius: float = 0.05,
                   **kwargs) -> pyvista.vtkActor:
        """
        Add an arrow to the plotter, drawn as a solid mesh rather than a line, so that it reads as an arrow from
        any viewing angle.

        :param start: The position of the arrow's tail in 3D space.
        :param direction: The direction the arrow points. Its magnitude sets the arrow's overall length.
        :param tip_length: The length of the arrow's head, as a fraction of the overall length.
        :param tip_radius: The radius of the arrow's head, as a fraction of the overall length.
        :param shaft_radius: The radius of the arrow's shaft, as a fraction of the overall length.
        :param kwargs: Additional keyword arguments to pass to the PyVista `Plotter.add_mesh` method.
        :return: The PyVista actor that was added to the plotter.
        """
        # PyVista 0.50 makes `start` and `direction` keyword-only, and warns about positional use
        # before then.
        pd = pyvista.Arrow(start=to_tuple3(start), direction=to_tuple3(direction),
                           tip_length=tip_length, tip_radius=tip_radius, shaft_radius=shaft_radius)
        kwargs.setdefault("color", "black")
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


class LineBuilder:
    def __init__(self):
        self.vertices = []
        self._skip = 1

    def add(self, points: Coords3):
        if self.vertices:
            if self._skip > 0:
                self._skip -= 1
            else:
                self.vertices.append(self.vertices[-1])

        self.vertices.append(to_tuple3(points))

    def skip(self):
        self._skip = 2

    def build(self) -> numpy.ndarray:
        return numpy.array(self.vertices, dtype=numpy.float64)
