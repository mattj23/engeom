"""
This module contains helper functions for working with PyVista.
"""

from __future__ import annotations

from typing import List, Any, Dict, Union, Iterable, Tuple

import numpy

from engeom.geom3 import Mesh, Curve3, Vector3, Point3, Iso3, SurfacePoint3, PointCloud
from engeom.metrology import Distance3
from .common import LabelPlace

PlotCoords = Union[Point3, Vector3, Iterable[float]]

try:
    import pyvista
except ImportError:
    pass
else:
    class PyvistaPlotterHelper:
        """
        A helper class for working with PyVista. It wraps around a PyVista `Plotter` object and provides direct methods
        for plotting some `engeom` entities.

        !!! example
            ```python
            import pyvista
            plotter = pyvista.Plotter()
            helper = PyvistaPlotterHelper(plotter)
            ```
        """

        def __init__(self, plotter: pyvista.Plotter):
            """
            Initialize the helper with a PyVista `Plotter` object.

            :param plotter: The PyVista `Plotter` object to wrap around.
            """
            self.plotter = plotter

        def add_points(self, *points, color: pyvista.ColorLike = "b", point_size: float = 5.0,
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
            return self.plotter.add_points(
                numpy.array([_tuplefy(p) for p in points], dtype=numpy.float64),
                color=color,
                point_size=point_size,
                render_points_as_spheres=render_points_as_spheres,
                **kwargs
            )

        def add_curves(
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
            return self.plotter.add_lines(
                vertices,
                color=color,
                width=width,
                label=label,
                name=name,
            )

        def add_mesh(self, mesh: Mesh, **kwargs) -> pyvista.vtkActor:
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
            data = pyvista.PolyData(mesh.vertices, faces)
            return self.plotter.add_mesh(data, **kwargs)

        def add_point_cloud(self, cloud: PointCloud, use_colors: bool = True, normal_arrow_size: float = 0.0, **kwargs):
            actors = []
            if normal_arrow_size >= 0.0 and cloud.normals is not None:
                arrow_color = kwargs.get("color", "gray")
                arrow_actor = self.plotter.add_arrows(cloud.points, cloud.normals, mag=normal_arrow_size,
                                                      color=arrow_color, reset_camera=False)
                actors.append(arrow_actor)

            if use_colors and cloud.colors is not None:
                kwargs.pop("color", None)  # Remove color if it's set, as colors will be used from the cloud
                kwargs["scalars"] = cloud.colors
                kwargs["rgba"] = True

            point_actor = self.plotter.add_points(cloud.points, **kwargs)
            actors.append(point_actor)

            return actors

        def distance(
                self,
                distance: Distance3,
                template: str = "{value:.3f}",
                label_place: LabelPlace = LabelPlace.Outside,
                label_offset: float | None = None,
                font_size: int = 12,
                scale_value: float = 1.0,
                font_family=None,
        ):
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
            """
            label_offset = label_offset or max(abs(distance.value), 1.0) * 3

            # The offset_dir is the direction from `a` to `b` projected so that it's parallel to the measurement
            # direction.
            offset_dir = distance.direction if distance.value >= 0 else -distance.direction

            # Rather than arrows, we'll use spheres to indicate the anchor points at the end of the leader lines
            spheres = [distance.a, distance.b]
            builder = LineBuilder()

            if label_place == LabelPlace.Inside:
                c = SurfacePoint3(*distance.center.point, *offset_dir)
                label_coords = c.at_distance(label_offset)

                builder.add(distance.a)
                builder.add(distance.a - offset_dir * label_offset * 0.25)
                builder.skip()

                builder.add(distance.b)
                builder.add(distance.b + offset_dir * label_offset * 0.25)

            elif label_place == LabelPlace.Outside:
                label_coords = distance.b + offset_dir * label_offset

                builder.add(distance.a)
                builder.add(distance.a - offset_dir * label_offset * 0.25)
                builder.skip()

                builder.add(distance.b)
                builder.add(label_coords)

            elif label_place == LabelPlace.OutsideRev:
                label_coords = distance.a - offset_dir * label_offset

                builder.add(distance.b)
                builder.add(distance.b + offset_dir * label_offset * 0.25)
                builder.skip()

                builder.add(distance.a)
                builder.add(label_coords)

            points = numpy.array([_tuplefy(p) for p in spheres], dtype=numpy.float64)
            self.plotter.add_points(points, color="black", point_size=4, render_points_as_spheres=True)

            lines = builder.build()
            self.plotter.add_lines(lines, color="black", width=1.5)

            value = distance.value * scale_value
            self.plotter.add_point_labels(
                [_tuplefy(label_coords)],
                [template.format(value=value)],
                show_points=False,
                background_color="white",
                font_family=font_family,
                # justification_vertical="center",
                # justification_horizontal="center",
                font_size=font_size,
                bold=False,
            )

        def coordinate_frame(self, iso, size: float = 1.0, line_width=3.0, label: str | None = None,
                             label_size: int = 12):
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
            """
            if not isinstance(iso, Iso3):
                if hasattr(iso, "as_numpy"):
                    iso = iso.as_numpy()

                if isinstance(iso, numpy.ndarray):
                    if iso.shape == (4, 4):
                        iso = Iso3(iso)
                    else:
                        raise ValueError("Invalid shape for iso: expected (4, 4), got {iso.shape}")
                else:
                    raise TypeError("Invalid type for iso: expected Iso3 or numpy.ndarray, got {type(iso)}")

            points = numpy.array([[0, 0, 0], [size, 0, 0], [0, size, 0], [0, 0, size]], dtype=numpy.float64)
            points = iso.transform_points(points)

            actors = [self.plotter.add_lines(points[[0, 1]], color="red", width=line_width),
                      self.plotter.add_lines(points[[0, 2]], color="green", width=line_width),
                      self.plotter.add_lines(points[[0, 3]], color="blue", width=line_width)]

            if label:
                actors.append(self.plotter.add_point_labels(
                    [points[0]],
                    [label],
                    show_points=False,
                    background_color="white",
                    font_family="courier",
                    font_size=label_size,
                    bold=False,
                ))

            return actors

        def label(self, point: PlotCoords, text: str, **kwargs):
            """
            Add a text label to the plotter.
            :param point: The position of the label in 3D space.
            :param text: The text to display in the label.
            :param kwargs: Additional keyword arguments to pass to the `pyvista.Label` constructor.
            """
            label = pyvista.Label(text=text, position=_tuplefy(point), **kwargs)
            self.plotter.add_actor(label)

        def arrow(self, start: PlotCoords, direction: PlotCoords,
                  tip_length: float = 0.25,
                  tip_radius: float = 0.1,
                  shaft_radius: float = 0.05,
                  **kwargs):
            pd = pyvista.Arrow(_tuplefy(start), _tuplefy(direction), tip_length=tip_length, tip_radius=tip_radius,
                               shaft_radius=shaft_radius)
            self.plotter.add_mesh(pd, **kwargs, color="black")


    def _cmap_extremes(item: Any) -> Dict[str, pyvista.ColorLike]:
        working = {}
        try:
            from matplotlib.colors import Colormap
        except ImportError:
            return working
        else:
            if isinstance(item, Colormap):
                over = getattr(item, "_rgba_over", None)
                under = getattr(item, "_rgba_under", None)
                if over is not None:
                    working["above_color"] = over
                if under is not None:
                    working["below_color"] = under
            return working


class LineBuilder:
    def __init__(self):
        self.vertices = []
        self._skip = 1

    def add(self, points: PlotCoords):
        if self.vertices:
            if self._skip > 0:
                self._skip -= 1
            else:
                self.vertices.append(self.vertices[-1])

        self.vertices.append(_tuplefy(points))

    def skip(self):
        self._skip = 2

    def build(self) -> numpy.ndarray:
        return numpy.array(self.vertices, dtype=numpy.float64)


def _tuplefy(item: PlotCoords) -> Tuple[float, float, float]:
    if isinstance(item, (Point3, Vector3)):
        return item.x, item.y, item.z
    else:
        x, y, z, *_ = item
        return x, y, z
