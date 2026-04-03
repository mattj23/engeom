from __future__ import annotations
from typing import List, Iterable, Tuple, Union
import numpy

from .common import LabelPlace
from engeom.geom2 import Curve2, Circle2, Aabb2, Point2, Vector2, SurfacePoint2, Arc2, Segment2
from engeom.geom3 import Vector3, Mesh, Point3, Iso3
from engeom.metrology import Distance2

PlotCoords = Union[Point2, Vector2, Iterable[float]]
PointLike = Union[Point2, Tuple[float, float], Point3, Tuple[float, float, float], Iterable[float]]

try:
    from matplotlib.pyplot import Axes, Circle
    from matplotlib.colors import ListedColormap
except ImportError:
    pass
else:

    class GomColorMap(ListedColormap):
        """
        A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software.

        You can use this to instantiate a color map, or you can use the `GOM_CMAP` object directly.
        """

        def __init__(self):
            colors = numpy.array(
                [
                    [1, 0, 160],
                    [1, 0, 255],
                    [0, 254, 255],
                    [0, 160, 0],
                    [0, 254, 0],
                    [255, 255, 0],
                    [255, 128, 0],
                    [255, 1, 0],
                ],
                dtype=numpy.float64,
            )
            colors /= 256.0
            colors = numpy.hstack((colors, numpy.ones((len(colors), 1))))
            super().__init__(colors)
            self.set_under("magenta")
            self.set_over("darkred")


    GOM_CMAP = GomColorMap()
    """
    A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software, already instantiated and
    available in the module.
    """


    def set_aspect_fill(ax: Axes):
        """
        Set the aspect ratio of a Matplotlib Axes (subplot) object to be 1:1 in x and y, while also having it expand
        to fill all available space.

        In comparison to the set_aspect('equal') method, this method will also expand the plot to prevent the overall
        figure from shrinking.  It does this by manually re-checking the x and y limits and adjusting whichever is the
        limiting value. Essentially, it will honor the larger of the two existing limits which were set before this
        function was called, and will only expand the limits on the other axis to fill the remaining space.

        Call this function after all visual elements have been added to the plot and any manual adjustments to the axis
        limits are performed. If you use fig.tight_layout(), call this function after that.
        :param ax: a Matplotlib Axes object
        :return: None
        """
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()

        bbox = ax.get_window_extent()
        width, height = bbox.width, bbox.height

        x_scale = width / (x1 - x0)
        y_scale = height / (y1 - y0)

        if y_scale > x_scale:
            y_range = y_scale / x_scale * (y1 - y0)
            y_mid = (y0 + y1) / 2
            ax.set_ylim(y_mid - y_range / 2, y_mid + y_range / 2)
        else:
            x_range = x_scale / y_scale * (x1 - x0)
            x_mid = (x0 + x1) / 2
            ax.set_xlim(x_mid - x_range / 2, x_mid + x_range / 2)


    class ViewPort:
        """
        A helper class that helps draw diagrams of 3d objects being rendered in parallel projection onto a 2d view,
        such as for generating diagrams and illustrations.
        """

        def __init__(self, view: Iso3, helper: MatplotlibAxesHelper):
            self.view = view
            self.helper = helper

        def dimension_arrow(
                self,
                point0: PointLike,
                point1: PointLike,
                label: str,
                leader_shift: PointLike | None = None,
                label_shift: PointLike | None = None,
                linewidth: float = 1.5,
                fontsize: int = 14,
                color: str = "black",
                label_position: float = 0.5,
        ):
            ldr_shift = _to_3d_point(leader_shift).coords if leader_shift is not None else Vector3.zero()
            lbl_shift = _to_3d_point(label_shift).coords if label_shift is not None else Vector3.zero()
            p0 = _to_3d_point(point0)
            p1 = _to_3d_point(point1)
            p0_end = p0 + ldr_shift
            p1_end = p1 + ldr_shift
            label_pos =  p0 + (p1 - p0) * label_position + ldr_shift + lbl_shift

            self.helper.arrow(self.view @ p0_end, self.view @ p1_end, color=color, linewidth=linewidth, arrow="<->")
            self.helper.text(label, self.view @ label_pos, fontsize=fontsize, ha="center", va="center", color=color,
                             bbox=dict(boxstyle="round", ec="black", fc="w", lw=linewidth))

            if ldr_shift.norm() > 0.0:
                self.helper.arrow(self.view @ p0, self.view @ p0_end, color=color, linewidth=linewidth * 0.75, arrow="-", linestyle="dotted")
                self.helper.arrow(self.view @ p1, self.view @ p1_end, color=color, linewidth=linewidth * 0.75, arrow="-", linestyle="dotted")


        def labeled_point(
                self,
                point: PointLike,
                label: str,
                offset_3d: PointLike | None = None,
                offset_2d: PointLike | None = None,
                fontsize: int = 14,
                color: str = "black",
                marker: str = "o",
                marker_size: float = 5.0,
                weight: str = "normal",
        ):
            """
            Generates a labeled point representation on a 3D or 2D space. The labeled point is defined
            by its position, label text, offsets for 3D and 2D rendering adjustments, and additional
            visual attributes.

            :param point: The coordinates of the point to be labeled, specified as a type compatible with PointLike.
            :param label: The label text to display near the point.
            :param offset_3d: The 3D offset value to position the label relative to the point in the 3D rendered space.
            :param offset_2d: The 2D offset value to position the label relative to the point in the 2D rendered space.
            :param fontsize: Font size of the label. Defaults to 14.
            :param color: The color of the label text. Defaults to "black".
            :param marker: The marker style used to visually indicate the point. Defaults to "o".
            :param marker_size: Size of the marker representing the point. Defaults to 5.0.
            :param weight: The font weight of the label text. Defaults to "normal".
            """
            p = _to_3d_point(point)
            o3 = Vector3.zero() if offset_3d is None else _to_3d_point(offset_3d).coords
            o2 = Vector2.zero() if offset_2d is None else Vector3(*_to_2d_point(offset_2d), 0)

            self.helper.points(self.view @ p, marker=marker, markersize=marker_size, color=color)
            self.helper.text(label, (self.view @ (p + o3)) + o2, fontsize=fontsize, ha="center", va="center",
                             color=color, weight=weight)

        def coordinate_system(
                self,
                cs: Iso3,
                length: float,
                linewidth: float = 2.0,
                label_offset: float = 0.1,
                fontsize: int = 14,
        ):
            """
            Draws a coordinate system in the 2d view using the provided view. The coordinate system is drawn as a set
            of red, green, and blue arrows representing the x, y, and z axes respectively. If one of the arrows
            aligns too closely to the view direction, it will be hidden.

            :param cs: The coordinate system to draw
            :param length: The length of the arrows representing the axes
            :param linewidth: The width of the arrows representing the axes
            :param label_offset: The additional scale factor for length to determine where the labels are placed.
            :param fontsize: The font size of the labels
            """

            def _visible(_v: Point3, _len: float) -> bool:
                return _v.with_z(0).coords.norm() > _len * 0.1

            o_view = self.view @ cs.origin
            to_draw = [("x", "red", Vector3.x_axis()),
                       ("y", "green", Vector3.y_axis()),
                       ("z", "blue", Vector3.z_axis())]

            for text, color, vector in to_draw:
                p_view: Point3 = self.view @ cs @ (Point3.origin() + vector * length)
                if _visible(p_view, length):
                    self.helper.arrow(o_view, p_view, color=color, linewidth=linewidth)
                    text_pos = self.view @ cs @ (Point3.origin() + vector * length * (1 + label_offset))
                    self.helper.text(f"${text}$", text_pos, color=color, fontsize=fontsize, ha="center", va="center")

        def mesh_outline(
                self, mesh:
                Mesh,
                visible_kwargs: dict | None = None,
                hidden_kwargs: dict | None = None,
                no_hidden=False,
                max_edge_len=10.0,
                corner_angle=None
        ):
            """
            Draws the outline of a mesh in the 2d view using the provided view and helper. The default parameters for
            the visible edges are {'color': 'black', 'linewidth': 1.0} and for the hidden edges are
            {'color': 'black', 'linewidth': 0.5, alpha=0.125}.

            :param mesh: The mesh to draw the outline of.
            :param visible_kwargs: Optional keyword arguments for the visible edges of the mesh.
            :param hidden_kwargs: Optional keyword arguments for the hidden edges of the mesh.
            :param no_hidden: If True, only the visible edges will be drawn.
            :param max_edge_len: Max length of edges for the outline, edges longer than this will be broken up
            :param corner_angle: Minimum angle between two adjacent faces for the common edge to be considered a corner and included in the outline
            """
            visible = TraceBuilder()
            hidden = TraceBuilder()

            visible_kwargs = visible_kwargs or {'color': 'black', 'linewidth': 1.0}
            hidden_kwargs = hidden_kwargs or {'color': 'black', 'linewidth': 0.5, 'alpha': 0.125}
            assert visible_kwargs is not None
            assert hidden_kwargs is not None

            points, edge_types = mesh.visual_outline(self.view.inverse() @ Vector3.z_axis(), max_edge_len, corner_angle)
            p0s = self.view.transform_points(points[:, :3])
            p1s = self.view.transform_points(points[:, 3:])
            for edge_type, (p0, p1) in zip(edge_types, zip(p0s, p1s)):
                if edge_type == 0:
                    visible.add_segment(p0, p1)
                elif not no_hidden:
                    hidden.add_segment(p0, p1)

            self.helper.ax.plot(*visible.xy, **visible_kwargs)
            if not no_hidden:
                self.helper.ax.plot(*hidden.xy, **hidden_kwargs)


    class MatplotlibAxesHelper:
        """
        A helper class for working with Matplotlib. It wraps around a Matplotlib `Axes` object and provides direct
        methods for plotting some `engeom` entities.  It also enforces the aspect ratio to be 1:1 and expands the
        subplot to fill its available space.

        !!! example
            ```python
            from matplotlib.pyplot import figure
            fig = figure()
            ax = fig.subplots()
            helper = MatplotlibAxesHelper(ax)
            ```
        """

        def __init__(self, ax: Axes, skip_aspect=False, hide_axes=False):
            """
            Initialize the helper with a Matplotlib `Axes` object.
            :param ax: The Matplotlib `Axes` object to wrap around.
            :param skip_aspect: Set this to true to skip enforcing the aspect ratio to be 1:1.
            :param hide_axes: Set this to true to hide the axes.
            """

            self.ax = ax
            if not skip_aspect:
                ax.set_aspect("equal", adjustable="datalim")

            if hide_axes:
                ax.axis("off")

        def get_3d_viewport(self, view: Iso3) -> ViewPort:
            """
            This method returns a ViewPort object that can be used to draw 3d objects in parallel projection onto the
            2d view, such as for generating diagrams and illustrations. The view should be an isometry describing a
            transformation from the 3d space into a 2d image plane where +X is to the right, +Y is up, and +Z is into
            the image plane.

            :param view: The isometry describing the transformation from 3d space into the 2d image plane.
            :return: A ViewPort object that can be used to draw 3d objects in parallel projection onto the 2d view.
            """
            return ViewPort(view, self)

        def set_bounds(self, box: Aabb2):
            """
            Set the bounds of a Matplotlib Axes object.
            :param box: an Aabb2 object
            """
            self.ax.set_xlim(box.min.x, box.max.x)
            self.ax.set_ylim(box.min.y, box.max.y)

        def arc(self, arc: Arc2, **kwargs):
            """

            :param arc:
            :return:
            """
            from matplotlib.patches import Arc

            # Arcs are drawn by matplotlib in a counter-clockwise direction, so if the sweep is negative we
            # need to swap the start and end angles.
            if arc.angle < 0.0:
                start_angle = arc.angle0 + arc.angle
                sweep_angle = -arc.angle
            else:
                start_angle = arc.angle0
                sweep_angle = arc.angle

            patch = Arc((arc.center.x, arc.center.y), 2 * arc.r, 2 * arc.r,
                        angle=numpy.rad2deg(start_angle),
                        theta1=0.0,
                        theta2=numpy.rad2deg(sweep_angle),
                        **kwargs)
            self.ax.add_patch(patch)

        def segment(self, seg: Segment2, **kwargs):
            """
            Plot a segment on a Matplotlib Axes object.
            :param seg: a Segment2 object
            :param kwargs: keyword arguments to pass to the plot function
            """
            self.ax.plot([seg.a.x, seg.b.x], [seg.a.y, seg.b.y], **kwargs)

        def plot_circle(self, *circle: Circle2 | Iterable[float], fill=False, **kwargs):
            """
            Plot a circle on a Matplotlib Axes object.
            :param circle: a Circle2 object
            :param kwargs: keyword arguments to pass to the plot function
            """
            from matplotlib.pyplot import Circle

            for cdata in circle:
                if isinstance(cdata, Circle2):
                    c = Circle((cdata.center.x, cdata.center.y), cdata.r, fill=fill, **kwargs)
                else:
                    x, y, r, *_ = cdata
                    c = Circle((x, y), r, fill=fill, **kwargs)
                self.ax.add_patch(c)

        def plot_curve(self, curve: Curve2, **kwargs):
            """
            Plot a curve on a Matplotlib Axes object.
            :param curve: a Curve2 object
            :param kwargs: keyword arguments to pass to the plot function
            """
            self.ax.plot(curve.points[:, 0], curve.points[:, 1], **kwargs)

        def fill_curve(self, curve: Curve2, **kwargs):
            """
            Fill a curve on a Matplotlib Axes object.
            :param curve: a Curve2 object (can be closed but doesn't need to be, will be closed automatically)
            :param kwargs: keyword arguments to pass to the inner Axes.fill function
            :return:
            """
            self.ax.fill(curve.points[:, 0], curve.points[:, 1], **kwargs)

        def distance(
                self,
                distance: Distance2,
                side_shift: float = 0,
                template: str = "{value:.3f}",
                fontsize: int = 10,
                label_place: LabelPlace = LabelPlace.Outside,
                label_offset: float | None = None,
                fontname: str | None = None,
                scale_value: float = 1.0,
        ):
            """
            Plot a `Distance2` object on a Matplotlib Axes, drawing the leader lines and adding a text label with the
            distance value.
            :param distance: The `Distance2` object to plot.
            :param side_shift: Shift the ends of the leader lines by this amount of data units. The direction of the
            shift is orthogonal to the distance direction, with positive values shifting to the right.
            :param template: The format string to use for the distance label. The default is "{value:.3f}".
            :param fontsize: The font size to use for the label.
            :param label_place: The placement of the label.
            :param label_offset: The distance offset to use for the label. Will have different meanings depending on
            the `label_place` parameter.
            :param fontname: The name of the font to use for the label.
            :param scale_value: A scaling factor to apply to the value before displaying it in the label. Use this to
            convert between different units of measurement without having to modify the actual value or the coordinate
            system.
            """
            pad_scale = self._font_height(12) * 1.5

            # The offset_dir is the direction from `a` to `b` projected so that it's parallel to the measurement
            # direction.
            offset_dir = distance.direction if distance.value >= 0 else -distance.direction
            center = SurfacePoint2(*distance.center.point, *offset_dir)
            center = center.shift_orthogonal(side_shift)
            leader_a = center.projection(distance.a)
            leader_b = center.projection(distance.b)

            if label_place == LabelPlace.Inside:
                label_offset = label_offset or 0.0
                label_coords = center.at_distance(label_offset)
                self.arrow(label_coords, leader_a)
                self.arrow(label_coords, leader_b)
            elif label_place == LabelPlace.Outside:
                label_offset = label_offset or pad_scale * 3
                label_coords = leader_b + offset_dir * label_offset
                self.arrow(leader_a - offset_dir * pad_scale, leader_a)
                self.arrow(label_coords, leader_b)
            elif label_place == LabelPlace.OutsideRev:
                label_offset = label_offset or pad_scale * 3
                label_coords = leader_a - offset_dir * label_offset
                self.arrow(leader_b + offset_dir * pad_scale, leader_b)
                self.arrow(label_coords, leader_a)

            # Do we need sideways leaders?
            self._line_if_needed(pad_scale, distance.a, leader_a)
            self._line_if_needed(pad_scale, distance.b, leader_b)

            kwargs = {"ha": "center", "va": "center", "fontsize": fontsize}
            if fontname is not None:
                kwargs["fontname"] = fontname

            value = distance.value * scale_value
            box_style = dict(boxstyle="round,pad=0.3", ec="black", fc="white")
            self.text(template.format(value=value), label_coords, bbox=box_style, **kwargs)

        def _line_if_needed(self, pad: float, actual: Point2, leader_end: Point2):
            half_pad = pad * 0.5
            v: Vector2 = leader_end - actual
            if v.norm() < half_pad:
                return
            work = SurfacePoint2(*actual, *v)
            t1 = work.scalar_projection(leader_end) + half_pad
            self.arrow(actual, work.at_distance(t1), arrow="-")

        def text(self, text: str, pos: PlotCoords, shift: PlotCoords | None = None, ha: str = "center",
                 va: str = "center", **kwargs):
            """
            Annotate a Matplotlib Axes object with text only, by default in the xy data plane.
            :param text: the text to annotate
            :param pos: the position of the annotation
            :param shift: an optional shift vector to apply to the position
            :param ha: horizontal alignment
            :param va: vertical alignment
            :param kwargs: keyword arguments to pass to the annotate function
            :return: the annotation object
            """
            xy = _tuplefy(pos)
            if shift is not None:
                shift = _tuplefy(shift)
                xy = (xy[0] + shift[0], xy[1] + shift[1])

            return self.ax.annotate(text, xy=xy, ha=ha, va=va, **kwargs)

        def points(self, *points: PlotCoords, marker="o", markersize=5.0, **kwargs):
            x, y = zip(*[_tuplefy(p) for p in points])
            return self.ax.plot(x, y, marker, markersize=markersize, **kwargs)

        def surface_points(self, *points: SurfacePoint2, arrow_len=1, marker="o", markersize="5", **kwargs):
            x, y = zip(*[(p.point.x, p.point.y) for p in points])
            color = kwargs.get("color", "black")
            kwargs["color"] = color
            self.ax.plot(x, y, marker, markersize=markersize, **kwargs)
            for p in points:
                p: SurfacePoint2
                self.arrow(p.point, p.at_distance(arrow_len), arrow="->", color=color)

        def labeled_arrow(self, start: PlotCoords, end: PlotCoords, text: str, fraction: float = 0.5,
                          shift: PlotCoords | None = None,
                          arrow="->", color="black", linewidth: float | None = None, linestyle="-",
                          **text_kwargs):
            """

            :param start:
            :param end:
            :param text:
            :param shift:
            :param fraction:
            :param arrow:
            :param color:
            :param linewidth:
            :param linestyle:
            :param text_kwargs: parameters to pass to the text function
            :return:
            """
            start = Point2(*_tuplefy(start))
            end = Point2(*_tuplefy(end))

            self.arrow(start, end, arrow=arrow, color=color, linewidth=linewidth, linestyle=linestyle)

            v: Vector2 = end - start
            position = start + v * fraction
            self.text(text, position, shift=shift, color=color, **text_kwargs)

        def arrow(self, start: PlotCoords, end: PlotCoords, arrow="->", color="black", linewidth: float | None = None,
                  linestyle="-", **arrow_kwargs):
            """
            Draw an arrow on a Matplotlib Axes object from `start` to `end`.
            :param start:
            :param end:
            :param arrow:
            :param color:
            :param linewidth:
            :param linestyle:
            """
            props = dict(
                arrowstyle=arrow,
                fc=color,
                ec=color,
                linewidth=linewidth,
                linestyle=linestyle,
                **arrow_kwargs
            )

            return self.ax.annotate("", xy=_tuplefy(end), xytext=_tuplefy(start), arrowprops=props)

        def _font_height(self, font_size: int) -> float:
            # Get the height of a font in data units
            fig_dpi = self.ax.figure.dpi
            font_height_inches = font_size * 1.0 / 72.0
            font_height_px = font_height_inches * fig_dpi

            px_per_data = self._get_scale()
            return font_height_px / px_per_data

        def _get_scale(self) -> float:
            # Get the scale of the plot in data units per pixel.
            x0, x1 = self.ax.get_xlim()
            y0, y1 = self.ax.get_ylim()

            bbox = self.ax.get_window_extent()
            width, height = bbox.width, bbox.height

            # Units are pixels per data unit
            x_scale = width / (x1 - x0)
            y_scale = height / (y1 - y0)

            return min(x_scale, y_scale)


def _tuplefy(item: PlotCoords) -> Tuple[float, float]:
    if isinstance(item, (Point2, Vector2)):
        return item.x, item.y
    else:
        x, y, *_ = item
        return x, y


class TraceBuilder:
    def __init__(self):
        self.xs = []
        self.ys = []
        self.c = []

    def bounds(self) -> Aabb2:
        xs = [x for x in self.xs if x is not None]
        ys = [y for y in self.ys if y is not None]
        return Aabb2(
            x_min=min(xs),
            x_max=max(xs),
            y_min=min(ys),
            y_max=max(ys),
        )

    def add_segment(self, *points: PointLike):
        self.add_points(*points)
        self.add_blank()

    def add_blank(self):
        self.xs.append(None)
        self.ys.append(None)
        self.c.append(None)

    def add_points(self, *points: PointLike):
        for x, y, *_ in points:
            self.xs.append(x)
            self.ys.append(y)

    def add_point_and_color(self, point: PointLike, color: float):
        self.xs.append(point[0])
        self.ys.append(point[1])
        self.c.append(color)

    def invert_y(self):
        self.ys = [-y if y is not None else None for y in self.ys]

    @property
    def kwargs(self):
        return dict(x=self.xs, y=self.ys)

    @property
    def xy(self):
        return self.xs, self.ys


def _to_3d_point(p: PointLike) -> Point3:
    if isinstance(p, Point3):
        return p
    elif isinstance(p, Point2):
        return Point3(p.x, p.y, 0)
    else:
        values = [float(v) for v in p]
        while len(values) < 3:
            values.append(0)
        return Point3(*values)


def _to_2d_point(p: PointLike) -> Point2:
    if isinstance(p, Point2):
        return p
    elif isinstance(p, Point3):
        return Point2(p.x, p.y)
    else:
        values = [float(v) for v in p]
        while len(values) < 2:
            values.append(0)
        return Point2(*values)
