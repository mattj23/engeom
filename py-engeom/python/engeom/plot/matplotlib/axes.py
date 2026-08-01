"""
The main helper wrapping a Matplotlib `Axes` object for drawing `engeom` entities.
"""

from __future__ import annotations

from typing import Iterable

import numpy
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, Polygon
from matplotlib.text import Annotation

from engeom.geom2 import Curve2, Circle2, Aabb2, Point2, Vector2, SurfacePoint2, Arc2, Segment2, Boundary2, \
    CubicSpline2
from engeom.geom3 import Iso3
from engeom.metrology import Distance2

from .._coerce import Coords2, to_tuple2
from .._common import LabelPlace, check_label_place
from .viewport import ViewPort3


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


class AxesHelper:
    """
    A helper class for working with Matplotlib. It wraps around a Matplotlib `Axes` object and provides direct
    methods for plotting some `engeom` entities.  It also enforces the aspect ratio to be 1:1 and expands the
    subplot to fill its available space.

    !!! example
        ```python
        from matplotlib.pyplot import figure
        fig = figure()
        ax = fig.subplots()
        helper = AxesHelper(ax)
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

    def viewport(self, view: Iso3) -> ViewPort3:
        """
        This method returns a ViewPort3 object that can be used to draw 3d objects in parallel projection onto the
        2d view, such as for generating diagrams and illustrations. The view should be an isometry describing a
        transformation from the 3d space into a 2d image plane where +X is to the right, +Y is up, and +Z is into
        the image plane.

        :param view: The isometry describing the transformation from 3d space into the 2d image plane.
        :return: A ViewPort3 object that can be used to draw 3d objects in parallel projection onto the 2d view.
        """
        return ViewPort3(view, self)

    def set_bounds(self, box: Aabb2) -> None:
        """
        Set the bounds of a Matplotlib Axes object.
        :param box: an Aabb2 object
        """
        self.ax.set_xlim(box.min.x, box.max.x)
        self.ax.set_ylim(box.min.y, box.max.y)

    def draw_arc(self, *arcs: Arc2, **kwargs) -> list[Arc]:
        """

        :param arcs:
        :return:
        """
        patches = []
        for arc in arcs:
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
            patches.append(self.ax.add_patch(patch))
        return patches

    def draw_segment(self, *segments: Segment2, **kwargs) -> list[Line2D]:
        """
        Plot a segment on a Matplotlib Axes object.
        :param segments: the Segment2 objects to draw
        :param kwargs: keyword arguments to pass to the plot function
        :return: one Line2D per segment
        """
        return [
            self.ax.plot([seg.a.x, seg.b.x], [seg.a.y, seg.b.y], **kwargs)[0]
            for seg in segments
        ]

    def draw_boundary(self, *boundaries: Boundary2, tol: float | None = None, **kwargs) -> list[Line2D]:
        lines = []
        for boundary in boundaries:
            if tol is None:
                use_tol = boundary.aabb.extent.norm() / 1000
            elif tol <= 0.0:
                raise ValueError(f"tol must be positive, got {tol}")
            else:
                use_tol = tol
            points = boundary.to_points(use_tol)
            lines.append(self.ax.plot(points[:, 0], points[:, 1], **kwargs)[0])
        return lines

    def draw_normals(self, *sources: Curve2 | Boundary2, count: int, length: float,
                     color: str | None = None, linewidth: float = 1.0, **kwargs) -> list[Annotation]:
        """
        Draw surface normals sampled at even arc-length intervals along one or more curves or
        boundaries, as arrows of a fixed length pointing outward from the geometry.

        This is mainly a debugging aid for confirming that a winding direction matches your
        intention, since the normal direction in 2D is derived from it.

        `Curve2` and `Boundary2` are both accepted, because both are arc-length parameterized and
        both report a surface point with a normal at a given length, which is all this needs.

        :param sources: the curves or boundaries to draw normals along.
        :param count: how many normals to draw along each source. Samples are evenly spaced by arc
            length from one end to the other, so on a closed source the first and last coincide.
        :param length: the length of each arrow, in data units.
        :param color: the color of the arrows. If None, the arrow default is used.
        :param linewidth: the width of the arrows.
        :param kwargs: additional keyword arguments to pass to `draw_arrow`.
        :return: one Annotation per normal drawn, across every source, in order.
        :raises TypeError: if a source is neither a `Curve2` nor a `Boundary2`.
        """
        if color is not None:
            kwargs["color"] = color
        kwargs["linewidth"] = linewidth

        arrows = []
        for source in sources:
            if not isinstance(source, (Curve2, Boundary2)):
                raise TypeError(
                    f"expected a Curve2 or Boundary2, got {type(source).__name__}"
                )
            for t in numpy.linspace(0, source.length(), count):
                station = source.at_length(t)
                arrows.append(
                    self.draw_arrow(station.point, station.surface_point.at_distance(length), **kwargs)
                )
        return arrows

    def draw_circle(self, *circles: Circle2 | Iterable[float], fill=False, **kwargs) -> list[Circle]:
        """
        Plot a circle on a Matplotlib Axes object.
        :param circles: the Circle2 objects, or (x, y, r) iterables, to draw
        :param kwargs: keyword arguments to pass to the plot function
        :return: one Circle patch per circle
        """
        patches = []
        for cdata in circles:
            if isinstance(cdata, Circle2):
                c = Circle((cdata.center.x, cdata.center.y), cdata.r, fill=fill, **kwargs)
            else:
                values = tuple(cdata)
                if len(values) != 3:
                    raise ValueError(
                        f"expected a Circle2 or an (x, y, r) iterable, got {len(values)} values"
                    )
                x, y, r = values
                c = Circle((x, y), r, fill=fill, **kwargs)
            patches.append(self.ax.add_patch(c))
        return patches

    def draw_curve(self, *curves: Curve2, **kwargs) -> list[Line2D]:
        """
        Plot a curve on a Matplotlib Axes object.
        :param curves: the Curve2 objects to draw
        :param kwargs: keyword arguments to pass to the plot function
        :return: one Line2D per curve
        """
        return [
            self.ax.plot(curve.points[:, 0], curve.points[:, 1], **kwargs)[0]
            for curve in curves
        ]

    def draw_spline(self, *splines: CubicSpline2, tol: float, **kwargs) -> list[Line2D]:
        """
        Plot a cubic Bezier spline by flattening it into an adaptive polyline whose linear
        segments deviate from the underlying curve by no more than `tol`.

        :param splines: the CubicSpline2 objects to draw
        :param tol: maximum Euclidean deviation between the drawn polyline and the spline.
            Smaller values produce a smoother-looking curve at the cost of more line segments.
        :param kwargs: keyword arguments to pass to `Axes.plot`
        :return: one Line2D per spline
        """
        lines = []
        for spline in splines:
            points = spline.polyline(tol)
            lines.append(self.ax.plot(points[:, 0], points[:, 1], **kwargs)[0])
        return lines

    def fill_curve(self, *curves: Curve2, **kwargs) -> list[Polygon]:
        """
        Fill a curve on a Matplotlib Axes object.
        :param curves: the Curve2 objects to fill (each can be closed but doesn't need to be, and
            will be closed automatically)
        :param kwargs: keyword arguments to pass to the inner Axes.fill function
        :return: one Polygon per curve
        """
        return [
            self.ax.fill(curve.points[:, 0], curve.points[:, 1], **kwargs)[0]
            for curve in curves
        ]

    def draw_distance(
            self,
            distance: Distance2,
            side_shift: float = 0,
            template: str = "{value:.3f}",
            fontsize: int = 10,
            label_place: LabelPlace = "outside",
            label_offset: float | None = None,
            fontname: str | None = None,
            scale_value: float = 1.0,
    ) -> list[Artist]:
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
        :raises ValueError: if `label_place` is not one of the valid placement tokens.
        """
        label_place = check_label_place(label_place)
        artists: list[Artist] = []
        pad_scale = self._font_height(12) * 1.5

        # The offset_dir is the direction from `a` to `b` projected so that it's parallel to the measurement
        # direction.
        offset_dir = distance.direction if distance.value >= 0 else -distance.direction
        center = SurfacePoint2(*distance.center.point, *offset_dir)
        center = center.shifted_orthogonal(side_shift)
        leader_a = center.projection(distance.a)
        leader_b = center.projection(distance.b)

        if label_place == "inside":
            label_offset = label_offset if label_offset is not None else 0.0
            label_coords = center.at_distance(label_offset)
            artists.append(self.draw_arrow(label_coords, leader_a))
            artists.append(self.draw_arrow(label_coords, leader_b))
        elif label_place == "outside":
            label_offset = label_offset if label_offset is not None else pad_scale * 3
            label_coords = leader_b + offset_dir * label_offset
            artists.append(self.draw_arrow(leader_a - offset_dir * pad_scale, leader_a))
            artists.append(self.draw_arrow(label_coords, leader_b))
        else:  # "outside_rev", the only remaining token after validation
            label_offset = label_offset if label_offset is not None else pad_scale * 3
            label_coords = leader_a - offset_dir * label_offset
            artists.append(self.draw_arrow(leader_b + offset_dir * pad_scale, leader_b))
            artists.append(self.draw_arrow(label_coords, leader_a))

        # Do we need sideways leaders?
        artists.extend(self._line_if_needed(pad_scale, distance.a, leader_a))
        artists.extend(self._line_if_needed(pad_scale, distance.b, leader_b))

        kwargs = {"ha": "center", "va": "center", "fontsize": fontsize}
        if fontname is not None:
            kwargs["fontname"] = fontname

        value = distance.value * scale_value
        box_style = dict(boxstyle="round,pad=0.3", ec="black", fc="white")
        artists.append(self.draw_text(template.format(value=value), label_coords, bbox=box_style, **kwargs))
        return artists

    def _line_if_needed(self, pad: float, actual: Point2, leader_end: Point2) -> list[Annotation]:
        half_pad = pad * 0.5
        v: Vector2 = leader_end - actual
        if v.norm() < half_pad:
            return []
        work = SurfacePoint2(*actual, *v)
        t1 = work.scalar_projection(leader_end) + half_pad
        return [self.draw_arrow(actual, work.at_distance(t1), arrow="-")]

    def draw_text(self, text: str, pos: Coords2, shift: Coords2 | None = None, ha: str = "center",
                  va: str = "center", **kwargs) -> Annotation:
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
        xy = to_tuple2(pos)
        if shift is not None:
            shift = to_tuple2(shift)
            xy = (xy[0] + shift[0], xy[1] + shift[1])

        return self.ax.annotate(text, xy=xy, ha=ha, va=va, **kwargs)

    def draw_point(self, *points: Coords2, marker="o", markersize=5.0, **kwargs) -> Line2D:
        # An empty call still produces a real (empty) artist, rather than failing in `zip`, so that
        # drawing a computed and possibly empty set of points needs no special case at the call site.
        coords = [to_tuple2(p) for p in points]
        x, y = zip(*coords) if coords else ((), ())
        return self.ax.plot(x, y, marker, markersize=markersize, **kwargs)[0]

    def draw_surface_point(self, *points: SurfacePoint2, arrow_len=1, marker="o", markersize=5.0,
                           **kwargs) -> list[Artist]:
        coords = [(p.point.x, p.point.y) for p in points]
        x, y = zip(*coords) if coords else ((), ())
        color = kwargs.get("color", "black")
        kwargs["color"] = color
        artists: list[Artist] = [self.ax.plot(x, y, marker, markersize=markersize, **kwargs)[0]]
        for p in points:
            p: SurfacePoint2
            artists.append(self.draw_arrow(p.point, p.at_distance(arrow_len), arrow="->", color=color))
        return artists

    def draw_labeled_arrow(self, start: Coords2, end: Coords2, text: str, fraction: float = 0.5,
                           shift: Coords2 | None = None,
                           arrow="->", color="black", linewidth: float | None = None, linestyle="-",
                           **text_kwargs) -> list[Artist]:
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
        start = Point2(*to_tuple2(start))
        end = Point2(*to_tuple2(end))

        drawn = self.draw_arrow(start, end, arrow=arrow, color=color, linewidth=linewidth, linestyle=linestyle)

        v: Vector2 = end - start
        position = start + v * fraction
        return [drawn, self.draw_text(text, position, shift=shift, color=color, **text_kwargs)]

    def draw_arrow(self, start: Coords2, end: Coords2, arrow="->", color="black", linewidth: float | None = None,
                   linestyle="-", **arrow_kwargs) -> Annotation:
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

        return self.ax.annotate("", xy=to_tuple2(end), xytext=to_tuple2(start), arrowprops=props)

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
