"""
The main helper wrapping a Matplotlib `Axes` object for drawing `engeom` entities.
"""

from __future__ import annotations

from typing import Iterable

import numpy
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, Polygon, Rectangle
from matplotlib.text import Annotation

from engeom.geom2 import Curve2, Circle2, Aabb2, Line2, Point2, Vector2, SurfacePoint2, Arc2, Segment2, Boundary2, \
    CubicSpline2
from engeom.geom3 import Iso3
from engeom.metrology import Distance2

from .._coerce import Coords2, to_tuple2
from .._common import LabelPlace, check_label_place
from ._style import merge_style
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
    Wraps a Matplotlib ``Axes`` and draws `engeom` entities onto it.

    Every method that adds an artist is prefixed ``draw_``, so ``helper.draw_<tab>`` enumerates
    everything that can be drawn; a method without that prefix configures or queries instead. Each
    ``draw_`` method returns the Matplotlib artist or artists it created, so anything the helper
    does not expose can still be reached by configuring the returned object directly.

    The wrapped ``Axes`` stays available as `ax`, which is the escape hatch for the rest of the
    Matplotlib API: ``helper.ax.legend()``, ``helper.ax.set_title(...)``, and so on.

    By default the aspect ratio is set to 1:1, since geometry drawn at unequal scales is misleading.

    !!! example
        ```python
        from matplotlib.pyplot import figure
        fig = figure()
        ax = fig.subplots()
        helper = AxesHelper(ax)
        ```
    """

    def __init__(self, ax: Axes, skip_aspect: bool = False, hide_axes: bool = False):
        """
        Wrap a Matplotlib ``Axes`` object.

        :param ax: the ``Axes`` object to wrap. It is not copied, so anything drawn through the
            helper appears on the caller's own figure.
        :param skip_aspect: set to True to leave the aspect ratio alone. By default the aspect is
            forced to 1:1 with ``adjustable="datalim"``, so that geometry is not visually distorted.
        :param hide_axes: set to True to hide the axis lines, ticks, and labels, which is usually
            what you want when the figure is a diagram rather than a plot.
        """

        self.ax = ax
        if not skip_aspect:
            ax.set_aspect("equal", adjustable="datalim")

        if hide_axes:
            ax.axis("off")

    def viewport(self, view: Iso3) -> ViewPort3:
        """
        Create a viewport for drawing 3D entities onto this 2D axes in parallel projection, such as
        for generating diagrams and illustrations of 3D geometry.

        The view isometry describes the transformation from 3D space into the 2D image plane, where
        +X is to the right, +Y is up, and +Z points into the image plane.

        :param view: the isometry transforming 3D space into the 2D image plane.
        :return: a `ViewPort3` bound to this helper, which draws onto the same axes.
        """
        return ViewPort3(view, self)

    def set_bounds(self, box: Aabb2) -> None:
        """
        Set both axis limits from a bounding box, in a single call.

        :param box: the box to fit the axes to. Its extents are applied directly, with no padding.
        """
        self.ax.set_xlim(box.min.x, box.max.x)
        self.ax.set_ylim(box.min.y, box.max.y)

    def draw_arc(self, *arcs: Arc2, color: str | None = None, linewidth: float | None = None,
                 linestyle: str | None = None, alpha: float | None = None, label: str | None = None,
                 **kwargs) -> list[Arc]:
        """
        Draw one or more circular arcs. Matplotlib's ``Arc`` is a stroke-only patch and cannot be
        filled.

        `Arc2` carries a signed sweep, while Matplotlib's ``Arc`` patch always sweeps
        counter-clockwise from its start angle. A negative sweep is therefore converted here into
        the equivalent positive sweep from the other end, so that clockwise arcs render correctly.

        :param arcs: the arcs to draw.
        :param color: the stroke color. If None, the Matplotlib default applies.
        :param linewidth: the stroke width in points. If None, the Matplotlib default applies.
        :param linestyle: the stroke style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the arc is fully opaque.
        :param label: the legend label. If None, the arc is not labeled.
        :param kwargs: any other keyword argument accepted by the Matplotlib ``Arc`` patch.
        :return: one ``Arc`` patch per arc, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
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

    def draw_segment(self, *segments: Segment2, color: str | None = None,
                     linewidth: float | None = None, linestyle: str | None = None,
                     alpha: float | None = None, label: str | None = None,
                     **kwargs) -> list[Line2D]:
        """
        Draw one or more line segments, each as its own separate line.

        :param segments: the segments to draw.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per segment, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
        return [
            self.ax.plot([seg.a.x, seg.b.x], [seg.a.y, seg.b.y], **kwargs)[0]
            for seg in segments
        ]

    def draw_line(self, *lines: Line2, t: tuple[float, float] | None = None,
                  color: str | None = None, linewidth: float | None = None,
                  linestyle: str | None = None, alpha: float | None = None,
                  label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more `Line2` entities, either as infinite lines spanning the axes or as finite
        pieces cut from a parameter range.

        A `Line2` has no endpoints, so by default each one is drawn as a Matplotlib ``AxLine``,
        which re-clips itself to the axes whenever the view changes. Such a line therefore stays
        correct if the limits are changed afterward, and it is deliberately excluded from
        autoscaling, so that a construction line passing far from the geometry does not drag the
        view out to meet it. Give `t` to draw a finite piece instead, which is an ordinary line that
        does participate in autoscaling.

        :param lines: the lines to draw.
        :param t: if given, the ``(t0, t1)`` parameter range to draw between, producing a finite
            segment from ``line.at(t0)`` to ``line.at(t1)``. The parameter is in units of the
            line's direction vector, which is arc length only if that vector is unit-length, so
            consider `Line2.normalized` first. If None, the line is drawn infinite.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per line, in the order given. Infinite lines return the ``AxLine``
            subclass of ``Line2D``.
        :raises ValueError: if `t` is given and does not hold exactly two values.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)

        if t is not None:
            values = tuple(t)
            if len(values) != 2:
                raise ValueError(f"t must hold exactly two values, got {len(values)}")
            t0, t1 = values
            drawn = []
            for line in lines:
                a, b = line.at(t0), line.at(t1)
                drawn.append(self.ax.plot([a.x, b.x], [a.y, b.y], **kwargs)[0])
            return drawn

        # An infinite line is a construction entity and should never influence the data limits, but
        # `axline` registers its two reference points like any other artist. Snapshot the limits and
        # put them back, so that a line through a far-off origin cannot rescale the plot.
        saved_bounds = self.ax.dataLim.frozen()
        saved_ignore = self.ax.ignore_existing_data_limits
        try:
            drawn = []
            for line in lines:
                origin, other = line.origin, line.at(1.0)
                drawn.append(self.ax.axline((origin.x, origin.y), (other.x, other.y), **kwargs))
            return drawn
        finally:
            self.ax.dataLim.set(saved_bounds)
            self.ax.ignore_existing_data_limits = saved_ignore

    def draw_boundary(self, *boundaries: Boundary2, tol: float | None = None,
                      color: str | None = None, linewidth: float | None = None,
                      linestyle: str | None = None, alpha: float | None = None,
                      label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more boundaries by flattening their arc segments into polylines.

        :param boundaries: the boundaries to draw.
        :param tol: the maximum deviation between the drawn polyline and the true boundary, in data
            units. Smaller values give smoother arcs at the cost of more points. If None, a
            tolerance of one thousandth of each boundary's own diagonal is used, so mixed-scale
            boundaries in one call are each flattened appropriately.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per boundary, in the order given.
        :raises ValueError: if `tol` is given and is not positive.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
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

    def draw_circle(self, *circles: Circle2 | Iterable[float], fill: bool = False,
                    color: str | None = None, edgecolor: str | None = None,
                    facecolor: str | None = None, linewidth: float | None = None,
                    linestyle: str | None = None, alpha: float | None = None,
                    label: str | None = None, **kwargs) -> list[Circle]:
        """
        Draw one or more circles as patches.

        A raw ``(x, y, r)`` iterable is accepted alongside `Circle2`, so circle data already held as
        an array row does not need converting first.

        :param circles: the circles to draw, each either a `Circle2` or an ``(x, y, r)`` iterable.
        :param fill: whether to fill the circles. Unlike stroking versus filling a curve, both are
            the same Matplotlib patch here, so this is a flag rather than a separate method.
        :param color: sets both the edge and face color at once. If None, the Matplotlib default
            applies. `edgecolor` and `facecolor` override this for the side they name.
        :param edgecolor: the outline color. If None, `color` or the Matplotlib default applies.
        :param facecolor: the fill color, which has no effect unless `fill` is True. If None,
            `color` or the Matplotlib default applies.
        :param linewidth: the outline width in points. If None, the Matplotlib default applies.
        :param linestyle: the outline style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the circle is fully opaque.
        :param label: the legend label. If None, the circle is not labeled.
        :param kwargs: any other keyword argument accepted by the Matplotlib ``Circle`` patch.
        :return: one ``Circle`` patch per circle, in the order given.
        :raises ValueError: if an iterable does not hold exactly three values.
        """
        kwargs = merge_style(kwargs, color=color, edgecolor=edgecolor, facecolor=facecolor,
                             linewidth=linewidth, linestyle=linestyle, alpha=alpha, label=label)
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

    def draw_aabb(self, *boxes: Aabb2, fill: bool = False, color: str | None = None,
                  edgecolor: str | None = None, facecolor: str | None = None,
                  linewidth: float | None = None, linestyle: str | None = None,
                  alpha: float | None = None, label: str | None = None,
                  **kwargs) -> list[Rectangle]:
        """
        Draw one or more axis-aligned bounding boxes as rectangle patches.

        To draw a box with some breathing room around the geometry it bounds, expand it first with
        `Aabb2.expand`, rather than looking for a margin argument here.

        :param boxes: the bounding boxes to draw.
        :param fill: whether to fill the boxes. Both stroking and filling are the same Matplotlib
            patch here, so this is a flag rather than a separate method.
        :param color: sets both the edge and face color at once. If None, the Matplotlib default
            applies. `edgecolor` and `facecolor` override this for the side they name.
        :param edgecolor: the outline color. If None, `color` or the Matplotlib default applies.
        :param facecolor: the fill color, which has no effect unless `fill` is True. If None,
            `color` or the Matplotlib default applies.
        :param linewidth: the outline width in points. If None, the Matplotlib default applies.
        :param linestyle: the outline style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the box is fully opaque.
        :param label: the legend label. If None, the box is not labeled.
        :param kwargs: any other keyword argument accepted by the Matplotlib ``Rectangle`` patch.
        :return: one ``Rectangle`` patch per box, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, edgecolor=edgecolor, facecolor=facecolor,
                             linewidth=linewidth, linestyle=linestyle, alpha=alpha, label=label)
        patches = []
        for box in boxes:
            extent = box.extent
            patch = Rectangle((box.min.x, box.min.y), extent.x, extent.y, fill=fill, **kwargs)
            patches.append(self.ax.add_patch(patch))
        return patches

    def draw_curve(self, *curves: Curve2, color: str | None = None, linewidth: float | None = None,
                   linestyle: str | None = None, alpha: float | None = None,
                   label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more curves as polylines, each as its own separate line.

        :param curves: the curves to draw.
        :param color: the line color. If None, the axes' color cycle supplies one, which means each
            curve in a multi-curve call takes the next cycle color rather than matching.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled. Note that a label given
            while drawing several curves is applied to every one of them, producing repeated legend
            entries; label them in separate calls instead.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per curve, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
        return [
            self.ax.plot(curve.points[:, 0], curve.points[:, 1], **kwargs)[0]
            for curve in curves
        ]

    def draw_spline(self, *splines: CubicSpline2, tol: float, color: str | None = None,
                    linewidth: float | None = None, linestyle: str | None = None,
                    alpha: float | None = None, label: str | None = None,
                    **kwargs) -> list[Line2D]:
        """
        Draw one or more cubic Bezier splines, flattening each into an adaptive polyline whose
        linear segments deviate from the true curve by no more than `tol`.

        :param splines: the splines to draw.
        :param tol: the maximum Euclidean deviation between the drawn polyline and the spline, in
            data units. Smaller values produce a smoother-looking curve at the cost of more line
            segments. Required, and keyword-only, since there is no scale-independent default.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per spline, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
        lines = []
        for spline in splines:
            points = spline.polyline(tol)
            lines.append(self.ax.plot(points[:, 0], points[:, 1], **kwargs)[0])
        return lines

    def fill_curve(self, *curves: Curve2, color: str | None = None, edgecolor: str | None = None,
                   facecolor: str | None = None, linewidth: float | None = None,
                   alpha: float | None = None, label: str | None = None,
                   **kwargs) -> list[Polygon]:
        """
        Fill the interior of one or more curves, rather than stroking their outline.

        This is a separate method from `draw_curve` rather than a flag on it, because filling goes
        through ``Axes.fill`` while stroking goes through ``Axes.plot``, and the two accept
        different keyword arguments.

        :param curves: the curves to fill. Each may be open, in which case it is closed
            automatically by joining its last point back to its first.
        :param color: sets both the edge and face color at once. If None, the axes' color cycle
            supplies one. `edgecolor` and `facecolor` override this for the side they name.
        :param edgecolor: the outline color of the filled region.
        :param facecolor: the interior color of the filled region.
        :param linewidth: the outline width in points. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. Filled regions are usually worth making partly
            transparent so that geometry drawn underneath stays visible.
        :param label: the legend label. If None, the region is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.fill``.
        :return: one ``Polygon`` per curve, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, edgecolor=edgecolor, facecolor=facecolor,
                             linewidth=linewidth, alpha=alpha, label=label)
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
        Draw a `Distance2` as a dimension: leader lines running out to the two anchor points, an
        arrow spanning between them, and a boxed text label carrying the measured value.

        Sizes that are not given in data units, such as the standoff between an anchor point and the
        start of its leader, are derived from the current font size and axis limits. That means the
        dimension is laid out for the axes as they are when this is called, so set the axis limits
        before drawing rather than after.

        :param distance: the measurement to draw.
        :param side_shift: offset the leader line ends by this many data units, orthogonally to the
            measurement direction, with positive values shifting to the right. Use this to keep
            several dimensions of the same feature from overlapping each other.
        :param template: the format string for the label. The measured value is substituted for the
            ``value`` key.
        :param fontsize: the font size of the label.
        :param label_place: where to put the label relative to the anchor points. See `LabelPlace`.
        :param label_offset: how far the label sits from its anchor, in data units. For ``"inside"``
            this is measured from the midpoint, and for the outside placements it is the standoff
            beyond the relevant leader. If None, a default derived from the font size is used.
        :param fontname: the name of the font to use for the label. If None, the Matplotlib default
            is used.
        :param scale_value: a factor applied to the value before it is formatted into the label.
            Use this to display a measurement in different units without altering the geometry or
            its coordinate system.
        :return: every artist drawn, in draw order: the leader and dimension arrows first, then any
            sideways leaders, and the label last.
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
        # A sideways leader is only drawn when the anchor point and the end of its leader are far
        # enough apart to be visually distinct; below that they would overlap into a smudge.
        half_pad = pad * 0.5
        v: Vector2 = leader_end - actual
        if v.norm() < half_pad:
            return []
        work = SurfacePoint2(*actual, *v)
        t1 = work.scalar_projection(leader_end) + half_pad
        return [self.draw_arrow(actual, work.at_distance(t1), arrow="-")]

    def draw_text(self, text: str, pos: Coords2, shift: Coords2 | None = None, ha: str = "center",
                  va: str = "center", color: str | None = None, fontsize: float | None = None,
                  alpha: float | None = None, **kwargs) -> Annotation:
        """
        Draw a text label at a position in data coordinates.

        :param text: the text to draw.
        :param pos: where to place it, as any 2D coordinate the helpers accept.
        :param shift: an optional offset added to `pos`, in data units. Useful for nudging a label
            clear of the feature it names without recomputing the position.
        :param ha: horizontal alignment of the text relative to its position.
        :param va: vertical alignment of the text relative to its position.
        :param color: the text color. If None, the Matplotlib default applies.
        :param fontsize: the font size in points. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the text is fully opaque.
        :param kwargs: any other keyword argument accepted by ``Axes.annotate``, such as ``bbox``,
            ``rotation``, ``fontweight``, or ``fontfamily``.
        :return: the ``Annotation`` that was drawn.
        """
        kwargs = merge_style(kwargs, color=color, fontsize=fontsize, alpha=alpha)
        xy = to_tuple2(pos)
        if shift is not None:
            shift = to_tuple2(shift)
            xy = (xy[0] + shift[0], xy[1] + shift[1])

        return self.ax.annotate(text, xy=xy, ha=ha, va=va, **kwargs)

    def draw_point(self, *points: Coords2, marker: str = "o", markersize: float = 5.0,
                   color: str | None = None, alpha: float | None = None, label: str | None = None,
                   **kwargs) -> Line2D:
        """
        Draw one or more discrete points as markers.

        Unlike the other varargs draw methods, all the points go into a single Matplotlib artist,
        so this returns one ``Line2D`` rather than a list of them.

        :param points: the points to draw, each as any 2D coordinate the helpers accept.
        :param marker: the Matplotlib marker style, such as ``"o"``, ``"x"``, or ``"^"``.
        :param markersize: the marker size in points.
        :param color: the marker color. If None, the axes' color cycle supplies one.
        :param alpha: the opacity from 0.0 to 1.0. If None, the markers are fully opaque.
        :param label: the legend label. If None, the markers are not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``, such as
            ``markeredgecolor`` or ``markerfacecolor``.
        :return: the single ``Line2D`` holding every point. Drawing no points still returns a valid
            but empty artist, so that a computed and possibly empty set needs no special case.
        """
        kwargs = merge_style(kwargs, color=color, alpha=alpha, label=label)
        coords = [to_tuple2(p) for p in points]
        x, y = zip(*coords) if coords else ((), ())
        return self.ax.plot(x, y, marker, markersize=markersize, **kwargs)[0]

    def draw_surface_point(self, *points: SurfacePoint2, arrow_len: float = 1, marker: str = "o",
                           markersize: float = 5.0, color: str = "black",
                           alpha: float | None = None, label: str | None = None,
                           **kwargs) -> list[Artist]:
        """
        Draw one or more surface points, each as a marker with an arrow showing its normal.

        :param points: the surface points to draw.
        :param arrow_len: the length of the normal arrows, in data units.
        :param marker: the Matplotlib marker style used for the point positions.
        :param markersize: the marker size in points.
        :param color: the color of both the markers and their normal arrows, so the two stay
            visually paired. Defaults to black rather than the color cycle, since a surface point
            and its arrow have to match to read as one glyph.
        :param alpha: the opacity from 0.0 to 1.0. If None, the markers are fully opaque.
        :param label: the legend label for the markers. If None, they are not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot`` for the markers.
        :return: every artist drawn: the single ``Line2D`` holding all the markers first, then one
            ``Annotation`` per normal arrow, in the order the points were given.
        """
        kwargs = merge_style(kwargs, color=color, alpha=alpha, label=label)
        coords = [(p.point.x, p.point.y) for p in points]
        x, y = zip(*coords) if coords else ((), ())
        artists: list[Artist] = [self.ax.plot(x, y, marker, markersize=markersize, **kwargs)[0]]
        for p in points:
            p: SurfacePoint2
            artists.append(self.draw_arrow(p.point, p.at_distance(arrow_len), arrow="->", color=color))
        return artists

    def draw_labeled_arrow(self, start: Coords2, end: Coords2, text: str, fraction: float = 0.5,
                           shift: Coords2 | None = None, arrow: str = "->", color: str = "black",
                           linewidth: float | None = None, linestyle: str = "-",
                           **text_kwargs) -> list[Artist]:
        """
        Draw an arrow with a text label positioned along it, for annotating a direction or a
        callout on a diagram.

        :param start: the tail of the arrow, as any 2D coordinate the helpers accept.
        :param end: the head of the arrow.
        :param text: the label text.
        :param fraction: where the label sits along the arrow, as a fraction from `start` to `end`.
            The default of 0.5 puts it at the midpoint; 0.0 and 1.0 put it at the ends.
        :param shift: an optional offset added to the label position, in data units, for nudging
            the text clear of the arrow itself.
        :param arrow: the Matplotlib arrow style. Use ``"-"`` for a plain line with no head.
        :param color: the color of both the arrow and the label, so the two read as one annotation.
        :param linewidth: the width of the arrow. If None, the Matplotlib default is used.
        :param linestyle: the line style of the arrow.
        :param text_kwargs: additional keyword arguments to pass to `draw_text`.
        :return: the arrow ``Annotation`` and the label ``Annotation``, in that order.
        """
        start = Point2(*to_tuple2(start))
        end = Point2(*to_tuple2(end))

        drawn = self.draw_arrow(start, end, arrow=arrow, color=color, linewidth=linewidth, linestyle=linestyle)

        v: Vector2 = end - start
        position = start + v * fraction
        return [drawn, self.draw_text(text, position, shift=shift, color=color, **text_kwargs)]

    def draw_arrow(self, start: Coords2, end: Coords2, arrow: str = "->", color: str = "black",
                   linewidth: float | None = None, linestyle: str = "-",
                   alpha: float | None = None, **arrow_kwargs) -> Annotation:
        """
        Draw an arrow from `start` to `end`.

        This is implemented as an empty annotation carrying an arrow patch, which is what lets the
        head keep a fixed size in points as the axes are zoomed, rather than scaling with the data.

        :param start: the tail of the arrow, as any 2D coordinate the helpers accept.
        :param end: the head of the arrow.
        :param arrow: the Matplotlib arrow style, such as ``"->"``, ``"<->"``, or ``"-"`` for a
            plain line with no head at all.
        :param color: the color of both the arrow's fill and its edge.
        :param linewidth: the width of the arrow. If None, the Matplotlib default is used.
        :param linestyle: the line style of the arrow, such as ``"dotted"``.
        :param alpha: the opacity from 0.0 to 1.0. If None, the arrow is fully opaque.
        :param arrow_kwargs: any other property accepted in the ``arrowprops`` dict passed to
            ``Axes.annotate``, such as ``connectionstyle``, ``shrinkA``, or ``mutation_scale``.
        :return: the ``Annotation`` that was drawn.
        """
        props = merge_style(
            dict(arrowstyle=arrow, fc=color, ec=color, linestyle=linestyle, **arrow_kwargs),
            linewidth=linewidth,
            alpha=alpha,
        )

        return self.ax.annotate("", xy=to_tuple2(end), xytext=to_tuple2(start), arrowprops=props)

    def _font_height(self, font_size: int) -> float:
        # The height of a font in data units, used to size dimension standoffs so that they stay
        # proportional to the text rather than to the geometry being measured.
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
