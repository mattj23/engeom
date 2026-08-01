"""
Accumulation of disconnected polylines into the single flat x/y sequence that Matplotlib's plotting
functions consume, using `None` as the break between runs.
"""

from __future__ import annotations

from engeom.geom2 import Aabb2

from .._coerce import PointLike, to_tuple2


class TraceBuilder:
    """
    Accumulates points into flat parallel lists suitable for handing to Matplotlib in a single call.

    Matplotlib draws a polyline through every point in one `Axes.plot` call, so drawing many
    disconnected runs would normally require one call per run. Inserting a `None` breaks the line
    instead, letting any number of disjoint segments go out as a single artist. This class manages
    that bookkeeping.

    The three lists are always the same length and are index-aligned:

    * `xs` and `ys` hold the coordinates, with `None` at a break between runs.
    * `c` holds a scalar per point for color mapping, or `None` where the point was added without
      one. It is only meaningful if every point was added through `add_point_and_color`.

    !!! example
        ```python
        builder = TraceBuilder()
        builder.add_segment(p0, p1)
        builder.add_segment(p2, p3)
        ax.plot(*builder.xy, color="black")
        ```
    """

    def __init__(self):
        self.xs = []
        self.ys = []
        self.c = []

    def bounds(self) -> Aabb2 | None:
        """
        Compute the axis-aligned bounding box of the points accumulated so far, ignoring the `None`
        entries that mark breaks between runs.

        :return: an `Aabb2` covering every point, or `None` if nothing has been added yet or if only
            breaks have been added.
        """
        xs = [x for x in self.xs if x is not None]
        ys = [y for y in self.ys if y is not None]
        if not xs or not ys:
            return None
        return Aabb2(
            x_min=min(xs),
            x_max=max(xs),
            y_min=min(ys),
            y_max=max(ys),
        )

    def add_segment(self, *points: PointLike):
        """
        Add a run of points followed by a break, so that the next thing added starts a new
        disconnected polyline.

        :param points: the points in the run, in order. Each may be any 2D or 3D coordinate the
            plotting helpers accept; a z coordinate is ignored.
        """
        self.add_points(*points)
        self.add_blank()

    def add_blank(self):
        """
        Add a break, terminating the current run so that the next point added is not connected to
        the previous one. Appends `None` to all three lists.
        """
        self.xs.append(None)
        self.ys.append(None)
        self.c.append(None)

    def add_points(self, *points: PointLike):
        """
        Add one or more points to the current run without terminating it.

        Points added this way carry no color, so `None` is appended to `c` for each one to keep the
        three lists index-aligned.

        :param points: the points to add. Each may be any 2D or 3D coordinate the plotting helpers
            accept; a z coordinate is ignored.
        """
        for point in points:
            x, y = to_tuple2(point)
            self.xs.append(x)
            self.ys.append(y)
            self.c.append(None)

    def add_point_and_color(self, point: PointLike, color: float):
        """
        Add a single point to the current run along with a scalar to color it by, for use with a
        color-mapped plotting call.

        :param point: the point to add. May be any 2D or 3D coordinate the plotting helpers accept;
            a z coordinate is ignored.
        :param color: the scalar value to associate with this point.
        """
        x, y = to_tuple2(point)
        self.xs.append(x)
        self.ys.append(y)
        self.c.append(color)

    def invert_y(self):
        """
        Negate every y coordinate accumulated so far, in place, leaving the `None` breaks intact.
        Useful when the source data uses a downward-positive y axis, as image and raster data does.
        """
        self.ys = [-y if y is not None else None for y in self.ys]

    @property
    def xy(self):
        """
        The accumulated coordinates as an `(xs, ys)` pair, ready to be unpacked into an `Axes.plot`
        call as ``ax.plot(*builder.xy, ...)``.

        :return: a tuple of the x list and the y list.
        """
        return self.xs, self.ys
