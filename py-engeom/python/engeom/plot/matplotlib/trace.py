"""
Accumulation of disconnected polylines into the single flat x/y sequence that Matplotlib's plotting
functions consume, using `None` as the break between runs.
"""

from __future__ import annotations

from engeom.geom2 import Aabb2

from .._coerce import PointLike


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
