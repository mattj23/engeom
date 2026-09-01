"""
Giving an unbounded entity something finite to draw.

A `Plane3` has no origin, in-plane orientation, or size, while a `Line3` has only a direction and a
point on it. Neither can be drawn without deciding how much to show. These helpers remove the need
for callers to work that out every time by resolving an extent from a convenient source—usually
the geometry already on the plotter—and clipping the entity to it.
"""

from __future__ import annotations

import numpy

from engeom.geom3 import Aabb3, Line3, Plane3, Point3

from .._common import plane_basis
from .._coerce import to_tuple3

__all__ = ["clip_line_to_aabb", "clip_plane_to_aabb", "resolve_extent"]

# The renderer bounds PyVista reports when there is nothing in it. This is a placeholder rather
# than a measurement, so an extent taken from an empty scene has to be refused rather than silently
# sizing later entities against a two-unit cube.
_EMPTY_BOUNDS = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)


def _box_of(value) -> Aabb3:
    """Resolve everything except `None` into a bounding box, without padding it."""
    if isinstance(value, Aabb3):
        return value

    # `Mesh3`, `Segment3`, `CurveGroup3` and `CubicSpline3` carry their box as a property, while a
    # `PointCloud3` computes one on demand, so both spellings are accepted.
    box = getattr(value, "aabb", None)
    if isinstance(box, Aabb3):
        return box
    if hasattr(value, "compute_aabb"):
        return value.compute_aabb()

    if numpy.ndim(value) == 2:
        return Aabb3.from_points(numpy.asarray(value, dtype=numpy.float64).reshape(-1, 3))

    raise TypeError(
        f"Cannot take an extent from {type(value).__name__}. Give an `Aabb3`, an entity with a "
        "bounding box, an (n, 3) array of points, or None to use the scene."
    )


def resolve_extent(value, plotter=None, pad: float = 0.1) -> Aabb3:
    """
    Work out the region of space an unbounded entity should be drawn across.

    :param value: what to take the extent from. `None` uses the bounds of everything already drawn
        into the plotter's active renderer, which is usually what is wanted: draw the part, then
        draw the plane through it. Otherwise an `Aabb3`, any entity carrying a bounding box
        (`Mesh3`, `Segment3`, `CurveGroup3`, `CubicSpline3`, `PointCloud3`), or an `(n, 3)` array
        of points.
    :param plotter: the plotter to read the scene bounds from. Only needed when `value` is `None`.
    :param pad: how far to grow the box, as a fraction of its diagonal, so that a plane drawn
        through a part protrudes past it rather than stopping flush with its surface.
    :return: the padded bounding box to draw across.
    :raises ValueError: if `value` is `None` and the scene has nothing bounded in it, or if the
        resolved box has no size to it.
    :raises TypeError: if `value` is not one of the accepted forms.
    """
    if value is None:
        if plotter is None:
            raise ValueError("An extent of None has to come from a plotter, but none was given")

        renderer = plotter.renderer
        bounds = tuple(float(v) for v in plotter.bounds)
        if not renderer.actors or bounds == _EMPTY_BOUNDS:
            raise ValueError(
                "There is nothing in the scene to take an extent from. Draw the geometry the "
                "entity belongs to first, or pass `extent=` explicitly."
            )
        box = Aabb3(*bounds[0:6:2], *bounds[1:6:2])
    else:
        box = _box_of(value)

    diagonal = float(numpy.linalg.norm(box.extent.as_numpy()))
    if diagonal <= 0.0:
        raise ValueError("The extent has no size to it, so there is nothing to draw across")

    return box.expand(diagonal * pad) if pad else box


def _corners(box: Aabb3) -> tuple[numpy.ndarray, numpy.ndarray]:
    lo = numpy.array(to_tuple3(box.min), dtype=numpy.float64)
    hi = numpy.array(to_tuple3(box.max), dtype=numpy.float64)
    return lo, hi


def _clip_to_slab(polygon: list, axis: int, value: float, keep_above: bool) -> list:
    """
    Clip a convex polygon against one of a box's six bounding planes, by Sutherland-Hodgman.

    Cutting a large in-plane square down against the six planes in turn is what makes the plane
    case free of special cases. Intersecting the box's twelve edges with the plane instead has to
    cope with a plane through a corner, which produces the same point three times, and with a plane
    lying in a face, where four edges give no crossing at all.
    """
    if not polygon:
        return []

    result = []
    for i, current in enumerate(polygon):
        following = polygon[(i + 1) % len(polygon)]
        current_in = current[axis] >= value if keep_above else current[axis] <= value
        following_in = following[axis] >= value if keep_above else following[axis] <= value

        if current_in:
            result.append(current)
        if current_in != following_in:
            # The two points fall on opposite sides, so they differ along this axis and the
            # denominator cannot vanish.
            t = (value - current[axis]) / (following[axis] - current[axis])
            result.append(current + t * (following - current))
    return result


def clip_plane_to_aabb(plane: Plane3, box: Aabb3) -> numpy.ndarray | None:
    """
    Find the convex polygon where a plane crosses a box.

    :param plane: the plane to clip.
    :param box: the box to clip it against.
    :return: an `(n, 3)` array of the polygon's vertices, wound consistently, or `None` if the
        plane misses the box entirely.
    """
    lo, hi = _corners(box)
    center = numpy.array(to_tuple3(plane.project_point(box.center)), dtype=numpy.float64)

    # Start from a square in the plane that is certainly larger than the box, so that clipping can
    # only ever cut it down to the true intersection.
    u, v = plane_basis(plane.normal)
    reach = float(numpy.linalg.norm(hi - lo))
    du = numpy.array(to_tuple3(u), dtype=numpy.float64) * reach
    dv = numpy.array(to_tuple3(v), dtype=numpy.float64) * reach
    polygon = [center + du + dv, center - du + dv, center - du - dv, center + du - dv]

    for axis in range(3):
        polygon = _clip_to_slab(polygon, axis, lo[axis], keep_above=True)
        polygon = _clip_to_slab(polygon, axis, hi[axis], keep_above=False)

    return numpy.array(polygon, dtype=numpy.float64) if len(polygon) >= 3 else None


def clip_line_to_aabb(line: Line3, box: Aabb3) -> numpy.ndarray | None:
    """
    Find the segment of a line that lies inside a box.

    :param line: the line to clip. Its direction need not be normalized.
    :param box: the box to clip it against.
    :return: a `(2, 3)` array of the segment's endpoints, or `None` if the line misses the box.
    :raises ValueError: if the line's direction has no length, so it points nowhere.
    """
    lo, hi = _corners(box)
    origin = numpy.array(to_tuple3(line.origin), dtype=numpy.float64)
    direction = numpy.array(to_tuple3(line.direction), dtype=numpy.float64)

    if not numpy.any(direction):
        raise ValueError("The line has no direction, so it cannot be clipped to an extent")

    # The usual slab test: keep the overlap of the parameter intervals over which the line is
    # inside each of the three pairs of parallel faces.
    near, far = -numpy.inf, numpy.inf
    for axis in range(3):
        if direction[axis] == 0.0:
            if not (lo[axis] <= origin[axis] <= hi[axis]):
                return None
            continue
        first = (lo[axis] - origin[axis]) / direction[axis]
        second = (hi[axis] - origin[axis]) / direction[axis]
        if first > second:
            first, second = second, first
        near = max(near, first)
        far = min(far, second)
        if near > far:
            return None

    return numpy.array([origin + near * direction, origin + far * direction], dtype=numpy.float64)
