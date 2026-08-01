"""
Coercion of the loose "point-like" arguments the plotting helpers accept into concrete `engeom`
types or plain tuples.

The helpers deliberately take anything that reads as a coordinate, so that callers can pass an
`engeom` primitive, a tuple, a list, or a numpy row interchangeably without converting first. This
module is where that permissiveness is implemented once, for every backend.
"""

from __future__ import annotations

from typing import Iterable, Tuple, Union

from engeom.geom2 import Point2, Vector2
from engeom.geom3 import Point3

# A 2D coordinate in any of the forms the helpers accept. Kept as `Union` rather than `X | Y`
# because a module-level alias is evaluated eagerly, and the package supports Python 3.8.
Coords2 = Union[Point2, Vector2, Iterable[float]]

# A coordinate of either dimension. Used where the helper promotes or truncates between 2D and 3D
# rather than requiring the caller to match the dimension of the target.
PointLike = Union[Point2, Tuple[float, float], Point3, Tuple[float, float, float], Iterable[float]]


def to_tuple2(item: Coords2) -> Tuple[float, float]:
    """
    Coerce a 2D coordinate into a plain ``(x, y)`` tuple.

    :param item: a `Point2`, a `Vector2`, or any iterable whose first two values are the x and y
        coordinates. Additional values are ignored.
    :return: the x and y coordinates as a tuple of floats.
    """
    if isinstance(item, (Point2, Vector2)):
        return item.x, item.y
    else:
        x, y, *_ = item
        return x, y


def to_point2(p: PointLike) -> Point2:
    """
    Coerce a coordinate of either dimension into a `Point2`, truncating the z coordinate of a 3D
    input and padding a short iterable with zeros.

    :param p: the coordinate to convert.
    :return: a `Point2`.
    """
    if isinstance(p, Point2):
        return p
    elif isinstance(p, Point3):
        return Point2(p.x, p.y)
    else:
        values = [float(v) for v in p]
        while len(values) < 2:
            values.append(0)
        return Point2(*values)


def to_point3(p: PointLike) -> Point3:
    """
    Coerce a coordinate of either dimension into a `Point3`, placing a 2D input on the z=0 plane and
    padding a short iterable with zeros.

    :param p: the coordinate to convert.
    :return: a `Point3`.
    """
    if isinstance(p, Point3):
        return p
    elif isinstance(p, Point2):
        return Point3(p.x, p.y, 0)
    else:
        values = [float(v) for v in p]
        while len(values) < 3:
            values.append(0)
        return Point3(*values)
