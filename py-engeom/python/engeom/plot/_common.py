"""
Backend-neutral pieces shared by every plotting helper.
"""

from __future__ import annotations

from typing import Literal, Tuple, get_args

from engeom.geom3 import Vector3

LabelPlace = Literal["outside", "inside", "outside_rev"]
"""
Where a measurement's label is placed relative to its two anchor points:

* ``"outside"`` places the label beyond the anchors, on the side of the second point.
* ``"inside"`` places the label between the two anchor points.
* ``"outside_rev"`` places the label beyond the anchors, on the side of the first point.
"""

LABEL_PLACES: Tuple[str, ...] = get_args(LabelPlace)
"""The valid `LabelPlace` tokens, derived from the alias so the two cannot drift apart."""


def check_label_place(value: str) -> str:
    """
    Validate a label placement token, so that an unknown value fails immediately with a message
    naming the alternatives rather than falling through the placement branches.

    :param value: the token to check.
    :return: the token unchanged, if it is valid.
    :raises ValueError: if the token is not one of the valid placements.
    """
    if value not in LABEL_PLACES:
        valid = ", ".join(repr(v) for v in LABEL_PLACES)
        raise ValueError(f"invalid label_place {value!r}, expected one of {valid}")
    return value


def plane_basis(normal: Vector3) -> Tuple[Vector3, Vector3]:
    """
    Build an orthonormal pair of vectors spanning the plane with the given normal.

    The pair is deterministic but otherwise arbitrary, since a plane normal alone does not fix a
    rotation about itself. The world axis least aligned with the normal is used as the seed, which
    keeps the construction numerically well conditioned for every possible normal.

    This is shared rather than per-backend because both backends need it for the same reason: an
    entity defined only by a center, a normal, and a size (a `Circle3`, a `Plane3`) has to be given
    an in-plane orientation before it can be turned into vertices.

    :param normal: the plane normal. Need not be unit length.
    :return: two orthonormal vectors, both perpendicular to `normal` and to each other.
    """
    n = normal.normalized()
    seed = min((Vector3.x_axis(), Vector3.y_axis(), Vector3.z_axis()), key=lambda a: abs(n.dot(a)))
    u = n.cross(seed).normalized()
    return u, n.cross(u)
