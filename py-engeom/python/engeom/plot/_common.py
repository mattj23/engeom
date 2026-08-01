"""
Backend-neutral pieces shared by every plotting helper.
"""

from __future__ import annotations

from typing import Literal, Tuple, get_args

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
