"""
Handling of the common Matplotlib styling arguments that the draw methods name explicitly.

The draw methods spell out the handful of styling arguments that get used constantly (color, line
width, and so on) as real keyword parameters, so that editors and type checkers can complete them,
while still forwarding an open ``**kwargs`` for everything else Matplotlib accepts. Naming them has
one consequence: an argument the caller did not supply arrives as `None`, and passing that straight
through would override Matplotlib's own default rather than defer to it. `merge_style` is where that
distinction is enforced.
"""

from __future__ import annotations

from typing import Any, Dict


def merge_style(kwargs: Dict[str, Any], **named: Any) -> Dict[str, Any]:
    """
    Merge explicitly named styling arguments into a keyword argument dict, dropping any that were
    left as `None`.

    Dropping them matters because Matplotlib's defaults are not all fixed values. Passing
    ``color=None`` to ``Axes.plot`` is not the same as omitting it, since omitting it draws the next
    color from the axes' color cycle.

    There is no risk of a collision between the two sources: an argument named in the signature is
    bound to that parameter and never reaches ``**kwargs``, so Python rejects a duplicate before
    this is called.

    :param kwargs: the open keyword arguments collected by the draw method.
    :param named: the explicitly named styling arguments, any of which may be `None`.
    :return: the `kwargs` dict, mutated in place and returned for convenience.
    """
    for key, value in named.items():
        if value is not None:
            kwargs[key] = value
    return kwargs
