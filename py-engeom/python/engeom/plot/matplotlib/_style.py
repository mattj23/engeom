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

from typing import Any, Dict, Optional, Union

# How a single element of a composite drawing is styled or suppressed. `False` hides the element,
# `True` or `None` draws it with the method's own defaults, and a dict of Matplotlib keyword
# arguments draws it with those merged over the defaults.
#
# This is a module-level alias, so it is evaluated eagerly and has to spell out `Union` and
# `Optional` rather than using `X | Y`, which the package's Python 3.8 floor does not accept at
# runtime. Annotations inside function signatures are exempt, since `from __future__ import
# annotations` keeps those unevaluated.
ElementStyle = Optional[Union[bool, Dict[str, Any]]]


def element_style(value: ElementStyle, defaults: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Resolve a per-element style argument from a composite draw method into the keyword arguments to
    draw that element with, or `None` if the caller asked for it to be left out.

    A supplied dict is merged *over* the defaults rather than replacing them, so that restyling one
    thing about an element does not silently discard the rest of its designed appearance. Pass the
    key explicitly to override a default; there is no way to unset one, since Matplotlib's own
    default is what a `None` value would mean and that is a different thing again.

    :param value: the caller's argument. `False` suppresses the element, `True` or `None` accepts
        the defaults, and a dict is merged over them.
    :param defaults: the element's default keyword arguments, which are not mutated.
    :return: the keyword arguments to draw with, or `None` if the element is suppressed.
    """
    if value is False:
        return None
    if value is None or value is True:
        return dict(defaults)
    return {**defaults, **value}


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
