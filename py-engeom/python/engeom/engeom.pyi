from __future__ import annotations
from typing import Literal

type ResampleEnum = Resample.Count | Resample.Spacing | Resample.MaxSpacing

VecDot = Literal["as_is", "abs", "clamp_pos"]
"""Controls how algorithms use a dot product between two direction vectors: ``"as_is"`` (raw value,
can be negative), ``"abs"`` (absolute value), or ``"clamp_pos"`` (clamped to zero from below)."""

DeviationMode = Literal["point", "plane"]
"""How deviation between a point and another geometry is measured: ``"point"`` (direct distance to
the closest point) or ``"plane"`` (distance to the tangent plane at the closest point)."""

SelectOp = Literal["add", "remove", "keep"]
"""A selection operation applied to a set: ``"add"``, ``"remove"``, or ``"keep"`` (keep only)."""


class Resample:
    """
    A common enum to describe different methods of resampling a point based geometry.
    """

    class Count:
        """ A resampling method that specifies a fixed number of entities. """
        def __init__(self, count: int):
            self.count = count

    class Spacing:
        """
        A resampling method that specifies a fixed spacing between entities. For some types of resampling operations,
        this may result in dead space near edges. Check with the specific operation documentation for details.
        """
        def __init__(self, spacing: float):
            self.spacing = spacing

    class MaxSpacing:
        """
        A resampling method which specifies a maximum permissible spacing between entities. Unlike `Spacing`, this will
        allow the operation to adjust spacing to remove potential dead space near edges, while also placing a limit on
        how far apart entities can be as a result.
        """
        def __init__(self, max_spacing: float):
            self.max_spacing = max_spacing