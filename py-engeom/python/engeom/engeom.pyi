from __future__ import annotations
from enum import Enum

type ResampleEnum = Resample_Count | Resample_Spacing | Resample_MaxSpacing

class DeviationMode(Enum):
    """
    Represents the different methods of calculating deviation between a point and another geometry.
    """
    Point = 0
    Plane = 1

class SelectOp(Enum):
    """
    A common enum to describe different types of selection operations.
    """
    Add=0
    Remove=1
    Keep=2

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