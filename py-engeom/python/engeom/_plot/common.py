from enum import Enum

class LabelPlace(Enum):
    """
    Represents the different locations where a label can be placed between its anchor points.
    """

    Outside = 1
    """ The label is placed outside the anchor points, on the side of the second point in the measurement. """

    Inside = 2
    """ The label is placed between the two anchor points. """

    OutsideRev = 3
    """ The label is placed outside the two anchor points, on the side of the first point in the measurement. """


