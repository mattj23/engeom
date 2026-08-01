"""
Color maps for use with the Matplotlib helpers.
"""

from __future__ import annotations

import numpy
from matplotlib.colors import ListedColormap


class GomColorMap(ListedColormap):
    """
    A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software.

    You can use this to instantiate a color map, or you can use the `GOM_CMAP` object directly.
    """

    def __init__(self):
        colors = numpy.array(
            [
                [1, 0, 160],
                [1, 0, 255],
                [0, 254, 255],
                [0, 160, 0],
                [0, 254, 0],
                [255, 255, 0],
                [255, 128, 0],
                [255, 1, 0],
            ],
            dtype=numpy.float64,
        )
        colors /= 256.0
        colors = numpy.hstack((colors, numpy.ones((len(colors), 1))))
        super().__init__(colors)
        self.set_under("magenta")
        self.set_over("darkred")


GOM_CMAP = GomColorMap()
"""
A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software, already instantiated and
available in the module.
"""
