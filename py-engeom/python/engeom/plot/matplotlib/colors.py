"""
Color maps and deviation scales for use with the Matplotlib helpers.
"""

from __future__ import annotations

import inspect

import numpy
from matplotlib.colors import Colormap, ListedColormap, Normalize
from numpy.typing import ArrayLike

# Matplotlib 3.11 added `bad`, `under`, and `over` constructor arguments and started raising a
# PendingDeprecationWarning from `set_under`/`set_over`. The package declares no matplotlib floor,
# so which path is available has to be detected rather than assumed.
_CTOR_TAKES_EXTREMES = "over" in inspect.signature(ListedColormap.__init__).parameters


class GomColorMap(ListedColormap):
    """
    A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software.

    Values beyond the mapped range are called out rather than clamped: magenta below and dark red
    above. That is the convention the GOM software uses, and it is what makes an out-of-tolerance
    deviation obvious instead of merging into the end color.

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

        if _CTOR_TAKES_EXTREMES:
            super().__init__(colors, under="magenta", over="darkred")
        else:
            super().__init__(colors)
            self.set_under("magenta")
            self.set_over("darkred")


GOM_CMAP = GomColorMap()
"""
A color map similar to the 8 discrete colors in the GOM/Zeiss Inspect software, already instantiated and
available in the module.
"""


def deviation_limit(values: ArrayLike, percentile: float = 100.0) -> float:
    """
    Choose a symmetric half-range for a deviation scale from the data itself.

    The result is the magnitude that `percentile` of the absolute values fall within, so the
    returned limit is always non-negative and is meant to be fed straight to `deviation_norm`. Lower
    the percentile to keep a handful of outliers from flattening the whole scale; at the default of
    100 the limit is simply the largest absolute deviation present.

    :param values: the deviation values, in any shape numpy will accept. NaNs are ignored.
    :param percentile: the percentile of absolute deviation the limit should cover, from 0 to 100.
    :return: the half-range to normalize about zero with.
    :raises ValueError: if `percentile` is outside 0 to 100, or if `values` holds no finite value.
    """
    if not 0.0 <= percentile <= 100.0:
        raise ValueError(f"percentile must be between 0 and 100, got {percentile}")

    magnitudes = numpy.abs(numpy.asarray(values, dtype=numpy.float64))
    magnitudes = magnitudes[numpy.isfinite(magnitudes)]
    if magnitudes.size == 0:
        raise ValueError("values holds no finite value to take a deviation limit from")

    return float(numpy.percentile(magnitudes, percentile))


def deviation_norm(limit: float) -> Normalize:
    """
    Build a normalization for a deviation scale, symmetric about zero.

    Symmetry is the whole point of this over calling ``Normalize`` directly. A deviation scale fit
    to the data's own minimum and maximum puts zero somewhere off center, which silently moves the
    color that reads as "nominal" and makes two plots of the same part incomparable.

    Used with `GOM_CMAP`, the range is divided into the map's eight discrete bands, and values
    outside it take the map's over and under colors.

    :param limit: the half-range, so the scale runs from ``-limit`` to ``+limit``. The sign is
        ignored, since a range is symmetric either way.
    :return: a ``Normalize`` centered on zero.
    :raises ValueError: if `limit` is zero, which would leave nothing to normalize across.
    """
    if limit == 0.0:
        raise ValueError("limit must be non-zero, or the scale has no range to divide")

    magnitude = abs(limit)
    return Normalize(vmin=-magnitude, vmax=magnitude)


def has_extremes(cmap: Colormap) -> tuple[bool, bool]:
    """
    Report whether a color map calls out values below and above its range with distinct colors,
    as `GOM_CMAP` does.

    Matplotlib exposes no direct way to ask this, so the end colors are compared against the map's
    own first and last colors; a map whose extreme color matches its end color is not calling
    anything out. This is what `AxesHelper.draw_colorbar` uses to decide whether to draw the
    pointed ends that mean "there is data beyond here".

    :param cmap: the color map to inspect.
    :return: whether the under color and the over color, in that order, are distinct from the ends.
    """
    under = not numpy.array_equal(numpy.asarray(cmap.get_under()), numpy.asarray(cmap(0.0)))
    over = not numpy.array_equal(numpy.asarray(cmap.get_over()), numpy.asarray(cmap(1.0)))
    return under, over


def extend_for(cmap: Colormap) -> str:
    """
    Return the Matplotlib colorbar ``extend`` token matching a color map's out-of-range colors.

    :param cmap: the color map to inspect.
    :return: one of ``"neither"``, ``"min"``, ``"max"``, or ``"both"``.
    """
    under, over = has_extremes(cmap)
    if under and over:
        return "both"
    if under:
        return "min"
    if over:
        return "max"
    return "neither"
