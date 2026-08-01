"""
Helpers for drawing `engeom` entities onto a Matplotlib `Axes` object.

`MatplotlibAxesHelper` wraps an `Axes` and draws 2D entities directly onto it. For 3D entities, ask
the helper for a `ViewPort` via `get_3d_viewport`, which renders them in parallel projection onto
the same axes.

Requires `matplotlib`, which is an optional dependency: `pip install engeom[matplotlib]`.
"""

try:
    # Probe for the dependency before importing the submodules, so that a missing matplotlib is
    # reported as itself rather than as a confusing failure deeper in the package. Probing
    # `pyplot` specifically because that is what the submodules import.
    import matplotlib.pyplot  # noqa: F401
except ImportError as _e:
    raise ImportError(
        "engeom.plot.matplotlib requires the `matplotlib` package, which is not installed. "
        "Install it with `pip install engeom[matplotlib]`."
    ) from _e

from .axes import MatplotlibAxesHelper
from .colors import GOM_CMAP, GomColorMap
from .trace import TraceBuilder
from .viewport import ViewPort

__all__ = [
    "GOM_CMAP",
    "GomColorMap",
    "MatplotlibAxesHelper",
    "TraceBuilder",
    "ViewPort",
]
