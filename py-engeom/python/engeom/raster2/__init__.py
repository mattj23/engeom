"""
This module provides a number of classes and functions for working with 2D raster (pixel) data.
"""

from ..engeom import _raster2

# Global import of all functions
for name in [n for n in dir(_raster2) if not n.startswith("_")]:
    globals()[name] = getattr(_raster2, name)
