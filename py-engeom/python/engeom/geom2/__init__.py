"""
This module contains a number of classes and functions for working with 2D geometry. These range from fundamental
primitives like points and vectors up to complex items like polylines.
"""

from ..engeom import _geom2

# Global import of all functions
for name in [n for n in dir(_geom2) if not n.startswith("_")]:
    globals()[name] = getattr(_geom2, name)
