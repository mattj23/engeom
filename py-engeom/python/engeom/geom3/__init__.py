"""
This module provides a number of classes and functions for working with 3D geometry. These range from fundamental
primitives like points and vectors up to complex items like meshes.
"""

from ..engeom import _geom3

# Global import of all functions
for name in [n for n in dir(_geom3) if not n.startswith("_")]:
    globals()[name] = getattr(_geom3, name)
