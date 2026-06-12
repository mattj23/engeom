"""
Tools for analyzing 2D airfoil cross-sections (new implementation).
"""

from ..engeom import _airfoil2

for name in [n for n in dir(_airfoil2) if not n.startswith("_")]:
    globals()[name] = getattr(_airfoil2, name)
