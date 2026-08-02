"""
Helpers for visualizing `engeom` entities with optional third-party plotting libraries.

Each backend lives in its own submodule, so the import statement is the declaration of which
optional dependency you need:

* `engeom.plot.matplotlib` requires `matplotlib` (`pip install engeom[matplotlib]`). It wraps a
  Matplotlib `Axes` object for drawing 2D entities, and provides a parallel-projection viewport for
  drawing 3D entities onto that same 2D axes, such as for generating diagrams and illustrations.

* `engeom.plot.pyvista` requires `pyvista` (`pip install engeom[pyvista]`). It wraps a PyVista
  `Plotter` object, which is built on VTK and provides both interactive and static 3D rendering.

Importing `engeom.plot` itself pulls in neither backend, so it is always safe. Importing a backend
submodule without its library installed raises an `ImportError` naming the missing package.
"""

from ._common import LabelPlace

__all__ = ["LabelPlace"]
