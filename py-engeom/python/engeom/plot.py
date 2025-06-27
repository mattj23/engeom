"""
This module contains tools to help with the plotting of geometry and other visual elements using optional plotting
libraries.

Currently, the following plotting libraries are supported:

* **PyVista**: A 3D visualization library that is built on top of VTK and provides tools for both interactive and
  static 3D rendering. The `PyvistaPlotterHelper` can wrap around a PyVista `Plotter` object to provide direct
  functionality for plotting some `engeom` entities.

* **Matplotlib**: A 2D plotting library that is widely used for creating static plots. The `MatplotlibAxesHelper` class
  can wrap around a Matplotlib `Axes` object to provide direct functionality for plotting some `engeom` entities, while
  also forcing the aspect ratio to be 1:1 and expanding the plot to fill all available space.
"""

from ._plot import LabelPlace

try:
    from ._plot.pyvista import PyvistaPlotterHelper
except ImportError:
    pass

try:
    from ._plot.matplotlib import GOM_CMAP, GomColorMap, MatplotlibAxesHelper, TraceBuilder
except ImportError:
    pass
