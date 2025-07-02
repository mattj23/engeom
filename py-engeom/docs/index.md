# Engeom (Python Bindings)

This site has the documentation for the Python bindings of the `engeom` library. 

The [`engeom` library](https://github.com/mattj23/engeom) is a Rust library for 2D and 3D engineering geometry, with a specific focus on metrology applications.  The Python bindings provide a way to use some of the library's functionality in Python, with a reasonably similar interface to the Rust library.

The bindings are built using `pyo3` and `maturin` and `engeom` is compiled as a self-contained extension library, so the python module is the only thing which needs to be installed on an interpreter for everything to work.  It is available for Python versions `>=3.8`, for Windows/Linux/macOS, and on most common 64-bit architectures.

## Quick Overview

The `engeom` Python library has a few general feature sets:

- 2D and 3D geometric primitives (points, vectors, planes, circles/arcs, line segments, transformation matrices, etc) which can be used to represent and perform basic geometric operations through a simple and consistent interface.
- More complicated entities like 3D triangle meshes, 2D and 3D polylines, and 3D point clouds, with common operations like intersections, projections, distance queries, transformations, etc.
- Some basic fitting and alignment algorithms for 2D and 3D data
- Helper functions for plotting and visualization to assist in use of the `matplotlib` and `pyvista` libraries.

## Installation

The Python bindings can be installed using `pip`:

```bash
pip install engeom
```



