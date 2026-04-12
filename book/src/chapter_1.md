# Introduction to Engeom

`engeom` is a Rust library for working with engineering geometry, with a particular focus on metrology and dimensional
inspection. It provides a comprehensive set of tools for constructing, measuring, and transforming 2D and 3D geometric
entities, and is designed to handle the full pipeline from raw measurement data through to final inspection results.

Python bindings are published to PyPI and expose most of the library's functionality as a compiled extension module:

```bash
pip install engeom
```

This book covers both the Rust API and the Python bindings together. Where the two interfaces differ, each is shown
separately; where they are equivalent the Rust example is given first.

## Design Priorities

Because `engeom` is a metrology library, algorithms and implementations are prioritized in the following order:

1. **Accuracy and correctness** - results must be numerically reliable for inspection use cases
2. **Speed** - performance matters, but never at the cost of correctness
3. **Memory usage** - efficient, but a secondary concern

## Capabilities

The library is organized around several interconnected domains:

**3D Geometry**
- Points, vectors, surface points, lines, planes, spheres, and meshes
- Measurements on point clouds and unstructured triangle meshes
- Levenberg-Marquardt best-fit and rigid-body alignment

**2D Geometry**
- Points, vectors, surface points, lines, circles, arcs, and polylines
- Curve fitting, alignment, and distance/angle measurement
- Specialized tools for airfoil cross-section construction and analysis

**Raster Fields (2D)**
- Depth-map style scalar fields on a regular grid
- Binning, filtering, smoothing, and in-painting operations

**1D Scalar Series**
- Spatial series sampled from 2D/3D surfaces or motion time-series
- Interpolation, smoothing, filtering, extrema detection, and curve fitting

**Cross-domain Transformations**
- 3D to 2D projections and mesh unrolling
- Sampling 3D/2D data onto 1D series
- Projecting deviations between domains

## Using the Rust Library

Add `engeom` to your `Cargo.toml`:

```toml
[dependencies]
engeom = "0.3"
```

The most commonly used types are re-exported directly from the crate root:

```rust
use engeom::{Point3, Vector3, SurfacePoint3, Plane3, Mesh};
use engeom::{Point2, Vector2, SurfacePoint2, Circle2, Curve2};
```

The library builds on [nalgebra](https://docs.rs/nalgebra/latest/nalgebra/) (via the `parry` physics libraries) for its core linear
algebra types. `Point2`, `Vector2`, `Point3`, `Vector3`, and the `Iso` isometry types are all nalgebra type aliases,
so the full nalgebra API is available on them wherever `engeom` is used.

## Using the Python Bindings

After installing the package, import from the `engeom.geom2` or `engeom.geom3` sub-modules:

```python
from engeom.geom2 import Point2, Vector2, SurfacePoint2, Circle2
from engeom.geom3 import Point3, Vector3, SurfacePoint3, Plane3, Mesh
```

> **Working with large datasets:** For operations on many points or vectors, prefer `numpy` arrays over Python lists
> of `Point`/`Vector` objects. Many `engeom` functions accept `numpy.ndarray` inputs directly, which is significantly
> faster than constructing individual Python objects.

## NumPy Integration

The Python bindings are designed to work seamlessly with [NumPy](https://numpy.org/), Python's standard library for
numerical array operations. Because `engeom` is frequently used to process large sets of measurement points, passing
data as `numpy.ndarray` objects is much more efficient than constructing individual `Point` or `Vector` instances in a
loop.

The following conventions apply throughout the library:

**Array shape for points and vectors.** When passing or receiving arrays of 2D or 3D points and vectors, arrays are
shaped `(n, 2)` or `(n, 3)` respectively. Each row is one point or vector, and the columns correspond to the `x`,
`y`, and `z` components.

```python
import numpy as np

# Ten random 3D points
points = np.random.rand(10, 3)  # shape (10, 3)

# Ten random 2D vectors
vectors = np.random.rand(10, 2)  # shape (10, 2)
```

**Read-only array access.** Some `engeom` objects expose their internal data directly as read-only `numpy.ndarray`
views. For example, `Mesh` provides its vertex positions and face index lists as arrays without copying the underlying
data. These arrays follow the same `(n, d)` shape convention, and floating-point data is always 64-bit (`float64`).

**Integer index arrays.** Where index lists are used - such as the face-vertex lists of a `Mesh` - the expected dtype
is an unsigned integer. Internally these are stored as Rust `u32` or `usize` values. If you encounter a construction
error when passing index data, check that the array dtype is an unsigned integer type (e.g. `numpy.uint32`) rather
than the default signed `int64`.

```python
import numpy as np

vertices = np.array([[0.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0]], dtype=np.float64)

faces = np.array([[0, 1, 2]], dtype=np.uint32)  # must be unsigned
```
