# engeom

A Rust library of metrology-focused tools for working with 2D and 3D geometry. The primary use case is
engineering applications such as GD&T and quality/dimensional inspection, covering every step of the process
starting from raw data formats such as meshes and point clouds.

Because this is a metrology-focused library, the order of priority for algorithms and methods is:

1. Accuracy/correctness
2. Speed
3. Memory usage

## What is included

* **3D geometry** (`geom3`): points, vectors, surface points, lines, planes, circles, spheres, cylinders, cones,
  triangle meshes, and point clouds, with Levenberg-Marquardt fitting and alignment (`align3`).
* **2D geometry** (`geom2`): points, vectors, surface points, lines, circles, arcs, and polyline curves, with
  fitting and alignment (`align2`).
* **Airfoil analysis** (`airfoil2`): construction and analysis of airfoil cross-sections, including camber
  lines, leading/trailing edge detection, inscribed circles, thickness, and orientation.
* **Raster fields** (`raster2`, `raster3`): 2D scalar fields for depth-map style data, with filtering,
  smoothing, in-painting, and region labeling.
* **1D scalar series** (`func1`): functions over a single variable domain, with interpolation, smoothing,
  curve fitting, minima/maxima detection, and partitioning.
* **Metrology** (`metrology`): dimension and tolerance types built on the geometry primitives.
* **Robust estimation** (`common::consensus`): a MAGSAC++ based consensus framework for fitting primitives to
  data containing outliers.
* **File I/O** (`io`): STL and PLY import/export (feature-gated), the `lptf3` laser profile sensor format, and
  tolerance-based curve, cloud and mesh compression via the [`tol-compress`](https://crates.io/crates/tol-compress)
  crate.
* **Sensor models** (`sensors`): pinhole camera and laser profile sensors.

Types that exist in both two and three dimensions carry a `2` or `3` suffix (`Circle2`/`Circle3`,
`Curve2`/`Curve3`, and so on). `Point2/3`, `Vector2/3`, and `Iso2/3` are aliases over `parry2d-f64`/
`parry3d-f64`, re-exported as `engeom::na`.

## Features

| Feature         | Description                                                    |
|-----------------|----------------------------------------------------------------|
| `stl`           | STL mesh import and export                                     |
| `ply`           | PLY mesh and point cloud import and export                     |
| `three_d`       | An interactive geometry viewer built on `three-d`, for debugging |
| `private_tests` | Gates integration tests that depend on non-public test assets  |

## Python bindings

Python bindings are published to PyPI:

```bash
pip install engeom
```

## Documentation

Combined documentation for the Rust and Python APIs is built with
[mdBook](https://rust-lang.github.io/mdBook/index.html) from the `book/` directory of the repository.

## License

Licensed under either of Apache License, Version 2.0 or MIT license at your option.
