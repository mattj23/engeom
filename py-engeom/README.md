# engeom

Python bindings for [engeom](https://github.com/mattj23/engeom), a Rust library of metrology-focused tools for
working with 2D and 3D geometry. The primary use case is engineering applications such as GD&T and
quality/dimensional inspection, covering every step of the process starting from raw data formats such as meshes
and point clouds.

```bash
pip install engeom
```

## What is included

* **3D geometry** (`engeom.geom3`): points, vectors, surface points, lines, planes, circles, spheres, cylinders,
  cones, triangle meshes, and point clouds, with Levenberg-Marquardt fitting and alignment
  (`engeom.align3`).
* **2D geometry** (`engeom.geom2`): points, vectors, surface points, lines, circles, arcs, and polyline curves,
  with fitting and alignment.
* **Airfoil analysis** (`engeom.airfoil2`): construction and analysis of airfoil cross-sections, including
  camber lines, leading/trailing edge detection, inscribed circles, thickness, and orientation.
* **Raster fields** (`engeom.raster2`, `engeom.raster3`): 2D scalar fields for depth-map style data, with
  filtering, smoothing, in-painting, and region labeling.
* **Metrology** (`engeom.metrology`): dimension and tolerance types built on the geometry primitives.
* **Sensor models** (`engeom.sensors`): pinhole camera and laser profile sensors.
* **Plotting helpers** (`engeom.plot`): convenience adapters for visualizing engeom types.

## Numpy at the boundary

Array-shaped data crosses the Rust/Python boundary as numpy arrays rather than lists of objects. Points and
vectors are shaped `(n, 2)` or `(n, 3)` with dtype `float64`, and face/index arrays use unsigned integer dtypes
such as `uint32`. This matters for both the API shape and for performance.

## Documentation

Combined documentation for the Rust and Python APIs is built with
[mdBook](https://rust-lang.github.io/mdBook/index.html) from the `book/` directory of the
[repository](https://github.com/mattj23/engeom).

## License

Licensed under either of Apache License, Version 2.0 or MIT license at your option.
