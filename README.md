# Engineering Geometry Library for Rust

The goal of this library is to provide a comprehensive set of metrology-focused tools for working with 2D and 3D geometry in Rust.  The primary use case is engineering applications such as GD&T and quality/dimensional inspection, and the library contains a wide set of tools to help work with every step of the process starting with raw data formats such as meshes and point clouds.

## Documentation

The documentation for this project is built with mkdocs-material.  To set up the documentation locally, use a virtual environment and run `mkdocs serve`.

```bash
pip install mkdocs-material

# Run mkdocs serve from the folder containing the mkdocs.yml file
mkdocs serve
```

##  General Principles

Because this is a metrology-focused library, the order of priority for algorithms and methods is:

1. Accuracy/correctness
2. Speed
3. Memory usage

Ultimately, the goal is that this library will contain support for:

* 3D Geometry
    * Measurements on points, point clouds, and unstructured meshes
    * Construction of geometric primitives such as surface points, lines, spheres, planes, and more
    * Levenberg-Marquardt fitting and alignment
    * Measurement of distances, angles, etc

* 2D Geometry:
    * Measurements on points and polylines
    * Construction of geometric primitives such as surface points, lines, circles
    * Levenberg-Marquardt fitting and alignment
    * Measurement of distances, angles, etc
    * Special tools for construction and analysis of airfoil cross-sections

* 2D Raster Fields:
    * Typically for applications like depth maps
    * Basic operations for binning, filtering, smoothing, in-painting, and other tools based on image processing

* 1D Scalar Series
    * Typical applications are for spatial series sampled off of 2D or 3D surfaces, or for time series sampled for motion
    * Represent series as a function over a single variable domain
    * Allow for operations such as interpolation, smoothing, filtering, minima/maxima detection, curve fitting, partitioning, etc

* Transformations between domains:
    * 3D to 2D projections
    * 3D mesh topology to flattened 2D topology
    * Transformation of 3D deviations to 2D raster fields
    * Sampling of 3D or 2D data to 1D scalar series
    * Projections of 2D data to 3D surfaces
    * Projections of 1D data to 2D or 3D points, lines, or other primitives






