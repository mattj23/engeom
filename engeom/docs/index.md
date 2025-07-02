# Overview

The `engeom` library is a set of tools for engineering geometry applications in Rust (with a reduced set of bindings in Python), with a focus on metrology, GD&T, and inspection.  

While there is significant overlap with other computational geometry applications (computer graphics, CAD, etc.), `engeom` was built on a set of base principles that favors metrology and fills the gap left by other libraries and approaches.

The `engeom` library puts its priorities in the following order:

1. Accuracy and correctness
2. Speed
3. Memory usage

Metrology differs from computer graphics in that accuracy and correctness are absolutely fundamental, whereas graphics applications can usually get away with only being visually correct. Speed is important insofar as it enables more complex algorithms and larger data sets.  Memory limitations, on the other hand, are usually addressed by purchasing more memory.

When compared to CAD applications, `engeom` is built around primitives mostly focused on real world measured data, leading to a focus on large discrete data sets with noise and uncertainty, and algorithms which operate on these.

## Dependencies

The `engeom` library is built on top of the `parry` libraries, both 2D and 3D, and specifically the `f64` versions of these libraries.

The `parry` libraries in turn are built on top of the `nalgebra` library, which provides the underlying linear algebra and vector/matrix types.  Thus, at its core, `engeom` uses 64 bit `nalgebra` types to build its geometric representations and primitives.


## Goals

The `engeom` library is still in its early stages, but a summary of the long-term goals for the project include:

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