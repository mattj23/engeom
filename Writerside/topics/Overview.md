# Overview

The `engeom` library is a set of tools for engineering geometry applications in Rust, with a focus on metrology, GD&T, and inspection.  While there will be a large amount of overlap with other computational geometry applications, such as computer graphics or CAD, `engeom` was built with a set of foundational principles that favors metrology and fills the gap left by other libraries and approaches.

First, `engeom` prioritizes accuracy/correctness, then speed, and finally memory usage.  Compared to computer graphics this ordering of concerns specifically serves metrology, where accuracy and correctness is absolutely fundamental, speed enables more complex algorithms and larger data sets, and memory problems are usually addressed by purchasing more memory.

Second, when compared to CAD applications, `engeom` is built around primitives mostly focused on real world measured data, leading to a focus on large discrete data sets with noise and uncertainty, and algorithms which operate on these.

## Dependencies

The `engeom` library is built on top of the `parry` libraries, both 2D and 3D, and specifically the `f64` versions of these libraries.

The `parry` libraries in turn are built on top of the `nalgebra` library, which provides the underlying linear algebra and vector/matrix types.  Thus, at its core, `engeom` uses 64 bit `nalgebra` types to build its geometric representations and primitives.


## Goals