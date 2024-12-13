# Curve (2D)

The `Curve2` struct and its associated functions and tools are one of the most important 2D constructs in the `engeom`
library.

A `Curve2` is a series of connected points and the line segments between them which may or may not be considered a
closed loop. In the closed case, the first and the last points are the same and the curve functions as a complex
polygon. In the open case, the curve is a series of connected line segments like a polyline.

In addition to having a 2D shape in the 2D plane, a `Curve2` is also a 1D manifold whose domain is its length from start
to end. Every point along that manifold has a corresponding location, direction, and normal in the 2D plane.

## General Features

A general, non-exhaustive list of features of the `Curve2` struct includes:

* Creation from a sequence of `Point2` values, a tolerance distance (below which points are considered to be the same),
  and the ability to force a curve to be closed or allow it to be detected by the location of the first and last points.

    * The option to ensure that on creation the curve is oriented so that its convex hull is clockwise or
      counter-clockwise.

    * The option to create the curve from surface points following the convention that the normal faces outward from a
      counter-clockwise oriented curve.

* The ability to efficiently find a position of interest on a `Curve2` by several different methods. A 'position' along
  the curve includes it's 2D point, it's direction, normal, total length from the start of the curve, and the index of
  the segment it is on.
    * By a total distance from the start of the curve
    * By a fraction of the total length of the curve
    * By the closest point in space to an arbitrary 2D point
    * By the index of a vertex on the curve

* Tools to simplify, refine, and resample the curve.

* Tools to take intersections between the curve and other 2D shapes.

* Tools to trim and subdivide the curve into smaller curves.

* Tools to iterate over the curve in various ways.