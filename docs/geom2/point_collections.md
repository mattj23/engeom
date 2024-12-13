# Point Collections (2D)

Although 2D point clouds may seem like a less common application than their 3D counterpart, the creation and
manipulation of 2D series of points is a frequent task in many higher level applications. This page describes a series
of tools that are available to make working with `Vec<Point2>` and `&[Point2]` values more convenient.

!!! Note

    For general information on points, vectors, and core spatial data types, see the [Points, Vectors, and Transformations](../common/core_space.md) page.

## Hull Algorithms

The module `engeom::geom2::hull` contains a number of tools for calculating hulls of collections of 2D points. The most
widely known algorithm is the [convex hull](https://en.wikipedia.org/wiki/Convex_hull), which finds the indices of
points which lie on the smallest convex set containing all the input points.

There are also a set of algorithms for finding a hull based on the ball pivot algorithm, which can be used to find a
more generally shaped enclosure (containing both convex and concave regions) around a set of points, effectively looking
for neighbors within a given radius.

### Convex Hull

### Ball Pivot Hull

## K-D Tree (2D)

The module `engeom::geom2::kd_tree2` contains an alias for a
2-dimensional [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) acceleration structure used to efficiently search for
neighbors between points in 2D space. Additionally, it contains a set of helper functions for building and working with
these trees.
