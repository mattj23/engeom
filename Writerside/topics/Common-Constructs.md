# Common Constructs

The `engeom` library is built on a foundation of common constructs, each of which may apply to many of the different geometric domains of the library.  This section provides a high level overview of those constructs and how they are used.

## Points and Vectors

Points and vectors in `engeom` are mostly convenience type declarations built on top of the `nalgebra` underlying types for `Point<f64>` and `SVector<f64>`, which are used to represent points and vectors in N-dimensional space.  These types are used throughout the library to represent positions, directions, and displacements.  For 2D and 3D geometry, there are the explicitly defined type aliases `Point2` and `Point3`, accessible at the top level of the library crate.

Extending the concept of the vector is the `UnitVec2` and `UnitVec3` types, which are aliases for the `nalgebra` unit type `Unit<Vector2>` and `Unit<Vector3>`.  These are used to represent unit vectors (vectors with a total length of exactly 1.0) in 2D and 3D space.  These types guarantee that the vector is a unit vector, and are required for use in certain cases where the algorithm requires a vector of unit length.  A unit vector can be converted to a regular vector by calling the `into_inner()` function on the `UnitVecN` type.

Finally, joining the concepts of point and vector is the `engeom` specific `SurfacePoint<D>` struct, which represents a point *and* a corresponding unit vector normal.  The aliases `SurfacePoint2` and `SurfacePoint3` are provided for 2D and 3D geometry.  These types are used to represent a point on the surface of a half-space of some sort, or a point that has an inherent direction.  They are also useful as an overdetermined representation for a plane or a ray.  There are a number of convenience functions for dealing with `SurfacePoint<D>` types of different dimensionality.
