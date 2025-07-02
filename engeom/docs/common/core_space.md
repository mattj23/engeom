# Points, Vectors, and Transformations

The core spatial concepts on which `engeom`'s 2D and 3D geometry tools are built are the representation of coordinates,
directions, and transformations in 2-dimensional and 3-dimensional space.

For that purpose, `engeom` uses convenience types declared on top of the `nalgebra` underlying types for `Point<f64, D>`
and `SVector<f64, D>`, which are used to represent points and vectors in `D`-dimensional space. These are aliased
individually for 2D and 3D versions of these types.

For transformations of points and vectors, `engeom` wraps the `Isometry<T, R, D>` type from `nalgebra` into the
2-dimensional `engeom::Iso2` and the 3-dimensional `engeom::Iso3`.

Algorithms and derived types which are built on these concepts are implemented over generic dimensions wherever
possible. The following sections go into more detail on such generic implementations available for both 2D and 3D
versions of these types.

## Vectors and Unit Vectors

Vectors (directions and distances in space) in `engeom` are represented by aliases on of the statically sized
`SVector<f64, D>` type from `nalgebra`, and so have access to all functions and traits implemented on that type.

For 2D vectors, use the `engeom::Vector2` type, and for 3D vectors use the `engeom::Vector3` type.

Basic vector operations are the same for both 2D and 3D, except for the presence of `.z` in the 3D versions.

```rust
use engeom::Vector2;

fn main() {
    // Creating a new 2D vector
    let v0 = Vector2::new(0.0, 1.0);
    let v1 = Vector2::new(1.0, 0.0);

    // Adding two vectors
    let v2 = v0 + v1;
    println!("v2: {:?}", v2); // Output: v2: [1.0, 1.0]

    // Subtracting two vectors
    let v3 = v0 - v1;
    println!("v3: {:?}", v3); // Output: v3: [-1.0, 1.0]

    // Scaling a vector
    let v4 = v0 * 2.0;
    println!("v4: {:?}", v4); // Output: v4: [0.0, 2.0]

    // Dot product of two vectors
    let dot_product = v0.dot(&v1);
    println!("Dot product: {}", dot_product); // Output: Dot product: 0.0

    // Magnitude (norm) of a vector
    println!("Magnitude: {}, {}", v0.magnitude(), v0.norm()); // Output: Magnitude: 1.0, 1.0

    // Normalizing a vector
    let v5 = v0.normalize();
    println!("v5: {:?}", v5); // Output: v5: [0.0, 1.0]

    // Accessing individual components
    println!("v0.x: {}, v0.y: {}", v0.x, v0.y); // Output: v0.x: 0.0, v0.y: 1.0
}
```

Unit vectors (vectors with a total length of exactly 1.0) are represented by the `engeom::UnitVec2` and
`engeom::UnitVec3` types, which are aliases for the `nalgebra` unit type `Unit<Vector2>` and `Unit<Vector3>`. These are
used to represent unit vectors in 2D and 3D space, which guarantee that the vector has a length of exactly 1.0 are
required for use in certain types and algorithms. A unit vector can be converted to a regular vector by calling the
`into_inner()` function, and created by using the `new_normalize`, `new_unchecked`, or `try_new` functions.

```rust
use engeom::{Vector2, UnitVec2};

fn main() {
    // Creating a new unit vector by normalizing a regular vector
    let v0 = Vector2::new(3.0, 4.0);
    let unit_v0 = UnitVec2::new_normalize(v0);
    println!("unit_v0: {:?}", unit_v0); // Output: unit_v0: [0.6, 0.8]

    // Creating a unit vector directly (unchecked, assumes the input is already normalized)
    let unit_v1 = UnitVec2::new_unchecked(Vector2::new(0.6, 0.8));
    println!("unit_v1: {:?}", unit_v1); // Output: unit_v1: [0.6, 0.8]

    // Attempting to create a unit vector with error checking
    if let Ok(unit_v2) = UnitVec2::try_new(Vector2::new(3.0, 4.0), 1e-6) {
        println!("unit_v2: {:?}", unit_v2); // Output: unit_v2: [0.6, 0.8]
    } else {
        println!("Failed to create unit vector");
    }

    // Converting a unit vector back to a regular vector
    let regular_v = unit_v0.into_inner();
    println!("regular_v: {:?}", regular_v); // Output: regular_v: [0.6, 0.8]

    // Using unit vectors in calculations
    let v1 = Vector2::new(1.0, 0.0);
    let dot_product = unit_v0.dot(&v1);
    println!("Dot product: {}", dot_product); // Output: Dot product: 0.6
}
```

## Points

Points (positions in space) in `engeom` are represented by aliases on of the statically sized `Point<f64, D>` type from
`nalgebra`. They are stored identically to vectors, but their type semantics are built around the concept of a position
in space, rather than a direction or distance.

For 2D points, use the `engeom::Point2` type, and for 3D points use the `engeom::Point3` type.

```rust
use engeom::Point2;
use engeom::common::points::dist;

fn main() {
    // Creating a new 2D point
    let p0 = Point2::new(1.0, 2.0);
    let p1 = Point2::new(3.0, 4.0);

    // Adding a vector to a point
    let v = engeom::Vector2::new(1.0, 1.0);
    let p2 = p0 + v;
    println!("p2: {:?}", p2); // Output: p2: [2.0, 3.0]

    // Subtracting a vector from a point
    let p3 = p1 - v;
    println!("p3: {:?}", p3); // Output: p3: [2.0, 3.0]

    // Subtracting two points to get a vector
    let v2 = p1 - p0;
    println!("v2: {:?}", v2); // Output: v2: [2.0, 2.0]

    // Distance between two points
    let d0 = dist(&p0, &p1);
    let d1 = (p0 - p1).magnitude();
    println!("Distance (dist): {}", d0); // Output: Distance (dist): 2.8284271247461903
    println!("Distance (magnitude): {}", d1); // Output: Distance (magnitude): 2.8284271247461903

    // Accessing individual components
    println!("p0.x: {}, p0.y: {}", p0.x, p0.y); // Output: p0.x: 1.0, p0.y: 2.0

    // Converting a point to a vector
    let v3 = p0.coords;
}
```

The `engeom::common::points` module contains a number of convenience functions for working with points with a generic
number of dimensions. Some of these functions include:

* `dist` - Calculate the distance between two points
* `mid_point` - Calculate the midpoint between two points
* `mean_point` - Calculate the mean point of a list of points
* `mean_point_weighted` - Calculate the weighted mean point of a list of points
* `evenly_spaced_points` - Generate a list of `n` evenly spaced points between two points, including the endpoints
* `evenly_spaced_points_between` - Generate a list of `n` evenly spaced points between two points, excluding the
  endpoints
* `fill_gaps` - Fill in gaps between points in a list of points
* `max_point_in_direction` - Find the point in a list of points that is the farthest in a given direction
* `linear_interpolation_error` - Calculate the error in linearly interpolating between two points rather than including
  a third point between them
* `ramer_dougnlas_peucker` - Simplify a list of points using the Ramer-Douglas-Peucker algorithm
* `transform_points` - Transform a list of points using an isometry

## Surface Points

Joining the concepts of point and vector is the `engeom` specific `SurfacePoint<D>` struct, which represents a point
*and* a corresponding unit vector normal. The aliases `SurfacePoint2` and `SurfacePoint3` are provided for 2D and 3D
geometry. These types are used to represent a point on the surface of a half-space of some sort, or a point that has an
inherent direction. They are also useful as an overdetermined representation for a plane or a ray. There are a number of
convenience functions for dealing with `SurfacePoint<D>` types of different dimensionality.

![`SurfacePoint` measurements](images/surface_point_meas.svg){width=400}
/// caption
The `SurfacePoint` consists of a point and a unit vector normal. The black point $sp$ is the `SurfacePoint`, and the
black attached arrow is its normal vector. The red point $p$ is another point in the same space. Two simple projection
measurements are show, the `scalar_projection` and the `planar_distance`.  
///

The following example demonstrates some of the functions available in both `SurfacePoint2` and `SurfacePoint3`, with
the 2D version used in the example for simplicity.

```rust
use engeom::{Point2, SurfacePoint2, Vector2};
use approx::assert_relative_eq;

fn main() {
    // Creating a new SurfacePoint2 by normalizing a vector
    let sp = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    println!("SurfacePoint: {:?}", sp);

    // Calculating the scalar projection of another point onto the line defined by the SurfacePoint2
    let other = Point2::new(-1.0, -1.0);
    let scalar_projection = sp.scalar_projection(&other);
    println!("Scalar projection: {}", scalar_projection); // Output: Scalar projection: -1.0

    // Finding the point at a certain distance along the normal from the SurfacePoint2
    let point_at_distance = sp.at_distance(2.0);
    println!("Point at distance: {:?}", point_at_distance); // Output: Point at distance: [0.0, 2.0]

    // Projecting another point onto the line defined by the SurfacePoint2
    let projection = sp.projection(&other);
    println!("Projection: {:?}", projection); // Output: Projection: [0.0, -1.0]

    // Reversing the normal of the SurfacePoint2
    let reversed_sp = sp.reversed();
    println!("Reversed SurfacePoint: {:?}", reversed_sp);

    // Calculating the planar distance between the SurfacePoint2 and another point
    let planar_distance = sp.planar_distance(&other);
    println!("Planar distance: {}", planar_distance); // Output: Planar distance: 1.0
}
```

Also, check the documentation for the `engeom::common::surface_point` module, which has additional tools for working
with `SurfacePoint<D>` types of different dimensionality.

## Transformations

Transformations in `engeom` are represented by the `engeom::Iso2` and `engeom::Iso3` types, which are straightforward
aliases for the underlying `Isometry<f64, R, D>` type from `nalgebra`. These types are used to represent a rigid
transformation of points and vectors in 2D and 3D space.

There are many options for creating a new isometry, but be aware that they differ significantly between 2D and 3D due to
the handling of the rotation component. The following examples demonstrates the creation of both types of isometries.

```rust
use engeom::{Iso2, Iso3, Vector2, Vector3};
use std::f64::consts::PI;

fn main() {
    // There is a shortcut for creating an identity isometry
    let identity = Iso2::identity();
    let identity = Iso3::identity();

    // 2D Isometries
    // =============================================================
    // Creating a new 2D isometry is done by providing a translation 
    // vector and a rotation angle
    let a = Iso2::new(Vector2::new(1.0, 2.0), PI / 2.0);
    let a = Iso2::new([1.0, 2.0].into(), PI / 2.0);

    // They can also be created as a translation or rotation only
    let b = Iso2::translation(1.0, 2.0);
    let c = Iso2::rotation(PI / 2.0);

    // 3D Isometries
    // =============================================================
    // Creating a new 2D isometry is done by providing a translation 
    // vector and an axis-angle rotation, in which the angle is 
    // specified by a vector whose direction is the axis of rotation
    // and whose magnitude is the angle of rotation
    let a = Iso3::new(Vector3::new(1.0, 2.0, 3.0),
                      Vector3::new(0.0, 0.0, PI / 2.0));

    let a = Iso3::new([1.0, 2.0, 3.0].into(),
                      [0.0, 0.0, PI / 2.0].into());

    // They can also be created as a translation or rotation only
    let b = crate::Iso3::translation(1.0, 2.0, 3.0);
    let c = crate::Iso3::rotation(Vector3::new(0.0, 0.0, PI / 2.0));
}
```

Transformations can be applied to points, vectors, and other isometries using multiplication, as well as through special
transformation methods. The following example demonstrates some of these capabilities, using the 2D version for
simplicity, although the syntax will be the same for both 2D and 3D.

```rust
use engeom::{Iso2, Point2, Vector2, SurfacePoint2};

fn main() {
    // Create a new 2D isometry that translates by (1.0, 2.0) and rotates by -90 degrees
    let iso = Iso2::new(Vector2::new(1.0, 2.0), -PI / 2.0);


    // When applying an isometry to a point, the point will be both translated and rotated
    // and the result will be a new point
    let p = Point2::new(1.0, 1.0);
    let pt: Point2 = iso * p;

    // When applying an isometry to a vector, the vector will only be rotated, and
    // the result will be a new vector
    let v = Vector2::new(1.0, 1.0);
    let vt: Vector2 = iso * v;

    // When applying an isometry to another isometry, the result will be a new isometry,
    // which is the composition of left isometry applied to the right isometry
    let i = Iso2::rotation(PI / 2.0);
    let it: Iso2 = iso * i;

    // When applying an isometry to a surface point, the point will be transformed 
    // and rotated as a point, and the normal will be rotated as a vector. The result
    // will be a new surface point
    let sp = SurfacePoint2::new_normalize(p, v);
    let spt: SurfacePoint2 = &iso * sp;
}
```
