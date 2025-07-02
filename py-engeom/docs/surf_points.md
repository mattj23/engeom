# Surface Points

## Overview

Surface points are a composite structure that consist of a point in space and a normal direction. Conceptually, they
come from metrology as a means of representing a point on the surface of an object along with the normal direction of
the surface at that point. However, they are also isomorphic with the concept of a ray or a parameterized line with a
direction of unit length, and can be used in that way as well.

## Common Features

This section defines features of the `SurfacePoint` types that are common to both 2D and 3D versions.

### Construction

Both the 2D and 3D versions need to be created by specifying the cartesian coordinates and the components of the
direction vector. The constructor will automatically scale the direction vector to unit length, so that does *not* need
to be done by the user.

```python
from engeom.geom2 import SurfacePoint2
from engeom.geom3 import SurfacePoint3

sp2 = SurfacePoint2(0, 0, 1, 1)
sp3 = SurfacePoint3(0, 0, 0, 1, 1, 1)
```

For convenience, they can also be built from points and vectors using the unpacking operator and the iterable feature of
the `Point` and `Vector` types.

```python
from engeom.geom2 import SurfacePoint2, Point2, Vector2

p = Point2(0, 0)
v = Vector2(1, 1)

sp = SurfacePoint2(*p, *v)
```

### Accessing the Point and Normal

The `SurfacePoint` types have properties that allow access to the point and normal components, giving you a standard
`Point` and `Vector` object.

```python
from engeom.geom3 import SurfacePoint3

sp = SurfacePoint3(0, 0, 0, 1, 1, 1)

print(sp.point)
# Point3(0, 0, 0)

print(sp.normal)
# Vector3(0.5773502691896258, 0.5773502691896258, 0.5773502691896258)
```

### Projection Distances

There are two distance measuring functions that can be used with `SurfacePoint` types. Both return scalar floating point
values.

* The `scalar_projection` function takes a point in space and returns the scalar projection of the point onto the
  normal vector of the `SurfacePoint`. This is the distance along the normal vector from the `SurfacePoint` to its
  closest approach with the point.

* The `planar_distance` function takes a point in space and returns the distance from the `SurfacePoint` to the test
  point projected into the plane defined by the normal vector of the `SurfacePoint`. This is also the closest distance
  from the test point to the ray defined by the `SurfacePoint`.

![`SurfacePoint` measurements](images/surface_point_meas.svg){width=400}
/// caption
The `SurfacePoint` consists of a point and a unit vector normal. The black point $sp$ is the `SurfacePoint`, and the
black attached arrow is its normal vector. The red point $p$ is another point in the same space. Two simple projection
measurements are show, the `scalar_projection` and the `planar_distance`.  
///

### Projection Points

There are two functions of the `SurfacePoint` types that return new `Point` objects. These are the `at_distance` and
the `projection` functions.

* The `at_distance` function takes a scalar distance and returns a new `Point` object that is that distance along the
  normal vector from the `SurfacePoint`.  This is the same as `sp.point + sp.normal * distance`.

* The `projection` function takes a point in space and returns the point on the `SurfacePoint`'s normal vector that is
  the closest approach to the test point.  This is the same as `sp.at_distance(sp.scalar_projection(p))`.

### Mutating Operations

The `SurfacePoint` types cannot naively participate in addition or subtraction operations because they are neither completely points nor completely vectors, and the results of addition/subtraction are not well-defined.  

However, the following operations are defined:

#### Direction Reverse

This is equivalent to multiplying the normal vector by -1.

```python
from engeom.geom3 import SurfacePoint3
a = SurfacePoint3(0, 0, 0, 1, 0, 0)
b = a.reversed()

print(b)
# SurfacePoint3(0, 0, 0, -1, 0, 0)
```

#### Multiplication/Division by Scalars

Multiplying or dividing a `SurfacePoint` by a scalar will multiply/divide the point's coordinates by the scalar.  The
magnitude of the normal vector will remain the same, but the direction will be inverted if the scalar is negative.

One common use of this is to change the units of the point while keeping the normal vector the same, for instance,
when converting from inches to millimeters, or from millimeters to meters.

```python
from engeom.geom3 import SurfacePoint3
a = SurfacePoint3(1, 2, 3, 1, 0, 0)
b = a * 2

print(b)
# SurfacePoint3(2, 4, 6, 1, 0, 0)
```

## 2D-Only Features

In two dimensions, all rotations take place in the cartesian XY plane and can be described by a single scalar angle. This means there can be a clear concept of orthogonality and rotations are simple to describe.

As a result, there are some additional features that are only available in `SurfacePoint2` objects.

### Normal Rotation

The normal vector of a `SurfacePoint2` can be rotated by a scalar angle. This is done using the `rot_normal` function.

```python
from engeom.geom2 import SurfacePoint2
from math import pi

a = SurfacePoint2(0, 0, 0, 1)
b = a.rot_normal(pi / 2)

    
print(b)
# SurfacePoint2(0, 0, -1, 0)
```

### Orthogonal Shift

This is a convenience feature which shifts the position of the `SurfacePoint2` by a scalar distance in the direction orthogonal to the normal vector. In keeping with the clockwise winding order convention, this is in the direction of the normal vector rotated by 90 degrees in the *clockwise* direction.

For example, if a surface point has a normal vector pointing in the positive Y direction, then a positive shift will move the point in the positive X direction, and a negative shift will move the point in the negative X direction.

```python
from engeom.geom2 import SurfacePoint2
a = SurfacePoint2(0, 0, 0, 1)
b = a.shift_orthogonal(2)

print(b)
# SurfacePoint2(2, 0, 0, 1)
```
