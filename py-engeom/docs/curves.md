# Curves

## Introduction

The `engeom` library has both a 2D and 3D curve type. This type represents a 1-dimensional manifold that consists of
a sequence of vertices in $\mathbb{R}^n$ space that are connected by line-segment edges, sometimes referred to as a
polyline.

Each curve entity consists of a single contiguous sequence of vertices with a clear start and end point. In the case of
the 2D `Curve2` type, the curve can also be "closed", meaning that the end point is connected to the start point to
form a closed loop and so operations on the manifold that cross the end point will wrap around to the start point,
and vice versa.

While the fundamental operations on a curve, such as distance queries and manifold traversal are the same or similar
for the 2D and 3D curve types, the 2D curve type has many more features that are made possible by the 2D plane.

## 2D Curves

The 2D curve type, `Curve2`, is a more feature-rich type than the 3D curve type, because the nature of the 2D plane
means that a 2D polyline is more conceptually related to a 3D `Mesh` object than it is to a 3D polyline. A `Curve2`
object can represent a boundary/surface with a clear sense of "inside" and "outside", and if it forms a closed loop
it can model a partitioning of the $\mathbb{R}^2$ plane into an interior and exterior region.

In 2D, a curve has concept of a surface normal direction, which is built from the concept inside/outside defined through
the winding order of the vertices. By convention, the segment from vertex $i$ to vertex $i+1$ defines the space to
its "right" as being outside the curve, and the space to its "left" as being inside the curve, resulting in a
counter-clockwise winding order defining a positive convex shape.

### Creation

Creation of a `Curve2` object is done by passing in a list of ordered vertices that define the curve. The vertices
will be interpreted in sequential order, so that each vertex will be connected to the next vertex in the list.

Adjacent vertices that are within a distance tolerance from each other will be de-duplicated, to prevent the curve from
having zero-length segments. Additionally, if the first and last point are within the distance tolerance, the curve
will be considered "closed" and algorithms which do manifold traversal will wrap between the first and last points.

However, because of the importance of winding order on 2D curves, the `Curve2` type must *also* be constructed with the
vertices in the specific order that matches the intended definition of "inside" vs "outside" for the curve.

While the space occupied by the curve might be exactly the same whether the vertices are constructed in forward or
reverse order, the inside vs outside will be opposite.

There are three ways that winding order can be specified during construction of a `Curve2` object:

1. You can prepare the list of vertices so that they are in the correct order and then pass them into the constructor.
   Counter-clockwise order of vertices for positive convex shapes/regions, and clockwise order for negative convex
   shapes/regions.

2. If you know that the curve is meant to represent a positive (convex) shape overall, you can set the `hull_ccw=True`
   flag. The constructor will build the convex hull of the vertices you provide and reverse their order if the hull
   sequence does not match the input sequence.

3. You may provide the constructor with an array of surface normals that correspond with the vertices. There must be one
   normal per vertex in the input list, and it must be pointing in the direction that you intend to be the "outside" of
   the surface. The constructor will reverse the input order if the majority of normals are not pointing in the same
   direction as the winding order would imply.

The input arguments for the constructor are:

| Argument       | Type                         | Description                                                                                                                                                                                                              |
|----------------|------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `vertices`     | `numpy.ndarray`              | A 2D array of shape `(n, 2)` where `n` is the number of vertices in the curve. The columns are the x and y components of the vertices                                                                                    |
| `normals`      | `numpy.ndarray` **OPTIONAL** | A 2D array of shape `(n, 2)` where `n` is the number of vertices in the curve. The columns are the x and y components of the normals at each vertex. Default is `None`.                                                  |
| `tol`          | `float` **OPTIONAL**         | A tolerance distance, below which points are considered to be the same and de-duplicated. Default is `1e-6`                                                                                                              |
| `force_closed` | `bool` **OPTIONAL**          | If `True`, the curve will be guaranteed to be closed. If the first and last point are more than `tol` distance apart, an additional vertex will be added to the end that overlaps with the beginning. Default is `False` |
| `hull_ccw`     | `bool` **OPTIONAL**          | If `True`, the vertices will be re-ordered to match the convex hull of the vertices. Default is `False`. This will be ignored if `normals` is not `None`                                                                 |

```python
import numpy
from engeom.geom2 import Curve2

# These are the corners of an open unit square
vertices = numpy.array([[0.0, 0], [1, 0], [1, 1], [0, 1]])

# Create a curve with these vertices and nothing else
c1 = Curve2(vertices)
print(c1)  # <Curve2 n=4, l=3 (open)>

# Force the curve to be closed
c2 = Curve2(vertices, force_closed=True)
print(c2)  # <Curve2 n=5, l=4 (closed)>
```

### Stations

Because a `Curve2` object is a 1D manifold, every unique position along the curve can be represented by a single scalar
value, which is the length from the start of the curve. Each unique position along a 2D curve has several geometric
properties which are useful.  These properties are bundled in a `CurveStation2` object, which is a lightweight data 
object that represents a single position on the manifold.

`CurveStation2` objects are not created directly, but are retrieved from a `Curve2` object through one of several different possible queries on the manifold.

The `CurveStation2` object has the following properties:

| Property           | Type            | Description                                                                                                                                      |
|--------------------|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| `.point`           | `Point2`        | The 2D position in space that corresponds with the station on the manifold.                                                                      |
| `.direction`       | `Vector2`       | The vector pointing in the direction of positive distance along the curve. Typically this is the vector from the last vertex to the next vertex. |
| `.normal`          | `Vector2`       | The vector pointing in the direction of the surface normal at the station. This is the `direction` vector rotated by $-90°$.                     |
| `.direction_point` | `SurfacePoint2` | A convenience surface point that combines the `point` position and `direction` vector.                                                           |
| `.surface_point`   | `SurfacePoint2` | A convenience surface point that combines the `point` position and `normal` vector.                                                              |
| `.index`           | `int`           | The index of the previous vertex on the curve, at or before the station.                                                                         | 
| `.length_along`    | `float`         | The distance along the curve from the start to the station. This is the manifold domain.                                                         |


### Querying




