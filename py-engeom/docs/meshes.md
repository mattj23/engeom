# Meshes

The `engeom` Python library provides a `Mesh` class which represents an unstructured triangle mesh in 3D space
(simplicial 2-complex).

A `Mesh` consists of a list of vertices and a list of faces. Vertices are a list of 3D points, and faces are a list of
unsigned integer triplets that refer to indices in the list of vertices. Each face has the indices of three vertices
that form a triangle, listed in counter-clockwise order.

!!! warning
    The `engeom` mesh module is still in development and is not yet stable. The API may change in the future.

## Creating a Mesh

There are currently two direct ways to create a new `Mesh` object.

1. From a `numpy` array of vertices and another `numpy` array of faces
2. Load from a STL file

### From `numpy` Arrays

```python
import numpy 
from engeom.geom3 import Mesh

vertices = numpy.array([
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 0],
], dtype=numpy.float64)

triangles = numpy.array([
    [0, 1, 2],
    [3, 4, 5],
], dtype=numpy.uint32)

mesh = Mesh(vertices, triangles)

print(mesh) # <Mesh 6 vertices, 2 faces>
```

!!! tip
    If you get an error `TypeError: argument 'faces': 'ndarray' object cannot be converted to 'PyArray<T, D>'`, make 
    sure to convert the faces array to an unsigned integer type, e.g. `numpy.uint32`.

There are two options which can be used during creation to alter how the vertices and faces are handled.

- `merge_duplicates`: If `True`, duplicate vertices will be merged and duplicate faces will be remoted. Default is
  `False`.
- `delete_degenerate`: If `True`, degenerate faces (faces with two or more vertices that are the same) will be removed. 
  Default is `False`.

```python 
import numpy
from engeom.geom3 import Mesh

vertices = numpy.array([
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 0],
], dtype=numpy.float64)

triangles = numpy.array([
    [0, 1, 2],
    [3, 4, 5],
], dtype=numpy.uint32)

mesh = Mesh(vertices, triangles, merge_duplicates=True)
print(mesh) # <Mesh 4 vertices, 2 faces>
```

### From STL File

To load a mesh from an STL file, use the static `load_stl` method.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")
```

The `load_stl` method also has the `merge_duplicates` and `delete_degenerate` options.

## Access Vertices and Faces

The vertices and faces of a mesh can be accessed using the `vertices` and `faces` properties.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")

print(mesh.vertices) # Numpy ndarray of shape (N, 3), numpy.float64 data type
print(mesh.faces) # Numpy ndarray of shape (M, 3), numpy.uint32 data type
```

## Bounding Volumes

Meshes have a bounding volume which can be accessed using the `aabb` property, yielding an `Aabb3` object.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")
print(mesh.aabb)
```

## Append, Clone, and Transform

The `Mesh` class supports appending the faces and vertices of another mesh.  After the operation, the mesh which calls
the `append` method will have the vertices and faces of the other mesh appended to it, while the other mesh remains
unchanged.

If the `remove_duplicates` flag and/or `delete_degenerate` flag is set to `True`, the vertices and faces of the mesh 
will be cleaned up after the append operation.

```python
from engeom.geom3 import Mesh

mesh1 = Mesh.load_stl("path/to/file1.stl")
mesh2 = Mesh.load_stl("path/to/file2.stl")

mesh1.append(mesh2)
```

The `cloned` method creates a deep copy of the mesh.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")

mesh2 = mesh.cloned()
```

The `transform_by` method applies a transformation matrix to the vertices of the mesh, leaving the faces unchanged. The
object is modified in place.

```python
from engeom.geom3 import Mesh, Iso3

t = Iso3.from_translation(1, 0, 0)
mesh = Mesh.load_stl("path/to/file.stl")

mesh.transform_by(t)
```

## Splitting and Sectioning

There are two operations which can be performed on a mesh using a `Plane3` object: `split` and `section`.

The `split` method takes a `Plane3` object and returns a tuple of two optional `Mesh` objects.

* The first element in the tuple is the part of the mesh which lies in the negative side of the plane, or `None` if
  there is no such part.

* The second element in the tuple is the part of the mesh which lies in the positive side of the plane, or `None` if
  there is no such part.

```python
from engeom.geom3 import Mesh, Plane3

mesh = Mesh.load_stl("path/to/file.stl")
plane = Plane3(1, 0, 0, 1)

mesh1, mesh2 = mesh.split(plane)
```

The `section` method takes a `Plane3` object and returns a list of `Curve3` objects which represent the continuous 
polylines which are the intersection of the mesh with the plane.

```python
from engeom.geom3 import Mesh, Plane3

mesh = Mesh.load_stl("path/to/file.stl")
plane = Plane3(1, 0, 0, 1)

curves = mesh.section(plane)

for curve in curves:
    print(curve)
```

## Splitting Patches

The `split_patches` method splits the mesh into connected components. The method returns a list of `Mesh` objects, each
representing a connected component of the mesh.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")

patches = mesh.separate_patches()
```

## Sampling

The `sample_poisson` method samples the mesh using a Poisson disk sampling algorithm. The method takes a `float`
parameter which represents the minimum distance between points.

```python
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")
result = mesh.sample_poisson(0.1)
```

In the above example, `result` will be a `numpy` array of shape `(N, 6)` where `N` is the number of points that resulted
from the sampling operation, and each row corresponds with an individual point. The first three columns of the 
resulting array are the $x$, $y$, and $z$ coordinates of the point, and the last three columns are the components of the 
normal direction of the surface at that point ($nx$, $ny$, $nz$).

The sampling will be done on the surfaces represented by the faces of the mesh, and the points will be distributed
approximately evenly over the surface. The points will coincide with the vertices of the mesh, but will instead be
random points which lie on the actual triangles.

## Measurements

There are a number of measurements which can be made on a mesh.

### Closest Point on Surface

The closest point on the surface of a `Mesh` to an arbitrary point can be found using the `surface_closest_to` method.
This will yield a `SurfacePoint3` object whose position is the closest point on the mesh, and whose normal is the 
normal of the face on which the point lies.

```python 
from engeom.geom3 import Mesh, Point3

mesh = Mesh.load_stl("path/to/file.stl")

closest = mesh.surface_closest_to(1, 2, 3)

# Don't forget about the unpacking operator if you are working 
# with points or other iterables
p = Point3(1, 2, 3)
cl = mesh.surface_closest_to(*p)
```

### Surface Deviation at a Single Point

Deviation from a surface is a common concept in metrology.  A test point is projected onto the closest face of a mesh,
and the "deviation" is measured as the distance between the test point and its projection.  The distance is signed so
that it is positive if the point lies in the direction of the face normal at the closest point, and negative 
otherwise.

The `measure_point_deviation` method calculates the deviation at the location of a single test point. The method returns
a `Length3` object, which is a metrology entity that represents a scalar distance measured between two positions in 3D 
along a specified direction.

There are two possible modes of computing the distance, specified using the `DeviationMode` enum.  The two modes are
essentially the same except for how they treat points which are beyond the edge of the closest face.

- `DeviationMode.Point`: The deviation is calculated as the direct distance from the test point to the closest point on
  the face.

- `DeviationMode.Plane`: The deviation is calculated as the distance from the test point to the plane of the face on
  which the closest point lies. This allows for points that are slightly beyond the edge of the closest face to have a
  deviation which would be the same as if the edge of the face extended to beyond the test point.

```python
from engeom import DeviationMode
from engeom.geom3 import Mesh

mesh = Mesh.load_stl("path/to/file.stl")

result = mesh.measure_point_deviation(1, 2, 3, DeviationMode.Point)
```

### Surface Deviation at Multiple Points

To calculate the deviation of a large number of points at once, use the `deviation` method, which will return a `numpy`
array with one scalar value of deviation for each test point. Each value will be positive if the test point is in the 
direction of the face normal at the closest point, or negative if it is in the opposite direction.

Additionally, like the `measure_point_deviation` method, there are two possible modes of computing the distance. The 
method is specified using the `DeviationMode` enum.  The two modes are essentially the same except for how they treat
points which are slightly beyond the edge of the closest face.

- `DeviationMode.Point`: The deviation is calculated as the direct distance from the test point to the closest point on
  the face.

- `DeviationMode.Plane`: The deviation is calculated as the distance from the test point to the plane of the face on 
  which the closest point lies. This allows for points that are slightly beyond the edge of the closest face to have a 
  deviation which would be the same as if the edge of the face extended to beyond the test point.


```python
import numpy
from engeom import DeviationMode
from engeom.geom3 import Mesh 

mesh = Mesh.load_stl("path/to/file.stl")

points = numpy.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
], dtype=numpy.float64)

values = mesh.deviation(points, DeviationMode.Plane)
```


## Face Selection and Filtering

Face selection and filtering allow for the selection/de-selection of faces on the mesh based on certain chained 
criteria.  Ultimately, the selection process will produce a list of face indices which can be used for other algorithms,
or to extract a new `Mesh` object constructed from copies of the selected faces.

!!! note
    This documentation is in progress.