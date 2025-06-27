# Bounding Volumes

There are currently two types of bounding volumes available in the `engeom` library. Both are axis-aligned bounding
boxes (AABBs), but one is for 2D (`Aabb2`) and the other is for 3D (`Aabb3`).

These bounding volumes are primarily used within the Rust language library as an internal mechanism to support data
structures which accelerate distance queries and intersection checks. However, they also have a number of useful
features and so are exposed to the Python API.

## Creating Bounding Volumes

Most commonly, bounding volumes are created internally by the library for geometries that are part of acceleration
structures or have clear spatial bounds, and accessed by retrieving them from those entities. However, through the
Python API, they can also be created directly.

```python
from engeom.geom2 import Aabb2
from engeom.geom3 import Aabb3

# Create a 2D AABB. The arguments are x_min, y_min, x_max, y_max.
box2 = Aabb2(-1, -2, 3, 4)

# Create a 3D AABB. The arguments are x_min, y_min, z_min, x_max, y_max, z_max.
box3 = Aabb3(-1, -2, -3, 3, 4, 5)
```

Two other convenient methods for creating bounding volumes are `from_points` and `at_point`.

The `from_points` function creates a bounding volume that contains all the points in the input list. Pass it a list of
points as a numpy array, and it will determine the bounds.

```python
import numpy
from engeom.geom2 import Aabb2
from engeom.geom3 import Aabb3

points = numpy.array([[0, 0, 0],
                      [1, 1, 1],
                      [2, 2, 2],
                      [3, 3, 3]]).astype(numpy.float64)

# Create boxes that contain all the points
box2 = Aabb2.from_points(points[:, :2])
box3 = Aabb3.from_points(points)
```

The `at_point` function creates a bounding volume centered at the specified point with the specified dimensions. The
arguments are the center point and the dimensions of the box. The 2D version takes 4 arguments (x, y, width, height),
and the 3D version takes 6 arguments (x, y, z, width, height, depth).

```python
from engeom.geom2 import Aabb2
from engeom.geom3 import Aabb3

# Create a 2D AABB centered at the point (1, 2) that is 3 units wide (x) 
# and 4 units tall (y).
box2 = Aabb2.at_point(1, 2, 3, 4)

# Create a 3D AABB centered at the point (1, 2, 3) that is 4 units wide (x),
# 5 units tall (y), and 6 units deep (z).
box3 = Aabb3.at_point(1, 2, 3, 4, 5, 6)
```

## Bounding Volume Properties

There are four properties of bounding volumes, both 2D and 3D, that can be accessed and yield point and vector values
related to the geometry of the volume.

| Property  | Type                | Description                                                      |
|-----------|---------------------|------------------------------------------------------------------|
| `.center` | `Point2`/`Point3`   | The center point of the bounding volume                          |
| `.min`    | `Point2`/`Point3`   | The minimum corner point of the bounding volume                  |
| `.max`    | `Point2`/`Point3`   | The maximum corner point of the bounding volume                  |
| `.extent` | `Vector2`/`Vector3` | The extents of the bounding volume (equivalent to `.max - .min`) |

```python
from engeom.geom3 import Aabb3

box = Aabb3(-1, -2, -3, 1, 2, 3)

# Get the center point of the box
print(box.center)  # Point3(0, 0, 0)
print(box.min)     # Point3(-1, -2, -3)
print(box.max)     # Point3(1, 2, 3)
print(box.extent)  # Vector3(2, 4, 6)
```

## Expand and Shrink

Bounding volumes can be expanded or shrunk by a specified amount, yielding a new bounding volume. The `expand` and
`shrink` methods take a single argument, which is the amount to expand or shrink the perimeter of the bounding volume.

The overall change in the size of the extents will be twice the amount specified. For example, if you expand  a 2D box 
by 1 unit, the width and height will each increase by 2 units.

```python
from engeom.geom2 import Aabb2
box = Aabb2(-1, -2, 1, 2)

print(box.extent)  # Vector2(2, 4)

expanded = box.expand(0.5)
print(expanded.extent)  # Vector2(3, 5)

shrunk = box.shrink(0.5) 
print(shrunk.extent)  # Vector2(1, 3)
```

## Intersection and Containment

!!! warning 
    Intersection and containment options are not yet bound to the Python API. This will be added in a future release.