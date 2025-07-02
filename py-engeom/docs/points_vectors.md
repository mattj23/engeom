# Points and Vectors

Points and vectors are a fundamental geometric primitive within the `engeom` library that allow for the representation
of positions and directions in 2D and 3D space. They are used as components in a number of other geometric entities, but
can also be used directly.

!!! note
    When working with large numbers of points or vectors in Python, it is much more efficient to use `numpy` arrays than it
    is to use Python lists of `Point2`/`Point3` or `Vector2`/`Vector3` objects. The `engeom` library provides a number of
    functions which can operate directly on `numpy.ndarray`s for clarity and speed.

The 2D point and vector types are located in the `geom2` module, while the 3D point and vector types are located in the
`geom3` module.

```python
from engeom.geom2 import Point2, Vector2
from engeom.geom3 import Point3, Vector3
```

## Creation

Points and vectors both are created by specifying `x`, `y`, and (in the case of the 3D versions) `z` components.

```python
from engeom.geom3 import Point3, Vector3

# Create a 3D point at (1, 2, 3)
p1 = Point3(1, 2, 3)

# Create a 3D vector with components (4, 5, 6)
v1 = Vector3(4, 5, 6)
```

For convenience, you can use Python's `*` unpacking operator to pass in a list or tuple of values.

```python
import numpy
from engeom.geom3 import Point3, Vector3

coords = numpy.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
points = [Point3(*row) for row in coords]

named = {
    "v0": (1, 2, 3),
    "v1": (4, 5, 6),
    "v2": (7, 8, 9)
}

vec0 = Vector3(*named["v0"])
vec1 = Vector3(*named["v1"])
vec2 = Vector3(*named["v2"])
```

Lastly, for convenience, a point's coordinates can be extracted in the form of a vector.

```python
from engeom.geom3 import Point3

p = Point3(1, 2, 3)

v = p.coords
print(v)  # Vector3(1, 2, 3)

```

## Accessing Components

The `x`, `y`, and `z` components of a point or vector can be accessed directly.

```python
from engeom.geom3 import Point3, Vector3

v0 = Vector3(1, 2, 3)

print(v0.x)  # 1
print(v0.y)  # 2
print(v0.z)  # 3
```

The components are also iterable, so you can unpack them back into a list, tuple, or into arguments for a function.

```python
from engeom.geom3 import Point3


def some_function(a: float, b: float, c: float) -> float:
    return a + b + c


p = Point3(1, 2, 3)

x, y, z = p
coords = list(p)
value = some_function(*p)

for component in p:
    print(component)
```

## Scaling Points and Vectors

Points and vectors can be multiplied and divided by scalars.

```python
from engeom.geom3 import Point3, Vector3

p0 = Point3(1, 2, 3)
v0 = Vector3(4, 5, 6)

p1 = p0 * 2
p2 = 3.0 * p0
v1 = v0 / 2
```

Both can also be inverted, which is equivalent to multiplying by `-1`.

```python
from engeom.geom3 import Point3, Vector3

p0 = Point3(1, 2, 3)
v0 = Vector3(4, 5, 6)

p1 = -p0
v1 = -v0
```

## Adding and Subtracting Points and Vectors

Points and vectors can be added to and subtracted from each other, however, the resulting type will vary depending on
the operation.

| Left   | Operand | Right  | Result  |
|--------|---------|--------|---------|
| Point  | `+`     | Vector | Point   |
| Point  | `-`     | Vector | Point   |
| Point  | `+`     | Point  | INVALID |
| Point  | `-`     | Point  | Vector  |
| Vector | `+`     | Vector | Vector  |
| Vector | `-`     | Vector | Vector  |
| Vector | `+`     | Point  | INVALID |
| Vector | `-`     | Point  | INVALID |

To summarize, a point cannot be added to anything, and a point can only be subtracted from another point. A vector can
be added to or subtracted from either a point or another vector.

## Vector Operations

The length of a vector can be calculated using the `norm()` function, and a normalized version of the vector can be
created using the `normalized()` function.

```python
from engeom.geom3 import Vector3

v = Vector3(1, 2, 3)

length = v.norm()
unit = v.normalized()
```

Vectors can be used to perform dot and cross products, and to measure the smallest angle between two vectors.

```python
from engeom.geom3 import Vector3

v0 = Vector3(1, 2, 3)
v1 = Vector3(4, 5, 6)

d = v0.dot(v1)

angle = v0.angle_to(v1)

# The cross product will be a new vector if these are `Vector3` objects, 
# or a scalar if they are `Vector2` objects.
c = v0.cross(v1)
```
