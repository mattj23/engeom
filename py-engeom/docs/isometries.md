# Isometries (Rigid-body Transformations)

Isometries are a class of transformations that preserve distances between points. They are also known as rigid-body
transformations, as they preserve the shape of the object being transformed.

In the underlying Rust `engeom` library, isometries are aliases for a native `nalgebra` struct. In two dimensions this
consists of a 2D translation and a unit complex for rotation, and in three dimensions it consists of a 3D translation
and a unit quaternion for rotation.

Isometries can be thought of as equivalent to transformation matrices, but with the limitation that certain matrices
are not valid isometries because they do not preserve distances. They can be composed together by multiplication
according to the rules of matrix multiplication, inverted according to the rules of matrix inversion, and they can be
multiplied against points or vectors to transform them.

## Creating Isometries

Two-dimensional isometries can be created directly from three scalar values: the x and y components of the translation,
and the angle of rotation in radians. Three-dimensional isometries are more complicated, and are easiest to create
through composition.

### Two-dimensional Isometries

```python
from math import pi
from engeom.geom2 import Iso2

# Shortcut for the identity transform
i0 = Iso2.identity()

# Create an isometry that translates by (1, 2) and rotates by pi/4 radians
i1 = Iso2(1, 2, pi / 4)
```

### Three-dimensional Isometries

```python
import numpy
from math import pi
from engeom.geom3 import Iso3

# Shortcut for the identity transform
i0 = Iso3.identity()

# Try to create an isometry directly from a 4x4 transformation matrix. Will throw
# an exception if the matrix is not a valid isometry.
m = numpy.array([[1, 0, 0, 1],
                 [0, 1, 0, 2],
                 [0, 0, 1, 3],
                 [0, 0, 0, 1]])
i1 = Iso3(m)

# Create an isometry that only translates by specifying the translation vector
i2 = Iso3.from_translation(1, 2, 3)

# Create an isometry that rotates by pi/4 radians around the x-axis. See the documentation
# for `from_rotation` for more information on the arguments.
i3 = Iso3.from_rotation(pi / 4, 1, 0, 0)
```

## Inverting Isometries

Isometries can, by definition, be inverted. The inverse of an isometry is an isometry that, when applied to the result
of the original isometry, returns the original input. An isometry multiplied by its inverse is the identity isometry.

```python
from math import pi
from engeom.geom3 import Iso3

i = Iso3.from_rotation(pi / 4, 1, 0, 0)

# Invert the isometry
i_inv = i.inverse()
```


## Composition

Isometries can be composed together by multiplying them together. The order of multiplication is important, as
isometries do not commute. The result of multiplying two isometries together is a new isometry that is equivalent to
apply the right hand isometry first, and then the left hand isometry.

The operator for isometry multiplication is the same as the matrix multiplication operator, `@`.

```python
from math import pi
from engeom.geom3 import Iso3

i1 = Iso3.from_rotation(pi / 4, 1, 0, 0)
i2 = Iso3.from_translation(1, 2, 3)

# Apply the rotation first, then the translation

i3 = i2 @ i1
```

## Transforming Primitives

### Vectors, Points, Surface Points, etc

Isometries can be applied to points, vectors, and other geometric primitives by multiplying them by the isometry using
the `@` operator.

* A point transformed by an isometry is a new point that is at a new position in space.
* A vector transformed by an isometry has been rotated, but its magnitude has not changed.
* A surface point transformed by an isometry is the result of transforming the point and the normal vector
  independently, and then re-constituting them into a new surface point. The position has been moved and the normal
  vector rotated but remains of unit magnitude.

```python
from math import pi
from engeom.geom2 import Iso2, Vector2, Point2, SurfacePoint2

p = Point2(1, 2)
v = Vector2(1, 2)
sp = SurfacePoint2(1, 2, 1, 0)

i = Iso2(1, 2, pi / 4)

p2 = i @ p
v2 = i @ v
sp2 = i @ sp
```

### Numpy Arrays

For efficient transformation of large numbers of points or vectors, both 2D and 3D isometries can be applied to `numpy`
arrays representing points or vectors according to the rules defined in the previous section.

```python
import numpy
from math import pi
from engeom.geom2 import Iso2

values = numpy.array([[1, 2],
                      [3, 4],
                      [5, 6],
                      [7, 8]])

i = Iso2(1, 2, pi / 4)

# Apply the isometry to the values as if they were points
new_points = i.transform_points(values)

# Apply the isometry to the values as if they were vectors
new_vectors = i.transform_vectors(values)
```