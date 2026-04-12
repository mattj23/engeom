# 2D and 3D Primitives

## Points and Vectors

Points and vectors are a fundamental geometric primitive within the `engeom` library that allow for the representation
of positions and directions in 2D and 3D space. They are used as components in a number of other geometric entities,
but can also be used directly.

In Rust, the point and vector types are type aliases for [nalgebra](https://nalgebra.org/) types, re-exported directly
from the `engeom` crate:

```rust
use engeom::{Point2, Vector2, UnitVec2};
use engeom::{Point3, Vector3, UnitVec3};
```

In Python, the 2D point and vector types are located in the `geom2` module, while the 3D types are in the `geom3`
module:

```python
from engeom.geom2 import Point2, Vector2
from engeom.geom3 import Point3, Vector3
```

> **Note (Python):** When working with large numbers of points or vectors in Python, it is much more efficient to use
> `numpy` arrays than lists of `Point2`/`Point3` or `Vector2`/`Vector3` objects. The `engeom` library provides a
> number of functions which can operate directly on `numpy.ndarray`s for clarity and speed.

### Creation

Points and vectors are created by specifying `x`, `y`, and (for 3D) `z` components.

**Rust:**

```rust
use engeom::{Point3, Vector3};

// Create a 3D point at (1, 2, 3)
let p1 = Point3::new(1.0, 2.0, 3.0);

// Create a 3D vector with components (4, 5, 6)
let v1 = Vector3::new(4.0, 5.0, 6.0);

// 2D equivalents
use engeom::{Point2, Vector2};
let p2 = Point2::new(1.0, 2.0);
let v2 = Vector2::new(3.0, 4.0);
```

You can also convert an array or slice into a point or vector using `from`:

```rust
use engeom::{Point3, Vector3};

let p = Point3::from([1.0, 2.0, 3.0]);
let v = Vector3::from([4.0, 5.0, 6.0]);
```

**Python:**

```python
from engeom.geom3 import Point3, Vector3

# Create a 3D point at (1, 2, 3)
p1 = Point3(1, 2, 3)

# Create a 3D vector with components (4, 5, 6)
v1 = Vector3(4, 5, 6)
```

Python's `*` unpacking operator can be used to pass in a list or tuple of values:

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
```

### Accessing Components

The `x`, `y`, and `z` components of a point or vector can be accessed directly in both languages.

**Rust:**

```rust
use engeom::Vector3;

let v0 = Vector3::new(1.0, 2.0, 3.0);

println!("{}", v0.x);  // 1
println!("{}", v0.y);  // 2
println!("{}", v0.z);  // 3
```

**Python:**

```python
from engeom.geom3 import Vector3

v0 = Vector3(1, 2, 3)

print(v0.x)  # 1
print(v0.y)  # 2
print(v0.z)  # 3
```

In Python, the components are also iterable, so you can unpack them back into a list, tuple, or function arguments:

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

### Coords

A point's position can be extracted in the form of a vector using the `coords` property.

**Rust:**

```rust
use engeom::Point3;

let p = Point3::new(1.0, 2.0, 3.0);
let v = p.coords;  // Vector3(1, 2, 3)
```

**Python:**

```python
from engeom.geom3 import Point3

p = Point3(1, 2, 3)
v = p.coords
print(v)  # Vector3(1, 2, 3)
```

### Scaling Points and Vectors

Points and vectors can be multiplied and divided by scalars.

**Rust:**

```rust
use engeom::{Point3, Vector3};

let p0 = Point3::new(1.0, 2.0, 3.0);
let v0 = Vector3::new(4.0, 5.0, 6.0);

let p1 = p0 * 2.0;
let v1 = v0 / 2.0;
```

Both can also be negated, which is equivalent to multiplying by `-1`:

```rust
use engeom::{Point3, Vector3};

let p0 = Point3::new(1.0, 2.0, 3.0);
let v0 = Vector3::new(4.0, 5.0, 6.0);

let p1 = -p0;
let v1 = -v0;
```

**Python:**

```python
from engeom.geom3 import Point3, Vector3

p0 = Point3(1, 2, 3)
v0 = Vector3(4, 5, 6)

p1 = p0 * 2
p2 = 3.0 * p0
v1 = v0 / 2

p3 = -p0
v2 = -v0
```

### Adding and Subtracting Points and Vectors

Points and vectors can be added to and subtracted from each other; however, the resulting type depends on the
operation:

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

To summarize: a point cannot be added to anything, and a point can only be subtracted from another point (yielding
the vector between them). A vector can be added to or subtracted from either a point or another vector.

**Rust:**

```rust
use engeom::{Point3, Vector3};

let p = Point3::new(1.0, 2.0, 3.0);
let v = Vector3::new(1.0, 0.0, 0.0);

// Translate a point by a vector
let p2 = p + v;

// Get the vector between two points
let p3 = Point3::new(4.0, 5.0, 6.0);
let diff: Vector3 = p3 - p;

// Add two vectors
let v2 = v + diff;
```

**Python:**

```python
from engeom.geom3 import Point3, Vector3

p = Point3(1, 2, 3)
v = Vector3(1, 0, 0)

# Translate a point by a vector
p2 = p + v

# Get the vector between two points
p3 = Point3(4, 5, 6)
diff = p3 - p  # Vector3

# Add two vectors
v2 = v + diff
```

### Vector Operations

The length of a vector can be calculated, and a normalized (unit-length) version can be obtained.

**Rust:**

```rust
use engeom::Vector3;

let v = Vector3::new(1.0, 2.0, 3.0);

let length = v.norm();
let unit = v.normalize();
```

For cases where you need to guarantee unit length at the type level, `UnitVec3` (or `UnitVec2`) wraps a vector and
enforces that it is always normalized:

```rust
use engeom::{UnitVec3, Vector3};

let unit = UnitVec3::new_normalize(Vector3::new(1.0, 2.0, 3.0));
let inner: &Vector3 = unit.as_ref();
```

Vectors support dot products, cross products, and measuring the smallest angle between two vectors:

```rust
use engeom::Vector3;

let v0 = Vector3::new(1.0, 2.0, 3.0);
let v1 = Vector3::new(4.0, 5.0, 6.0);

let d = v0.dot(&v1);

// angle_between returns the angle in radians
let angle = v0.angle(&v1);

// cross product produces a new Vector3 (or a scalar for Vector2)
let c = v0.cross(&v1);
```

**Python:**

```python
from engeom.geom3 import Vector3

v = Vector3(1, 2, 3)

length = v.norm()
unit = v.normalized()

v0 = Vector3(1, 2, 3)
v1 = Vector3(4, 5, 6)

d = v0.dot(v1)

angle = v0.angle_to(v1)

# The cross product will be a new vector if these are `Vector3` objects,
# or a scalar if they are `Vector2` objects.
c = v0.cross(v1)
```

## Surface Points

A surface point is a composite structure consisting of a point in space and a unit-length normal direction. They originate from metrology as a way to represent a point on the surface of an object together with the surface normal at that location. However, they are also isomorphic with a ray or a parameterized line with a unit direction vector, and can be used in that way as well.

In Rust the types are `SurfacePoint2` and `SurfacePoint3`, generic over dimension internally but exposed as concrete aliases:

```rust
use engeom::{SurfacePoint2, SurfacePoint3};
```

In Python they live in the respective geometry modules:

```python
from engeom.geom2 import SurfacePoint2
from engeom.geom3 import SurfacePoint3
```

### Construction

Both types are constructed by providing the cartesian coordinates of the point followed by the components of the normal vector. The normal is automatically normalized, so it does not need to be a unit vector beforehand.

**Rust:**

```rust
use engeom::{Point2, Point3, SurfacePoint2, SurfacePoint3, Vector2, Vector3};

// 2D: point at origin, normal pointing in the +Y direction
let sp2 = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));

// 3D: point at origin, normal pointing along (1, 1, 1)
let sp3 = SurfacePoint3::new_normalize(Point3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
```

If you already have a `UnitVec2` or `UnitVec3` you can use `new` to skip the normalization step:

```rust
use engeom::{Point3, SurfacePoint3, UnitVec3, Vector3};

let normal = UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0));
let sp = SurfacePoint3::new(Point3::new(1.0, 2.0, 3.0), normal);
```

**Python:**

```python
from engeom.geom2 import SurfacePoint2
from engeom.geom3 import SurfacePoint3

sp2 = SurfacePoint2(0, 0, 1, 1)      # x, y, nx, ny
sp3 = SurfacePoint3(0, 0, 0, 1, 1, 1)  # x, y, z, nx, ny, nz
```

The `*` unpacking operator can be combined with `Point` and `Vector` iterability for convenience:

```python
from engeom.geom2 import SurfacePoint2, Point2, Vector2

p = Point2(0, 0)
v = Vector2(1, 1)

sp = SurfacePoint2(*p, *v)
```

### Accessing the Point and Normal

The point and normal fields are directly accessible as public members.

**Rust:**

```rust
use engeom::{Point3, SurfacePoint3, Vector3};

let sp = SurfacePoint3::new_normalize(Point3::new(1.0, 2.0, 3.0), Vector3::new(0.0, 0.0, 1.0));

let p = sp.point;           // Point3
let n = sp.normal;          // UnitVec3
let n_vec = n.into_inner(); // Vector3
```

**Python:**

```python
from engeom.geom3 import SurfacePoint3

sp = SurfacePoint3(0, 0, 0, 1, 1, 1)

print(sp.point)   # Point3(0, 0, 0)
print(sp.normal)  # Vector3(0.577..., 0.577..., 0.577...)
```

### Projection Distances

Two scalar distance measurements are defined on surface points:

- `scalar_projection`: the signed distance from the surface point to a test point, measured along the normal direction. Positive means the test point is on the same side as the normal; negative means the opposite side.
- `planar_distance`: the perpendicular distance from the test point to the line/ray defined by the surface point. This is always non-negative.

**Rust:**

```rust
use engeom::{Point3, SurfacePoint3, Vector3};

let sp = SurfacePoint3::new_normalize(Point3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0));
let p  = Point3::new(3.0, 4.0, 5.0);

let along  = sp.scalar_projection(&p);  // 5.0  (distance along normal)
let lateral = sp.planar_distance(&p);   // 5.0  (distance in XY plane, sqrt(3²+4²))
```

**Python:**

```python
from engeom.geom3 import SurfacePoint3

sp = SurfacePoint3(0, 0, 0, 0, 0, 1)
p  = (3, 4, 5)  # or a Point3

along   = sp.scalar_projection(p)
lateral = sp.planar_distance(p)
```

### Projection Points

Two functions return new `Point` objects derived from projecting a test point onto the surface point's line:

- `at_distance(d)`: returns the point that is distance `d` along the normal from the surface point's position. Equivalent to `sp.point + sp.normal * d`.
- `projection(p)`: returns the closest point on the surface point's line to a test point. Equivalent to `sp.at_distance(sp.scalar_projection(p))`.

**Rust:**

```rust
use engeom::{Point3, SurfacePoint3, Vector3};

let sp = SurfacePoint3::new_normalize(Point3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0));

let ahead = sp.at_distance(3.0);       // Point3(0, 3, 0)
let proj  = sp.projection(&Point3::new(5.0, 2.0, 0.0)); // Point3(0, 2, 0)
```

**Python:**

```python
from engeom.geom3 import SurfacePoint3

sp = SurfacePoint3(0, 0, 0, 0, 1, 0)

ahead = sp.at_distance(3.0)
proj  = sp.projection((5, 2, 0))
```

### Mutating Operations

Surface points cannot participate in plain addition or subtraction because the semantics are ambiguous. Instead, several dedicated operations are provided.

#### Reversing the Normal

Returns a new surface point at the same position with the normal vector flipped.

**Rust:**

```rust
use engeom::{Point3, SurfacePoint3, Vector3};

let sp  = SurfacePoint3::new_normalize(Point3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
let rev = sp.reversed();
```

**Python:**

```python
from engeom.geom3 import SurfacePoint3

a = SurfacePoint3(0, 0, 0, 1, 0, 0)
b = a.reversed()
# SurfacePoint3(0, 0, 0, -1, 0, 0)
```

#### Shifting Along the Normal

Returns a new surface point displaced from the original along the normal direction, keeping the same normal.

**Rust:**

```rust
use engeom::{Point3, SurfacePoint3, Vector3};

let sp      = SurfacePoint3::new_normalize(Point3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0));
let shifted = sp.new_shifted(2.0); // point moves to (0, 2, 0), normal unchanged
```

#### Scaling (Python only)

In Python, multiplying or dividing a surface point by a scalar scales the point coordinates while keeping the normal direction unchanged. If the scalar is negative the normal is also flipped. A common use case is unit conversion.

```python
from engeom.geom3 import SurfacePoint3

a = SurfacePoint3(1, 2, 3, 1, 0, 0)
b = a * 2
# SurfacePoint3(2, 4, 6, 1, 0, 0)
```

## 2D-Only Surface Point Features

In two dimensions all rotations are described by a single scalar angle, which makes it practical to define additional operations on `SurfacePoint2`.

### Normal Rotation

The normal vector can be rotated by an arbitrary angle (in radians).

**Rust:**

```rust
use std::f64::consts::PI;
use engeom::{Point2, SurfacePoint2, Vector2};

let sp      = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
let rotated = sp.rot_normal(PI / 2.0);
// Normal now points in the -X direction
```

**Python:**

```python
from engeom.geom2 import SurfacePoint2
from math import pi

a = SurfacePoint2(0, 0, 0, 1)
b = a.rot_normal(pi / 2)
# SurfacePoint2(0, 0, -1, 0)
```

### Orthogonal Shift

Shifts the position of the surface point in the direction orthogonal to its normal. Following the clockwise winding order convention, this direction is the normal rotated 90 degrees clockwise. For a normal pointing in +Y, a positive shift moves the point in the +X direction.

**Rust:**

```rust
use engeom::{Point2, SurfacePoint2, Vector2};

let sp      = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
let shifted = sp.shift_orthogonal(2.0);
// point moves to (2, 0), normal unchanged
```

**Python:**

```python
from engeom.geom2 import SurfacePoint2

a = SurfacePoint2(0, 0, 0, 1)
b = a.shift_orthogonal(2)
# SurfacePoint2(2, 0, 0, 1)
```
