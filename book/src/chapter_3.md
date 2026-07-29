# Isometries

An isometry is a transformation that preserves distances between points. Informally, think of it as moving and
rotating a rigid object in space: no matter how you reposition it, the distances between any two points on the object
remain the same. Isometries are also called **rigid-body transformations** for exactly this reason.

In `engeom`, an isometry encodes two things:

- A **rotation**: how much to rotate, and around which axis
- A **translation**: how far to move after rotating

Isometries are equivalent to homogeneous transformation matrices, but with the constraint that the matrix must
represent a valid rigid-body transformation (no shear, no scale). In practice this means they are more numerically
stable than raw matrices because the rotation component is stored as a unit quaternion (3D) or unit complex number
(2D), which is easy to keep normalized.

## Points vs Vectors Under Transformation

A key property of isometries, shared with transformation matrices, is that they treat **points** and **vectors**
differently:

- **Points** are transformed by both the rotation and the translation. A point represents a location in space, so
  moving the coordinate frame moves the point.
- **Vectors** are transformed by rotation only. A vector represents a direction or displacement, not a location, so
  translating the coordinate frame leaves it unchanged.

This matches the behavior of homogeneous coordinates: points have an implicit \\( w = 1 \\) component so translation
applies, while vectors have \\( w = 0 \\) so it does not.

## The Coordinate Frame Interpretation

The most important mental model for isometries in a metrology context is the **coordinate frame** interpretation.
An isometry \\( T \\) describes a local coordinate frame embedded in world space. Applying \\( T \\) to a point answers the
question: *"a point at this location in the local frame: where is it in the world frame?"*

The inverse \\( T^{-1} \\) answers the reverse question: *"a point at this world-space location: where is it in the local
frame?"*

This forward/inverse duality is used constantly in inspection workflows. For example, a datum reference frame might
be described by an isometry \\( T \\). Measured points arrive in world space and need to be expressed in the datum frame
for tolerance evaluation: apply \\( T^{-1} \\). Nominal geometry defined in the datum frame needs to be placed back in world
space for visualization: apply \\( T \\).

## Types

In Rust, `Iso2` and `Iso3` are type aliases for the corresponding nalgebra isometry types:

```rust
use engeom::{Iso2, Iso3};
```

In Python:

```python
from engeom.geom2 import Iso2
from engeom.geom3 import Iso3
```

## Creating Isometries

### Identity

The identity isometry applies no rotation and no translation. It is the starting point for building up
transformations by composition.

**Rust:**

```rust
use engeom::{Iso2, Iso3};

let id2 = Iso2::identity();
let id3 = Iso3::identity();
```

**Python:**

```python
from engeom.geom2 import Iso2
from engeom.geom3 import Iso3

id2 = Iso2.identity()
id3 = Iso3.identity()
```

### 2D Isometries

A 2D isometry is fully specified by an x translation, a y translation, and a rotation angle in radians.

The `from_*` constructors below come from the `IsoExtensions2` trait (`engeom::geom2`), which must be in scope to use
them. `from_translation` and `from_rotation` forward to nalgebra's `Iso2::translation` and `Iso2::rotation`; they exist
so that the whole constructor family follows engeom's `from_<description>` naming convention, and because nalgebra uses
the name `translation` for both a constructor and the isometry's `translation` field.

**Rust:**

```rust
use engeom::geom2::IsoExtensions2;
use engeom::{Iso2, Vector2};
use std::f64::consts::PI;

// Translate by (1, 2) and rotate by π/4 radians
let iso = Iso2::new(Vector2::new(1.0, 2.0), PI / 4.0);

// Translation only
let t = Iso2::from_translation(1.0, 2.0);

// Rotation only
let r = Iso2::from_rotation(PI / 4.0);

// From a 3x3 homogeneous matrix (row-major), from the `IsoExtensions2` trait
let array = [1.0, 0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 1.0];
let m = Iso2::from_array(&array).unwrap();
```

**Python:**

```python
import numpy as np
from math import pi
from engeom.geom2 import Iso2

# Translate by (1, 2) and rotate by π/4 radians
iso = Iso2(1, 2, pi / 4)

# Identity
i0 = Iso2.identity()

# From a 3x3 homogeneous matrix (raises an exception if not a valid isometry)
m = np.array([[1, 0, 1],
              [0, 1, 2],
              [0, 0, 1]], dtype=np.float64)
i1 = Iso2.from_array(m)
```

### 2D Isometries from an Arbitrary Rotation Point

`Iso2::from_rotation` always rotates around the origin. To rotate around an arbitrary point, such as the center of a
fitted circle, use `from_rotation_about`, provided by the `IsoExtensions2` trait (`engeom::geom2`). In 2D the axis of
rotation is always perpendicular to the plane, so unlike the 3D `from_rot_axis` (which needs a `Line3` to describe an
arbitrary direction), it is fully specified by a single point.

**Rust:**

```rust
use engeom::geom2::IsoExtensions2;
use engeom::{Iso2, Point2};
use std::f64::consts::PI;

let center = Point2::new(1.0, 1.0);
let iso = Iso2::from_rotation_about(&center, PI / 2.0);
```

**Python:**

```python
from math import pi
from engeom.geom2 import Iso2, Point2

center = Point2(1, 1)
iso = Iso2.from_rotation_about(center, pi / 2)
```

### 2D Isometries from a Basis Vector

A 2D isometry's rotation is a single degree of freedom, so unlike the 3D `from_basis_xy`/`from_basis_xz`/etc. family
(which needs two vectors to pin down all three axes), a single vector is enough to fully determine the frame: the
second axis is always a fixed 90 degree turn away from the first. `from_basis_x` treats the vector as the local
x-axis; `from_basis_y` treats it as the local y-axis.

**Rust:**

```rust
use engeom::geom2::IsoExtensions2;
use engeom::{Iso2, Point2, Vector2};

let e0 = Vector2::new(1.0, 1.0); // will become the X axis
let origin = Some(Point2::new(1.0, 2.0));

let frame = Iso2::from_basis_x(&e0, origin).unwrap();
```

**Python:**

```python
from engeom.geom2 import Iso2, Vector2, Point2

e0 = Vector2(1, 1)  # will become the X axis
origin = Point2(1, 2)

frame = Iso2.from_basis_x(e0, origin=origin)
```

### 3D Isometries

3D isometries are more involved because rotations in three dimensions require three degrees of freedom. Several
construction paths are available.

**Rust:**

Rotations are encoded as an axis-angle vector, whose direction is the rotation axis and whose magnitude is the
rotation angle in radians. `Iso3::new` (from nalgebra) takes a translation vector and an axis-angle vector together,
while `from_translation` and `from_rotation` (from the `IsoExtensions3` trait) each build one component on its own.
The latter two forward to nalgebra's `Iso3::translation` and `Iso3::rotation`, under names that follow engeom's
`from_<description>` convention.

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

// Translate by (1, 2, 3) with no rotation
let t = Iso3::from_translation(1.0, 2.0, 3.0);

// Rotate by π/4 around the x-axis (axis-angle encoding)
let r = Iso3::from_rotation(&(Vector3::x() * (PI / 4.0)));

// Combine translation and rotation directly
let iso = Iso3::new(Vector3::new(1.0, 2.0, 3.0), Vector3::x() * (PI / 4.0));
```

The `IsoExtensions3` trait (from `engeom::geom3`) provides convenience constructors for the common case of a pure
rotation around a principal axis:

```rust
use engeom::geom3::IsoExtensions3;
use engeom::Iso3;
use std::f64::consts::PI;

let rx = Iso3::from_rx(PI / 4.0); // rotate around X
let ry = Iso3::from_ry(PI / 4.0); // rotate around Y
let rz = Iso3::from_rz(PI / 4.0); // rotate around Z
```

**Python:**

```python
import numpy as np
from math import pi
from engeom.geom3 import Iso3

# Identity
i0 = Iso3.identity()

# Translation only
i1 = Iso3.from_translation(1, 2, 3)

# Rotation only: angle in radians, then the axis components
i2 = Iso3.from_rotation(pi / 4, 1, 0, 0)  # rotate π/4 around X

# From a 4x4 homogeneous matrix (raises an exception if not a valid isometry)
m = np.array([[1, 0, 0, 1],
              [0, 1, 0, 2],
              [0, 0, 1, 3],
              [0, 0, 0, 1]], dtype=np.float64)
i3 = Iso3(m)
```

### 3D Isometries from an Arbitrary Rotation Axis

`from_rx`, `from_ry`, and `from_rz` only rotate around the principal axes as they pass through the origin. To rotate
around an arbitrary axis, such as one derived from a fitted cylinder or a measured hinge line, use `from_rot_axis`,
which takes a `Line3` and an angle. The line's direction does not need to be normalized, and the axis does not need
to pass through the origin.

**Rust:**

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Line3, Point3, Vector3};
use std::f64::consts::PI;

let axis = Line3::new(Point3::new(1.0, 1.0, 0.0), Vector3::z());
let iso = Iso3::from_rot_axis(&axis, PI / 2.0).unwrap();
```

**Python:**

```python
from math import pi
from engeom.geom3 import Iso3, Line3

# A line through (1, 1, 0) pointing along the z-axis
axis = Line3(1, 1, 0, 0, 0, 1)
iso = Iso3.from_rot_axis(axis, pi / 2)
```

### 3D Isometries from Basis Vectors

A common need in metrology is to construct a coordinate frame from two measured directions: for example, the normal
of a fitted plane and the direction of a measured edge. A family of methods exists for this purpose.

Each method takes two vectors that define two of the three axes, plus an optional origin point. The third axis is
computed by cross product to form a right-handed coordinate system. The naming convention is `from_basis_xy` where
the first letter is the primary axis (set directly from the first argument after normalizing) and the second letter
is the secondary axis (orthogonalized against the primary). All variants fail if the two input vectors are parallel
(i.e. the cross product is near zero).

The full set of variants is: `from_basis_xy`, `from_basis_xz`, `from_basis_yx`, `from_basis_yz`, `from_basis_zx`,
`from_basis_zy`.

**Rust:**

In Rust these are provided by the `IsoExtensions3` trait and return a `Result`.

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Point3, Vector3};

// Build a frame where e0 becomes X and e1 (after orthogonalization) becomes Y
let e0 = Vector3::new(1.0, 1.0, 0.0); // will become the X axis
let e1 = Vector3::new(0.0, 1.0, 0.0); // Y axis derived from this
let origin = Some(Point3::new(1.0, 2.0, 3.0));

let frame = Iso3::from_basis_xy(&e0, &e1, origin).unwrap();
```

**Python:**

In Python the methods are static methods on `Iso3`. They raise a `ValueError` if the
vectors are parallel. The `origin` argument is an optional `Point3`.

```python
from engeom.geom3 import Iso3, Vector3, Point3

e0 = Vector3(1, 1, 0)  # will become the X axis
e1 = Vector3(0, 1, 0)  # Y axis derived from this
origin = Point3(1, 2, 3)

frame = Iso3.from_basis_xy(e0, e1, origin=origin)

# origin is optional; omit it for a frame centred at the world origin
frame_no_origin = Iso3.from_basis_xy(e0, e1)
```

## Inverting Isometries

The inverse of an isometry undoes the transformation. If \\( T \\) maps points from local space into world space, then
\\( T^{-1} \\) maps points from world space back into local space. An isometry multiplied by its inverse is the identity.

**Rust:**

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

let iso = Iso3::from_rotation(&(Vector3::x() * (PI / 4.0)));
let inv = iso.inverse();
```

**Python:**

```python
from math import pi
from engeom.geom3 import Iso3

iso = Iso3.from_rotation(pi / 4, 1, 0, 0)
inv = iso.inverse()
```

## Flipping an Isometry 180 Degrees

Sometimes a coordinate frame needs to be turned around without moving its origin, for example when a measured
feature's direction convention is the opposite of what a downstream calculation expects. This is a 180-degree
**rotation**, not a mirror/reflection: the origin stays put, and the axes that reverse do so in pairs so that the
result is still a proper rigid-body transformation.

In 3D, `IsoExtensions3` provides `flipped_around_x`, `flipped_around_y`, and `flipped_around_z`, each of which keeps
one axis fixed and reverses the other two.

In 2D there is only one such operation: `IsoExtensions2::flipped` keeps the origin fixed and reverses both the x-axis
and y-axis together, since the only axis a proper 180-degree rotation can turn around in the plane is the implicit
axis perpendicular to it. There is no 2D equivalent of `flipped_around_x`/`flipped_around_y`, because reversing only
one in-plane axis is a reflection (determinant -1), which `Iso2` cannot represent at all.

**Rust:**

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

let iso3 = Iso3::from_rotation(&(Vector3::x() * (PI / 4.0)));
let flipped3 = iso3.flipped_around_z();
```

```rust
use engeom::geom2::IsoExtensions2;
use engeom::Iso2;
use std::f64::consts::PI;

let iso2 = Iso2::from_rotation(PI / 4.0);
let flipped2 = iso2.flipped();
```

**Python:**

```python
from math import pi
from engeom.geom3 import Iso3

iso3 = Iso3.from_rotation(pi / 4, 1, 0, 0)
flipped3 = iso3.flipped_around_z()
```

```python
from math import pi
from engeom.geom2 import Iso2

iso2 = Iso2(0, 0, pi / 4)
flipped2 = iso2.flipped()
```

## Composing Isometries

Isometries can be composed by multiplication. The result is a single isometry equivalent to applying the right-hand
operand first and then the left-hand operand: the same convention as matrix multiplication.

Importantly, isometry multiplication is **not commutative**: rotating then translating produces a different result
than translating then rotating.

**Rust:**

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

let r = Iso3::from_rotation(&(Vector3::x() * (PI / 4.0)));
let t = Iso3::from_translation(1.0, 2.0, 3.0);

// Rotate first, then translate
let rt = t * r;

// Translate first, then rotate
let tr = r * t;
```

**Python:**

Python uses the `@` operator for isometry composition (matching numpy's matrix-multiply convention):

```python
from math import pi
from engeom.geom3 import Iso3

r = Iso3.from_rotation(pi / 4, 1, 0, 0)
t = Iso3.from_translation(1, 2, 3)

# Rotate first, then translate
rt = t @ r

# Translate first, then rotate
tr = r @ t
```

## Applying Isometries to Primitives

### Points, Vectors, and Surface Points

An isometry is applied to a geometric primitive by multiplication. In Rust the operator is `*`; in Python it is `@`.

Recall that points are affected by both rotation and translation, while vectors are only rotated.

**Rust:**

```rust
use engeom::{Iso2, Point2, SurfacePoint2, Vector2};
use std::f64::consts::PI;

let iso = Iso2::new(Vector2::new(1.0, 2.0), PI / 4.0);

let p  = Point2::new(1.0, 0.0);
let v  = Vector2::new(1.0, 0.0);
let sp = SurfacePoint2::new_normalize(Point2::new(1.0, 0.0), Vector2::new(1.0, 0.0));

let p2  = iso * p;   // rotated and translated
let v2  = iso * v;   // rotated only
let sp2 = iso * sp;  // point rotated+translated, normal rotated only
```

**Python:**

```python
from math import pi
from engeom.geom2 import Iso2, Point2, Vector2, SurfacePoint2

iso = Iso2(1, 2, pi / 4)

p  = Point2(1, 0)
v  = Vector2(1, 0)
sp = SurfacePoint2(1, 0, 1, 0)

p2  = iso @ p
v2  = iso @ v
sp2 = iso @ sp
```

### NumPy Arrays (Python only)

For large batches of points or vectors, isometries can be applied directly to `(n, 2)` or `(n, 3)` numpy arrays
using `transform_points` and `transform_vectors`. This is significantly faster than transforming individual Python
objects in a loop.

```python
import numpy as np
from math import pi
from engeom.geom2 import Iso2

iso = Iso2(1, 2, pi / 4)

values = np.array([[1.0, 0.0],
                   [2.0, 1.0],
                   [3.0, 2.0]])

# Transform as points (rotation + translation)
new_points = iso.transform_points(values)

# Transform as vectors (rotation only)
new_vectors = iso.transform_vectors(values)
```
