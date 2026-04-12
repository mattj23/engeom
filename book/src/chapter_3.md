# Isometries

An isometry is a transformation that preserves distances between points. Informally, think of it as moving and
rotating a rigid object in space: no matter how you reposition it, the distances between any two points on the object
remain the same. Isometries are also called **rigid-body transformations** for exactly this reason.

In `engeom`, an isometry encodes two things:

- A **rotation**: how much to rotate, and around which axis
- A **translation**: how far to move after rotating

Isometries are equivalent to 4x4 homogeneous transformation matrices, but with the constraint that the matrix must
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

**Rust:**

```rust
use engeom::{Iso2, Vector2};
use std::f64::consts::PI;

// Translate by (1, 2) and rotate by π/4 radians
let iso = Iso2::new(Vector2::new(1.0, 2.0), PI / 4.0);

// Translation only
let t = Iso2::translation(1.0, 2.0);

// Rotation only
let r = Iso2::rotation(PI / 4.0);
```

**Python:**

```python
from math import pi
from engeom.geom2 import Iso2

# Translate by (1, 2) and rotate by π/4 radians
iso = Iso2(1, 2, pi / 4)

# Identity
i0 = Iso2.identity()
```

### 3D Isometries

3D isometries are more involved because rotations in three dimensions require three degrees of freedom. Several
construction paths are available.

**Rust:**

In nalgebra, `Iso3::new` takes a translation vector and an axis-angle vector. The axis-angle vector's direction is
the rotation axis and its magnitude is the rotation angle in radians.

```rust
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

// Translate by (1, 2, 3) with no rotation
let t = Iso3::translation(1.0, 2.0, 3.0);

// Rotate by π/4 around the x-axis (axis-angle encoding)
let r = Iso3::rotation(Vector3::x() * (PI / 4.0));

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

In Rust these are provided by the `IsoExtensions3` trait and are prefixed with `try_`, returning a `Result`.

```rust
use engeom::geom3::IsoExtensions3;
use engeom::{Iso3, Point3, Vector3};

// Build a frame where e0 becomes X and e1 (after orthogonalization) becomes Y
let e0 = Vector3::new(1.0, 1.0, 0.0); // will become the X axis
let e1 = Vector3::new(0.0, 1.0, 0.0); // Y axis derived from this
let origin = Some(Point3::new(1.0, 2.0, 3.0));

let frame = Iso3::try_from_basis_xy(&e0, &e1, origin).unwrap();
```

**Python:**

In Python the methods are static methods on `Iso3`, without the `try_` prefix. They raise a `ValueError` if the
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
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

let iso = Iso3::rotation(Vector3::x() * (PI / 4.0));
let inv = iso.inverse();
```

**Python:**

```python
from math import pi
from engeom.geom3 import Iso3

iso = Iso3.from_rotation(pi / 4, 1, 0, 0)
inv = iso.inverse()
```

## Composing Isometries

Isometries can be composed by multiplication. The result is a single isometry equivalent to applying the right-hand
operand first and then the left-hand operand: the same convention as matrix multiplication.

Importantly, isometry multiplication is **not commutative**: rotating then translating produces a different result
than translating then rotating.

**Rust:**

```rust
use engeom::{Iso3, Vector3};
use std::f64::consts::PI;

let r = Iso3::rotation(Vector3::x() * (PI / 4.0));
let t = Iso3::translation(1.0, 2.0, 3.0);

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
