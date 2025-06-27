# SVD Basis

One of the applications
of [singular value decomposition (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition) is to find a set of
orthonormal basis vectors for a larger set of vectors. This is closely related to the concept of principal component
analysis (PCA) in statistics.

If given a set of points in $\mathbb{R}^n$, the singular value decomposition can be used to find a set of $n$
orthonormal basis vectors describing the primary directions of the data. The first basis vector will be the direction
which accounts for the most variance in the data, and the last basis vector will be the direction which accounts for the
least variance.

These basis vectors have a number of different useful applications, the most obvious of which is to use them to
create an isometry that maps the original data to a new coordinate system where the basis vectors are the axes. In the
case of a large number of co-planar points in 3-space, for example, this could be used to efficiently find a
transformation to map them to the XY plane and back. In two or three dimensions, when coupled with careful resampling,
it can be used to best fit lines or planes, or to rotate a set of points or geometry to align them with the cartesian
axes.

## Finding the Basis

To use the SVD basis in `engeom`, import either the 2D or 3D version of the class from the appropriate `geom` module.
The class will need to be constructed from a numpy array of shape `(n, m)` where `n` is the number of points and `m` is
the dimensionality of the points.

In the following example, we'll use 10000 random points. X and Z will be random, but Y will be a linear function of
X and Z with some noise added. We'll then demonstrate both the 2D and 3D versions of the SVD basis.

```python
import numpy
from engeom.geom2 import SvdBasis2
from engeom.geom3 import SvdBasis3

# Create 10000 random points in 3-space
xs = numpy.random.rand(10000)
ys = xs + numpy.random.rand(10000) * 0.1
zs = numpy.random.rand(10000)
points = numpy.vstack((xs, ys, zs)).T

# Create the 2D SVD basis
b2 = SvdBasis2(points[:, :2])
print(b2.largest())
# Will print something like: Vector2(0.704, 0.709)

# Create the 3D SVD basis
b3 = SvdBasis3(points)
print(b3.largest())
# Will print something like: Vector3(0.704, 0.709, 0.0)
```

The `SvdBasis` constructors also have an optional `weights` argument. To use it, pass in a numpy array of the same
length as the number of points, where each element is a floating point weight associated with the point at the same
index. The weights will be used to scale the power of each point in the SVD calculation. If the weights are not
specified, it's the same as passing in an array of ones.

## Using the Basis

Once created, the `SvdBasis` objects have a number of methods related to the basis result.

- `largest()` returns the largest basis vector.
- `smallest()` returns the smallest basis vector.
- `rank(tol)` returns the rank of the basis, which is the number of singular values greater than a tolerance `tol`. See
  the function documentation for more details on the meaning.
- `basis_variances()` returns the variances associated with the basis vectors. The first element is the variance of the
  largest basis vector, and so on.
- `basis_stdevs()` returns the standard deviations associated with the basis vectors. The first element is the standard
  deviation of the largest basis vector, and so on.

Finally, the `to_iso2()` and `to_iso3()` methods can be used to create an isometry that maps the original data in the
world coordinate system to the new coordinate system defined by the basis vectors. This is useful for transforming the
data to a new coordinate system where the basis vectors are the axes.

```python 
import numpy
from engeom.geom3 import SvdBasis3

# Create four points at the corners of a tall rectangle
points = numpy.array([[0, 0, 0],
                      [1, 1, 0],
                      [1, 1, 5],
                      [0, 0, 5]]).astype(numpy.float64)

b = SvdBasis3(points)
iso = b.to_iso3()

moved = iso.transform_points(points)
print(moved)

# Will print something like:
# [[ 2.50000000e+00 -7.07106781e-01  3.10862447e-15]
#  [ 2.50000000e+00  7.07106781e-01  2.66453526e-15]
#  [-2.50000000e+00  7.07106781e-01 -3.10862447e-15]
#  [-2.50000000e+00 -7.07106781e-01 -2.22044605e-15]]
```

In the above example, the original points formed the corners of a rectangle 5 units tall in the Z direction, $\sqrt{2}$
units wide, and rotated 45 degrees around the Z axis. 

After the basis was found and the points transformed by the isometry, the rectangle is now centered at the origin,
the long direction goes from -2.5 to +2.5, and the short direction goes from -0.707 to +0.707. The Z values are all
essentially zero.

By inverting the isometry, you can move points from the basis coordinate system back to the world coordinate system.
