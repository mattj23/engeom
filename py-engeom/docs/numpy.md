# Use of Numpy

Numpy is an extension library for the Python programming language which offers support for large, multidimensional arrays and matrices and a number of linear algebra operations that can be performed on them.  It's effectively Python's standard go-to library for linear algebra, and is used heavily by other Python libraries in the scientific computing and data analysis spaces.

The underlying `engeom` Rust library uses `nalgebra` for dense linear algebra operations and `faer` for sparse matrix operations. To provide a consistent interface for the Python bindings, the Python `engeom` library assumes that many users will already be using `numpy` for bulk operations on arrays of points, vectors, face lists, etc, and so it makes sense to map many of the `engeom` features which would operate on points/vectors to operate seamlessly on `numpy.ndarray` objects.

To that end, there are a number of conventions used in `engeom` regarding the use of `numpy` arrays:

1. When mapping large arrays of points/vectors, `numpy.ndarray` objects should be sized as `(n, 2)` or `(n, 3)` for 2D and 3D points/vectors, respectively.  Rows represent individual points/vectors, and columns represent the `x`, `y`, and `z`, components of the points/vectors.

2. Some `engeom` objects will provide access to their data as read-only `numpy.ndarray` objects. For instance, the `Mesh` object will expose its vertex and face data as `numpy.ndarray`s.  The shaping of these arrays will be consistent with the above convention, and floating point data will be 64-bit.

3. Where index-based lists, such as the face vertices lists of the `Mesh` objects are used, these will typically be expected to be unsigned integer types.  Internally, most are stored as Rust's `u32` or `usize` types. If you get an error while constructing objects, such as the `Mesh`, make sure that the data you are passing in is of the correct type.