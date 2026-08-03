# tol-compress

Tolerance-bounded compression of practical spatial coordinate data.

The `tol-compress` library performs aggressive, practical compression of common spatial structures like point clouds, meshes, 1-simplices, vector and unit-vector fields in two and three dimensions.  Its name comes from the guaranteed recovery tolerance it uses to determine the amount of compression possible.

This library is similar to `cloudini` or `lepcc`, but supports more than three dimensional point clouds and attempts more aggressive compression.

## Theory

The `tol-compress` library takes advantage of two specific features of real-world measured spatial data:

1. Real world measurement systems have known precisions below which differences in position do not represent meaningful information, and users of measurement data have knowledge about the smallest differences in values that have relevance to their use case.

2. Values produced by the measurement of physical objects rarely span more than a few orders of magnitude more than the smallest meaningful precision for the measurement system and/or application.

To leverage this, the `tol-compress` library will attempt to reorient and/or cluster spatial positions, then encode dimensions separately as an unsigned integer value between double precision bounds using the smallest possible bit width to guarantee that all positions are recovered within a user-supplied tolerance.

Theoretically, the tolerance based compression will do the best when the ratio of the range of positions (max - min) to the acceptable recovery tolerance is lower (five orders of magnitude, for example, allows encoding as 16-bit values).  As the ratio increases the tolerance based compression will approach the inefficiency of regular single and double precision based formats.

## Feature List

This is the working list of features that I want to examine adding to the library:

- For point data, use SVD or PCA to reorient the points prior to encoding
- Allow the bit width of each dimension to be different (I believe the `engeom::io` module already does this)
- Allow for bit widths that result in partial bytes, for both the encoding of floating points _and_ for the indices used in simplex connectivity maps.
- Examine spatial clustering, since the worst conditioned cases would involve small clusters of points spread far apart.
- `n-1` dimension encoding for unit vectors, allowing an angle tolerance to be specified 
- Pre-built full formats for meshes and point clouds in 2 and 3 dimensions, including with domain attributes like in `engeom`'s mesh entities.
- Pre-built formats for 1-simplices in 2 and 3 dimensions
- Minimal dependencies, maybe just `nalgebra` to start with





