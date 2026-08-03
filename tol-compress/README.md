# tol-compress

Tolerance-bounded compression of practical spatial coordinate data.

The `tol-compress` library performs aggressive, practical compression of common spatial structures like point clouds, meshes, 1-simplices, vector and unit-vector fields in two and three dimensions.  Its name comes from the guaranteed recovery tolerance it uses to determine the amount of compression possible.

## Differences from other libraries

This library is similar to `cloudini` or `lepcc`, but supports more than three dimensional point clouds and attempts more aggressive compression.

Compared to Google's Draco, which to my knowledge is the state of the art for point clouds and triangles:

- Draco is aimed towards the 3d visual asset use case, while `tol-compress` is focused on suitability for storing measured and/or engineering data.
- The point of the `tol-compress` library is to give an end-to-end recovery tolerance guarantee, so the library will take an allowable tolerance from the user and determine the appropriate bit-width on its own, rather than requiring the caller to figure out what it must be.
- All in all, the `tol-compress` library is better described as a simple file format than a compression codec. It forgoes some of the more aggressive optimizations in the name of simplicity.
<!-- TODO: build this section out more once we figure out what things we can actually borrow effectively -->

## Theory

The `tol-compress` library takes advantage of two specific features of real-world measured spatial data:

1. Real world measurement systems have known precisions below which differences in position do not represent meaningful information, and users of measurement data have knowledge about the smallest differences in values that have relevance to their use case.

2. Values produced by the measurement of physical objects rarely span more than a few orders of magnitude more than the smallest meaningful precision for the measurement system and/or application.

Put into plain language, sensing systems are only so accurate and storing precision well below their noise floor or below what makes any difference for your use case is just wasted space.  Furthermore, most measurement systems can't reliably measure a meter long object with better than micron precision, or a kilometer long structure with millimeter precision...and if they could it wouldn't be all that meaningful anyway (to understand why, consider that even the coefficient of thermal expansion for most common materials is between `1e-5` and `1e-4`).

That said, the reason that 64-bit floating points (aka `double` or `f64`) are still the go-to for these applications is because floating point precision varies by _where_ on the number line a value is, and the further you are from the origin the worse the precision is.  Using 64 bits essentially allows you to not worry about where points are while still maintaining the integrity of the value, but it makes for wasteful storage.

To leverage this, the `tol-compress` library will attempt to reorient and/or cluster spatial positions, then encode dimensions separately as an unsigned integer value between double precision bounds using the smallest possible bit width to guarantee that all positions are recovered within a user-supplied tolerance.

Theoretically, the tolerance based compression will do the best when the ratio of the range of positions (max - min) to the acceptable recovery tolerance is lower (five orders of magnitude, for example, allows encoding as 16-bit values, taking 25% of the space of a `double`).  As the ratio of range to tolerance increases, the tolerance based compression will approach the storage efficiency of regular single and double precision based formats, rather than exceeding them.

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





