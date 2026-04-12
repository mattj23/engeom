# Alignments

One of the core features of the Engeom library is the ability to align entities in 2D and 3D. Alignments are typically optimization problems, in which a solver attempts to minimize the difference between two sets of geometric entities by adjusting an isometry.

Engeom's alignment tools are mostly built around Levenberg-Marquardt solving of difference residuals between entities that can be distance matched.  For example, the alignment a set of points to another set of points, or a set of points to simplices like lines or triangle meshes in 2D or 3D space.

## Common Concepts

In order to provide both consistency and flexibility for alignments, the typical parameterization of problems is conceptually the same for both 2D and 3D. There are also a few common concepts to introduce up front.

First, alignments take place between *test entities* and *target entities*.  What the entities are will vary with the different tools, but all will be composed of some sort of geometry. That geometry will be defined with (x, y) or (x, y, z) coordinates of some form or another.

There are no particular guarantees about where the geometry is in space, if they're near the origin and/or if the test and target entities are close to each other. Most types of geometry in the Engeom library can be transformed, but some geometric entities are large and making unnecessary copies of them is a waste of CPU and memory.

To give more control over alignments, Engeom uses three distinct isometries when performing alignments.

1. A local origin, \\( L \\), which is defined in the same coordinate system as the test geometry.
2. An alignment isometry, \\( A \\), which is created from the parameters being optimized by the alignment algorithm. This isometry is constantly changing as the algorithm runs.
3. A working offset isometry, \\( O \\), which moves the test geometry after the alignment isometry has been applied. It is used to provide a starting position, to negate the effect of the local origin, or some combination of the two.

The _total_ isometry provided by the alignment parameters is calculated through matrix multiplication:

\\[ O \times A \times L^{-1} \\]

Conceptually, this means that the test geometry is first moved by \\( L^{-1} \\) so that it sits at the world origin with the same relationship it originally had to the local origin. Then the transform \\( A \\) from the alignment parameters is applied.  Finally, the geometry is moved by the working offset transform \\( O \\) to its final position.

There are two special cases to think about:

- If the transforms \\( L \\) and \\( O \\) are equal to the identity matrix, the alignment reduces to translation and rotation about the origin.
- If the transforms are equal such that \\( L = O \\), the result is that the alignment transform behaves as if it's acting directly at the local origin \\( L \\), meaning that the local origin is the center of rotation and translations happen in the directions of the X, Y, and Z axes of \\( L \\).
