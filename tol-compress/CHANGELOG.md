# Changelog

This crate releases on its own cadence, under `tol-compress-vX.Y.Z` tags.

Versions below 1.0 follow Cargo's convention that the minor position carries breaking changes, so
`0.1` and `0.2` are not compatible and Cargo will not substitute one for the other.

## 0.2.0

### The file format changed, and 0.1 files cannot be read

The container's format version is now **2**. A file written by 0.1 is refused with
`Error::UnsupportedVersion(1)` rather than misread, and **this release contains no reader for it**.
If you are holding 0.1 files, decode them with 0.1 before upgrading.

The point block gained partitions and moved partition corners off `f64`, which is what forced the
break. It now records one anchor box for the whole block and stores each partition's corners as
integer codes on that anchor's lattice. An `f64` corner places a value to a relative precision of
1e-16 when the lattice it anchors resolves to about `2 * tol`, so the remaining bits were waste. 
This didn't make much of a difference before partitioning blocks, but when implementing that feature
the cost of the header had to be paid back by the block and so it started to matter.  With the
change a 3D partition header went from 56 bytes to about 20.

The reserved multiple-partition field is now in use, and so is the per-partition transform flag,
though nothing sets it yet; see below. Compression and per-item attributes remain reserved and are
still rejected if set.

### `Effort`

Rather than start to accumulate a zoo of different settings and flags that are unique to each of
the different encoding features, I decided to set up a single "effort" concept, which will be
reused across the different techniques.  The goal here was:

- Prevent a user from having to remember what all of the different flags and options would do
  across every different situation.
- Simplify the public facing API by taking away clutter.
- Stabilize the public facing API by mapping new features into the existing "effort" levels, or
  giving us the ability to shuffle features around based on what we learn and how they change
  with improvements.

The mechanism is now in the enum `Effort::{Quick, Balanced, Thorough}`, defaulting to `Balanced`.

Three rules hold to try to make this safe despite the complexity it abstracts away:

- It never changes what comes back. Semantic settings such as `mesh::VertexOrder` stay separate.
- It never reaches the reader. The level is not recorded in the file.
- It never costs size. A higher level never produces a larger file than a lower one.

(Right now, `Effort::Quick` writes a single partition, which is the shape 0.1 wrote, though its
bytes are different.)

### Partitioning

A point block may now hold several partitions, each with its own box and widths, so a set whose
points bunch up in places no longer charges every point for the extent of the whole.

Partitions are contiguous runs of the sequence the caller supplied. **Nothing is reordered**, which
is what lets this apply to every container: a polyline's order is the curve itself, a cloud's order
is documented as preserved, and a mesh's vertex order is the first-use order its index coding needs.

Measured on large scans that still arrive in acquisition order, cutting that order into partitions is
worth 19% of a turbine stator's coordinates and 31% of an unthinned phone scan's. Small or
spatially incoherent inputs gain little or nothing, and are written as a single partition.

### Partition transforms, carried but idle

A partition can be written in a frame of its own rather than on world axes, which for a scanned
surface means the tangent frame of a patch, where one axis collapses to how far the patch bends away
from its own plane. The flag names it, the writer emits a rotation when handed one, and the reader
applies it. **Nothing chooses a frame**, so no file written by this release has the flag set.

Only the rotation is stored, 38 bits in 3D and 12 in 2D. The frame's origin is the anchor's centre
and its corners sit on the anchor rotated into the same frame, both of which the decoder derives.
Rotation precision does not enter the tolerance budget: the encoder quantizes the rotation first and
then works in the frame that defines, so a coarse rotation costs a slightly looser bounding box but
doesn't affect the accuracy of the points.

It's unused right now because of stability rather than size. Two estimators were built and measured, a
covariance eigen fit and DiTO-14, and they came out similar: -4.9% of coordinates and -1.9% of a
whole file on a 370k point turbine stator, -14.3% and -3.8% on an unthinned phone scan. But writing
a file, reading it back and writing it again stopped settling. Without frames that cycle closes in
two to four passes with every point inside about 1.5 tolerances of where it started; with frames it
did not close at all on some inputs and points wandered at `sqrt(passes)`, past thirty tolerances and
still climbing. No individual write breaks its guarantee, and the behaviour note further down still
holds for what this release actually produces, but I was hoping to not have a format where the 
store-and-restore cycle could walk data around indefinitely.

Two causes were found and fixed, one of which was a bug in the axis-aligned path as well: the 
outward corner snapping compared exactly, so a corner landing on a lattice point could grow its box 
by a whole step depending on which side the arithmetic fell.  There's at least one third cause which
I didn't find.

I'm not sure where to go from here on this, so I'm leaving transformation mechanism in place while
I decide what's worth trying to squeeze out of the space that's left.  I also may revisit the idea
of the store/restore error settling, because right now I've only tested it on cases where you load
and save the same data without modification, and that's probably not a very realistic case, since
if you haven't modified a file you wouldn't be saving it again.  At some point I'll test what happens
to modified data, and then decide if keeping the settling guarantee is actually valuable.

### Breaking API changes

- `tol_compress::WriteOptions` is gone from the crate root. There are now three of them, one per
  container, and a bare name that silently meant the mesh one was misleading. Use
  `mesh::WriteOptions`, `cloud::WriteOptions` or `polyline::WriteOptions`.
- `mesh::WriteOptions` is `#[non_exhaustive]` and carries a new `effort` field, so it can no longer
  be built with a struct literal from outside the crate. Use `WriteOptions::new()` with the `with_`
  methods, or `Default`.
- `write_points` still has its old signature but no longer produces the same bytes, since it now
  partitions by default. `write_points_with(.., Effort::Quick)` is the nearest equivalent to the old
  behaviour in shape, though not byte for byte.

### Added

- `effort` module and `Effort`.
- `segment` module: `plan`, `Plan`, `Run`, `stride`, `corner_codes`, `decode_corners`, `LOOKBACK`.
- `transform` module: `Rotation` with `from_axes`, `from_direction`, `to_local`, `to_world`, plus
  `supported`, `bits`, `rotated_anchor` and `centre_of`. Nothing in the crate calls the builders, so
  this is here for the reader, for anyone implementing the format, and for whoever picks the
  estimator back up.
- `points::write_points_with` and `points::partition_bytes`.
- `quantize::width_fits`, the format's width rule as a single predicate, so a caller stepping widths
  upward and `bits_for_tol` searching from scratch cannot disagree.
- `cloud::WriteOptions`, `polyline::WriteOptions`, and `write_to_with` on both.
- `write_file_with` on all three containers. Previously the options were reachable only when writing
  to a `Write`, so a caller using a path could not set them at all.

### Behaviour worth knowing

**Re-encoding is not idempotent.** Where partitions go is decided from where the points are, and the
quantizer moves them, so a file read back and rewritten is not byte-identical to itself. It reaches a
short closed orbit instead, usually a fixed point and otherwise a pair, within a few rewrites.
Because a repeat means the decoded points repeated exactly, error cannot accumulate through a
pipeline that stores and restores. Every file in such an orbit is within tolerance of its input.

**Renumbering a mesh no longer always pays.** `reorder` orders vertices for the index block, by
traversing connectivity. The points block now wants them ordered for space. On a mesh whose
connectivity is a poor guide to position, the traversal can scramble a spatially coherent input and
cost the points block more than the index block gains. It still pays on the great majority of meshes,
and `write_to_preserving_order` is worth measuring against if you have one that it does not.

## 0.1.0

Initial release. Tolerance-bounded encoding of points, meshes, polylines and point clouds in 2D and
3D, with exact per-axis bit widths, exact index coding at whichever of two codings measures smaller,
and connectivity-derived vertex and face renumbering.
