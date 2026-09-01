# tol-compress

Tolerance-bounded compression of practical spatial coordinate data, no dependencies.

For point clouds, meshes, and polylines, you give it a tolerance in the same units as the coordinates, and it finds the smallest representation that guarantees every coordinate comes back within that distance of where it started.  On real scanned geometry that is roughly 60 to 70% smaller than the same data stored as `f32` coordinates with `u32` indices.

```toml
[dependencies]
tol-compress = "0.2"
```

---

## Quick start

```rust
use tol_compress::{Mesh3, mesh};

fn main() -> Result<(), tol_compress::Error> {
    // `tol-compress` has no concept of units, so you just need to be consistent between
    // your geometry and your tolerance. If these positions are in meters, then 1e-6 is
    // a micron.  If they were millimeters, 1e-3 would be a micron.

    //  We'll assume these numbers are in meters
    let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = vec![[0u32, 1, 2]];
    let path = std::env::temp_dir().join("bracket.tcmesh");

    // Every vertex is guaranteed to come back within one micron of where it went in.
    let mesh = Mesh3::new(points, faces).named("bracket");
    mesh::write_one_file(&path, &mesh, 1e-6)?;

    let back = mesh::read_one_file(&path)?;
    assert_eq!(back.name.as_deref(), Some("bracket"));
    Ok(())
}
```

Point clouds and polylines work the same way, in two or three dimensions:

```rust
use tol_compress::{Cloud3, Polyline2};

let cloud = Cloud3::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
let outline = Polyline2::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], true); // closed
```

Every container is a *collection*, so a file may hold one item or many, each optionally named and carrying its own metadata.  The `write_one_*` and `read_one_*` functions are conveniences over that rather than a different format; reading one item out of a file that holds several is an error rather than a silent truncation.

## Practical compression ratios

Measured on test fixtures, `vs f32` is the same geometry stored as single-precision coordinates with `u32` indices, which is roughly what something like a binary single-point precision `.ply` file would be.

| mesh | vertices | faces | tolerance | file | bytes/vertex | vs f32 |
|---|--:|--:|--:|--:|--:|--:|
| Stanford bunny (res 4) | 453 | 948 | 10 um | 4,717 | 10.4 | **-72%** |
| scanner output chunk | 8,091 | 15,539 | 1 um | 90,581 | 11.2 | **-68%** |
| Turbine stator | 369,763 | 734,117 | 1 um | 4,803,640 | 13.0 | **-64%** |
| Pixel 8 Pro, unthinned | 3,742,403 | 7,483,443 | 1 um | 50,445,803 | 13.5 | **-63%** |

The stator is a 127 mm part scanned on a system that usually gets about 10-15 micron repeatability. I compressed it with a 1um tolerance, and it shows how the ratio in the theory section starts to come into play.  Five orders of magnitude between extent and tolerance needs 17 bits per axis, so coordinates cost more per vertex than on the smaller parts even though the file is still 64% smaller than `f32`.

The Pixel 8 Pro scan is especially interesting, because it was polygonized by Zeiss without the standard baked-in thinning settings, so it has an extremely dense triangulation.

The two large meshes are not in this repository.  `engeom`'s `tc_size_row` example takes a path and a tolerance and prints the row, so figures like these can be taken from files that stay wherever they live.

## Differences from other libraries

This library is similar to `cloudini` or `lepcc`, but covers more than point clouds: meshes and polylines are first class kinds of their own, in two dimensions as well as three, and it attempts more aggressive compression.

Compared to Google's Draco, which to my knowledge is the state of the art for point clouds and triangles:

- Draco is aimed towards the 3d visual asset use case, while `tol-compress` is focused on suitability for storing measured and/or engineering data.  This shows up in places like the Edgebreaker algorithm requiring manifold surfaces, and being optimized for smooth meshes.
- The point of the `tol-compress` library is to give an end-to-end recovery tolerance guarantee, so the library will take an allowable tolerance from the user and determine the appropriate bit-width on its own, rather than requiring the caller to figure out what it must be.
- All in all, the `tol-compress` library is better described as a simple file format than a compression codec. It forgoes some of the more aggressive optimizations in the name of simplicity, and in particular it has no entropy coding stage at all, which keeps decoding cheap and keeps the format something you could (hopefully) implement from a written description.

## Theory

The `tol-compress` library takes advantage of two specific features of real-world measured spatial data:

1. Real world measurement systems have known precisions below which differences in position do not represent meaningful information, and users of measurement data have knowledge about the smallest differences in values that have relevance to their use case.

2. Values produced by the measurement of physical objects rarely span more than a few orders of magnitude more than the smallest meaningful precision for the measurement system and/or application.

Put into plain language, sensing systems are only so accurate and storing precision well below their noise floor or below what makes any difference for your use case is just wasted space.  Furthermore, most measurement systems can't reliably measure a meter long object with sub-micron precision, or a kilometer long structure with millimeter precision...and if they could it wouldn't be all that meaningful anyway (to understand why, consider that even the coefficient of thermal expansion for most common materials is between `1e-5` and `1e-4`).

That said, the reason that 64-bit floating points (aka `double` or `f64`) are still the go-to for working with and persisting this type of data is because floating point precision varies by _where_ on the number line a value is, and the further you are from the origin the worse the precision is.  Using 64 bits essentially allows you to not worry about where points are while still maintaining the integrity of the value, but it makes for incredibly wasteful storage.

To leverage this, the `tol-compress` library encodes each dimension separately as an unsigned integer value between double precision bounds, using the smallest bit width that guarantees all positions are recovered within a user-supplied tolerance.

Theoretically, the tolerance based compression will do the best when the ratio of the range of positions (max - min) to the acceptable recovery tolerance is lower (just under 5 orders of magnitude, for example, allows encoding as 16-bit values, taking 25% of the space of a `double`).  As the ratio of range to tolerance increases, the tolerance based compression will approach the storage efficiency of regular single and double precision based formats, rather than exceeding them.

## What's in the file

This section should have enough information to get you started, but the module documentation on [docs.rs](https://docs.rs/tol-compress) carries the byte-level detail, and a full specification will follow once the format is complete.

A file begins with a **12 byte container header**: the magic `TOLC`, a format version, a "kind" byte (mesh, polyline or cloud, in 2D or 3D), a reserved compression byte, flags, and an item count N.  Optional file level **metadata** follows, a sorted key-value map of `Bool`, `I64`, `F64`, `Text` and `Bytes` values.  Metadata is stored and never interpreted, so a caller can record whatever their domain needs without the format learning anything about it.

Then there are N **items**, each with a short preamble (flags, optional name, optional metadata, a reserved attribute count) followed by the geometry blocks its kind calls for:

- A **points block** records one anchor box as `f64` with the bit width per axis the whole set needs, then one or more partitions.  Each partition stores its own corners as codes on the anchor's lattice rather than as `f64`, then its own bit width per axis, then the codes packed back to back with no padding between them.  Widths are exact rather than rounded up to whole bytes, and an axis whose values are all identical costs no bits at all.  Partitioning lets a set whose points bunch up in places avoid charging every point for the extent of the whole, and the encoder cuts the sequence where it already is rather than moving points around, so point order survives.

  Corners are codes rather than `f64` because an `f64` places a value to a relative precision of 1e-16 when the lattice it anchors resolves to about `2 x tol`.  At a typical metrology resolution that takes a 3D partition header from 56 bytes to about 20, and since a cut has to earn its header back, that is also what sets how finely the encoder can afford to cut.
- An **index block**, for meshes, stores triangles either as absolute indices at an exact width or as distances back from a running high-water mark, whichever measures smaller.  The second is much smaller on real meshes.  Both are exact; connectivity is never approximated.

The tolerance itself is deliberately not stored.  It is derivable from the bounds and widths already present, and the derived figure is the guarantee actually achieved.

> [!WARNING]
> **Vertex and face order is not preserved by default!** Unless you tell it not to, the encoder renumbers both so that indices compress, which is where most of the saving on a mesh comes from, and the renumbering is derived entirely from connectivity so nothing about it has to be stored.  The mesh that comes back has the same geometry: the same triangles over the same positions, but consistently renumbered.  If you hold data keyed by vertex index outside the file, either use `write_to_preserving_order` to disable the re-indexing, or reorder your own data first with the `reorder` module.

### Forward compatibility

Several fields are **reserved and rejected if set**: compression and per-item attributes.  A reader from this version meeting a file that uses one of them fails with a specific error rather than misreading it, and files written today stay readable as those features arrive.  That is what those bytes were spent on.

The point block's per-partition transform flag is a step further along: both ends of it are implemented and tested, but no encoder chooses a frame, so no file written by this version sets it.  See below.

The point block layout changed once, at format version 2, when partitions arrived and their corners moved off `f64`.  A version 1 file is refused by name rather than misread, and there is no version 1 writer.

## Status

Working and measured:

- Meshes, polylines and point clouds, in 2D and 3D, as named collections with metadata.
- Exact per-axis bit widths, including zero-width degenerate axes.
- Face and vertex renumbering with high-water-mark index coding, plus per-block adaptive widths.
- A dimension-generic point block: `points::write_points` and `read_points` handle any dimension from 1 to 255 and are tested through 4D.  Containers are a different matter, since the file's kind byte names only `Cloud2`, `Cloud3`, `Polyline2`, `Polyline3` and `Mesh3`, so a cloud above three dimensions has nothing to be written as and `cloud::write_to` panics on one.  Reaching higher dimensions means either using the block encoder directly or spending a new kind byte.
- Partitioning of a point sequence into runs that each carry their own box and widths, cutting the caller's order in place so that nothing is reordered.  Worth a third of the coordinates on a structured surface and rather more on scattered clusters, and nothing at all where the sequence has no locality.
- A single `Effort` setting, `Quick`, `Balanced` or `Thorough`, standing in for how hard the encoder searches rather than a flag per technique.  It never changes what comes back, never reaches the reader, and a higher level never produces a larger file.

Not yet implemented, with the format space reserved:

- Per-element attributes: scalars, labels, colors, vectors, and unit vectors with an angle tolerance.
- Partitioning that is allowed to reorder points rather than only cut the order it was given.

Built and measured, but not switched on:

- **Per-partition transforms**, which write a partition in a frame of its own so that a patch of scanned surface is measured in its own tangent frame instead of against the world axes.  The format, the writer and the reader are all done; what is missing is the estimator that picks a frame.  Two were tried, a covariance eigen fit and DiTO-14, and they measured level with each other at -4.9% of the coordinates and -1.9% of a whole mesh file on a turbine stator, -14.3% and -3.8% on an unthinned phone scan.  The problem is that re-encoding stopped settling: without frames, writing a file and reading it back a few times converges with every point within about 1.5 tolerances of where it started, and with frames some inputs never converged and drifted past thirty.  Every individual write still met its tolerance, but I had thought that re-encoding stabilization was worth trying to guarantee and this violated it.  I'll revisit this after I've had time to think it through and try some other stuff.

Considered and declined after testing/measurement, recorded here for posterity:

- I looked at trying to reclaim some of the bits for encoding point blocks, since for an N-dimensional space I currently shrink the tolerance given to each axis so that in cartesian space the total magnitude sums to the user-specified tolerance.  That means that usually the total tolerance stays well under the limit, since each axis assumed the other would completely fill its range.  I tried doing a greedy reclaimation of bits, allowing an axis to shrink its bit width so long as the other axes had enough room left in their tolerance that the total still stayed under the guarantee.  It worked, but it took very little (usually less than 1%) off of the file size for a lot of complexity, so I didn't think it was worth it.

- I originally planned to have an optional compression dependency that would let you activate compression on the blocks.  The `tol-compress` crate is an overhaul and expansion of a similar but much simpler tool in the `engeom` library, which in turn was based on a tiny quantized mesh format I had made for embedding small meshes into `engeom`'s test suite. For both of those implementations, adding a layer of compression on top of it shrunk the files by a decent amount, and seemed worth trying. What I didn't realize was that `flate2` worked well because I had only been quantizing to the nearest whole byte.  Once we tried adding compression to the `tol-compress` outputs, we saw very minimal size reduction (in one case `gzip` actually made it bigger). If you find yourself with some data that does benefit from compression, and/or you're trying to squeeze as much space out as you can, you can compress the file as a whole.

- My understanding of Google's Draco is that it does some coordinate prediction from well formed neighbor triangles, encoding the actual point as a residual which can then get away with a much smaller bit width.  We did some rough experimentation and it didn't show much benefit, maybe between +1% and -8% on real meshes, and -14% with exception coding, against the 30 to 45% that index coding gives on the same files. The residual distribution ends up with a heavy tail that per-block widths can't really exploit.

## License

Dual licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
