# tol-compress

Tolerance-bounded compression of practical spatial coordinate data, no dependencies.

For point clouds, meshes, and polylines (for now), you give it a tolerance in the same units as the coordinates, and it finds the smallest representation that guarantees every coordinate comes back within that distance of where it started.  On real scanned geometry that is roughly 60 to 70% smaller than the same data stored as `f32` coordinates with `u32` indices, and smaller than a gzipped 16-bit quantized format even without being compressed at all.

```toml
[dependencies]
tol-compress = "0.1"
```

## Quick start

```rust
use std::path::Path;
use tol_compress::{Mesh3, mesh};

fn main() -> Result<(), tol_compress::Error> {
    let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = vec![[0u32, 1, 2]];

    // Every vertex is guaranteed to come back within one micron of where it went in.
    let mesh = Mesh3::new(points, faces).named("bracket");
    mesh::write_one_file(Path::new("bracket.tcmesh"), &mesh, 1e-6)?;

    let back = mesh::read_one_file(Path::new("bracket.tcmesh"))?;
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
| Stanford bunny (res 4) | 453 | 948 | 10 um | 4,704 | 10.4 | **-72%** |
| scanner output chunk | 8,091 | 15,539 | 1 um | 92,335 | 11.4 | **-67%** |
| Turbine stator | 369,763 | 734,117 | 1 um | 5,235,443 | 14.2 | **-60%** |

The stator is a 127 mm part scanned on a system that usually gets about 10-15 micron repeatability. I compressed it with a 1um tolerance, and it and shows how the ratio in the theory section starts come into play.  Five orders of magnitude between extent and tolerance needs 17 bits per axis, so coordinates cost more per vertex than on the smaller parts even though the file is still 60% smaller than `f32`.

Against a gzipped 16-bit quantized format holding the same meshes, uncompressed `tol-compress` output is 9 to 26% smaller.

## Differences from other libraries

This library is similar to `cloudini` or `lepcc`, but supports more than three dimensional point clouds and attempts more aggressive compression.

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

## What is in a file

This section should have enough information to get you started.  The module documentation on [docs.rs](https://docs.rs/tol-compress) carries the byte-level detail, and a full specification will follow once the format is complete.

A file begins with a **12 byte container header**: the magic `TOLC`, a format version, a "kind" byte (mesh, polyline or cloud, in 2D or 3D), a reserved compression byte, flags, and an item count.  Optional file level **metadata** follows, a sorted key-value map of `Bool`, `I64`, `F64`, `Text` and `Bytes` values.  Metadata is stored and never interpreted, so a caller can record whatever their domain needs without the format learning anything about it.

Then that many **items**, each with a short preamble (flags, optional name, optional metadata, a reserved attribute count) followed by the geometry blocks its kind calls for:

- A **points block** records the bounding box as `f64`, then a bit width per axis chosen so the round-trip error meets the tolerance, then the codes packed back to back with no padding between them.  Widths are exact rather than rounded up to whole bytes, and an axis whose values are all identical costs no bits at all.
- An **index block**, for meshes, stores triangles either as absolute indices at an exact width or as distances back from a running high-water mark, whichever measures smaller.  The second is much smaller on real meshes.  Both are exact; connectivity is never approximated.

The tolerance itself is deliberately not stored.  It is derivable from the bounds and widths already present, and the derived figure is the guarantee actually achieved.

> [!WARNING]
> **Vertex and face order is not preserved by default!** Unless you tell it not to, the encoder renumbers both so that indices compress, which is where most of the saving on a mesh comes from, and the renumbering is derived entirely from connectivity so nothing about it has to be stored.  The mesh that comes back is the same mesh: the same triangles over the same positions, consistently renumbered.  If you hold data keyed by vertex index outside the file, either use `write_to_preserving_order` or reorder your own data first with the `reorder` module.

### Forward compatibility

Several fields are **reserved and rejected if set**: compression, per-item attributes, a point block transform, and multiple partitions.  A reader from this version meeting a file that uses one of them fails with a specific error rather than misreading it, and files written today stay readable as those features arrive.  That is what those bytes were spent on.

## Status

Working and measured:

- Meshes, polylines and point clouds, in 2D and 3D, as named collections with metadata.
- Exact per-axis bit widths, including zero-width degenerate axes.
- Face and vertex renumbering with high-water-mark index coding, plus per-block adaptive widths.
- Point clouds of more than three dimensions, since the point encoder is dimension-generic.

Not yet implemented, with the format space reserved:

- Per-element attributes: scalars, labels, colors, vectors, and unit vectors with an angle tolerance.
- Reorientation of point sets before encoding, and spatial partitioning of scattered clusters.

Considered and declined after testing/measurement, recorded here for posterity:

- I originally planned to have an optional compression dependency that would let you activate compression on the blocks.  The `tol-compress` crate is an overhaul and expansion of a similar but much simpler tool in the `engeom` library, which in turn was based on a tiny quantized mesh format I had made for embedding small meshes into `engeom`'s test suite. For both of those implementations, adding a layer of compression on top of it shrunk the files by a decent amount, and seemed worth trying. What I didn't realize was that `flate2` worked well because I had only been quantizing to the nearest whole byte.  Once we tried adding compression to the `tol-compress` outputs, we saw very minimal size reduction (in one case `gzip` actually made it bigger). If you find yourself with some data that does benefit from compression, and/or you're trying to squeeze as much space out as you can, you can compress the file as a whole.

- My understanding of Google's Draco is that it does some coordinate prediction from well formed neighbor triangles, encoding the actual point as a residual which can then get away with a much smaller bit width.  We did some rough experimentation and it didn't show much benefit, maybe between +1% and -8% on real meshes, and -14% with exception coding, against the 30 to 45% that index coding gives on the same files. The residual distribution ends up with a heavy tail that per-block widths can't really exploit.

## License

Dual licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
