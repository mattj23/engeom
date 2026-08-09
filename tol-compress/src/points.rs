//! Encoding and decoding blocks of points in any number of dimensions.
//!
//! This is the core of the format. A block records where its points live and how finely each axis
//! was quantized, then packs the codes back to back with no padding between them.
//!
//! # Block layout
//!
//! ```text
//! u8       dimension N, checked against the caller's expectation
//! u32      partition count; a count of zero ends the block here
//! N x f64  anchor lower corner
//! N x f64  anchor upper corner
//! N x u8   anchor bit width per axis
//!   per partition:
//!     u32                point count, never zero
//!     u8                 transform flag, 0 for identity and 1 for a rotation
//!     N x u8             bit width per axis
//!     bit-packed, padded to a byte boundary:
//!       rotation         only when the flag is 1; 38 bits in 3D and 12 in 2D
//!       2N codes         this partition's lower then upper corner, at the anchor's widths
//!       count x N codes  the points, in axis order, at this partition's widths
//! ```
//!
//! # Partitions
//!
//! A block may hold several partitions, each with its own box and its own widths, which is how a
//! set whose points bunch up in places avoids charging every point for the full extent of the
//! whole. [`crate::segment`] decides where they go, and it cuts the sequence **where it already
//! is** rather than moving points around, so points come back out in the order they went in.
//!
//! Partitions are written in sequence order and the decoder concatenates them, so nothing about the
//! arrangement is stored beyond the boxes themselves. [`Effort::Quick`] writes a single partition
//! over the whole set.
//!
//! A partition's corners are **not** stored as `f64`. The block records one anchor box that encloses
//! everything, along with the widths the set as a whole would need, and each partition's corners are
//! codes on that lattice. An `f64` corner places a value to a relative precision of 1e-16 when the
//! lattice it anchors resolves to about `2 * tol`, so the other 47 bits an axis were buying nothing.
//! At a typical metrology resolution this takes a 3D partition header from 56 bytes to about 20,
//! and since a cut has to earn its header back, that is also what sets how finely
//! [`crate::segment`] can afford to cut.
//!
//! Corners are rounded **outward** onto the anchor lattice, so the box a partition declares still
//! encloses every point in it. A point outside its own box would be clamped to the boundary by the
//! quantizer, which is an error no bit width can bound.
//!
//! # Transformed partitions
//!
//! A partition may be written in a frame of its own rather than in world axes, which is what the
//! transform flag says. Only the rotation is stored: its origin is the anchor's centre and its
//! corners sit on the anchor rotated into that same frame, both of which the decoder derives.
//! [`crate::transform`] holds the encoding and the reasoning, including why the rotation's own
//! precision does not enter the tolerance budget.
//!
//! **Nothing chooses a frame at the moment**, so the flag is always written clear and no file this
//! crate produces carries a rotation. Both directions are implemented and tested all the same, so
//! picking the work back up means writing an estimator rather than rebuilding the format.
//! [`crate::transform`] records what was tried, what it was worth, and why it is idle.
//!
//! The flag is `u8` and only 0 and 1 are defined. Anything else is refused rather than ignored,
//! which is what keeps the byte available for a translation or a scale later.
//!
//! # Tolerance budget
//!
//! `tol` bounds the Euclidean distance between an original point and its recovered position. Each
//! axis is therefore given `tol / sqrt(N)`, so that `N` independent per-axis errors at their worst
//! combine in quadrature to exactly `tol`.

use crate::bits::{BitReader, BitWriter, read_payload};
use crate::bounds::Bounds;
use crate::effort::Effort;
use crate::error::{Error, Result};
use crate::quantize::{Quantizer, bits_for_tol};
use crate::raw::{MAX_PREALLOC, read_f64, read_u8, read_u32, write_f64, write_u8, write_u32};
use crate::segment;
use crate::transform::{self, Rotation};
use std::io::{Read, Write};

/// The number of bytes one partition of `count` points at `widths` occupies, against a block whose
/// anchor is stored at `anchor_widths`.
///
/// A partition always takes a whole number of bytes: its header is byte oriented and
/// [`BitWriter`] pads the payload out to a boundary, so there is no fractional part to carry
/// between partitions.
///
/// This exists for the same reason [`crate::blocks::encoded_bits`] does. Anything deciding where to
/// cut a point set into partitions has to score its candidates against what the writer will
/// actually emit.
pub fn partition_bytes<const N: usize>(
    count: usize,
    widths: &[u8; N],
    anchor_widths: &[u8; N],
    transform: Option<&Rotation<N>>,
) -> u64 {
    // Count, transform flag, and a width byte per axis.
    let header = 4 + 1 + N as u64;

    // The rotation, the corners and the points share one bitstream, so the lot is padded to a byte
    // once at the end rather than three times.
    let rotation = match transform {
        Some(_) => u64::from(transform::bits(N)),
        None => 0,
    };
    let corners: u64 = anchor_widths.iter().map(|&w| 2 * u64::from(w)).sum();
    let per_point: u64 = widths.iter().map(|&w| u64::from(w)).sum();

    header + (rotation + corners + count as u64 * per_point).div_ceil(8)
}

/// The narrowest per-axis widths that recover any point in `bounds` to within `tol`.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis spans a range too wide to meet `tol`.
pub(crate) fn widths_for<const N: usize>(bounds: &Bounds<N>, tol: f64) -> Result<[u8; N]> {
    // Splitting the budget in quadrature: N axes each within tol/sqrt(N) put the point within tol.
    let axis_tol = tol / (N as f64).sqrt();
    let extents = bounds.extents();

    let mut widths = [0u8; N];
    for (i, w) in widths.iter_mut().enumerate() {
        *w = bits_for_tol(extents[i], axis_tol)?;
    }

    Ok(widths)
}

/// Write a block of points, choosing the narrowest per-axis widths that meet `tol`.
///
/// `tol` is the largest acceptable Euclidean distance between any original point and the position
/// recovered from the block.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis spans a range too wide to meet `tol`, and
/// [`Error::Malformed`] if any coordinate is not finite.
///
/// # Panics
///
/// Panics if `N` is zero or above 255.
pub fn write_points<W: Write, const N: usize>(
    writer: &mut W,
    points: &[[f64; N]],
    tol: f64,
) -> Result<()> {
    write_points_with(writer, points, tol, Effort::default())
}

/// Write a block of points at a given search effort. See [`write_points`].
///
/// [`Effort::Quick`] writes the set as a single partition, which is what this block has always
/// been. The levels above it hand the sequence to [`crate::segment`] first, which cuts it into
/// partitions wherever that measures smaller. Point order is untouched at every level.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis spans a range too wide to meet `tol`, and
/// [`Error::Malformed`] if any coordinate is not finite.
///
/// # Panics
///
/// Panics if `N` is zero or above 255.
pub fn write_points_with<W: Write, const N: usize>(
    writer: &mut W,
    points: &[[f64; N]],
    tol: f64,
    effort: Effort,
) -> Result<()> {
    assert!((1..=255).contains(&N), "dimension {N} is out of range");

    write_u8(writer, N as u8)?;

    let plan = segment::plan(points, tol, effort)?;

    // An empty set has no anchor box, and inventing one would be a lie a decoder might believe.
    // A partition count of zero ends the block before any of it is written.
    let count = u32::try_from(plan.runs.len())
        .map_err(|_| Error::Malformed("point block holds more partitions than a u32 can count"))?;
    write_u32(writer, count)?;
    if count == 0 {
        return Ok(());
    }

    for i in 0..N {
        write_f64(writer, plan.anchor.mins[i])?;
    }
    for i in 0..N {
        write_f64(writer, plan.anchor.maxs[i])?;
    }
    for &w in plan.anchor_widths.iter() {
        write_u8(writer, w)?;
    }

    for run in &plan.runs {
        write_partition(
            writer,
            &points[run.start..run.end],
            &run.bounds,
            &run.widths,
            &plan.anchor,
            &plan.anchor_widths,
            run.transform.as_ref(),
        )?;
    }

    Ok(())
}

/// The lattice a partition's corners are placed on, which is the anchor turned into the partition's
/// own frame when it has one.
///
/// Both sides of the format compute this the same way from the same inputs, which is what lets a
/// rotated partition reuse the corner coding unchanged.
fn corner_lattice<const N: usize>(
    anchor: &Bounds<N>,
    transform: Option<&Rotation<N>>,
) -> Bounds<N> {
    match transform {
        Some(rotation) => transform::rotated_anchor(anchor, rotation),
        None => *anchor,
    }
}

/// Write one non-empty partition against a box and a set of widths already decided for it.
///
/// The bounds and widths are handed in rather than derived here because the planner had to price
/// this partition before committing to it, and a price computed from a different pass over the data
/// than the bytes are is not a price at all. The bounds must already sit on the corner lattice,
/// which is what [`crate::segment::plan`] guarantees.
///
/// When `transform` is set, `bounds` and `widths` describe the points **in that frame**, not in
/// world axes, and the points are turned into it on the way out.
fn write_partition<W: Write, const N: usize>(
    writer: &mut W,
    points: &[[f64; N]],
    bounds: &Bounds<N>,
    widths: &[u8; N],
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
    transform: Option<&Rotation<N>>,
) -> Result<()> {
    debug_assert!(
        !points.is_empty(),
        "an empty partition has no bounds to write"
    );
    debug_assert!(
        transform.is_none() || transform::supported(N),
        "a frame was fitted at a dimension that cannot store one"
    );

    let count = u32::try_from(points.len())
        .map_err(|_| Error::Malformed("partition holds more points than a u32 can index"))?;
    write_u32(writer, count)?;
    write_u8(writer, u8::from(transform.is_some()))?;
    for &w in widths.iter() {
        write_u8(writer, w)?;
    }

    // Recomputed from the bounds rather than carried alongside them, so the codes written and the
    // box they were derived from cannot come from different passes. The bounds are already on the
    // lattice, so this round-trips them rather than moving them.
    let lattice = corner_lattice(anchor, transform);
    let (lo, hi) = segment::corner_codes(bounds, &lattice, anchor_widths);
    debug_assert_eq!(
        segment::decode_corners(&lo, &hi, &lattice, anchor_widths),
        *bounds,
        "a planned box must already sit on the corner lattice"
    );

    // Rebuilt from the bounds and widths the header just recorded, which is the same lattice the
    // decoder will reconstruct from them.
    let quants: Vec<Quantizer> = (0..N)
        .map(|i| Quantizer::new(bounds.mins[i], bounds.maxs[i] - bounds.mins[i], widths[i]))
        .collect();
    let centre = transform::centre_of(anchor);

    let mut bw = BitWriter::new(&mut *writer);
    if let Some(rotation) = transform {
        rotation.write(&mut bw)?;
    }
    for i in 0..N {
        bw.write_bits(lo[i], anchor_widths[i])?;
    }
    for i in 0..N {
        bw.write_bits(hi[i], anchor_widths[i])?;
    }
    for p in points {
        // The frame is applied before quantization, so the codes describe local coordinates and the
        // widths chosen for them are the ones that were priced.
        let local = match transform {
            Some(rotation) => rotation.to_local(p, &centre),
            None => *p,
        };
        for (i, q) in quants.iter().enumerate() {
            bw.write_bits(q.encode(local[i]), q.bits())?;
        }
    }
    bw.finish()?;

    Ok(())
}

/// Read a block of points written by [`write_points`].
///
/// # Errors
///
/// [`Error::Malformed`] if the block was written at a different dimension than `N`, or if it uses
/// a feature this version does not implement. I/O errors propagate for truncated input.
pub fn read_points<R: Read, const N: usize>(reader: &mut R) -> Result<Vec<[f64; N]>> {
    assert!((1..=255).contains(&N), "dimension {N} is out of range");

    let dimension = read_u8(reader)?;
    if usize::from(dimension) != N {
        return Err(Error::Malformed(
            "point block was written at a different dimension than requested",
        ));
    }

    let partitions = read_u32(reader)?;
    if partitions == 0 {
        return Ok(Vec::new());
    }

    let mut mins = [0.0f64; N];
    let mut maxs = [0.0f64; N];
    for m in mins.iter_mut() {
        *m = read_f64(reader)?;
    }
    for m in maxs.iter_mut() {
        *m = read_f64(reader)?;
    }
    for i in 0..N {
        if !mins[i].is_finite() || !maxs[i].is_finite() || maxs[i] < mins[i] {
            return Err(Error::Malformed("point block has an invalid anchor box"));
        }
    }
    let anchor = Bounds { mins, maxs };

    let mut anchor_widths = [0u8; N];
    for w in anchor_widths.iter_mut() {
        *w = read_u8(reader)?;
        if *w > 53 {
            return Err(Error::Malformed(
                "point block declares an impossible anchor width",
            ));
        }
    }

    let mut out = Vec::new();
    for _ in 0..partitions {
        read_partition(reader, &anchor, &anchor_widths, &mut out)?;
    }

    Ok(out)
}

fn read_partition<R: Read, const N: usize>(
    reader: &mut R,
    anchor: &Bounds<N>,
    anchor_widths: &[u8; N],
    out: &mut Vec<[f64; N]>,
) -> Result<()> {
    let count = read_u32(reader)? as usize;
    if count == 0 {
        return Err(Error::Malformed("point block holds an empty partition"));
    }

    // Only identity and rotation are defined. An unknown value is refused rather than ignored,
    // since ignoring one would place every point in the partition somewhere plausible and wrong.
    let flag = read_u8(reader)?;
    let transformed = match flag {
        0 => false,
        1 => true,
        _ => {
            return Err(Error::Malformed(
                "point block carries a transform this version does not define",
            ));
        }
    };
    if transformed && !transform::supported(N) {
        return Err(Error::Malformed(
            "point block carries a transform at a dimension that has no encoding for one",
        ));
    }

    let mut widths = [0u8; N];
    for w in widths.iter_mut() {
        *w = read_u8(reader)?;
        if *w > 53 {
            return Err(Error::Malformed(
                "point block declares an impossible bit width",
            ));
        }
    }

    out.reserve(count.min(MAX_PREALLOC));

    // The reader must not read past the payload, because a byte-oriented header follows it. The
    // flag, the widths and the count fully determine how long it is, rotation and corners included.
    let rotation_bits = if transformed {
        u64::from(transform::bits(N))
    } else {
        0
    };
    let corners: u64 = anchor_widths.iter().map(|&w| 2 * u64::from(w)).sum();
    let per_point: u64 = widths.iter().map(|&w| u64::from(w)).sum();
    let payload = (rotation_bits + corners + count as u64 * per_point).div_ceil(8) as usize;

    let bytes = read_payload(reader, payload)?;
    let mut br = BitReader::new(&bytes);

    let transform = if transformed {
        Some(Rotation::<N>::read(&mut br)?)
    } else {
        None
    };
    let lattice = corner_lattice(anchor, transform.as_ref());

    let mut lo = [0u64; N];
    let mut hi = [0u64; N];
    for (i, code) in lo.iter_mut().enumerate() {
        *code = br.read_bits(anchor_widths[i])?;
    }
    for (i, code) in hi.iter_mut().enumerate() {
        *code = br.read_bits(anchor_widths[i])?;
    }

    let bounds = crate::segment::decode_corners(&lo, &hi, &lattice, anchor_widths);
    for i in 0..N {
        if bounds.maxs[i] < bounds.mins[i] {
            return Err(Error::Malformed("partition has an inverted box"));
        }
    }

    let quants: Vec<Quantizer> = (0..N)
        .map(|i| Quantizer::new(bounds.mins[i], bounds.maxs[i] - bounds.mins[i], widths[i]))
        .collect();
    let centre = transform::centre_of(anchor);

    for _ in 0..count {
        let mut p = [0.0f64; N];
        for (i, q) in quants.iter().enumerate() {
            p[i] = q.decode(br.read_bits(q.bits())?);
        }
        // The codes are local coordinates when a frame is set, so the point comes back out of it.
        if let Some(rotation) = &transform {
            p = rotation.to_world(&p, &centre);
        }
        out.push(p);
    }
    br.finish()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testgen::Rng;
    use std::io::Cursor;

    /// Euclidean distance, which is what `tol` actually bounds.
    fn distance<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
        (0..N).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    fn round_trip<const N: usize>(points: &[[f64; N]], tol: f64) -> Vec<[f64; N]> {
        let mut buf = Vec::new();
        write_points(&mut buf, points, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered: Vec<[f64; N]> = read_points(&mut cursor).unwrap();

        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder did not consume the whole block"
        );
        assert_eq!(recovered.len(), points.len());
        recovered
    }

    fn assert_within<const N: usize>(points: &[[f64; N]], recovered: &[[f64; N]], tol: f64) {
        for (i, (o, r)) in points.iter().zip(recovered.iter()).enumerate() {
            let d = distance(o, r);
            assert!(d <= tol, "point {i}: error {d} exceeds tolerance {tol}");
        }
    }

    #[test]
    fn round_trips_in_two_three_and_four_dimensions() {
        let mut rng = Rng::new(42);

        for tol in [1e-2, 1e-4, 1e-6] {
            let p2: Vec<[f64; 2]> = rng.points(2000, -50.0, 50.0);
            assert_within(&p2, &round_trip(&p2, tol), tol);

            let p3: Vec<[f64; 3]> = rng.points(2000, -50.0, 50.0);
            assert_within(&p3, &round_trip(&p3, tol), tol);

            let p4: Vec<[f64; 4]> = rng.points(2000, -50.0, 50.0);
            assert_within(&p4, &round_trip(&p4, tol), tol);
        }
    }

    /// The scale of set the previous implementation was tested at, kept so the comparison stays
    /// honest.
    #[test]
    fn round_trips_a_large_set() {
        let tol = 1e-4;
        let mut rng = Rng::new(42);
        let points: Vec<[f64; 3]> = rng.points(100_000, 0.0, 100.0);

        assert_within(&points, &round_trip(&points, tol), tol);
    }

    /// Axes with wildly different spans must get wildly different widths, which is the point of
    /// sizing each one independently.
    #[test]
    fn axes_are_sized_independently() {
        let mut rng = Rng::new(3);
        let points: Vec<[f64; 3]> = (0..1000)
            .map(|_| {
                [
                    rng.in_range(0.0, 1000.0),
                    rng.in_range(0.0, 1.0),
                    rng.in_range(0.0, 0.001),
                ]
            })
            .collect();

        let tol = 1e-4;
        assert_within(&points, &round_trip(&points, tol), tol);

        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();

        // The anchor widths, which describe the set as a whole, sit after the dimension, the
        // partition count and the anchor box.
        let offset = 1 + 4 + 8 * 2 * 3;
        let widths = &buf[offset..offset + 3];
        assert!(
            widths[0] > widths[1] && widths[1] > widths[2],
            "expected decreasing widths for decreasing spans, got {widths:?}"
        );
    }

    /// A flat axis costs nothing at all. The previous implementation charged two bytes per point
    /// for it.
    #[test]
    fn a_degenerate_axis_is_free() {
        let mut rng = Rng::new(9);
        let points: Vec<[f64; 3]> = (0..1000)
            .map(|_| [rng.in_range(-10.0, 10.0), rng.in_range(-10.0, 10.0), 4.25])
            .collect();

        let tol = 1e-3;
        let recovered = round_trip(&points, tol);
        assert_within(&points, &recovered, tol);

        // Exactly recovered, not merely within tolerance: there is only one value to return.
        for r in &recovered {
            assert_eq!(r[2], 4.25);
        }

        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();
        let offset = 1 + 4 + 8 * 2 * 3;
        assert_eq!(buf[offset + 2], 0, "flat axis should have a width of zero");
    }

    #[test]
    fn handles_a_single_point() {
        let points = [[1.5, -2.5, 3.5]];
        let recovered = round_trip(&points, 1e-6);

        // A single point defines a zero-extent box on every axis, so it comes back exactly.
        assert_eq!(recovered[0], points[0]);
    }

    #[test]
    fn handles_identical_points() {
        let points = vec![[7.0, 8.0, 9.0]; 500];
        let recovered = round_trip(&points, 1e-6);

        for r in &recovered {
            assert_eq!(*r, [7.0, 8.0, 9.0]);
        }

        // Every axis is degenerate, so there are no corner codes and no payload, and only the
        // headers remain.
        let mut buf = Vec::new();
        write_points(&mut buf, &points, 1e-6).unwrap();
        assert_eq!(buf.len(), 1 + 4 + 8 * 2 * 3 + 3 + 4 + 1 + 3);
    }

    #[test]
    fn handles_an_empty_set() {
        let points: Vec<[f64; 3]> = Vec::new();
        let recovered = round_trip(&points, 1e-6);

        assert!(recovered.is_empty());

        let mut buf = Vec::new();
        write_points(&mut buf, &points, 1e-6).unwrap();
        assert_eq!(
            buf.len(),
            1 + 4,
            "an empty block is a dimension and a partition count of zero"
        );
    }

    #[test]
    fn encoded_size_matches_the_width_arithmetic() {
        let mut rng = Rng::new(11);
        let count = 1000;
        let points: Vec<[f64; 3]> = rng.points(count, 0.0, 100.0);
        let tol = 1e-3;

        // Held to one partition so the arithmetic can be written out in full.
        let mut buf = Vec::new();
        write_points_with(&mut buf, &points, tol, Effort::Quick).unwrap();

        let anchor_at = 1 + 4 + 8 * 2 * 3;
        let header = anchor_at + 3 + 4 + 1 + 3;
        let widths: Vec<u32> = buf[header - 3..header]
            .iter()
            .map(|&b| u32::from(b))
            .collect();
        let anchor: Vec<u32> = buf[anchor_at..anchor_at + 3]
            .iter()
            .map(|&b| u32::from(b))
            .collect();

        // A lone partition's corners are the anchor's own, and they share the payload's bitstream.
        let corner_bits = 2 * anchor.iter().sum::<u32>();
        let payload_bits = count as u32 * widths.iter().sum::<u32>();
        let expected = header + (corner_bits + payload_bits).div_ceil(8) as usize;

        assert_eq!(buf.len(), expected);
    }

    /// The accounting a planner prices candidate partitions with has to match what the writer
    /// emits, or the plan is optimizing a cost the file does not have. Counts either side of a
    /// payload byte boundary are where a `div_ceil` gone the wrong way would show up.
    #[test]
    fn partition_bytes_matches_the_writer() {
        let mut rng = Rng::new(17);
        let tol = 1e-3;

        let cases = [
            (1usize, [100.0, 100.0, 100.0]),
            (7, [100.0, 1.0, 0.001]),
            (63, [1000.0, 10.0, 0.0]),
            (64, [1000.0, 10.0, 0.0]),
            (65, [1000.0, 10.0, 0.0]),
            (1000, [100.0, 100.0, 100.0]),
            (1001, [0.0, 0.0, 0.0]),
        ];

        for (count, spans) in cases {
            let points: Vec<[f64; 3]> = (0..count)
                .map(|_| std::array::from_fn(|i| rng.next_f64() * spans[i]))
                .collect();

            // Held to one partition, so the block is a header, an anchor and a single partition and
            // the arithmetic can be written out in full.
            let mut buf = Vec::new();
            write_points_with(&mut buf, &points, tol, Effort::Quick).unwrap();

            // The anchor widths sit after the dimension, the partition count and the anchor box.
            let anchor_at = 1 + 4 + 8 * 2 * 3;
            let anchor: [u8; 3] = buf[anchor_at..anchor_at + 3].try_into().unwrap();
            let widths: [u8; 3] = buf[anchor_at + 3 + 4 + 1..anchor_at + 3 + 4 + 4]
                .try_into()
                .unwrap();

            // Nothing chooses a frame at the moment, so every partition is written on world axes.
            assert_eq!(
                buf[anchor_at + 3 + 4],
                0,
                "the transform flag should be clear"
            );

            assert_eq!(
                buf.len() as u64,
                1 + 4 + 8 * 2 * 3 + 3 + partition_bytes(count, &widths, &anchor, None),
                "count {count}, spans {spans:?}, widths {widths:?}, anchor {anchor:?}"
            );
        }

        // An empty set has no anchor to record, so the block ends after a partition count of zero.
        let mut buf = Vec::new();
        let empty: Vec<[f64; 3]> = Vec::new();
        write_points(&mut buf, &empty, tol).unwrap();
        assert_eq!(buf.len(), 1 + 4);
    }

    /// The property the whole partitioning design rests on. A block may come back in several
    /// pieces, but it comes back as one sequence in the order it went in, because the encoder cuts
    /// the caller's order rather than rearranging it.
    ///
    /// Checked with points that are individually identifiable, so a partition written or read out
    /// of turn shows up as a mismatch rather than as a plausible-looking set.
    #[test]
    fn points_come_back_in_the_order_they_went_in() {
        // Tight clusters visited one after another, which is what makes the encoder cut at all.
        let mut rng = Rng::new(23);
        let mut points: Vec<[f64; 3]> = Vec::new();
        for k in 0..10 {
            let centre = [k as f64 * 500.0, 0.0, 0.0];
            for i in 0..600 {
                // The third coordinate counts, so every point says where it belongs.
                points.push([
                    centre[0] + rng.gaussian(0.5),
                    centre[1] + rng.gaussian(0.5),
                    (k * 600 + i) as f64 * 1e-3,
                ]);
            }
        }

        let tol = 1e-6;
        let recovered = round_trip(&points, tol);

        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();
        let partitions = u32::from_le_bytes(buf[1..5].try_into().unwrap());
        assert!(
            partitions > 1,
            "expected several partitions, got {partitions}"
        );

        assert_within(&points, &recovered, tol);
        for (i, r) in recovered.iter().enumerate() {
            let expected = i as f64 * 1e-3;
            assert!(
                (r[2] - expected).abs() <= tol,
                "point {i} came back carrying position {}, not {expected}",
                r[2]
            );
        }
    }

    /// Partitioning is an optimization, so it has to be one: the same data at [`Effort::Quick`] is
    /// the single partition the format wrote before, and it must not be the smaller of the two.
    #[test]
    fn partitioning_never_costs_bytes() {
        let mut rng = Rng::new(29);
        let cases: Vec<Vec<[f64; 3]>> = vec![
            // Clustered, which partitioning is for.
            (0..5000)
                .map(|i| {
                    let c = (i / 500) as f64 * 400.0;
                    [c + rng.gaussian(0.5), rng.gaussian(0.5), rng.gaussian(0.5)]
                })
                .collect(),
            // A chain, which it also helps.
            (0..5000)
                .map(|i| {
                    let t = i as f64 * 0.01;
                    [50.0 * t.cos(), 50.0 * t.sin(), 4.0 * t]
                })
                .collect(),
            // And scattered, where it should decline to cut and match Quick byte for byte.
            rng.points(5000, -50.0, 50.0),
        ];

        for (i, points) in cases.iter().enumerate() {
            let tol = 1e-3;
            let mut sizes = Vec::new();
            for effort in [Effort::Quick, Effort::Balanced, Effort::Thorough] {
                let mut buf = Vec::new();
                write_points_with(&mut buf, points, tol, effort).unwrap();
                assert_within(points, &read_points(&mut Cursor::new(&buf)).unwrap(), tol);
                sizes.push(buf.len());
            }

            assert!(
                sizes[0] >= sizes[1] && sizes[1] >= sizes[2],
                "case {i}: sizes went up with effort, {sizes:?}"
            );
        }
    }

    /// Re-encoding a block must come back to somewhere it has already been, rather than drifting.
    ///
    /// It is not a fixed point and cannot be. Where the partitions go is decided from where the
    /// points are, the quantizer moves them, and a set read back out therefore plans slightly
    /// differently than the one that went in. Some inputs settle on one file and some alternate
    /// between two.
    ///
    /// What matters is that the orbit is short and closed. A repeat means the decoded points
    /// repeated exactly, so rewriting a file over and over revisits a handful of representations
    /// instead of walking away from the geometry a tolerance at a time. Nothing about correctness
    /// rests on it, since every one of those files is within tolerance of its input; what rests on
    /// it is that error cannot accumulate through a pipeline that stores and restores.
    /// A partition count is attacker-controlled like any other length field, and it sits in front of
    /// a loop rather than an allocation, so the failure it must not produce is a hang.
    #[test]
    fn an_absurd_partition_count_does_not_allocate() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, u32::MAX).unwrap();

        let result: Result<Vec<[f64; 3]>> = read_points(&mut Cursor::new(&buf));
        assert!(result.is_err(), "should run out of input, not spin");
    }

    /// Against the previous whole-byte scheme, which charged 3 bytes per axis here.
    #[test]
    fn beats_whole_byte_widths() {
        let mut rng = Rng::new(13);
        let points: Vec<[f64; 3]> = rng.points(50_000, 0.0, 100.0);
        let tol = 1e-4;

        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();

        // The old scheme: 3 bytes per axis for this range and tolerance, 9 bytes per point.
        let old = points.len() * 9;
        let ratio = buf.len() as f64 / old as f64;
        assert!(
            ratio < 0.85,
            "expected a clear improvement on whole-byte widths, got {ratio}"
        );
    }

    #[test]
    fn reading_at_the_wrong_dimension_is_an_error() {
        let points: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0]];
        let mut buf = Vec::new();
        write_points(&mut buf, &points, 1e-6).unwrap();

        let mut cursor = Cursor::new(&buf);
        let wrong: Result<Vec<[f64; 2]>> = read_points(&mut cursor);
        assert!(matches!(wrong, Err(Error::Malformed(_))));
    }

    #[test]
    fn a_range_too_wide_for_the_tolerance_is_an_error() {
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut buf = Vec::new();

        assert!(matches!(
            write_points(&mut buf, &points, 1e-18),
            Err(Error::ToleranceNotRepresentable { .. })
        ));
    }

    #[test]
    fn non_finite_coordinates_are_rejected() {
        let points = [[0.0, 0.0, 0.0], [f64::NAN, 1.0, 1.0]];
        let mut buf = Vec::new();

        assert!(matches!(
            write_points(&mut buf, &points, 1e-3),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn truncated_input_is_an_error() {
        let mut rng = Rng::new(5);
        let points: Vec<[f64; 3]> = rng.points(100, 0.0, 10.0);
        let mut buf = Vec::new();
        write_points(&mut buf, &points, 1e-4).unwrap();

        for cut in [0, 1, 4, 20, buf.len() / 2, buf.len() - 1] {
            let mut cursor = Cursor::new(&buf[..cut]);
            let result: Result<Vec<[f64; 3]>> = read_points(&mut cursor);
            assert!(result.is_err(), "truncating to {cut} bytes should fail");
        }
    }

    /// A length field is attacker-controlled in any format. Claiming four billion points must not
    /// let the decoder allocate for them before a single one has arrived.
    #[test]
    fn an_absurd_length_field_does_not_allocate() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, 1).unwrap();
        for _ in 0..3 {
            write_f64(&mut buf, 0.0).unwrap();
        }
        for _ in 0..3 {
            write_f64(&mut buf, 1.0).unwrap();
        }
        for _ in 0..3 {
            write_u8(&mut buf, 16).unwrap();
        }
        write_u32(&mut buf, u32::MAX).unwrap();
        write_u8(&mut buf, 0).unwrap();
        for _ in 0..3 {
            write_u8(&mut buf, 16).unwrap();
        }

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[f64; 3]>> = read_points(&mut cursor);
        assert!(result.is_err(), "should run out of input, not memory");
    }

    /// Re-encoding a block must come back to somewhere it has already been, rather than drifting.
    ///
    /// It is not a fixed point and cannot be. Where the partitions go is decided from where the
    /// points are, the quantizer moves them, and a set read back out therefore plans slightly
    /// differently than the one that went in. Some inputs settle on one file and some alternate
    /// between two.
    ///
    /// What matters is that the orbit is short and closed. A repeat means the decoded points
    /// repeated exactly, so rewriting a file over and over revisits a handful of representations
    /// instead of walking away from the geometry a tolerance at a time. Nothing about correctness
    /// rests on it, since every one of those files is within tolerance of its input; what rests on
    /// it is that error cannot accumulate through a pipeline that stores and restores.
    #[test]
    fn re_encoding_reaches_a_cycle() {
        let mut rng = Rng::new(31);
        let cases: Vec<Vec<[f64; 3]>> = vec![
            (0..4000)
                .map(|i| {
                    let c = (i / 400) as f64 * 400.0;
                    [c + rng.gaussian(0.5), rng.gaussian(0.5), rng.gaussian(0.5)]
                })
                .collect(),
            (0..4000)
                .map(|i| {
                    let t = i as f64 * 0.01;
                    [50.0 * t.cos(), 50.0 * t.sin(), 4.0 * t]
                })
                .collect(),
            rng.points(4000, -50.0, 50.0),
        ];

        for (i, points) in cases.iter().enumerate() {
            let tol = 1e-3;
            let mut current = points.clone();
            let mut seen: Vec<Vec<u8>> = Vec::new();
            let mut closed = None;

            for round in 0..8 {
                let mut buf = Vec::new();
                write_points(&mut buf, &current, tol).unwrap();
                if seen.contains(&buf) {
                    closed = Some(round);
                    break;
                }
                seen.push(buf.clone());
                current = read_points(&mut Cursor::new(&buf)).unwrap();
            }

            // Measured across the corpus, the orbit closes by the fourth rewrite at the latest,
            // usually on a fixed point and otherwise on a pair.
            assert!(
                matches!(closed, Some(r) if r <= 4),
                "case {i}: orbit closed at {closed:?}"
            );
        }
    }

    /// Write a block of one partition in a frame chosen by the caller, which is what an encoder
    /// would do if anything chose frames.
    ///
    /// This is the only thing that exercises the transform path in either direction. Nothing in the
    /// crate emits a rotation on its own, so without a test that puts one in the stream by hand the
    /// reader's half would never run at all. See [`crate::transform`] for why it is idle.
    fn write_oriented<const N: usize>(
        points: &[[f64; N]],
        rotation: &Rotation<N>,
        tol: f64,
    ) -> Vec<u8> {
        let anchor = Bounds::from_points(points).unwrap();
        let anchor_widths = widths_for(&anchor, tol).unwrap();

        let centre = transform::centre_of(&anchor);
        let lattice = transform::rotated_anchor(&anchor, rotation);
        let local: Vec<[f64; N]> = points
            .iter()
            .map(|p| rotation.to_local(p, &centre))
            .collect();

        // The box has to be snapped onto the lattice its corners are coded on before the writer
        // sees it, which is the same thing `segment::finish_run` does for an ordinary run.
        let tight = Bounds::from_points(&local).unwrap();
        let (lo, hi) = segment::corner_codes(&tight, &lattice, &anchor_widths);
        let bounds = segment::decode_corners(&lo, &hi, &lattice, &anchor_widths);
        let widths = widths_for(&bounds, tol).unwrap();

        let mut buf = Vec::new();
        write_u8(&mut buf, N as u8).unwrap();
        write_u32(&mut buf, 1).unwrap();
        for i in 0..N {
            write_f64(&mut buf, anchor.mins[i]).unwrap();
        }
        for i in 0..N {
            write_f64(&mut buf, anchor.maxs[i]).unwrap();
        }
        for &w in anchor_widths.iter() {
            write_u8(&mut buf, w).unwrap();
        }
        write_partition(
            &mut buf,
            points,
            &bounds,
            &widths,
            &anchor,
            &anchor_widths,
            Some(rotation),
        )
        .unwrap();

        buf
    }

    /// A partition written in a frame has to come back within tolerance, in both dimensions that
    /// have a rotation encoding.
    #[test]
    fn a_transformed_partition_round_trips() {
        let mut rng = Rng::new(4_242);
        let tol = 1e-3;

        // An oblique slab, which is the shape a frame is for: thin on one axis of its own and thick
        // on all three of the world's.
        let frame = Rotation::<3>::from_axes([1.0, 1.0, 0.2], [-0.3, 0.1, 1.0]).unwrap();
        let centre = [0.0; 3];
        let solid: Vec<[f64; 3]> = (0..2_000)
            .map(|_| {
                let flat = [
                    rng.in_range(-20.0, 20.0),
                    rng.in_range(-15.0, 15.0),
                    rng.in_range(-0.05, 0.05),
                ];
                frame.to_world(&flat, &centre)
            })
            .collect();

        let buf = write_oriented(&solid, &frame, tol);
        assert_eq!(buf[1 + 4 + 8 * 2 * 3 + 3 + 4], 1, "the flag must be set");

        let mut cursor = Cursor::new(&buf);
        let back: Vec<[f64; 3]> = read_points(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "the rotation's bits have to be accounted for in the payload length"
        );
        assert_within(&solid, &back, tol);

        // And the frame has to have paid for itself, or there would be no reason to carry it.
        let mut plain = Vec::new();
        write_points_with(&mut plain, &solid, tol, Effort::Quick).unwrap();
        assert!(
            buf.len() < plain.len(),
            "oriented {} bytes against {} on world axes",
            buf.len(),
            plain.len()
        );

        let angled = Rotation::<2>::from_direction([2.0, 1.0]);
        let origin = [0.0; 2];
        let plane: Vec<[f64; 2]> = (0..2_000)
            .map(|_| {
                let flat = [rng.in_range(-20.0, 20.0), rng.in_range(-0.05, 0.05)];
                angled.to_world(&flat, &origin)
            })
            .collect();

        let buf = write_oriented(&plane, &angled, tol);
        let mut cursor = Cursor::new(&buf);
        let back: Vec<[f64; 2]> = read_points(&mut cursor).unwrap();
        assert_eq!(cursor.position() as usize, buf.len());
        assert_within(&plane, &back, tol);
    }

    /// A block header holding one partition, up to and including its transform flag, for tests that
    /// need to hand the reader something it must refuse.
    fn block_up_to_the_flag<const N: usize>(flag: u8) -> Vec<u8> {
        let mut buf = Vec::new();
        write_u8(&mut buf, N as u8).unwrap();
        write_u32(&mut buf, 1).unwrap();
        for _ in 0..N {
            write_f64(&mut buf, 0.0).unwrap();
        }
        for _ in 0..N {
            write_f64(&mut buf, 1.0).unwrap();
        }
        for _ in 0..N {
            write_u8(&mut buf, 16).unwrap();
        }
        write_u32(&mut buf, 1).unwrap();
        write_u8(&mut buf, flag).unwrap();
        buf
    }

    /// Only identity and rotation are defined. A flag this version does not know has to be refused,
    /// since carrying on regardless would put every point in the partition somewhere plausible and
    /// wrong, and the byte is reserved for a translation or a scale later.
    #[test]
    fn an_undefined_transform_is_refused_rather_than_ignored() {
        for flag in [2u8, 3, 17, 255] {
            let buf = block_up_to_the_flag::<3>(flag);
            let mut cursor = Cursor::new(&buf);
            let result: Result<Vec<[f64; 3]>> = read_points(&mut cursor);
            assert!(
                matches!(result, Err(Error::Malformed(_))),
                "flag {flag} was not refused"
            );
        }
    }

    /// Rotations are defined in two and three dimensions only, so a transform flag anywhere else
    /// names an encoding that does not exist and cannot be guessed at.
    #[test]
    fn a_transform_outside_two_and_three_dimensions_is_refused() {
        let buf = block_up_to_the_flag::<4>(1);
        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[f64; 4]>> = read_points(&mut cursor);
        assert!(matches!(result, Err(Error::Malformed(_))));

        // And the same block without the flag set reads as far as its payload, which is what says
        // the refusal above is about the transform rather than about the dimension.
        let buf = block_up_to_the_flag::<4>(0);
        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[f64; 4]>> = read_points(&mut cursor);
        assert!(matches!(
            result,
            Err(Error::Malformed(_)) | Err(Error::Io(_))
        ));
    }
}
