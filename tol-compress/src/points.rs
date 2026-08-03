//! Encoding and decoding blocks of points in any number of dimensions.
//!
//! This is the core of the format. A block records where its points live and how finely each axis
//! was quantized, then packs the codes back to back with no padding between them.
//!
//! # Block layout
//!
//! ```text
//! u8       dimension N, checked against the caller's expectation
//! u32      partition count
//!   per partition:
//!     u32                point count; a count of zero ends the partition here
//!     u8                 transform flag, 0 for identity
//!     N x f64            lower corner
//!     N x f64            upper corner
//!     N x u8             bit width per axis
//!     bit-packed codes   N codes per point, in axis order, padded to a byte boundary
//! ```
//!
//! Two fields are placeholders for later work and are currently fixed. The partition count is
//! always 1, and the transform flag is always 0. Both are read back generically, so the decoder
//! written today will read the files a partitioning or reorienting encoder produces later.
//!
//! # Tolerance budget
//!
//! `tol` bounds the Euclidean distance between an original point and its recovered position. Each
//! axis is therefore given `tol / sqrt(N)`, so that `N` independent per-axis errors at their worst
//! combine in quadrature to exactly `tol`.

use crate::bits::{BitReader, BitWriter};
use crate::bounds::Bounds;
use crate::error::{Error, Result};
use crate::quantize::Quantizer;
use crate::raw::{read_f64, read_u8, read_u32, write_f64, write_u8, write_u32};
use std::io::{Read, Write};

/// Ceiling on how much a decoder will pre-allocate from a length field before seeing the data.
///
/// A corrupt or hostile stream can claim four billion points. Growing the vector as points
/// actually arrive costs a few reallocations on genuinely large blocks and bounds the damage on
/// bad ones.
const MAX_PREALLOC: usize = 1 << 16;

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
    assert!((1..=255).contains(&N), "dimension {N} is out of range");

    write_u8(writer, N as u8)?;
    // One partition for now. Spatial partitioning is a later increment, and the decoder already
    // loops, so adding it will not change the format.
    write_u32(writer, 1)?;

    write_partition(writer, points, tol)
}

fn write_partition<W: Write, const N: usize>(
    writer: &mut W,
    points: &[[f64; N]],
    tol: f64,
) -> Result<()> {
    let count = u32::try_from(points.len())
        .map_err(|_| Error::Malformed("partition holds more points than a u32 can index"))?;
    write_u32(writer, count)?;

    // An empty partition ends after its count. There is no meaningful bounding box to record and
    // inventing one would be a lie a decoder might believe.
    if points.is_empty() {
        return Ok(());
    }

    write_u8(writer, 0)?;

    let bounds = Bounds::from_points(points)?;
    for i in 0..N {
        write_f64(writer, bounds.mins[i])?;
    }
    for i in 0..N {
        write_f64(writer, bounds.maxs[i])?;
    }

    // Splitting the budget in quadrature: N axes each within tol/sqrt(N) put the point within tol.
    let axis_tol = tol / (N as f64).sqrt();
    let extents = bounds.extents();
    let mut quants = Vec::with_capacity(N);
    for (i, &extent) in extents.iter().enumerate() {
        quants.push(Quantizer::from_tol(bounds.mins[i], extent, axis_tol)?);
    }

    for q in &quants {
        write_u8(writer, q.bits())?;
    }

    let mut bw = BitWriter::new(&mut *writer);
    for p in points {
        for (i, q) in quants.iter().enumerate() {
            bw.write_bits(q.encode(p[i]), q.bits())?;
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
    let mut out = Vec::new();
    for _ in 0..partitions {
        read_partition(reader, &mut out)?;
    }

    Ok(out)
}

fn read_partition<R: Read, const N: usize>(reader: &mut R, out: &mut Vec<[f64; N]>) -> Result<()> {
    let count = read_u32(reader)? as usize;
    if count == 0 {
        return Ok(());
    }

    let transform = read_u8(reader)?;
    if transform != 0 {
        return Err(Error::Malformed(
            "point block carries a transform, which this version cannot apply",
        ));
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
            return Err(Error::Malformed("point block has invalid bounds"));
        }
    }

    let mut quants = Vec::with_capacity(N);
    for i in 0..N {
        let bits = read_u8(reader)?;
        if bits > 53 {
            return Err(Error::Malformed(
                "point block declares an impossible bit width",
            ));
        }
        quants.push(Quantizer::new(mins[i], maxs[i] - mins[i], bits));
    }

    out.reserve(count.min(MAX_PREALLOC));

    let mut br = BitReader::new(&mut *reader);
    for _ in 0..count {
        let mut p = [0.0f64; N];
        for (i, q) in quants.iter().enumerate() {
            p[i] = q.decode(br.read_bits(q.bits())?);
        }
        out.push(p);
    }
    br.align();

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

        // Widths sit after the dimension, partition count, point count, transform flag and bounds.
        let offset = 1 + 4 + 4 + 1 + 8 * 2 * 3;
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
        let offset = 1 + 4 + 4 + 1 + 8 * 2 * 3;
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

        // Every axis is degenerate, so the payload is empty and only the header remains.
        let mut buf = Vec::new();
        write_points(&mut buf, &points, 1e-6).unwrap();
        assert_eq!(buf.len(), 1 + 4 + 4 + 1 + 8 * 2 * 3 + 3);
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
            1 + 4 + 4,
            "an empty block should carry no bounds"
        );
    }

    #[test]
    fn encoded_size_matches_the_width_arithmetic() {
        let mut rng = Rng::new(11);
        let count = 1000;
        let points: Vec<[f64; 3]> = rng.points(count, 0.0, 100.0);
        let tol = 1e-3;

        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();

        let header = 1 + 4 + 4 + 1 + 8 * 2 * 3 + 3;
        let widths: Vec<u32> = buf[header - 3..header]
            .iter()
            .map(|&b| u32::from(b))
            .collect();
        let payload_bits = count as u32 * widths.iter().sum::<u32>();
        let expected = header + payload_bits.div_ceil(8) as usize;

        assert_eq!(buf.len(), expected);
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
        write_u32(&mut buf, u32::MAX).unwrap();
        write_u8(&mut buf, 0).unwrap();
        for _ in 0..3 {
            write_f64(&mut buf, 0.0).unwrap();
        }
        for _ in 0..3 {
            write_f64(&mut buf, 1.0).unwrap();
        }
        for _ in 0..3 {
            write_u8(&mut buf, 16).unwrap();
        }

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[f64; 3]>> = read_points(&mut cursor);
        assert!(result.is_err(), "should run out of input, not memory");
    }

    #[test]
    fn a_declared_transform_is_refused_rather_than_ignored() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, 1).unwrap();
        write_u32(&mut buf, 1).unwrap();
        write_u8(&mut buf, 1).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[f64; 3]>> = read_points(&mut cursor);
        assert!(matches!(result, Err(Error::Malformed(_))));
    }
}
