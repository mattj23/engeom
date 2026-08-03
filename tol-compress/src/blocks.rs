//! Per-block adaptive bit widths, the PFor / SIMD-BP128 technique without the SIMD.
//!
//! A block of fixed-width codes has to be as wide as its widest value. Chopping the stream into
//! short blocks and giving each one its own width means a single large value inflates only the 64
//! values around it instead of the whole stream. It needs no probability model and no entropy
//! coder, so it stays inside this crate's format-not-codec decision, and it keeps random access at
//! block granularity.
//!
//! # What this does and does not buy
//!
//! It only pays on data whose local range is much smaller than its global range, which means data
//! that has already been through a predictor. Measured on real meshes:
//!
//! - On raw quantized coordinates it is worth **0 to 2.5%**, because the bounding box is by
//!   construction the tightest one, so nearly every block spans nearly the full range.
//! - On high-water-mark index deltas it is worth about 7% beyond what a single flat width
//!   over the whole stream would give (engine blade, -26.7% flat against -33.8% blocked).
//!
//! So this is infrastructure for [`crate::indices`] mode 1 and for the scalar and label attribute
//! streams that come later, not something to use on the points block.
//!
//! # Stream layout
//!
//! ```text
//! u32    payload byte length
//!        per block of up to 64 values:
//!          6 bits         width w, 0 to 63
//!          count x w      the block's values, LSB-first as everywhere else
//!        padded to a byte boundary
//! ```
//!
//! The number of blocks is implied by the value count, which the caller already knows from the
//! header of whatever block this stream belongs to; only the final block may be short. A width of
//! 0 stays legal and free, so a run of 64 equal values costs six bits.
//!
//! The byte length is redundant with the widths in the stream, and it earns its four bytes twice
//! over: it lets the payload be read out in one piece before decoding, which is what
//! [`crate::bits::BitReader`] needs, and it lets a reader skip the stream without parsing it. A
//! length that disagrees with the blocks is caught rather than ignored.

use crate::bits::{BitReader, BitWriter, read_payload};
use crate::error::{Error, Result};
use crate::raw::{MAX_PREALLOC, read_u32, write_u32};
use std::io::{Read, Write};

/// Values per block.
///
/// Larger blocks amortize the width header better and smaller ones adapt faster. Measured on index
/// deltas the curve is shallow and monotonic (16 gives -32.0%, 64 gives -33.8%, 256 gives -34.2%),
/// so this sits where the remaining tenths stop being worth the coarser random-access granularity.
/// It is part of the format and cannot change without a version bump.
pub const BLOCK: usize = 64;

/// Bits used by a block's width header.
const WIDTH_BITS: u8 = 6;

/// The widest code this scheme can express, set by the six-bit width header.
///
/// Every current and planned use fits: indices are at most 32 bits, and quantized coordinates and
/// scalars at most 53.
pub const MAX_WIDTH: u8 = 63;

/// The narrowest width that can hold every value in a block.
///
/// # Errors
///
/// [`Error::Malformed`] if some value needs all 64 bits, which this scheme cannot express.
fn block_width(block: &[u64]) -> Result<u8> {
    let max = block.iter().copied().max().unwrap_or(0);
    let width = (64 - max.leading_zeros()) as u8;

    if width > MAX_WIDTH {
        return Err(Error::Malformed(
            "block holds a value needing all 64 bits, which block coding cannot express",
        ));
    }

    Ok(width)
}

/// The exact number of bits [`write_stream`] emits, excluding the length prefix and the final
/// padding to a byte boundary.
///
/// Useful for deciding whether block coding is worth using for a given stream at all, which is a
/// choice the caller makes and records.
///
/// # Errors
///
/// [`Error::Malformed`] if any value needs all 64 bits.
pub fn encoded_bits(values: &[u64]) -> Result<u64> {
    let mut total = 0u64;

    for chunk in values.chunks(BLOCK) {
        let width = block_width(chunk)?;
        total += u64::from(WIDTH_BITS) + u64::from(width) * chunk.len() as u64;
    }

    Ok(total)
}

/// Write values as a length-prefixed, block-coded stream.
///
/// # Errors
///
/// [`Error::Malformed`] if any value needs all 64 bits, or if the encoded payload is longer than a
/// `u32` can measure.
pub fn write_stream<W: Write>(writer: &mut W, values: &[u64]) -> Result<()> {
    let bytes = encoded_bits(values)?.div_ceil(8);
    let length = u32::try_from(bytes)
        .map_err(|_| Error::Malformed("block-coded stream is longer than a u32 can measure"))?;
    write_u32(writer, length)?;

    let mut bw = BitWriter::new(&mut *writer);
    for chunk in values.chunks(BLOCK) {
        // Recomputed rather than carried over from `encoded_bits`, so the width the header declares
        // and the width the values are packed at cannot come from different passes over the data.
        let width = block_width(chunk)?;
        bw.write_bits(u64::from(width), WIDTH_BITS)?;

        for &value in chunk {
            bw.write_bits(value, width)?;
        }
    }
    bw.finish()?;

    Ok(())
}

/// Read a stream written by [`write_stream`].
///
/// `count` is how many values the stream holds, which the caller knows from the header of the block
/// this stream belongs to.
///
/// # Errors
///
/// [`Error::Malformed`] if the declared length disagrees with the blocks it contains. I/O errors
/// propagate for truncated input.
pub fn read_stream<R: Read>(reader: &mut R, count: usize) -> Result<Vec<u64>> {
    let length = read_u32(reader)? as usize;
    let payload = read_payload(reader, length)?;

    let mut out = Vec::with_capacity(count.min(MAX_PREALLOC));
    let mut br = BitReader::new(&payload);

    let mut left = count;
    while left > 0 {
        let n = left.min(BLOCK);

        // Six bits cannot express a width above 63, so no range check is possible or needed here.
        let width = br.read_bits(WIDTH_BITS)? as u8;
        for _ in 0..n {
            out.push(br.read_bits(width)?);
        }

        left -= n;
    }

    // Catches a length prefix that does not match the blocks, which would otherwise leave the rest
    // of the containing stream misaligned.
    br.finish()?;

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testgen::Rng;
    use std::io::Cursor;

    fn round_trip(values: &[u64]) -> Vec<u64> {
        let mut buf = Vec::new();
        write_stream(&mut buf, values).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_stream(&mut cursor, values.len()).unwrap();

        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder did not consume the whole stream"
        );
        assert_eq!(back, values, "block coding must be exact");
        back
    }

    /// Every width the scheme can express, at the values that sit on its boundaries.
    #[test]
    fn round_trips_at_every_width() {
        for width in 0..=MAX_WIDTH {
            let values: Vec<u64> = if width == 0 {
                vec![0; 100]
            } else {
                let max = (1u64 << width) - 1;
                vec![0, 1, max / 2, max, max, 1, 0]
            };

            round_trip(&values);
        }
    }

    /// Block boundaries are where an off-by-one in the chunking would hide.
    #[test]
    fn round_trips_around_every_block_boundary() {
        let mut rng = Rng::new(4_001);
        for count in [0usize, 1, 2, 63, 64, 65, 127, 128, 129, 255, 256, 257, 1000] {
            let values: Vec<u64> = (0..count).map(|_| rng.next_u64() & 0xFFFF).collect();
            round_trip(&values);
        }
    }

    #[test]
    fn an_empty_stream_is_just_a_length() {
        let mut buf = Vec::new();
        write_stream(&mut buf, &[]).unwrap();

        assert_eq!(buf, vec![0, 0, 0, 0]);
        assert!(read_stream(&mut Cursor::new(&buf), 0).unwrap().is_empty());
    }

    /// The whole point of the scheme: blocks of very different magnitude each get their own width
    /// rather than every block paying for the widest.
    #[test]
    fn blocks_are_sized_independently() {
        let mut values = Vec::new();
        values.extend(std::iter::repeat_n(3u64, BLOCK)); // 2 bits
        values.extend(std::iter::repeat_n(1_000_000u64, BLOCK)); // 20 bits
        values.extend(std::iter::repeat_n(0u64, BLOCK)); // 0 bits

        round_trip(&values);

        // Three width headers, then 2 bits, 20 bits and nothing at all per value.
        let expected_bits = (6 + 2 * BLOCK) + (6 + 20 * BLOCK) + 6;
        assert_eq!(encoded_bits(&values).unwrap(), expected_bits as u64);

        // A single width over the whole stream would charge 20 bits for all 192 values.
        let flat = 20 * 3 * BLOCK;
        assert!(
            (expected_bits as f64) < 0.6 * flat as f64,
            "expected a clear win over a flat width, got {expected_bits} against {flat}"
        );
    }

    /// One large value must inflate its own block and leave its neighbours alone. This is what a
    /// single flat width over the stream cannot do, and it is why the tail of a residual
    /// distribution is survivable.
    #[test]
    fn an_outlier_inflates_only_its_own_block() {
        let mut values = vec![1u64; BLOCK * 3];
        values[BLOCK + 10] = u64::from(u32::MAX);

        round_trip(&values);

        // Two blocks at 1 bit, one at 32.
        let expected = (6 + BLOCK) + (6 + 32 * BLOCK) + (6 + BLOCK);
        assert_eq!(encoded_bits(&values).unwrap(), expected as u64);
    }

    /// A run of identical values costs only its width headers, which is what makes degenerate data
    /// nearly free rather than merely cheap.
    #[test]
    fn a_constant_run_costs_only_its_headers() {
        let values = vec![0u64; BLOCK * 10];

        round_trip(&values);
        assert_eq!(encoded_bits(&values).unwrap(), 6 * 10);

        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();
        assert_eq!(buf.len(), 4 + (60f64 / 8.0).ceil() as usize);
    }

    #[test]
    fn encoded_size_matches_the_block_arithmetic() {
        let mut rng = Rng::new(4_002);
        let values: Vec<u64> = (0..1000).map(|_| rng.next_u64() & 0x3FF).collect();

        let bits = encoded_bits(&values).unwrap();
        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();

        assert_eq!(buf.len() as u64, 4 + bits.div_ceil(8));
    }

    /// The scheme tops out at 63 bits because the width header is six bits wide. A value needing
    /// all 64 has to be refused rather than silently truncated.
    #[test]
    fn a_value_needing_all_64_bits_is_refused() {
        let values = vec![1u64, u64::MAX];

        assert!(matches!(encoded_bits(&values), Err(Error::Malformed(_))));

        let mut buf = Vec::new();
        assert!(matches!(
            write_stream(&mut buf, &values),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn the_widest_expressible_value_is_accepted() {
        let values = vec![(1u64 << MAX_WIDTH) - 1, 0];
        round_trip(&values);
    }

    #[test]
    fn truncated_input_is_an_error() {
        let mut rng = Rng::new(4_003);
        let values: Vec<u64> = (0..500).map(|_| rng.next_u64() & 0xFFFF).collect();

        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();

        for cut in [0, 1, 3, 4, 5, 20, buf.len() / 2, buf.len() - 1] {
            assert!(
                read_stream(&mut Cursor::new(&buf[..cut]), values.len()).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
    }

    /// A count larger than the payload holds must fail, not return whatever the padding decodes to.
    #[test]
    fn a_count_beyond_the_payload_is_an_error() {
        let values = vec![7u64; 10];
        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();

        assert!(read_stream(&mut Cursor::new(&buf), 1000).is_err());
    }

    /// A count smaller than the payload holds leaves bytes unread, which would desynchronize
    /// everything after this stream.
    #[test]
    fn a_count_short_of_the_payload_is_an_error() {
        let values: Vec<u64> = (0..200u64).collect();
        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();

        assert!(matches!(
            read_stream(&mut Cursor::new(&buf), 64),
            Err(Error::Malformed(_))
        ));
    }

    /// A length prefix that disagrees with the blocks must be caught rather than trusted.
    #[test]
    fn a_corrupt_length_prefix_is_refused() {
        let values = vec![5u64; 100];
        let mut buf = Vec::new();
        write_stream(&mut buf, &values).unwrap();

        let real = u32::from_le_bytes(buf[..4].try_into().unwrap());
        buf[..4].copy_from_slice(&(real + 4).to_le_bytes());

        assert!(read_stream(&mut Cursor::new(&buf), 100).is_err());
    }

    /// A length field is attacker-controlled in any format.
    #[test]
    fn an_absurd_length_does_not_allocate() {
        let mut buf = Vec::new();
        write_u32(&mut buf, u32::MAX).unwrap();
        buf.extend_from_slice(&[0u8; 16]);

        assert!(
            read_stream(&mut Cursor::new(&buf), 10).is_err(),
            "should run out of input, not memory"
        );
    }
}
