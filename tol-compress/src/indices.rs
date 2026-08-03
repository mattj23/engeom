//! Encoding and decoding blocks of simplex connectivity.
//!
//! A mesh stores three indices per triangle and roughly twice as many triangles as vertices, so
//! the index block is typically the larger part of a mesh file, not the coordinates. It is also
//! the more wasteful of the two under whole-byte widths: addressing 100,000 vertices needs 17
//! bits, which the previous scheme rounded up to 24.
//!
//! # Block layout
//!
//! ```text
//! u8       arity N, checked against the caller's expectation
//! u32      simplex count; a count of zero ends the block here
//! u8       bit width per index
//! bit-packed codes   N indices per simplex, in order, padded to a byte boundary
//! ```
//!
//! The width is recorded in the stream rather than recomputed by the reader from a vertex count it
//! was told separately. The extra byte buys the guarantee that an encoder and decoder cannot
//! silently disagree about how to parse the payload.

use crate::bits::{BitReader, BitWriter};
use crate::error::{Error, Result};
use crate::raw::{MAX_PREALLOC, read_u8, read_u32, write_u8, write_u32};
use std::io::{Read, Write};

/// The number of bits needed to address `count` distinct items, that is, to hold any value in
/// `0..count`.
///
/// Returns 0 when `count` is 0 or 1: a mesh with a single vertex needs no bits to say which vertex
/// a face refers to.
pub fn bits_for_count(count: u32) -> u8 {
    if count <= 1 {
        0
    } else {
        // The largest index is `count - 1`, so that is what has to fit.
        (32 - (count - 1).leading_zeros()) as u8
    }
}

/// Write a block of simplices, using the narrowest width that can address `index_limit` items.
///
/// `index_limit` is the number of items the indices refer into, normally the vertex count. Every
/// index must be below it.
///
/// # Errors
///
/// [`Error::Malformed`] if any index is at or above `index_limit`, or if there are more simplices
/// than a `u32` can count.
///
/// # Panics
///
/// Panics if `N` is zero or above 255.
pub fn write_indices<W: Write, const N: usize>(
    writer: &mut W,
    indices: &[[u32; N]],
    index_limit: u32,
) -> Result<()> {
    assert!((1..=255).contains(&N), "arity {N} is out of range");

    write_u8(writer, N as u8)?;

    let count = u32::try_from(indices.len())
        .map_err(|_| Error::Malformed("block holds more simplices than a u32 can count"))?;
    write_u32(writer, count)?;

    if indices.is_empty() {
        return Ok(());
    }

    // Checked here rather than left to the bit writer's width assertion, because an index past the
    // end of the vertex array is bad data rather than a programming error, and deserves an error
    // a caller can handle.
    for simplex in indices {
        for &i in simplex {
            if i >= index_limit {
                return Err(Error::Malformed(
                    "simplex refers to an index at or beyond the item count",
                ));
            }
        }
    }

    let bits = bits_for_count(index_limit);
    write_u8(writer, bits)?;

    let mut bw = BitWriter::new(&mut *writer);
    for simplex in indices {
        for &i in simplex {
            bw.write_bits(u64::from(i), bits)?;
        }
    }
    bw.finish()?;

    Ok(())
}

/// Read a block of simplices written by [`write_indices`].
///
/// `index_limit` is the number of items the indices are expected to refer into. Every decoded
/// index is checked against it, so a corrupt block fails here rather than further downstream where
/// an out-of-range index would index past the end of a vertex array.
///
/// # Errors
///
/// [`Error::Malformed`] if the block was written at a different arity than `N`, declares an
/// impossible width, or contains an index at or beyond `index_limit`. I/O errors propagate for
/// truncated input.
pub fn read_indices<R: Read, const N: usize>(
    reader: &mut R,
    index_limit: u32,
) -> Result<Vec<[u32; N]>> {
    assert!((1..=255).contains(&N), "arity {N} is out of range");

    let arity = read_u8(reader)?;
    if usize::from(arity) != N {
        return Err(Error::Malformed(
            "index block was written at a different arity than requested",
        ));
    }

    let count = read_u32(reader)? as usize;
    if count == 0 {
        return Ok(Vec::new());
    }

    let bits = read_u8(reader)?;
    if bits > 32 {
        return Err(Error::Malformed(
            "index block declares a width wider than a u32",
        ));
    }

    let mut out: Vec<[u32; N]> = Vec::with_capacity(count.min(MAX_PREALLOC));

    let mut br = BitReader::new(&mut *reader);
    for _ in 0..count {
        let mut simplex = [0u32; N];
        for slot in simplex.iter_mut() {
            let value = br.read_bits(bits)? as u32;
            if value >= index_limit {
                return Err(Error::Malformed(
                    "simplex refers to an index at or beyond the item count",
                ));
            }
            *slot = value;
        }
        out.push(simplex);
    }
    br.align();

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testgen::Rng;
    use std::io::Cursor;

    fn round_trip<const N: usize>(indices: &[[u32; N]], limit: u32) -> Vec<[u32; N]> {
        let mut buf = Vec::new();
        write_indices(&mut buf, indices, limit).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered: Vec<[u32; N]> = read_indices(&mut cursor, limit).unwrap();

        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder did not consume the whole block"
        );
        assert_eq!(recovered, indices, "indices must survive exactly");
        recovered
    }

    #[test]
    fn width_matches_a_hand_computed_table() {
        assert_eq!(bits_for_count(0), 0);
        assert_eq!(bits_for_count(1), 0);
        assert_eq!(bits_for_count(2), 1);
        assert_eq!(bits_for_count(3), 2);
        assert_eq!(bits_for_count(4), 2);
        assert_eq!(bits_for_count(5), 3);
        assert_eq!(bits_for_count(256), 8);
        assert_eq!(bits_for_count(257), 9);
        assert_eq!(bits_for_count(65_536), 16);
        assert_eq!(bits_for_count(65_537), 17);
        assert_eq!(bits_for_count(u32::MAX), 32);
    }

    /// The width has to hold `count - 1`, not `count`. Getting this wrong costs a whole bit at
    /// every power of two, which is exactly where meshes tend to sit.
    #[test]
    fn width_is_exact_at_every_power_of_two() {
        for exp in 1..32u32 {
            let count = 1u32 << exp;
            assert_eq!(bits_for_count(count), exp as u8, "count {count}");
            assert_eq!(
                bits_for_count(count + 1),
                exp as u8 + 1,
                "count {count} + 1"
            );
        }
    }

    #[test]
    fn round_trips_triangles() {
        let mut rng = Rng::new(21);
        let limit = 50_000u32;
        let faces: Vec<[u32; 3]> = (0..20_000)
            .map(|_| {
                [
                    (rng.next_u64() % limit as u64) as u32,
                    (rng.next_u64() % limit as u64) as u32,
                    (rng.next_u64() % limit as u64) as u32,
                ]
            })
            .collect();

        round_trip(&faces, limit);
    }

    #[test]
    fn round_trips_edges() {
        let mut rng = Rng::new(22);
        let limit = 1000u32;
        let edges: Vec<[u32; 2]> = (0..500)
            .map(|_| {
                [
                    (rng.next_u64() % limit as u64) as u32,
                    (rng.next_u64() % limit as u64) as u32,
                ]
            })
            .collect();

        round_trip(&edges, limit);
    }

    /// Either side of every byte boundary, where a whole-byte scheme would jump and this one
    /// should not.
    #[test]
    fn round_trips_across_width_boundaries() {
        for limit in [
            2u32, 3, 255, 256, 257, 1023, 1024, 1025, 65_535, 65_536, 65_537, 1_048_576,
        ] {
            let faces: Vec<[u32; 3]> = vec![
                [0, limit - 1, limit / 2],
                [limit - 1, 0, 1],
                [limit / 3, limit / 2, limit - 1],
            ];
            round_trip(&faces, limit);
        }
    }

    #[test]
    fn a_single_item_needs_no_bits_at_all() {
        let faces = vec![[0u32, 0, 0]; 100];
        round_trip(&faces, 1);

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 1).unwrap();

        // Arity, count, width, and nothing else: every index is known to be zero.
        assert_eq!(buf.len(), 1 + 4 + 1);
    }

    #[test]
    fn handles_an_empty_block() {
        let empty: Vec<[u32; 3]> = Vec::new();
        round_trip(&empty, 100);

        let mut buf = Vec::new();
        write_indices(&mut buf, &empty, 100).unwrap();
        assert_eq!(buf.len(), 1 + 4);
    }

    #[test]
    fn encoded_size_matches_the_width_arithmetic() {
        let limit = 100_000u32;
        let count = 5_000usize;
        let faces: Vec<[u32; 3]> = (0..count).map(|i| [i as u32, 0, 1]).collect();

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();

        let bits = u32::from(bits_for_count(limit));
        assert_eq!(bits, 17);

        let payload_bits = count as u32 * 3 * bits;
        assert_eq!(buf.len(), 1 + 4 + 1 + payload_bits.div_ceil(8) as usize);
    }

    /// The headline case for this block: a 100,000 vertex mesh needs 17 bits per index, where the
    /// previous scheme charged 3 whole bytes.
    #[test]
    fn beats_whole_byte_widths() {
        let limit = 100_000u32;
        let count = 200_000usize;
        let faces: Vec<[u32; 3]> = (0..count).map(|i| [(i % 100_000) as u32, 0, 1]).collect();

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();

        let old = count * 3 * 3;
        let ratio = buf.len() as f64 / old as f64;
        assert!(
            ratio < 0.72,
            "expected roughly 17/24 of the old size, got {ratio}"
        );
    }

    /// Small meshes suffered most under the old two-byte floor, which charged 16 bits per index
    /// regardless of how few vertices there were.
    #[test]
    fn small_meshes_gain_the_most() {
        let limit = 200u32;
        let count = 400usize;
        let faces: Vec<[u32; 3]> = (0..count).map(|i| [(i % 200) as u32, 0, 1]).collect();

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();

        // 8 bits per index against the old floor of 16.
        assert_eq!(bits_for_count(limit), 8);
        let old = count * 3 * 2;
        assert!((buf.len() as f64 / old as f64) < 0.55);
    }

    #[test]
    fn an_out_of_range_index_is_refused_on_write() {
        let faces = [[0u32, 1, 10]];
        let mut buf = Vec::new();

        assert!(matches!(
            write_indices(&mut buf, &faces, 10),
            Err(Error::Malformed(_))
        ));
    }

    /// A block whose indices point past the end of the vertex array must fail here, not by
    /// indexing out of bounds in whatever consumes it.
    #[test]
    fn an_out_of_range_index_is_refused_on_read() {
        let faces = [[0u32, 1, 999]];
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 1000).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[u32; 3]>> = read_indices(&mut cursor, 500);
        assert!(matches!(result, Err(Error::Malformed(_))));
    }

    #[test]
    fn reading_at_the_wrong_arity_is_an_error() {
        let faces = [[0u32, 1, 2]];
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 10).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[u32; 2]>> = read_indices(&mut cursor, 10);
        assert!(matches!(result, Err(Error::Malformed(_))));
    }

    #[test]
    fn truncated_input_is_an_error() {
        let faces: Vec<[u32; 3]> = (0..100).map(|i| [i, i + 1, i + 2]).collect();
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 1000).unwrap();

        for cut in [0, 1, 4, 6, buf.len() / 2, buf.len() - 1] {
            let mut cursor = Cursor::new(&buf[..cut]);
            let result: Result<Vec<[u32; 3]>> = read_indices(&mut cursor, 1000);
            assert!(result.is_err(), "truncating to {cut} bytes should fail");
        }
    }

    #[test]
    fn an_absurd_length_field_does_not_allocate() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, u32::MAX).unwrap();
        write_u8(&mut buf, 17).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[u32; 3]>> = read_indices(&mut cursor, 100_000);
        assert!(result.is_err(), "should run out of input, not memory");
    }

    #[test]
    fn an_impossible_width_is_refused() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, 1).unwrap();
        write_u8(&mut buf, 40).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[u32; 3]>> = read_indices(&mut cursor, 100);
        assert!(matches!(result, Err(Error::Malformed(_))));
    }
}
