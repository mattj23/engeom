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
//! u32      simplex count; a count of zero ends the block here, before the mode byte
//! u8       coding mode: 0 absolute, 1 high-water-mark delta
//!   mode 0:  u8 bit width per index, then count x N packed codes
//!   mode 1:  a block-coded stream of count x N deltas, see [`crate::blocks`]
//! ```
//!
//! In mode 0 the width is recorded in the stream rather than recomputed by the reader from a vertex
//! count it was told separately. The extra byte buys the guarantee that an encoder and decoder
//! cannot silently disagree about how to parse the payload.
//!
//! # High-water-mark coding
//!
//! When the simplices are in **first-use order**, meaning every index is at most the number of
//! vertices introduced before it, a corner can be coded as the distance back from that running
//! high-water mark rather than as an absolute index. A delta of zero means "the next new vertex"
//! and advances the mark. On a mesh whose faces have been ordered so neighbours follow one another,
//! those distances are small and tightly clustered, and the block coder packs them at a fraction of
//! the width an absolute index needs. Measured on real meshes at -25% to -50% against mode 0.
//!
//! [`crate::reorder`] produces that order. This module never reorders anything itself, because the
//! reordering has to move the points and every per-vertex attribute in lockstep and only the caller
//! owns those.
//!
//! # Choosing a mode
//!
//! [`write_indices`] tries mode 1 whenever the input happens to be in first-use order and keeps it
//! only when it actually comes out smaller, so **no input encodes larger than it used to**. Small
//! meshes are where this matters: the block coder's per-block width headers and length prefix are
//! not free, and on a few hundred triangles they can outweigh what the deltas save.
//!
//! An out-of-range index is unrepresentable in mode 1 rather than merely rejected, since a delta can
//! only reach back to vertices that have already been introduced.

use crate::bits::{BitReader, BitWriter, read_payload};
use crate::blocks;
use crate::error::{Error, Result};
use crate::raw::{MAX_PREALLOC, read_u8, read_u32, write_u8, write_u32};
use std::io::{Read, Write};

/// Indices stored as absolute values at one fixed width.
pub const MODE_ABSOLUTE: u8 = 0;

/// Indices stored as block-coded distances back from a running high-water mark.
pub const MODE_HIGH_WATER: u8 = 1;

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

/// Distances back from a running high-water mark, if `indices` is in first-use order.
///
/// Returns `None` the moment it is not, which is the ordinary outcome for a caller who kept their
/// own vertex order.
fn high_water_deltas<const N: usize>(indices: &[[u32; N]]) -> Option<Vec<u64>> {
    let mut deltas = Vec::with_capacity(indices.len() * N);
    let mut high_water = 0u32;

    for simplex in indices {
        for &index in simplex {
            if index > high_water {
                return None;
            }

            deltas.push(u64::from(high_water - index));
            if index == high_water {
                high_water += 1;
            }
        }
    }

    Some(deltas)
}

/// Write a block of simplices, choosing whichever coding comes out smaller.
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

    let width = bits_for_count(index_limit);
    let absolute_bytes = (indices.len() as u64 * N as u64 * u64::from(width)).div_ceil(8);

    // Only worth keeping if it actually wins, so a mesh too small for the block coder's overheads
    // falls back rather than growing.
    if let Some(deltas) = high_water_deltas(indices) {
        let coded_bytes = 4 + blocks::encoded_bits(&deltas)?.div_ceil(8);

        if coded_bytes < absolute_bytes {
            write_u8(writer, MODE_HIGH_WATER)?;
            return blocks::write_stream(writer, &deltas);
        }
    }

    write_u8(writer, MODE_ABSOLUTE)?;
    write_u8(writer, width)?;

    let mut bw = BitWriter::new(&mut *writer);
    for simplex in indices {
        for &i in simplex {
            bw.write_bits(u64::from(i), width)?;
        }
    }
    bw.finish()?;

    Ok(())
}

/// Read a block of simplices written by [`write_indices`], in either coding mode.
///
/// `index_limit` is the number of items the indices are expected to refer into. Every decoded
/// index is checked against it, so a corrupt block fails here rather than further downstream where
/// an out-of-range index would index past the end of a vertex array.
///
/// # Errors
///
/// [`Error::Malformed`] if the block was written at a different arity than `N`, declares a coding
/// mode or width this version does not know, or contains an index at or beyond `index_limit`. I/O
/// errors propagate for truncated input.
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

    match read_u8(reader)? {
        MODE_ABSOLUTE => read_absolute(reader, count, index_limit),
        MODE_HIGH_WATER => read_high_water(reader, count, index_limit),
        _ => Err(Error::Malformed(
            "index block declares a coding mode this version does not know",
        )),
    }
}

fn read_absolute<R: Read, const N: usize>(
    reader: &mut R,
    count: usize,
    index_limit: u32,
) -> Result<Vec<[u32; N]>> {
    let width = read_u8(reader)?;
    if width > 32 {
        return Err(Error::Malformed(
            "index block declares a width wider than a u32",
        ));
    }

    let mut out: Vec<[u32; N]> = Vec::with_capacity(count.min(MAX_PREALLOC));

    // The reader must not read past the payload, because a byte-oriented header follows it. The
    // width, arity and count fully determine how long it is.
    let payload = (count as u64 * N as u64 * u64::from(width)).div_ceil(8) as usize;

    let bytes = read_payload(reader, payload)?;
    let mut br = BitReader::new(&bytes);
    for _ in 0..count {
        let mut simplex = [0u32; N];
        for slot in simplex.iter_mut() {
            let value = br.read_bits(width)? as u32;
            if value >= index_limit {
                return Err(Error::Malformed(
                    "simplex refers to an index at or beyond the item count",
                ));
            }
            *slot = value;
        }
        out.push(simplex);
    }
    br.finish()?;

    Ok(out)
}

fn read_high_water<R: Read, const N: usize>(
    reader: &mut R,
    count: usize,
    index_limit: u32,
) -> Result<Vec<[u32; N]>> {
    let total = count
        .checked_mul(N)
        .ok_or(Error::Malformed("index block declares an impossible size"))?;

    let deltas = blocks::read_stream(reader, total)?;

    let mut out: Vec<[u32; N]> = Vec::with_capacity(count.min(MAX_PREALLOC));
    let mut high_water = 0u32;
    let mut deltas = deltas.into_iter();

    for _ in 0..count {
        let mut simplex = [0u32; N];

        for slot in simplex.iter_mut() {
            let delta = deltas
                .next()
                .ok_or(Error::Malformed("index block is short of deltas"))?;

            // Reaching back further than the mark would name a vertex that has not been introduced,
            // which no encoder can produce and no decoder can resolve.
            if delta > u64::from(high_water) {
                return Err(Error::Malformed(
                    "index delta reaches back past the high water mark",
                ));
            }

            let value = high_water - delta as u32;
            if value >= index_limit {
                return Err(Error::Malformed(
                    "simplex refers to an index at or beyond the item count",
                ));
            }

            if delta == 0 {
                high_water += 1;
            }
            *slot = value;
        }

        out.push(simplex);
    }

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

    /// Which coding a block came out in, read straight from the mode byte.
    fn mode_of(buf: &[u8]) -> u8 {
        buf[5]
    }

    /// Faces in first-use order over `vertex_count` vertices, as a reorderer would produce.
    fn first_use_faces(vertex_count: u32) -> Vec<[u32; 3]> {
        let mut faces = Vec::new();
        for v in 2..vertex_count {
            faces.push([v - 2, v - 1, v]);
        }
        faces
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

        // Arity, count, mode, width, and nothing else: every index is known to be zero.
        assert_eq!(buf.len(), 1 + 4 + 1 + 1);
        assert_eq!(mode_of(&buf), MODE_ABSOLUTE);
    }

    #[test]
    fn handles_an_empty_block() {
        let empty: Vec<[u32; 3]> = Vec::new();
        round_trip(&empty, 100);

        let mut buf = Vec::new();
        write_indices(&mut buf, &empty, 100).unwrap();

        // An empty block ends before the mode byte, so container framing does not move.
        assert_eq!(buf.len(), 1 + 4);
    }

    #[test]
    fn encoded_size_matches_the_width_arithmetic() {
        let limit = 100_000u32;
        let count = 5_000usize;
        // Descending, so the very first corner is past the high water mark and mode 1 is not
        // available. Ascending indices would quietly qualify for it and measure something else.
        let faces: Vec<[u32; 3]> = (0..count).map(|i| [(count - i) as u32, 0, 1]).collect();

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();
        assert_eq!(mode_of(&buf), MODE_ABSOLUTE, "not in first-use order");

        let bits = u32::from(bits_for_count(limit));
        assert_eq!(bits, 17);

        let payload_bits = count as u32 * 3 * bits;
        assert_eq!(buf.len(), 1 + 4 + 1 + 1 + payload_bits.div_ceil(8) as usize);
    }

    /// The headline case for absolute coding: a 100,000 vertex mesh needs 17 bits per index, where
    /// the previous scheme charged 3 whole bytes.
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

    // ---------------------------------------------------------------------------------------
    // High-water-mark coding
    // ---------------------------------------------------------------------------------------

    #[test]
    fn first_use_order_takes_the_high_water_mode() {
        let faces = first_use_faces(20_000);
        let limit = 20_000u32;

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();

        assert_eq!(mode_of(&buf), MODE_HIGH_WATER);
        round_trip(&faces, limit);
    }

    /// The point of the mode: on a mesh whose corners stay near the high water mark, the deltas
    /// pack far narrower than an absolute index can.
    #[test]
    fn high_water_coding_is_much_smaller() {
        let limit = 20_000u32;
        let faces = first_use_faces(limit);

        let mut coded = Vec::new();
        write_indices(&mut coded, &faces, limit).unwrap();
        assert_eq!(mode_of(&coded), MODE_HIGH_WATER);

        // What mode 0 would have cost for the same block.
        let absolute = 1
            + 4
            + 1
            + 1
            + (faces.len() as u64 * 3 * u64::from(bits_for_count(limit))).div_ceil(8) as usize;

        let ratio = coded.len() as f64 / absolute as f64;
        assert!(
            ratio < 0.55,
            "expected a large gain from delta coding, got {ratio}"
        );
    }

    /// Not being in first-use order is not an error, it just means the block falls back.
    #[test]
    fn out_of_order_indices_fall_back_to_absolute() {
        // The very first corner names a vertex that has not been introduced.
        let faces = vec![[5u32, 0, 1], [0, 1, 2]];

        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 10).unwrap();

        assert_eq!(mode_of(&buf), MODE_ABSOLUTE);
        round_trip(&faces, 10);
    }

    /// The mode is chosen on measured size, so a block too small for the block coder's overheads
    /// keeps the plain coding rather than growing.
    #[test]
    fn a_tiny_block_does_not_grow() {
        let faces = first_use_faces(6);
        let limit = 6u32;

        let mut chosen = Vec::new();
        write_indices(&mut chosen, &faces, limit).unwrap();

        let absolute = 1
            + 4
            + 1
            + 1
            + (faces.len() as u64 * 3 * u64::from(bits_for_count(limit))).div_ceil(8) as usize;

        assert!(
            chosen.len() <= absolute,
            "picked a {} byte encoding over a {} byte one",
            chosen.len(),
            absolute
        );
        round_trip(&faces, limit);
    }

    /// Every corpus mesh, reordered, must round-trip its connectivity exactly and never come out
    /// larger than absolute coding would have made it.
    #[test]
    fn no_reordered_corpus_mesh_encodes_larger() {
        for case in crate::corpus::all() {
            if case.faces.is_empty() {
                continue;
            }

            let limit = case.points.len() as u32;
            let r = crate::reorder::optimize(&case.faces, case.points.len()).unwrap();

            let mut buf = Vec::new();
            write_indices(&mut buf, &r.faces, limit).unwrap();

            let absolute = 1
                + 4
                + 1
                + 1
                + (r.faces.len() as u64 * 3 * u64::from(bits_for_count(limit))).div_ceil(8)
                    as usize;

            assert!(
                buf.len() <= absolute,
                "{}: {} bytes against {} for absolute coding",
                case.name,
                buf.len(),
                absolute
            );

            let back: Vec<[u32; 3]> = read_indices(&mut Cursor::new(&buf), limit).unwrap();
            assert_eq!(back, r.faces, "{}", case.name);
        }
    }

    /// A structurally valid mode 1 block whose first corner reaches back before any vertex has
    /// been introduced. This is the check that makes an out-of-range index unrepresentable rather
    /// than merely rejected.
    #[test]
    fn a_delta_past_the_high_water_mark_is_refused() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, 1).unwrap();
        write_u8(&mut buf, MODE_HIGH_WATER).unwrap();
        blocks::write_stream(&mut buf, &[5, 0, 0]).unwrap();

        assert!(matches!(
            read_indices::<_, 3>(&mut Cursor::new(&buf), 100),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn an_unknown_coding_mode_is_refused() {
        let faces = [[0u32, 1, 2]];
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 3).unwrap();

        buf[5] = 7;

        assert!(matches!(
            read_indices::<_, 3>(&mut Cursor::new(&buf), 3),
            Err(Error::Malformed(_))
        ));
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

    /// The same, for a block that chose high-water coding: a mesh decoded against a shorter vertex
    /// array must be refused rather than resolving to a vertex that is not there.
    #[test]
    fn a_high_water_block_is_checked_against_the_vertex_count() {
        let faces = first_use_faces(2000);
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 2000).unwrap();
        assert_eq!(mode_of(&buf), MODE_HIGH_WATER);

        let result: Result<Vec<[u32; 3]>> = read_indices(&mut Cursor::new(&buf), 100);
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

    /// The same, for a block in high-water coding, whose payload is framed differently.
    #[test]
    fn a_truncated_high_water_block_is_an_error() {
        let faces = first_use_faces(2000);
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, 2000).unwrap();
        assert_eq!(mode_of(&buf), MODE_HIGH_WATER);

        for cut in [6, 8, 10, 40, buf.len() / 2, buf.len() - 1] {
            let result: Result<Vec<[u32; 3]>> = read_indices(&mut Cursor::new(&buf[..cut]), 2000);
            assert!(result.is_err(), "truncating to {cut} bytes should fail");
        }
    }

    #[test]
    fn an_absurd_length_field_does_not_allocate() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 3).unwrap();
        write_u32(&mut buf, u32::MAX).unwrap();
        write_u8(&mut buf, MODE_ABSOLUTE).unwrap();
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
        write_u8(&mut buf, MODE_ABSOLUTE).unwrap();
        write_u8(&mut buf, 40).unwrap();

        let mut cursor = Cursor::new(&buf);
        let result: Result<Vec<[u32; 3]>> = read_indices(&mut cursor, 100);
        assert!(matches!(result, Err(Error::Malformed(_))));
    }
}
