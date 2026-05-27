//! This module has the tools for reading and writing indices used in the tolerance compression
//! scheme.

use std::io::{Read, Write};
use crate::io::tol_compress::core::read_u32;

/// Detect the total number of _whole_ bytes needed to store indices up to a maximum count value of
/// `total_items`. For example, a single byte can store up to 256 items, two bytes can store 66536,
/// etc.
///
/// # Arguments
///
/// * `total_items`: The total number of items to be stored, assuming that the indices run from 0
///   up to the last item.
///
/// returns: u8
pub fn bytes_for_count(total_items: u32) -> u8 {
    let mut bytes: u8 = 2;
    while (1u64 << (bytes as u64 * 8)) <= total_items as u64 {
        bytes += 1;
    }
    bytes
}

pub fn write_indices<W: Write, const N: usize>(
    writer: &mut W,
    indices: &[[u32; N]],
    max_index: u32,
) -> crate::Result<()> {
    // Find the bytes needed to represent the highest valued index.
    let bytes = bytes_for_count(max_index);

    // Write the number of indices
    writer.write_all(&(indices.len() as u32).to_le_bytes())?;

    // Now write the indices
    for i in indices {
        writer.write_all(&i[0].to_le_bytes()[..bytes as usize])?;
        writer.write_all(&i[1].to_le_bytes()[..bytes as usize])?;
        writer.write_all(&i[2].to_le_bytes()[..bytes as usize])?;
    }

    Ok(())
}

pub fn read_indices<R: Read, const N: usize>(
    reader: &mut R,
    max_index: u32,
) -> crate::Result<Vec<[u32; N]>> {
    let bytes = bytes_for_count(max_index) as usize;
    let n = read_u32(reader)? as usize;
    let mut indices = Vec::with_capacity(n);
    for _ in 0..n {
        let mut entry = [0u32; N];
        for e in entry.iter_mut() {
            let mut buf = [0u8; 4];
            reader.read_exact(&mut buf[..bytes])?;
            *e = u32::from_le_bytes(buf);
        }
        indices.push(entry);
    }
    Ok(indices)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // bytes_for_count thresholds: 2 bytes for max_index < 65536,
    //                             3 bytes for max_index < 16_777_216,
    //                             4 bytes for max_index <= u32::MAX

    #[test]
    fn round_trip_indices_2_bytes() {
        // max_index=1000 → bytes_for_count returns 2 (65536 > 1000)
        let indices: Vec<[u32; 3]> = vec![[0, 500, 1000], [1, 999, 500], [250, 750, 1000]];
        let mut buf = Vec::new();
        write_indices(&mut buf, &indices, 1000).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_indices::<_, 3>(&mut cursor, 1000).unwrap();
        assert_eq!(indices, recovered);
    }

    #[test]
    fn round_trip_indices_3_bytes() {
        // max_index=100_000 → bytes_for_count returns 3 (16_777_216 > 100_000)
        let indices: Vec<[u32; 3]> = vec![
            [0, 50_000, 100_000],
            [1, 99_999, 50_000],
            [25_000, 75_000, 100_000],
        ];
        let mut buf = Vec::new();
        write_indices(&mut buf, &indices, 100_000).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_indices::<_, 3>(&mut cursor, 100_000).unwrap();
        assert_eq!(indices, recovered);
    }

    #[test]
    fn round_trip_indices_4_bytes() {
        // max_index=20_000_000 → bytes_for_count returns 4 (4_294_967_296 > 20_000_000)
        let indices: Vec<[u32; 3]> = vec![
            [0, 10_000_000, 20_000_000],
            [1, 19_999_999, 10_000_000],
            [5_000_000, 15_000_000, 20_000_000],
        ];
        let mut buf = Vec::new();
        write_indices(&mut buf, &indices, 20_000_000).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_indices::<_, 3>(&mut cursor, 20_000_000).unwrap();
        assert_eq!(indices, recovered);
    }
}
