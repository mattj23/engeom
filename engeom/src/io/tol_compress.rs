//! This module contains common tools for serialization that uses the practical tolerance-compressed
//! method of storing spatial coordinates.
//!
//! This method attempts to take advantage of two specific features of practical measurement data:
//!
//! 1. Real world measurement systems have known precisions below which differences in position
//!    do not represent meaningful information, and users of measurement data have knowledge about
//!    where increasingly small differences in values cease to be relevant to their use case.
//! 2. Values produced by the measurement of physical objects rarely span more than a few orders
//!    of magnitude more than the smallest meaningful precision for the measurement system and/or
//!    the application.
//!
//! If a user can supply a tolerance that represents the largest acceptable round-trip accuracy of
//! any value in their dataset, the storage algorithm can find a representation scheme that uses
//! the smallest amount of bytes necessary to guarantee every position is recovered within that
//! tolerance.

use crate::geom3::Aabb3;
use crate::na::{Quaternion, Translation3, UnitQuaternion};
use crate::{Iso3, Point3, Result};
use std::io::{Read, Write};

pub(crate) fn bytes_for_count(total_items: u32) -> u8 {
    let mut bytes: u8 = 2;
    while (1u64 << (bytes as u64 * 8)) <= total_items as u64 {
        bytes += 1;
    }
    bytes
}

pub(crate) fn write_indices<W: Write, const N: usize>(
    writer: &mut W,
    indices: &[[u32; N]],
    max_index: u32,
) -> Result<()> {
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

pub(crate) fn read_indices<R: Read, const N: usize>(
    reader: &mut R,
    max_index: u32,
) -> Result<Vec<[u32; N]>> {
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

pub(crate) fn write_tc_points3<W: Write>(
    writer: &mut W,
    points: &[Point3],
    tol: f64,
) -> Result<()> {
    // In the future, this can be replaced with a more complicated partitioning algorithm if it
    // proves to actually produce real improvements.  For now, though, we'll only use one box
    // and a single identity transform.

    // Write the total number of partitions as an u32
    writer.write_all(&1u32.to_le_bytes())?;

    // Now we would iterate through the partitions, but we only have one in this implementation
    // -----------------------------------------------------------------------------------------
    // Write the number of points as an u32
    writer.write_all(&(points.len() as u32).to_le_bytes())?;

    // Make and write the recovery information
    let iso = Iso3::identity();
    write_iso3(writer, &iso)?;
    let bounds = Aabb3::from_points_ref(points);
    write_aabb3(writer, &bounds)?;

    // Now we need to figure out the number of bytes for each of the dimensions
    let d_tol = tol / 3f64.sqrt();
    let x_bytes = bytes_for_tol(bounds.extents().x, d_tol);
    let y_bytes = bytes_for_tol(bounds.extents().y, d_tol);
    let z_bytes = bytes_for_tol(bounds.extents().z, d_tol);

    // Now we write the byte counts
    writer.write_all(&x_bytes.to_le_bytes())?;
    writer.write_all(&y_bytes.to_le_bytes())?;
    writer.write_all(&z_bytes.to_le_bytes())?;

    // Finally, we will convert and write the actual vertices
    for p in points {
        let t = &iso * p;
        write_var_width(writer, (t.x - bounds.mins.x) / bounds.extents().x, x_bytes)?;
        write_var_width(writer, (t.y - bounds.mins.y) / bounds.extents().y, y_bytes)?;
        write_var_width(writer, (t.z - bounds.mins.z) / bounds.extents().z, z_bytes)?;
    }

    Ok(())
}

pub(crate) fn read_tc_points3<R: Read>(reader: &mut R) -> Result<Vec<Point3>> {
    let n_partitions = read_u32(reader)? as usize;
    let mut all_points = Vec::new();

    for _ in 0..n_partitions {
        let n_points = read_u32(reader)? as usize;
        let iso = read_iso3(reader)?;
        let bounds = read_aabb3(reader)?;
        let x_bytes = read_u8(reader)?;
        let y_bytes = read_u8(reader)?;
        let z_bytes = read_u8(reader)?;

        let inv_iso = iso.inverse();
        let ext = bounds.extents();

        for _ in 0..n_points {
            let fx = read_var_width(reader, x_bytes)?;
            let fy = read_var_width(reader, y_bytes)?;
            let fz = read_var_width(reader, z_bytes)?;
            let t = Point3::new(
                fx * ext.x + bounds.mins.x,
                fy * ext.y + bounds.mins.y,
                fz * ext.z + bounds.mins.z,
            );
            all_points.push(&inv_iso * t);
        }
    }

    Ok(all_points)
}

fn write_iso3<W: Write>(writer: &mut W, iso: &Iso3) -> Result<()> {
    writer.write_all(&iso.rotation.i.to_le_bytes())?;
    writer.write_all(&iso.rotation.j.to_le_bytes())?;
    writer.write_all(&iso.rotation.k.to_le_bytes())?;
    writer.write_all(&iso.rotation.w.to_le_bytes())?;
    writer.write_all(&iso.translation.x.to_le_bytes())?;
    writer.write_all(&iso.translation.y.to_le_bytes())?;
    writer.write_all(&iso.translation.z.to_le_bytes())?;

    Ok(())
}

fn write_aabb3<W: Write>(writer: &mut W, aabb3: &Aabb3) -> Result<()> {
    writer.write_all(&aabb3.mins.x.to_le_bytes())?;
    writer.write_all(&aabb3.mins.y.to_le_bytes())?;
    writer.write_all(&aabb3.mins.z.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.x.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.y.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.z.to_le_bytes())?;
    Ok(())
}

fn bytes_for_tol(range: f64, tol: f64) -> u8 {
    let mut bytes = 2;
    while range / 2f64.powi((bytes as i32) * 8) > tol {
        bytes += 1;
    }

    bytes
}

fn write_var_width<W: Write>(writer: &mut W, fraction: f64, bytes: u8) -> Result<()> {
    // Using max_int = 2^(bytes*8) - 1 ensures u fits in `bytes` bytes (avoids overflow at
    // fraction = 1.0 which would produce 2^(bytes*8) and silently truncate to 0).
    let max_int = 2u64.pow((bytes as u32) * 8) - 1;
    let u = (fraction * max_int as f64).round() as u64;

    for i in 0..bytes {
        let byte = ((u >> (i * 8)) & 0xFF) as u8;
        writer.write_all(&[byte])?;
    }

    Ok(())
}

fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_f64<R: Read>(reader: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

fn read_iso3<R: Read>(reader: &mut R) -> Result<Iso3> {
    let i = read_f64(reader)?;
    let j = read_f64(reader)?;
    let k = read_f64(reader)?;
    let w = read_f64(reader)?;
    let tx = read_f64(reader)?;
    let ty = read_f64(reader)?;
    let tz = read_f64(reader)?;
    let rotation = UnitQuaternion::new_unchecked(Quaternion::new(w, i, j, k));
    let translation = Translation3::new(tx, ty, tz);
    Ok(Iso3::from_parts(translation, rotation))
}

fn read_aabb3<R: Read>(reader: &mut R) -> Result<Aabb3> {
    let min_x = read_f64(reader)?;
    let min_y = read_f64(reader)?;
    let min_z = read_f64(reader)?;
    let max_x = read_f64(reader)?;
    let max_y = read_f64(reader)?;
    let max_z = read_f64(reader)?;
    Ok(Aabb3::new(
        Point3::new(min_x, min_y, min_z),
        Point3::new(max_x, max_y, max_z),
    ))
}

fn read_var_width<R: Read>(reader: &mut R, bytes: u8) -> Result<f64> {
    let mut u: u64 = 0;
    for i in 0..bytes {
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf)?;
        u |= (buf[0] as u64) << (i * 8);
    }
    let max_int = 2u64.pow((bytes as u32) * 8) - 1;
    Ok(u as f64 / max_int as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry;
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

    #[test]
    fn round_trip_points3() {
        let tol = 1e-4;
        let mut rg = RandomGeometry::from_seed(42);
        let points: Vec<Point3> = (0..100000).map(|_| rg.point3(100.0)).collect();

        let mut buf = Vec::new();
        write_tc_points3(&mut buf, &points, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_points3(&mut cursor).unwrap();

        assert_eq!(points.len(), recovered.len());
        for (original, recovered) in points.iter().zip(recovered.iter()) {
            let dist = (original - recovered).norm();
            assert!(dist <= tol, "point error {dist} exceeds tolerance {tol}");
        }
    }
}
