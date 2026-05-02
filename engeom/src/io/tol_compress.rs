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

use std::io::Write;
use crate::{Iso3, Point3, Result};
use crate::geom3::Aabb3;

pub(crate) fn write_tc_points3<W: Write>(writer: &mut W, points: &[Point3], tol: f64) -> Result<()> {
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
    let max_val = 2f64.powi((bytes as i32) * 8);
    let u_frac = fraction * max_val;
    let u = u_frac.round() as u64;

    for i in 0..bytes {
        let byte = ((u >> (i * 8)) & 0xFF) as u8;
        writer.write_all(&[byte])?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry;
}
