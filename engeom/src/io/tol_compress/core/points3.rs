//! Common tolerance compression tools for 3D geometry

use crate::geom3::Aabb3;
use crate::io::tol_compress::core::{
    bytes_for_tol, read_f64, read_u8, read_u32, read_var_width, write_var_width,
};
use crate::na::{Quaternion, Translation3, UnitQuaternion};
use crate::{Iso3, Point3};
use std::io::{Read, Write};

pub fn write_tc_points3<W: Write>(
    writer: &mut W,
    points: &[Point3],
    tol: f64,
) -> crate::Result<()> {
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

pub fn read_tc_points3<R: Read>(reader: &mut R) -> crate::Result<Vec<Point3>> {
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

fn write_iso3<W: Write>(writer: &mut W, iso: &Iso3) -> crate::Result<()> {
    writer.write_all(&iso.rotation.i.to_le_bytes())?;
    writer.write_all(&iso.rotation.j.to_le_bytes())?;
    writer.write_all(&iso.rotation.k.to_le_bytes())?;
    writer.write_all(&iso.rotation.w.to_le_bytes())?;
    writer.write_all(&iso.translation.x.to_le_bytes())?;
    writer.write_all(&iso.translation.y.to_le_bytes())?;
    writer.write_all(&iso.translation.z.to_le_bytes())?;

    Ok(())
}

fn write_aabb3<W: Write>(writer: &mut W, aabb3: &Aabb3) -> crate::Result<()> {
    writer.write_all(&aabb3.mins.x.to_le_bytes())?;
    writer.write_all(&aabb3.mins.y.to_le_bytes())?;
    writer.write_all(&aabb3.mins.z.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.x.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.y.to_le_bytes())?;
    writer.write_all(&aabb3.maxs.z.to_le_bytes())?;
    Ok(())
}

fn read_iso3<R: Read>(reader: &mut R) -> crate::Result<Iso3> {
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

fn read_aabb3<R: Read>(reader: &mut R) -> crate::Result<Aabb3> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry3;
    use std::io::Cursor;

    // bytes_for_count thresholds: 2 bytes for max_index < 65536,
    //                             3 bytes for max_index < 16_777_216,
    //                             4 bytes for max_index <= u32::MAX

    #[test]
    fn round_trip_points3() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(42);
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
