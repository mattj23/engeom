//! Common tolerance compression tools for 2D geometry

use crate::geom2::Aabb2;
use crate::io::tol_compress::core::{
    bytes_for_tol, read_f64, read_u8, read_u32, read_var_width, write_var_width,
};
use crate::na::UnitComplex;
use crate::{Iso2, Point2};
use parry2d_f64::na::Translation2;
use std::io::{Read, Write};

pub fn write_tc_points2<W: Write>(
    writer: &mut W,
    points: &[Point2],
    tol: f64,
) -> crate::Result<()> {
    // Write the total number of partitions as an u32
    writer.write_all(&1u32.to_le_bytes())?;

    // Write the number of points as an u32
    writer.write_all(&(points.len() as u32).to_le_bytes())?;

    // Make and write the recovery information
    let iso = Iso2::identity();
    write_iso2(writer, &iso)?;
    let bounds = Aabb2::from_points(points.iter().cloned());
    write_aabb2(writer, &bounds)?;

    // Figure out the number of bytes for each dimension
    let d_tol = tol / 2f64.sqrt();
    let x_bytes = bytes_for_tol(bounds.extents().x, d_tol);
    let y_bytes = bytes_for_tol(bounds.extents().y, d_tol);

    writer.write_all(&x_bytes.to_le_bytes())?;
    writer.write_all(&y_bytes.to_le_bytes())?;

    for p in points {
        let t = &iso * p;
        write_var_width(writer, (t.x - bounds.mins.x) / bounds.extents().x, x_bytes)?;
        write_var_width(writer, (t.y - bounds.mins.y) / bounds.extents().y, y_bytes)?;
    }

    Ok(())
}

pub fn read_tc_points2<R: Read>(reader: &mut R) -> crate::Result<Vec<Point2>> {
    let n_partitions = read_u32(reader)? as usize;
    let mut all_points = Vec::new();

    for _ in 0..n_partitions {
        let n_points = read_u32(reader)? as usize;
        let iso = read_iso2(reader)?;
        let bounds = read_aabb2(reader)?;
        let x_bytes = read_u8(reader)?;
        let y_bytes = read_u8(reader)?;

        let inv_iso = iso.inverse();
        let ext = bounds.extents();

        for _ in 0..n_points {
            let fx = read_var_width(reader, x_bytes)?;
            let fy = read_var_width(reader, y_bytes)?;
            let t = Point2::new(fx * ext.x + bounds.mins.x, fy * ext.y + bounds.mins.y);
            all_points.push(&inv_iso * t);
        }
    }

    Ok(all_points)
}

fn write_iso2<W: Write>(writer: &mut W, iso: &Iso2) -> crate::Result<()> {
    writer.write_all(&iso.rotation.angle().to_le_bytes())?;
    writer.write_all(&iso.translation.x.to_le_bytes())?;
    writer.write_all(&iso.translation.y.to_le_bytes())?;
    Ok(())
}

fn write_aabb2<W: Write>(writer: &mut W, aabb2: &Aabb2) -> crate::Result<()> {
    writer.write_all(&aabb2.mins.x.to_le_bytes())?;
    writer.write_all(&aabb2.mins.y.to_le_bytes())?;
    writer.write_all(&aabb2.maxs.x.to_le_bytes())?;
    writer.write_all(&aabb2.maxs.y.to_le_bytes())?;
    Ok(())
}

fn read_iso2<R: Read>(reader: &mut R) -> crate::Result<Iso2> {
    let angle = read_f64(reader)?;
    let tx = read_f64(reader)?;
    let ty = read_f64(reader)?;
    let rotation = UnitComplex::new(angle);
    let translation = Translation2::new(tx, ty);
    Ok(Iso2::from_parts(translation, rotation))
}

fn read_aabb2<R: Read>(reader: &mut R) -> crate::Result<Aabb2> {
    let min_x = read_f64(reader)?;
    let min_y = read_f64(reader)?;
    let max_x = read_f64(reader)?;
    let max_y = read_f64(reader)?;
    Ok(Aabb2::new(
        Point2::new(min_x, min_y),
        Point2::new(max_x, max_y),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry3;
    use std::io::Cursor;

    #[test]
    fn round_trip_points2() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(42);
        let points: Vec<Point2> = (0..100000)
            .map(|_| Point2::new(rg.f64_sym(100.0), rg.f64_sym(100.0)))
            .collect();

        let mut buf = Vec::new();
        write_tc_points2(&mut buf, &points, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_points2(&mut cursor).unwrap();

        assert_eq!(points.len(), recovered.len());
        for (original, recovered) in points.iter().zip(recovered.iter()) {
            let dist = (original - recovered).norm();
            assert!(dist <= tol, "point error {dist} exceeds tolerance {tol}");
        }
    }
}
