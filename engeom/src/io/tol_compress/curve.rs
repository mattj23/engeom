//! Curve (2d and 3d) serialization for the practical tol-compression format. This uses the
//! algorithms in the `engeom::io::tol_compress::core` module and its submodules to dynamically
//! adjust the byte width and store curve information in an efficient binary representation.
//!
//! The recommended format extensions are `.tccurve2` and `.tccurve3` for 2d and 3d curves,
//! respectively.

use crate::Result;
use crate::geom2::Curve2;
use crate::geom3::Curve3;
use crate::io::tol_compress::core::{
    read_f64, read_tc_points2, read_tc_points3, write_tc_points2, write_tc_points3,
};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

const MAGIC2: &[u8; 8] = b"TCCURVE2";
const MAGIC3: &[u8; 8] = b"TCCURVE3";

/// Serialize a 2D curve into the tccurve2 format, writing to any [`Write`] sink.
///
/// The file stores the curve's reconstruction tolerance and closed/open state, followed by the
/// vertex positions encoded as variable-width integers quantized within the point bounding box.
/// The `tol` parameter sets the maximum acceptable round-trip position error for any vertex.
/// Smaller values produce more accurate output at the cost of more bytes per vertex.
///
/// Use [`write_tc_curve2_file`] for the common case of writing directly to a file path.
pub fn write_tc_curve2_to<W: Write>(writer: &mut W, curve: &Curve2, tol: f64) -> Result<()> {
    writer.write_all(MAGIC2)?;
    writer.write_all(&curve.tol().to_le_bytes())?;
    writer.write_all(&[curve.is_closed() as u8])?;
    write_tc_points2(writer, curve.points(), tol)?;
    Ok(())
}

/// Deserialize a 2D curve from a tccurve2-format byte stream.
///
/// Returns an error if the magic bytes do not match or the data is malformed. The recovered
/// vertex positions are guaranteed to be within the tolerance that was supplied at write time,
/// and the curve's closed/open state is faithfully restored.
///
/// Use [`read_tc_curve2_file`] for the common case of reading from a file path.
pub fn read_tc_curve2_from<R: Read>(reader: &mut R) -> Result<Curve2> {
    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC2 {
        return Err("Not a tccurve2 file: invalid magic bytes".into());
    }
    let curve_tol = read_f64(reader)?;
    let mut closed_byte = [0u8; 1];
    reader.read_exact(&mut closed_byte)?;
    let is_closed = closed_byte[0] != 0;
    let points = read_tc_points2(reader)?;
    Curve2::from_points(&points, curve_tol, is_closed)
}

/// Write a 2D curve to a tccurve2 file at the given path. See [`write_tc_curve2_to`] for format details.
pub fn write_tc_curve2_file(path: &Path, curve: &Curve2, tol: f64) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_tc_curve2_to(&mut writer, curve, tol)
}

/// Read a 2D curve from a tccurve2 file at the given path. See [`read_tc_curve2_from`] for format details.
pub fn read_tc_curve2_file(path: &Path) -> Result<Curve2> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    read_tc_curve2_from(&mut reader)
}

/// Serialize a 3D curve into the tccurve3 format, writing to any [`Write`] sink.
///
/// The file stores the curve's reconstruction tolerance followed by the vertex positions encoded
/// as variable-width integers quantized within the point bounding box. The `tol` parameter sets
/// the maximum acceptable round-trip position error for any vertex. Smaller values produce more
/// accurate output at the cost of more bytes per vertex.
///
/// Use [`write_tc_curve3_file`] for the common case of writing directly to a file path.
pub fn write_tc_curve3_to<W: Write>(writer: &mut W, curve: &Curve3, tol: f64) -> Result<()> {
    writer.write_all(MAGIC3)?;
    writer.write_all(&curve.tol().to_le_bytes())?;
    write_tc_points3(writer, curve.vertices(), tol)?;
    Ok(())
}

/// Deserialize a 3D curve from a tccurve3-format byte stream.
///
/// Returns an error if the magic bytes do not match or the data is malformed. The recovered
/// vertex positions are guaranteed to be within the tolerance that was supplied at write time.
///
/// Use [`read_tc_curve3_file`] for the common case of reading from a file path.
pub fn read_tc_curve3_from<R: Read>(reader: &mut R) -> Result<Curve3> {
    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC3 {
        return Err("Not a tccurve3 file: invalid magic bytes".into());
    }
    let curve_tol = read_f64(reader)?;
    let points = read_tc_points3(reader)?;
    Curve3::from_points(&points, curve_tol)
}

/// Write a 3D curve to a tccurve3 file at the given path. See [`write_tc_curve3_to`] for format details.
pub fn write_tc_curve3_file(path: &Path, curve: &Curve3, tol: f64) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_tc_curve3_to(&mut writer, curve, tol)
}

/// Read a 3D curve from a tccurve3 file at the given path. See [`read_tc_curve3_from`] for format details.
pub fn read_tc_curve3_file(path: &Path) -> Result<Curve3> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    read_tc_curve3_from(&mut reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Point2;
    use crate::geom3::Point3;
    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;
    use std::io::Cursor;

    fn make_curve2(rg: &mut RandomGeometry3, n: usize, tol: f64) -> Curve2 {
        let points: Vec<Point2> = (0..n)
            .map(|_| Point2::new(rg.f64_sym(100.0), rg.f64_sym(100.0)))
            .collect();
        Curve2::from_points(&points, tol, false).unwrap()
    }

    fn make_closed_curve2(rg: &mut RandomGeometry3, n: usize, tol: f64) -> Curve2 {
        let mut points: Vec<Point2> = (0..n)
            .map(|_| Point2::new(rg.f64_sym(100.0), rg.f64_sym(100.0)))
            .collect();
        points.push(points[0]);
        Curve2::from_points(&points, tol, true).unwrap()
    }

    fn make_curve3(rg: &mut RandomGeometry3, n: usize, tol: f64) -> Curve3 {
        let points: Vec<Point3> = (0..n)
            .map(|_| Point3::new(rg.f64_sym(100.0), rg.f64_sym(100.0), rg.f64_sym(100.0)))
            .collect();
        Curve3::from_points(&points, tol).unwrap()
    }

    fn check_curve2_round_trip(original: &Curve2, recovered: &Curve2, tol: f64) {
        assert_eq!(original.points().len(), recovered.points().len());
        assert_eq!(original.is_closed(), recovered.is_closed());
        for (a, b) in original.points().iter().zip(recovered.points().iter()) {
            assert_relative_eq!(a, b, epsilon = tol);
        }
    }

    fn check_curve3_round_trip(original: &Curve3, recovered: &Curve3, tol: f64) {
        assert_eq!(original.vertices().len(), recovered.vertices().len());
        for (a, b) in original.vertices().iter().zip(recovered.vertices().iter()) {
            assert_relative_eq!(a, b, epsilon = tol);
        }
    }

    #[test]
    fn round_trip_curve2_bytes() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(42);
        let curve = make_curve2(&mut rg, 1000, 1e-6);

        let mut buf = Vec::new();
        write_tc_curve2_to(&mut buf, &curve, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_curve2_from(&mut cursor).unwrap();
        check_curve2_round_trip(&curve, &recovered, tol);
    }

    #[test]
    fn round_trip_curve2_closed_bytes() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(7);
        let curve = make_closed_curve2(&mut rg, 999, 1e-6);
        assert!(curve.is_closed());

        let mut buf = Vec::new();
        write_tc_curve2_to(&mut buf, &curve, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_curve2_from(&mut cursor).unwrap();
        check_curve2_round_trip(&curve, &recovered, tol);
    }

    #[test]
    fn round_trip_curve2_file() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(13);
        let curve = make_curve2(&mut rg, 500, 1e-6);

        let path = std::env::temp_dir().join("test_curve2_round_trip.tccurve2");
        write_tc_curve2_file(&path, &curve, tol).unwrap();
        let recovered = read_tc_curve2_file(&path).unwrap();
        check_curve2_round_trip(&curve, &recovered, tol);
    }

    #[test]
    fn round_trip_curve3_bytes() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(42);
        let curve = make_curve3(&mut rg, 1000, 1e-6);

        let mut buf = Vec::new();
        write_tc_curve3_to(&mut buf, &curve, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_curve3_from(&mut cursor).unwrap();
        check_curve3_round_trip(&curve, &recovered, tol);
    }

    #[test]
    fn round_trip_curve3_file() {
        let tol = 1e-4;
        let mut rg = RandomGeometry3::from_seed(99);
        let curve = make_curve3(&mut rg, 500, 1e-6);

        let path = std::env::temp_dir().join("test_curve3_round_trip.tccurve3");
        write_tc_curve3_file(&path, &curve, tol).unwrap();
        let recovered = read_tc_curve3_file(&path).unwrap();
        check_curve3_round_trip(&curve, &recovered, tol);
    }

    #[test]
    fn invalid_magic_curve2_rejected() {
        let bad = b"NOPE0000\x00\x00";
        let mut cursor = Cursor::new(bad);
        assert!(read_tc_curve2_from(&mut cursor).is_err());
    }

    #[test]
    fn invalid_magic_curve3_rejected() {
        let bad = b"NOPE0000\x00\x00";
        let mut cursor = Cursor::new(bad);
        assert!(read_tc_curve3_from(&mut cursor).is_err());
    }
}
