//! Curve (2d and 3d) serialization for the practical tol-compression format.
//!
//! This is a thin adapter over the [`tol_compress`] crate: it converts between engeom's [`Curve2`]
//! and [`Curve3`] and that crate's dimension-generic polyline containers, and adds nothing of its
//! own. The `tol` parameter throughout is the maximum acceptable round-trip position error for any
//! vertex, in the same units as the coordinates.
//!
//! The recommended format extensions are `.tccurve2` and `.tccurve3` for 2d and 3d curves,
//! respectively.
//!
//! A curve carries a chordal reconstruction tolerance that is not a property of the stored points,
//! so it travels as item metadata under [`CHORD_TOL`] rather than as a field of the format. The
//! `tol-compress` crate stores it and never interprets it, which is the point: a general-purpose
//! geometry format should not have to learn what engeom means by a chord tolerance.

use crate::Result;
use crate::geom2::Curve2;
use crate::geom3::Curve3;
use std::io::{Read, Write};
use std::path::Path;
use tol_compress::{Polyline2, Polyline3, polyline};

/// The metadata key a curve's chordal reconstruction tolerance is stored under.
///
/// Namespaced because the metadata map is shared with any other writer of these files, and an
/// unprefixed `chord_tol` would be an obvious thing for a second one to also claim.
pub const CHORD_TOL: &str = "engeom.chord_tol";

/// Pull the chordal tolerance back out of an item's metadata.
///
/// Missing or wrongly typed is an error rather than a default. Substituting a plausible value would
/// silently change how the recovered curve resamples, which is a worse failure than refusing a file
/// engeom did not write.
fn chord_tol(metadata: &tol_compress::Metadata, what: &str) -> Result<f64> {
    metadata
        .get(CHORD_TOL)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| {
            format!("{what} has no `{CHORD_TOL}` metadata; was it written by engeom?").into()
        })
}

/// Serialize a 2D curve into the tccurve2 format, writing to any [`Write`] sink.
///
/// The vertex positions are quantized to the narrowest bit width that keeps every one of them
/// within `tol` of where it started; the curve's chordal tolerance and closed/open state travel
/// alongside them. Smaller `tol` values produce more accurate output at the cost of more bytes per
/// vertex.
///
/// Use [`write_tc_curve2_file`] for the common case of writing directly to a file path.
pub fn write_tc_curve2_to<W: Write>(writer: &mut W, curve: &Curve2, tol: f64) -> Result<()> {
    let points = curve.points().iter().map(|p| [p.x, p.y]).collect();
    let item = Polyline2::new(points, curve.is_closed()).with_meta(CHORD_TOL, curve.tol());
    polyline::write_one_to(writer, &item, tol)?;
    Ok(())
}

/// Deserialize a 2D curve from a tccurve2-format byte stream.
///
/// The recovered vertex positions are guaranteed to be within the tolerance that was supplied at
/// write time, and the curve's chordal tolerance and closed/open state are faithfully restored.
///
/// Use [`read_tc_curve2_file`] for the common case of reading from a file path.
pub fn read_tc_curve2_from<R: Read>(reader: &mut R) -> Result<Curve2> {
    let item = polyline::read_one_from::<R, 2>(reader)?;
    let tol = chord_tol(&item.metadata, "tccurve2 file")?;
    let points: Vec<_> = item
        .points
        .iter()
        .map(|p| crate::Point2::new(p[0], p[1]))
        .collect();
    Curve2::from_points(&points, tol, item.closed)
}

/// Write a 2D curve to a tccurve2 file at the given path. See [`write_tc_curve2_to`] for format details.
pub fn write_tc_curve2_file(path: &Path, curve: &Curve2, tol: f64) -> Result<()> {
    let points = curve.points().iter().map(|p| [p.x, p.y]).collect();
    let item = Polyline2::new(points, curve.is_closed()).with_meta(CHORD_TOL, curve.tol());
    polyline::write_one_file(path, &item, tol)?;
    Ok(())
}

/// Read a 2D curve from a tccurve2 file at the given path. See [`read_tc_curve2_from`] for format details.
pub fn read_tc_curve2_file(path: &Path) -> Result<Curve2> {
    let item = polyline::read_one_file::<2>(path)?;
    let tol = chord_tol(&item.metadata, "tccurve2 file")?;
    let points: Vec<_> = item
        .points
        .iter()
        .map(|p| crate::Point2::new(p[0], p[1]))
        .collect();
    Curve2::from_points(&points, tol, item.closed)
}

/// Serialize a 3D curve into the tccurve3 format, writing to any [`Write`] sink.
///
/// The vertex positions are quantized to the narrowest bit width that keeps every one of them
/// within `tol` of where it started, and the curve's chordal tolerance travels alongside them.
/// Smaller `tol` values produce more accurate output at the cost of more bytes per vertex.
///
/// Use [`write_tc_curve3_file`] for the common case of writing directly to a file path.
pub fn write_tc_curve3_to<W: Write>(writer: &mut W, curve: &Curve3, tol: f64) -> Result<()> {
    polyline::write_one_to(writer, &curve3_item(curve), tol)?;
    Ok(())
}

/// Deserialize a 3D curve from a tccurve3-format byte stream.
///
/// The recovered vertex positions are guaranteed to be within the tolerance that was supplied at
/// write time.
///
/// Use [`read_tc_curve3_file`] for the common case of reading from a file path.
pub fn read_tc_curve3_from<R: Read>(reader: &mut R) -> Result<Curve3> {
    curve3_from(polyline::read_one_from::<R, 3>(reader)?)
}

/// Write a 3D curve to a tccurve3 file at the given path. See [`write_tc_curve3_to`] for format details.
pub fn write_tc_curve3_file(path: &Path, curve: &Curve3, tol: f64) -> Result<()> {
    polyline::write_one_file(path, &curve3_item(curve), tol)?;
    Ok(())
}

/// Read a 3D curve from a tccurve3 file at the given path. See [`read_tc_curve3_from`] for format details.
pub fn read_tc_curve3_file(path: &Path) -> Result<Curve3> {
    curve3_from(polyline::read_one_file::<3>(path)?)
}

/// A 3D curve as a storable polyline. `Curve3` has no notion of closure, so the flag is always
/// false and a file claiming otherwise is refused on read rather than quietly flattened.
fn curve3_item(curve: &Curve3) -> Polyline3 {
    let points = curve.vertices().iter().map(|p| [p.x, p.y, p.z]).collect();
    Polyline3::new(points, false).with_meta(CHORD_TOL, curve.tol())
}

/// The inverse of [`curve3_item`].
fn curve3_from(item: Polyline3) -> Result<Curve3> {
    if item.closed {
        return Err("tccurve3 file is marked closed, which a Curve3 cannot represent".into());
    }
    let tol = chord_tol(&item.metadata, "tccurve3 file")?;
    let points: Vec<_> = item
        .points
        .iter()
        .map(|p| crate::Point3::new(p[0], p[1], p[2]))
        .collect();
    Curve3::from_points(&points, tol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom2::Point2;
    use crate::geom3::Point3;
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
        assert_eq!(original.tol(), recovered.tol());
        for (a, b) in original.points().iter().zip(recovered.points().iter()) {
            assert_relative_eq!(a, b, epsilon = tol);
        }
    }

    fn check_curve3_round_trip(original: &Curve3, recovered: &Curve3, tol: f64) {
        assert_eq!(original.vertices().len(), recovered.vertices().len());
        assert_eq!(original.tol(), recovered.tol());
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

    /// A tccurve2 written by something other than engeom has no chord tolerance in it. Guessing one
    /// would change how the curve resamples without any sign that it happened.
    #[test]
    fn missing_chord_tol_is_rejected() {
        let item = Polyline2::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], false);
        let mut buf = Vec::new();
        polyline::write_one_to(&mut buf, &item, 1e-4).unwrap();

        // Matched rather than unwrapped: `Curve2` is not `Debug`, so `unwrap_err` will not compile.
        let mut cursor = Cursor::new(&buf);
        let err = match read_tc_curve2_from(&mut cursor) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("a file with no chord tolerance must not read back as a curve"),
        };
        assert!(err.contains(CHORD_TOL), "unhelpful error: {err}");
    }

    /// Reading a 2d curve out of a 3d file, or the reverse, must fail rather than reinterpret the
    /// coordinates. The container records its kind, so this is the crate's job, but engeom is the
    /// one that would suffer if it ever stopped being true.
    #[test]
    fn wrong_dimension_is_rejected() {
        let mut rg = RandomGeometry3::from_seed(5);
        let curve = make_curve3(&mut rg, 50, 1e-6);

        let mut buf = Vec::new();
        write_tc_curve3_to(&mut buf, &curve, 1e-4).unwrap();

        let mut cursor = Cursor::new(&buf);
        assert!(read_tc_curve2_from(&mut cursor).is_err());
    }

    /// The committed airfoil fixture is read by the airfoil test suite and will outlive any memory
    /// of how it was converted, so its contents are asserted here rather than assumed.
    #[test]
    fn airfoil_fixture_reads_back() {
        let curve = crate::tests::airfoil_curve();
        assert_eq!(curve.points().len(), 544);
        assert!(curve.is_closed());
        assert_relative_eq!(curve.tol(), 1e-4);
    }
}
