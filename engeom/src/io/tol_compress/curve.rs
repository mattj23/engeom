//! Curve (2d and 3d) serialization for the practical tol-compression format.
//!
//! This is a thin adapter over the [`tol_compress`] crate: it converts between engeom's [`Curve2`]
//! and [`Curve3`] and that crate's 2D and 3D polyline containers, and adds nothing of its
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
//!
//! # One file, one or many curves
//!
//! These files are ordered collections, which is what lets a whole [`crate::geom2::CurveGroup2`] or
//! [`crate::geom3::CurveGroup3`] live in one of them: a planar section is naturally several loops
//! and strands, and splitting them across files would lose both their grouping and their order.
//!
//! The single-curve functions are a convenience over that rather than a separate format, and they
//! **refuse a file holding several** rather than silently returning the first. The collection
//! functions accept any count, including one, so reading a file written by either path works.
//!
//! Every curve keeps its own chord tolerance and closed state, so a collection may freely mix
//! closed loops with open strands and curves built at different tolerances.

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
    polyline::write_one_to(writer, &curve2_item(curve), tol)?;
    Ok(())
}

/// Deserialize a 2D curve from a tccurve2-format byte stream.
///
/// The recovered vertex positions are guaranteed to be within the tolerance that was supplied at
/// write time, and the curve's chordal tolerance and closed/open state are faithfully restored.
///
/// Use [`read_tc_curve2_file`] for the common case of reading from a file path.
pub fn read_tc_curve2_from<R: Read>(reader: &mut R) -> Result<Curve2> {
    curve2_from(polyline::read_one_from::<R, 2>(reader)?)
}

/// Write a 2D curve to a tccurve2 file at the given path. See [`write_tc_curve2_to`] for format details.
pub fn write_tc_curve2_file(path: &Path, curve: &Curve2, tol: f64) -> Result<()> {
    polyline::write_one_file(path, &curve2_item(curve), tol)?;
    Ok(())
}

/// Read a 2D curve from a tccurve2 file at the given path. See [`read_tc_curve2_from`] for format details.
pub fn read_tc_curve2_file(path: &Path) -> Result<Curve2> {
    curve2_from(polyline::read_one_file::<2>(path)?)
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

// ================================================================================================
// Collections
// ================================================================================================

/// Serialize any number of 2D curves into a single tccurve2 stream, every one of them at the same
/// storage tolerance.
///
/// Order is preserved, and each curve keeps its own closed state and chord tolerance. This is what
/// a [`crate::geom2::CurveGroup2`] is written with.
///
/// Use [`write_tc_curves2_file`] for the common case of writing directly to a file path.
pub fn write_tc_curves2_to<W: Write>(writer: &mut W, curves: &[Curve2], tol: f64) -> Result<()> {
    let items: Vec<Polyline2> = curves.iter().map(curve2_item).collect();
    polyline::write_to(writer, &items, tol)?;
    Ok(())
}

/// Deserialize any number of 2D curves from a tccurve2-format byte stream, in the order they were
/// written.
///
/// A file holding a single curve reads back as a one-element vector, so this accepts anything
/// [`write_tc_curve2_to`] produces as well.
///
/// Use [`read_tc_curves2_file`] for the common case of reading from a file path.
pub fn read_tc_curves2_from<R: Read>(reader: &mut R) -> Result<Vec<Curve2>> {
    polyline::read_from::<R, 2>(reader)?
        .into_iter()
        .map(curve2_from)
        .collect()
}

/// Write any number of 2D curves to a tccurve2 file. See [`write_tc_curves2_to`].
pub fn write_tc_curves2_file(path: &Path, curves: &[Curve2], tol: f64) -> Result<()> {
    let items: Vec<Polyline2> = curves.iter().map(curve2_item).collect();
    polyline::write_file(path, &items, tol)?;
    Ok(())
}

/// Read any number of 2D curves from a tccurve2 file. See [`read_tc_curves2_from`].
pub fn read_tc_curves2_file(path: &Path) -> Result<Vec<Curve2>> {
    polyline::read_file::<2>(path)?
        .into_iter()
        .map(curve2_from)
        .collect()
}

/// Serialize any number of 3D curves into a single tccurve3 stream, every one of them at the same
/// storage tolerance.
///
/// Order is preserved, and each curve keeps its own chord tolerance. This is what a
/// [`crate::geom3::CurveGroup3`] is written with.
///
/// Use [`write_tc_curves3_file`] for the common case of writing directly to a file path.
pub fn write_tc_curves3_to<W: Write>(writer: &mut W, curves: &[Curve3], tol: f64) -> Result<()> {
    let items: Vec<Polyline3> = curves.iter().map(curve3_item).collect();
    polyline::write_to(writer, &items, tol)?;
    Ok(())
}

/// Deserialize any number of 3D curves from a tccurve3-format byte stream, in the order they were
/// written.
///
/// Use [`read_tc_curves3_file`] for the common case of reading from a file path.
pub fn read_tc_curves3_from<R: Read>(reader: &mut R) -> Result<Vec<Curve3>> {
    polyline::read_from::<R, 3>(reader)?
        .into_iter()
        .map(curve3_from)
        .collect()
}

/// Write any number of 3D curves to a tccurve3 file. See [`write_tc_curves3_to`].
pub fn write_tc_curves3_file(path: &Path, curves: &[Curve3], tol: f64) -> Result<()> {
    let items: Vec<Polyline3> = curves.iter().map(curve3_item).collect();
    polyline::write_file(path, &items, tol)?;
    Ok(())
}

/// Read any number of 3D curves from a tccurve3 file. See [`read_tc_curves3_from`].
pub fn read_tc_curves3_file(path: &Path) -> Result<Vec<Curve3>> {
    polyline::read_file::<3>(path)?
        .into_iter()
        .map(curve3_from)
        .collect()
}

/// A 2D curve as a storable polyline, carrying its closed state and chord tolerance.
fn curve2_item(curve: &Curve2) -> Polyline2 {
    let points = curve.points().iter().map(|p| [p.x, p.y]).collect();
    Polyline2::new(points, curve.is_closed()).with_meta(CHORD_TOL, curve.tol())
}

/// The inverse of [`curve2_item`].
fn curve2_from(item: Polyline2) -> Result<Curve2> {
    let tol = chord_tol(&item.metadata, "tccurve2 file")?;
    let points: Vec<_> = item
        .points
        .iter()
        .map(|p| crate::Point2::new(p[0], p[1]))
        .collect();
    Curve2::from_points(&points, tol, item.closed)
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

    // ============================================================================================
    // Collections
    // ============================================================================================

    /// A collection has to keep its order and let its members differ from each other, which is what
    /// makes it able to hold a whole curve group rather than a uniform batch.
    #[test]
    fn a_collection_of_curve2_round_trips_with_mixed_members() {
        let tol = 1e-5;
        let mut rg = RandomGeometry3::from_seed(101);
        let curves = vec![
            make_closed_curve2(&mut rg, 200, 1e-4),
            make_curve2(&mut rg, 50, 1e-6),
            make_closed_curve2(&mut rg, 75, 1e-3),
        ];

        let mut buf = Vec::new();
        write_tc_curves2_to(&mut buf, &curves, tol).unwrap();
        let back = read_tc_curves2_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.len(), curves.len());
        for (a, b) in curves.iter().zip(back.iter()) {
            check_curve2_round_trip(a, b, tol);
        }

        // The per-curve chord tolerances really are independent, not collapsed to a common one.
        assert_relative_eq!(back[0].tol(), 1e-4);
        assert_relative_eq!(back[1].tol(), 1e-6);
        assert_relative_eq!(back[2].tol(), 1e-3);
    }

    #[test]
    fn a_collection_of_curve3_round_trips() {
        let tol = 1e-5;
        let mut rg = RandomGeometry3::from_seed(202);
        let curves = vec![
            make_curve3(&mut rg, 300, 1e-6),
            make_curve3(&mut rg, 40, 1e-4),
        ];

        let mut buf = Vec::new();
        write_tc_curves3_to(&mut buf, &curves, tol).unwrap();
        let back = read_tc_curves3_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.len(), curves.len());
        for (a, b) in curves.iter().zip(back.iter()) {
            check_curve3_round_trip(a, b, tol);
        }
    }

    #[test]
    fn a_collection_round_trips_through_a_real_file() {
        let tol = 1e-5;
        let mut rg = RandomGeometry3::from_seed(303);
        let curves = vec![
            make_closed_curve2(&mut rg, 120, 1e-4),
            make_curve2(&mut rg, 60, 1e-4),
        ];

        let path = std::env::temp_dir().join("engeom_tc_curves2_collection.tccurve2");
        write_tc_curves2_file(&path, &curves, tol).unwrap();
        let back = read_tc_curves2_file(&path).unwrap();

        assert_eq!(back.len(), 2);
        for (a, b) in curves.iter().zip(back.iter()) {
            check_curve2_round_trip(a, b, tol);
        }
        let _ = std::fs::remove_file(&path);
    }

    /// The two paths write the same format, so a file from either has to be readable by the
    /// collection reader. Only the single reader is picky, and only about count.
    #[test]
    fn the_collection_reader_accepts_a_single_curve_file() {
        let mut rg = RandomGeometry3::from_seed(404);
        let curve = make_curve2(&mut rg, 30, 1e-6);

        let mut buf = Vec::new();
        write_tc_curve2_to(&mut buf, &curve, 1e-5).unwrap();

        let back = read_tc_curves2_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(back.len(), 1);
        check_curve2_round_trip(&curve, &back[0], 1e-5);
    }

    /// ...and the single reader refuses a collection rather than returning its first member, which
    /// would discard the rest with no sign that it happened.
    #[test]
    fn the_single_reader_refuses_a_collection() {
        let mut rg = RandomGeometry3::from_seed(505);
        let curves = vec![
            make_curve2(&mut rg, 20, 1e-6),
            make_curve2(&mut rg, 20, 1e-6),
        ];

        let mut buf = Vec::new();
        write_tc_curves2_to(&mut buf, &curves, 1e-5).unwrap();

        let refused = read_tc_curve2_from(&mut Cursor::new(&buf)).is_err();
        assert!(refused, "a two-curve file must not read back as one curve");
    }

    /// An empty collection is a valid file. It is the group constructors, not the format, which
    /// decide that a body needs at least one curve.
    #[test]
    fn an_empty_collection_is_a_valid_file() {
        let mut buf = Vec::new();
        write_tc_curves2_to(&mut buf, &[], 1e-5).unwrap();

        let back = read_tc_curves2_from(&mut Cursor::new(&buf)).unwrap();
        assert!(back.is_empty());
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
