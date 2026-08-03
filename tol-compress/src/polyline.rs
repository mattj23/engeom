//! Ordered polyline containers, conventionally `.tccurve2` and `.tccurve3`.
//!
//! Unlike trying to store generic 1-simplices, polylines are chains, so the ordering of the points
//! replaces the need for a block of indices. There is a flag to say if the last vertex joins back
//! to the first, but that's the only real specialization here over just storing points in general.
//!
//! # Collections
//!
//! A polyline file is an ordered collection, and cross sections, contours and toolpaths tend to
//! arrive in groups. Entries may carry a name; order is meaningful either way. See
//! [`crate::find_by_name`].
//!
//! ```
//! use tol_compress::{Polyline2, find_by_name, polyline};
//!
//! let sections = vec![
//!     Polyline2::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]], true).named("root"),
//!     Polyline2::new(vec![[0.0, 0.5], [1.0, 0.5]], false).named("tip"),
//! ];
//!
//! let mut buf = Vec::new();
//! polyline::write_to(&mut buf, &sections, 1e-5)?;
//!
//! let back: Vec<Polyline2> = polyline::read_from(&mut buf.as_slice())?;
//! assert_eq!(back.len(), 2);
//! assert!(back[0].closed);
//! assert_eq!(find_by_name(&back, "tip").map(|p| p.points.len()), Some(2));
//! # Ok::<(), tol_compress::Error>(())
//! ```

use crate::container::{self, Kind, Named, item};
use crate::error::{Error, Result};
use crate::metadata::Metadata;
use crate::points::{read_points, write_points};
use crate::raw::MAX_PREALLOC;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// An ordered polyline in `N` dimensions, with an optional name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polyline<const N: usize> {
    /// Optional identifier, preserved through a round trip. Order is meaningful regardless.
    pub name: Option<String>,
    /// Vertices, in order.
    pub points: Vec<[f64; N]>,
    /// Whether the last vertex joins back to the first.
    pub closed: bool,
    /// Caller metadata, empty when the polyline carries none. Stored, never interpreted.
    ///
    /// This is where a chord tolerance belongs, under a namespaced key.
    pub metadata: Metadata,
}

/// A polyline in the plane.
pub type Polyline2 = Polyline<2>;

/// A polyline in space.
pub type Polyline3 = Polyline<3>;

impl<const N: usize> Polyline<N> {
    /// A polyline with no name and no metadata.
    pub fn new(points: Vec<[f64; N]>, closed: bool) -> Self {
        Self {
            name: None,
            points,
            closed,
            metadata: Metadata::new(),
        }
    }

    /// The same polyline carrying a name.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The same polyline carrying one metadata entry.
    pub fn with_meta(mut self, key: impl Into<String>, value: impl Into<crate::Value>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }
}

impl<const N: usize> Named for Polyline<N> {
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

/// The container kind for a polyline of this dimension.
///
/// # Panics
///
/// Panics outside 2 and 3 dimensions. Higher-dimensional polylines have no assigned kind byte;
/// adding one later needs a new kind, not a format version, since unknown kinds already fail.
fn kind_for<const N: usize>() -> Kind {
    match N {
        2 => Kind::Polyline2,
        3 => Kind::Polyline3,
        _ => panic!("polyline containers exist only in 2 and 3 dimensions, not {N}"),
    }
}

/// Write a collection of polylines, every item at the same storage tolerance.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis is too wide to meet `tol`, and
/// [`Error::Malformed`] for a non-finite coordinate.
pub fn write_to<W: Write, const N: usize>(
    writer: &mut W,
    polylines: &[Polyline<N>],
    tol: f64,
) -> Result<()> {
    write_to_with_meta(writer, polylines, tol, &Metadata::new())
}

/// Write a collection of polylines with file-level metadata attached. See [`write_to`].
pub fn write_to_with_meta<W: Write, const N: usize>(
    writer: &mut W,
    polylines: &[Polyline<N>],
    tol: f64,
    file_metadata: &Metadata,
) -> Result<()> {
    let count = u32::try_from(polylines.len())
        .map_err(|_| Error::Malformed("container holds more items than a u32 can count"))?;
    container::write_header(writer, kind_for::<N>(), count, file_metadata)?;

    for polyline in polylines {
        write_item(writer, polyline, tol)?;
    }

    Ok(())
}

/// Write a single polyline as a one-item collection.
pub fn write_one_to<W: Write, const N: usize>(
    writer: &mut W,
    polyline: &Polyline<N>,
    tol: f64,
) -> Result<()> {
    write_to(writer, std::slice::from_ref(polyline), tol)
}

/// Read a collection of polylines.
///
/// # Errors
///
/// [`Error::Malformed`] if the file holds a different dimension than `N`.
pub fn read_from<R: Read, const N: usize>(reader: &mut R) -> Result<Vec<Polyline<N>>> {
    let header = container::read_header(reader, kind_for::<N>())?;

    let mut out = Vec::with_capacity((header.count as usize).min(MAX_PREALLOC));
    for _ in 0..header.count {
        out.push(read_item(reader)?);
    }

    Ok(out)
}

/// Read a container that holds exactly one polyline.
///
/// # Errors
///
/// [`Error::NotASingleItem`] if the container holds any other number, because returning the first
/// of several would silently discard data.
pub fn read_one_from<R: Read, const N: usize>(reader: &mut R) -> Result<Polyline<N>> {
    let header = container::read_header(reader, kind_for::<N>())?;
    if header.count != 1 {
        return Err(Error::NotASingleItem {
            found: header.count,
        });
    }
    read_item(reader)
}

/// Write a collection of polylines to a file. See [`write_to`].
pub fn write_file<const N: usize>(path: &Path, polylines: &[Polyline<N>], tol: f64) -> Result<()> {
    // Buffering is required, not tidiness: the bit reader and writer move a byte at a time, so an
    // unbuffered File would mean a syscall per byte.
    let mut writer = BufWriter::new(File::create(path)?);
    write_to(&mut writer, polylines, tol)?;
    writer.flush()?;
    Ok(())
}

/// Write a single polyline to a file as a one-item collection.
pub fn write_one_file<const N: usize>(path: &Path, polyline: &Polyline<N>, tol: f64) -> Result<()> {
    write_file(path, std::slice::from_ref(polyline), tol)
}

/// Read a collection of polylines from a file.
pub fn read_file<const N: usize>(path: &Path) -> Result<Vec<Polyline<N>>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_from(&mut reader)
}

/// Read a file holding exactly one polyline. See [`read_one_from`].
pub fn read_one_file<const N: usize>(path: &Path) -> Result<Polyline<N>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_one_from(&mut reader)
}

fn write_item<W: Write, const N: usize>(
    writer: &mut W,
    polyline: &Polyline<N>,
    tol: f64,
) -> Result<()> {
    item::write_preamble(
        writer,
        polyline.name.as_deref(),
        &polyline.metadata,
        polyline.closed,
    )?;
    write_points(writer, &polyline.points, tol)?;

    Ok(())
}

fn read_item<R: Read, const N: usize>(reader: &mut R) -> Result<Polyline<N>> {
    let preamble = item::read_preamble(reader, true)?;
    let points: Vec<[f64; N]> = read_points(reader)?;

    Ok(Polyline {
        name: preamble.name,
        points,
        closed: preamble.closed,
        metadata: preamble.metadata,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::find_by_name;
    use crate::metadata::Value;
    use crate::testgen::Rng;
    use std::io::Cursor;

    struct TempPath(std::path::PathBuf);

    impl TempPath {
        fn new(tag: &str) -> Self {
            let mut rng = Rng::new(std::process::id() as u64);
            let name = format!("tol-compress-{tag}-{:x}.tccurve", rng.next_u64());
            Self(std::env::temp_dir().join(name))
        }
    }

    impl Drop for TempPath {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn circle(n: usize, radius: f64) -> Vec<[f64; 2]> {
        (0..n)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / n as f64;
                [radius * t.cos(), radius * t.sin()]
            })
            .collect()
    }

    fn helix(n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = 6.0 * std::f64::consts::TAU * i as f64 / n as f64;
                [10.0 * t.cos(), 10.0 * t.sin(), 0.5 * t]
            })
            .collect()
    }

    fn distance<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
        (0..N).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    fn assert_matches<const N: usize>(
        original: &Polyline<N>,
        recovered: &Polyline<N>,
        tol: f64,
        what: &str,
    ) {
        assert_eq!(recovered.name, original.name, "{what}: name");
        assert_eq!(recovered.closed, original.closed, "{what}: closed");
        assert_eq!(recovered.metadata, original.metadata, "{what}: metadata");
        assert_eq!(
            recovered.points.len(),
            original.points.len(),
            "{what}: point count"
        );

        for (i, (o, r)) in original
            .points
            .iter()
            .zip(recovered.points.iter())
            .enumerate()
        {
            let d = distance(o, r);
            assert!(d <= tol, "{what}: vertex {i} recovered {d} away, tol {tol}");
        }
    }

    fn round_trip<const N: usize>(polyline: &Polyline<N>, tol: f64) -> Polyline<N> {
        let mut buf = Vec::new();
        write_one_to(&mut buf, polyline, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_one_from(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder left bytes unread"
        );
        back
    }

    #[test]
    fn open_and_closed_both_round_trip() {
        let tol = 1e-5;

        let closed = Polyline2::new(circle(200, 25.0), true);
        assert_matches(&closed, &round_trip(&closed, tol), tol, "closed");

        let open = Polyline2::new(circle(200, 25.0), false);
        assert_matches(&open, &round_trip(&open, tol), tol, "open");

        // The flag is the only difference, so the encodings must differ by nothing else.
        let mut a = Vec::new();
        let mut b = Vec::new();
        write_one_to(&mut a, &closed, tol).unwrap();
        write_one_to(&mut b, &open, tol).unwrap();
        assert_eq!(a.len(), b.len());
        assert_ne!(a, b);
    }

    #[test]
    fn three_dimensional_polylines_round_trip() {
        let tol = 1e-4;
        let p = Polyline3::new(helix(2000), false).named("helix");
        assert_matches(&p, &round_trip(&p, tol), tol, "helix");
    }

    /// Metadata is opt-in: a polyline that carries none pays nothing for the mechanism.
    #[test]
    fn metadata_is_free_when_unused() {
        let tol = 1e-5;
        let bare = Polyline2::new(circle(64, 10.0), true);
        let tagged = Polyline2::new(circle(64, 10.0), true).with_meta("engeom.chord_tol", 1e-4);

        let mut a = Vec::new();
        let mut b = Vec::new();
        write_one_to(&mut a, &bare, tol).unwrap();
        write_one_to(&mut b, &tagged, tol).unwrap();

        assert!(bare.metadata.is_empty());
        assert!(
            b.len() > a.len(),
            "metadata must cost something when present"
        );

        assert!(round_trip(&bare, tol).metadata.is_empty());
        assert_eq!(
            round_trip(&tagged, tol).metadata["engeom.chord_tol"].as_f64(),
            Some(1e-4)
        );
    }

    /// The storage tolerance and whatever a caller records about curve approximation are entirely
    /// separate. Storing points more finely than the curve approximates its source is normal, and
    /// the crate must not conflate the two or interpret the latter at all.
    #[test]
    fn storage_tolerance_and_caller_metadata_are_independent() {
        let p = Polyline2::new(circle(128, 35.0), true).with_meta("engeom.chord_tol", 1e-4);
        let back = round_trip(&p, 1e-6);

        assert_eq!(back.metadata["engeom.chord_tol"].as_f64(), Some(1e-4));
        for (o, r) in p.points.iter().zip(back.points.iter()) {
            assert!(distance(o, r) <= 1e-6);
        }
    }

    /// A value this crate would consider nonsense is still a value. Judging it belongs to whoever
    /// wrote it, which is exactly what makes the mechanism reusable.
    #[test]
    fn metadata_values_are_not_second_guessed() {
        let p = Polyline2::new(circle(8, 1.0), true)
            .with_meta("weird.nan", f64::NAN)
            .with_meta("weird.negative", -1.0)
            .with_meta("weird.zero", 0.0);

        let back = round_trip(&p, 1e-4);
        assert!(back.metadata["weird.nan"].as_f64().unwrap().is_nan());
        assert_eq!(back.metadata["weird.negative"].as_f64(), Some(-1.0));
        assert_eq!(back.metadata["weird.zero"].as_f64(), Some(0.0));
    }

    #[test]
    fn file_level_metadata_round_trips() {
        let mut file_meta = Metadata::new();
        file_meta.insert("units".into(), "mm".into());
        file_meta.insert("scanner".into(), Value::I64(7));

        let items = vec![Polyline2::new(circle(16, 2.0), true)];
        let mut buf = Vec::new();
        write_to_with_meta(&mut buf, &items, 1e-4, &file_meta).unwrap();

        let header = container::read_header(&mut Cursor::new(&buf), Kind::Polyline2).unwrap();
        assert_eq!(header.metadata, file_meta);

        // And the items still decode past it.
        let back: Vec<Polyline2> = read_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(back.len(), 1);
    }

    #[test]
    fn a_collection_preserves_order_and_mixed_naming() {
        let tol = 1e-5;
        let items: Vec<Polyline2> = (0..12)
            .map(|i| {
                let p = Polyline2::new(circle(32, 1.0 + i as f64), i % 2 == 0);
                if i % 3 == 0 {
                    p.named(format!("section-{i}"))
                } else {
                    p
                }
            })
            .collect();

        let mut buf = Vec::new();
        write_to(&mut buf, &items, tol).unwrap();
        let back: Vec<Polyline2> = read_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.len(), 12);
        for (i, (o, r)) in items.iter().zip(back.iter()).enumerate() {
            assert_matches(o, r, tol, &format!("item {i}"));
        }
        assert_eq!(back[0].name.as_deref(), Some("section-0"));
        assert_eq!(back[1].name, None);
        assert_eq!(back[3].name.as_deref(), Some("section-3"));
    }

    #[test]
    fn find_by_name_locates_an_entry() {
        let items = vec![
            Polyline2::new(circle(8, 1.0), true),
            Polyline2::new(circle(8, 2.0), true).named("wanted"),
            Polyline2::new(circle(8, 3.0), true).named("other"),
        ];

        assert!(find_by_name(&items, "wanted").is_some());
        assert_eq!(find_by_name(&items, "other").unwrap().points.len(), 8);
        assert!(find_by_name(&items, "absent").is_none());
    }

    /// Duplicate names are permitted by the format, so the documented first-match behaviour is
    /// what callers rely on.
    #[test]
    fn duplicate_names_resolve_to_the_first() {
        let items = vec![
            Polyline2::new(circle(4, 1.0), true).named("dup"),
            Polyline2::new(circle(16, 2.0), true).named("dup"),
        ];

        let found = find_by_name(&items, "dup").unwrap();
        assert_eq!(found.points.len(), 4, "should be the first of the two");
    }

    #[test]
    fn reading_at_the_wrong_dimension_is_rejected() {
        let mut buf = Vec::new();
        write_one_to(&mut buf, &Polyline2::new(circle(8, 1.0), true), 1e-4).unwrap();

        let wrong: Result<Vec<Polyline3>> = read_from(&mut Cursor::new(&buf));
        assert!(matches!(wrong, Err(Error::Malformed(_))));
    }

    #[test]
    fn degenerate_sizes_round_trip() {
        let tol = 1e-6;

        let empty = Polyline2::new(Vec::new(), false);
        assert_matches(&empty, &round_trip(&empty, tol), tol, "empty");

        let single = Polyline2::new(vec![[1.5, -2.5]], false);
        assert_matches(&single, &round_trip(&single, tol), tol, "single");

        // A closed polyline of two points is geometrically degenerate but structurally valid, and
        // the format has no business refusing it.
        let two = Polyline2::new(vec![[0.0, 0.0], [1.0, 0.0]], true);
        assert_matches(&two, &round_trip(&two, tol), tol, "two closed");
    }

    #[test]
    fn read_one_refuses_a_multi_item_file() {
        let p = Polyline2::new(circle(8, 1.0), true);
        let mut buf = Vec::new();
        write_to(&mut buf, &[p.clone(), p.clone()], 1e-4).unwrap();

        let back: Result<Polyline2> = read_one_from(&mut Cursor::new(&buf));
        assert!(matches!(back, Err(Error::NotASingleItem { found: 2 })));
    }

    #[test]
    fn round_trips_through_a_real_file() {
        let tol = 1e-5;
        let path = TempPath::new("polyline");

        let items = vec![
            Polyline2::new(circle(100, 5.0), true).named("outer"),
            Polyline2::new(circle(50, 2.0), true)
                .named("inner")
                .with_meta("engeom.chord_tol", 1e-3),
        ];

        write_file(&path.0, &items, tol).unwrap();
        let back: Vec<Polyline2> = read_file(&path.0).unwrap();

        assert_eq!(back.len(), 2);
        for (o, r) in items.iter().zip(back.iter()) {
            assert_matches(o, r, tol, "file");
        }
        assert_eq!(
            find_by_name(&back, "inner").unwrap().metadata["engeom.chord_tol"].as_f64(),
            Some(1e-3)
        );

        // And the single-item path.
        write_one_file(&path.0, &items[0], tol).unwrap();
        let one: Polyline2 = read_one_file(&path.0).unwrap();
        assert_matches(&items[0], &one, tol, "single file");
    }

    #[test]
    fn truncated_input_is_an_error() {
        let p = Polyline3::new(helix(200), false)
            .named("truncate me")
            .with_meta("engeom.chord_tol", 1e-3);
        let mut buf = Vec::new();
        write_one_to(&mut buf, &p, 1e-4).unwrap();

        for cut in [0, 5, 11, 12, 20, 30, buf.len() / 2, buf.len() - 1] {
            let back: Result<Polyline3> = read_one_from(&mut Cursor::new(&buf[..cut]));
            assert!(back.is_err(), "truncating to {cut} bytes should fail");
        }
    }

    #[test]
    #[should_panic(expected = "only in 2 and 3 dimensions")]
    fn an_unsupported_dimension_panics() {
        let mut buf = Vec::new();
        let _ = write_one_to(&mut buf, &Polyline::<4>::new(Vec::new(), false), 1e-3);
    }
}
