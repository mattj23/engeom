//! Row-organized point containers, conventionally `.tcrpf3`.
//!
//! A rasterizing sensor produces a cloud whose points are grouped into ordered rows, such as a
//! snapshot sensor's raster rows, a laser line's profiles, or a time-of-flight frame's scanlines.
//! Preserving this grouping makes the points meshable without surface reconstruction. Converting
//! the points to an unordered cloud discards the grouping.
//!
//! This format is more general than a raster. Rows can contain different numbers of points,
//! ordinals can have gaps where rows were dropped, and column indices are optional. For example, a
//! laser-line sweep identifies each point's profile but has no meaningful index across that
//! profile. The format guarantees that points remain grouped and ordered and that each row's
//! position in the sweep is recorded. See [`crate::rows`] for the grouping block.
//!
//! This module does not interpret the geometry. Record row pitch, sensor frame, handedness, and
//! units in [`Metadata`]; the crate assigns no meaning to a row.
//!
//! ```
//! use tol_compress::{PointRow3, RowPoints3, row_points};
//!
//! let scan = RowPoints3::new(vec![
//!     PointRow3::new(0, vec![[0.0, 0.0, 1.0], [0.1, 0.0, 1.0]]),
//!     PointRow3::new(1, vec![[0.0, 0.1, 1.0], [0.1, 0.1, 1.0]]),
//! ])
//! .named("exposure 500")
//! .with_meta("engeom.row_pitch", 0.1);
//!
//! let mut buf = Vec::new();
//! row_points::write_one_to(&mut buf, &scan, 1e-4)?;
//!
//! let back = row_points::read_one_from(&mut buf.as_slice())?;
//! assert_eq!(back.rows.len(), 2);
//! assert_eq!(back.rows[1].ordinal, 1);
//! # Ok::<(), tol_compress::Error>(())
//! ```

use crate::container::{self, Kind, Named, item};
use crate::effort::Effort;
use crate::error::{Error, Result};
use crate::metadata::Metadata;
use crate::points::{read_points, write_points_with};
use crate::raw::MAX_PREALLOC;
use crate::rows::{RowLayout, read_rows, write_rows};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// One row of points in x, y, z space.
///
/// The `ordinal` records the row's position in the sensor sweep, independent of its position in the
/// file. Dropping a row therefore leaves a gap. Rows must appear in increasing ordinal order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointRow3 {
    /// Where this row sits in the sensor's sweep.
    pub ordinal: u32,
    /// The row's points, in the order the sensor produced them.
    pub points: Vec<[f64; 3]>,
    /// Each point's strictly increasing column index. This is `None` when the source has no
    /// meaningful index across a row. Either every row has column indices or none does.
    pub columns: Option<Vec<u32>>,
}

impl PointRow3 {
    /// A row with no column indices.
    pub fn new(ordinal: u32, points: Vec<[f64; 3]>) -> Self {
        Self {
            ordinal,
            points,
            columns: None,
        }
    }

    /// The same row carrying a column index for every point.
    pub fn with_columns(mut self, columns: Vec<u32>) -> Self {
        self.columns = Some(columns);
        self
    }
}

/// Points in x, y, z space grouped into ordered rows, with an optional name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RowPoints3 {
    /// Optional identifier, preserved through a round trip.
    pub name: Option<String>,
    /// The rows, in increasing ordinal order.
    pub rows: Vec<PointRow3>,
    /// Caller metadata, empty when the item carries none. Stored, never interpreted.
    pub metadata: Metadata,
}

impl RowPoints3 {
    /// A row-organized point set with no name and no metadata.
    pub fn new(rows: Vec<PointRow3>) -> Self {
        Self {
            name: None,
            rows,
            metadata: Metadata::new(),
        }
    }

    /// The same item carrying a name.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The same item carrying one metadata entry.
    pub fn with_meta(mut self, key: impl Into<String>, value: impl Into<crate::Value>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }

    /// How many points the rows hold in total.
    pub fn point_count(&self) -> usize {
        self.rows.iter().map(|row| row.points.len()).sum()
    }

    /// Validate all format requirements without writing data.
    ///
    /// [`write_to`] and related functions call this method automatically. Call it directly only to
    /// validate before writing.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Malformed`] with a distinct message for an empty row, non-increasing
    /// ordinals, inconsistent column presence between rows, a column count that differs from the
    /// row's point count, non-increasing columns within a row, or counts too large for `u32`.
    pub fn validate(&self) -> Result<()> {
        self.layout()?.validate()
    }

    /// The grouping, as the row block stores it.
    fn layout(&self) -> Result<RowLayout> {
        let has_columns = self.rows.first().map(|row| row.columns.is_some());

        let mut ordinals = Vec::with_capacity(self.rows.len());
        let mut counts = Vec::with_capacity(self.rows.len());
        let mut columns = Vec::new();

        for row in &self.rows {
            if row.columns.is_some() != has_columns.unwrap_or(false) {
                return Err(Error::Malformed(
                    "some rows carry column indices and others do not",
                ));
            }

            let count = u32::try_from(row.points.len())
                .map_err(|_| Error::Malformed("a row holds more points than a u32 can count"))?;

            if let Some(row_columns) = &row.columns {
                if row_columns.len() != row.points.len() {
                    return Err(Error::Malformed(
                        "a row holds a different number of column indices than points",
                    ));
                }
                columns.extend_from_slice(row_columns);
            }

            ordinals.push(row.ordinal);
            counts.push(count);
        }

        let layout = RowLayout::new(ordinals, counts);
        Ok(if has_columns == Some(true) {
            layout.with_columns(columns)
        } else {
            layout
        })
    }
}

impl Named for RowPoints3 {
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

/// Everything [`write_to_with`] can be told beyond the geometry and the tolerance.
///
/// This type is non-exhaustive so later settings do not break callers. Build it with
/// [`WriteOptions::new`] and the `with_` methods, or use [`Default`].
#[non_exhaustive]
#[derive(Debug, Clone, Default)]
pub struct WriteOptions {
    /// File-level metadata. Stored, never interpreted.
    pub metadata: Metadata,
    /// How thoroughly to search for a smaller file. This affects file size and encoding time only.
    pub effort: Effort,
}

impl WriteOptions {
    /// The defaults: no metadata, [`Effort::Balanced`].
    pub fn new() -> Self {
        Self::default()
    }

    /// The same options carrying file-level metadata.
    pub fn with_metadata(mut self, metadata: Metadata) -> Self {
        self.metadata = metadata;
        self
    }

    /// The same options at a different search effort.
    pub fn with_effort(mut self, effort: Effort) -> Self {
        self.effort = effort;
        self
    }
}

/// Write a collection of row-organized point sets, every item at the same storage tolerance.
///
/// # Errors
///
/// Everything [`RowPoints3::validate`] can return, plus [`Error::ToleranceNotRepresentable`] if any
/// axis is too wide to meet `tol`, and [`Error::Malformed`] for a non-finite coordinate.
pub fn write_to<W: Write>(writer: &mut W, items: &[RowPoints3], tol: f64) -> Result<()> {
    write_to_with(writer, items, tol, &WriteOptions::default())
}

/// Write a collection with file-level metadata attached. See [`write_to_with`].
pub fn write_to_with_meta<W: Write>(
    writer: &mut W,
    items: &[RowPoints3],
    tol: f64,
    file_metadata: &Metadata,
) -> Result<()> {
    let options = WriteOptions::new().with_metadata(file_metadata.clone());
    write_to_with(writer, items, tol, &options)
}

/// Write a collection with full control over metadata and effort. See [`write_to`].
pub fn write_to_with<W: Write>(
    writer: &mut W,
    items: &[RowPoints3],
    tol: f64,
    options: &WriteOptions,
) -> Result<()> {
    let count = u32::try_from(items.len())
        .map_err(|_| Error::Malformed("container holds more items than a u32 can count"))?;
    container::write_header(writer, Kind::RowPoints3, count, &options.metadata)?;

    for item in items {
        write_item(writer, item, tol, options.effort)?;
    }

    Ok(())
}

/// Write a single item as a one-item collection.
pub fn write_one_to<W: Write>(writer: &mut W, item: &RowPoints3, tol: f64) -> Result<()> {
    write_to(writer, std::slice::from_ref(item), tol)
}

/// Read a collection of row-organized point sets.
///
/// # Errors
///
/// [`Error::Malformed`] if the file holds a different kind, or if its row block and points block
/// disagree about how many points there are.
pub fn read_from<R: Read>(reader: &mut R) -> Result<Vec<RowPoints3>> {
    let header = container::read_header(reader, Kind::RowPoints3)?;

    let mut out = Vec::with_capacity((header.count as usize).min(MAX_PREALLOC));
    for _ in 0..header.count {
        out.push(read_item(reader)?);
    }

    Ok(out)
}

/// Read a container that holds exactly one item.
///
/// # Errors
///
/// [`Error::NotASingleItem`] if the container holds any other number.
pub fn read_one_from<R: Read>(reader: &mut R) -> Result<RowPoints3> {
    let header = container::read_header(reader, Kind::RowPoints3)?;
    if header.count != 1 {
        return Err(Error::NotASingleItem {
            found: header.count,
        });
    }
    read_item(reader)
}

/// Write a collection to a file. See [`write_to`].
pub fn write_file(path: &Path, items: &[RowPoints3], tol: f64) -> Result<()> {
    write_file_with(path, items, tol, &WriteOptions::default())
}

/// Write a collection to a file with full control over metadata and effort.
pub fn write_file_with(
    path: &Path,
    items: &[RowPoints3],
    tol: f64,
    options: &WriteOptions,
) -> Result<()> {
    // Buffering is required because the bit reader and writer process one byte at a time. An
    // unbuffered file would cause one system call per byte.
    let mut writer = BufWriter::new(File::create(path)?);
    write_to_with(&mut writer, items, tol, options)?;
    writer.flush()?;
    Ok(())
}

/// Write a single item to a file as a one-item collection.
pub fn write_one_file(path: &Path, item: &RowPoints3, tol: f64) -> Result<()> {
    write_file(path, std::slice::from_ref(item), tol)
}

/// Read a collection from a file.
pub fn read_file(path: &Path) -> Result<Vec<RowPoints3>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_from(&mut reader)
}

/// Read a file holding exactly one item. See [`read_one_from`].
pub fn read_one_file(path: &Path) -> Result<RowPoints3> {
    let mut reader = BufReader::new(File::open(path)?);
    read_one_from(&mut reader)
}

fn write_item<W: Write>(writer: &mut W, scan: &RowPoints3, tol: f64, effort: Effort) -> Result<()> {
    let layout = scan.layout()?;
    item::write_preamble(writer, scan.name.as_deref(), &scan.metadata, false)?;
    write_rows(writer, &layout)?;

    // One points block covers all rows, allowing the partitioner to cut according to geometry
    // without being constrained by row boundaries.
    let points = scan
        .rows
        .iter()
        .flat_map(|row| row.points.iter().copied())
        .collect::<Vec<_>>();
    write_points_with(writer, &points, tol, effort)?;

    Ok(())
}

fn read_item<R: Read>(reader: &mut R) -> Result<RowPoints3> {
    let preamble = item::read_preamble(reader, false)?;
    let layout = read_rows(reader)?;
    let points: Vec<[f64; 3]> = read_points(reader)?;

    if layout.point_count() != points.len() as u64 {
        return Err(Error::Malformed(
            "row block and points block disagree about how many points there are",
        ));
    }

    let mut rows = Vec::with_capacity(layout.row_count().min(MAX_PREALLOC));
    let mut at = 0usize;
    for (i, &count) in layout.counts.iter().enumerate() {
        let count = count as usize;
        let row_points = points[at..at + count].to_vec();
        let columns = layout
            .columns
            .as_ref()
            .map(|columns| columns[at..at + count].to_vec());

        rows.push(PointRow3 {
            ordinal: layout.ordinals[i],
            points: row_points,
            columns,
        });
        at += count;
    }

    Ok(RowPoints3 {
        name: preamble.name,
        rows,
        metadata: preamble.metadata,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::probe;
    use crate::metadata::Value;
    use crate::testgen::Rng;
    use crate::{Cloud3, cloud};
    use std::io::Cursor;

    const TOL: f64 = 1e-4;

    fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    fn assert_matches(original: &RowPoints3, recovered: &RowPoints3, tol: f64, what: &str) {
        assert_eq!(recovered.name, original.name, "{what}: name");
        assert_eq!(recovered.metadata, original.metadata, "{what}: metadata");
        assert_eq!(
            recovered.rows.len(),
            original.rows.len(),
            "{what}: row count"
        );

        for (i, (o, r)) in original.rows.iter().zip(recovered.rows.iter()).enumerate() {
            assert_eq!(r.ordinal, o.ordinal, "{what}: row {i} ordinal");
            assert_eq!(r.columns, o.columns, "{what}: row {i} columns");
            assert_eq!(
                r.points.len(),
                o.points.len(),
                "{what}: row {i} point count"
            );

            for (j, (op, rp)) in o.points.iter().zip(r.points.iter()).enumerate() {
                let d = distance(op, rp);
                assert!(d <= tol, "{what}: row {i} point {j} recovered {d} away");
            }
        }
    }

    fn round_trip(scan: &RowPoints3, tol: f64) -> RowPoints3 {
        let mut buf = Vec::new();
        write_one_to(&mut buf, scan, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_one_from(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder left bytes unread"
        );
        back
    }

    /// Model the row geometry of a real snapshot sensor. Each row represents a plane through the
    /// camera, so its y-coordinate varies with z instead of forming a constant-y line.
    fn jittered_grid(rows: usize, per_row: usize, seed: u64) -> RowPoints3 {
        let mut rng = Rng::new(seed);
        let mut out = Vec::with_capacity(rows);

        for i in 0..rows {
            let mut points = Vec::with_capacity(per_row);
            for j in 0..per_row {
                let z = 10.0 + 2.0 * rng.next_f64();
                let x = j as f64 * 0.05 + 0.002 * rng.next_f64();
                let y = i as f64 * 0.1 + 0.3 * z;
                points.push([x, y, z]);
            }
            out.push(PointRow3::new(i as u32, points));
        }

        RowPoints3::new(out)
    }

    #[test]
    fn an_empty_item_round_trips() {
        let scan = RowPoints3::new(Vec::new()).named("nothing");
        assert_matches(&scan, &round_trip(&scan, TOL), TOL, "empty");
    }

    #[test]
    fn a_single_row_round_trips() {
        let scan = RowPoints3::new(vec![PointRow3::new(
            3,
            vec![[0.0, 0.0, 0.0], [1.0, 0.5, 0.25]],
        )]);
        assert_matches(&scan, &round_trip(&scan, TOL), TOL, "single row");
    }

    #[test]
    fn a_dense_jittered_grid_round_trips_within_tolerance() {
        let scan = jittered_grid(120, 200, 31);
        assert_matches(&scan, &round_trip(&scan, TOL), TOL, "grid");
    }

    #[test]
    fn ordinal_gaps_survive_the_round_trip() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(0, vec![[0.0, 0.0, 0.0]]),
            PointRow3::new(4, vec![[0.0, 0.4, 0.0]]),
            PointRow3::new(900, vec![[0.0, 90.0, 0.0]]),
        ]);
        let back = round_trip(&scan, TOL);
        assert_eq!(
            back.rows.iter().map(|r| r.ordinal).collect::<Vec<_>>(),
            vec![0, 4, 900]
        );
    }

    #[test]
    fn columns_round_trip_when_present() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(0, vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]).with_columns(vec![0, 7]),
            PointRow3::new(1, vec![[0.0, 0.1, 0.0]]).with_columns(vec![3]),
        ]);
        assert_matches(&scan, &round_trip(&scan, TOL), TOL, "columns");
    }

    /// Row grouping should add only a fraction of a percent to the cost of storing the same points
    /// as an unordered cloud.
    #[test]
    fn the_grouping_costs_almost_nothing_against_a_cloud() {
        let scan = jittered_grid(120, 200, 5);
        let mut rows_buf = Vec::new();
        write_one_to(&mut rows_buf, &scan, TOL).unwrap();

        let points = scan
            .rows
            .iter()
            .flat_map(|row| row.points.iter().copied())
            .collect::<Vec<_>>();
        let mut cloud_buf = Vec::new();
        cloud::write_one_to(&mut cloud_buf, &Cloud3::new(points), TOL).unwrap();

        let overhead = (rows_buf.len() as f64 - cloud_buf.len() as f64) / cloud_buf.len() as f64;
        assert!(
            overhead < 0.01,
            "row grouping cost {:.3}% over the same points as a cloud",
            overhead * 100.0
        );
    }

    #[test]
    fn a_row_with_no_points_is_refused() {
        let scan = RowPoints3::new(vec![PointRow3::new(0, Vec::new())]);
        assert!(matches!(
            write_one_to(&mut Vec::new(), &scan, TOL),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn ordinals_that_do_not_increase_are_refused() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(2, vec![[0.0; 3]]),
            PointRow3::new(2, vec![[0.0; 3]]),
        ]);
        assert!(matches!(
            write_one_to(&mut Vec::new(), &scan, TOL),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn mixing_rows_with_and_without_columns_is_refused() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(0, vec![[0.0; 3]]).with_columns(vec![0]),
            PointRow3::new(1, vec![[0.0; 3]]),
        ]);
        assert!(matches!(
            write_one_to(&mut Vec::new(), &scan, TOL),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn a_column_count_that_disagrees_with_a_row_is_refused() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(0, vec![[0.0; 3], [1.0; 3]]).with_columns(vec![0]),
        ]);
        assert!(matches!(
            write_one_to(&mut Vec::new(), &scan, TOL),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn columns_that_do_not_increase_within_a_row_are_refused() {
        let scan = RowPoints3::new(vec![
            PointRow3::new(0, vec![[0.0; 3], [1.0; 3]]).with_columns(vec![4, 4]),
        ]);
        assert!(matches!(
            write_one_to(&mut Vec::new(), &scan, TOL),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn probe_reports_the_kind_without_decoding_geometry() {
        let scan = jittered_grid(4, 8, 2);
        let mut buf = Vec::new();
        write_one_to(&mut buf, &scan, TOL).unwrap();

        let header = probe(&mut buf.as_slice()).unwrap();
        assert_eq!(header.kind, Kind::RowPoints3);
        assert_eq!(header.kind.extension(), "tcrpf3");
        assert_eq!(header.count, 1);
    }

    #[test]
    fn reading_row_points_as_a_cloud_is_rejected() {
        let scan = jittered_grid(4, 8, 3);
        let mut buf = Vec::new();
        write_one_to(&mut buf, &scan, TOL).unwrap();

        assert!(matches!(
            cloud::read_one_from::<_, 3>(&mut buf.as_slice()),
            Err(Error::Malformed(_))
        ));
    }

    /// The closed flag applies only to polylines. Reject it here instead of accepting an undefined
    /// flag bit.
    #[test]
    fn a_closed_flag_is_rejected() {
        let scan = jittered_grid(2, 4, 4);
        let mut buf = Vec::new();
        write_one_to(&mut buf, &scan, TOL).unwrap();

        // The item preamble's flag byte follows the 12-byte container header.
        buf[12] |= item::CLOSED;
        assert!(matches!(
            read_one_from(&mut buf.as_slice()),
            Err(Error::Malformed(_))
        ));
    }

    /// Row counts and point data are stored independently. Reject a file when the two blocks
    /// disagree instead of slicing rows from the wrong points.
    #[test]
    fn a_row_and_point_count_mismatch_is_refused() {
        let scan = RowPoints3::new(vec![PointRow3::new(0, vec![[0.0; 3], [1.0; 3]])]);
        let layout = RowLayout::new(vec![0], vec![3]);

        let mut buf = Vec::new();
        container::write_header(&mut buf, Kind::RowPoints3, 1, &Metadata::new()).unwrap();
        item::write_preamble(&mut buf, None, &Metadata::new(), false).unwrap();
        write_rows(&mut buf, &layout).unwrap();
        let points = scan.rows[0].points.clone();
        write_points_with(&mut buf, &points, TOL, Effort::default()).unwrap();

        assert!(matches!(
            read_one_from(&mut buf.as_slice()),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn metadata_round_trips() {
        let scan = jittered_grid(6, 10, 9)
            .named("pass 1")
            .with_meta("engeom.row_pitch", 0.1)
            .with_meta("engeom.source_frame", "gocator-3x00 left-handed, x negated");

        let back = round_trip(&scan, TOL);
        assert_eq!(back.metadata, scan.metadata);
        assert_eq!(
            back.metadata.get("engeom.row_pitch"),
            Some(&Value::F64(0.1))
        );
    }

    #[test]
    fn a_collection_keeps_its_order_and_names() {
        let items = vec![
            jittered_grid(4, 6, 11).named("500"),
            jittered_grid(5, 6, 12).named("1000"),
            jittered_grid(3, 6, 13),
        ];

        let mut buf = Vec::new();
        write_to(&mut buf, &items, TOL).unwrap();
        let back = read_from(&mut buf.as_slice()).unwrap();

        assert_eq!(back.len(), 3);
        for (o, r) in items.iter().zip(back.iter()) {
            assert_matches(o, r, TOL, "collection");
        }
        assert!(crate::find_by_name(&back, "1000").is_some());
    }

    #[test]
    fn read_one_refuses_a_multi_item_file() {
        let items = vec![jittered_grid(3, 4, 14), jittered_grid(3, 4, 15)];
        let mut buf = Vec::new();
        write_to(&mut buf, &items, TOL).unwrap();

        assert!(matches!(
            read_one_from(&mut buf.as_slice()),
            Err(Error::NotASingleItem { found: 2 })
        ));
    }

    #[test]
    fn a_file_round_trips_through_the_disk() {
        let scan = jittered_grid(20, 40, 16).named("disk");
        let dir = std::env::temp_dir().join("tol-compress-row-points");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("scan.tcrpf3");

        write_one_file(&path, &scan, TOL).unwrap();
        let back = read_one_file(&path).unwrap();
        assert_matches(&scan, &back, TOL, "file");

        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn truncation_is_an_error_rather_than_a_partial_read() {
        let scan = jittered_grid(8, 12, 17);
        let mut buf = Vec::new();
        write_one_to(&mut buf, &scan, TOL).unwrap();

        for cut in [8, 16, 24, buf.len() / 2, buf.len() - 1] {
            assert!(
                read_one_from(&mut &buf[..cut]).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
    }

    /// Effort changes the partitioner's search depth without changing decoded data. A higher effort
    /// must not produce a larger file than the least expensive setting.
    #[test]
    fn effort_never_costs_bytes() {
        let scan = jittered_grid(40, 60, 18);

        let mut sizes = Vec::new();
        for effort in [Effort::Quick, Effort::Balanced, Effort::Thorough] {
            let mut buf = Vec::new();
            let options = WriteOptions::new().with_effort(effort);
            write_to_with(&mut buf, std::slice::from_ref(&scan), TOL, &options).unwrap();

            let back = read_one_from(&mut buf.as_slice()).unwrap();
            assert_matches(&scan, &back, TOL, "effort");
            sizes.push(buf.len());
        }

        assert!(sizes[1] <= sizes[0], "balanced {sizes:?}");
        assert!(sizes[2] <= sizes[1], "thorough {sizes:?}");
    }
}
