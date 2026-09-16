//! The row block: how points are grouped into ordered rows.
//!
//! A rasterizing sensor produces points organized into rows. This organization makes the points
//! meshable because two neighboring rows can be triangulated into a strip without surface
//! reconstruction. Geometry remains encoded through [`crate::points`]; this block records only the
//! integer grouping information.
//!
//! # Block layout
//!
//! ```text
//! u32     row count R; a count of zero ends the block here, before the flag byte
//! u8      flags: bit0 has columns
//! stream  R ordinal codes: code[0] = ordinal[0], code[i] = ordinal[i] - ordinal[i-1] - 1
//! stream  R point counts, each at least 1
//! stream  P column codes, only with bit0 set, where P is the sum of the counts:
//!         per row, the first code is column[0], then column[j] - column[j-1] - 1
//! ```
//!
//! Each `stream` is a [`crate::blocks`] stream. Consecutive rows and dense-grid columns encode as
//! runs of zeros, and a block of 64 zeros requires six bits. On a full raster, the complete row
//! block uses a fraction of a bit per point.
//!
//! # Ordinals and columns
//!
//! The **ordinal** records a row's position in the sensor sweep, independent of its position in the
//! file. A mesher uses it to determine the distance between rows. Removing an empty or low-quality
//! row leaves an ordinal gap and does not pull neighboring rows together. The delta subtracts one
//! so consecutive ordinals encode as zero and a gap uses only the bits required for its size.
//!
//! **Columns** are optional because some sensors do not provide them. A laser-line sweep identifies
//! the profile containing each point but has no meaningful index across the profile. A
//! time-of-flight or structured-light sensor provides both indices. Column indices let a reader
//! rasterize the points along the other axis, while rows remain meshable without them.
//!
//! Skipping an ordinal represents an empty row. The format rejects a zero point count to maintain
//! one representation for "nothing was measured here."

use crate::blocks;
use crate::error::{Error, Result};
use crate::raw::{MAX_PREALLOC, read_u8, read_u32, write_u8, write_u32};
use std::io::{Read, Write};

/// Column indices are present for every row.
pub const HAS_COLUMNS: u8 = 1 << 0;

/// Every flag bit this version understands.
const KNOWN_FLAGS: u8 = HAS_COLUMNS;

/// How a run of points is divided into rows.
///
/// `ordinals` and `counts` contain one parallel entry per row. When present, `columns` contains one
/// entry per point in row order.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RowLayout {
    /// Each row's position in the sensor's sweep, strictly increasing.
    pub ordinals: Vec<u32>,
    /// How many points each row holds, each at least 1.
    pub counts: Vec<u32>,
    /// Each point's column index, in row order, strictly increasing within a row. `None` when the
    /// source has no meaningful column index.
    pub columns: Option<Vec<u32>>,
}

impl RowLayout {
    /// A layout with no column indices.
    pub fn new(ordinals: Vec<u32>, counts: Vec<u32>) -> Self {
        Self {
            ordinals,
            counts,
            columns: None,
        }
    }

    /// The same layout carrying a column index for every point.
    pub fn with_columns(mut self, columns: Vec<u32>) -> Self {
        self.columns = Some(columns);
        self
    }

    /// How many rows the layout describes.
    pub fn row_count(&self) -> usize {
        self.ordinals.len()
    }

    /// How many points the rows hold in total.
    pub fn point_count(&self) -> u64 {
        self.counts.iter().map(|&c| u64::from(c)).sum()
    }

    /// Validate all format requirements for a layout.
    ///
    /// [`write_rows`] calls this method automatically. Call it directly only to validate before
    /// writing.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Malformed`] with a distinct message for mismatched ordinal and count lengths,
    /// an empty row, non-increasing ordinals, a column-vector length that differs from the point
    /// count, or non-increasing columns within a row.
    pub fn validate(&self) -> Result<()> {
        if self.ordinals.len() != self.counts.len() {
            return Err(Error::Malformed(
                "row layout has a different number of ordinals than counts",
            ));
        }

        for (i, &ordinal) in self.ordinals.iter().enumerate() {
            if i > 0 && ordinal <= self.ordinals[i - 1] {
                return Err(Error::Malformed("row ordinals do not strictly increase"));
            }
        }

        if self.counts.contains(&0) {
            return Err(Error::Malformed(
                "a row holds no points; express an empty row by skipping its ordinal",
            ));
        }

        if let Some(columns) = &self.columns {
            if u64::try_from(columns.len()) != Ok(self.point_count()) {
                return Err(Error::Malformed(
                    "row layout holds a different number of column indices than points",
                ));
            }

            let mut at = 0usize;
            for &count in &self.counts {
                let row = &columns[at..at + count as usize];
                for j in 1..row.len() {
                    if row[j] <= row[j - 1] {
                        return Err(Error::Malformed(
                            "column indices do not strictly increase within a row",
                        ));
                    }
                }
                at += count as usize;
            }
        }

        Ok(())
    }
}

/// Write a row block.
///
/// # Errors
///
/// Everything [`RowLayout::validate`] can return, plus [`Error::Malformed`] if there are more rows
/// than a `u32` can count.
pub fn write_rows<W: Write>(writer: &mut W, layout: &RowLayout) -> Result<()> {
    layout.validate()?;

    let count = u32::try_from(layout.ordinals.len())
        .map_err(|_| Error::Malformed("block holds more rows than a u32 can count"))?;
    write_u32(writer, count)?;

    if count == 0 {
        return Ok(());
    }

    let flags = if layout.columns.is_some() {
        HAS_COLUMNS
    } else {
        0
    };
    write_u8(writer, flags)?;

    let mut ordinal_codes = Vec::with_capacity(layout.ordinals.len());
    for (i, &ordinal) in layout.ordinals.iter().enumerate() {
        let code = if i == 0 {
            u64::from(ordinal)
        } else {
            u64::from(ordinal - layout.ordinals[i - 1] - 1)
        };
        ordinal_codes.push(code);
    }
    blocks::write_stream(writer, &ordinal_codes)?;

    let counts = layout
        .counts
        .iter()
        .map(|&c| u64::from(c))
        .collect::<Vec<_>>();
    blocks::write_stream(writer, &counts)?;

    if let Some(columns) = &layout.columns {
        let mut column_codes = Vec::with_capacity(columns.len());
        let mut at = 0usize;
        for &count in &layout.counts {
            let row = &columns[at..at + count as usize];
            column_codes.push(u64::from(row[0]));
            for j in 1..row.len() {
                column_codes.push(u64::from(row[j] - row[j - 1] - 1));
            }
            at += count as usize;
        }
        blocks::write_stream(writer, &column_codes)?;
    }

    Ok(())
}

/// Read a row block written by [`write_rows`].
///
/// # Errors
///
/// [`Error::Malformed`] for unknown flag bits, a row holding no points, an ordinal or column index
/// that overflows a `u32`, or a total point count that overflows a `u32`. I/O errors propagate for
/// truncated input.
pub fn read_rows<R: Read>(reader: &mut R) -> Result<RowLayout> {
    let count = read_u32(reader)? as usize;
    if count == 0 {
        return Ok(RowLayout::default());
    }

    let flags = read_u8(reader)?;
    if flags & !KNOWN_FLAGS != 0 {
        return Err(Error::Malformed("row block sets unknown flag bits"));
    }

    let ordinal_codes = blocks::read_stream(reader, count)?;
    let mut ordinals = Vec::with_capacity(count.min(MAX_PREALLOC));
    let mut previous: Option<u32> = None;
    for code in ordinal_codes {
        let ordinal = match previous {
            None => u32::try_from(code).map_err(|_| Error::Malformed("row ordinal overflows"))?,
            Some(prev) => {
                let step = u32::try_from(code)
                    .ok()
                    .and_then(|c| c.checked_add(1))
                    .ok_or(Error::Malformed("row ordinal overflows"))?;
                prev.checked_add(step)
                    .ok_or(Error::Malformed("row ordinal overflows"))?
            }
        };
        previous = Some(ordinal);
        ordinals.push(ordinal);
    }

    let count_codes = blocks::read_stream(reader, count)?;
    let mut counts = Vec::with_capacity(count.min(MAX_PREALLOC));
    let mut total = 0u64;
    for code in count_codes {
        let n = u32::try_from(code).map_err(|_| Error::Malformed("row point count overflows"))?;
        if n == 0 {
            return Err(Error::Malformed("row block declares a row with no points"));
        }
        total += u64::from(n);
        counts.push(n);
    }

    let total = usize::try_from(total)
        .ok()
        .filter(|&t| u32::try_from(t).is_ok())
        .ok_or(Error::Malformed(
            "row block declares more points than a u32 can count",
        ))?;

    let columns = if flags & HAS_COLUMNS != 0 {
        let column_codes = blocks::read_stream(reader, total)?;
        let mut columns: Vec<u32> = Vec::with_capacity(total.min(MAX_PREALLOC));
        let mut at = 0usize;
        for &n in &counts {
            for j in 0..n as usize {
                let code = column_codes[at + j];
                let column = if j == 0 {
                    u32::try_from(code).map_err(|_| Error::Malformed("column index overflows"))?
                } else {
                    let step = u32::try_from(code)
                        .ok()
                        .and_then(|c| c.checked_add(1))
                        .ok_or(Error::Malformed("column index overflows"))?;
                    columns[at + j - 1]
                        .checked_add(step)
                        .ok_or(Error::Malformed("column index overflows"))?
                };
                columns.push(column);
            }
            at += n as usize;
        }
        Some(columns)
    } else {
        None
    };

    Ok(RowLayout {
        ordinals,
        counts,
        columns,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn round_trip(layout: &RowLayout) -> RowLayout {
        let mut buf = Vec::new();
        write_rows(&mut buf, layout).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_rows(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder left bytes unread"
        );
        back
    }

    fn dense(rows: u32, per_row: u32) -> RowLayout {
        let ordinals = (0..rows).collect::<Vec<_>>();
        let counts = vec![per_row; rows as usize];
        let columns = (0..rows).flat_map(|_| 0..per_row).collect::<Vec<_>>();
        RowLayout::new(ordinals, counts).with_columns(columns)
    }

    #[test]
    fn an_empty_layout_is_four_bytes() {
        let mut buf = Vec::new();
        write_rows(&mut buf, &RowLayout::default()).unwrap();
        assert_eq!(buf, vec![0, 0, 0, 0]);

        assert_eq!(
            read_rows(&mut Cursor::new(&buf)).unwrap(),
            RowLayout::default()
        );
    }

    #[test]
    fn rows_round_trip_without_columns() {
        let layout = RowLayout::new(vec![0, 1, 2, 5, 9], vec![3, 1, 7, 2, 4]);
        assert_eq!(round_trip(&layout), layout);
    }

    #[test]
    fn rows_round_trip_with_columns() {
        let layout =
            RowLayout::new(vec![4, 5, 7], vec![2, 3, 1]).with_columns(vec![0, 9, 1, 2, 40, 17]);
        assert_eq!(round_trip(&layout), layout);
    }

    /// Delta coding makes dense-raster grouping nearly free compared with storing the same points
    /// as an unordered cloud.
    #[test]
    fn a_dense_raster_costs_only_width_headers() {
        let layout = dense(200, 640);
        let mut buf = Vec::new();
        write_rows(&mut buf, &layout).unwrap();

        let points = layout.point_count();
        let bits_per_point = (buf.len() as f64 * 8.0) / points as f64;
        assert!(
            bits_per_point < 0.2,
            "{bits_per_point} bits per point for a dense raster"
        );
        assert_eq!(round_trip(&layout), layout);
    }

    #[test]
    fn ordinal_gaps_round_trip_to_the_top_of_the_range() {
        let layout = RowLayout::new(vec![0, 1_000_000, u32::MAX], vec![1, 1, 1]);
        assert_eq!(round_trip(&layout), layout);
    }

    #[test]
    fn a_single_row_round_trips() {
        let layout = RowLayout::new(vec![7], vec![4]).with_columns(vec![10, 11, 12, 13]);
        assert_eq!(round_trip(&layout), layout);
    }

    #[test]
    fn a_row_with_no_points_is_refused() {
        let layout = RowLayout::new(vec![0, 1], vec![3, 0]);
        assert!(matches!(
            write_rows(&mut Vec::new(), &layout),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn ordinals_that_do_not_increase_are_refused() {
        for ordinals in [vec![0, 0], vec![5, 4]] {
            let layout = RowLayout::new(ordinals, vec![1, 1]);
            assert!(matches!(
                write_rows(&mut Vec::new(), &layout),
                Err(Error::Malformed(_))
            ));
        }
    }

    #[test]
    fn a_column_count_that_disagrees_with_the_points_is_refused() {
        let layout = RowLayout::new(vec![0, 1], vec![2, 2]).with_columns(vec![0, 1, 0]);
        assert!(matches!(
            write_rows(&mut Vec::new(), &layout),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn columns_that_do_not_increase_within_a_row_are_refused() {
        let layout = RowLayout::new(vec![0], vec![3]).with_columns(vec![0, 5, 5]);
        assert!(matches!(
            write_rows(&mut Vec::new(), &layout),
            Err(Error::Malformed(_))
        ));
    }

    /// Column indices restart for each row, so a row can begin left of the previous row's endpoint.
    #[test]
    fn columns_may_restart_lower_in_the_next_row() {
        let layout = RowLayout::new(vec![0, 1], vec![2, 2]).with_columns(vec![10, 11, 0, 1]);
        assert_eq!(round_trip(&layout), layout);
    }

    #[test]
    fn a_mismatched_ordinal_and_count_length_is_refused() {
        let layout = RowLayout::new(vec![0, 1], vec![1]);
        assert!(matches!(
            write_rows(&mut Vec::new(), &layout),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn unknown_flag_bits_are_rejected() {
        let layout = RowLayout::new(vec![0], vec![1]);
        let mut buf = Vec::new();
        write_rows(&mut buf, &layout).unwrap();
        buf[4] = 1 << 5;

        assert!(matches!(
            read_rows(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    /// The `u32` row count is attacker-controlled. A claim of four billion rows must fail while
    /// reading the following stream without allocating storage for all claimed rows.
    #[test]
    fn an_absurd_row_count_does_not_allocate() {
        let mut buf = Vec::new();
        write_u32(&mut buf, u32::MAX).unwrap();
        write_u8(&mut buf, 0).unwrap();
        blocks::write_stream(&mut buf, &[0, 0, 0]).unwrap();

        assert!(read_rows(&mut Cursor::new(&buf)).is_err());
    }

    #[test]
    fn truncation_at_every_offset_is_an_error() {
        let layout = dense(4, 6);
        let mut buf = Vec::new();
        write_rows(&mut buf, &layout).unwrap();

        for cut in 1..buf.len() {
            assert!(
                read_rows(&mut Cursor::new(&buf[..cut])).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
        assert!(read_rows(&mut Cursor::new(&buf)).is_ok());
    }
}
