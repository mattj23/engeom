//! Row-organized point serialization for the practical tol-compression format.
//!
//! This adapter converts between engeom's [`RowPointsScan3`] and the [`tol_compress`] crate's
//! `RowPoints3` container. The `tol` parameter is the maximum permitted round-trip position error
//! for any point, in the coordinate units.
//!
//! The recommended format extension is `.tcrpf3`.
//!
//! # What this format is for
//!
//! A rastering or sweeping sensor produces points grouped into rows. This grouping permits meshing
//! by triangulating neighboring rows into strips without surface reconstruction. Writing the points
//! as a point cloud discards the grouping. This format preserves it, and
//! [`load_tc_row_points_mesh_data`] uses the same strip mesher as the LPTF3 loader.
//!
//! # What travels as metadata
//!
//! The format stores only geometry and grouping. Information that engeom needs to interpret a scan
//! is stored in the item's metadata map, which `tol-compress` preserves without interpretation:
//!
//! - [`ROW_PITCH`], required: the nominal distance between consecutive row ordinals. A row's strip
//!   coordinate is its ordinal multiplied by this value. It is not derived from a point because a
//!   snapshot-sensor raster row is a plane through the camera, and its points differ in position by
//!   millimeters.
//! - [`COL_PITCH`], optional: the nominal spacing between points along a row, recorded for callers
//!   who need it.
//! - [`ALONG_AXIS`], optional, `"x"` (the default) or `"y"`: the world coordinate that varies along
//!   a row. A sensor whose strips run along y sets this value. Transposing coordinates during
//!   conversion would reflect the points out of the sensor's frame.
//! - [`SOURCE_FRAME`], optional free text: the transformation used to place points in a right-handed
//!   frame. Many sensors report in a left-handed frame, and a converter that flips an axis should
//!   record that operation here.
//!
//! # Winding
//!
//! The mesher gives faces a +z normal when rows ascend in the strip coordinate and the frame is
//! right-handed with +z toward the sensor. A converter from a left-handed source must negate an
//! axis and number its ordinals in ascending sweep order. Otherwise, the resulting mesh is inside
//! out.
//!
//! # Attributes
//!
//! The format stores geometry only. Writing a scan with additional attributes returns an error that
//! names them, preventing silent data loss before the format gains attribute support.

use crate::io::Lptf3Load;
use crate::io::row_scan::{
    RowAxis, RowScan3, compute_row_scan_mesh, smooth_row_scan, thin_row_scan,
};
use crate::{MeshData3, Point3, PointCloud3, Result};
use std::io::{Read, Write};
use std::path::Path;
use tol_compress::{Metadata, PointRow3 as TcRow, RowPoints3 as TcScan, row_points as tc};

/// The metadata key a scan's nominal row pitch is stored under.
///
/// The key is namespaced because other file writers share the metadata map and could also use an
/// unprefixed `row_pitch` key.
pub const ROW_PITCH: &str = "engeom.row_pitch";

/// The metadata key a scan's nominal column pitch is stored under, when it has one.
pub const COL_PITCH: &str = "engeom.col_pitch";

/// The metadata key naming which world coordinate varies along a row, `"x"` or `"y"`.
pub const ALONG_AXIS: &str = "engeom.along_axis";

/// The metadata key describing the frame the points came from and what was done to it.
pub const SOURCE_FRAME: &str = "engeom.source_frame";

/// One row of a scan, as engeom holds it.
///
/// The `ordinal` records the row's position in the sensor sweep, independent of its position in the
/// file. Dropping an empty row therefore leaves a gap without pulling neighboring rows together.
/// Rows must appear in increasing ordinal order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointRow3 {
    /// Where this row sits in the sensor's sweep.
    pub ordinal: u32,
    /// The row's points.
    pub points: Vec<Point3>,
    /// Each point's index across the row, strictly increasing, when the sensor has one. Either
    /// every row carries these or none does.
    pub columns: Option<Vec<u32>>,
}

impl PointRow3 {
    /// A row with no column indices.
    pub fn new(ordinal: u32, points: Vec<Point3>) -> Self {
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

/// A scan whose points are grouped into the rows the sensor produced them in.
#[derive(Debug, Clone, PartialEq)]
pub struct RowPointsScan3 {
    /// The rows, in increasing ordinal order.
    pub rows: Vec<PointRow3>,
    /// The nominal distance between consecutive row ordinals.
    pub row_pitch: f64,
    /// The nominal spacing between points along a row, when it is known.
    pub col_pitch: Option<f64>,
    /// Which world coordinate varies along a row.
    pub along: RowAxis,
    /// Optional identifier, preserved through a round trip.
    pub name: Option<String>,
    /// Additional caller data. It is stored without interpretation and merged with this module's
    /// keys when the scan is written.
    pub metadata: Metadata,
}

impl RowPointsScan3 {
    /// A scan of rows running along world x, with no name and no extra metadata.
    pub fn new(rows: Vec<PointRow3>, row_pitch: f64) -> Self {
        Self {
            rows,
            row_pitch,
            col_pitch: None,
            along: RowAxis::X,
            name: None,
            metadata: Metadata::new(),
        }
    }

    /// The same scan carrying a nominal column pitch.
    pub fn with_col_pitch(mut self, col_pitch: f64) -> Self {
        self.col_pitch = Some(col_pitch);
        self
    }

    /// The same scan with its rows running along a different world axis.
    pub fn with_along_axis(mut self, along: RowAxis) -> Self {
        self.along = along;
        self
    }

    /// The same scan carrying a name.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The same scan carrying one extra metadata entry.
    pub fn with_meta(
        mut self,
        key: impl Into<String>,
        value: impl Into<tol_compress::Value>,
    ) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }

    /// How many points the rows hold in total.
    pub fn point_count(&self) -> usize {
        self.rows.iter().map(|r| r.points.len()).sum()
    }

    /// Convert the scan to the shared row intermediate for thinning, meshing, or flattening.
    ///
    /// This conversion sorts each row by the along coordinate instead of trusting file order. Each
    /// row's strip coordinate is its ordinal multiplied by the row pitch.
    ///
    /// # Errors
    ///
    /// An error if the row pitch is not a positive distance.
    pub fn to_row_scan(&self) -> Result<RowScan3> {
        let rows = self
            .rows
            .iter()
            .map(|r| r.points.clone())
            .collect::<Vec<_>>();
        let row_y = self
            .rows
            .iter()
            .map(|r| r.ordinal as f64 * self.row_pitch)
            .collect::<Vec<_>>();

        RowScan3::new(rows, row_y, self.row_pitch, 1, self.along)
    }
}

/// Serialize a scan into the tcrpf3 format, writing to any [`Write`] sink.
///
/// Point positions are quantized to the narrowest per-axis bit width that keeps every point within
/// `tol` of its original position. Row grouping and ordinals are stored exactly. Smaller `tol`
/// values increase accuracy and require more bytes per point.
///
/// Use [`write_tc_row_points_file`] for the common case of writing directly to a file path.
///
/// # Errors
///
/// An error if a row holds no points, if the ordinals do not strictly increase, if some rows carry
/// column indices and others do not, or if the row pitch is not a positive distance.
pub fn write_tc_row_points_to<W: Write>(
    writer: &mut W,
    scan: &RowPointsScan3,
    tol: f64,
) -> Result<()> {
    tc::write_one_to(writer, &to_container(scan)?, tol)?;
    Ok(())
}

/// Deserialize a scan from a tcrpf3-format byte stream.
///
/// Recovered positions are guaranteed to remain within the tolerance supplied during writing, and
/// the grouping is exact.
///
/// Use [`read_tc_row_points_file`] for the common case of reading from a file path.
pub fn read_tc_row_points_from<R: Read>(reader: &mut R) -> Result<RowPointsScan3> {
    from_container(tc::read_one_from(reader)?)
}

/// Write a scan to a tcrpf3 file at the given path. See [`write_tc_row_points_to`].
pub fn write_tc_row_points_file(path: &Path, scan: &RowPointsScan3, tol: f64) -> Result<()> {
    tc::write_one_file(path, &to_container(scan)?, tol)?;
    Ok(())
}

/// Read a scan from a tcrpf3 file at the given path. See [`read_tc_row_points_from`].
pub fn read_tc_row_points_file(path: &Path) -> Result<RowPointsScan3> {
    from_container(tc::read_one_file(path)?)
}

// ================================================================================================
// Collections
// ================================================================================================

/// Serialize any number of scans into a single tcrpf3 stream, every one at the same storage
/// tolerance.
///
/// Order is preserved, and each scan retains its row pitch, axis, and metadata. Use this function
/// for a series of captures of one scene, such as captures at different exposures.
pub fn write_tc_row_points_all_to<W: Write>(
    writer: &mut W,
    scans: &[RowPointsScan3],
    tol: f64,
) -> Result<()> {
    let items = scans.iter().map(to_container).collect::<Result<Vec<_>>>()?;
    tc::write_to(writer, &items, tol)?;
    Ok(())
}

/// Deserialize any number of scans from a tcrpf3-format byte stream, in the order they were
/// written. See [`write_tc_row_points_all_to`].
pub fn read_tc_row_points_all_from<R: Read>(reader: &mut R) -> Result<Vec<RowPointsScan3>> {
    tc::read_from(reader)?
        .into_iter()
        .map(from_container)
        .collect()
}

/// Write any number of scans to a tcrpf3 file. See [`write_tc_row_points_all_to`].
pub fn write_tc_row_points_all_file(path: &Path, scans: &[RowPointsScan3], tol: f64) -> Result<()> {
    let items = scans.iter().map(to_container).collect::<Result<Vec<_>>>()?;
    tc::write_file(path, &items, tol)?;
    Ok(())
}

/// Read any number of scans from a tcrpf3 file. See [`read_tc_row_points_all_from`].
pub fn read_tc_row_points_all_file(path: &Path) -> Result<Vec<RowPointsScan3>> {
    tc::read_file(path)?
        .into_iter()
        .map(from_container)
        .collect()
}

// ================================================================================================
// Loading
// ================================================================================================

/// An uncertainty model for row-organized points that receives the complete point.
///
/// Unlike `Lptf3UncertaintyModel`, which accepts only x and z, this model receives all coordinates.
/// A laser profiler's uncertainty does not depend on the profile's position along the sweep, while
/// a snapshot sensor's uncertainty varies across its full field of view.
pub trait RowPointsUncertaintyModel {
    /// The predicted standard deviation of the measurement at `p`.
    fn value(&self, p: &Point3) -> f64;
}

/// Read a tcrpf3 file and return a [`PointCloud3`], applying the requested thinning.
///
/// The load modes have the same meaning as for LPTF3:
///   - [`Lptf3Load::All`]: every point in the file.
///   - [`Lptf3Load::TakeEveryN`]: approximately one point per `n` row pitches in both directions,
///     producing approximately square spacing in the along/across plane.
///   - [`Lptf3Load::SmoothSample`]: the same thinning, with each retained point moved to a
///     Gaussian-weighted local mean of its full-resolution neighborhood. Row structure enables the
///     smoothing operation to use points that thinning will discard.
///
/// The format does not store the sensor's color channel, so the returned cloud contains geometry
/// only.
pub fn load_tc_row_points(file_path: &Path, load: Lptf3Load) -> Result<PointCloud3> {
    let scan = loaded_row_scan(file_path, load)?;
    let points = scan.rows().iter().flatten().copied().collect::<Vec<_>>();
    Ok(PointCloud3::new(points))
}

/// Read a tcrpf3 file and mesh it by triangulating between adjacent rows, returning a
/// [`MeshData3`].
///
/// See [`load_tc_row_points`] for what the load modes do, and
/// [`crate::io::row_scan::compute_row_scan_mesh`] for the triangulation and for why some points do
/// not make it into the mesh.
///
/// # Attributes
///
/// If supplied, the uncertainty model is evaluated at every point in the finished mesh, and the
/// results are stored as `point_stdev`. Evaluation occurs after the mesher removes points that
/// belong to no face, keeping the attribute aligned with retained points.
pub fn load_tc_row_points_mesh_data(
    file_path: &Path,
    load: Lptf3Load,
    uncertainty: Option<&dyn RowPointsUncertaintyModel>,
) -> Result<MeshData3> {
    let scan = loaded_row_scan(file_path, load)?;
    let mut mesh = compute_row_scan_mesh(&scan)?;

    if let Some(m) = uncertainty {
        let stdev = mesh.points().iter().map(|p| m.value(p)).collect::<Vec<_>>();
        mesh.set_point_stdev(Some(stdev))?;
    }

    Ok(mesh)
}

/// Read a tcrpf3 file into the shared row intermediate, applying the requested thinning.
fn loaded_row_scan(file_path: &Path, load: Lptf3Load) -> Result<RowScan3> {
    let scan = read_tc_row_points_file(file_path)?.to_row_scan()?;

    match load {
        Lptf3Load::All => Ok(scan),
        Lptf3Load::TakeEveryN(n) => thin_row_scan(&scan, n),
        Lptf3Load::SmoothSample(params) => smooth_row_scan(&scan, params),
    }
}

// ================================================================================================
// Conversion
// ================================================================================================

/// A scan as the storable container, with the keys this module owns merged into its metadata.
fn to_container(scan: &RowPointsScan3) -> Result<TcScan> {
    if !scan.row_pitch.is_finite() || scan.row_pitch <= 0.0 {
        return Err(format!("row pitch {} is not a positive distance", scan.row_pitch).into());
    }

    let mut metadata = scan.metadata.clone();
    metadata.insert(ROW_PITCH.into(), scan.row_pitch.into());
    if let Some(col_pitch) = scan.col_pitch {
        metadata.insert(COL_PITCH.into(), col_pitch.into());
    }
    metadata.insert(ALONG_AXIS.into(), axis_name(scan.along).into());

    let rows = scan
        .rows
        .iter()
        .map(|row| {
            let points = row.points.iter().map(|p| [p.x, p.y, p.z]).collect();
            let out = TcRow::new(row.ordinal, points);
            match &row.columns {
                Some(c) => out.with_columns(c.clone()),
                None => out,
            }
        })
        .collect();

    let mut item = TcScan::new(rows);
    item.name = scan.name.clone();
    item.metadata = metadata;

    // Validate here so a caller that constructs a scan manually receives the error from this
    // conversion function instead of the writer.
    item.validate()?;

    Ok(item)
}

/// The inverse of [`to_container`].
fn from_container(item: TcScan) -> Result<RowPointsScan3> {
    let row_pitch = item
        .metadata
        .get(ROW_PITCH)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| {
            format!("tcrpf3 file has no `{ROW_PITCH}` metadata; was it written by engeom?")
        })?;

    if !row_pitch.is_finite() || row_pitch <= 0.0 {
        return Err(format!("tcrpf3 file records a row pitch of {row_pitch}").into());
    }

    let col_pitch = item.metadata.get(COL_PITCH).and_then(|v| v.as_f64());

    let along = match item.metadata.get(ALONG_AXIS).and_then(|v| v.as_text()) {
        None | Some("x") => RowAxis::X,
        Some("y") => RowAxis::Y,
        Some(other) => {
            return Err(format!(
                "tcrpf3 file records `{ALONG_AXIS}` as {other:?}, which is not \"x\" or \"y\""
            )
            .into());
        }
    };

    // Promote this module's keys to fields and remove them from the caller-visible map so metadata
    // cannot diverge from those fields during a round trip.
    let mut metadata = item.metadata;
    for key in [ROW_PITCH, COL_PITCH, ALONG_AXIS] {
        metadata.remove(key);
    }

    let rows = item
        .rows
        .into_iter()
        .map(|row| PointRow3 {
            ordinal: row.ordinal,
            points: row
                .points
                .iter()
                .map(|p| Point3::new(p[0], p[1], p[2]))
                .collect(),
            columns: row.columns,
        })
        .collect();

    Ok(RowPointsScan3 {
        rows,
        row_pitch,
        col_pitch,
        along,
        name: item.name,
        metadata,
    })
}

fn axis_name(axis: RowAxis) -> &'static str {
    match axis {
        RowAxis::X => "x",
        RowAxis::Y => "y",
    }
}

// ================================================================================================
// Tests
// ================================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::mesh::algorithms::normals::compute_face_normal;
    use approx::assert_relative_eq;
    use std::io::Cursor;

    const TOL: f64 = 1e-6;
    const ROW_PITCH_MM: f64 = 0.1;
    const COL_PITCH_MM: f64 = 0.05;

    /// A synthetic snapshot scan, with the row geometry such a sensor really produces.
    ///
    /// A raster row is a plane through the camera rather than a constant-y line, so along a row the
    /// y-coordinate varies with depth: here by `slope` millimeters per millimeter of z, which is
    /// the same order as the several millimeters observed in the Gocator data. A mesher that took a
    /// row's strip coordinate from one of its points rather than from its ordinal would see the
    /// rows interleaved rather than stacked.
    fn snapshot_scan(rows: u32, cols: u32, slope: f64) -> RowPointsScan3 {
        let out = (0..rows)
            .map(|i| {
                let points = (0..cols)
                    .map(|j| {
                        // A gentle bowl, so the surface is not flat and the depth really varies.
                        let x = j as f64 * COL_PITCH_MM;
                        let nominal_y = i as f64 * ROW_PITCH_MM;
                        let z = 10.0 + 0.5 * (x - 1.0).powi(2);
                        Point3::new(x, nominal_y + slope * z, z)
                    })
                    .collect();
                PointRow3::new(i, points).with_columns((0..cols).collect())
            })
            .collect();

        RowPointsScan3::new(out, ROW_PITCH_MM)
            .with_col_pitch(COL_PITCH_MM)
            .named("synthetic")
    }

    /// The same scan with its rows running along world y instead of world x, which is what a
    /// column-major sensor produces. Built by swapping the roles of x and y in the construction
    /// rather than by transposing a built scan, since a coordinate swap is a reflection.
    fn transposed_scan(rows: u32, cols: u32) -> RowPointsScan3 {
        let out = (0..rows)
            .map(|i| {
                let points = (0..cols)
                    .map(|j| {
                        let y = j as f64 * COL_PITCH_MM;
                        let x = i as f64 * ROW_PITCH_MM;
                        let z = 10.0 + 0.5 * (y - 1.0).powi(2);
                        Point3::new(x, y, z)
                    })
                    .collect();
                PointRow3::new(i, points)
            })
            .collect();

        RowPointsScan3::new(out, ROW_PITCH_MM).with_along_axis(RowAxis::Y)
    }

    fn check_round_trip(original: &RowPointsScan3, back: &RowPointsScan3, tol: f64) {
        assert_eq!(back.rows.len(), original.rows.len());
        assert_eq!(back.name, original.name);
        assert_relative_eq!(back.row_pitch, original.row_pitch);
        assert_eq!(back.col_pitch, original.col_pitch);
        assert_eq!(back.along, original.along);
        assert_eq!(back.metadata, original.metadata);

        for (a, b) in original.rows.iter().zip(back.rows.iter()) {
            assert_eq!(b.ordinal, a.ordinal);
            assert_eq!(b.columns, a.columns);
            assert_eq!(b.points.len(), a.points.len());
            for (p, q) in a.points.iter().zip(b.points.iter()) {
                assert_relative_eq!(p, q, epsilon = tol);
            }
        }
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "engeom-tcrpf3-{}-{}.tcrpf3",
            name,
            std::process::id()
        ))
    }

    #[test]
    fn round_trips_through_bytes() {
        let scan = snapshot_scan(12, 30, 0.3);

        let mut buf = Vec::new();
        write_tc_row_points_to(&mut buf, &scan, TOL).unwrap();
        let back = read_tc_row_points_from(&mut Cursor::new(&buf)).unwrap();

        check_round_trip(&scan, &back, TOL);
    }

    #[test]
    fn round_trips_through_a_file() {
        let scan = snapshot_scan(8, 20, 0.2);
        let path = temp_path("file");

        write_tc_row_points_file(&path, &scan, TOL).unwrap();
        let back = read_tc_row_points_file(&path).unwrap();

        check_round_trip(&scan, &back, TOL);
        let _ = std::fs::remove_file(&path);
    }

    /// The keys this module owns are promoted to fields on read, so a caller's own metadata comes
    /// back unchanged, without engeom's bookkeeping mixed into it.
    #[test]
    fn caller_metadata_survives_and_stays_separate() {
        let scan = snapshot_scan(4, 10, 0.1)
            .with_meta(SOURCE_FRAME, "gocator-3x00 left-handed, x negated")
            .with_meta("engeom.exposure_us", 500i64);

        let mut buf = Vec::new();
        write_tc_row_points_to(&mut buf, &scan, TOL).unwrap();
        let back = read_tc_row_points_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.metadata.len(), 2);
        assert!(!back.metadata.contains_key(ROW_PITCH));
        assert_eq!(
            back.metadata.get(SOURCE_FRAME).and_then(|v| v.as_text()),
            Some("gocator-3x00 left-handed, x negated")
        );
        check_round_trip(&scan, &back, TOL);
    }

    /// Ordinals are the row's place in the sweep, so a gap has to survive: it is what tells the
    /// mesher not to join the two rows across it.
    #[test]
    fn ordinal_gaps_survive() {
        let mut scan = snapshot_scan(5, 8, 0.0);
        scan.rows[3].ordinal = 40;
        scan.rows[4].ordinal = 41;

        let mut buf = Vec::new();
        write_tc_row_points_to(&mut buf, &scan, TOL).unwrap();
        let back = read_tc_row_points_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(
            back.rows.iter().map(|r| r.ordinal).collect::<Vec<_>>(),
            vec![0, 1, 2, 40, 41]
        );
    }

    /// A tcrpf3 written by something other than engeom has no row pitch in it, and guessing one
    /// would silently change where the mesher thinks every row is.
    #[test]
    fn a_file_without_a_row_pitch_is_rejected() {
        let item = TcScan::new(vec![TcRow::new(0, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])]);
        let mut buf = Vec::new();
        tc::write_one_to(&mut buf, &item, TOL).unwrap();

        let err = read_tc_row_points_from(&mut Cursor::new(&buf)).unwrap_err();
        assert!(err.to_string().contains(ROW_PITCH), "unhelpful: {err}");
    }

    #[test]
    fn an_unknown_along_axis_is_rejected() {
        let item = TcScan::new(vec![TcRow::new(0, vec![[0.0, 0.0, 0.0]])])
            .with_meta(ROW_PITCH, 0.1)
            .with_meta(ALONG_AXIS, "z");
        let mut buf = Vec::new();
        tc::write_one_to(&mut buf, &item, TOL).unwrap();

        assert!(read_tc_row_points_from(&mut Cursor::new(&buf)).is_err());
    }

    #[test]
    fn a_non_positive_row_pitch_is_refused_on_write() {
        let mut scan = snapshot_scan(3, 5, 0.0);
        scan.row_pitch = 0.0;

        assert!(write_tc_row_points_to(&mut Vec::new(), &scan, TOL).is_err());
    }

    #[test]
    fn ordinals_that_do_not_increase_are_refused_on_write() {
        let mut scan = snapshot_scan(3, 5, 0.0);
        scan.rows[2].ordinal = 1;

        assert!(write_tc_row_points_to(&mut Vec::new(), &scan, TOL).is_err());
    }

    // ============================================================================================
    // Meshing
    // ============================================================================================

    /// The point of the format. A snapshot sensor's rows are planes through the camera, so their y
    /// varies by millimeters with depth while the rows themselves are one-tenth of a millimeter
    /// apart. The strip coordinate has to come from the ordinal, and this is the test that fails if
    /// it ever goes back to coming from a point.
    #[test]
    fn meshes_a_scan_whose_rows_are_not_constant_y() {
        let scan = snapshot_scan(24, 60, 0.3);
        let path = temp_path("mesh");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let mesh = load_tc_row_points_mesh_data(&path, Lptf3Load::All, None).unwrap();

        // 23 strips of 59 quads, two triangles each, and nothing rejected on a smooth surface.
        assert_eq!(mesh.face_count(), 23 * 59 * 2);
        assert_eq!(mesh.point_count(), 24 * 60);

        for f in mesh.faces() {
            let p = f.map(|i| mesh.points()[i as usize]);
            let n = compute_face_normal(&p).expect("a face is degenerate");
            assert!(n.z > 0.0, "a face points away from the sensor: {n:?}");
        }

        let _ = std::fs::remove_file(&path);
    }

    /// A column-major sensor's strips run along y. Flattening to (y, x) is a mirror, so the winding
    /// has to be put back or the whole mesh comes out inside out.
    #[test]
    fn a_scan_along_y_meshes_with_its_normals_the_same_way_up() {
        let scan = transposed_scan(20, 50);
        let path = temp_path("along-y");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let mesh = load_tc_row_points_mesh_data(&path, Lptf3Load::All, None).unwrap();
        assert_eq!(mesh.face_count(), 19 * 49 * 2);

        for f in mesh.faces() {
            let p = f.map(|i| mesh.points()[i as usize]);
            let n = compute_face_normal(&p).expect("a face is degenerate");
            assert!(n.z > 0.0, "a face points away from the sensor: {n:?}");
        }

        let _ = std::fs::remove_file(&path);
    }

    /// A gap in the ordinals is a gap in the scan, and the mesher must not stitch across it.
    #[test]
    fn rows_separated_by_an_ordinal_gap_are_not_joined() {
        let mut scan = snapshot_scan(6, 20, 0.0);
        for row in scan.rows.iter_mut().skip(3) {
            row.ordinal += 50;
        }

        let path = temp_path("gap");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();
        let mesh = load_tc_row_points_mesh_data(&path, Lptf3Load::All, None).unwrap();

        // Two groups of three rows: two strips each, rather than the five a contiguous scan gives.
        assert_eq!(mesh.face_count(), 4 * 19 * 2);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn loads_every_point_as_a_cloud() {
        let scan = snapshot_scan(10, 25, 0.25);
        let path = temp_path("cloud");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let cloud = load_tc_row_points(&path, Lptf3Load::All).unwrap();
        assert_eq!(cloud.points().len(), 250);

        let _ = std::fs::remove_file(&path);
    }

    /// Thinning is meant to be roughly square in the along/across plane rather than one-directional,
    /// so taking every fourth row should also drop points along the rows.
    #[test]
    fn thinning_drops_points_in_both_directions() {
        let scan = snapshot_scan(40, 100, 0.2);
        let path = temp_path("thin");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let full = load_tc_row_points(&path, Lptf3Load::All).unwrap();
        let thin = load_tc_row_points(&path, Lptf3Load::TakeEveryN(4)).unwrap();

        // One row in four, and one point in every four row-pitches along a row whose points are
        // half a row pitch apart: about 1/4 x 1/8 of the total.
        assert!(thin.points().len() < full.points().len() / 20);
        assert!(!thin.points().is_empty());

        let _ = std::fs::remove_file(&path);
    }

    /// The smoothing filter is what the row structure buys that a cloud cannot: it uses the points
    /// that are about to be discarded. It moves points along z only.
    #[test]
    fn smoothing_moves_points_along_z_only() {
        let scan = snapshot_scan(40, 100, 0.2);
        let path = temp_path("smooth");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let params = crate::io::Lptf3DsParams::new(4, 1.5, 1.0, 1.0);
        let thin = load_tc_row_points(&path, Lptf3Load::TakeEveryN(4)).unwrap();
        let smooth = load_tc_row_points(&path, Lptf3Load::SmoothSample(params)).unwrap();

        assert_eq!(smooth.points().len(), thin.points().len());
        for (a, b) in thin.points().iter().zip(smooth.points().iter()) {
            assert_relative_eq!(a.x, b.x, epsilon = 1e-12);
            assert_relative_eq!(a.y, b.y, epsilon = 1e-12);
        }

        let _ = std::fs::remove_file(&path);
    }

    /// An uncertainty model for these points sees the whole point, not just x and z, because a
    /// snapshot sensor's uncertainty varies across its whole field of view.
    #[test]
    fn the_uncertainty_model_sees_the_whole_point() {
        struct Ramp;
        impl RowPointsUncertaintyModel for Ramp {
            fn value(&self, p: &Point3) -> f64 {
                0.001 + 0.0001 * p.y
            }
        }

        let scan = snapshot_scan(12, 30, 0.2);
        let path = temp_path("stdev");
        write_tc_row_points_file(&path, &scan, TOL).unwrap();

        let mesh = load_tc_row_points_mesh_data(&path, Lptf3Load::All, Some(&Ramp)).unwrap();
        let stdev = mesh.point_stdev().expect("a model was supplied");

        for (p, s) in mesh.points().iter().zip(stdev.iter()) {
            assert_relative_eq!(*s, 0.001 + 0.0001 * p.y, epsilon = 1e-12);
        }

        let _ = std::fs::remove_file(&path);
    }

    // ============================================================================================
    // Collections
    // ============================================================================================

    /// A capture series is several scans of one scene, so a collection has to keep its order and
    /// let its members differ from each other.
    #[test]
    fn a_collection_round_trips_with_mixed_members() {
        let scans = vec![
            snapshot_scan(6, 12, 0.1).named("500"),
            snapshot_scan(8, 10, 0.3).named("1000"),
            transposed_scan(5, 9),
        ];

        let mut buf = Vec::new();
        write_tc_row_points_all_to(&mut buf, &scans, TOL).unwrap();
        let back = read_tc_row_points_all_from(&mut Cursor::new(&buf)).unwrap();

        assert_eq!(back.len(), 3);
        for (a, b) in scans.iter().zip(back.iter()) {
            check_round_trip(a, b, TOL);
        }
        assert_eq!(back[2].along, RowAxis::Y);
    }

    /// The two write paths produce the same format, so the collection reader has to accept a file
    /// from either. Only the single reader is picky, and only about count.
    #[test]
    fn the_single_reader_refuses_a_collection() {
        let scans = vec![snapshot_scan(4, 8, 0.0), snapshot_scan(4, 8, 0.0)];
        let mut buf = Vec::new();
        write_tc_row_points_all_to(&mut buf, &scans, TOL).unwrap();

        assert!(read_tc_row_points_from(&mut Cursor::new(&buf)).is_err());
        assert_eq!(
            read_tc_row_points_all_from(&mut Cursor::new(&buf))
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn an_empty_collection_is_a_valid_file() {
        let mut buf = Vec::new();
        write_tc_row_points_all_to(&mut buf, &[], TOL).unwrap();

        assert!(
            read_tc_row_points_all_from(&mut Cursor::new(&buf))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_tcmesh_is_not_read_as_row_points() {
        let bad = b"NOPE0000\x00\x00";
        assert!(read_tc_row_points_from(&mut Cursor::new(bad)).is_err());
    }
}
