//! The row-structured intermediate that sits between a sensor file and a triangle mesh.
//!
//! Sweeping sensors produce points grouped into rows, such as a laser profiler's frames, a snapshot
//! profiler's raster rows, or a time-of-flight camera's scanlines. This grouping makes the points
//! meshable without general surface reconstruction because two neighboring rows can be triangulated
//! into a strip. [`RowScan3`] represents the grouping, and [`compute_row_scan_mesh`] provides the
//! strip mesher shared by every format that can create one.
//!
//! # A row is not a constant-coordinate line
//!
//! A laser-profiler row lies at one y-coordinate. A snapshot-sensor raster row represents a plane
//! through the camera, so its y-coordinate varies by millimeters with depth across the row. The
//! strip coordinate used for flattening is therefore the **nominal** value in [`RowScan3::row_y`],
//! usually the row's sweep ordinal multiplied by the row pitch. It is never derived from a point in
//! the row.
//!
//! # Along and across
//!
//! By default, rows run along world x and the sweep advances along world y, as with a laser
//! profiler. A sensor whose strips run in the other direction uses [`RowAxis::Y`]. Transposing
//! coordinates during conversion would reflect the points out of the sensor's frame. The mesher
//! flattens coordinates to `(along, across)` and reverses face winding for the mirrored case, so
//! surface normals remain on the same side for either axis.

use crate::common::triangulation::parallel_row2::{StripRowPoint, build_parallel_row_strip};
use crate::geom3::mesh::MeshData3;
use crate::io::lptf3::{Lptf3DsParams, adjust_by_gwm, gaussian_weight};
use crate::{Point3, Result};
use rayon::prelude::*;

/// The maximum edge ratio for candidate faces evaluated in the flattened 2D space of two joined
/// rows.
pub(crate) const STRIP_EDGE_RATIO: f64 = 2.0;

/// The maximum edge ratio measured on the 3D points. This rejects faces that appear valid in the
/// flattened space but span a scan-depth discontinuity.
pub(crate) const WORLD_EDGE_RATIO: f64 = 5.0;

/// Which world coordinate varies along a row.
///
/// The other coordinate is the sweep direction measured by [`RowScan3::row_y`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowAxis {
    /// Rows run along world x, and the sweep advances along world y. This is the usual convention
    /// produced by a laser profile sensor.
    X,
    /// Rows run along world y, and the sweep advances along world x. This represents a sensor whose
    /// raster is transposed relative to the usual convention.
    Y,
}

impl RowAxis {
    /// The coordinate that varies along a row.
    pub fn along(self, p: &Point3) -> f64 {
        match self {
            RowAxis::X => p.x,
            RowAxis::Y => p.y,
        }
    }

    /// The coordinate the sweep advances in.
    pub fn across(self, p: &Point3) -> f64 {
        match self {
            RowAxis::X => p.y,
            RowAxis::Y => p.x,
        }
    }
}

/// Points grouped into the rows a sweeping sensor produced, in sweep order.
///
/// Each row is sorted in ascending order by the along coordinate, and `row_y` gives its nominal
/// strip coordinate. When `colors` is `Some`, it has the same row and column structure as the
/// points, so flattened colors and points remain aligned element by element.
pub struct RowScan3 {
    rows: Vec<Vec<Point3>>,
    row_y: Vec<f64>,
    row_spacing: f64,
    take_every: u32,
    along: RowAxis,
    colors: Option<Vec<Vec<u8>>>,
}

impl RowScan3 {
    /// Create a scan and sort its rows by the along coordinate.
    ///
    /// `row_spacing` is the nominal distance between consecutive sensor rows before thinning.
    /// `take_every` is the thinning factor already applied; together, they define the actual row
    /// separation.
    ///
    /// # Errors
    ///
    /// An error if `row_y` has a different length than `rows`, or if `row_spacing` is not finite
    /// and positive.
    pub fn new(
        rows: Vec<Vec<Point3>>,
        row_y: Vec<f64>,
        row_spacing: f64,
        take_every: u32,
        along: RowAxis,
    ) -> Result<Self> {
        let mut scan = Self::new_unchecked(rows, row_y, row_spacing, take_every, along)?;
        for row in scan.rows.iter_mut() {
            row.sort_by(|a, b| along.along(a).total_cmp(&along.along(b)));
        }
        Ok(scan)
    }

    /// Create a scan without sorting rows that are already ordered by the along coordinate, as most
    /// sensors provide them.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`RowScan3::new`]. The caller must ensure correct sort order;
    /// downstream operations do not check it.
    pub fn new_unchecked(
        rows: Vec<Vec<Point3>>,
        row_y: Vec<f64>,
        row_spacing: f64,
        take_every: u32,
        along: RowAxis,
    ) -> Result<Self> {
        if rows.len() != row_y.len() {
            return Err(format!(
                "a row scan needs one strip coordinate per row, but has {} rows and {} coordinates",
                rows.len(),
                row_y.len()
            )
            .into());
        }

        if !row_spacing.is_finite() || row_spacing <= 0.0 {
            return Err(format!("row spacing {row_spacing} is not a positive distance").into());
        }

        Ok(Self {
            rows,
            row_y,
            row_spacing,
            take_every: take_every.max(1),
            along,
            colors: None,
        })
    }

    /// The same scan carrying the sensor's single-channel color/intensity value per point.
    ///
    /// # Errors
    ///
    /// An error if the color rows do not match the point rows element for element.
    pub fn with_colors(mut self, colors: Option<Vec<Vec<u8>>>) -> Result<Self> {
        if let Some(c) = &colors
            && (c.len() != self.rows.len()
                || c.iter()
                    .zip(self.rows.iter())
                    .any(|(a, b)| a.len() != b.len()))
        {
            return Err("the color rows do not match the point rows".into());
        }
        self.colors = colors;
        Ok(self)
    }

    /// The rows of points, in sweep order, each sorted in the along coordinate.
    pub fn rows(&self) -> &[Vec<Point3>] {
        &self.rows
    }

    /// The nominal strip coordinate of each row.
    pub fn row_y(&self) -> &[f64] {
        &self.row_y
    }

    /// The sensor's nominal distance between consecutive rows, before thinning.
    pub fn row_spacing(&self) -> f64 {
        self.row_spacing
    }

    /// The thinning factor already applied to the rows.
    pub fn take_every(&self) -> u32 {
        self.take_every
    }

    /// Which world coordinate varies along a row.
    pub fn along_axis(&self) -> RowAxis {
        self.along
    }

    /// The sensor's color/intensity channel, if it has one.
    pub fn colors(&self) -> Option<&Vec<Vec<u8>>> {
        self.colors.as_ref()
    }

    /// How many points the rows hold in total.
    pub fn point_count(&self) -> usize {
        self.rows.iter().map(|r| r.len()).sum()
    }

    /// The largest strip coordinate gap two rows may span and still be triangulated together.
    ///
    /// This is twice the actual row spacing, allowing meshing across one dropped row while
    /// preserving larger scan breaks.
    pub fn max_strip_spacing(&self) -> f64 {
        self.take_every as f64 * self.row_spacing * 2.0
    }
}

/// Build a triangle mesh by triangulating between each adjacent pair of rows.
///
/// Rows that are too far apart in the strip coordinate are not joined. Faces spanning a depth
/// discontinuity are rejected, and points that belong to no face are removed. When present, the
/// sensor's color channel is transferred to the mesh. The caller must attach other point data after
/// meshing because the mesh point buffer is a subset of the scan point buffer.
///
/// # Orphan points
///
/// **Points that belong to no face are discarded.** The edge criteria reject triangles around a
/// dropout or depth discontinuity, and the first and last scan rows can have nothing to join. A real
/// scan therefore leaves some measured points unconnected. These points are measurements but do not
/// form part of the requested surface mesh.
pub fn compute_row_scan_mesh(scan: &RowScan3) -> Result<MeshData3> {
    let (mut points, strip_rows) = flatten_rows(scan);
    let mut colors = flatten_colors(scan);
    let mut faces = build_faces(scan, &points, &strip_rows)?;

    drop_orphan_points(&mut points, &mut faces, &mut colors);

    let mut mesh = MeshData3::new(points, faces)?;
    if colors.is_some() {
        mesh.set_point_colors(colors)?;
    }

    Ok(mesh)
}

/// Copy row-structured points into one buffer and return per-row strip points containing each
/// point's buffer index.
fn flatten_rows(scan: &RowScan3) -> (Vec<Point3>, Vec<Vec<StripRowPoint>>) {
    let mut points = Vec::with_capacity(scan.point_count());
    let mut strip_rows = Vec::with_capacity(scan.rows.len());

    for row in scan.rows.iter() {
        let mut strip_row = Vec::with_capacity(row.len());
        for p in row.iter() {
            strip_row.push(StripRowPoint::new(scan.along.along(p), points.len() as u32));
            points.push(*p);
        }
        strip_rows.push(strip_row);
    }

    (points, strip_rows)
}

/// Flatten the sensor's color or intensity channel into one buffer aligned with the point buffer,
/// expanding each 8-bit value to gray.
fn flatten_colors(scan: &RowScan3) -> Option<Vec<[u8; 3]>> {
    scan.colors
        .as_ref()
        .map(|color_rows| color_rows.iter().flatten().map(|&c| [c, c, c]).collect())
}

/// Triangulate each adjacent row pair and discard faces spanning gaps larger than the allowed row
/// spacing.
fn build_faces(
    scan: &RowScan3,
    points: &[Point3],
    strip_rows: &[Vec<StripRowPoint>],
) -> Result<Vec<[u32; 3]>> {
    let max_spacing = scan.max_strip_spacing();
    let mut faces = Vec::new();

    // A scan with zero or one row has nothing to triangulate between, and the subtraction below
    // must not be allowed to underflow.
    for row_i in 0..strip_rows.len().saturating_sub(1) {
        if strip_rows[row_i].is_empty() || strip_rows[row_i + 1].is_empty() {
            continue; // Skip empty rows
        }

        // Use each row's nominal strip coordinate. A snapshot-sensor row is a plane through the
        // camera rather than a constant-y line, so its points differ in y by millimeters.
        let y0 = scan.row_y[row_i];
        let y1 = scan.row_y[row_i + 1];

        // Skip triangulation when the rows are too far apart.
        if (y1 - y0).abs() > max_spacing {
            continue;
        }

        let row0 = &strip_rows[row_i];
        let row1 = &strip_rows[row_i + 1];

        // Build the strip triangulation between the two rows.
        let r = build_parallel_row_strip(row0, y0, row1, y1, STRIP_EDGE_RATIO)?;
        for [i0, i1, i2] in r {
            // Check the edge ratio on the actual 3D points, which the flattened triangulation
            // could not see.
            let pa = points[i0 as usize];
            let pb = points[i1 as usize];
            let pc = points[i2 as usize];

            let ea = (pa - pb).norm();
            let eb = (pb - pc).norm();
            let ec = (pc - pa).norm();

            // TODO: decide whether to reject faces which stand exactly on end.
            //
            // The sensor quantizes x, so two points in the same row can land on the same x value
            // at different depths. A face joining that pair to a point in the next row is exactly
            // vertical, with a normal that has no z component at all, and it represents a step in
            // the surface which the sensor could not actually resolve rather than measured
            // geometry. On the sample scan these are about 0.05% of the faces.
            //
            // Rejecting these faces would leave points on opposite sides of the step unconnected.
            // This might represent the unresolved surface more accurately, or it might create holes
            // along every part edge. Resolve this choice before relying on the loader for measurement.
            let edge_ratio = ea.max(eb).max(ec) / max_spacing;
            if edge_ratio < WORLD_EDGE_RATIO {
                // Flattening to (y, x) reflects the usual (x, y) coordinates. Reverse the strip
                // winding to restore the original orientation.
                match scan.along {
                    RowAxis::X => faces.push([i1, i0, i2]),
                    RowAxis::Y => faces.push([i0, i1, i2]),
                }
            }
        }
    }

    Ok(faces)
}

/// Remove every point that no face references, renumbering faces and subsetting the color
/// buffer to match. Returns the number of points removed.
///
/// A scan always produces some orphan points because edge criteria reject faces around dropouts and
/// depth discontinuities, and the first and last rows can have nothing to join. These measured
/// points do not form part of the requested surface. Retaining them would produce undefined point
/// normals and bounds that include geometry absent from the triangulation. Load the scan as a point
/// cloud when every measured point matters.
fn drop_orphan_points(
    points: &mut Vec<Point3>,
    faces: &mut [[u32; 3]],
    colors: &mut Option<Vec<[u8; 3]>>,
) -> usize {
    let mut used = vec![false; points.len()];
    for f in faces.iter() {
        for i in f.iter() {
            used[*i as usize] = true;
        }
    }

    let removed = used.iter().filter(|u| !**u).count();
    if removed == 0 {
        return 0;
    }

    // Build the map from old index to new, then compact everything indexed by a point.
    let mut remap = vec![u32::MAX; points.len()];
    let mut next = 0u32;
    for (old, u) in used.iter().enumerate() {
        if *u {
            remap[old] = next;
            next += 1;
        }
    }

    *points = points
        .iter()
        .zip(used.iter())
        .filter(|(_, u)| **u)
        .map(|(p, _)| *p)
        .collect();

    if let Some(c) = colors {
        *c = c
            .iter()
            .zip(used.iter())
            .filter(|(_, u)| **u)
            .map(|(v, _)| *v)
            .collect();
    }

    for f in faces.iter_mut() {
        for i in f.iter_mut() {
            *i = remap[*i as usize];
        }
    }

    removed
}

/// Which points of which rows survive a thinning by `take_every`.
///
/// Rows are kept when they are at least `take_every` row spacings past the last one kept, and
/// within a kept row points are kept at roughly that same spacing along the row, so the result is
/// approximately square in the along/across plane rather than thinned in one direction only. This
/// is the same rule the LPTF3 loader applies, expressed over rows already in memory.
fn thinning_plan(scan: &RowScan3, take_every: u32) -> Vec<(usize, Vec<usize>)> {
    let spacing = take_every as f64 * scan.row_spacing;
    let mut plan = Vec::new();
    let mut last_y: Option<f64> = None;

    for (row_i, row) in scan.rows.iter().enumerate() {
        if row.is_empty() {
            continue;
        }
        if let Some(y) = last_y
            && (scan.row_y[row_i] - y).abs() < spacing
        {
            continue;
        }
        last_y = Some(scan.row_y[row_i]);

        let mut take = Vec::new();
        let mut last_along: Option<f64> = None;
        for (col_i, p) in row.iter().enumerate() {
            let along = scan.along.along(p);
            if let Some(a) = last_along
                && along - a < spacing
            {
                continue;
            }
            last_along = Some(along);
            take.push(col_i);
        }

        plan.push((row_i, take));
    }

    plan
}

/// A thinned copy of the scan, keeping roughly one point per `take_every` row spacings in both
/// directions.
pub fn thin_row_scan(scan: &RowScan3, take_every: u32) -> Result<RowScan3> {
    let plan = thinning_plan(scan, take_every);
    let rows = plan
        .iter()
        .map(|(row_i, take)| take.iter().map(|&col_i| scan.rows[*row_i][col_i]).collect())
        .collect::<Vec<_>>();

    scan_from_plan(scan, &plan, take_every, rows)
}

/// A thinned copy of the scan whose surviving points have been pulled onto a gaussian-weighted
/// local mean of their full-resolution neighbourhood.
///
/// The smoothing is what a bare thinning cannot do: it uses the points that are about to be thrown
/// away, which no longer exist once the scan is a cloud. Each surviving point moves along z only,
/// by at most `params.max_move`.
///
/// # Errors
///
/// An error if `params.take_every` is below 2, since smoothing without thinning is not what this
/// filter does.
pub fn smooth_row_scan(scan: &RowScan3, params: Lptf3DsParams) -> Result<RowScan3> {
    if params.take_every < 2 {
        return Err("take_every must be at least 2".into());
    }

    let plan = thinning_plan(scan, params.take_every);

    // Measure the neighborhood in rows of the full-resolution scan, as the LPTF3
    // downsampling filter does, so the two produce the same result on the same geometry.
    let look_rows = (params.take_every as f64 * params.look_scale.abs()).ceil() as i32;
    let look_dist = look_rows as f64 * scan.row_spacing * 1.25;
    let weight_sigma = params.weight_scale * look_dist;

    let smoothed = plan
        .par_iter()
        .map(|(row_i, take)| {
            let row = &scan.rows[*row_i];
            take.iter()
                .map(|&col_i| {
                    let p = row[col_i];
                    let mut samples = Vec::new();

                    for check_i in (*row_i as i32 - look_rows)..=(*row_i as i32 + look_rows) {
                        if check_i < 0 || check_i >= scan.rows.len() as i32 {
                            continue;
                        }
                        let check_row = &scan.rows[check_i as usize];

                        // The rows are sorted along the axis, so the window can be found rather
                        // than scanned for.
                        let target = scan.along.along(&p) - look_dist;
                        let start = check_row
                            .binary_search_by(|a| scan.along.along(a).total_cmp(&target))
                            .unwrap_or_else(|i| i);

                        for check_p in check_row.iter().skip(start) {
                            let d = (check_p - p).norm();
                            if d <= look_dist {
                                samples.push((*check_p, gaussian_weight(d, weight_sigma)));
                            }
                            if scan.along.along(check_p) > scan.along.along(&p) + look_dist {
                                break;
                            }
                        }
                    }

                    adjust_by_gwm(&p, &samples, params.max_move)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    scan_from_plan(scan, &plan, params.take_every, smoothed)
}

/// The scan a thinning plan describes, given the rows the caller built from that same plan.
fn scan_from_plan(
    scan: &RowScan3,
    plan: &[(usize, Vec<usize>)],
    take_every: u32,
    rows: Vec<Vec<Point3>>,
) -> Result<RowScan3> {
    let row_y = plan.iter().map(|(i, _)| scan.row_y[*i]).collect::<Vec<_>>();

    let colors = scan.colors.as_ref().map(|c| {
        plan.iter()
            .map(|(row_i, take)| take.iter().map(|&col_i| c[*row_i][col_i]).collect())
            .collect::<Vec<Vec<u8>>>()
    });

    // The rows are filled by the caller from the same plan, so they are already in order.
    RowScan3::new_unchecked(rows, row_y, scan.row_spacing, take_every, scan.along)?
        .with_colors(colors)
}
