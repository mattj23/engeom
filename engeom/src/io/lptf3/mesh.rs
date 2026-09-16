//! Construction of a triangle mesh from the row-structured point data in a LPTF3 file.
//!
//! The mesh is built by triangulating between adjacent scan rows rather than by any general
//! surface reconstruction technique, which is only possible because the file format preserves the
//! sensor's original row structure. This is both much faster and much more predictable than a
//! general method operating on an unordered point cloud.
//!
//! This module contains only loading logic. Triangulation is provided by
//! [`compute_row_scan_mesh`](crate::io::row_scan::compute_row_scan_mesh), shared with every other
//! format that can produce a [`RowScan3`](crate::io::row_scan::RowScan3).

use crate::Result;
use crate::geom3::mesh::MeshData3;
use crate::io::Lptf3Load;
use crate::io::lptf3::{Lptf3UncertaintyModel, get_downfilter_point_rows, get_loader_point_rows};
use crate::io::row_scan::compute_row_scan_mesh;
use std::path::Path;

pub fn load_lptf3_mesh_data_core(
    file_path: &Path,
    load: Lptf3Load,
    model: Option<&dyn Lptf3UncertaintyModel>,
) -> Result<MeshData3> {
    let rows = match load {
        Lptf3Load::All => get_loader_point_rows(file_path, None),
        Lptf3Load::TakeEveryN(n) => get_loader_point_rows(file_path, Some(n)),
        Lptf3Load::SmoothSample(params) => get_downfilter_point_rows(file_path, params),
    }?;

    let mut mesh = compute_row_scan_mesh(&rows)?;
    attach_uncertainty(&mut mesh, model)?;

    Ok(mesh)
}

/// Attach modeled measurement uncertainty when a model is supplied.
///
/// Evaluate the model after meshing because the mesher removes points that belong to no face. This
/// keeps the attribute aligned with retained points.
fn attach_uncertainty(
    mesh: &mut MeshData3,
    model: Option<&dyn Lptf3UncertaintyModel>,
) -> Result<()> {
    let Some(m) = model else {
        return Ok(());
    };

    let stdev = mesh
        .points()
        .iter()
        .map(|p| m.value(p.x, p.z))
        .collect::<Vec<_>>();
    mesh.set_point_stdev(Some(stdev))
}

// ===============================================================================================
// Tests
// ===============================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::load_lptf3_mesh_data;
    use approx::assert_relative_eq;
    use std::path::PathBuf;

    /// All of the synthetic files below use a 1 mm resolution in every axis, so a raw coordinate
    /// value maps directly onto a millimetre and the expected geometry can be read off the input.
    const RES_NM: u32 = 1_000_000;

    /// A row of a synthetic scan, given as the raw (x, z, color) values of its points in ascending
    /// order by x, together with the index of the frame it belongs to. The frame index sets the y
    /// position of the row, so a gap in the indices is a gap in the scan.
    struct Row {
        frame: u32,
        points: Vec<(i16, i16, u8)>,
    }

    impl Row {
        /// A flat row of `n` points, one per millimetre, lying at z = 0.
        fn flat(frame: u32, n: i16) -> Self {
            Self {
                frame,
                points: (0..n).map(|x| (x, 0, x as u8)).collect(),
            }
        }
    }

    /// Encode a synthetic LPTF3 file, in the 16-bit coordinate format described at the top of the
    /// `lptf3` module. Written by hand rather than by a writer under test, so that a change to the
    /// reader cannot quietly agree with a matching change to the writer.
    fn encode(y_translation_nm: u32, has_color: bool, rows: &[Row]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(b"LPTF3");
        out.extend_from_slice(&1u16.to_le_bytes()); // version
        out.extend_from_slice(&(if has_color { 0x0002u16 } else { 0u16 }).to_le_bytes()); // flags
        out.push(0); // motion type: fixed y translation
        out.extend_from_slice(&y_translation_nm.to_le_bytes());

        for row in rows {
            out.extend_from_slice(&row.frame.to_le_bytes());
            out.extend_from_slice(&(row.points.len() as u32).to_le_bytes());
            out.extend_from_slice(&0i32.to_le_bytes()); // x offset
            out.extend_from_slice(&0i32.to_le_bytes()); // z offset
            out.extend_from_slice(&RES_NM.to_le_bytes()); // x resolution
            out.extend_from_slice(&RES_NM.to_le_bytes()); // z resolution

            for (x, z, c) in row.points.iter() {
                out.extend_from_slice(&x.to_le_bytes());
                out.extend_from_slice(&z.to_le_bytes());
                if has_color {
                    out.push(*c);
                }
            }
        }

        out
    }

    /// A synthetic file on disk, removed when the test finishes with it. The loader opens paths
    /// rather than readers, so these tests cannot stay in memory the way the PLY ones do.
    struct TempFile {
        path: PathBuf,
    }

    impl TempFile {
        fn new(name: &str, bytes: &[u8]) -> Self {
            let path = std::env::temp_dir().join(format!(
                "engeom-test-{}-{}.lptf3",
                name,
                std::process::id()
            ));
            std::fs::write(&path, bytes).expect("failed to write the synthetic scan");
            Self { path }
        }

        fn path(&self) -> &Path {
            &self.path
        }
    }

    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.path);
        }
    }

    /// A model whose predicted deviation grows linearly with x, so a per-point value can be
    /// checked against the point it was computed from.
    struct RampModel;

    impl Lptf3UncertaintyModel for RampModel {
        fn value(&self, x: f64, _z: f64) -> f64 {
            0.001 + 0.0001 * x
        }
    }

    /// A model which returns a value that is not a standard deviation, to prove the container
    /// rejects it rather than storing it.
    struct NegativeModel;

    impl Lptf3UncertaintyModel for NegativeModel {
        fn value(&self, _x: f64, _z: f64) -> f64 {
            -1.0
        }
    }

    fn flat_grid(rows: u32, cols: i16) -> Vec<Row> {
        (0..rows).map(|f| Row::flat(f, cols)).collect()
    }

    /// A flat scan of consecutive rows should mesh into a fully connected grid, with the points
    /// laid out in scan order.
    #[test]
    fn builds_a_grid_mesh_from_a_flat_scan() -> Result<()> {
        let file = TempFile::new("grid", &encode(RES_NM, false, &flat_grid(4, 5)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert_eq!(mesh.point_count(), 20);

        // Points arrive in scan order: row by row, and within a row in ascending x.
        for row in 0..4 {
            for col in 0..5 {
                let p = mesh.points()[row * 5 + col];
                assert_relative_eq!(p.x, col as f64, epsilon = 1.0e-12);
                assert_relative_eq!(p.y, row as f64, epsilon = 1.0e-12);
                assert_relative_eq!(p.z, 0.0, epsilon = 1.0e-12);
            }
        }

        // Three strips of four quads, each quad split into two triangles.
        assert_eq!(mesh.face_count(), 24);

        // Nothing was supplied to attach, so nothing should have been invented.
        assert!(!mesh.has_attrs());

        Ok(())
    }

    /// The color channel is a per-point value the scan carries, and it used to be dropped on the
    /// floor by the mesh pathway even though the point cloud pathway kept it.
    #[test]
    fn carries_the_color_channel_onto_the_mesh() -> Result<()> {
        let file = TempFile::new("color", &encode(RES_NM, true, &flat_grid(3, 4)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        let colors = mesh
            .point_colors()
            .expect("the file declares a color channel");
        assert_eq!(colors.len(), mesh.point_count());

        // The channel is a single 8-bit value expanded to gray, and `Row::flat` sets it to the
        // point's x index, so it should track the x coordinate of the point it belongs to.
        for (p, c) in mesh.points().iter().zip(colors.iter()) {
            let expected = p.x as u8;
            assert_eq!(*c, [expected, expected, expected]);
        }

        Ok(())
    }

    /// A file with no color channel gets no colors, rather than a black array.
    #[test]
    fn a_scan_without_color_gets_no_color_attribute() -> Result<()> {
        let file = TempFile::new("nocolor", &encode(RES_NM, false, &flat_grid(3, 4)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.point_colors().is_none());

        Ok(())
    }

    /// The whole reason `point_stdev` is a first class field: the uncertainty model's output used
    /// to be returned as a loose parallel array, or discarded entirely.
    #[test]
    fn stores_the_modeled_uncertainty_as_point_stdev() -> Result<()> {
        let file = TempFile::new("stdev", &encode(RES_NM, false, &flat_grid(3, 4)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, Some(&RampModel))?;

        let stdev = mesh.point_stdev().expect("a model was supplied");
        assert_eq!(stdev.len(), mesh.point_count());

        for (p, s) in mesh.points().iter().zip(stdev.iter()) {
            assert_relative_eq!(*s, 0.001 + 0.0001 * p.x, epsilon = 1.0e-12);
        }

        Ok(())
    }

    /// Without a model there is nothing to store, and no reason to invent a value.
    #[test]
    fn no_model_means_no_standard_deviations() -> Result<()> {
        let file = TempFile::new("nostdev", &encode(RES_NM, false, &flat_grid(3, 4)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.point_stdev().is_none());

        Ok(())
    }

    /// A model is client code, and a negative value is not a standard deviation. The container's
    /// validation should surface that rather than let it through.
    #[test]
    fn a_model_returning_an_invalid_deviation_is_an_error() {
        let file = TempFile::new("badmodel", &encode(RES_NM, false, &flat_grid(3, 4)));
        let result = load_lptf3_mesh_data(file.path(), Lptf3Load::All, Some(&NegativeModel));

        assert!(result.is_err());
    }

    /// A single row has no neighbor to triangulate against, so every one of its points is an
    /// orphan and the mesh comes back empty. It must not underflow the strip loop on the way.
    #[test]
    fn a_single_row_scan_loads_as_an_empty_mesh() -> Result<()> {
        let file = TempFile::new("onerow", &encode(RES_NM, false, &flat_grid(1, 6)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.is_empty());

        Ok(())
    }

    /// Points which no face reaches are dropped, and everything indexed by a point is compacted
    /// along with them. The rows here are a connected pair followed by an isolated one.
    #[test]
    fn points_which_end_up_in_no_face_are_dropped() -> Result<()> {
        // Frames 0 and 1 are adjacent and will mesh together. Frame 6 is five row spacings away,
        // past the limit, so its points connect to nothing.
        let rows = vec![Row::flat(0, 5), Row::flat(1, 5), Row::flat(6, 5)];
        let file = TempFile::new("orphan", &encode(RES_NM, true, &rows));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, Some(&RampModel))?;

        // Only the ten points of the first two rows survive.
        assert_eq!(mesh.point_count(), 10);
        assert_eq!(mesh.face_count(), 8);
        assert!(mesh.points().iter().all(|p| p.y < 1.5));

        // The color and standard deviation arrays are compacted to match, not left at the
        // original length or misaligned against the surviving points.
        let colors = mesh
            .point_colors()
            .expect("the file declares a color channel");
        let stdev = mesh.point_stdev().expect("a model was supplied");
        assert_eq!(colors.len(), 10);
        assert_eq!(stdev.len(), 10);

        for ((p, c), s) in mesh.points().iter().zip(colors.iter()).zip(stdev.iter()) {
            assert_eq!(c[0], p.x as u8);
            assert_relative_eq!(*s, 0.001 + 0.0001 * p.x, epsilon = 1.0e-12);
        }

        Ok(())
    }

    /// A file with no frames at all is an empty mesh, which `MeshData3` permits.
    #[test]
    fn an_empty_scan_loads_as_an_empty_mesh() -> Result<()> {
        let file = TempFile::new("empty", &encode(RES_NM, false, &[]));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.is_empty());

        Ok(())
    }

    /// Rows separated by more than twice the nominal row spacing are a gap in the scan, and must
    /// not be bridged with triangles that were never measured. With nothing to join, both rows
    /// are orphaned and the mesh is empty.
    #[test]
    fn rows_further_apart_than_the_spacing_allows_are_not_joined() -> Result<()> {
        let rows = vec![Row::flat(0, 5), Row::flat(5, 5)];
        let file = TempFile::new("gap", &encode(RES_NM, false, &rows));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert_eq!(mesh.face_count(), 0);
        assert!(mesh.is_empty());

        Ok(())
    }

    /// A step in z larger than the world edge ratio permits is a depth discontinuity, and the two
    /// sides of it should not be stitched together.
    #[test]
    fn a_depth_discontinuity_is_not_stitched_across() -> Result<()> {
        // Two rows of five points each, where the second row sits 20 mm below the first. Every
        // candidate face would need an edge of at least 20 mm, well past the limit.
        let rows = vec![
            Row {
                frame: 0,
                points: (0..5).map(|x| (x, 0, 0)).collect(),
            },
            Row {
                frame: 1,
                points: (0..5).map(|x| (x, -20, 0)).collect(),
            },
        ];
        let file = TempFile::new("step", &encode(RES_NM, false, &rows));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert_eq!(mesh.face_count(), 0);
        assert!(mesh.is_empty());

        Ok(())
    }

    /// Every face must index a point which exists, and must not be degenerate. `MeshData3::new`
    /// checks the first, this checks the second.
    #[test]
    fn no_face_is_degenerate() -> Result<()> {
        let file = TempFile::new("degen", &encode(RES_NM, false, &flat_grid(5, 7)));
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.face_count() > 0);
        for f in mesh.faces() {
            assert_ne!(f[0], f[1]);
            assert_ne!(f[1], f[2]);
            assert_ne!(f[2], f[0]);
        }

        Ok(())
    }

    // ===========================================================================================
    // Point data loading
    // ===========================================================================================

    /// The point pathway keeps every measured point, unlike the mesh pathway which drops the ones
    /// that end up in no face.
    #[test]
    fn loads_every_point_of_a_flat_scan() -> Result<()> {
        use crate::io::load_lptf3;

        let file = TempFile::new("cloud-grid", &encode(RES_NM, false, &flat_grid(4, 5)));
        let cloud = load_lptf3(file.path(), Lptf3Load::All)?;

        assert_eq!(cloud.point_count(), 20);

        // Points arrive in scan order: row by row, and within a row in ascending x.
        for row in 0..4 {
            for col in 0..5 {
                let p = cloud.points()[row * 5 + col];
                assert_relative_eq!(p.x, col as f64, epsilon = 1.0e-12);
                assert_relative_eq!(p.y, row as f64, epsilon = 1.0e-12);
                assert_relative_eq!(p.z, 0.0, epsilon = 1.0e-12);
            }
        }

        // Nothing was supplied to attach, so nothing should have been invented.
        assert!(cloud.attrs().is_empty());

        Ok(())
    }

    #[test]
    fn carries_the_color_channel_onto_the_cloud() -> Result<()> {
        use crate::io::load_lptf3;

        let file = TempFile::new("cloud-color", &encode(RES_NM, true, &flat_grid(3, 4)));
        let cloud = load_lptf3(file.path(), Lptf3Load::All)?;

        let colors = cloud
            .point_colors()
            .expect("the file declares a color channel");
        assert_eq!(colors.len(), cloud.point_count());

        // The synthetic rows set the channel to the x index, expanded to a gray triple.
        for row in 0..3 {
            for col in 0..4 {
                let c = colors[row * 4 + col];
                assert_eq!(c, [col as u8; 3]);
            }
        }

        Ok(())
    }

    /// A scan with no color channel must not invent one.
    #[test]
    fn a_scan_without_color_produces_a_bare_cloud() -> Result<()> {
        use crate::io::load_lptf3;

        let file = TempFile::new("cloud-nocolor", &encode(RES_NM, false, &flat_grid(2, 3)));
        let cloud = load_lptf3(file.path(), Lptf3Load::All)?;

        assert!(cloud.point_colors().is_none());
        assert!(cloud.attrs().is_empty());

        Ok(())
    }

    /// The mesh pathway drops the points which end up in no face, so its point buffer is a strict
    /// subset of what the point pathway returns for the same file.
    #[test]
    fn the_mesh_pathway_keeps_no_more_points_than_the_point_pathway() -> Result<()> {
        use crate::io::load_lptf3;

        let file = TempFile::new("cloud-vs-mesh", &encode(RES_NM, false, &flat_grid(4, 5)));

        let cloud = load_lptf3(file.path(), Lptf3Load::All)?;
        let mesh = load_lptf3_mesh_data(file.path(), Lptf3Load::All, None)?;

        assert!(mesh.point_count() <= cloud.point_count());

        Ok(())
    }
}

// ===============================================================================================
// Tests against a real scan
// ===============================================================================================

/// These exercise the loader against `laser-sample.lptf3`, a real laser profile scan rather than a
/// synthetic one. The synthetic tests above pin the behavior that can be reasoned about exactly;
/// these pin the behavior on data with genuine sensor noise, dropouts, and a non-trivial surface.
///
/// The fixture was produced from a proprietary scan by keeping only the points lying within 0.1 mm
/// of a mesh covering the region cleared for publication, which left 87,422 of the original
/// 17,297,108 points across 457 of the original 12,012 frames. The surviving points are the raw
/// values as they were measured, so the noise characteristics are real. Frame indices were
/// preserved, so the gaps left by the removal appear to the loader as gaps in the scan, which is
/// itself worth having in a fixture.
#[cfg(test)]
mod real_scan_tests {
    use super::*;
    use crate::geom3::mesh::algorithms::normals::compute_face_normal;
    use crate::io::row_scan::WORLD_EDGE_RATIO;
    use crate::io::{DiffTanModel, Lptf3DsParams, load_lptf3, load_lptf3_mesh_data};
    use crate::tests::get_test_file_path;
    use approx::assert_relative_eq;
    use std::path::PathBuf;

    /// The parameters of the sensor which produced the fixture, taken from the private alignment
    /// test that was built against this same data.
    fn sensor_model() -> DiffTanModel {
        DiffTanModel::new(230.0, 80.0, 0.00037716)
    }

    fn scan_path() -> PathBuf {
        get_test_file_path("laser-sample.lptf3")
    }

    /// A downsampling configuration matching the private test's, at a coarser move limit.
    fn smooth_params() -> Lptf3DsParams {
        Lptf3DsParams::new(8, 1.0, 1.0, 0.25)
    }

    #[test]
    fn loads_the_real_scan_at_full_resolution() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;

        // 87,422 points are in the file; the 1,805 which end up in no face are dropped.
        assert_eq!(mesh.point_count(), 85_617);
        assert_eq!(mesh.face_count(), 141_554);

        // Every point sits on one of the 457 surviving scan rows.
        let rows = mesh
            .points()
            .iter()
            .map(|p| (p.y * 1.0e6).round() as i64)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(rows.len(), 457);

        Ok(())
    }

    /// The mesh pathway walks the same rows in the same order as the point cloud pathway, then
    /// drops the points which end up in no face. Its point buffer must therefore be a subsequence
    /// of the cloud's: the same values, in the same order, with gaps.
    ///
    /// This is what makes the two comparable, and it is the property that would break if the
    /// flattening or the orphan compaction ever reordered anything.
    #[test]
    fn the_mesh_points_are_a_subsequence_of_the_cloud_points() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;
        let cloud = load_lptf3(&scan_path(), Lptf3Load::All)?;

        assert_eq!(cloud.points().len(), 87_422);
        assert!(mesh.point_count() < cloud.points().len());

        let mut cloud_i = 0;
        for m in mesh.points() {
            // Advance through the cloud until this mesh point is found. Anything skipped was an
            // orphan; running off the end means the order does not match.
            while cloud_i < cloud.points().len() && cloud.points()[cloud_i] != *m {
                cloud_i += 1;
            }
            assert!(
                cloud_i < cloud.points().len(),
                "mesh point {:?} is not in the cloud in order",
                m
            );
            cloud_i += 1;
        }

        Ok(())
    }

    /// The scan carries a laser return intensity per point, which the mesh pathway used to drop.
    #[test]
    fn the_real_scan_carries_its_intensity_channel() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;

        let colors = mesh.point_colors().expect("the scan has a color channel");
        assert_eq!(colors.len(), mesh.point_count());

        // Every channel is the same 8-bit value expanded to gray.
        assert!(colors.iter().all(|c| c[0] == c[1] && c[1] == c[2]));

        // The values must span a real range rather than being a constant, which is what a
        // misaligned read of the point payload would most likely produce.
        let lo = colors.iter().map(|c| c[0]).min().unwrap();
        let hi = colors.iter().map(|c| c[0]).max().unwrap();
        assert_eq!(lo, 5);
        assert_eq!(hi, 255);

        Ok(())
    }

    /// The uncertainty model is the reason `point_stdev` is a first class field. Against the real
    /// sensor parameters it should produce deviations on the order of a couple of microns.
    #[test]
    fn the_sensor_model_lands_on_the_real_scan() -> Result<()> {
        let model = sensor_model();
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, Some(&model))?;

        let stdev = mesh.point_stdev().expect("a model was supplied");
        assert_eq!(stdev.len(), mesh.point_count());

        for (p, s) in mesh.points().iter().zip(stdev.iter()) {
            assert_relative_eq!(*s, model.value(p.x, p.z), epsilon = 1.0e-15);
        }

        // The model depends only on z, and this patch spans about 8 mm of depth, so the deviation
        // varies by only a few percent across it. Roughly 2.6 to 2.8 microns.
        let lo = stdev.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = stdev.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(lo, 0.0026286, epsilon = 1.0e-6);
        assert_relative_eq!(hi, 0.0028173, epsilon = 1.0e-6);

        Ok(())
    }

    /// Taking every Nth row is the cheap downsample, and should thin the scan without moving any
    /// point that survives it.
    #[test]
    fn taking_every_eighth_row_thins_the_scan() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::TakeEveryN(8), None)?;

        assert_eq!(mesh.point_count(), 2_128);
        assert_eq!(mesh.face_count(), 3_998);

        // Still a mesh with attributes intact, not just a bag of points.
        assert_eq!(mesh.point_colors().unwrap().len(), mesh.point_count());

        Ok(())
    }

    /// The smoothing downsample selects the same points as `TakeEveryN`, then moves each one along
    /// z onto a gaussian-weighted fit of its neighborhood. The x and y coordinates are untouched
    /// and the move is bounded by `max_move`.
    #[test]
    fn smoothing_moves_points_along_z_only() -> Result<()> {
        let params = smooth_params();
        let plain =
            load_lptf3_mesh_data(&scan_path(), Lptf3Load::TakeEveryN(params.take_every), None)?;
        let smooth = load_lptf3_mesh_data(&scan_path(), Lptf3Load::SmoothSample(params), None)?;

        assert_eq!(smooth.point_count(), plain.point_count());
        assert_eq!(smooth.face_count(), 3_998);

        let mut moved = 0;
        for (a, b) in plain.points().iter().zip(smooth.points().iter()) {
            assert_relative_eq!(a.x, b.x, epsilon = 1.0e-12);
            assert_relative_eq!(a.y, b.y, epsilon = 1.0e-12);

            let dz = (b.z - a.z).abs();
            assert!(
                dz <= params.max_move + 1.0e-12,
                "a point moved {} mm, past the {} mm limit",
                dz,
                params.max_move
            );
            if dz > 1.0e-9 {
                moved += 1;
            }
        }

        // The smoothing has to actually do something, or this test proves nothing.
        assert!(
            moved > plain.point_count() / 2,
            "only {} of {} points moved",
            moved,
            plain.point_count()
        );

        Ok(())
    }

    /// Real scan data has dropouts and depth steps, so the edge criteria have work to do. Whatever
    /// they let through must still be a valid mesh.
    #[test]
    fn every_face_of_the_real_scan_is_valid() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;

        for f in mesh.faces() {
            assert!(f.iter().all(|i| (*i as usize) < mesh.point_count()));
            assert_ne!(f[0], f[1]);
            assert_ne!(f[1], f[2]);
            assert_ne!(f[2], f[0]);

            // No face may span more than the world edge ratio allows, measured against the
            // nominal row spacing of the full-resolution scan.
            let p = f.map(|i| mesh.points()[i as usize]);
            let longest = (p[0] - p[1])
                .norm()
                .max((p[1] - p[2]).norm())
                .max((p[2] - p[0]).norm());
            assert!(longest < WORLD_EDGE_RATIO * 2.0 * 0.015);
        }

        Ok(())
    }

    /// Every face of a scanned surface must face back toward the sensor, which looks along -z. If
    /// the winding were reversed, or mixed between the two triangles of a quad, this would not
    /// hold. This is the cheapest end to end check that the triangulation is coherent.
    #[test]
    fn the_real_scan_has_coherent_face_winding() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;

        let mut vertical = 0;
        for f in mesh.faces() {
            let p = f.map(|i| mesh.points()[i as usize]);
            let n = compute_face_normal(&p).expect("a face of the scan is degenerate");

            assert!(n.z >= 0.0, "a face points away from the sensor: {:?}", n);
            if n.z == 0.0 {
                vertical += 1;
            }
        }

        // A handful of faces stand exactly on end, because the sensor quantizes x and two points
        // in a row can land on the same value at different depths. Those are legitimate, but they
        // should be a rounding error in the total rather than a structural feature.
        assert!(
            vertical * 1000 < mesh.face_count(),
            "{} of {} faces are exactly vertical",
            vertical,
            mesh.face_count()
        );

        Ok(())
    }

    /// The loader drops the points which end up in no face, so point normals are computable on a
    /// scan mesh without any further cleanup. Before that, `compute_point_normals` refused the
    /// whole mesh, because a point belonging to no face has no defined normal.
    #[test]
    fn the_real_scan_has_no_orphan_points() -> Result<()> {
        let mesh = load_lptf3_mesh_data(&scan_path(), Lptf3Load::All, None)?;

        let mut used = vec![false; mesh.point_count()];
        for f in mesh.faces() {
            for i in f {
                used[*i as usize] = true;
            }
        }
        assert!(used.iter().all(|u| *u));

        let normals = mesh.compute_point_normals()?;
        assert_eq!(normals.len(), mesh.point_count());

        // The sensor looks along -z, so no point of the measured surface may face away from it.
        assert!(normals.iter().all(|n| n.z >= 0.0));

        // A handful sit exactly on end, on points whose every incident face is one of the
        // vertical ones described in `build_faces`.
        let flat = normals.iter().filter(|n| n.z == 0.0).count();
        assert!(
            flat * 1000 < normals.len(),
            "{} of {} point normals are exactly vertical",
            flat,
            normals.len()
        );

        Ok(())
    }
}
