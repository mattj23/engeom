use crate::common::IndexMask;
use crate::io::lptf3::{Lptf3Loader, Lptf3UncertaintyModel};
use crate::io::{Lptf3DsParams, Lptf3Load, load_lptf3_mesh, write_mesh_stl};
use crate::{Mesh, PointCloud, Result};
use std::path::Path;
use stl_io::write_stl;

pub fn load_lptf3_comprehensive(
    file_path: &Path,
    uncertainty_model: &dyn Lptf3UncertaintyModel,
    bad_edge_count: usize,
    emitter_z: f64,
    detector_y: f64,
    detector_z: f64,
) -> Result<(PointCloud, Vec<f64>)> {
    let base_params = Lptf3Load::SmoothSample(Lptf3DsParams::new(8, 1.5, 1.0, 1.0));
    let half_mesh = load_lptf3_mesh(file_path, base_params)?;
    let mesh = Mesh::try_from(&half_mesh)?;

    let mut loader = Lptf3Loader::new(file_path, None, false)?;
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let mut colors = Vec::new();
    let mut uncertainties = Vec::new();

    while let Some(frame) = loader.get_next_frame_points()? {
        if frame.points.len() < 2 {
            continue;
        }

        let mut edge_mask = IndexMask::new(frame.points.len(), false);
        edge_mask.set(0, true);
        edge_mask.set(frame.points.len() - 1, true);

        for i0 in 0..frame.points.len() - 2 {
            let i1 = i0 + 1;
            let p0 = frame.points[i0].as_point2();
            let p1 = frame.points[i1].as_point2();

            if (p1 - p0).norm() > loader.y_translation * 100.0 {
                edge_mask.set(i0, true);
                edge_mask.set(i1, true);
            }
        }

        // Pass to mark neighbors of bad edges
        let mut xs = edge_mask
            .to_indices()
            .iter()
            .map(|&i| frame.points[i].x)
            .collect::<Vec<_>>();

        for i in 0..edge_mask.len() {
            if edge_mask.get(i) {
                continue;
            }
            let x_dist = xs
                .iter()
                .map(|x| (frame.points[i].x - *x).abs())
                .min_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap());
            if let Some(d) = x_dist {
                if d.abs() < loader.y_translation * bad_edge_count as f64 {
                    edge_mask.set(i, true);
                }
            }
        }

        edge_mask.not_mut();
        for i in edge_mask.to_indices() {
            let p = frame.points[i].at_y(frame.y_pos);
            let mp = mesh.surf_closest_to(&p);
            let u = uncertainty_model.value(p.x, p.z);
            uncertainties.push(u);
            points.push(p);
            normals.push(mp.sp.normal);

            if let Some(color) = frame.points[i].color {
                colors.push([color; 3]);
            }
        }
    }

    let c = if loader.has_color { Some(colors) } else { None };

    let cloud = PointCloud::try_new(points, Some(normals), c)?;

    Ok((cloud, uncertainties))
}
