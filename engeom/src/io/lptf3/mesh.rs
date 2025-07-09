use crate::Point3;
use crate::Result;
use crate::common::triangulation::parallel_row2::{StripRowPoint, build_parallel_row_strip};
use crate::geom3::mesh::HalfEdgeMesh;
use crate::io::Lptf3Load;
use crate::io::lptf3::{Lptf3UncertaintyModel, get_downfilter_point_rows, get_loader_point_rows};
use std::path::Path;

pub fn load_lptf3_mesh_core(
    file_path: &Path,
    load: Lptf3Load,
    model: Option<&dyn Lptf3UncertaintyModel>,
) -> Result<(HalfEdgeMesh, Option<Vec<f64>>)> {
    let (take_every, y_shift, point_rows) = match load {
        Lptf3Load::All => get_loader_point_rows(file_path, None),
        Lptf3Load::TakeEveryN(n) => get_loader_point_rows(file_path, Some(n)),
        Lptf3Load::SmoothSample(params) => get_downfilter_point_rows(file_path, params),
    }?;

    // Set the edge ratios for the strip and world triangulation
    let strip_r = 2.0; // The maximum edge ratio for the strip triangulation.
    let world_r = 5.0; // The maximum edge ratio for world triangulation.

    let max_spacing = take_every as f64 * y_shift * 2.0;
    let mut uncert = Vec::new();

    // First build the mesh vertices and the corresponding rows of strip row points
    let mut mesh = HalfEdgeMesh::new();
    let mut strip_rows = Vec::new();
    for row in point_rows.iter() {
        let mut strip_row = Vec::new();
        for p in row.iter() {
            if let Some(m) = model {
                uncert.push(m.value(p.x, p.z));
            }

            let ih = mesh
                .add_vertex(p.coords)
                .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
            strip_row.push(StripRowPoint::new(p.x, ih));
        }
        strip_rows.push(strip_row);
    }

    // Now iterate through the strip rows and build the mesh
    for row_i in 0..strip_rows.len() - 1 {
        if point_rows[row_i].is_empty() || point_rows[row_i + 1].is_empty() {
            continue; // Skip empty rows
        }

        let y0 = point_rows[row_i][0].y;
        let y1 = point_rows[row_i + 1][0].y;

        // If the rows are too far apart, skip the triangulation
        if (y1 - y0).abs() > max_spacing {
            continue;
        }

        let row0 = &strip_rows[row_i];
        let row1 = &strip_rows[row_i + 1];

        // Build the strip triangulation between the two rows
        let r = build_parallel_row_strip(row0, y0, row1, y1, strip_r)?;
        for (i0, i1, i2) in r {
            // Check the edge ratio on actual points
            let pa: Point3 = mesh
                .point(i0)
                .map_err(|e| format!("Failed to get point {}: {:?}", i0, e))?
                .into();
            let pb: Point3 = mesh
                .point(i1)
                .map_err(|e| format!("Failed to get point {}: {:?}", i1, e))?
                .into();
            let pc: Point3 = mesh
                .point(i2)
                .map_err(|e| format!("Failed to get point {}: {:?}", i2, e))?
                .into();
            let ea = (pa - pb).norm();
            let eb = (pb - pc).norm();
            let ec = (pc - pa).norm();

            let edge_ratio = ea.max(eb).max(ec) / max_spacing;
            if edge_ratio < world_r {
                mesh.add_tri_face(i1, i0, i2)
                    .map_err(|e| format!("Failed to add face: {:?}", e))?;
            }
        }
    }

    if model.is_some() {
        Ok((mesh, Some(uncert)))
    } else {
        Ok((mesh, None))
    }
}
