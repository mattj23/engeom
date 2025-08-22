use crate::common::{DiscreteDomain, IndexMask};
use crate::io::lptf3::{Lptf3Loader, Lptf3UncertaintyModel};
use crate::io::{Lptf3DsParams, Lptf3Load, load_lptf3_mesh, write_mesh_stl};
use crate::sensors::LaserProfileGeom;
use crate::{Mesh, Point3, PointCloud, Result, SurfacePoint3, UnitVec3};
use parry3d_f64::query::{Ray, RayCast};
use rayon::prelude::*;
use std::path::Path;

struct ProcessedFrame {
    y: f64,
    points: Vec<Point3>,
    normals: Vec<UnitVec3>,
    colors: Vec<u8>,
}

impl ProcessedFrame {
    fn new(y: f64, points: Vec<Point3>, normals: Vec<UnitVec3>, colors: Vec<u8>) -> Self {
        Self {
            y,
            points,
            normals,
            colors,
        }
    }
}

pub fn load_lptf3_comprehensive(
    file_path: &Path,
    uncertainty_model: &dyn Lptf3UncertaintyModel,
    bad_edge_count: usize,
    ray_check: Option<(&LaserProfileGeom, f64)>,
) -> Result<(PointCloud, Vec<f64>)> {
    let base_params = Lptf3Load::SmoothSample(Lptf3DsParams::new(8, 1.5, 1.0, 1.0));
    let half_mesh = load_lptf3_mesh(file_path, base_params)?;
    let mesh = Mesh::try_from(&half_mesh)?;

    let mut loader = Lptf3Loader::new(file_path, None, false)?;

    let mut frames = Vec::new();
    while let Some(frame) = loader.get_next_frame_points()? {
        if frame.points.len() < 2 {
            continue;
        }
        frames.push(frame);
    }

    let mut results = frames
        .par_iter()
        .map(|frame| {
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
            let xs = edge_mask
                .to_indices()
                .iter()
                .map(|&i| frame.points[i].x)
                .collect::<Vec<_>>();
            let xs = DiscreteDomain::try_from(xs).expect("Failed to create discrete domain");
            let mut ci = xs
                .closest_index(frame.points[0].x)
                .expect("Failed to find closest index");

            for i in 0..edge_mask.len() {
                if edge_mask.get(i) {
                    continue;
                }

                let d0 = (frame.points[i].x - xs.values()[ci]).abs();
                let d1 = if ci < xs.values().len() - 1 {
                    (frame.points[i].x - xs.values()[ci + 1]).abs()
                } else {
                    f64::INFINITY
                };
                if d1 < d0 {
                    ci += 1;
                }
                let d = d0.min(d1);

                if d < loader.y_translation * bad_edge_count as f64 {
                    edge_mask.set(i, true);
                    continue;
                }
            }

            // Figure out the position of the scanner head

            let mut points = Vec::new();
            let mut normals = Vec::new();
            let mut colors = Vec::new();
            edge_mask.not_mut();

            for i in edge_mask.to_indices() {
                let p = frame.points[i].at_y(frame.y_pos);

                if let Some((rc, offset)) = ray_check {
                    let detector = Point3::new(0.0, frame.y_pos + rc.detector_y, rc.detector_z);
                    let ray_point = SurfacePoint3::new_normalize(p, detector - p).shift(offset);
                    let ray = Ray::from(&ray_point);
                    let intersect = mesh.tri_mesh().cast_local_ray(&ray, f64::INFINITY, false);
                    if intersect.is_some() {
                        continue;
                    }
                }

                let mp = mesh.surf_closest_to(&p);
                points.push(p);
                normals.push(mp.sp.normal);

                if let Some(color) = frame.points[i].color {
                    colors.push(color);
                }
            }

            ProcessedFrame::new(frame.y_pos, points, normals, colors)
        })
        .collect::<Vec<_>>();

    results.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let mut colors = Vec::new();
    let mut uncertainties = Vec::new();

    for frame in results {
        points.extend(frame.points);
        normals.extend(frame.normals);
        colors.extend(frame.colors.iter().map(|&c| [c; 3]));
    }

    for p in points.iter() {
        uncertainties.push(uncertainty_model.value(p.x, p.z));
    }

    let c = if loader.has_color { Some(colors) } else { None };

    let cloud = PointCloud::try_new(points, Some(normals), c)?;

    Ok((cloud, uncertainties))
}
