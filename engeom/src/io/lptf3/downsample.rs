use crate::io::lptf3::{Lptf3Loader, expand_colors};
use crate::{Point3, PointCloud, Result, SurfacePoint3, UnitVec3, Vector3};
use rayon::prelude::*;
use std::path::Path;

#[derive(Copy, Clone, Debug)]
pub struct Lptf3DsParams {
    pub take_every: u32,
    pub look_scale: f64,
    pub weight_sigma: f64,
    pub max_move: f64,
}

impl Lptf3DsParams {
    pub fn new(take_every: u32, look_scale: f64, weight_sigma: f64, max_move: f64) -> Self {
        Self {
            take_every,
            look_scale,
            weight_sigma,
            max_move,
        }
    }
}

pub fn load_lptf3_downfilter(file_path: &Path, params: Lptf3DsParams) -> Result<PointCloud> {
    let downsampled = load_downsample_filter_lptf3(file_path, params)?;
    let final_points = downsampled.rows.into_iter().flatten().collect::<Vec<_>>();

    let c = if let Some(colors) = downsampled.colors {
        let final_colors = colors.into_iter().flatten().collect::<Vec<_>>();
        Some(expand_colors(&final_colors))
    } else {
        None
    };
    PointCloud::try_new(final_points, None, c)
}

pub struct Lptf3Downsampled {
    pub rows: Vec<Vec<Point3>>,
    pub colors: Option<Vec<Vec<u8>>>,
    pub y_translation: f64,
}

pub fn load_downsample_filter_lptf3(
    file_path: &Path,
    params: Lptf3DsParams,
) -> Result<Lptf3Downsampled> {
    if params.take_every < 2 {
        return Err("take_every must be at least 2".into());
    }

    let mut loader = Lptf3Loader::new(file_path, Some(params.take_every), true)?;

    // Prepare the full point cloud into a set of rows
    // =========================================================================================
    // At the end of this step, `all_points` will contain a vector of vectors, where each inner
    // vector is the 3d point data for a single row of points, sorted in ascending order by x
    // coordinate.  The `row_data` will contain the indices of the points that are destined for the
    // final point cloud. The color vector will have been filled with the color values in the order
    // of the points matching with a flattening of `all_points`.
    let mut all_points = Vec::new();
    let mut all_colors = Vec::new();
    let mut row_data = Vec::new();

    while let Some(full) = loader.get_next_frame_points()? {
        let mut row = Vec::new();
        let mut c_row = Vec::new();
        for p in full.points.iter() {
            row.push(p.at_y(full.y_pos));
            c_row.push(p.color.unwrap_or(0));
        }
        all_points.push(row);
        all_colors.push(c_row);
        row_data.push(full.to_take)
    }

    // Sample the final point cloud
    // =========================================================================================
    // We will iterate through each row of points and for each index in the `row_data` element
    // we will find all points within a sampling distance of the point at that index.  We will
    // perform a gaussian weighted SVD on those points and then correct the working point's z value
    // to intersect with the plane.

    // The number of rows to look forward and backwards when sampling the point cloud.
    let look_rows = (params.take_every as f64 * params.look_scale.abs()).ceil() as i32;
    let look_dist = look_rows as f64 * loader.y_translation * 1.25;
    let weight_sigma = params.weight_sigma * look_dist;

    let working_indices = (0..row_data.len())
        .filter(|&i| !row_data[i].is_empty())
        .collect::<Vec<_>>();

    let mut combined = working_indices
        .par_iter()
        .map(|&row_i| {
            let to_take = &row_data[row_i];

            let mut row_points = Vec::new();
            let mut row_colors = Vec::new();
            for col_i in to_take.iter() {
                // The point that we're working on
                let p = all_points[row_i][*col_i];

                let mut samples = Vec::new();
                for check_i in (row_i as i32 - look_rows)..=(row_i as i32 + look_rows) {
                    if check_i < 0 || check_i >= all_points.len() as i32 {
                        continue; // Skip rows that are out of bounds
                    }
                    let check_row = &all_points[check_i as usize];

                    // Binary search for the first point that is p.x - look_dist
                    let target = p.x - look_dist;
                    let start = check_row
                        .binary_search_by(|a| a.x.total_cmp(&target))
                        .unwrap_or_else(|i| i);

                    for check_p in check_row.iter().skip(start) {
                        let d = (check_p - p).norm();
                        if d <= look_dist {
                            samples.push((*check_p, gaussian_weight(d, weight_sigma)));
                        }
                        if check_p.x > p.x + look_dist {
                            // If the point is beyond the look distance, we can stop checking this row
                            break;
                        }
                    }
                }

                row_points.push(adjust_by_gwm(&p, &samples, params.max_move));
                row_colors.push(all_colors[row_i][*col_i]);
            }

            (row_i, row_points, row_colors)
        })
        .collect::<Vec<_>>();

    // Get everything back in order
    combined.sort_by(|(row_i1, _, _), (row_i2, _, _)| row_i1.cmp(row_i2));
    let mut final_rows = Vec::with_capacity(combined.len());
    let mut final_row_colors = Vec::with_capacity(combined.len());
    for (_, row_points, row_colors) in combined {
        final_rows.push(row_points);
        final_row_colors.push(row_colors);
    }

    let c = if loader.has_color {
        Some(final_row_colors)
    } else {
        None
    };
    Ok(Lptf3Downsampled {
        rows: final_rows,
        colors: c,
        y_translation: loader.y_translation,
    })
}

fn adjust_by_gwm(p: &Point3, samples: &[(Point3, f64)], max_move: f64) -> Point3 {
    if samples.len() < 3 {
        return *p;
    }
    let (points, weights) = samples
        .iter()
        .map(|(pnt, wt)| (*pnt, *wt))
        .unzip::<Point3, f64, Vec<_>, Vec<_>>();

    let sp = SurfacePoint3::new(*p, UnitVec3::new_unchecked(Vector3::z()));

    let mut ts = points
        .iter()
        .map(|pnt| sp.scalar_projection(pnt))
        .collect::<Vec<_>>();
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let total_weight = weights.iter().sum::<f64>();
    let weighted_mean = ts
        .iter()
        .zip(weights.iter())
        .map(|(t, w)| t * w)
        .sum::<f64>()
        / total_weight;

    let weighted_mean = weighted_mean.clamp(-max_move, max_move);
    sp.at_distance(weighted_mean)
}

fn gaussian_weight(x: f64, sigma: f64) -> f64 {
    (-0.5 * (x.powi(2) / sigma.powi(2))).exp()
}
