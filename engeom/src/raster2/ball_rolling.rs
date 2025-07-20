//! Implementation of a ball rolling background algorithm for scalar raster data.  This algorithm
//! is based on the widely known algorithm first introduced in a 1983 paper by Stanley Sternberg.

use crate::na::DMatrix;
use crate::raster2::{Point2IIndexAccess, ScalarRaster, SizeForIndex};
use crate::{Point2, Result};
use rayon::prelude::*;

/// This function performs a ball rolling algorithm to compute a background for a raster of
/// scalar values. The algorithm was first introduced by Stanley Sternberg in 1983 and is widely
/// used in image processing for background estimation.
///
/// In the case that the raster of scalar data represents a height field, the ball rolling
/// operation has a physical interpretation. Each point in the resulting scalar field represents
/// the lowest height that the surface of the ball reached over that point, and the combined set
/// of points is the hull of the ball over the entire surface.
///
/// If subtracted from the original raster, the result is a distance map, showing how close the
/// ball's surface got to the original point. Areas where the ball contacted will have a value of
/// zero, while scratches, grooves, and pits will have negative values that correspond with
/// actual physical distances.
///
/// # Arguments
///
/// * `raster`:
/// * `radius`:
///
/// returns: Result<ScalarRaster, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn ball_rolling_background(raster: &ScalarRaster, radius: f64) -> Result<ScalarRaster> {
    let ball = ball_matrix(radius, raster.px_size);
    let mut matrix = raster.to_matrix();
    for p in matrix.iter_indices() {
        if !raster.mask.get_point(p) {
            matrix.set_at(p, f64::NEG_INFINITY)?
        }
    }
    let tail_n = (ball.nrows() - 1) as i32 / 2;

    let xi_work = (-tail_n..(raster.width() as i32 + tail_n)).collect::<Vec<_>>();

    let collected = xi_work
        .par_chunks(xi_work.len() / rayon::current_num_threads())
        .map(|slice| {
            let mut rolled = DMatrix::zeros(matrix.nrows(), matrix.ncols());
            rolled.fill(f64::INFINITY);

            for xi in slice {
                for yi in -tail_n..(raster.height() as i32 + tail_n) {
                    let (min_mi, min_ki, count_ki) = wv(yi, matrix.nrows(), tail_n);
                    let (min_mj, min_kj, count_kj) = wv(*xi, matrix.ncols(), tail_n);

                    let mut contact_height = f64::INFINITY;

                    for j in 0..count_kj {
                        for i in 0..count_ki {
                            let ki = min_ki + i;
                            let kj = min_kj + j;
                            let mi = min_mi + i;
                            let mj = min_mj + j;
                            contact_height = contact_height.min(ball[(ki, kj)] - matrix[(mi, mj)]);
                        }
                    }

                    if !contact_height.is_finite() {
                        continue;
                    }

                    for j in 0..count_kj {
                        for i in 0..count_ki {
                            let ki = min_ki + i;
                            let kj = min_kj + j;
                            let mi = min_mi + i;
                            let mj = min_mj + j;
                            rolled[(mi, mj)] =
                                rolled[(mi, mj)].min(ball[(ki, kj)] - contact_height);
                        }
                    }
                }
            }

            rolled
        })
        .collect::<Vec<_>>();

    let mut rolled = DMatrix::zeros(matrix.nrows(), matrix.ncols());
    rolled.fill(f64::INFINITY);

    for r in collected {
        for j in 0..r.nrows() {
            for i in 0..r.ncols() {
                rolled[(j, i)] = rolled[(j, i)].min(r[(j, i)]);
            }
        }
    }

    for p in rolled.iter_indices() {
        if !raster.mask.get_point(p) {
            rolled.set_at(p, f64::NAN)?;
        }
    }

    Ok(ScalarRaster::from_matrix(
        &rolled,
        raster.px_size,
        raster.min_z,
        raster.max_z,
    ))
}

fn ball_matrix(radius: f64, px_size: f64) -> DMatrix<f64> {
    let px_dia = (2.0 * radius / px_size).ceil() as usize | 1;

    let mut matrix = DMatrix::zeros(px_dia, px_dia);
    let c = Point2::new(radius, radius);

    for p in matrix.iter_indices() {
        let a = Point2::new(p.x as f64, p.y as f64) * px_size;
        let b = (a - c).norm();
        let v = if b >= radius {
            f64::INFINITY
        } else {
            radius - (radius.powi(2) - b.powi(2)).sqrt()
        };
        matrix.set_at(p, v).unwrap()
    }

    matrix
}

fn wv(a: i32, count: usize, tail_n: i32) -> (usize, usize, usize) {
    let min_mi = (a - tail_n).max(0);
    let min_ki = min_mi + tail_n - a;
    let max_mi = (a + tail_n).min(count as i32 - 1);
    let count_ki = max_mi - min_mi + 1;

    (min_mi as usize, min_ki as usize, count_ki as usize)
}
