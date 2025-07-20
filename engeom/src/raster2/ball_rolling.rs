//! Implementation of a ball rolling background algorithm for scalar raster data.  This algorithm
//! is based on the widely known algorithm first introduced in a 1983 paper by Stanley Sternberg.

use crate::na::DMatrix;
use crate::raster2::{Point2IIndexAccess, ScalarRaster, SizeForIndex};
use crate::{Point2, Result};

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
    let matrix = raster.to_matrix();
    let ball = ball_matrix(radius, raster.px_size);

    let mut rolled = DMatrix::zeros(matrix.nrows(), matrix.ncols());
    rolled.fill(f64::INFINITY);

    let tail_n = (ball.nrows() - 1) as i32 / 2;
    for p in matrix.iter_indices() {
        if !matrix.get_at(p).unwrap().is_finite() {
            continue;
        }

        let (min_mi, min_ki, count_ki) = window_vals(p.y as usize, matrix.nrows(), tail_n);
        let (min_mj, min_kj, count_kj) = window_vals(p.x as usize, matrix.ncols(), tail_n);

        let mut contact_height = f64::INFINITY;

        for j in 0..count_kj {
            for i in 0..count_ki {
                let ki = min_ki + i;
                let kj = min_kj + j;
                let mi = min_mi + i;
                let mj = min_mj + j;
                if ball[(ki, kj)].is_finite() && matrix[(mi, mj)].is_finite() {
                    contact_height = contact_height.min(ball[(ki, kj)] - matrix[(mi, mj)]);
                }
            }
        }

        // The total amount we can sink has been captured, we can now take a bite out of the
        // output matrix
        for j in 0..count_kj {
            for i in 0..count_ki {
                let ki = min_ki + i;
                let kj = min_kj + j;
                let mi = min_mi + i;
                let mj = min_mj + j;
                if ball[(ki, kj)].is_finite() && matrix[(mi, mj)].is_finite() {
                    rolled[(mi, mj)] = rolled[(mi, mj)].min(ball[(ki, kj)] - contact_height); // - matrix[(mi, mj)]);
                }
            }
        }
    }

    let mut result = DMatrix::zeros(matrix.nrows(), matrix.ncols());
    result.fill(f64::INFINITY);
    for p in rolled.iter_indices() {
        let ov = rolled.get_at(p).unwrap();
        if ov.is_finite() {
            let sv = matrix.get_at(p).unwrap();

            result.set_at(p, sv - ov).unwrap();
        }
    }

    Ok(ScalarRaster::from_matrix(
        &result,
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

fn window_vals(a: usize, count: usize, tail_n: i32) -> (usize, usize, usize) {
    let min_mi = (a as i32 - tail_n).max(0) as usize;
    let min_ki = min_mi + (tail_n as usize) - a;
    let max_mi = (a + tail_n as usize).min(count - 1);
    let count_ki = max_mi - min_mi + 1;

    (min_mi, min_ki, count_ki)
}
