//! Implementation of a ball rolling background algorithm for scalar raster data.  This algorithm
//! is based on the widely known algorithm first introduced in a 1983 paper by Stanley Sternberg.

use crate::na::DMatrix;
use crate::raster2::roi::RoiOverlay;
use crate::raster2::{
    Point2I, Point2IIndexAccess, RasterRoi, ScalarImage, ScalarRaster, SizeForIndex, Vector2I,
};
use crate::{Point2, Result};
use colorgrad::preset::turbo;
use imageproc::definitions::Image;
use std::path::Path;
use crate::image::Luma;

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
    let ball = ball_buffer(radius, &raster);
    ball.render_with_cmap(
        &Path::new("D:/temp/k/ball.png"),
        &turbo(),
        Some((-2.0, 2.0)),
    )?;

    // The mask will be the same as the original raster, so we only need to create a buffer to
    // hold the hull values
    // let mut hull = ScalarImage::new(raster.width(), raster.height());
    let mut hull = ScalarImage::new(raster.width(), raster.height());
    hull.fill(u16::MAX);

    let tail_n = (ball.width() - 1) as i32 / 2;

    for xi in -tail_n..(hull.width() as i32 + tail_n) {
        for yi in -tail_n..(hull.height() as i32 + tail_n) {
            let mut contact_height = u16::MAX;

            let (min_mi, min_ki, count_ki) = window_vals(yi, hull.height(), tail_n);
            let (min_mj, min_kj, count_kj) = window_vals(xi, hull.width(), tail_n);

            for j in 0..count_kj {
                for i in 0..count_ki {
                    let ki = min_ki + i;
                    let kj = min_kj + j;
                    let mi = min_mi + i;
                    let mj = min_mj + j;

                    if ball.mask.buffer.get_pixel(kj, ki)[0] == 255
                        && raster.mask.buffer.get_pixel(mj, mi)[0] == 255
                    {
                        let b = ball.values.get_pixel(kj, ki)[0];
                        let m = raster.values.get_pixel(mj, mi)[0];
                        contact_height = contact_height.min(b - m)
                    }
                }
            }

            if contact_height == u16::MAX {
                continue; // No contact, skip this point
            }

            for j in 0..count_kj {
                for i in 0..count_ki {
                    let ki = min_ki + i;
                    let kj = min_kj + j;
                    let mi = min_mi + i;
                    let mj = min_mj + j;

                    if ball.mask.buffer.get_pixel(kj, ki)[0] == 255
                        && raster.mask.buffer.get_pixel(mj, mi)[0] == 255
                    {
                        let b = ball.values.get_pixel(kj, ki)[0];
                        let v = hull.get_pixel(mj, mi)[0];
                        hull.put_pixel(mj, mi, Luma([v.min(b - contact_height)]));
                    }
                }
            }
        }
    }

    ScalarRaster::try_new(hull, raster.mask.clone(), raster.px_size, raster.min_z, raster.max_z)
}

fn ball_buffer(radius: f64, target: &ScalarRaster) -> ScalarRaster {
    let px_dia = (2.0 * radius / target.px_size).ceil() as usize | 1;

    let mut matrix = DMatrix::zeros(px_dia, px_dia);
    let c = Point2::new(radius, radius);

    for p in matrix.iter_indices() {
        let a = Point2::new(p.x as f64, p.y as f64) * target.px_size;
        let b = (a - c).norm();
        let v = if b >= radius {
            f64::INFINITY
        } else {
            radius - (radius.powi(2) - b.powi(2)).sqrt()
        };
        matrix.set_at(p, v).unwrap()
    }

    ScalarRaster::from_matrix(&matrix, target.px_size, target.min_z, target.max_z)
}

fn window_vals(a: i32, count: u32, tail_n: i32) -> (u32, u32, u32) {
    let min_mi = (a - tail_n).max(0);
    let min_ki = min_mi + tail_n - a;
    let max_mi = (a + tail_n).min(count as i32 - 1);
    let count_ki = max_mi - min_mi + 1;

    (min_mi as u32, min_ki as u32, count_ki as u32)
}
