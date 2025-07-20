//! Implementation of a ball rolling background algorithm for scalar raster data.  This algorithm
//! is based on the widely known algorithm first introduced in a 1983 paper by Stanley Sternberg.

use std::path::Path;
use colorgrad::preset::turbo;
use crate::na::DMatrix;
use crate::raster2::roi::RoiOverlay;
use crate::raster2::{Point2I, Point2IIndexAccess, RasterRoi, ScalarImage, ScalarRaster, SizeForIndex, Vector2I};
use crate::{Point2, Result};
use imageproc::definitions::Image;

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
    ball.render_with_cmap(&Path::new("D:/temp/k/ball.png"), &turbo(), Some((-2.0, 2.0)))?;

    // The mask will be the same as the original raster, so we only need to create a buffer to
    // hold the hull values
    // let mut hull = ScalarImage::new(raster.width(), raster.height());
    let mut hull = ScalarRaster::filled_like(&raster, u16::MAX);

    let tail_n = (ball.width() - 1) as i32 / 2;
    let tail_v = Vector2I::new(tail_n, tail_n);
    let full_roi = RasterRoi::new(
        Point2I::origin(),
        Point2I::new(raster.width() as i32, raster.height() as i32),
    );

    for p in ExpandedIter::from_raster(&raster, tail_n) {
        let roi = RasterRoi::new(p - tail_v, p + tail_v);
        let overlay = RoiOverlay::new(roi, full_roi);

        let mut contact_height = u16::MAX;

        for q in overlay.iter_intersection_a() {
            if let (Some(b), Some(m)) = (ball.u_at(q.local), raster.u_at(q.parent)) {
                contact_height = contact_height.min(b - m);
            }
        }

        for q in overlay.iter_intersection_a() {
            if let (Some(b), Some(v)) = (ball.u_at(q.local), hull.u_at(q.parent)) {
                hull.set_u_at(q.parent, Some(v.min(b - contact_height)))?;
            }
        }
    }

    Ok(hull)
}

struct ExpandedIter {
    min_x: i32,
    min_y: i32,
    max_x: i32,
    max_y: i32,
    current_x: i32,
    current_y: i32,
}

impl ExpandedIter {
    fn from_raster(raster: &ScalarRaster, expand: i32) -> Self {
        Self::new(
            -expand,
            -expand,
            raster.width() as i32 + expand - 1,
            raster.height() as i32 + expand - 1,
        )
    }
    fn new(min_x: i32, min_y: i32, max_x: i32, max_y: i32) -> Self {
        Self {
            min_x,
            min_y,
            max_x,
            max_y,
            current_x: min_x,
            current_y: min_y,
        }
    }
}

impl Iterator for ExpandedIter {
    type Item = Point2I;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_y > self.max_y {
            return None;
        }

        let point = Point2I::new(self.current_x, self.current_y);

        // Move to the next point
        self.current_x += 1;
        if self.current_x > self.max_x {
            self.current_x = self.min_x;
            self.current_y += 1;
        }

        Some(point)
    }
}

fn indices(
    i: usize,
    j: usize,
    min_mi: usize,
    min_ki: usize,
    min_mj: usize,
    min_kj: usize,
) -> (usize, usize, usize, usize) {
    let ki = min_ki + i;
    let kj = min_kj + j;
    let mi = min_mi + i;
    let mj = min_mj + j;

    (ki, kj, mi, mj)
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

fn window_vals(a: i32, count: usize, tail_n: i32) -> (usize, usize, usize) {
    let min_mi = (a - tail_n).max(0);
    let min_ki = min_mi + tail_n - a;
    let max_mi = (a + tail_n).min(count as i32 - 1);
    let count_ki = max_mi - min_mi + 1;

    (min_mi as usize, min_ki as usize, count_ki as usize)
}
