//! This module has a RANSAC implementation for creating lines

use crate::Line2;
use crate::Result;
use crate::common::PCoords;
use crate::common::ransac_tools::ransac_indices;
use itertools::Itertools;

impl Line2 {
    pub fn try_new_ransac(
        points: &[impl PCoords<2>],
        threshold: f64,
        iterations: Option<usize>,
    ) -> Result<Line2> {
        let iterations = iterations.unwrap_or(500);

        let mut best_line = None;
        let mut best_inlier_count = 0;
        let indices = ransac_indices::<2>(iterations, points.len())?;
        for i in indices {
            let line = Line2::from_points(&points[i[0]], &points[i[1]]);
            let mut count = 0;
            for p in points {
                if line.distance_to(p) < threshold {
                    count += 1;
                }
            }
            if count > best_inlier_count {
                best_inlier_count = count;
                best_line = Some(line);
            }
        }

        best_line.ok_or("No good candidate found".into())
    }
}
