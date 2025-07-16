//! This module contains tools for performing Poisson disk sampling on points in 2D and 3D.

use crate::common::kd_tree::{KdTreeSearch, PartialKdTree};
use crate::common::{IndexMask, PCoords, voxel_downsample};

/// Performs Poisson disk sampling on a set of points in D-dimensional space, returning a mask
/// indicating which points are retained.
///
/// Internally, this function first performs a voxel downsampling of the points at 1/2 the size
/// of the Poisson disk radius. This is a quick operation that guarantees that at most there will
/// be about 16 points within the radius of any point in the pre-sampled set. This pre-sampling
/// reduces the size of the KD tree that needs to be built, and prevents the issue with Kiddo not
/// returning all points within the radius when the number of points is too large from breaking
/// the sampling process.
///
/// # Arguments
///
/// * `points`: A slice of points implementing the `PCoords` trait for the specified dimension `D`.
/// * `radius`: The radius of the Poisson disk sampling. This value should be positive and non-zero.
///
/// returns: IndexMask
pub fn sample_poisson_disk_all<const D: usize>(
    points: &[impl PCoords<D>],
    radius: f64,
) -> IndexMask {
    let pre_mask = voxel_downsample(points, radius * 0.25);
    let partial_tree = PartialKdTree::new(points, &pre_mask);

    let mut skip_mask = pre_mask.clone();
    let mut final_mask = IndexMask::new(points.len(), false);

    for i in pre_mask.iter_true() {
        if !skip_mask.get(i) {
            continue;
        }

        final_mask.set(i, true);

        let within = partial_tree.within(&points[i], radius);
        for w in within {
            skip_mask.set(w.0, false);
        }
    }

    final_mask
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point2;
    use crate::common::kd_tree::KdTree;
    use crate::na::Point;
    use rand;
    use rand::Rng;

    #[test]
    fn stress_test_poisson_disk() {
        let n = 5000;
        let mx = 10.0;
        let r = 0.2;

        for _ in 0..100 {
            let points = random_points(n, mx);
            // let mut indices = index_vec(None, n);
            // indices.shuffle(&mut rand::rng());

            let keep = sample_poisson_disk_all(&points, r);
            let at_least = (mx * mx) / (r * r) * 0.25;
            assert!(keep.len() > at_least as usize);

            // Brute force check that each point only has one point (itself) within the radius
            let kept = keep.clone_indices_of(&points).unwrap();
            let tree = KdTree::new(&kept);
            for (i, &p) in kept.iter().enumerate() {
                let within = tree.within(&p, r);
                assert_eq!(within.len(), 1);
                assert_eq!(within[0].0, i);
            }
        }
    }

    fn random_points(n: usize, mx: f64) -> Vec<Point<f64, 2>> {
        let mut rng = rand::rng();
        (0..n)
            .map(|_| Point2::new(rng.random_range(0.0..mx), rng.random_range(0.0..mx)))
            .collect()
    }
}
