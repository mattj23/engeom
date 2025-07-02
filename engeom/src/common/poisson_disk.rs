//! This module contains tools for performing Poisson disk sampling on points in 2D and 3D.

use crate::common::kd_tree::{KdTree, KdTreeSearch};
use parry3d_f64::na::Point;

pub fn sample_poisson_disk<const D: usize>(
    all_points: &[Point<f64, D>],
    working_indices: &[usize],
    radius: f64,
) -> Vec<usize> {
    /*
       We're going to do this work with a mask array and a kd-tree.

       First we'll reduce the `all_points` array to just the points at the `working_indices`,
       producing the `working_points` array. The point at `working_points[m]` will be the same as
       the point at `all_points[working_indices[m]]`, where `m` goes from 0 ->
       `working_indices.len()`. The outside index `i` is `working_indices[m]`.

       The mask will be the same length as the `working_indices` array. Each mask element at index
       `m` corresponds with the point at `working_indices[m]`. If the mask element is true, then
       point index `working_indices[m]` (the point at `working_points[m]`) is a valid candidate for
       inclusion in the final set.

       The kd-tree will be built with the `working_points` array. When we query the kd-tree, the
       results will be indices into the `working_points` array.

       We will iterate over the mask array. We skip any mask elements that are false. For each
       mask element that is true, we query the kd-tree for all points within `radius` of the point
       and set the mask elements corresponding to those points to false. After visiting a true
       mask element, we can add the index `working_indices[m]` to the result set.
    */

    let mut results = Vec::new();

    let working_points = working_indices
        .iter()
        .map(|i| all_points[*i])
        .collect::<Vec<_>>();
    let mut mask = vec![true; working_indices.len()];
    let tree = KdTree::new(&working_points);

    for (m, &i) in working_indices.iter().enumerate() {
        if !mask[m] {
            continue;
        }
        results.push(i);
        let within = tree.within(&working_points[m], radius);
        for w in within {
            mask[w.0] = false;
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point2;
    use crate::common::indices::index_vec;
    use rand;
    use rand::Rng;
    use rand::prelude::SliceRandom;

    #[test]
    fn stress_test_poisson_disk() {
        let n = 5000;
        let mx = 10.0;
        let r = 0.2;

        for _ in 0..100 {
            let points = random_points(n, mx);
            let mut indices = index_vec(None, n);
            indices.shuffle(&mut rand::rng());

            let keep = sample_poisson_disk(&points, &indices, r);
            let at_least = (mx * mx) / (r * r) * 0.25;
            assert!(keep.len() > at_least as usize);

            // Brute force check that each point only has one point (itself) within the radius
            let kept = keep.iter().map(|i| points[*i]).collect::<Vec<_>>();
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
