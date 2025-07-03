//! This module has tools for estimating normals from point clouds.

use crate::common::kd_tree::KdTreeSearch;
use crate::{KdTree3, Point3, SvdBasis3, UnitVec3, Vector3};
use rayon::prelude::*;

/// Estimates of normals from a point cloud, including confidence values.
pub struct NormalEstimates {
    pub normals: Vec<UnitVec3>,
    pub confidence: Vec<f64>,
}

///
///
/// # Arguments
///
/// * `points`:
/// * `must_match`:
/// * `tree`:
/// * `radius`:
///
/// returns: NormalEstimates
///
/// # Examples
///
/// ```
///
/// ```
pub fn estimate_by_neighborhood(
    points: &[Point3],
    must_match: &[Vector3],
    tree: &KdTree3,
    radius: f64,
) -> NormalEstimates {
    // We'll do this with a parallel iterator, so we can use rayon.
    let indices = (0..points.len()).collect::<Vec<_>>();
    let mut combined = indices
        .par_iter()
        .map(|&i| {
            let neighbors = tree
                .within(&points[i], radius)
                .iter()
                .map(|(j, _)| *j)
                .collect::<Vec<_>>();
            let (n, c) = svd_normal(&neighbors, points);

            // The normal is pointing in the wrong direction, flip it.
            let checked_n = if n.dot(&must_match[i]) < 0.0 { -n } else { n };

            (i, checked_n, c)
        })
        .collect::<Vec<_>>();

    // Put everything in order by index.
    combined.sort_by(|a, b| a.0.cmp(&b.0));

    let mut normals = Vec::with_capacity(points.len());
    let mut confidence = Vec::with_capacity(points.len());

    for (_, n, c) in combined {
        normals.push(n);
        confidence.push(c);
    }

    NormalEstimates {
        normals,
        confidence,
    }
}

fn svd_normal(neighbors: &[usize], points: &[Point3]) -> (UnitVec3, f64) {
    if neighbors.len() < 3 {
        return (UnitVec3::new_unchecked(Vector3::new(0.0, 0.0, 1.0)), 0.0);
    }
    let working = neighbors.iter().map(|&i| points[i]).collect::<Vec<_>>();

    let svd = SvdBasis3::from_points(&working, None);
    let st_dev = svd.basis_stdevs();
    // TODO: I don't remember where this formula comes from
    let certainty = (st_dev[1] / st_dev[0]) * (st_dev[0] - st_dev[2]) / st_dev[0];
    (UnitVec3::new_normalize(svd.basis[2]), certainty)
}
