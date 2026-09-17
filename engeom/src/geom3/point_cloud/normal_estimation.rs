//! Estimating surface normals at the points of a cloud by fitting a plane to each neighborhood.
//!
//! # Axes first, directions second
//!
//! A plane fit through a neighborhood gives an *axis* without a direction. The two unit vectors
//! along the plane's smallest principal direction fit the data equally well, and local geometry
//! cannot identify which one points out of the surface. [`estimate_axes_by_neighborhood`] performs
//! the fit and leaves the sign arbitrary. The functions in [`super::orientation`] choose a sign as
//! a separate step.
//!
//! [`estimate_by_neighborhood`] is the two steps together for the case where the caller already
//! knows roughly which way each normal should face, such as a scan whose sensor position is known.

use crate::common::kd_tree::KdTreeSearch;
use crate::{KdTree3, Point3, SvdBasis3, UnitVec3, Vector3};
use rayon::prelude::*;

/// Estimates of normals from a point cloud, including confidence values.
pub struct NormalEstimates {
    /// The estimated normal at each point, in the order of the points they came from.
    pub normals: Vec<UnitVec3>,

    /// How plane-like each neighborhood was. Values are normally between zero and one. A fully
    /// degenerate neighborhood in which all principal deviations are zero can produce `NaN`. See
    /// [`estimate_axes_by_neighborhood`] for the definition and interpretation.
    pub confidence: Vec<f64>,
}

/// Fit a plane to the neighbors of each point within `radius` and report the plane's normal axis.
///
/// The direction of each result is arbitrary, since a plane fit cannot resolve it. Pass the result
/// through one of the functions in [`super::orientation`] to settle that.
///
/// # The confidence value
///
/// The confidence is `(s1 / s0) * (s0 - s2) / s0`, where `s0 >= s1 >= s2` are the standard
/// deviations of the neighborhood along its three principal directions. It is the product of two
/// conditions required for a meaningful normal:
///
/// - `s1 / s0` is near one when the neighborhood spreads out in two directions and near zero when
///   it falls along a line. A line has no plane through it, so the smallest principal direction is
///   determined by noise.
/// - `(s0 - s2) / s0` is near one when the neighborhood is flat and near zero when it fills space
///   in all three directions and does not define a plane.
///
/// The value is therefore low on edges, corners, sparse regions, and neighborhoods that include
/// two surfaces. It is dimensionless and is intended for relative ordering, not as a probability.
/// If all three deviations are zero, the divisions produce `NaN`; callers should reject non-finite
/// confidence values.
///
/// # Arguments
///
/// * `points`: the points to estimate at, which the tree must have been built over
/// * `tree`: a k-d tree over `points`
/// * `radius`: the neighborhood radius, in the units of the points. It must include several points
///   without crossing a feature.
///
/// returns: NormalEstimates
pub fn estimate_axes_by_neighborhood(
    points: &[Point3],
    tree: &KdTree3,
    radius: f64,
) -> NormalEstimates {
    // Rayon preserves the order of an indexed parallel iterator, so the results line up with the
    // points without a later sort.
    let combined: Vec<(UnitVec3, f64)> = (0..points.len())
        .into_par_iter()
        .map(|i| {
            let neighbors = tree
                .within(&points[i], radius)
                .iter()
                .map(|(j, _)| *j)
                .collect::<Vec<_>>();
            svd_normal(&neighbors, points)
        })
        .collect();

    let mut normals = Vec::with_capacity(points.len());
    let mut confidence = Vec::with_capacity(points.len());

    for (n, c) in combined {
        normals.push(n);
        confidence.push(c);
    }

    NormalEstimates {
        normals,
        confidence,
    }
}

/// Fit a plane to the neighbors of each point within `radius`, taking the direction of each normal
/// from a vector the caller supplies.
///
/// Each estimated axis is flipped where needed so that it makes an angle of less than ninety
/// degrees with the matching entry of `must_match`. The usual source of those vectors is the
/// direction back toward the sensor that measured the point, for which
/// [`super::orientation::must_match_toward`] builds the array.
///
/// # Arguments
///
/// * `points`: the points to estimate at, which the tree must have been built over
/// * `must_match`: one vector per point, giving the side each normal should end up on. A zero
///   vector leaves the sign of that normal arbitrary.
/// * `tree`: a k-d tree over `points`
/// * `radius`: the neighborhood radius, in the units of the points
///
/// returns: NormalEstimates
pub fn estimate_by_neighborhood(
    points: &[Point3],
    must_match: &[Vector3],
    tree: &KdTree3,
    radius: f64,
) -> NormalEstimates {
    let mut estimates = estimate_axes_by_neighborhood(points, tree, radius);

    for (n, m) in estimates.normals.iter_mut().zip(must_match.iter()) {
        if n.dot(m) < 0.0 {
            *n = -*n;
        }
    }

    estimates
}

/// The smallest principal axis of a neighborhood, with a measure of how plane-like it was.
///
/// Fewer than three points cannot define a plane. These neighborhoods receive `+Z` at zero
/// confidence. Callers that cannot accept a placeholder normal should filter by confidence.
fn svd_normal(neighbors: &[usize], points: &[Point3]) -> (UnitVec3, f64) {
    if neighbors.len() < 3 {
        return (UnitVec3::new_unchecked(Vector3::new(0.0, 0.0, 1.0)), 0.0);
    }
    let working = neighbors.iter().map(|&i| points[i]).collect::<Vec<_>>();

    let svd = SvdBasis3::from_points(&working, None).unwrap();
    let st_dev = svd.basis_stdevs();
    let certainty = (st_dev[1] / st_dev[0]) * (st_dev[0] - st_dev[2]) / st_dev[0];
    (UnitVec3::new_normalize(svd.basis[2]), certainty)
}
