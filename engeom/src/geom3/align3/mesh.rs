//! This module has some common abstractions and tools for aligning meshes

use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::{dist, mean_point};
use crate::common::vec_f64::mean_and_stdev;
use crate::geom3::mesh::MeshSurfPoint;
use crate::{Iso3, KdTree3, Mesh, SvdBasis3, To2D, TransformBy};
use parry2d_f64::transformation::convex_hull;
use std::num::NonZero;

#[derive(Debug, Clone, Copy)]
pub struct GAPParams {
    pub sample_spacing: f64,
    pub max_neighbor_angle: f64,
    pub out_of_plane_ratio: f64,
    pub centroid_ratio: f64,
    pub filter_distances: Option<f64>,
}

impl GAPParams {
    /// Creates a new set of parameters for the mesh sampling algorithm.
    ///
    /// # Arguments
    ///
    /// * `sample_spacing`: The spacing to use for the Poisson disk sampling of the test mesh. This
    ///   will also be used as a base value for `out_of_plane_ratio` and `centroid_ratio`. Smaller
    ///   values mean a more dense sampling, but also more points to check and more influence of smaller
    ///   features in the mesh.
    /// * `max_neighbor_angle`: The maximum permissible angle between the normals of any test point and
    ///   its pre-filtered closest 6 neighbors. Leave this large if you have a noisy mesh. A default
    ///   value of `PI / 3.0` (60 degrees) is a good starting point.
    /// * `out_of_plane_ratio`: The maximum permissible out-of-plane distance of any point in the
    ///   neighborhood after the SVD basis is computed. This is a ratio of the sample spacing, so a
    ///   value of `0.5` means that the maximum out-of-plane distance is half of the sample spacing.
    ///   Smaller values enforce more flatness, while larger values allow for more curvature. A
    ///   reasonable starting point of 0.1 to 0.05 will work for many engineering shapes.
    /// * `centroid_ratio`: The maximum permissible distance of the centroid of the neighborhood's 2D
    ///   convex hull points to the 2D projection of the test point. Smaller values require test point
    ///   to be closer to the center of the neighborhood, away from edges.  A value of `0.5` to `1.0`
    ///   is a good starting point.
    /// * `filter_distances`: An optional value that, if provided, will result in a filtering operation
    ///   of the final selected candidates based on the distance between them and the reference mesh.
    ///   A value of `Some(3.0)` will filter out candidates that are more than 3 standard deviations
    ///   above the mean candidate distance to the reference mesh. This can be used to remove outlying
    ///   areas of the test mesh as the two meshes begin to converge towards alignment.
    ///
    /// returns: MshSmParams
    pub fn new(
        sample_spacing: f64,
        max_neighbor_angle: f64,
        out_of_plane_ratio: f64,
        centroid_ratio: f64,
        filter_distances: Option<f64>,
    ) -> Self {
        Self {
            sample_spacing,
            max_neighbor_angle,
            out_of_plane_ratio,
            centroid_ratio,
            filter_distances,
        }
    }

    /// Creates a new set of default parameters for the mesh sampling algorithm, requiring only the
    /// sample spacing to be specified.
    ///
    /// # Arguments
    ///
    /// * `sample_spacing`: the spacing to use for the Poisson disk sampling of the test mesh.
    ///
    /// returns: MshSmParams
    pub fn defaults(sample_spacing: f64) -> Self {
        Self {
            sample_spacing,
            max_neighbor_angle: std::f64::consts::PI / 3.0,
            out_of_plane_ratio: 0.05,
            centroid_ratio: 1.0,
            filter_distances: Some(3.0),
        }
    }
}

/// A sampling algorithm that finds a set of ideal alignment points on a test mesh which can be
/// used to align it with a reference mesh.  Pay close attention to the parameters.
///
/// The method will begin with a Poisson disk sampling of the test mesh, and then identify ideal
/// points based on the local neighborhood of each point. Points which are in neighborhoods of low
/// curvature and away from edges are retained, and then the entire neighborhood is projected onto
/// the reference mesh and the same criteria are applied to the projected neighborhood.
///
/// 1. The 6 nearest neighbors to the point in the test mesh (not including the original point) are
///    found, and they must all be within 2x the sample spacing distance of the original point.
/// 2. The angle between the normals of the original point and each neighbor must be less than the
///    `max_neighbor_angle` parameter.
/// 3. A SVD basis is computed from the neighborhood points, and the maximum out-of-plane distance
///    of each point must be less than the `out_of_plane_ratio` parameter times the sample spacing.
/// 4. The centroid of the neighborhood's 2D convex hull points must be within
///    `centroid_ratio * sample_spacing` of the 2D test point.
///
/// If all of these criteria are met, the neighborhood is projected onto the reference mesh using
/// the provided `Iso3` transform, and the same criteria are applied to the projected points.
/// Additionally, the test point and its corresponding projected point must have normals facing in
/// the same direction.
///
/// At the very end, if a sigma value is provided in `filter_distances`, the mean and standard
/// deviation of the distance from each test point candidate to its corresponding projected point
/// are computed, and any candidates with a distance more than `sigma` standard deviations from
/// the mean distance are filtered out.
///
/// # Arguments
///
/// * `test_mesh`: The mesh which will be sampled for alignment points: the resulting points will be
///   on the surface of this mesh.
/// * `ref_mesh`: The mesh which is used as a reference for the alignment. The resulting points
///   will be ideal points to use when aligning the test mesh to this reference mesh.
/// * `iso`: An initial `Iso3` transform that will be applied to the test mesh points before
///   projecting them onto the reference mesh. This would represent an initial guess of the
///   alignment that is to follow.
/// * `params`: The parameters for the sampling algorithm: see the `MshSmParams` struct for details.
///
/// returns: Vec<MeshSurfPoint, Global>
pub fn generate_alignment_points(
    test_mesh: &Mesh,
    ref_mesh: &Mesh,
    iso: &Iso3,
    params: &GAPParams,
) -> Vec<MeshSurfPoint> {
    // We start with a Poisson disk sampling of the test mesh to get a set of points that are
    // well distributed across the surface and spaced at a roughly known distance.
    let all_points = test_mesh.sample_poisson(params.sample_spacing);
    let tree = KdTree3::new(&all_points);

    // Now we're going to iterate through the points and find ones which meet the criteria for
    // being paired with the reference mesh.
    let mut candidates = Vec::new();
    for (i, pnt) in all_points.iter().enumerate() {
        // Find the nearest 7 to the point in the reference mesh.
        let nearest = tree.nearest(pnt, NonZero::new(7).unwrap());

        // Prepare a vec with the neighbor points
        let mut neighbors = Vec::new();
        for (idx, _) in nearest.iter() {
            if *idx != i {
                neighbors.push(all_points[*idx])
            }
        }

        // We'll execute the sample validity check, and if it passes, we'll add the point to
        // the mask
        let (ok, d) = smpl_check(pnt, &neighbors, ref_mesh, iso, &params);
        if ok {
            candidates.push((d, pnt))
        }
    }

    // Lastly, we'll filter out candidates more than 3 standard deviations beyond the mean distance
    // to the reference mesh.
    if let Some(sigma) = params.filter_distances {
        if candidates.len() > 10 {
            let distances = candidates.iter().map(|(d, _)| *d).collect::<Vec<_>>();
            if let Some((mean, stdev)) = mean_and_stdev(&distances) {
                let threshold = mean + sigma * stdev;
                candidates.retain(|(d, _)| *d < threshold);
            }
        }
    }

    candidates.into_iter().map(|(_, p)| *p).collect()
}

fn smpl_check(
    check: &MeshSurfPoint,
    neighbors: &[MeshSurfPoint],
    reference: &Mesh,
    iso: &Iso3,
    params: &GAPParams,
) -> (bool, f64) {
    // Actual points check
    if !sac_check(check, neighbors, params) {
        return (false, f64::MAX);
    }

    // If the points on the test mesh pass, we project the points to the reference mesh and
    // run the same check.
    let moved = iso * check.sp;
    let check_ref = reference.surf_closest_to(&moved.point);

    // Normals must be facing the same direction
    if check_ref.sp.normal.dot(&moved.normal) < 0.0 {
        return (false, f64::MAX);
    }

    let neighbors_ref = neighbors
        .iter()
        .map(|sp| reference.surf_closest_to(&(iso * sp.sp.point)))
        .collect::<Vec<_>>();

    // The minimum spacing to the check_ref point should be max_spacing
    for sp in &neighbors_ref {
        if dist(sp, &check_ref) < params.sample_spacing {
            return (false, f64::MAX);
        }
    }

    (
        sac_check(&check_ref, &neighbors_ref, params),
        dist(&moved, &check_ref),
    )
}

pub fn sac_check(
    check_point: &MeshSurfPoint,
    neighbors: &[MeshSurfPoint],
    params: &GAPParams,
) -> bool {
    if neighbors.len() < 5 {
        return false;
    }

    for n in neighbors.iter() {
        if dist(n, check_point) > params.sample_spacing * 2.0
            || n.sp.normal.angle(&check_point.sp.normal) > params.max_neighbor_angle
            || check_point.sp.scalar_projection(&n.sp.point).abs() > params.sample_spacing
        {
            return false;
        }
    }

    // We'll turn the neighbors into a collection of simple points and then compute a
    // basis to perform the final checks
    let mut points = neighbors.iter().map(|n| n.sp.point).collect::<Vec<_>>();
    points.push(check_point.sp.point);

    let basis = SvdBasis3::from_points(&points, None);
    let iso = Iso3::from(&basis);

    // Now we'll move the points to the basis coordinates
    let points = (&points).transform_by(&iso);

    if points
        .iter()
        .map(|p| p.z.abs())
        .any(|z| z > params.out_of_plane_ratio * params.sample_spacing)
    {
        return false;
    }

    let check2 = (iso * check_point.sp.point).to_2d();
    let points2 = points.to_2d();

    let centroid = mean_point(&convex_hull(&points2));
    if dist(&centroid, &check2) > params.sample_spacing * params.centroid_ratio {
        return false;
    }

    true
}
