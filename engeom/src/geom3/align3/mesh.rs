//! This module contains support functionality for mesh-to-mesh alignment, including the core
//! building blocks consumed by other mesh alignment routine in the `align3` module tree.  All
//! public items are re-exported directly from `align3`, so callers interact with them through that
//! parent module.  If you are implementing your own alignment routines, you may find that you want
//! to use some functionality contained here.
//!
//! # What this module provides
//!
//! ## Smart alignment point sampling
//!
//! [`generate_alignment_points`] is the primary function for smart alignment point sampling, but
//! pay close attention to its documentation. It takes a test (moving) mesh and a reference (fixed)
//! mesh and does point sampling that tries to filter for quality sample points on the test mesh
//! which are away from edges, corners, and high curvature, and within an expected distribution
//! of distances.  See the [`GAPParams`] configuration struct.
//!
//! ## Alignment mesh container
//!
//! [`AlignmentMesh`] wraps a `&Mesh` together with optional per-vertex uncertainty values, an
//! initial transform, and a list of [`MeshWeight`] providers.  It is used by the multi-mesh bundle
//! adjustment solver (`multi_mesh`) where each entity needs its own configuration.
//!
//! ## Residual weighting
//!
//! The [`MeshWeight`] trait and its two built-in implementations let callers scale alignment
//! residuals by region:
//!
//! - [`FaceIndexWeight`] — applies a weight to points whose face index falls within a given
//!   [`IndexMask`], leaving all other points at weight `1.0`.
//! - [`NearMeshWeight`] — applies a weight to points that are close to (and similarly oriented as)
//!   a second reference mesh, based on distance and normal angle thresholds.

use crate::common::IndexMask;
use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::{dist, mean_point};
use crate::common::vec_f64::mean_and_stdev;
use crate::geom3::mesh::MeshSurfPoint;
use crate::{Iso3, KdTree3, Mesh, SelectOp, Selection, SvdBasis3, To2D, TransformBy};
use parry2d_f64::transformation::convex_hull;
use std::f64::consts::PI;

// ===============================================================================================
// Alignment mesh container
// ===============================================================================================

/// A container structure that holds all the information necessary to align this mesh against a
/// reference. This provides a unified interface for all additional options used to refine the
/// alignment process, such as the uncertainty of the mesh vertex points, an initial alignment,
/// and methods of applying weights to the sample points.
pub struct AlignmentMesh<'a> {
    pub mesh: &'a Mesh,
    pub uncertainty: Option<&'a [f64]>,
    pub initial: Option<&'a Iso3>,
    pub weights: Option<&'a [Box<dyn MeshWeight + Sync>]>,
}

impl<'a> AlignmentMesh<'a> {
    /// Creates a new `AlignmentMesh` instance.
    ///
    /// # Arguments
    ///
    /// * `mesh`: The mesh to align.
    /// * `uncertainty`: Optional uncertainty values for the mesh vertices, should be in the form
    ///   of a slice of f64 values the same length as the number of vertices in the mesh. The values
    ///   should represent standard deviations of distance the vertex would be from the current
    ///   position upon repeated measurements. A normal distribution is assumed for the sake of
    ///   calculating relative probabilities.
    /// * `initial`: Optional initial transformation for the mesh.  If not specified, the identity
    ///   transformation will be used.
    /// * `weights`: An optional list of weight providing entities that will be used to calculate
    ///   weights of the alignment points _once_ upon initialization. These weights will be combined
    ///   and will then scale the residual calculated at the associated alignment point.
    pub fn new(
        mesh: &'a Mesh,
        uncertainty: Option<&'a [f64]>,
        initial: Option<&'a Iso3>,
        weights: Option<&'a [Box<dyn MeshWeight + Sync>]>,
    ) -> Self {
        Self {
            mesh,
            uncertainty,
            initial,
            weights,
        }
    }

    pub fn transform(&self) -> Iso3 {
        *self.initial.unwrap_or(&Iso3::identity())
    }
}
// ===============================================================================================
// Residual weighting
// ===============================================================================================

/// This is a trait for a generic mesh weight providing entity. When given a `MeshSurfPoint`, it
/// should return a weight that will be applied to the residual at that point during the
/// alignment process. This allows for flexible weighting strategies, such as proximity to another
/// mesh, or being in a specific set of faces, or having a specific direction, etc.
pub trait MeshWeight {
    /// Returns the weight of the mesh point.
    ///
    /// # Arguments
    ///
    /// * `point`: The mesh point for which to compute the weight.
    ///
    /// returns: f64
    fn weight(&self, point: &MeshSurfPoint) -> f64;
}

#[derive(Clone)]
pub struct FaceIndexWeight {
    mask: IndexMask,
    weight: f64,
}

impl FaceIndexWeight {
    /// Creates a new FaceIndexWeight instance.
    ///
    /// # Arguments
    ///
    /// * `mask`: The mask containing the face indices to apply the weight to.
    /// * `weight`: The weight to apply to the points in the specified faces.
    pub fn new(mask: IndexMask, weight: f64) -> Self {
        Self { mask, weight }
    }

    pub fn to_boxed_trait(self) -> Box<dyn MeshWeight + Sync> {
        Box::new(self)
    }
}

impl MeshWeight for FaceIndexWeight {
    fn weight(&self, point: &MeshSurfPoint) -> f64 {
        // If the point's face index is in the mask, return the weight, otherwise return 0.0
        if self.mask.get(point.face_index as usize) {
            self.weight
        } else {
            1.0
        }
    }
}

#[derive(Clone)]
pub struct NearMeshWeight {
    mesh: Mesh,
    weight: f64,
    max_dist: f64,
    max_angle: f64,
}

impl NearMeshWeight {
    /// Creates a new NearMeshWeight instance.
    ///
    /// # Arguments
    ///
    /// * `mesh`: The mesh to compute the weight against.
    /// * `weight`: The weight to apply to the mesh points.
    /// * `max_dist`: The maximum distance to consider for the weight.
    /// * `max_angle`: The maximum angle between normals to consider for the weight.
    pub fn new(mesh: Mesh, weight: f64, max_dist: f64, max_angle: f64) -> Self {
        Self {
            mesh,
            weight,
            max_dist,
            max_angle,
        }
    }

    pub fn to_boxed_trait(self) -> Box<dyn MeshWeight + Sync> {
        Box::new(self)
    }
}

impl MeshWeight for NearMeshWeight {
    fn weight(&self, point: &MeshSurfPoint) -> f64 {
        // Find the nearest point in the mesh to the given point
        let nearest = self.mesh.surf_closest_to(&point.sp.point);
        let dist = dist(&nearest.sp, &point.sp);

        // If the distance is greater than the maximum distance, return 1.0
        if dist > self.max_dist {
            return 1.0;
        }

        // Check the angle between the normals
        if nearest.sp.normal.angle(&point.sp.normal) > self.max_angle {
            return 1.0;
        }

        // Return the weight if all conditions are met
        self.weight
    }
}

// ===============================================================================================
// Alignment point sampling
// ===============================================================================================

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
            max_neighbor_angle: PI / 3.0,
            out_of_plane_ratio: 0.05,
            centroid_ratio: 1.0,
            filter_distances: Some(3.0),
        }
    }
}

pub fn simple_alignment_points(
    test_mesh: &Mesh,
    ref_mesh: &Mesh,
    spacing: f64,
) -> Vec<MeshSurfPoint> {
    let overlap = test_mesh
        .face_select(Selection::None)
        .faces_overlap(ref_mesh, PI / 4.0, 2.0, SelectOp::Add)
        .take_mask();

    if overlap.count_true() == 0 {
        Vec::new()
    } else {
        test_mesh.sample_poisson(spacing, Some(&overlap))
    }
}

/// A sampling algorithm that finds a set of ideal alignment points on a test mesh which can be
/// used to align it with a reference mesh.  Pay close attention to the parameters.
///
/// The method will begin with a Poisson disk sampling of the test mesh and then identify ideal
/// points based on the local neighborhood of each point. Points which are in neighborhoods of low
/// curvature and away from edges are retained, and then the entire neighborhood is projected onto
/// the reference mesh, and the same criteria are applied to the projected neighborhood.
///
/// 1. The 6 nearest neighbors to the point in the test mesh (not including the original point) are
///    found, and they must all be within 2x the sample spacing distance of the original point.
/// 2. The angle between the normals of the original point and each neighbor must be less than the
///    `max_neighbor_angle` parameter.
/// 3. An SVD basis is computed from the neighborhood points, and the maximum out-of-plane distance
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
    let all_points = test_mesh.sample_poisson(params.sample_spacing, None);
    let tree = KdTree3::new(&all_points).expect("KD tree build failed");

    // Now we're going to iterate through the points and find ones which meet the criteria for
    // being paired with the reference mesh.
    let mut candidates = Vec::new();
    for (i, pnt) in all_points.iter().enumerate() {
        // Find the nearest 7 to the point in the reference mesh.
        let nearest = tree.nearest(pnt, 7);

        // Prepare a vec with the neighbor points
        let mut neighbors = Vec::new();
        for (idx, _) in nearest.iter() {
            if *idx != i {
                neighbors.push(all_points[*idx])
            }
        }

        // We'll execute the sample validity check, and if it passes, we'll add the point to
        // the mask
        let (ok, d) = smpl_check(pnt, &neighbors, ref_mesh, iso, params);
        if ok {
            candidates.push((d, pnt))
        }
    }

    // Lastly, we'll filter out candidates with distances more than a specified number of standard
    // deviations from the mean distance.
    if let Some(sigma) = params.filter_distances {
        if candidates.len() > 10 {
            let distances = candidates.iter().map(|(d, _)| *d).collect::<Vec<_>>();
            if let Some((mean, st_dev)) = mean_and_stdev(&distances) {
                let threshold = mean + sigma * st_dev;
                candidates.retain(|(d, _)| *d < threshold);
            }
        }
    }

    candidates.into_iter().map(|(_, p)| *p).collect()
}

/// Sample validity check, which looks at the conditions around the sample check point based on
/// its neighbors. Internally, this uses the `sac_check` acceptability test on a test point and
///  its neighbors.
///
/// First, it does this on the check point and its neighbors on the test mesh, then it projects
/// those points onto the reference mesh and performs the same check.
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

/// Sample acceptability check, which examines a bunch of requirements around a check point using
/// only its neighboring points. Checks include the minimum number of neighbors, their spacing
/// around the check point, and how they project into a local basis.
pub fn sac_check(
    check_point: &MeshSurfPoint,
    neighbors: &[MeshSurfPoint],
    params: &GAPParams,
) -> bool {
    // We require some minimum number of neighbors to do the downstream checks, we're using 5 for now
    if neighbors.len() < 5 {
        return false;
    }

    // We check that all the neighbors are within 2x the sample spacing distance of the check point,
    // that the angle between the normals of the check point and each neighbor is less than the
    // `max_neighbor_angle` parameter, and that the scalar projection of each neighbor onto the
    // check point is less than the sample spacing (preventing large out-of-plane neighbors).
    for n in neighbors.iter() {
        if dist(n, check_point) > params.sample_spacing * 2.0
            || n.sp.normal.angle(&check_point.sp.normal) > params.max_neighbor_angle
            || check_point.sp.scalar_projection(&n.sp.point).abs() > params.sample_spacing
        {
            return false;
        }
    }

    // We'll turn the neighbors into a collection of simple points and then compute a basis
    // coordinate system where we'll perform the final checks.
    let mut points = neighbors.iter().map(|n| n.sp.point).collect::<Vec<_>>();
    points.push(check_point.sp.point);

    // If we can't compute a basis, the points are too far apart or too noisy, and we can simply
    // decide that this sample fails the check.
    let Some(basis) = SvdBasis3::from_points(&points, None) else {
        return false;
    };

    // Now we'll move the points to the basis coordinate system. The x-y plane should roughly
    // represent the surface of the mesh around the check point, and the z-axis will roughly
    // be the surface normal direction. The origin is the centroid of the neighbors + the
    // original check point.
    let iso = Iso3::from(&basis);
    let points = (&points).transform_by(&iso);

    // Check that all the points are no more than out_of_plane_ratio * sample_spacing from the
    // x-y plane.
    if points
        .iter()
        .map(|p| p.z.abs())
        .any(|z| z > params.out_of_plane_ratio * params.sample_spacing)
    {
        return false;
    }

    // Finally, we're going to project everything the 2d plane and check that the centroid of
    // the convex hull of the collection of all points is within centroid_ratio * sample_spacing
    // of the 2d projection of the check point.
    let check2 = (iso * check_point.sp.point).to_2d();
    let points2 = points.to_2d();

    let centroid = mean_point(&convex_hull(&points2));
    if dist(&centroid, &check2) > params.sample_spacing * params.centroid_ratio {
        return false;
    }

    true
}
