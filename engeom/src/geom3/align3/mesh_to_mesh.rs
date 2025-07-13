use crate::common::DistMode;
use crate::geom3::Align3;
use crate::geom3::align3::mesh::generate_alignment_points;
use crate::geom3::align3::{GAPParams, points_to_mesh};
use crate::{Iso3, Mesh, Result};

/// Perform an iterative alignment of one mesh to another. Each iteration is a full
/// Levenberg-Marquardt optimization of the alignment of a set of specific points sampled from the
/// moving mesh, with the set of points being updated each iteration based on where the previous
/// alignment ended.
///
/// The alignment points are chosen using the `sample_alignment_candidates` method of the `Mesh`
/// struct, which performs an extremely selective sampling of points on the surface of the moving
/// mesh at areas of low curvature, away from corners and edges, and with projections onto the
/// reference mesh that have the same restrictive properties. Lastly, points with a distance more
/// than 3 standard deviations from the mean deviation from the reference mesh are pruned from the
/// sample set.
///
/// This results in a set of alignment points that are selectively distributed in areas of the
/// moving mesh that are most likely to match closely to the reference mesh as the two move, in an
/// attempt to pre-select away areas which are most likely to not overlap between the two meshes.
///
/// With each iteration, the two meshes should become more and more aligned, and the alignment
/// candidates will better reflect the areas of overlap between the two meshes. As this process
/// continues, the average residual distance between the moving mesh and the reference mesh will
/// converge and then slightly increase as more and more points are included in the comparison.
/// This plateau/climb in average residual is the condition used to terminate the iterations.
///
/// # Arguments
///
/// * `moving`: the mesh which will be moved into alignment
/// * `reference`: the mesh which is stationary
/// * `sample_spacing`: a Poisson disk sampling spacing used to pick the initial candidates for the
///   alignment points. The spacing will then be used to derive a set of physical limits for the
///   quality of the alignment candidates. See the `sample_alignment_candidates` method of the
///   `Mesh` struct for more details.
/// * `initial`: an initial guess for the alignment transform
/// * `mode`: the distance mode to use for the alignment. This can be either `DistMode::ToPoint` or
///   `DistMode::ToPlane`. The two are identical except for cases where an alignment point drifts
///   off the edge of the closest triangle on the reference mesh.
/// * `max_iter`: if the alignment does not converge after this many iterations, the function will
///   return an error. This is to prevent infinite loops in case the alignment never converges.
///
/// returns: Result<Alignment<Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn mesh_to_mesh_iterative(
    moving: &Mesh,
    reference: &Mesh,
    initial: &Iso3,
    mode: DistMode,
    max_iter: usize,
    params: &GAPParams,
) -> Result<Align3> {
    let mut last_residual = f64::MAX;
    let mut iter = 0;

    loop {
        let test_points = generate_alignment_points(moving, reference, initial, params);
        if test_points.len() < 5 {
            return Err(format!(
                "Failed on iteration {iter}, not enough alignment candidate \
            points were found to align the meshes."
            )
            .into());
        }
        let points = test_points
            .into_iter()
            .map(|p| p.sp.point)
            .collect::<Vec<_>>();

        let result = points_to_mesh(&points, reference, initial, mode)?;
        let avg = result.avg_residual();
        let improvement = (last_residual - avg) / last_residual;
        // println!("Iteration {iter}, improvement: {improvement:.6}, avg residual: {}",
        //          result.avg_residual());
        if improvement < 0.01 {
            return Ok(result);
        }
        last_residual = avg;
        iter += 1;

        if iter >= max_iter {
            return Err(format!("Failed to converge after {max_iter} iterations.").into());
        }
    }
}
