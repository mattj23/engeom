//! Sampling points from the surface of a triangle mesh.
//!
//! This module provides three sampling strategies. Each is available as a free function over raw
//! point and face buffers and as methods on both `MeshData3` and `Mesh3`. None uses the BVH: every
//! sample is computed from the three vertices of its face, and Poisson thinning operates on the
//! samples themselves.
//!
//! - **dense**: a barycentric grid on every face, guaranteeing that no point on the surface is
//!   farther than `max_spacing` from a sample. The samples are not evenly distributed and can be
//!   arbitrarily close together where faces are small.
//! - **poisson**: a dense sampling thinned so that no two samples are closer than `radius`. This is
//!   the sampler to use when a roughly even spacing matters, such as for alignment.
//! - **uniform**: `n` samples drawn at random with probability proportional to area. This is the
//!   only sampler whose output is not deterministic.
//!
//! A fourth sampler, one point per occupied voxel of a regular grid, lives in
//! [`voxelize`](super::voxelize).
//!
//! # Two tiers of result
//!
//! Each sampler comes in two forms. The `sample_surface_*` methods return [`MeshSurfPoint`]s that
//! retain the face index and barycentric coordinates of every sample, allowing callers to trace
//! samples back to the mesh. The `sample_*` methods return a [`PointCloud3`] containing only
//! positions and normals, which is more convenient for most callers and works directly with the
//! point-cloud tools.

use super::{Mesh3, MeshData3, MeshSurfPoint};
use crate::common::IndexMask;
use crate::common::barycentric::{barycentric_grid, barycentric_point};
use crate::common::points::dist;
use crate::common::poisson_disk::sample_poisson_disk_all;
use crate::geom3::point_cloud::PointCloud3;
use crate::{Point3, Result, SurfacePoint3};
use parry3d_f64::shape::Triangle;

// ===============================================================================================
// Public entry points
// ===============================================================================================

impl MeshData3 {
    /// Sample the surface densely, keeping the face index and barycentric coordinates of every
    /// sample.
    ///
    /// Every face is covered by a barycentric grid whose pitch is at most `max_spacing`, ensuring
    /// that no point on the surface is farther than that from a sample. A face whose edges are all
    /// shorter than `max_spacing` contributes only its centroid. There is no lower bound on the
    /// distance between samples, so density varies with face size; use the Poisson sampler when
    /// even spacing matters.
    ///
    /// Degenerate faces, meaning those with no computable normal, are skipped.
    ///
    /// # Arguments
    ///
    /// * `max_spacing`: the maximum distance from any point on the surface to a sample, which
    ///   must be finite and positive
    /// * `face_mask`: optional mask whose length must match the face count. When provided, sampling
    ///   is restricted to faces set to true.
    ///
    /// returns: `Result<Vec<MeshSurfPoint>>`
    pub fn sample_surface_dense(
        &self,
        max_spacing: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<Vec<MeshSurfPoint>> {
        compute_dense_surface_points(self.points(), self.faces(), max_spacing, face_mask)
    }

    /// Sample the surface with a Poisson disk spacing, keeping the face index and barycentric
    /// coordinates of every sample.
    ///
    /// A dense sampling at half the radius is thinned so that no two samples are closer than
    /// `radius`. The thinning is deterministic for a given mesh, so repeated calls give the same
    /// result.
    ///
    /// # Arguments
    ///
    /// * `radius`: the minimum distance between any two samples, which must be finite and
    ///   positive
    /// * `face_mask`: optional mask whose length must match the face count. When provided, sampling
    ///   is restricted to faces set to true.
    ///
    /// returns: `Result<Vec<MeshSurfPoint>>`
    pub fn sample_surface_poisson(
        &self,
        radius: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<Vec<MeshSurfPoint>> {
        compute_poisson_surface_points(self.points(), self.faces(), radius, face_mask)
    }

    /// Sample the surface densely, returning a cloud of positions and normals.
    ///
    /// See [`MeshData3::sample_surface_dense`] for the sampling behavior. This form discards the
    /// link back to the faces.
    ///
    /// # Arguments
    ///
    /// * `max_spacing`: the maximum distance from any point on the surface to a sample, which
    ///   must be finite and positive
    /// * `face_mask`: optional mask whose length must match the face count. When provided, sampling
    ///   is restricted to faces set to true.
    ///
    /// returns: `Result<PointCloud3>`
    pub fn sample_dense(
        &self,
        max_spacing: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        Ok(to_cloud(
            &self.sample_surface_dense(max_spacing, face_mask)?,
        ))
    }

    /// Sample the surface with a Poisson disk spacing, returning a cloud of positions and normals.
    ///
    /// See [`MeshData3::sample_surface_poisson`] for the sampling behavior. This form discards the
    /// link back to the faces.
    ///
    /// # Arguments
    ///
    /// * `radius`: the minimum distance between any two samples, which must be finite and
    ///   positive
    /// * `face_mask`: optional mask whose length must match the face count. When provided, sampling
    ///   is restricted to faces set to true.
    ///
    /// returns: `Result<PointCloud3>`
    pub fn sample_poisson(
        &self,
        radius: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        Ok(to_cloud(&self.sample_surface_poisson(radius, face_mask)?))
    }

    /// Draw `n` random samples from the surface with sampling probability proportional to area,
    /// returning a cloud of positions and normals.
    ///
    /// This uses the thread's random number generator and is not deterministic. Degenerate faces
    /// have zero area and are never selected.
    ///
    /// # Arguments
    ///
    /// * `n`: the number of samples to draw
    ///
    /// returns: `PointCloud3`
    pub fn sample_uniform(&self, n: usize) -> PointCloud3 {
        to_cloud(&compute_uniform_surface_points(
            self.points(),
            self.faces(),
            n,
        ))
    }
}

impl Mesh3 {
    /// Sample the surface densely, keeping the face index and barycentric coordinates of every
    /// sample.
    ///
    /// This delegates to [`MeshData3::sample_surface_dense`].
    pub fn sample_surface_dense(
        &self,
        max_spacing: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<Vec<MeshSurfPoint>> {
        compute_dense_surface_points(self.points(), self.faces(), max_spacing, face_mask)
    }

    /// Sample the surface with a Poisson disk spacing, keeping the face index and barycentric
    /// coordinates of every sample.
    ///
    /// This delegates to [`MeshData3::sample_surface_poisson`].
    pub fn sample_surface_poisson(
        &self,
        radius: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<Vec<MeshSurfPoint>> {
        compute_poisson_surface_points(self.points(), self.faces(), radius, face_mask)
    }

    /// Sample the surface densely, returning a cloud of positions and normals.
    ///
    /// This delegates to [`MeshData3::sample_dense`].
    pub fn sample_dense(
        &self,
        max_spacing: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        Ok(to_cloud(
            &self.sample_surface_dense(max_spacing, face_mask)?,
        ))
    }

    /// Sample the surface with a Poisson disk spacing, returning a cloud of positions and normals.
    ///
    /// This delegates to [`MeshData3::sample_poisson`].
    pub fn sample_poisson(
        &self,
        radius: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        Ok(to_cloud(&self.sample_surface_poisson(radius, face_mask)?))
    }

    /// Draw `n` random samples from the surface with sampling probability proportional to area,
    /// returning a cloud of positions and normals.
    ///
    /// This delegates to [`MeshData3::sample_uniform`].
    pub fn sample_uniform(&self, n: usize) -> PointCloud3 {
        to_cloud(&compute_uniform_surface_points(
            self.points(),
            self.faces(),
            n,
        ))
    }
}

// ===============================================================================================
// Free functions over raw buffers
// ===============================================================================================

/// Densely sample the given triangles using a barycentric grid on every face with a pitch of at
/// most `max_spacing`, ensuring that no point on the surface is farther than that from a sample.
///
/// This is the primitive behind [`MeshData3::sample_surface_dense`]; see there for the guarantees.
///
/// # Arguments
///
/// * `points`: the vertex buffer
/// * `faces`: triangles as indices into `points`
/// * `max_spacing`: the maximum distance from any point on the surface to a sample, which must be
///   finite and positive
/// * `face_mask`: optional mask whose length must match `faces`. When provided, sampling is
///   restricted to faces set to true.
///
/// returns: `Result<Vec<MeshSurfPoint>>`
pub fn compute_dense_surface_points(
    points: &[Point3],
    faces: &[[u32; 3]],
    max_spacing: f64,
    face_mask: Option<&IndexMask>,
) -> Result<Vec<MeshSurfPoint>> {
    validate(faces, max_spacing, "Maximum spacing", face_mask)?;

    let mut sampled = Vec::with_capacity(faces.len());

    for (face_i, vert) in faces.iter().enumerate() {
        if let Some(face_mask) = face_mask
            && !face_mask.get(face_i)
        {
            continue;
        }

        let tri = triangle(points, vert);
        let Some(normal) = tri.normal() else {
            continue;
        };
        let face_index = face_i as u32;

        let mut push = |bc: [f64; 3]| {
            let at = barycentric_point(&tri.a, &tri.b, &tri.c, bc);
            let sp = SurfacePoint3::new(at, normal);
            sampled.push(MeshSurfPoint { face_index, bc, sp });
        };

        if dist(&tri.a, &tri.b) < max_spacing
            && dist(&tri.a, &tri.c) < max_spacing
            && dist(&tri.b, &tri.c) < max_spacing
        {
            // If every pair of vertices is closer than the maximum spacing, sample only the face
            // centroid. The centroid of an equally sized neighboring triangle should be within
            // the maximum spacing of this triangle's centroid.
            push([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]);
        } else {
            for bc in barycentric_grid(&tri.a, &tri.b, &tri.c, max_spacing) {
                push(bc);
            }
        }
    }

    Ok(sampled)
}

/// Sample the given triangles with a Poisson disk spacing by thinning a dense sampling at half the
/// radius until no two samples are closer than `radius`.
///
/// This is the primitive behind [`MeshData3::sample_surface_poisson`].
///
/// # Arguments
///
/// * `points`: the vertex buffer
/// * `faces`: triangles as indices into `points`
/// * `radius`: the minimum distance between any two samples, which must be finite and positive
/// * `face_mask`: optional mask whose length must match `faces`. When provided, sampling is
///   restricted to faces set to true.
///
/// returns: `Result<Vec<MeshSurfPoint>>`
pub fn compute_poisson_surface_points(
    points: &[Point3],
    faces: &[[u32; 3]],
    radius: f64,
    face_mask: Option<&IndexMask>,
) -> Result<Vec<MeshSurfPoint>> {
    validate(faces, radius, "Radius", face_mask)?;

    let starting = compute_dense_surface_points(points, faces, radius * 0.5, face_mask)?;
    let mask = sample_poisson_disk_all(&starting, radius);
    Ok(mask.iter_true().map(|i| starting[i]).collect())
}

/// Draw `n` random samples from the given triangles with sampling probability proportional to area.
///
/// This is the primitive behind [`MeshData3::sample_uniform`]. It uses the thread's random number
/// generator and is not deterministic. Degenerate faces have zero area and are never selected, so
/// an input with no usable area yields an empty result.
///
/// # Arguments
///
/// * `points`: the vertex buffer
/// * `faces`: triangles as indices into `points`
/// * `n`: the number of samples to draw
///
/// returns: `Vec<MeshSurfPoint>`
pub fn compute_uniform_surface_points(
    points: &[Point3],
    faces: &[[u32; 3]],
    n: usize,
) -> Vec<MeshSurfPoint> {
    // Only faces with a usable normal participate. Each is recorded with the cumulative area up to
    // and including that face, allowing a uniform draw over the total area to map to a face by
    // binary search.
    let mut candidates = Vec::with_capacity(faces.len());
    let mut total_area = 0.0;
    for (face_i, vert) in faces.iter().enumerate() {
        let tri = triangle(points, vert);
        let area = tri.area();
        if let Some(normal) = tri.normal()
            && area > 0.0
        {
            total_area += area;
            candidates.push((face_i as u32, total_area, tri, normal));
        }
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        let r = rand::random::<f64>() * total_area;
        let i = candidates
            .partition_point(|(_, cumulative, _, _)| *cumulative <= r)
            .min(candidates.len() - 1);
        let (face_index, _, tri, normal) = &candidates[i];

        // Square-root parameterization gives a uniform distribution over the triangle's area.
        let r1 = rand::random::<f64>();
        let r2 = rand::random::<f64>();
        let s = r1.sqrt();
        let bc = [1.0 - s, s * (1.0 - r2), s * r2];
        let at = barycentric_point(&tri.a, &tri.b, &tri.c, bc);
        result.push(MeshSurfPoint {
            face_index: *face_index,
            bc,
            sp: SurfacePoint3::new(at, *normal),
        });
    }

    result
}

// ===============================================================================================
// Helpers
// ===============================================================================================

/// Validate the arguments shared by the dense and Poisson samplers.
///
/// A non-positive spacing would reach `barycentric_grid` as a divisor and request a grid of
/// unbounded size. A mask shorter than the face buffer would be indexed out of bounds. Both are
/// therefore rejected before sampling starts.
///
/// # Arguments
///
/// * `faces`: the face buffer the mask is checked against
/// * `spacing`: the sampler's spacing argument, which must be finite and positive
/// * `label`: how to name that argument in the error message
/// * `face_mask`: the optional mask to check
///
/// returns: `Result<()>`
fn validate(
    faces: &[[u32; 3]],
    spacing: f64,
    label: &str,
    face_mask: Option<&IndexMask>,
) -> Result<()> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(format!("{label} must be finite and positive, got {spacing}").into());
    }

    if let Some(mask) = face_mask
        && mask.len() != faces.len()
    {
        return Err(format!(
            "A face mask of length {} does not match a mesh with {} faces",
            mask.len(),
            faces.len()
        )
        .into());
    }

    Ok(())
}

fn triangle(points: &[Point3], vert: &[u32; 3]) -> Triangle {
    Triangle::new(
        points[vert[0] as usize],
        points[vert[1] as usize],
        points[vert[2] as usize],
    )
}

/// Build the cloud form from the surface form, omitting the link back to the faces.
fn to_cloud(samples: &[MeshSurfPoint]) -> PointCloud3 {
    let sps: Vec<SurfacePoint3> = samples.iter().map(|m| m.sp).collect();
    PointCloud3::from_surface_points(&sps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KdTree3;
    use crate::common::kd_tree::*;

    #[test]
    fn check_kiddo_bug() {
        let mesh = Mesh3::create_sphere(100.0, 0.011).unwrap();
        let r = 5.0;
        let sampled = mesh.sample_poisson(r, None).unwrap();

        let tree = KdTree3::try_new(sampled.points()).expect("Tree construction failed");
        for p in sampled.points() {
            let neighbors = tree.within(p, r);
            assert_eq!(neighbors.len(), 1, "Missed duplicate");
        }
    }

    /// A non-positive spacing reaches `barycentric_grid` as a divisor, where it would request a
    /// grid of unbounded size rather than simply producing a wrong answer.
    #[test]
    fn a_nonsense_spacing_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(mesh.sample_dense(bad, None).is_err(), "{bad} was accepted");
            assert!(
                mesh.sample_poisson(bad, None).is_err(),
                "{bad} was accepted"
            );
        }
    }

    /// A mask shorter than the face buffer used to be indexed out of bounds, panicking partway
    /// through the sampling rather than reporting a bad argument.
    #[test]
    fn a_mask_of_the_wrong_length_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        for len in [mesh.faces().len() - 1, mesh.faces().len() + 1] {
            let wrong = IndexMask::new(len, true);
            assert!(mesh.sample_dense(0.25, Some(&wrong)).is_err());
            assert!(mesh.sample_poisson(0.25, Some(&wrong)).is_err());
        }
    }

    #[test]
    fn sample_poisson_index_mask_restricts_to_masked_faces() {
        let mesh = Mesh3::create_sphere(100.0, 1.1).unwrap();
        let n_faces = mesh.faces().len();

        // Mask only the first half of the faces
        let masked_indices: Vec<usize> = (0..n_faces / 2).collect();
        let mask = IndexMask::try_from_indices(&masked_indices, n_faces).unwrap();

        let sampled = mesh.sample_surface_poisson(5.0, Some(&mask)).unwrap();

        assert!(!sampled.is_empty(), "Expected at least one sample");
        for mp in &sampled {
            assert!(
                mask.get(mp.face_index as usize),
                "face_index {} is not in the mask",
                mp.face_index
            );
        }
    }

    #[test]
    fn both_containers_sample_identically() {
        let data = MeshData3::create_sphere(20.0, 0.5).unwrap();
        let mesh = Mesh3::from_data(data.clone(), false).unwrap();

        let from_data = data.sample_surface_dense(2.0, None).unwrap();
        let from_mesh = mesh.sample_surface_dense(2.0, None).unwrap();
        assert_eq!(from_data.len(), from_mesh.len());
        for (a, b) in from_data.iter().zip(&from_mesh) {
            assert_eq!(a.face_index, b.face_index);
            assert_eq!(a.bc, b.bc);
            assert_eq!(a.sp.point, b.sp.point);
            assert_eq!(a.sp.normal, b.sp.normal);
        }

        let from_data = data.sample_poisson(3.0, None).unwrap();
        let from_mesh = mesh.sample_poisson(3.0, None).unwrap();
        assert_eq!(from_data.points(), from_mesh.points());
        assert_eq!(from_data.point_normals(), from_mesh.point_normals());
    }

    #[test]
    fn the_cloud_tier_carries_the_surface_tier_positions_and_normals() {
        let mesh = Mesh3::create_box(10.0, 6.0, 4.0, false);

        let surface = mesh.sample_surface_dense(1.0, None).unwrap();
        let cloud = mesh.sample_dense(1.0, None).unwrap();

        assert_eq!(cloud.point_count(), surface.len());
        let normals = cloud.point_normals().expect("the sample carries normals");
        for (i, mp) in surface.iter().enumerate() {
            assert_eq!(cloud.points()[i], mp.sp.point);
            assert_eq!(normals[i], mp.sp.normal);
        }
        assert!(cloud.point_colors().is_none());
        assert!(cloud.point_stdev().is_none());
    }

    #[test]
    fn dense_sampling_skips_degenerate_faces() {
        // A square with one of its two faces collapsed onto a line.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 1, 3]];

        let sampled = compute_dense_surface_points(&points, &faces, 0.25, None).unwrap();
        assert!(!sampled.is_empty());
        assert!(sampled.iter().all(|mp| mp.face_index == 0));

        let sampled = compute_uniform_surface_points(&points, &faces, 50);
        assert_eq!(sampled.len(), 50);
        assert!(sampled.iter().all(|mp| mp.face_index == 0));

        let only_degenerate = compute_uniform_surface_points(&points, &faces[1..], 10);
        assert!(only_degenerate.is_empty());
    }

    #[test]
    fn uniform_samples_lie_on_the_surface_with_the_face_normal() {
        let mesh = Mesh3::create_box(10.0, 6.0, 4.0, false);
        let cloud = mesh.sample_uniform(500);
        assert_eq!(cloud.point_count(), 500);

        let normals = cloud.point_normals().unwrap();
        for (p, n) in cloud.points().iter().zip(normals) {
            let closest = mesh.surface_closest_to(p);
            assert!(dist(p, &closest.sp.point) < 1e-9);
            assert!((closest.sp.normal.dot(n) - 1.0).abs() < 1e-9);
        }
    }
}
