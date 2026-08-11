//! Mapping points onto a regular grid of voxels, and the two ways of using that grid to make a
//! large point set smaller.
//!
//! Both operations start the same way, by dividing each coordinate by the cell size and flooring it
//! to an integer key, so that two points share a cell exactly when their keys match. What they do
//! afterwards is different, and the difference matters more than the shared name suggests.
//!
//! [`compute_first_per_voxel_mask`] **thins**: it keeps one original point per occupied cell and
//! discards the rest, reporting the survivors as a mask. Nothing is averaged, so every surviving
//! point is a measurement that was actually taken, carrying its original noise.
//!
//! [`compute_voxel_groups`] **partitions**: it reports every point in every occupied cell, which is
//! what a caller needs when it intends to combine them. `PointCloud3::reduce_by_voxel` builds on it
//! to **reduce**, replacing each cell with a new averaged point that was never measured but is
//! correspondingly less noisy.
//!
//! They use different algorithms on purpose. Thinning only needs to answer "have I seen this voxel
//! before", so a hash set does it in one linear pass. Grouping needs the members of each voxel
//! together, and sorting by key gets that with two allocations rather than one per voxel, and gets
//! a deterministic output order at the same time, since hash iteration order is not stable.

use crate::Result;
use crate::common::IndexMask;
use crate::common::PCoords;
use crate::na::SVector;
use faer::prelude::default;
use parry3d_f64::utils::hashmap::HashMap;

/// Thin a set of points by keeping the first one encountered in each voxel, reported as a mask over
/// the input.
///
/// One linear pass and a hash set, which makes this the cheapest thing in this module and a good
/// pre-filter for anything more expensive. Poisson disk sampling uses it that way, to shrink the
/// point set before paying for a spatial tree.
///
/// # Which point survives, and why the name says so
///
/// The survivor is the **lowest-indexed** point in each cell, which is what makes the result
/// deterministic. It is worth knowing that this is not a neutral choice. Real scan data arrives in
/// acquisition order, so "lowest index" usually means "earliest in the scan", and the survivors are
/// therefore biased toward one side of each cell rather than scattered within it. Where that bias
/// matters, [`compute_voxel_groups`] and a centroid give an unbiased cell representative instead.
///
/// This also guarantees nothing about the minimum spacing of the result. Two survivors in adjacent
/// cells can be arbitrarily close to each other; the guarantee is one point per cell, not a minimum
/// distance. Poisson disk sampling is the tool for that use case.
///
/// You _can_ safely assume that if two neighboring voxels are occupied the _maximum_ distance
/// between the associated points will be `√(2v² + (2v)²)`, or roughly `2.5v`, where `v` is
/// the voxel edge size.  So if you have points that were originally more dense than than voxel edge
/// spacing, you know that the distance between neighboring points in the output mask will be no
/// smaller than the floating point epsilon and no larger than 2.5x the voxel edge size.
///
/// # Arguments
///
/// * `points`: A slice of points implementing the `PCoords` trait for the specified dimension `D`.
/// * `voxel_size`: The size of the voxel grid to use. This value should be positive and non-zero.
///
/// returns: IndexMask
pub fn compute_first_per_voxel_mask<const D: usize>(
    points: &[impl PCoords<D>],
    voxel_size: f64,
) -> IndexMask {
    let mut voxel_map = HashMap::with_hasher(default());
    let mut mask = IndexMask::new(points.len(), false);

    for (i, xyz) in points.iter().enumerate() {
        let mut key: SVector<i32, D> = SVector::zeros();
        for (d, &coord) in xyz.coords().iter().enumerate() {
            key[d] = (coord / voxel_size).floor() as i32;
        }

        if !voxel_map.contains_key(&key) {
            voxel_map.insert(key, i);
            mask.set(i, true);
        }
    }

    mask
}

/// The points of a set, partitioned by which voxel of a regular grid they fall into.
///
/// Groups are non-empty, and both the groups and the indices within each group are in ascending
/// order, so the same input always produces the same output.
pub struct VoxelGroups {
    /// Point indices, arranged so that the members of each voxel are contiguous.
    ///
    /// These are `u32` rather than `usize` because this is the largest thing the partition
    /// allocates, at one entry per input point, and a cloud large enough to exhaust a `u32` index
    /// space would be consuming well over a hundred gigabytes in its coordinates alone.
    /// `compute_voxel_groups` rejects such an input rather than truncating.
    indices: Vec<u32>,

    /// Where each group starts in `indices`, with a trailing entry holding `indices.len()` so that
    /// group `i` is always `indices[offsets[i]..offsets[i + 1]]`.
    offsets: Vec<u32>,
}

impl VoxelGroups {
    /// The number of occupied voxels. Empty voxels are not represented.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// The total number of points across all groups, which is the size of the input.
    pub fn point_count(&self) -> usize {
        self.indices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The point indices in group `i`, ascending.
    pub fn group(&self, i: usize) -> &[u32] {
        &self.indices[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }

    /// Every group in order, each as a slice of ascending point indices.
    pub fn iter(&self) -> impl Iterator<Item = &[u32]> {
        (0..self.len()).map(|i| self.group(i))
    }
}

/// Partition a set of points by which voxel of a regular grid of `voxel_size` they fall into.
///
/// Only occupied voxels are reported, so the result has one group per distinct voxel and every
/// group is non-empty. See the module documentation for why this is separate from
/// [`compute_first_per_voxel_mask`] rather than sharing its implementation: thinning only needs to
/// know whether a cell has been seen, while grouping needs its members gathered together.
///
/// # Arguments
///
/// * `points`: A slice of points implementing the `PCoords` trait for the specified dimension `D`.
/// * `voxel_size`: The edge length of the grid cells. Must be finite and greater than zero.
///
/// returns: `Result<VoxelGroups>`
///
/// # A note on very large coordinates
///
/// Voxel keys are `i32`. A coordinate more than about `2e9 * voxel_size` from the origin saturates
/// the cast and collides with everything else out that far. That is far outside anything this
/// library measures, and detecting it would cost a branch per coordinate, so it is documented
/// rather than checked.
pub fn compute_voxel_groups<const D: usize>(
    points: &[impl PCoords<D>],
    voxel_size: f64,
) -> Result<VoxelGroups> {
    if !voxel_size.is_finite() || voxel_size <= 0.0 {
        return Err(format!("Voxel size must be finite and positive, got {voxel_size}").into());
    }

    if points.len() > u32::MAX as usize {
        return Err(format!(
            "Voxel grouping indexes points with u32 and cannot take {} of them",
            points.len()
        )
        .into());
    }

    if points.is_empty() {
        return Ok(VoxelGroups {
            indices: Vec::new(),
            offsets: vec![0],
        });
    }

    // Key every point, then sort by (key, index). Including the index makes the order total, so the
    // result does not depend on the sort being stable.
    let mut keyed: Vec<([i32; D], u32)> = points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut key = [0i32; D];
            for (d, &coord) in p.coords().iter().enumerate() {
                key[d] = (coord / voxel_size).floor() as i32;
            }
            (key, i as u32)
        })
        .collect();

    keyed.sort_unstable();

    let mut indices: Vec<u32> = Vec::with_capacity(keyed.len());
    let mut offsets: Vec<u32> = vec![0];

    let mut current = keyed[0].0;
    for (key, i) in keyed {
        if key != current {
            offsets.push(indices.len() as u32);
            current = key;
        }
        indices.push(i);
    }
    offsets.push(indices.len() as u32);

    Ok(VoxelGroups { indices, offsets })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Aabb2;
    use crate::{Point2, Point3};

    #[test]
    fn downsample_grid() -> Result<()> {
        let mut points = Vec::new();
        let s0 = 0.01;
        let s1 = s0 * 10.0;

        for i in 0..1000 {
            for j in 0..1000 {
                points.push(Point2::new(
                    i as f64 * s0 + s0 / 2.0,
                    j as f64 * s0 + s0 / 2.0,
                ));
            }
        }

        let mask = compute_first_per_voxel_mask(&points, 0.1);
        let thinned = mask.clone_indices_of(&points)?;

        for i in 0..100 {
            for j in 0..100 {
                let x = i as f64 * s1;
                let y = j as f64 * s1;
                let aabb = Aabb2::new(Point2::new(x, y), Point2::new(x + s1, y + s1));

                let count = thinned
                    .iter()
                    .filter(|p| aabb.contains_local_point(p))
                    .count();

                assert_eq!(count, 1);
            }
        }

        Ok(())
    }

    // ===========================================================================================
    // compute_voxel_groups
    // ===========================================================================================

    #[test]
    fn voxel_groups_partition_every_point_exactly_once() {
        let points: Vec<Point3> = (0..500)
            .map(|i| {
                let t = i as f64 * 0.37;
                Point3::new(t.sin() * 10.0, t.cos() * 10.0, t * 0.1)
            })
            .collect();

        let groups = compute_voxel_groups(&points, 1.0).expect("grouping failed");

        let mut seen = vec![0usize; points.len()];
        for g in groups.iter() {
            assert!(!g.is_empty(), "a reported group was empty");
            for &i in g {
                seen[i as usize] += 1;
            }
        }

        assert!(
            seen.iter().all(|&c| c == 1),
            "every point must appear in exactly one group"
        );
    }

    #[test]
    fn voxel_groups_put_co_located_points_together() {
        // Three clusters, each well inside its own voxel, at a spacing which cannot merge them.
        let points = vec![
            Point3::new(0.1, 0.1, 0.1),
            Point3::new(0.2, 0.2, 0.2),
            Point3::new(5.1, 0.1, 0.1),
            Point3::new(5.2, 0.2, 0.2),
            Point3::new(5.3, 0.3, 0.3),
            Point3::new(10.1, 0.1, 0.1),
        ];

        let groups = compute_voxel_groups(&points, 1.0).expect("grouping failed");

        let sizes: Vec<usize> = groups.iter().map(|g| g.len()).collect();
        assert_eq!(sizes, vec![2, 3, 1]);

        assert_eq!(groups.group(0), &[0, 1]);
        assert_eq!(groups.group(1), &[2, 3, 4]);
        assert_eq!(groups.group(2), &[5]);
    }

    /// Groups and their contents are ascending regardless of input order, which is what makes a
    /// reduction built on this reproducible run to run.
    #[test]
    fn voxel_groups_are_ordered_deterministically() {
        let points = vec![
            Point3::new(5.2, 0.0, 0.0),
            Point3::new(0.1, 0.0, 0.0),
            Point3::new(5.1, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
        ];

        let groups = compute_voxel_groups(&points, 1.0).expect("grouping failed");

        assert_eq!(groups.len(), 2);
        assert_eq!(groups.group(0), &[1, 3]);
        assert_eq!(groups.group(1), &[0, 2]);
    }

    #[test]
    fn voxel_groups_handle_negative_coordinates() {
        let points = vec![
            Point3::new(-0.5, -0.5, -0.5),
            Point3::new(-1.5, -0.5, -0.5),
            Point3::new(-0.9, -0.9, -0.9),
        ];

        let groups = compute_voxel_groups(&points, 1.0).expect("grouping failed");

        assert_eq!(groups.len(), 2);
        assert_eq!(groups.group(0), &[1]);
        assert_eq!(groups.group(1), &[0, 2]);
    }

    #[test]
    fn voxel_groups_reject_a_nonsense_size() {
        let points = vec![Point3::origin()];
        assert!(compute_voxel_groups(&points, 0.0).is_err());
        assert!(compute_voxel_groups(&points, -1.0).is_err());
        assert!(compute_voxel_groups(&points, f64::NAN).is_err());
        assert!(compute_voxel_groups(&points, f64::INFINITY).is_err());
    }

    #[test]
    fn voxel_groups_of_nothing_is_nothing() {
        let points: Vec<Point3> = Vec::new();
        let groups = compute_voxel_groups(&points, 1.0).expect("grouping failed");
        assert_eq!(groups.len(), 0);
        assert!(groups.is_empty());
    }

    /// The first index of each group is the lowest, which is the same point
    /// `compute_first_per_voxel_mask` keeps. The two are independent implementations, one hashing
    /// and one sorting, so this pins them to the same answer.
    #[test]
    fn voxel_groups_agree_with_the_thinning_mask_on_which_point_comes_first() {
        let points: Vec<Point3> = (0..300)
            .map(|i| {
                let t = i as f64 * 0.21;
                Point3::new(t.sin() * 4.0, t.cos() * 4.0, (t * 0.5).sin() * 4.0)
            })
            .collect();

        let groups = compute_voxel_groups(&points, 0.75).expect("grouping failed");
        let mask = compute_first_per_voxel_mask(&points, 0.75);

        let mut from_groups: Vec<usize> = groups.iter().map(|g| g[0] as usize).collect();
        from_groups.sort_unstable();

        assert_eq!(from_groups, mask.to_indices());
    }
}
