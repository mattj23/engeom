//! This module provides functionality for voxel downsampling of points

use crate::common::IndexMask;
use crate::common::PCoords;
use crate::na::SVector;
use parry3d_f64::utils::hashmap::HashMap;

/// Perform a voxel downsample of the given set of points, returning a mask that indicates which
/// points are retained after the operation.
///
/// A voxel downsampling is a relatively fast operation that converts the point coordinates to
/// integers representing a grid the size of the voxel spacing, and then retains only the first
/// point encountered in each voxel.
///
/// While this doesn't guarantee anything resembling a uniform distribution, it does serve as a
/// fast and effective pre-filtering step for more complicated downsampling methods, such as the
/// Poisson disk sampling method which has to construct a spatial tree.
///
/// # Arguments
///
/// * `points`: A slice of points implementing the `PCoords` trait for the specified dimension `D`.
/// * `voxel_size`: The size of the voxel grid to use for downsampling. This value should be
///   positive and non-zero.
///
/// returns: IndexMask
pub fn voxel_downsample<const D: usize>(points: &[impl PCoords<D>], voxel_size: f64) -> IndexMask {
    let mut voxel_map = HashMap::new();
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Aabb2;
    use crate::{Point2, Result};

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

        let mask = voxel_downsample(&points, 0.1);
        let thinned = mask.clone_indices_of(&points)?;

        for i in 0..100 {
            for j in 0..100 {
                let x = i as f64 * s1;
                let y = j as f64 * s1;
                let aabb = Aabb2::new(Point2::new(x, y), Point2::new(x + s1, y + s1));

                let count = thinned
                    .iter()
                    .filter(|p| aabb.contains_local_point(*p))
                    .count();

                assert_eq!(count, 1);
            }
        }

        Ok(())
    }
}
