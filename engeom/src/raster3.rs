//! This module contains tools for working with 3D voxel grids and the surfaces extracted from
//! them.

mod marching_cubes;
mod mc_tables;
mod sparse_grid;

pub use marching_cubes::{IsoSurfaceStats3, extract_isosurface};
pub use sparse_grid::{
    BLOCK_EDGE, BLOCK_VOXELS, Block3, BlockKey3, SdfVoxel, SparseGrid3, VoxelKey3, VoxelValue,
};

use std::collections::HashSet;

/// This function takes a set of coordinates in a 3D grid and returns a list of clusters of
/// connected voxel coordinates.
///
/// # Arguments
///
/// * `indices`:
///
/// returns: Vec<Vec<(u32, u32, u32), Global>, Global>
///
/// # Examples
///
/// ```
///
/// ```
pub fn clusters_from_sparse(mut indices: HashSet<(i32, i32, i32)>) -> Vec<Vec<(i32, i32, i32)>> {
    let mut results = Vec::new();

    while !indices.is_empty() {
        let mut working = Vec::new();
        let mut to_visit = Vec::new();

        to_visit.push(pop_index(&mut indices));

        while let Some(current) = to_visit.pop() {
            working.push(current);

            for x in -1..=1 {
                for y in -1..=1 {
                    for z in -1..=1 {
                        if x == 0 && y == 0 && z == 0 {
                            continue;
                        }

                        let neighbor = (current.0 + x, current.1 + y, current.2 + z);
                        if indices.remove(&neighbor) {
                            to_visit.push(neighbor);
                        }
                    }
                }
            }
        }

        results.push(working.into_iter().collect());
    }

    results
}

fn pop_index(indices: &mut HashSet<(i32, i32, i32)>) -> (i32, i32, i32) {
    let result = *indices.iter().next().unwrap();
    indices.remove(&result);
    result
}
