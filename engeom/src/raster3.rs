//! This module contains tools for working with 3D voxel grids.

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

        while !to_visit.is_empty() {
            let current = to_visit.pop().unwrap();
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
    let result = indices.iter().next().unwrap().clone();
    indices.remove(&result);
    result
}
