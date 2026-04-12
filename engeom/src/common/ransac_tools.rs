//! This module is for collecting tools that are useful in RANSAC

use crate::Result;
use rand::{RngExt, rng};
use std::collections::HashSet;

/// Generate `group_count` groups of `N` indices from a list of `index_count` total possibilities.
/// The indices in a single group will be unique.
///
/// # Arguments
///
/// * `group_count`: the total number of groups to generate
/// * `index_count`: the total number of possible indices
///
/// returns: Result<Vec<[usize; N], Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn ransac_indices<const N: usize>(
    group_count: usize,
    index_count: usize,
) -> Result<Vec<[usize; N]>> {
    if index_count < N {
        return Err("Not enough indices".into());
    }

    // TODO: this method relies on there being a lot more indices than the size of each group.
    // TODO: Make a more reasonable implementation that adjusts based on the size

    let mut result = Vec::with_capacity(group_count);
    let mut working = HashSet::with_capacity(N);
    let mut rng = rng();

    for _ in 0..group_count {
        working.clear();
        let mut indices = [0; N];
        for i in 0..N {
            loop {
                let j = rng.random_range(0..index_count);
                if !working.contains(&j) {
                    working.insert(j);
                    indices[i] = j;
                    break;
                }
            }
        }
        result.push(indices);
    }

    Ok(result)
}
