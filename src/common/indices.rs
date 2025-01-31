//! This module should have tools for working with indices

pub fn index_vec(indices: Option<&[usize]>, len: usize) -> Vec<usize> {
    if let Some(items) = indices {
        items.to_vec()
    } else {
        (0..len).collect::<Vec<_>>()
    }
}

/// Identify chains of indices. Starting from a set of index pairs, assemble consecutive pairs into
/// chains of contiguous indices while preserving the order of each pair.  The input is a slice of
/// u32 pairs, and the output is a vector of vectors of u32. Each inner vector represents a chain
/// of contiguous indices.
///
/// # Arguments
///
/// * `indices`: original pairs of indices
///
/// returns: Vec<Vec<u32, Global>, Global>
///
/// # Examples
///
/// ```
///
/// ```
pub fn chained_indices(indices: &[[u32; 2]]) -> Vec<Vec<u32>> {
    let mut pairs = index_vec(None, indices.len());

    let mut chains: Vec<Vec<u32>> = Vec::new();
    let mut working = Vec::new();
    let mut forward = true;

    while !pairs.is_empty() {
        // If working is empty, start a new chain with the first pair
        if working.is_empty() {
            let i = pairs.pop().unwrap();
            working.push(indices[i][0]);
            working.push(indices[i][1]);
            forward = true;
        }

        if forward {
            let last = *working.last().unwrap();
            if let Some((k, i)) = chain_candidates(&pairs, indices, last, true) {
                working.push(indices[i][1]);
                pairs.swap_remove(k);
            } else {
                forward = false;
            }
        } else {
            let first = *working.first().unwrap();
            if let Some((k, i)) = chain_candidates(&pairs, indices, first, false) {
                working.insert(0, indices[i][0]);
                pairs.swap_remove(k);
            } else {
                chains.push(working.clone());
                working.clear();
            }
        }
    }

    if !working.is_empty() {
        chains.push(working);
    }

    chains
}

fn chain_candidates(
    pairs: &[usize],
    indices: &[[u32; 2]],
    last: u32,
    forward: bool,
) -> Option<(usize, usize)> {
    let mut candidates = Vec::new();
    let j = if forward { 0 } else { 1 };

    for (k, &i) in pairs.iter().enumerate() {
        if indices[i][j] == last {
            candidates.push((k, i));
        }
    }

    if candidates.len() == 1 {
        Some(candidates[0])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chain_candidates_forward() {
        let indices = [[0, 1], [1, 2], [2, 3]];
        let pairs = [1, 2];

        let result = chain_candidates(&pairs, &indices, 1, true);

        assert_eq!(result, Some((0, 1)));
    }

    #[test]
    fn chain_candidates_backwards() {
        let indices = [[0, 1], [1, 2], [2, 3]];
        let pairs = [0, 1];

        let result = chain_candidates(&pairs, &indices, 2, false);

        assert_eq!(result, Some((1, 1)));
    }

    #[test]
    fn test_chained_indices() {
        let indices = [[0, 1], [1, 2], [2, 3], [4, 5], [5, 6], [7, 8]];

        let mut chains = chained_indices(&indices);
        chains.sort_by(|a, b| a[0].cmp(&b[0]));

        assert_eq!(chains.len(), 3);
        assert_eq!(chains[0], vec![0, 1, 2, 3]);
        assert_eq!(chains[1], vec![4, 5, 6]);
        assert_eq!(chains[2], vec![7, 8]);
    }
}
