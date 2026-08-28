//! Choosing the static body of a multi-body alignment.
//!
//! The multi-body adjustments hold one body fixed, and when the caller does not choose it these
//! helpers rank the candidates: the body most broadly referenced by the others makes the best
//! static reference. The dimension-specific part, counting how many sample points one body
//! contributes against another, comes in as a closure.

use crate::na::DMatrix;
use rayon::prelude::*;

/// The pairwise correspondence counts of a set of bodies, as a symmetric matrix with a zero
/// diagonal.
///
/// Each `i, j` entry is `pair_count(i, j)`: the number of sample points in body `j` which are a
/// good match for body `i`. Only the upper triangle is evaluated (in parallel), and the count is
/// mirrored, so `pair_count` must be symmetric in meaning even though it is called with `i < j`.
pub(crate) fn correspondence_matrix<F>(count: usize, pair_count: F) -> DMatrix<f64>
where
    F: Fn(usize, usize) -> f64 + Sync,
{
    let mut matrix = DMatrix::<f64>::zeros(count, count);

    let mut work_list = Vec::new();
    for i in 0..count {
        for j in (i + 1)..count {
            work_list.push((i, j));
        }
    }

    let collected = work_list
        .par_iter()
        .map(|&(i, j)| (i, j, pair_count(i, j)))
        .collect::<Vec<_>>();

    for (i, j, count) in collected {
        matrix[(i, j)] = count;
        matrix[(j, i)] = count;
    }

    matrix
}

/// Orders the bodies of a correspondence matrix by how broadly the others reference them,
/// most-referenced first. The head of the returned order is the best candidate for the static
/// body.
///
/// The row with the highest column sum has the most points referencing it, but a raw count would
/// let two heavily overlapping bodies inflate each other without either being a good static
/// reference. Instead each cell is scaled to the range 0..1 and square-rooted, which grants
/// diminishing returns to a large count from any single body and favors bodies referenced by many
/// others.
pub(crate) fn reference_priority(matrix: DMatrix<f64>) -> Vec<usize> {
    let max = matrix.max();
    let mut corr = matrix;
    if max > 0.0 {
        corr /= max;
    }
    corr.apply(|x| *x = x.sqrt());

    let mut pairs = corr
        .column_sum()
        .iter()
        .enumerate()
        .map(|(i, x)| (i, *x))
        .collect::<Vec<_>>();
    pairs.sort_by(|a, b| b.1.total_cmp(&a.1));
    pairs.iter().map(|(i, _)| *i).collect()
}
