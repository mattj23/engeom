//! This module should have tools for working with indices

pub fn index_vec(indices: Option<&[usize]>, len: usize) -> Vec<usize> {
    if let Some(items) = indices {
        items.to_vec()
    } else {
        (0..len).collect::<Vec<_>>()
    }
}
