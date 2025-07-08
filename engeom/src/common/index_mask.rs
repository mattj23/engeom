//! This module contains tools for working with indices as a mask of boolean values (may
//! eventually be implemented with bitvectors depending on real world performance).

use crate::Result;

#[derive(Clone)]
pub struct IndexMask {
    mask: Vec<bool>,
}

impl IndexMask {

    /// Create a new IndexMask with the specified length and initial value.
    ///
    /// # Arguments
    ///
    /// * `len`: the length of the mask
    /// * `value`: the initial value for each index in the mask
    ///
    /// returns: IndexMask
    pub fn new(len: usize, value: bool) -> Self {
        IndexMask {
            mask: vec![value; len],
        }
    }

    /// Get the index values stored in the mask as a vector of usize.
    pub fn to_indices(&self) -> Vec<usize> {
        self.mask
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v { Some(i) } else { None })
            .collect()
    }

    /// Modify the mask so that all values are the opposite of their current value.
    pub fn flip(&mut self) {
        for value in &mut self.mask {
            *value = !*value;
        }
    }

    /// Set the value at the specified index
    pub fn set(&mut self, index: usize, value: bool) {
        if index < self.mask.len() {
            self.mask[index] = value;
        }
    }

    /// Get the value at the specified index.
    pub fn get(&self, index: usize) -> bool {
        self.mask[index]
    }

    /// Get the length of the mask.
    pub fn len(&self) -> usize {
        self.mask.len()
    }
}
