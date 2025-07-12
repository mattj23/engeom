//! This module contains tools for working with indices as a mask of boolean values (may
//! eventually be implemented with bitvectors depending on real world performance).
use bitvec::prelude::*;

#[derive(Clone)]
pub struct IndexMask {
    mask: BitVec,
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
        let mask = BitVec::repeat(value, len);
        IndexMask {
            mask
        }
    }

    /// Get the index values stored in the mask as a vector of usize.
    pub fn to_indices(&self) -> Vec<usize> {
        let mut indices = Vec::new();
        for i in 0..self.len() {
            if self.get(i) {
                indices.push(i);
            }
        }

        indices
    }

    /// Modify the mask so that all values are the opposite of their current value.
    pub fn flip(&mut self) {
        for mut b in self.mask.iter_mut() {
            *b = !*b;
        }
    }

    /// Set the value at the specified index
    pub fn set(&mut self, index: usize, value: bool) {
        self.mask.set(index, value);
    }

    /// Get the value at the specified index.
    pub fn get(&self, index: usize) -> bool {
        self.mask[index]
    }

    #[allow(clippy::len_without_is_empty)]
    /// Get the length of the mask.
    pub fn len(&self) -> usize {
        self.mask.len()
    }

}
