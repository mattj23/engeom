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
        IndexMask { mask }
    }

    /// Get the index values stored in the mask as a vector of usize.
    pub fn to_indices(&self) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.mask.len() / 16);

        for (i, bit) in self.mask.iter().enumerate() {
            if *bit {
                indices.push(i);
            }
        }

        indices
    }

    /// Modify the mask so that all values are the opposite of their current value.
    pub fn flip(&mut self) {
        for mut u in self.mask.as_raw_mut_slice() {
            let flipped = !*u;
            *u = flipped;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(test)]
    mod stress_tests {
        use super::*;
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        #[test]
        fn stress_set_and_get() {
            let len = 10_000;
            let mut rng = StdRng::seed_from_u64(42);
            let mut mask = IndexMask::new(len, false);
            let mut vec = vec![false; len];

            for _ in 0..100_000 {
                let idx = rng.random_range(0..len);
                let val = rng.random_bool(0.5);
                mask.set(idx, val);
                vec[idx] = val;
            }

            for i in 0..len {
                assert_eq!(mask.get(i), vec[i], "Mismatch at index {}", i);
            }
        }

        #[test]
        fn stress_flip() {
            let len = 5_000;
            let mut rng = StdRng::seed_from_u64(123);
            let mut mask = IndexMask::new(len, false);
            let mut vec = vec![false; len];

            for _ in 0..10 {
                for i in 0..len {
                    let val = rng.random_bool(0.5);
                    mask.set(i, val);
                    vec[i] = val;
                }
                mask.flip();
                for i in 0..len {
                    vec[i] = !vec[i];
                }
                for i in 0..len {
                    assert_eq!(mask.get(i), vec[i], "Mismatch after flip at index {}", i);
                }
            }
        }

        #[test]
        fn stress_to_indices() {
            let len = 2_000;
            let mut rng = StdRng::seed_from_u64(999);
            let mut mask = IndexMask::new(len, false);
            let mut vec = vec![false; len];

            for i in 0..len {
                let val = rng.random_bool(0.3);
                mask.set(i, val);
                vec[i] = val;
            }

            let mask_indices = mask.to_indices();
            let vec_indices: Vec<usize> = vec
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v { Some(i) } else { None })
                .collect();

            assert_eq!(mask_indices, vec_indices, "Indices mismatch");
        }
    }
}
