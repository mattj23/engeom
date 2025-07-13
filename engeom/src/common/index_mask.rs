//! This module contains tools for working with indices as a mask of boolean values (may
//! eventually be implemented with bitvectors depending on real world performance).

use std::ops::Not;
use crate::Result;
use bitvec::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
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

    pub fn try_from_indices(indices: &[usize], len: usize) -> Result<Self> {
        let mut mask = BitVec::repeat(false, len);
        for &index in indices {
            if index >= len {
                return Err(format!("Index {} is out of bounds for length {}", index, len).into());
            }
            mask.set(index, true);
        }
        Ok(IndexMask { mask })
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
    pub fn not_mut(&mut self) {
        for u in self.mask.as_raw_mut_slice() {
            *u = !*u;
        }
    }
    
    pub fn not(&self) -> Self {
        let mut new_mask = self.clone();
        new_mask.not_mut();
        new_mask
    }
    

    pub fn or_mut(&mut self, other: &IndexMask) -> Result<()> {
        if self.mask.len() != other.mask.len() {
            return Err("Masks must be of the same length".into());
        }

        let self_mask = self.mask.as_raw_mut_slice();
        let other_mask = other.mask.as_raw_slice();

        for (a, b) in self_mask.iter_mut().zip(other_mask.iter()) {
            *a |= *b;
        }

        Ok(())
    }
    
    pub fn or(&self, other: &IndexMask) -> Result<Self> {
        let mut new_mask = self.clone();
        new_mask.or_mut(other)?;
        Ok(new_mask)
    }

    pub fn and_mut(&mut self, other: &IndexMask) -> Result<()> {
        if self.mask.len() != other.mask.len() {
            return Err("Masks must be of the same length".into());
        }

        let self_mask = self.mask.as_raw_mut_slice();
        let other_mask = other.mask.as_raw_slice();

        for (a, b) in self_mask.iter_mut().zip(other_mask.iter()) {
            *a &= *b;
        }

        Ok(())
    }
    
    pub fn and(&self, other: &IndexMask) -> Result<Self> {
        let mut new_mask = self.clone();
        new_mask.and_mut(other)?;
        Ok(new_mask)
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

    pub fn iter_true(&self) -> MaskTrueIterator {
        MaskTrueIterator {
            mask: self,
            current: 0,
        }
    }

    pub fn clone_indices_of<T: Clone>(&self, items: &[T]) -> Result<Vec<T>> {
        if items.len() != self.len() {
            return Err(format!(
                "Items length {} does not match mask length {}",
                items.len(),
                self.len()
            )
            .into());
        }

        let mut result = Vec::with_capacity(self.to_indices().len());
        for index in self.iter_true() {
            result.push(items[index].clone());
        }
        Ok(result)
    }
}

pub struct MaskTrueIterator<'a> {
    mask: &'a IndexMask,
    current: usize,
}

impl<'a> Iterator for MaskTrueIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        while self.current < self.mask.len() {
            if self.mask.get(self.current) {
                let index = self.current;
                self.current += 1;
                return Some(index);
            }
            self.current += 1;
        }
        None
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
                mask.not_mut();
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

        #[test]
        fn stress_to_true_iter() {
            let len = 10_000;
            let mut rng = StdRng::seed_from_u64(345);
            let mut mask = IndexMask::new(len, false);
            let mut vec = vec![false; len];

            for i in 0..len {
                let val = rng.random_bool(0.3);
                mask.set(i, val);
                vec[i] = val;
            }

            let mut mask_indices = Vec::new();
            for index in mask.iter_true() {
                mask_indices.push(index);
            }

            let vec_indices: Vec<usize> = vec
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| if v { Some(i) } else { None })
                .collect();

            assert_eq!(mask_indices, vec_indices, "Indices mismatch");
        }

        #[test]
        fn stress_and() {
            let len = 100_000;
            let mut rng = StdRng::seed_from_u64(345);

            for _ in 0..100 {
                let mut mask0 = IndexMask::new(len, false);
                let mut mask1 = IndexMask::new(len, false);

                let mut vec0 = vec![false; len];
                let mut vec1 = vec![false; len];

                for i in 0..len {
                    let val0 = rng.random_bool(0.3);
                    vec0[i] = val0;
                    mask0.set(i, val0);

                    let val1 = rng.random_bool(0.3);
                    vec1[i] = val1;
                    mask1.set(i, val1);
                }

                let expected = vec0
                    .iter()
                    .zip(vec1.iter())
                    .map(|(&v0, &v1)| v0 && v1)
                    .collect::<Vec<bool>>();

                mask0.and_mut(&mask1).unwrap();

                for i in 0..len {
                    assert_eq!(
                        mask0.get(i),
                        expected[i],
                        "Mismatch after OR at index {}",
                        i
                    );
                }
            }
        }

        #[test]
        fn stress_or() {
            let len = 100_000;
            let mut rng = StdRng::seed_from_u64(345);

            for _ in 0..100 {
                let mut mask0 = IndexMask::new(len, false);
                let mut mask1 = IndexMask::new(len, false);

                let mut vec0 = vec![false; len];
                let mut vec1 = vec![false; len];

                for i in 0..len {
                    let val0 = rng.random_bool(0.3);
                    vec0[i] = val0;
                    mask0.set(i, val0);

                    let val1 = rng.random_bool(0.3);
                    vec1[i] = val1;
                    mask1.set(i, val1);
                }

                let expected = vec0
                    .iter()
                    .zip(vec1.iter())
                    .map(|(&v0, &v1)| v0 || v1)
                    .collect::<Vec<bool>>();

                mask0.or_mut(&mask1).unwrap();

                for i in 0..len {
                    assert_eq!(
                        mask0.get(i),
                        expected[i],
                        "Mismatch after OR at index {}",
                        i
                    );
                }
            }
        }
    }
}
