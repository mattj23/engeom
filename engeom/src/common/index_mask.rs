//! This module contains tools for working with indices as a mask of boolean values (may
//! eventually be implemented with bitvectors depending on real world performance).

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

    pub fn count_true(&self) -> usize {
        self.iter_true().count()
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

    // ==========================================================================================
    // NOT Operations
    // ==========================================================================================

    /// Modify the mask so that all values are the opposite of their current value.
    pub fn not_mut(&mut self) {
        for u in self.mask.as_raw_mut_slice() {
            *u = !*u;
        }
    }

    /// Return a new mask that is the opposite/inversion of the current mask.
    pub fn not(&self) -> Self {
        let mut new_mask = self.clone();
        new_mask.not_mut();
        new_mask
    }

    // ==========================================================================================
    // OR (Union) Operations
    // ==========================================================================================

    /// Creates and returns a new mask that is the result of a bitwise OR operation between this
    /// mask and another mask. This is the equivalent of a set union operation, as it will
    /// combine the indices that are true in either mask.
    ///
    /// This will return an error if the two masks are not of the same length.
    ///
    /// # Arguments
    ///
    /// * `other`: The other `IndexMask` to OR with this one.
    ///
    /// returns: Result<IndexMask, Box<dyn Error, Global>>
    pub fn or(&self, other: &IndexMask) -> Result<Self> {
        let mut new_mask = self.clone();
        new_mask.or_mut(other)?;
        Ok(new_mask)
    }

    /// Performs a bitwise OR operation on the mask with another mask, modifying the current mask.
    /// This is the equivalent of a set union operation, as it will combine the indices that are
    /// true in either mask.
    ///
    /// This will return an error if the two masks are not of the same length.
    ///
    /// # Arguments
    ///
    /// * `other`: The other `IndexMask` to OR with this one.
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    pub fn or_mut(&mut self, other: &IndexMask) -> Result<()> {
        if self.mask.len() != other.mask.len() {
            return Err("Masks must be of the same length".into());
        }

        self.or_mut_unchecked(other);
        Ok(())
    }

    /// Performs a bitwise OR operation on the mask with another mask, modifying the current mask.
    /// This does not check for length, so it zips the two masks directly. It will not panic if the
    /// lengths are different, but the result will not be what's expected. ONLY USE THIS IF YOU ARE
    /// SURE THAT THE MASKS ARE OF THE SAME LENGTH.
    ///
    /// Because the index mask internally uses a bit vector backed by `usize`, don't expect this
    /// method to combine only up to the length of the shorter mask. If this mask is longer than
    /// `other`, the extra bits at the end of the last 64 bit chunk will contaminate their
    /// counterparts in this mask.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: ()
    pub fn or_mut_unchecked(&mut self, other: &IndexMask) {
        let self_mask = self.mask.as_raw_mut_slice();
        let other_mask = other.mask.as_raw_slice();

        for (a, b) in self_mask.iter_mut().zip(other_mask.iter()) {
            *a |= *b;
        }
    }

    // ==========================================================================================
    // AND (Intersection) Operations
    // ==========================================================================================

    /// Creates and returns a new mask that is the result of a bitwise AND operation between this
    /// mask and another mask. This is the equivalent of a set intersection operation, as it will
    /// only keep the indices that are true in both masks.
    ///
    /// Will return an error if the two masks are not of the same length.
    ///
    /// # Arguments
    ///
    /// * `other`: the other `IndexMask` to AND with this one.
    ///
    /// returns: Result<IndexMask, Box<dyn Error, Global>>
    pub fn and(&self, other: &IndexMask) -> Result<Self> {
        let mut new_mask = self.clone();
        new_mask.and_mut(other)?;
        Ok(new_mask)
    }
    pub fn and_mut(&mut self, other: &IndexMask) -> Result<()> {
        if self.mask.len() != other.mask.len() {
            return Err("Masks must be of the same length".into());
        }

        self.and_mut_unchecked(other);
        Ok(())
    }

    /// Performs a bitwise AND operation on the mask with another mask, modifying the current mask.
    /// This does not check for length, so it zips the two masks directly. It will not panic if the
    /// lengths are different, but the result will not be what's expected. ONLY USE THIS IF YOU ARE
    /// SURE THAT THE MASKS ARE OF THE SAME LENGTH.
    ///
    /// Because the index mask internally uses a bit vector backed by `usize`, don't expect this
    /// method to combine only up to the length of the shorter mask. If this mask is longer than
    /// `other`, the extra bits at the end of the last 64 bit chunk will contaminate their
    /// counterparts in this mask.
    ///
    /// # Arguments
    ///
    /// * `other`: a reference to another `IndexMask` to AND with this one.
    ///
    /// returns: ()
    pub fn and_mut_unchecked(&mut self, other: &IndexMask) {
        let self_mask = self.mask.as_raw_mut_slice();
        let other_mask = other.mask.as_raw_slice();

        for (a, b) in self_mask.iter_mut().zip(other_mask.iter()) {
            *a &= *b;
        }
    }

    // ==========================================================================================
    // AND NOT (Difference) Operations
    // ==========================================================================================

    /// Creates and returns a new mask that is the result of a bitwise AND NOT operation between
    /// this mask and another mask. This is the equivalent of a set difference operation, as it will
    /// remove the elements of this mask that are also true in the other mask.
    ///
    /// This is more efficient than using `this.and(other.not())` because it avoids creating a new
    /// mask for the NOT operation, instead doing the AND NOT operation directly.
    ///
    /// This will return an error if the two masks are not of the same length.
    ///
    /// # Arguments
    ///
    /// * `other`: the other `IndexMask` to AND NOT with this one.
    ///
    /// returns: Result<IndexMask, Box<dyn Error, Global>>
    pub fn and_not(&self, other: &IndexMask) -> Result<Self> {
        let mut new_mask = self.clone();
        new_mask.and_not_mut(other)?;
        Ok(new_mask)
    }

    /// Performs a bitwise AND NOT operation on this mask using another mask, modifying the current
    /// mask in place. This is the equivalent of a set difference operation, as it will
    /// remove the elements of this mask that are also true in the other mask.
    ///
    /// This is more efficient than using `this.and(other.not())` because it avoids creating a new
    /// mask for the NOT operation, instead doing the AND NOT operation directly.
    ///
    /// This will return an error if the two masks are not of the same length.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    pub fn and_not_mut(&mut self, other: &IndexMask) -> Result<()> {
        if self.mask.len() != other.mask.len() {
            return Err("Masks must be of the same length".into());
        }

        let self_mask = self.mask.as_raw_mut_slice();
        let other_mask = other.mask.as_raw_slice();

        for (a, b) in self_mask.iter_mut().zip(other_mask.iter()) {
            *a &= !*b;
        }
        Ok(())
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
