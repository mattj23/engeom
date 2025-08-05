//! This is a very simple averaging struct that will compute the average or weighted average of a
//! set of values using two scalars: a sum and a count.

use num_traits::{One, Zero};
use std::ops::{Add, Div};

/// This is a very simple, low memory struct to compute the average or weighted average of a set of
/// values using a sum and a count. The only reason that you would use this instead of two
/// variables is when you're trying to average several different things at once, and you want to
/// halve the number of variables you need to keep track of.
pub struct Averager<T>
where
    T: Add + Div + Zero + One + Copy,
{
    sum: T,
    count: T,
}

impl<T> Default for Averager<T>
where
    T: Add<Output = T> + Div<Output = T> + Zero + One + Copy,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Averager<T>
where
    T: Add<Output = T> + Div<Output = T> + Zero + One + Copy,
{
    pub fn new() -> Self {
        Self {
            sum: T::zero(),
            count: T::zero(),
        }
    }

    pub fn add(&mut self, value: T) {
        self.sum = self.sum + value;
        self.count = self.count + T::one();
    }

    pub fn add_weight(&mut self, value: T, weight: T) {
        self.sum = self.sum + (value * weight);
        self.count = self.count + weight;
    }

    pub fn average(&self) -> Option<T> {
        if self.count.is_zero() {
            return None;
        }
        Some(self.sum / self.count)
    }
}
