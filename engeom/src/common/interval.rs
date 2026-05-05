//! This module contains an abstraction which represents an interval on a continuous scalar domain,
//! such as the interval [0, 1] on the real number line.  Intervals allow for the testing of
//! intersections between ranges of values, and provide tools for gathering, merging, and separating
//! regions.

mod angle_domain;
mod infinite_domain;
mod merge;
mod merge_domain;

pub use angle_domain::AngleInterval;
pub use infinite_domain::Interval;
pub use merge::MergeDomain;
pub use merge_domain::*;

pub type IntervalMerge = MergeDomain<Interval>;

pub trait IntervalOps: Copy {
    fn min(&self) -> f64;
    fn max(&self) -> f64;
    fn contains_value(&self, x: f64) -> bool;
    fn contains_other(&self, other: Self) -> bool;
    fn extent(&self) -> f64;
    fn overlaps(&self, other: Self) -> bool;
    fn intersection(&self, other: Self) -> (Option<Self>, Option<Self>);
    fn clamp_value(&self, x: f64) -> f64;
    fn center(&self) -> f64;
    fn is_empty(&self) -> bool;
    fn new_containing(&self, other: &Self) -> Self;
    fn offset(&self, x: f64) -> Self;

    fn expanded(&self, half_width: f64) -> Self;
    fn dilated(&self, half_width: f64) -> Self;

    fn new_full() -> Self;
    fn wraps() -> bool;
}
