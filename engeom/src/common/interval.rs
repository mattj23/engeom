//! This module contains an abstraction which represents an interval on a continuous scalar domain,
//! such as the interval [0, 1] on the real number line.  Intervals allow for the testing of
//! intersections between ranges of values, and provide tools for gathering, merging, and separating
//! regions.

mod merge_domain;
mod infinite_domain;
mod angle_domain;

pub use merge_domain::*;
pub use infinite_domain::Interval;
pub use angle_domain::AngleInterval;

