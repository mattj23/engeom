//! This module contains common tools for serialization that uses the practical tolerance-compressed
//! method of storing spatial coordinates.
//!
//! This method attempts to take advantage of two specific features of practical measurement data:
//!
//! 1. Real world measurement systems have known precisions below which differences in position
//!    do not represent meaningful information, and users of measurement data have knowledge about
//!    where increasingly small differences in values cease to be relevant to their use case.
//! 2. Values produced by the measurement of physical objects rarely span more than a few orders
//!    of magnitude more than the smallest meaningful precision for the measurement system and/or
//!    the application.
//!
//! If a user can supply a tolerance that represents the largest acceptable round-trip accuracy of
//! any value in their dataset, the storage algorithm can find a representation scheme that uses
//! the smallest amount of bytes necessary to guarantee every position is recovered within that
//! tolerance.

mod common;
mod common3;
mod indices;
mod common2;

use crate::io::tol_compress::common::read_u32;

pub use indices::*;
pub use common2::*;
pub use common3::*;