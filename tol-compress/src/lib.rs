//! Tolerance-bounded compression of spatial coordinate data.
//!
//! This crate stores spatial coordinates using a practical tolerance-compressed method, which
//! takes advantage of two specific features of real world measurement data:
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

pub mod bits;
pub mod bounds;
pub mod container;
pub mod error;
pub mod indices;
pub mod mesh;
pub mod metadata;
pub mod points;
pub mod polyline;
pub mod quantize;
mod raw;

/// Deliberately awkward test geometry, behind the non-default `corpus` feature.
#[cfg(any(test, feature = "corpus"))]
pub mod corpus;

#[cfg(any(test, feature = "corpus"))]
pub mod testgen;

pub use bounds::Bounds;
pub use container::{Header, Kind, probe};
pub use error::{Error, Result};
pub use indices::{bits_for_count, read_indices, write_indices};
pub use mesh::Mesh3;
pub use metadata::{Metadata, Value};
pub use points::{read_points, write_points};
pub use polyline::{Polyline, Polyline2, Polyline3};
pub use quantize::{Quantizer, bits_for_tol};
