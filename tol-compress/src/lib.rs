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
//!
//! # Where to start
//!
//! Most callers want one of the three container modules, which read and write whole files:
//!
//! | module | holds | extension |
//! |---|---|---|
//! | [`mesh`] | triangle meshes in 3D | `.tcmesh` |
//! | [`polyline`] | ordered polylines in 2D or 3D, open or closed | `.tccurve2`, `.tccurve3` |
//! | [`cloud`] | unordered point sets in 2D or 3D | `.tccloud2`, `.tccloud3` |
//!
//! ```no_run
//! use std::path::Path;
//! use tol_compress::{Mesh3, mesh};
//!
//! let mesh = Mesh3::new(
//!     vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
//!     vec![[0, 1, 2]],
//! );
//!
//! // Every vertex is guaranteed back within one micron of where it started.
//! mesh::write_one_file(Path::new("part.tcmesh"), &mesh, 1e-6)?;
//! let back = mesh::read_one_file(Path::new("part.tcmesh"))?;
//! # Ok::<(), tol_compress::Error>(())
//! ```
//!
//! Every container is an ordered **collection**, so a file may hold one item or many. The
//! `write_one_*` and `read_one_*` functions are conveniences over that rather than a separate
//! format, and reading one item from a file holding several is an error rather than a silent
//! truncation. [`probe`] identifies a file without decoding its geometry.
//!
//! Items carry an optional name and an optional [`Metadata`] map, which is stored and never
//! interpreted, so a caller can record whatever their domain needs without this crate learning
//! anything about it.
//!
//! # Underneath
//!
//! The containers are built from block encoders that are useful on their own if you are storing
//! geometry inside some other format:
//!
//! - [`points`] encodes coordinates in any number of dimensions against a tolerance, and
//!   [`quantize`] and [`bounds`] are where the tolerance guarantee is actually decided.
//! - [`indices`] encodes simplex connectivity, exactly, at whichever of its two codings is smaller.
//! - [`reorder`] renumbers a mesh so its indices compress. [`mesh`] applies this by default; see
//!   its documentation if you hold per-vertex data outside the file.
//! - [`blocks`] and [`bits`] are the packing layer everything else sits on.
//!
//! # Errors and forward compatibility
//!
//! Everything fallible returns [`Error`], which implements [`std::error::Error`], so a caller using
//! a boxed-error alias absorbs it through `?` with no glue.
//!
//! Several fields in the format are reserved for features not yet implemented: compression,
//! per-item attributes, a point-block transform, and multiple partitions. A reader from this
//! version meeting a file that uses one of them **fails with a specific error rather than
//! misreading it**, and files written today stay readable as those features arrive.

#![deny(missing_docs)]
#![deny(rustdoc::broken_intra_doc_links)]

pub mod bits;
pub mod blocks;
pub mod bounds;
pub mod cloud;
pub mod container;
pub mod error;
pub mod indices;
pub mod mesh;
pub mod metadata;
pub mod points;
pub mod polyline;
pub mod quantize;
mod raw;
pub mod reorder;

/// Compiles and runs the README's examples as doctests, so the first thing a reader copies cannot
/// silently rot. Not part of the public API and not rendered anywhere.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

/// Deliberately awkward test geometry, behind the non-default `corpus` feature.
#[cfg(any(test, feature = "corpus"))]
pub mod corpus;

#[cfg(any(test, feature = "corpus"))]
pub mod testgen;

pub use blocks::BLOCK;
pub use bounds::Bounds;
pub use cloud::{Cloud, Cloud2, Cloud3};
pub use container::{Header, Kind, Named, find_by_name, probe};
pub use error::{Error, Result};
pub use indices::{bits_for_count, read_indices, write_indices};
pub use mesh::{Mesh3, VertexOrder, WriteOptions};
pub use metadata::{Metadata, Value};
pub use points::{read_points, write_points};
pub use polyline::{Polyline, Polyline2, Polyline3};
pub use quantize::{Quantizer, bits_for_tol};
pub use reorder::{Reordering, permute};
