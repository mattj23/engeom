//! Serialization that uses the practical tolerance-compressed method of storing spatial
//! coordinates: the caller supplies the largest round-trip error they are willing to accept, and
//! the storage layer finds the smallest representation that guarantees it.
//!
//! The method and the on-disk format live in the standalone [`tol_compress`] crate, which has no
//! dependencies and knows nothing about engeom. What is here is the adapter layer that converts
//! between engeom's geometry types and that crate's containers, plus the two engeom-specific
//! decisions that adapter has to make:
//!
//! - A [`Curve2`](crate::Curve2) or [`Curve3`](crate::Curve3) carries a chordal reconstruction
//!   tolerance, which is not a property of the stored points. It travels as item metadata under
//!   [`curve::CHORD_TOL`], a key the format stores and never interprets.
//! - The format stores geometry only, so writing a mesh carrying attributes is an error naming
//!   them rather than a silent loss. See the [`mesh`] module.

pub mod curve;
pub mod mesh;
