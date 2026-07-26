//! This module holds the implementations which are shared between the different mesh
//! representations.
//!
//! A triangle mesh shows up in this library in more than one form: `MeshData3` as a plain editable
//! container, the accelerated representation built on `parry3d` for spatial queries, and the
//! half-edge structure for topological editing. Many useful operations depend on nothing but the
//! point and face buffers, and every one of those should exist exactly once, here, rather than being
//! reimplemented per representation and left to drift apart.
//!
//! # Why free functions
//!
//! Everything in this module takes raw buffers rather than a mesh type:
//!
//! ```ignore
//! pub fn compute_point_normals(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<UnitVec3>>
//! ```
//!
//! That makes the algorithms independent of any particular container, testable without constructing
//! one, and usable directly by a caller that already has buffers in hand. A trait which lets these
//! be called as methods on the mesh types is expected eventually, and it will be a thin delegation
//! layer over these functions rather than a second home for the logic.

pub mod normals;
pub mod offset;

pub use normals::compute_point_normals;
pub use offset::{OffsetOpts, compute_normal_displaced_points, compute_face_offset_points};
