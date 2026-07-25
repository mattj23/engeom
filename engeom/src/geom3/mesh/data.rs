//! This module contains `MeshData3`, a plain container for triangle mesh data and its associated
//! per-element attributes.
//!
//! Unlike `Mesh`, which is built around a `parry3d` `TriMesh` and pays for a BVH build on every
//! construction, `MeshData3` holds nothing but the point and face buffers and whatever attributes
//! are attached to them. It is the type to reach for when editing, subsetting, or serializing mesh
//! data, and is converted into an accelerated representation only when spatial queries are needed.

mod attribute_set;
mod attributes;

pub use attribute_set::MeshAttrSet3;
pub use attributes::MeshAttr3;
