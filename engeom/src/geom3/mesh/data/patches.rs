//! Connected-patch labeling and filtering on the unaccelerated container.
//!
//! These methods mirror the methods with the same names on `Mesh3`. They are available on
//! `MeshData3` because patch operations only walk the face list and do not need a bounding volume
//! hierarchy. Building a `Mesh3` only to filter patches creates a hierarchy over geometry that is
//! about to be discarded, followed by another hierarchy over the retained geometry.

use crate::common::IndexMask;
use crate::geom3::mesh::nav_structure::MeshNav;
use crate::geom3::mesh::patches::{PatchFilter, PatchLabels};
use crate::{MeshData3, Result};

impl MeshData3 {
    /// Build a navigation structure for traversing this mesh by edge and face.
    ///
    /// Retain this structure when performing multiple structural queries. Each convenience method
    /// below creates and discards its own structure.
    pub fn compute_nav(&self) -> MeshNav<'_> {
        MeshNav::new(self)
    }

    /// Label each face with its connected patch.
    ///
    /// See [`Mesh3::compute_patch_labels`](crate::Mesh3::compute_patch_labels), which this mirrors.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask that restricts participation to selected faces, as if the mesh
    ///   had been pruned to those faces
    ///
    /// returns: `Result<PatchLabels>`
    pub fn compute_patch_labels(&self, mask: Option<&IndexMask>) -> Result<PatchLabels> {
        self.compute_nav().patch_labels(mask)
    }

    /// Build a face mask that selects the connected patches accepted by a filter.
    ///
    /// Use this method to inspect which faces [`MeshData3::remove_small_patches`] would discard or
    /// to combine the selection with other criteria before extraction.
    ///
    /// # Arguments
    ///
    /// * `filter`: the criteria that determine which patches to keep
    ///
    /// returns: `Result<IndexMask>` over this mesh's faces
    pub fn patch_mask(&self, filter: &PatchFilter) -> Result<IndexMask> {
        self.compute_nav().patch_mask(filter, None)
    }

    /// Discard connected patches rejected by a filter and return the retained mesh data.
    ///
    /// If the filter accepts every patch, this method clones the mesh without rebuilding it. All
    /// index mappings remain unchanged, and the no-op case is inexpensive.
    ///
    /// # Arguments
    ///
    /// * `filter`: the criteria that determine which patches to keep
    ///
    /// returns: `Result<MeshData3>`, failing if the filter would discard every face
    pub fn remove_small_patches(&self, filter: &PatchFilter) -> Result<Self> {
        let mask = self.patch_mask(filter)?;
        let kept = mask.count_true();

        if kept == 0 {
            return Err(
                "Every patch was discarded by the filter, which would leave an empty mesh".into(),
            );
        }

        if kept == self.faces().len() {
            return Ok(self.clone());
        }

        self.extract_subset_faces(&mask)
    }
}
