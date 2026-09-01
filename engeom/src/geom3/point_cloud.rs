//! This module contains the 3D point cloud types.
//!
//! [`PointCloud3`] is the owning container: a buffer of points plus a validated `PointAttrSet3`.
//! It is cheap to construct, cheap to edit, and is what file I/O reads and writes.
//!
//! [`CloudIndex3`] is the accelerated view over one. It borrows a cloud and owns a k-d tree built
//! from it, and every spatial query lives on it rather than on the cloud.
//!
//! # Why the tree is borrowed rather than owned by the point cloud
//!
//! `Mesh3` keeps its bounding volume hierarchy inside itself, but that's what `parry` chose and
//! it's complexity that `parry` manages internally.  When I've built accelerators or query
//! structures, I've done it with borrowed views that get built on demand.  That's how `MeshEdges`
//! and `MeshNav` work, and there's a lot to recommend of that approach too.
//!
//! I'm not 100% sold either way, but I decided to try doing the point cloud kd tree as a borrowed
//! view and see how it goes in use.  My main two thoughts were:
//!
//! - Building the view is not free compared to _not_ building it.  That's to say it's still an
//!   operation that takes time (a test, albiet on older hardware, showed building a tree over 2M
//!   points took 258ms), and that's a lot to accept when you aren't going to query anything.  This
//!   is an argument against having the point cloud own its tree, but it could also be settled by
//!   using the `MeshData3`/`Mesh3` split that the mesh features use.
//!
//! - A k-d tree can't be rigidly transformed; you have to rebuild it.  So if the cloud owned its
//!   own tree it would have to be rebuilt on evey mutation, which is what happens in `Mesh3`. The
//!   difference is that `Mesh3` is kind of just a wrapper over `parry::TriMesh`, and so the
//!   complexity isn't managed by me, it's managed by the authors of `parry`.  Having the borrowed
//!   view means that the borrow checker prevents mistakes from slipping through the cracks.
//!
//! That said, there are two reasons I don't like it.  First, it's different from how `Mesh3` and
//! `MeshData3` work, so it's going to be confusing to any user.  Second, I'm anticipating there's
//! going to be some friction from the borrow checker when working with it in the real world.  I'm
//! setting things up this way now so that, if in the course of putting some city miles on the
//! abstraction it turns out to have been a bad choice, I know to revert it and not try this again.

mod data;
mod normal_estimation;

use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::dist;
use crate::{KdTree3, Mesh3, Point3, Result};

pub use data::{PointCloud3, VOXEL_COHERENCE_ATTR, VOXEL_COUNT_ATTR};
pub use normal_estimation::NormalEstimates;

/// Finding which points of one entity are also covered by another, by testing whether the closest
/// point in each direction agrees.
pub trait PointCloudOverlap<TOther> {
    fn overlap_by_reciprocity(&self, other: &TOther, max_distance: f64) -> Vec<usize>;
}

/// A k-d tree built over a [`PointCloud3`], borrowed from it for the index's lifetime.
///
/// Build one with [`PointCloud3::compute_index`]. The borrow is what makes staleness impossible:
/// the cloud cannot be mutated while this exists, so the tree can never describe a point set that
/// has moved out from under it.
///
/// Throwing one away and building another is a reasonable thing to do when ownership gets awkward.
/// See the module documentation for the measurements which say so.
pub struct CloudIndex3<'a> {
    cloud: &'a PointCloud3,
    tree: TreeRef<'a>,
}

/// The tree an index queries, which it either owns or has been handed.
///
/// Owning is the normal case. Borrowing exists for a caller which keeps a tree alive across many
/// queries and rebuilds it on its own schedule; see [`PointCloud3::index_with_tree_unchecked`].
enum TreeRef<'a> {
    Owned(KdTree3),
    Borrowed(&'a KdTree3),
}

impl TreeRef<'_> {
    fn get(&self) -> &KdTree3 {
        match self {
            TreeRef::Owned(tree) => tree,
            TreeRef::Borrowed(tree) => tree,
        }
    }
}

impl<'a> CloudIndex3<'a> {
    /// Build an index over `cloud`. Prefer [`PointCloud3::compute_index`], which calls this.
    pub fn try_new(cloud: &'a PointCloud3) -> Result<Self> {
        let tree = KdTree3::try_new(cloud.points())?;
        Ok(Self {
            cloud,
            tree: TreeRef::Owned(tree),
        })
    }

    /// Pair a cloud with a tree the caller has already built over it, without checking that the two
    /// match. Reached through [`PointCloud3::index_with_tree_unchecked`], which is where the
    /// obligation is documented.
    pub(crate) fn with_tree_unchecked(cloud: &'a PointCloud3, tree: &'a KdTree3) -> Self {
        Self {
            cloud,
            tree: TreeRef::Borrowed(tree),
        }
    }

    /// The cloud this index was built over.
    pub fn cloud(&self) -> &'a PointCloud3 {
        self.cloud
    }

    /// The underlying k-d tree, for the query methods on [`KdTreeSearch`].
    pub fn tree(&self) -> &KdTree3 {
        self.tree.get()
    }

    /// The points of the underlying cloud. Query results index into this.
    pub fn points(&self) -> &[Point3] {
        self.cloud.points()
    }

    /// Estimate a normal at every point by fitting a plane to the neighbors within `radius`.
    ///
    /// Each estimated normal is flipped where needed to agree with the matching entry of
    /// `must_match`, since a plane fit gives an axis rather than a direction and cannot resolve the
    /// sign on its own. Points with fewer than three neighbors get `+Z` at zero confidence.
    ///
    /// The returned confidence is a measure of how plane-like the neighborhood was, so it is low on
    /// edges and corners and in sparse regions.
    pub fn estimate_normals(
        &self,
        must_match: &[crate::Vector3],
        radius: f64,
    ) -> Result<NormalEstimates> {
        if must_match.len() != self.points().len() {
            return Err(format!(
                "must_match has {} entries but the cloud has {} points",
                must_match.len(),
                self.points().len()
            )
            .into());
        }

        Ok(normal_estimation::estimate_by_neighborhood(
            self.points(),
            must_match,
            self.tree.get(),
            radius,
        ))
    }
}

impl KdTreeSearch<3> for CloudIndex3<'_> {
    fn nearest_one(&self, point: &impl crate::common::PCoords<3>) -> (usize, f64) {
        self.tree.get().nearest_one(point)
    }

    fn nearest(&self, point: &impl crate::common::PCoords<3>, count: usize) -> Vec<(usize, f64)> {
        self.tree.get().nearest(point, count)
    }

    fn within(&self, point: &impl crate::common::PCoords<3>, radius: f64) -> Vec<(usize, f64)> {
        self.tree.get().within(point, radius)
    }

    fn len(&self) -> usize {
        self.tree.get().len()
    }
}

impl PointCloudOverlap<Mesh3> for CloudIndex3<'_> {
    /// Find the indices of points in this cloud which overlap a mesh, by looking for reciprocity in
    /// the closest point in each direction.
    ///
    /// For each point `p_this` in this cloud, find the closest point on the mesh, `p_other`, then
    /// find the closest point to `p_other` back in this cloud, `p_recip`. Where the two clouds
    /// genuinely cover the same surface `p_recip` is `p_this`, so `p_this` counts as overlapping
    /// when the two are within `max_distance`.
    fn overlap_by_reciprocity(&self, mesh: &Mesh3, max_distance: f64) -> Vec<usize> {
        let mut result = Vec::new();
        for (i, p_this) in self.points().iter().enumerate() {
            let p_other = mesh.point_closest_to(p_this);

            let (k, _) = self.tree.get().nearest_one(&p_other);
            let p_recip = self.points()[k];

            if dist(p_this, &p_recip) < max_distance {
                result.push(i);
            }
        }

        result
    }
}

impl PointCloudOverlap<CloudIndex3<'_>> for CloudIndex3<'_> {
    /// Find the indices of points in this cloud which overlap another cloud, by looking for
    /// reciprocity in the closest point in each direction.
    ///
    /// See the mesh implementation for the reciprocity argument; the only difference is that
    /// `p_other` comes from a nearest-neighbor query into `other` rather than a surface projection.
    fn overlap_by_reciprocity(&self, other: &CloudIndex3, max_distance: f64) -> Vec<usize> {
        let mut result = Vec::new();
        for (i, p_this) in self.points().iter().enumerate() {
            let (j, _) = other.tree.get().nearest_one(p_this);
            let p_other = other.points()[j];

            let (k, _) = self.tree.get().nearest_one(&p_other);
            let p_recip = self.points()[k];

            if dist(p_this, &p_recip) < max_distance {
                result.push(i);
            }
        }

        result
    }
}
