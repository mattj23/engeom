//! This module holds the results of splitting a mesh into connected patches, and the summary
//! statistics used to rank them.
//!
//! The labeling itself lives on `MeshNav`, because it needs the edge/face adjacency. What lives
//! here is the *representation* of the result. A patch decomposition used to come back as a
//! `Vec<IndexMask>`, one mask per patch, which costs `patches * faces` bits. That is fine when a
//! mesh splits into a handful of pieces and quietly ruinous on the workload that actually needs
//! it: a dirty scan with thousands of small flying patches, where the masks alone can outweigh
//! the mesh by orders of magnitude. `PatchLabels` stores the same information as one `u32` per
//! face regardless of how many patches there are, and can still hand out masks on request for the
//! callers that want them.

use crate::common::IndexMask;
use crate::common::points::triangle_area;
use crate::{Mesh3, Result};

/// A per-face patch decomposition of a mesh.
///
/// Every face carries the index of the patch it belongs to, or [`PatchLabels::NO_PATCH`] if it was
/// excluded by the mask the labeling ran under. Patch indices are contiguous from `0` to
/// `count() - 1`.
///
/// Patches are numbered in order of their lowest-indexed face, so the decomposition is
/// deterministic for a given mesh and mask.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatchLabels {
    labels: Vec<u32>,
    count: usize,
}

impl PatchLabels {
    /// The label carried by a face which is not part of any patch, because the mask the labeling
    /// ran under excluded it.
    pub const NO_PATCH: u32 = u32::MAX;

    /// Build the labeling directly from a label buffer and a patch count.
    ///
    /// This is for the code doing the flood fill; ordinary callers get a `PatchLabels` from
    /// `MeshNav::patch_labels`.
    ///
    /// # Arguments
    ///
    /// * `labels`: one label per face, or `NO_PATCH` for an excluded face
    /// * `count`: the number of distinct patches
    ///
    /// returns: `Result<PatchLabels>`, failing if any label is outside `0..count`
    pub fn new(labels: Vec<u32>, count: usize) -> Result<Self> {
        for (i, &label) in labels.iter().enumerate() {
            if label != Self::NO_PATCH && label as usize >= count {
                return Err(format!(
                    "Face {} carries patch label {}, which is out of range for {} patches",
                    i, label, count
                )
                .into());
            }
        }

        Ok(Self { labels, count })
    }

    /// The patch label of every face, in face order. Faces excluded by the mask carry
    /// [`PatchLabels::NO_PATCH`].
    pub fn labels(&self) -> &[u32] {
        &self.labels
    }

    /// The number of distinct patches found.
    pub fn count(&self) -> usize {
        self.count
    }

    /// The number of faces the labeling covers, which is the face count of the mesh it ran on.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether the labeling covers no faces at all.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// The patch a face belongs to, or `None` if it was excluded by the mask.
    ///
    /// # Arguments
    ///
    /// * `face`: the index of the face
    ///
    /// returns: `Option<usize>`
    pub fn label_of(&self, face: usize) -> Option<usize> {
        match self.labels.get(face) {
            Some(&label) if label != Self::NO_PATCH => Some(label as usize),
            _ => None,
        }
    }

    /// Count the faces in each patch, indexed by patch label.
    ///
    /// returns: `Vec<usize>` of length `count()`
    pub fn face_counts(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.count];
        for &label in self.labels.iter() {
            if label != Self::NO_PATCH {
                counts[label as usize] += 1;
            }
        }
        counts
    }

    /// Build a face mask selecting a single patch.
    ///
    /// # Arguments
    ///
    /// * `patch`: the patch label
    ///
    /// returns: `Result<IndexMask>`, failing if the label is out of range
    pub fn mask_of(&self, patch: usize) -> Result<IndexMask> {
        if patch >= self.count {
            return Err(format!(
                "Patch {} is out of range for a labeling with {} patches",
                patch, self.count
            )
            .into());
        }

        let mut mask = IndexMask::new(self.labels.len(), false);
        let target = patch as u32;
        for (i, &label) in self.labels.iter().enumerate() {
            if label == target {
                mask.set(i, true);
            }
        }

        Ok(mask)
    }

    /// Build a face mask selecting every patch whose label passes a predicate.
    ///
    /// This is the cheap way to act on a decomposition, because it walks the label buffer once and
    /// allocates a single mask no matter how many patches survive.
    ///
    /// # Arguments
    ///
    /// * `keep`: called with each patch label, returning whether that patch should be selected
    ///
    /// returns: `IndexMask` over the faces
    pub fn mask_where(&self, keep: impl Fn(usize) -> bool) -> IndexMask {
        // Decide once per patch rather than once per face, since a predicate may be doing real
        // work and there are typically far fewer patches than faces.
        let keep: Vec<bool> = (0..self.count).map(&keep).collect();

        let mut mask = IndexMask::new(self.labels.len(), false);
        for (i, &label) in self.labels.iter().enumerate() {
            if label != Self::NO_PATCH && keep[label as usize] {
                mask.set(i, true);
            }
        }

        mask
    }

    /// Expand the labeling into one face mask per patch.
    ///
    /// Prefer `mask_where` where it will do. This allocates `count()` masks of `len()` bits each,
    /// which is the cost the labeling exists to avoid.
    ///
    /// returns: `Vec<IndexMask>` of length `count()`
    pub fn to_masks(&self) -> Vec<IndexMask> {
        let mut masks = vec![IndexMask::new(self.labels.len(), false); self.count];
        for (i, &label) in self.labels.iter().enumerate() {
            if label != Self::NO_PATCH {
                masks[label as usize].set(i, true);
            }
        }
        masks
    }

    /// Compute the summary statistics of every patch in a single pass over the faces.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh the labeling was computed on, whose face count must match
    ///
    /// returns: `Result<Vec<PatchStats>>` indexed by patch label
    pub fn compute_stats(&self, mesh: &Mesh3) -> Result<Vec<PatchStats>> {
        let faces = mesh.faces();
        if faces.len() != self.labels.len() {
            return Err(format!(
                "A labeling over {} faces does not match a mesh with {} faces",
                self.labels.len(),
                faces.len()
            )
            .into());
        }

        let points = mesh.points();
        let mut stats = vec![PatchStats::empty(); self.count];

        for (i, &label) in self.labels.iter().enumerate() {
            if label == Self::NO_PATCH {
                continue;
            }

            let face = &faces[i];
            let entry = &mut stats[label as usize];
            let p = [
                &points[face[0] as usize],
                &points[face[1] as usize],
                &points[face[2] as usize],
            ];

            entry.face_count += 1;
            entry.area += triangle_area(p[0], p[1], p[2]);
            for point in p {
                for axis in 0..3 {
                    entry.min[axis] = entry.min[axis].min(point[axis]);
                    entry.max[axis] = entry.max[axis].max(point[axis]);
                }
            }
        }

        Ok(stats)
    }
}

/// Summary measures of a single connected patch, used to decide whether it is worth keeping.
///
/// The three measures answer different questions and a filter usually wants more than one of them.
/// `face_count` is what tessellation density gives you and says nothing about physical size;
/// `area` is the honest measure of how much surface a patch represents; `aabb_diagonal` catches
/// the patch that covers little area but spans a long distance, such as a sliver of stray points
/// dragged out by a bad scan line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PatchStats {
    /// The number of faces in the patch.
    pub face_count: usize,

    /// The summed area of the patch's faces.
    pub area: f64,

    /// The lower corner of the patch's axis-aligned bounding box.
    pub min: [f64; 3],

    /// The upper corner of the patch's axis-aligned bounding box.
    pub max: [f64; 3],
}

impl PatchStats {
    /// The identity for accumulation: no faces, no area, and an inverted bounding box which any
    /// real point will correct.
    fn empty() -> Self {
        Self {
            face_count: 0,
            area: 0.0,
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        }
    }

    /// The length of the diagonal of the patch's axis-aligned bounding box, which is zero for an
    /// empty patch.
    pub fn aabb_diagonal(&self) -> f64 {
        if self.face_count == 0 {
            return 0.0;
        }

        let mut total = 0.0;
        for axis in 0..3 {
            let d = self.max[axis] - self.min[axis];
            total += d * d;
        }
        total.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Iso3;
    use approx::assert_relative_eq;

    /// Two boxes far enough apart to be separate patches, with exactly known area and extent. The
    /// larger box is built first, so it takes patch 0 under the lowest-face numbering rule.
    fn two_boxes() -> Mesh3 {
        let mut mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let mut other = Mesh3::create_box(1.0, 1.0, 1.0, false);
        other.transform_in_place(&Iso3::translation(100.0, 0.0, 0.0));
        mesh.append_in_place(&other).unwrap();
        mesh
    }

    #[test]
    fn stats_of_known_boxes() {
        let mesh = two_boxes();
        let labels = mesh.compute_patch_labels(None).unwrap();
        assert_eq!(labels.count(), 2);

        let stats = labels.compute_stats(&mesh).unwrap();

        // A closed box tessellates as two triangles per side.
        assert_eq!(stats[0].face_count, 12);
        assert_eq!(stats[1].face_count, 12);

        // Surface area of a rectangular prism is 2(lw + lh + wh).
        assert_relative_eq!(stats[0].area, 2.0 * (8.0 + 12.0 + 24.0), epsilon = 1e-10);
        assert_relative_eq!(stats[1].area, 6.0, epsilon = 1e-10);

        assert_relative_eq!(stats[0].aabb_diagonal(), 56.0_f64.sqrt(), epsilon = 1e-10);
        assert_relative_eq!(stats[1].aabb_diagonal(), 3.0_f64.sqrt(), epsilon = 1e-10);

        // The small box sits 100 out along x, so its bounding box must have moved with it.
        assert_relative_eq!(stats[1].min[0], 99.5, epsilon = 1e-10);
        assert_relative_eq!(stats[1].max[0], 100.5, epsilon = 1e-10);
    }

    #[test]
    fn face_counts_match_stats() {
        let mesh = two_boxes();
        let labels = mesh.compute_patch_labels(None).unwrap();

        let counts = labels.face_counts();
        let stats = labels.compute_stats(&mesh).unwrap();

        assert_eq!(counts.len(), labels.count());
        for (count, stat) in counts.iter().zip(stats.iter()) {
            assert_eq!(*count, stat.face_count);
        }
        assert_eq!(counts.iter().sum::<usize>(), mesh.faces().len());
    }

    /// `mask_where` is the allocation-light path and must select exactly what picking the same
    /// patches out of `to_masks` would have.
    #[test]
    fn mask_where_matches_expanded_masks() {
        let mesh = two_boxes();
        let labels = mesh.compute_patch_labels(None).unwrap();
        let masks = labels.to_masks();

        for keep in [0usize, 1] {
            let selected = labels.mask_where(|p| p == keep);
            assert_eq!(selected, masks[keep]);
            assert_eq!(selected, labels.mask_of(keep).unwrap());
        }

        let all = labels.mask_where(|_| true);
        assert_eq!(all.count_true(), mesh.faces().len());

        let none = labels.mask_where(|_| false);
        assert_eq!(none.count_true(), 0);
    }

    #[test]
    fn mask_of_rejects_an_out_of_range_patch() {
        let mesh = two_boxes();
        let labels = mesh.compute_patch_labels(None).unwrap();

        assert!(labels.mask_of(labels.count()).is_err());
    }

    #[test]
    fn stats_reject_a_mesh_of_the_wrong_size() {
        let mesh = two_boxes();
        let labels = mesh.compute_patch_labels(None).unwrap();
        let smaller = Mesh3::create_box(1.0, 1.0, 1.0, false);

        assert!(labels.compute_stats(&smaller).is_err());
    }

    #[test]
    fn new_rejects_a_label_beyond_the_count() {
        assert!(PatchLabels::new(vec![0, 1, 2], 2).is_err());
        assert!(PatchLabels::new(vec![0, 1, PatchLabels::NO_PATCH], 2).is_ok());
    }

    /// An excluded face carries no label and contributes to no patch's statistics.
    #[test]
    fn masked_faces_are_excluded_from_stats() {
        let mesh = two_boxes();

        // Keep only the large box, which is the first twelve faces.
        let mut mask = IndexMask::new(mesh.faces().len(), false);
        for i in 0..12 {
            mask.set(i, true);
        }

        let labels = mesh.compute_patch_labels(Some(&mask)).unwrap();
        assert_eq!(labels.count(), 1);

        let stats = labels.compute_stats(&mesh).unwrap();
        assert_eq!(stats[0].face_count, 12);
        assert_relative_eq!(stats[0].area, 88.0, epsilon = 1e-10);

        for face in 12..mesh.faces().len() {
            assert_eq!(labels.label_of(face), None);
        }
    }
}
