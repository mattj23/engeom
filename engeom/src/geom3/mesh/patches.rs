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

/// Which connected patches of a mesh are worth keeping.
///
/// Every criterion is optional and they combine as an intersection: a patch survives only if it
/// passes all of the ones which are set. The default keeps everything.
///
/// Marked `#[non_exhaustive]` so new criteria will not break your code, which means a struct
/// literal only works inside this crate. Build one with the chained setters instead:
///
/// ```
/// use engeom::geom3::mesh::PatchFilter;
///
/// let filter = PatchFilter::default()
///     .with_min_faces(50)
///     .with_min_area_fraction(0.01);
/// ```
///
/// The three threshold criteria measure genuinely different things and scan data usually wants
/// more than one. `min_faces` is cheap but tracks tessellation density rather than physical size,
/// so it will keep a dense speckle and discard a coarsely tessellated real feature. `min_area` is
/// the measure of how much surface a patch represents. `min_aabb_diagonal` catches the patch which
/// covers almost no area but spans a long distance, such as the sliver a bad scan line drags out.
/// Prefer the two physical measures when the tessellation density is not uniform.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PatchFilter {
    /// Discard a patch with fewer than this many faces.
    pub min_faces: Option<usize>,

    /// Discard a patch whose summed face area is below this, in the mesh's own units.
    pub min_area: Option<f64>,

    /// Discard a patch whose axis-aligned bounding box diagonal is below this, in the mesh's own
    /// units.
    ///
    /// This is the criterion to reach for when the absolute size of the junk is known but its
    /// tessellation is not, which is the usual situation with a flying patch: it is some small
    /// object a long way from everything else, and its diameter is the thing you can predict.
    pub min_aabb_diagonal: Option<f64>,

    /// Discard a patch whose area is below this fraction of the largest patch's area.
    ///
    /// Unlike the absolute thresholds this needs no knowledge of the part's size or units, which
    /// makes it the right default when the same filter runs across parts of different scales. A
    /// value of `0.01` keeps anything within two orders of magnitude of the main body.
    pub min_area_fraction: Option<f64>,

    /// Keep only the largest patches by area, discarding the rest.
    ///
    /// Ranked by area rather than face count, so a coarsely tessellated body still outranks a
    /// finely tessellated speck. `Some(1)` is the common "throw away everything but the part"
    /// case, available as `PatchFilter::keep_largest()`.
    pub keep_largest_n: Option<usize>,
}

impl PatchFilter {
    /// Keep only the single largest patch by area and discard everything else.
    pub fn keep_largest() -> Self {
        Self {
            keep_largest_n: Some(1),
            ..Default::default()
        }
    }

    /// Set the minimum face count, returning the modified filter.
    pub fn with_min_faces(mut self, min_faces: usize) -> Self {
        self.min_faces = Some(min_faces);
        self
    }

    /// Set the minimum area, returning the modified filter.
    pub fn with_min_area(mut self, min_area: f64) -> Self {
        self.min_area = Some(min_area);
        self
    }

    /// Set the minimum bounding box diagonal, returning the modified filter.
    pub fn with_min_aabb_diagonal(mut self, min_aabb_diagonal: f64) -> Self {
        self.min_aabb_diagonal = Some(min_aabb_diagonal);
        self
    }

    /// Set the minimum area as a fraction of the largest patch's, returning the modified filter.
    pub fn with_min_area_fraction(mut self, min_area_fraction: f64) -> Self {
        self.min_area_fraction = Some(min_area_fraction);
        self
    }

    /// Set how many of the largest patches by area to keep, returning the modified filter.
    pub fn with_keep_largest_n(mut self, keep_largest_n: usize) -> Self {
        self.keep_largest_n = Some(keep_largest_n);
        self
    }

    /// Decide which patches survive, given their statistics.
    ///
    /// This is the whole of the filter's meaning, separated from the mesh so the semantics can be
    /// reasoned about and tested on their own.
    ///
    /// # Arguments
    ///
    /// * `stats`: the statistics of each patch, indexed by patch label
    ///
    /// returns: `Vec<bool>` of length `stats.len()`, `true` where the patch survives
    pub fn keep_flags(&self, stats: &[PatchStats]) -> Vec<bool> {
        let mut keep = vec![true; stats.len()];

        for (i, s) in stats.iter().enumerate() {
            if let Some(min) = self.min_faces
                && s.face_count < min
            {
                keep[i] = false;
            }

            if let Some(min) = self.min_area
                && s.area < min
            {
                keep[i] = false;
            }

            if let Some(min) = self.min_aabb_diagonal
                && s.aabb_diagonal() < min
            {
                keep[i] = false;
            }
        }

        if let Some(fraction) = self.min_area_fraction {
            let largest = stats.iter().map(|s| s.area).fold(0.0, f64::max);
            let threshold = largest * fraction;
            for (i, s) in stats.iter().enumerate() {
                if s.area < threshold {
                    keep[i] = false;
                }
            }
        }

        if let Some(n) = self.keep_largest_n {
            // Ties break by patch index, so the ranking is deterministic on a mesh with several
            // patches of identical area, which a fixture of repeated primitives will have.
            let mut order: Vec<usize> = (0..stats.len()).collect();
            order.sort_by(|&a, &b| {
                stats[b]
                    .area
                    .partial_cmp(&stats[a].area)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.cmp(&b))
            });

            for &i in order.iter().skip(n) {
                keep[i] = false;
            }
        }

        keep
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

    // ---- PatchFilter -----------------------------------------------------------------------

    /// Stand-in statistics, so the filter's semantics can be checked without building meshes.
    fn stat(face_count: usize, area: f64, span: f64) -> PatchStats {
        PatchStats {
            face_count,
            area,
            min: [0.0; 3],
            max: [span, 0.0, 0.0],
        }
    }

    #[test]
    fn default_filter_keeps_everything() {
        let stats = [stat(1, 0.001, 0.01), stat(100_000, 5000.0, 300.0)];
        assert_eq!(PatchFilter::default().keep_flags(&stats), vec![true, true]);
    }

    #[test]
    fn each_threshold_acts_on_its_own_measure() {
        // A dense speck and a coarse but physically large patch, which is the case that separates
        // the face count criterion from the two physical ones.
        let stats = [stat(5000, 0.5, 1.0), stat(20, 400.0, 120.0)];

        let by_faces = PatchFilter {
            min_faces: Some(100),
            ..Default::default()
        };
        assert_eq!(by_faces.keep_flags(&stats), vec![true, false]);

        let by_area = PatchFilter {
            min_area: Some(10.0),
            ..Default::default()
        };
        assert_eq!(by_area.keep_flags(&stats), vec![false, true]);

        let by_span = PatchFilter {
            min_aabb_diagonal: Some(10.0),
            ..Default::default()
        };
        assert_eq!(by_span.keep_flags(&stats), vec![false, true]);
    }

    /// Criteria intersect, so a patch must clear every one that is set.
    #[test]
    fn criteria_intersect() {
        let stats = [stat(5000, 0.5, 1.0), stat(20, 400.0, 120.0)];
        let both = PatchFilter {
            min_faces: Some(100),
            min_area: Some(10.0),
            ..Default::default()
        };

        assert_eq!(both.keep_flags(&stats), vec![false, false]);
    }

    #[test]
    fn area_fraction_is_relative_to_the_largest_patch() {
        let stats = [
            stat(10, 1000.0, 100.0),
            stat(10, 5.0, 5.0),
            stat(10, 0.5, 1.0),
        ];
        let filter = PatchFilter {
            min_area_fraction: Some(0.001),
            ..Default::default()
        };

        // The threshold sits at 1.0, so only the smallest falls below it.
        assert_eq!(filter.keep_flags(&stats), vec![true, true, false]);
    }

    /// The rank cut goes by area, so a coarsely tessellated body outranks a finely tessellated
    /// speck even though the speck has far more faces.
    #[test]
    fn keep_largest_ranks_by_area_not_face_count() {
        let stats = [stat(9000, 2.0, 2.0), stat(30, 900.0, 110.0)];

        assert_eq!(
            PatchFilter::keep_largest().keep_flags(&stats),
            vec![false, true]
        );
    }

    #[test]
    fn keep_largest_n_takes_the_top_n() {
        let stats = [
            stat(10, 5.0, 5.0),
            stat(10, 100.0, 50.0),
            stat(10, 50.0, 25.0),
            stat(10, 1.0, 1.0),
        ];
        let filter = PatchFilter {
            keep_largest_n: Some(2),
            ..Default::default()
        };

        assert_eq!(filter.keep_flags(&stats), vec![false, true, true, false]);
    }

    /// Equal areas must break by index rather than by whatever order the sort happens to produce,
    /// or a fixture of repeated primitives would give different answers between runs.
    #[test]
    fn keep_largest_breaks_ties_deterministically() {
        let stats = [stat(10, 7.0, 3.0), stat(10, 7.0, 3.0), stat(10, 7.0, 3.0)];
        let filter = PatchFilter {
            keep_largest_n: Some(2),
            ..Default::default()
        };

        for _ in 0..8 {
            assert_eq!(filter.keep_flags(&stats), vec![true, true, false]);
        }
    }

    #[test]
    fn keep_largest_n_beyond_the_patch_count_keeps_everything() {
        let stats = [stat(10, 5.0, 5.0), stat(10, 1.0, 1.0)];
        let filter = PatchFilter {
            keep_largest_n: Some(10),
            ..Default::default()
        };

        assert_eq!(filter.keep_flags(&stats), vec![true, true]);
    }

    #[test]
    fn an_empty_patch_list_filters_to_nothing() {
        assert!(PatchFilter::keep_largest().keep_flags(&[]).is_empty());
    }

    // ---- Applying the filter to a mesh -----------------------------------------------------

    /// A body with a small flyer parked well away from it, which is the shape of the problem on
    /// real scan data.
    fn body_with_flyer() -> Mesh3 {
        let mut mesh = Mesh3::create_box(20.0, 20.0, 20.0, false);
        let mut flyer = Mesh3::create_box(0.5, 0.5, 0.5, false);
        flyer.transform_in_place(&Iso3::translation(300.0, 0.0, 0.0));
        mesh.append_in_place(&flyer).unwrap();
        mesh
    }

    #[test]
    fn remove_small_patches_drops_the_flyer() {
        let mesh = body_with_flyer();
        assert_eq!(mesh.compute_patch_labels(None).unwrap().count(), 2);

        let cleaned = mesh
            .remove_small_patches(&PatchFilter::keep_largest())
            .unwrap();

        assert_eq!(cleaned.compute_patch_labels(None).unwrap().count(), 1);
        assert_eq!(cleaned.faces().len(), 12);

        // The body is 20 on a side, so its diagonal is sqrt(3) * 20. With the flyer still attached
        // the bounding box would span more than 300.
        let stats = cleaned
            .compute_patch_labels(None)
            .unwrap()
            .compute_stats(&cleaned)
            .unwrap();
        assert_relative_eq!(
            stats[0].aabb_diagonal(),
            (3.0_f64).sqrt() * 20.0,
            epsilon = 1e-9
        );
    }

    /// Each physical criterion can pick the flyer out on its own.
    #[test]
    fn absolute_thresholds_drop_the_flyer() {
        let mesh = body_with_flyer();

        for filter in [
            PatchFilter {
                min_area: Some(10.0),
                ..Default::default()
            },
            PatchFilter {
                min_aabb_diagonal: Some(5.0),
                ..Default::default()
            },
            PatchFilter {
                min_area_fraction: Some(0.01),
                ..Default::default()
            },
        ] {
            let cleaned = mesh.remove_small_patches(&filter).unwrap();
            assert_eq!(
                cleaned.faces().len(),
                12,
                "filter {:?} kept the wrong set",
                filter
            );
        }
    }

    /// Dropping whole faces keeps an index mapping back to the original, so attributes come
    /// through. This is the reason patch filtering does not need an attribute-loss opt-in.
    #[test]
    fn attributes_survive_patch_removal() {
        let mut mesh = body_with_flyer();
        let labels: Vec<u32> = (0..mesh.faces().len() as u32).collect();
        mesh.set_face_labels(Some(labels)).unwrap();

        let cleaned = mesh
            .remove_small_patches(&PatchFilter::keep_largest())
            .unwrap();

        let kept = cleaned.face_labels().expect("face labels should survive");
        assert_eq!(kept.len(), 12);

        // The body was appended first, so it holds the first twelve original face indices.
        assert_eq!(kept, &(0..12u32).collect::<Vec<_>>()[..]);
    }

    /// A filter with nothing to do returns the mesh unchanged rather than rebuilding it.
    #[test]
    fn a_filter_that_keeps_everything_is_a_no_op() {
        let mesh = body_with_flyer();
        let same = mesh.remove_small_patches(&PatchFilter::default()).unwrap();

        assert_eq!(same.faces().len(), mesh.faces().len());
        assert_eq!(same.points().len(), mesh.points().len());
    }

    #[test]
    fn discarding_every_patch_is_an_error() {
        let mesh = body_with_flyer();
        let filter = PatchFilter {
            min_faces: Some(1_000_000),
            ..Default::default()
        };

        assert!(mesh.remove_small_patches(&filter).is_err());
    }

    #[test]
    fn patch_mask_reports_what_would_be_kept() {
        let mesh = body_with_flyer();
        let mask = mesh.patch_mask(&PatchFilter::keep_largest()).unwrap();

        assert_eq!(mask.count_true(), 12);
        assert_eq!(mask.len(), mesh.faces().len());
        for face in 0..12 {
            assert!(mask.get(face));
        }
        for face in 12..mesh.faces().len() {
            assert!(!mask.get(face));
        }
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
