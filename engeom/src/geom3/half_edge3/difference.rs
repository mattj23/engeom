//! Brute-force measurement of how much two meshes differ, which is what every tolerance claim in
//! the rest of `half_edge3` gets checked against.
//!
//! Everything here is deliberately naive and shares no machinery with the operations it verifies.
//! It's all test-only, and it's hella slow on purpose.
//!
//! # Why does this exist?
//!
//! These checks came out of the mesh decimation work when I was forced to confront the deeper and
//! more annoying question of "what does _difference_ between shapes actually mean?".  Originally it
//! lived in that submodule, but got moved once I realized I was going to have to re-tread this same
//! ground on every mesh editing operation.
//!
//! In the big picture, this exists because any operation in which I claim a tolerance has to be
//! checked by _something_, and that something needs to not share any of its internal machinery
//! (obvious in retrospect, but naturally I learned it the hard way). To that end, everything in
//! here is deliberately naive so that it can be obviously correct.  That's the idea, at least.
//!
//! # What _does_ "difference" even mean?
//!
//! Colloquially, the idea of two shapes being different seems pretty intuitive; it's a normal part
//! of the everyday human experience, and we have (reasonably or otherwise) very strong spatial
//! processing and intuition.
//!
//! As a computational formality, it's not that simple.  The idea splits into two, completely
//! separate questions that have to be answered independently:
//!
//! 1. Does shape B **represent** shape A?
//! 2. Does shape B **_not_ misrepresent** shape A?
//!
//! Let's dip into mathematical formalism for a moment, and define "shape" in the above context to
//! mean "a _set_ of points in Rⁿ space". Most likely they're all part of a group of 1 or more
//! continuous manifolds, but we don't actually need that to be true.
//!
//! Question #1 asks whether every point in A is acceptably close to a point in B.  If you copied A
//! and then removed some part of the copy, you'd answer "no" to this question.
//!
//! Question #2 asks whether any point in B is _not_ acceptably close to a point in A. If you copied
//! A and then added something new onto the copy, you'd answer "no" to this question.
//!
//! The two questions are independent: you can easily come up with scenarios where you'd answer
//! every one of the four possible combinations of "yes" and "no".
//!
//! And, unfortunately, when it comes to handling measurement data, answering "no" to either
//! question is violating the integrity of the measurement. A "no" to the first question means
//! you've modified or removed some of the values that _were_ measured.  A "no" to the second
//! question means you've claimed the measurement has stuff it originally didn't.
//!
//! # Approximation of the difference check
//!
//! Meshes can be thought of as "the set of all points contained by the union of all triangle faces
//! in the mesh" (assuming that you define the face as also containing its edges and vertices).
//! I believe that _rigorously_ testing all points in a triangle set requires, at the very least:
//!
//! - Distance bounding volume tools for the triangle faces (basically the ones defined in Geuziec's
//!   _Surface Simplification Inside a Tolerance Volume_)
//! - Tools for creating 2D boundary outlines of the intersection of those bounding volumes with
//!   arbitrary 3D planes, consisting of arc and line segment edges, and the tools for doing union
//!   and subtraction operations on them.
//!
//! That's theoretically possible, but a lot of work and I have no confidence in getting it right
//! without mistakes.  It's also unclear to me how computationally expensive it would be.
//!
//! As an alternative, this module chooses to _estimate_ the test by using point-to-mesh distance
//! queries, which are well understood and already supported by this library. To approximate "the
//! set of all points" representing by a mesh, I use dense sampling to generate points across the
//! faces of the triangles.
//!
//! The two questions from the section above resolve to a pair of reciprocal distance queries. For
//! each direction, the entire set of densely sampled points are queried against the other mesh,
//! and the worst case value is the Hausdorff distance estimate.
//!
//! # The checks in this module
//!
//! The main check in the module, [`ShapeDifference::dense`], is the acceptance criterion for any
//! operation that is trying to give a tolerance guarantee. It is the two-sided Hausdorff distance
//! between the surfaces, approximated by sampling as described above.  Assertions about tolerances
//! should be written against it.
//!
//! [`ShapeDifference::old_vertices_to_new`] and [`ShapeDifference::new_vertices_to_old`] are weaker,
//! utility measures that are done _without_ dense sampling. Instead, they just use the vertices
//! of the meshes. As a check, they're much closer to what's in the literature and what's commonly
//! implemented out in the wild.  You can assert against them if you're not trying to give a full
//! shape tolerance guarantee.
//!
//! [`ShapeDifference::spread`] is a helpful tool for characterization when developing a method, the
//! idea is to use the common machinery to get a distribution of distances rather than just the
//! worst case.
//!
//! [`ShapeDifference::boundary`] and [`ShapeDifference::volume_ratio`] catch two things the Hausdorff
//! distance is blind to: an outline that creeps inward while staying inside tolerance, and a part
//! that shrinks systematically. Both are applicable to metrology and neither shows up in a
//! deviation number.

use crate::geom3::mesh::SurfaceDeviation;
use crate::geom3::mesh::measurement::{derived_spacing, sample_distances};
use crate::{Mesh3, Point3, Result};

/// How many points each boundary segment is sampled at when comparing outlines.
///
/// The endpoints alone would miss an outline which bows between its vertices, which is the same
/// vertices-are-not-the-surface mistake in one dimension lower.
const BOUNDARY_SAMPLES_PER_SEGMENT: usize = 4;

/// The multiples of the threshold that [`Spread::over`] counts samples against.
pub(crate) const OVER_MULTIPLES: [f64; 4] = [1.0, 2.0, 5.0, 10.0];

/// The distribution of dense sample distances, pooled over both directions.
///
/// Pooled rather than per-direction because a piece of surface being wrong is equally a defect
/// whichever way it was measured, and the pooled fraction over a threshold is directly the "how
/// much of this mesh is out of tolerance" number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Spread {
    /// The median sample distance.
    pub p50: f64,

    /// The 90th percentile sample distance.
    pub p90: f64,

    /// The 99th percentile sample distance.
    pub p99: f64,

    /// The 99.9th percentile sample distance.
    ///
    /// Here because `p99` and `max` together cannot tell a thin tail of outliers apart from a
    /// tenth of the mesh having moved, and that distinction is the whole question when judging a
    /// method with no bound.
    pub p999: f64,

    /// The largest sample distance, which is the Hausdorff estimate itself.
    pub max: f64,

    /// The fraction of samples further than the threshold, and than 2x, 5x and 10x it.
    ///
    /// Indexed by [`OVER_MULTIPLES`]. For a method with a real bound every entry is zero by
    /// construction. For a best-effort method the shape of these four numbers is what says whether
    /// exceeding the threshold is a thin tail or a wholesale failure to respect it.
    pub over: [f64; 4],
}

impl std::fmt::Display for Spread {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "p50 {:.6} p90 {:.6} p99 {:.6} p99.9 {:.6} max {:.6}, \
             over {:.3}% / {:.3}% / {:.3}% / {:.3}% at 1x / 2x / 5x / 10x",
            self.p50,
            self.p90,
            self.p99,
            self.p999,
            self.max,
            self.over[0] * 100.0,
            self.over[1] * 100.0,
            self.over[2] * 100.0,
            self.over[3] * 100.0
        )
    }
}

/// A complete independent measurement of what an editing operation did to a mesh.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ShapeDifference {
    /// The worst distance from an original vertex to the new surface. _Not for use on methods that
    /// are trying to give strong tolerance guarantees._
    pub old_vertices_to_new: f64,

    /// The worst distance from a new vertex to the original surface. _Not for use on methods that
    /// are trying to give strong tolerance guarantees._
    pub new_vertices_to_old: f64,

    /// Both surfaces sampled densely and each measured against the other.
    pub dense: SurfaceDeviation,

    /// The distribution behind `dense`, present only when it was asked for.
    pub spread: Option<Spread>,

    /// The two-sided distance between the two meshes' boundary outlines.
    ///
    /// `None` when neither mesh has a boundary, since two watertight meshes have no outlines to
    /// compare. Infinite when one has a boundary and the other does not, which is a sign that
    /// something went horribly, horribly wrong.
    pub boundary: Option<f64>,

    /// The enclosed volume after divided by the volume before, or `None` unless both meshes are
    /// closed.
    ///
    /// A decimator that stays inside tolerance while consistently cutting corners off is biased,
    /// and surface deviation can't directly see it. Gueziec's Section 6 is about preserving this.
    pub volume_ratio: Option<f64>,
}

impl ShapeDifference {
    /// Measure `new` against `old`, without the sample distribution.
    ///
    /// # Arguments
    ///
    /// * `old`: the mesh as it was before the operation
    /// * `new`: the mesh the operation produced
    ///
    /// returns: `Result<ShapeDifference>`
    pub(crate) fn new(old: &Mesh3, new: &Mesh3) -> Result<Self> {
        Self::build(old, new, None)
    }

    /// Measure `new` against `old`, including the distribution of sample distances.
    ///
    /// Costs a second dense sampling pass and retains one `f64` per sample, so this is for
    /// characterizing a method rather than for asserting on one.
    ///
    /// # Arguments
    ///
    /// * `old`: the mesh as it was before the operation
    /// * `new`: the mesh the operation produced
    /// * `threshold`: the distance `Spread::over` counts samples against, normally whatever
    ///   tolerance the operation was asked for
    ///
    /// returns: `Result<ShapeDifference>`
    pub(crate) fn with_spread(old: &Mesh3, new: &Mesh3, threshold: f64) -> Result<Self> {
        Self::build(old, new, Some(threshold))
    }

    fn build(old: &Mesh3, new: &Mesh3, threshold: Option<f64>) -> Result<Self> {
        let dense = old.measure_surface_deviation(new, None)?;

        let spread = match threshold {
            Some(t) => {
                let mut all = sample_distances(old, new, derived_spacing(old)?, "reference")?;
                all.extend(sample_distances(new, old, derived_spacing(new)?, "test")?);
                Some(Spread::from_samples(&mut all, t))
            }
            None => None,
        };

        Ok(Self {
            old_vertices_to_new: worst_vertex_distance(old, new),
            new_vertices_to_old: worst_vertex_distance(new, old),
            dense,
            spread,
            boundary: boundary_deviation(old, new),
            volume_ratio: volume_ratio(old, new),
        })
    }

    /// The symmetric Hausdorff estimate, which is what a tolerance claim means.
    pub(crate) fn hausdorff(&self) -> f64 {
        self.dense.hausdorff()
    }

    /// Panic unless both dense directions are within `tol`, reporting everything measured.
    pub(crate) fn assert_within(&self, tol: f64) {
        assert!(
            self.hausdorff() <= tol,
            "surfaces deviate by more than the tolerance of {}\n  {}",
            tol,
            self
        );
    }
}

impl std::fmt::Display for ShapeDifference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "dense: {}", self.dense)?;
        write!(
            f,
            "\n  vertices only (diagnostic): old->new {:.6}, new->old {:.6}",
            self.old_vertices_to_new, self.new_vertices_to_old
        )?;
        if let Some(s) = self.spread {
            write!(f, "\n  spread: {}", s)?;
        }
        if let Some(b) = self.boundary {
            write!(f, "\n  boundary: {:.6}", b)?;
        }
        if let Some(v) = self.volume_ratio {
            write!(f, "\n  volume ratio: {:.6}", v)?;
        }
        Ok(())
    }
}

impl Spread {
    /// Build a spread from unsorted sample distances, sorting them in place.
    fn from_samples(samples: &mut [f64], threshold: f64) -> Self {
        samples.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        if samples.is_empty() {
            return Self {
                p50: 0.0,
                p90: 0.0,
                p99: 0.0,
                p999: 0.0,
                max: 0.0,
                over: [0.0; 4],
            };
        }

        let at = |q: f64| samples[(((samples.len() - 1) as f64) * q).round() as usize];
        let n = samples.len() as f64;
        let fraction_over =
            |m: f64| samples.iter().filter(|d| **d > threshold * m).count() as f64 / n;

        Self {
            p50: at(0.5),
            p90: at(0.9),
            p99: at(0.99),
            p999: at(0.999),
            max: at(1.0),
            over: std::array::from_fn(|i| fraction_over(OVER_MULTIPLES[i])),
        }
    }
}

/// The worst distance from any vertex of `from` to the surface of `to`.
///
/// This is the weak measure kept only for comparison. See the module documentation.
fn worst_vertex_distance(from: &Mesh3, to: &Mesh3) -> f64 {
    from.points()
        .iter()
        .map(|p| to.distance_closest_to(p))
        .fold(0.0f64, f64::max)
}

/// The two-sided distance between two meshes' boundary outlines.
fn boundary_deviation(old: &Mesh3, new: &Mesh3) -> Option<f64> {
    let a = boundary_segments(old);
    let b = boundary_segments(new);

    match (a.is_empty(), b.is_empty()) {
        (true, true) => None,
        (true, false) | (false, true) => Some(f64::INFINITY),
        (false, false) => Some(one_sided_boundary(&a, &b).max(one_sided_boundary(&b, &a))),
    }
}

/// The boundary edges of a mesh, as pairs of positions.
fn boundary_segments(mesh: &Mesh3) -> Vec<[Point3; 2]> {
    let points = mesh.points();
    mesh.compute_nav()
        .boundary_edges(None)
        .into_iter()
        .map(|e| [points[e[0] as usize], points[e[1] as usize]])
        .collect()
}

/// The worst distance from any point sampled along `from` to the nearest segment of `to`.
fn one_sided_boundary(from: &[[Point3; 2]], to: &[[Point3; 2]]) -> f64 {
    let mut worst = 0.0f64;

    for seg in from.iter() {
        for i in 0..BOUNDARY_SAMPLES_PER_SEGMENT {
            let t = i as f64 / (BOUNDARY_SAMPLES_PER_SEGMENT - 1) as f64;
            let p = seg[0] + (seg[1] - seg[0]) * t;

            let nearest = to
                .iter()
                .map(|s| point_to_segment(&p, &s[0], &s[1]))
                .fold(f64::INFINITY, f64::min);

            worst = worst.max(nearest);
        }
    }

    worst
}

/// The distance from a point to a line segment.
fn point_to_segment(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = b - a;
    let denom = ab.norm_squared();
    let t = if denom <= 0.0 {
        0.0
    } else {
        ((p - a).dot(&ab) / denom).clamp(0.0, 1.0)
    };

    (p - (a + ab * t)).norm()
}

/// The enclosed volume ratio, or `None` unless both meshes are closed.
fn volume_ratio(old: &Mesh3, new: &Mesh3) -> Option<f64> {
    if !boundary_segments(old).is_empty() || !boundary_segments(new).is_empty() {
        return None;
    }

    let before = signed_volume(old);
    if before.abs() <= f64::MIN_POSITIVE {
        return None;
    }

    Some(signed_volume(new) / before)
}

/// The signed volume enclosed by a closed mesh, by the divergence theorem.
///
/// Each face contributes the signed volume of the tetrahedron it forms with the origin, which sums
/// to the enclosed volume when the surface is closed and consistently wound.
fn signed_volume(mesh: &Mesh3) -> f64 {
    let points = mesh.points();
    let total: f64 = mesh
        .faces()
        .iter()
        .map(|f| {
            let a = points[f[0] as usize].coords;
            let b = points[f[1] as usize].coords;
            let c = points[f[2] as usize].coords;
            a.dot(&b.cross(&c))
        })
        .sum();

    total / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A mesh compared against itself has to come back at zero on every measure, or the harness is
    /// reporting its own noise as deviation.
    #[test]
    fn a_mesh_matches_itself_exactly() {
        let mesh = crate::tests::engine_blade();
        let g = ShapeDifference::new(&mesh, &mesh).unwrap();

        // Not exactly zero: projecting a point onto a triangle it already sits on leaves a
        // residue around 1e-13, so the statement is "nothing beyond arithmetic noise".
        assert!(g.hausdorff() < 1.0e-9, "{}", g);
        assert!(g.old_vertices_to_new < 1.0e-9, "{}", g);
        assert!(g.new_vertices_to_old < 1.0e-9, "{}", g);
        g.assert_within(1.0e-9);
    }

    /// The reason the whole harness exists, as an executable statement: two surfaces can agree at
    /// every single vertex and still deviate substantially in between.
    ///
    /// Two triangulations of the same saddle quad share all four vertices, so both vertex measures
    /// are exactly zero while the surfaces cross in the middle. A check built on vertices would call
    /// these meshes identical.
    #[test]
    fn vertex_measures_miss_what_the_dense_measure_catches() {
        let points = vec![
            Point3::new(-1.0, -1.0, 1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(-1.0, 1.0, -1.0),
        ];

        let one = Mesh3::new(points.clone(), vec![[0, 1, 2], [0, 2, 3]], false);
        let two = Mesh3::new(points, vec![[0, 1, 3], [1, 2, 3]], false);

        let g = ShapeDifference::new(&one, &two).unwrap();

        assert!(g.old_vertices_to_new < 1.0e-9, "the vertices are shared");
        assert!(g.new_vertices_to_old < 1.0e-9, "the vertices are shared");
        assert!(g.hausdorff() > 0.1, "the surfaces genuinely differ: {}", g);
    }

    /// A shape moved bodily by a known amount should measure that amount, which is the calibration
    /// that says the units and the directions are right.
    #[test]
    fn a_known_offset_measures_as_itself() {
        let a = Mesh3::create_box(4.0, 4.0, 4.0, false);
        let b = a.transform_copy(&crate::Iso3::translation(0.25, 0.0, 0.0));

        let g = ShapeDifference::new(&a, &b).unwrap();

        assert!(
            (g.hausdorff() - 0.25).abs() < 1.0e-9,
            "expected 0.25, got {}",
            g
        );
    }

    /// Both meshes are closed, so the volume ratio has to be present and, for a rigid motion, one.
    #[test]
    fn a_rigid_motion_preserves_the_volume() {
        let a = Mesh3::create_box(4.0, 4.0, 4.0, false);
        let b = a.transform_copy(&crate::Iso3::translation(0.25, 0.0, 0.0));

        let ratio = ShapeDifference::new(&a, &b).unwrap().volume_ratio.unwrap();
        assert!((ratio - 1.0).abs() < 1.0e-9, "got {}", ratio);
    }

    /// A mesh with a boundary has no volume to speak of, and two watertight meshes have no outline
    /// to compare. Each measure should decline rather than invent a number.
    #[test]
    fn the_optional_measures_decline_when_they_do_not_apply() {
        let closed = Mesh3::create_box(4.0, 4.0, 4.0, false);
        let g = ShapeDifference::new(&closed, &closed).unwrap();
        assert!(g.boundary.is_none(), "a closed mesh has no outline");
        assert!(g.volume_ratio.is_some());

        let open = crate::tests::stanford_bun_4();
        let g = ShapeDifference::new(&open, &open).unwrap();
        assert!(g.boundary.is_some(), "an open mesh has an outline");
        assert!(g.boundary.unwrap() < 1.0e-9, "matched against itself");
        assert!(g.volume_ratio.is_none(), "an open mesh encloses nothing");
    }

    /// The boundary measure has to see an outline move even when the surface has not, since that is
    /// the failure it was added for.
    #[test]
    fn the_boundary_measure_sees_an_outline_move() {
        let square = |w: f64| {
            Mesh3::new(
                vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(w, 0.0, 0.0),
                    Point3::new(w, 1.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                ],
                vec![[0, 1, 2], [0, 2, 3]],
                false,
            )
        };

        let g = ShapeDifference::new(&square(1.0), &square(0.7)).unwrap();
        let b = g.boundary.unwrap();
        assert!((b - 0.3).abs() < 1.0e-9, "expected 0.3, got {}", b);
    }

    /// The spread is only computed on request, and its maximum has to agree with the Hausdorff
    /// value the criterion reports from a separate pass.
    #[test]
    fn the_spread_agrees_with_the_criterion() {
        let a = Mesh3::create_box(4.0, 4.0, 4.0, false);
        let b = a.transform_copy(&crate::Iso3::translation(0.25, 0.0, 0.0));

        assert!(ShapeDifference::new(&a, &b).unwrap().spread.is_none());

        let g = ShapeDifference::with_spread(&a, &b, 0.1).unwrap();
        let s = g.spread.unwrap();

        assert!((s.max - g.hausdorff()).abs() < 1.0e-9, "{}", g);
        assert!(s.p50 <= s.p90 && s.p90 <= s.p99 && s.p99 <= s.p999 && s.p999 <= s.max);

        // A rigid 0.25 offset in x slides four of the box's six faces along their own planes, so
        // only the two perpendicular ones move at all and about a third of the samples are over
        // the 0.1 threshold. Over 2x it (0.2) they still are, and over 5x it (0.5) none are, which
        // brackets the real distance from both sides.
        assert!(
            s.over[0] > 0.3,
            "a 0.25 offset is over a 0.1 threshold: {}",
            g
        );
        assert!(s.over[1] > 0.3, "and over 2x it: {}", g);
        assert_eq!(s.over[2], 0.0, "but nothing is over 5x it: {}", g);
        assert!(
            s.over.windows(2).all(|w| w[0] >= w[1]),
            "the exceedance fractions must fall as the multiple rises: {}",
            g
        );
    }
}
