//! The best-effort decimation path: original vertices held inside the tolerance, and nothing else.
//!
//! This is the counterpart to [`super::guarantee`], and the difference between them is the whole
//! reason both exist. That one proves a two-sided bound over the surfaces and pays for it in
//! conservatism. This one deliberately does not, because a bound is not what every use of mesh
//! decimation needs, and paying for one you are not going to use is just a smaller mesh you did
//! not get.
//!
//! Reached through [`HalfEdgeMesh3::decimate_best_effort`].
//!
//! # What decides a collapse
//!
//! Currently, there's just one test used to allow/prevent a collapse. Every face carries the
//! original vertices it has absorbed, and a collapse hands the star's absorbed points to whichever
//! surviving face they land nearest; if that would put any of them further than the tolerance, the
//! collapse is refused. See [`Absorbed`].
//!
//! The positions never change, so nothing accumulates and nothing compounds, but this is obviously
//! a very limited subset of the original shape, and the test does not (currently) run in the
//! opposite direction.
//!
//! # Current state of the tolerance
//!
//! Tolerance controls the face count, but what happens with the _deviation_ seems to depend on
//! two things: how loose the tolerance is, and what kind of surface it is applied to. Only the
//! distribution shows either, so that is what is recorded here rather than the Hausdorff figure
//! alone. Distances are pooled over both directions, as multiples of the stated tolerance.
//!
//! <div class="warning">
//!
//! This is supposed to be a best-effort decimation tool, and it is.  I also spent most of my time
//! working on the decimator in [`super::guarantee`], so this one may not be fully baked. See the
//! results in this section for an idea of what you can expect for now.  I plan to spend more time
//! improving it in the future.
//!
//! </div>
//!
//! ## Engine blade, 43586 faces, 70 x 46 x 127 mm
//!
//! This mesh was a scan of an actual metal engine blade (from ebay, most likely from an obsolete
//! engine) that had most likley been rejected for repair.  It was already heavily smoothed and
//! decimated, so there's not as much for the decimator to work with.  It's also one of the main
//! test and sample meshes in the library.
//!
//! ```text
//! tol      faces      p50    p90    p99  p99.9     max   over 1x  over 2x  over 5x  over 10x
//! 0.001    85.9%    0.00x  2.85x 12.58x 24.78x   87.6x   19.712%  13.265%   5.415%    1.699%
//! 0.005    56.2%    0.28x  1.33x  3.52x  6.55x   17.3x   15.420%   4.492%   0.293%    0.016%
//! 0.010    36.7%    0.31x  0.94x  2.08x  3.53x    9.0x    8.613%   1.148%   0.015%    0.000%
//! 0.050     9.6%    0.26x  0.65x  1.00x  1.38x    3.7x    1.016%   0.009%   0.000%    0.000%
//! 0.100     5.1%    0.23x  0.60x  0.89x  1.13x    2.8x    0.318%   0.000%   0.000%    0.000%
//! 0.500     1.1%    0.21x  0.57x  0.83x  0.96x    4.2x    0.039%   0.002%   0.000%    0.000%
//! ```
//!
//! It looks like there's a bifurcation into two separate regimes:
//!
//! - **Loose tolerance -> thin tail.** At `0.05` and above, 99% of the surface is inside the
//!   number and essentially nothing is more than 2x. The mesh is broadly right with a few bad
//!   spots, which is what "best effort" is supposed to mean.
//!
//! - **Tight tolerance -> whole distribution.** At `0.001`, a fifth of the surface is outside the
//!   tolerance, 5% is past 5x it and 1.7% past 10x.
//!
//! The cause is a pair of floors in absolute terms which the tolerance cannot get under, about
//! `0.013` on p99 and `0.087` on the worst point. A vertex-level check cannot see what a collapse
//! does *between* the vertices, and tightening the number does not make it look. Above the floors
//! the whole distribution tracks the tolerance; below them only the median still does.
//!
//! ## Stanford bunny, 16127 faces, 0.2 x 0.2 x 0.1
//!
//! ```text
//! tol      faces      p50    p90    p99  p99.9     max   over 1x  over 2x  over 5x  over 10x
//! 0.0002   24.0%    0.29x  0.78x  1.69x  2.96x    6.2x    5.016%   0.536%   0.003%    0.000%
//! 0.0005    9.6%    0.26x  0.64x  1.04x  1.58x    3.7x    1.187%   0.022%   0.000%    0.000%
//! 0.0010    5.1%    0.24x  0.60x  0.86x  1.07x    2.2x    0.204%   0.000%   0.000%    0.000%
//! 0.0050    1.8%    0.25x  0.62x  0.84x  0.93x    2.1x    0.012%   0.000%   0.000%    0.000%
//! 0.0100    1.5%    0.23x  0.60x  0.80x  0.86x    2.5x    0.024%   0.005%   0.000%    0.000%
//! 0.0500    1.3%    0.10x  0.28x  0.37x  0.41x    0.4x    0.000%   0.000%   0.000%    0.000%
//! ```
//!
//! ## ATOS 5 reference scan, 734117 faces, 127 x 104 x 56 mm
//!
//! A structured light scan of a production part, from the private multi-align test data. This is
//! closer to the test case the module is aimed at, with the caveat that the Zeiss software was
//! run with polygonization defaults that do some decimation as part of the meshing process, so
//! the data represents a mesh with some triangle optimization already.
//!
//! ```text
//! tol      faces      p50    p90    p99  p99.9     max   over 1x  over 2x  over 5x  over 10x
//! 0.001    33.1%    0.30x  0.95x  2.76x  6.25x   57.2x    9.064%   2.090%   0.206%    0.017%
//! 0.005     8.9%    0.26x  0.65x  1.05x  1.75x    8.7x    1.285%   0.054%   0.001%    0.000%
//! 0.010     4.6%    0.25x  0.61x  0.90x  1.20x    6.1x    0.379%   0.004%   0.000%    0.000%
//! 0.050     1.4%    0.23x  0.57x  0.83x  0.94x    3.5x    0.018%   0.000%   0.000%    0.000%
//! 0.500     0.8%    0.18x  0.52x  0.78x  0.90x    4.1x    0.002%   0.000%   0.000%    0.000%
//! ```
//!
//! The guaranteed path on the same mesh reaches 67.7% / 26.7% / 16.3% / 4.9% / 1.4% of its faces
//! at those five tolerances. So this path is delivering **two to four times fewer
//! faces** while, from `0.005` upward, leaving under 1.3% of the surface outside the number and
//! essentially nothing past twice it.
//!
//! # Side-by-side comparison
//!
//! Engine blade, 43586 faces, same stated tolerance. `fwd` is original to result and `rev` is
//! result to original, both dense over the surfaces, from `difference`, which neither
//! gate consults. From the `best_effort_profile` harness in this module's tests.
//!
//! ```text
//!              ------------- guaranteed -------------    ------------ best effort ------------
//! tol          faces    time       fwd       rev         faces    time       fwd       rev
//! 0.001        99.8%    3.1s  0.000871  0.000926         85.8%    2.0s  0.084581  0.087650
//! 0.005        93.2%    3.2s  0.004687  0.004708         56.2%    2.7s  0.086538  0.086209
//! 0.010        77.9%    5.2s  0.009512  0.009467         36.7%    4.0s  0.077803  0.090397
//! 0.050        28.4%    8.1s  0.044231  0.042316          9.6%    4.2s  0.123640  0.180947
//! 0.100        17.5%    8.0s  0.082765  0.080127          5.1%    4.6s  0.180536  0.303376
//! 0.500         5.4%    7.6s  0.321812  0.306397          1.1%    5.6s  0.569687  0.893837
//! ```
//!
//! **Two to five times fewer faces, in about two thirds of the time.**

use super::{
    AlumPolyMesh, CollapseGate, DecimateReport, DecimateStats, DecimateTarget, Placement, Quadric,
    QuadricKind, StarFace, Structural, candidate_positions, run_decimation, set_boundary_locks,
};
use crate::common::barycentric::{barycentric_point, closest_barycentric};
use crate::geom3::half_edge3::HalfEdgeMesh3;
use crate::{Point3, Result, Vector3};
use alum::{FH, HH, Handle, HasIterators, HasTopology, VH};
use std::f64::consts::PI;

/// How best-effort decimation should behave.
///
/// Marked `#[non_exhaustive]`, so build one with [`BestEffortOpts::to_tolerance`],
/// [`BestEffortOpts::to_face_count`], or [`BestEffortOpts::to_ratio`] and the chained setters.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BestEffortOpts {
    /// How far an original *vertex* may end up from the result, in the mesh's own units.
    ///
    /// Exact on the thing it names: no original vertex ends up further than this from the result.
    /// **It is not a bound on how far the surface moves**, and how close it comes to being one
    /// depends on where you set it. Loose, around half the mesh's edge length or more, and 99% of
    /// the surface lands inside it with a worst point under 4x. Tight, and it stops meaning very
    /// much: at a fiftieth of that, a fifth of the surface ends up outside it.
    ///
    /// The module documentation has the distribution behind those two sentences, and is worth
    /// reading before choosing a number. Use [`super::DecimateOpts`] if the result is going to be
    /// measured.
    pub tol: f64,

    /// The largest angle, in radians, a face may be turned through by a collapse.
    ///
    /// Unlike [`tol`](BestEffortOpts::tol) this one is exact. It matters more here than on the
    /// guaranteed path, because it is the only test standing between an aggressive quadric and a
    /// folded surface. The default of `PI / 4` matches [`super::DecimateOpts`].
    pub max_normal_deviation: f64,

    /// The lowest triangle quality a collapse may produce, on a scale where an equilateral triangle
    /// is 1 and a degenerate one is 0. Also exact. The default of `0.01` matches
    /// [`super::DecimateOpts`].
    pub min_aspect: f64,

    /// Whether vertices on a boundary are held fixed.
    ///
    /// Worth more here than on the guaranteed path, because the first failure mode in the module
    /// documentation is specifically about outlines: a quadric charges nothing for sliding along a
    /// plane, and a flat boundary is exactly the case where that is free. Defaults to `true`.
    ///
    /// Each run decides this afresh: a run with it `false` frees a boundary an earlier run on the
    /// same mesh locked.
    pub lock_boundary: bool,

    /// Where the surviving vertex of a collapse is placed.
    pub placement: Placement,

    /// When to stop.
    ///
    /// [`DecimateTarget::FaceCount`] and [`DecimateTarget::Ratio`] are the usual way to drive this
    /// path, since a face budget is normally what a display mesh or a fixturing mesh is actually
    /// constrained by. The tolerance then acts as a floor which stops the count target doing
    /// something absurd.
    pub target: DecimateTarget,

    /// Which quadric formulation to accumulate.
    ///
    /// [`QuadricKind::Probabilistic`] adds an uncertainty term to the residual which the area
    /// weight does not account for, so the estimate reads high with it and the same `tol` decimates
    /// less. That is the right direction to be wrong in, but it does mean the two kinds are not
    /// interchangeable at a fixed tolerance.
    pub quadric: QuadricKind,
}

impl Default for BestEffortOpts {
    fn default() -> Self {
        Self {
            tol: 0.0,
            max_normal_deviation: PI / 4.0,
            min_aspect: 0.01,
            lock_boundary: true,
            placement: Placement::Optimal,
            target: DecimateTarget::ToTolerance,
            quadric: QuadricKind::Triangle,
        }
    }
}

impl BestEffortOpts {
    /// Decimate as far as the given estimated deviation allows.
    pub fn to_tolerance(tol: f64) -> Self {
        Self {
            tol,
            target: DecimateTarget::ToTolerance,
            ..Default::default()
        }
    }

    /// Decimate towards a face count, still limited by the given estimated deviation.
    pub fn to_face_count(tol: f64, faces: usize) -> Self {
        Self {
            tol,
            target: DecimateTarget::FaceCount(faces),
            ..Default::default()
        }
    }

    /// Decimate towards a fraction of the starting face count, still limited by the estimate.
    pub fn to_ratio(tol: f64, ratio: f64) -> Self {
        Self {
            tol,
            target: DecimateTarget::Ratio(ratio),
            ..Default::default()
        }
    }

    /// Set the maximum angle a face may be turned through, returning the modified options.
    pub fn with_max_normal_deviation(mut self, value: f64) -> Self {
        self.max_normal_deviation = value;
        self
    }

    /// Set the minimum triangle quality, returning the modified options.
    pub fn with_min_aspect(mut self, value: f64) -> Self {
        self.min_aspect = value;
        self
    }

    /// Set whether boundary vertices are locked, returning the modified options.
    pub fn with_lock_boundary(mut self, value: bool) -> Self {
        self.lock_boundary = value;
        self
    }

    /// Set where the surviving vertex is placed, returning the modified options.
    pub fn with_placement(mut self, value: Placement) -> Self {
        self.placement = value;
        self
    }

    /// Set which quadric formulation to use, returning the modified options.
    pub fn with_quadric(mut self, value: QuadricKind) -> Self {
        self.quadric = value;
        self
    }

    fn validate(&self) -> Result<()> {
        if self.tol.is_nan() || self.tol <= 0.0 {
            return Err(format!(
                "A best-effort decimation tolerance must be positive, got {}",
                self.tol
            )
            .into());
        }

        if self.max_normal_deviation.is_nan() || self.max_normal_deviation <= 0.0 {
            return Err("The maximum normal deviation must be positive".into());
        }

        self.target.validate()
    }
}

/// What a collapse would do, once it has been found acceptable.
pub(crate) struct Accepted {
    /// Where the merged vertex goes.
    position: Vector3,

    /// The quadric the merged vertex inherits.
    quadric: Quadric,

    /// The estimate at `position`, which is the queue key. Ordering by the same quantity the veto
    /// uses means the cheapest collapse by the queue's reckoning is also the one furthest from
    /// being refused.
    cost: f64,

    /// The estimate the merged vertex takes on, held to the running maximum along the merge chain.
    error: f64,

    /// Where each original vertex the star had absorbed ends up, as `(vertex, face)`.
    ///
    /// Worked out during the check rather than after it, because the check *is* this computation:
    /// finding the nearest surviving face of every absorbed point both measures the collapse and
    /// says who should own the point next. See [`Absorbed`].
    reassigned: Vec<(u32, FH)>,

    /// Every face the star had, which is the set whose ownership the collapse invalidates. Carried
    /// here because the collapse driver hands [`CollapseGate::commit`] the surviving vertex and
    /// nothing else, and by then the star is gone.
    cleared: Vec<FH>,
}

/// A decimator which both orders and decides collapses by the quadric error, in the normalized
/// form described in the module documentation.
///
/// Everything except that decision is shared with [`super::guarantee::ToleranceVolumeDecimator`]:
/// the queue, the collapse driver, the structural vetoes, the candidate positions, and the pinning
/// rule for locked vertices.
pub struct BestEffortDecimator {
    opts: BestEffortOpts,

    /// The accumulated quadric per vertex, indexed by vertex handle.
    ///
    /// Only [`Placement::Optimal`] reads these, to site the merged vertex. Nothing here judges a
    /// collapse by them; see the module documentation for the measurement that says why.
    quadrics: Vec<Quadric>,

    /// The running maximum estimate per vertex. A plain vector rather than a mesh property, since
    /// the quadrics it is derived from do not persist across runs either.
    error: Vec<f64>,

    /// Which vertices may not move, indexed by vertex handle. See
    /// [`candidate_positions`].
    pinned: Vec<bool>,

    /// Which original vertices each surviving face has absorbed.
    absorbed: Absorbed,

    /// The structural vetoes and the work counters, shared with the guaranteed path.
    structural: Structural,
}

/// The original vertices, and which face of the mesh as it stands now has absorbed each one.
///
/// # What this is for
///
/// Every other measurement in this module compares the star against the other star, which is
/// already-decimated geometry after the first collapse, so per-collapse errors compose over a run
/// in a way nothing local can see. This one does not: the positions here are the mesh as the run
/// found it and they never change, so the distance from an original vertex to the face which owns
/// it is measured from the truth every time it is measured, with nothing accumulated and nothing
/// compounding.
///
/// That turns a collapse test into a real invariant rather than an estimate. After every collapse
/// every original vertex is within the tolerance of the mesh, because a collapse which would
/// have put one further away was refused.
///
/// # Why it stays local
///
/// A point is attached to exactly one face, and only faces of the star move. A point owned by a
/// face outside the star therefore has the same distance it had before, so a collapse only has
/// to look at the points its own star holds. That is what buys the invariant without a global
/// query, and it is why the ownership has to be a partition rather than a per-face list of
/// whatever happens to be nearby.
///
/// Distance to the owning face is an upper bound on distance to the surface rather than the true
/// value, since some other face may be nearer. It is only ever wrong in the conservative
/// direction.
///
/// # What it does not say
///
/// This is a check in one direction, and only at vertices. It says every original vertex is
/// covered; it says nothing about anything going on _between_ its vertices, and nothing at all
/// about the reverse direction, where the new surface bridges a hole or spans a concavity the
/// original did not have. The `Kobbelt, Campagna & Seidel` quote in [`super`]'s documentation is
/// the argument that the first of those does not matter, and the paragraph after it is why I don't
/// buy that argument here. Use [`super::guarantee`] when you actually need a guaranteed bound.
struct Absorbed {
    /// Every original vertex position, indexed as the mesh indexed them at the start of the run.
    ///
    /// "Original" means the mesh this run was handed. A mesh which has been decimated before
    /// arrives with its previous run's output here, the same limitation the quadrics have, and for
    /// the same reason: neither survives in the mesh's own storage.
    points: Vec<Point3>,

    /// The original vertices each face owns, indexed by face handle.
    ///
    /// Sized once at construction. A collapse only ever deletes faces, so no index here is ever out
    /// of range, and the entries of deleted faces are emptied as their contents are handed on.
    owned: Vec<Vec<u32>>,
}

impl Absorbed {
    /// Attach every original vertex to one of its own incident faces, at distance zero.
    fn new(mesh: &AlumPolyMesh) -> Result<Self> {
        let points = mesh.points();
        let points = points
            .try_borrow()
            .map_err(|_| "Failed to borrow points while seeding the absorbed vertices")?;

        let mut owned = vec![Vec::new(); mesh.num_faces()];

        for v in mesh.vertices() {
            if let Some(f) = mesh.vf_ccw_iter(v).next() {
                owned[f.index() as usize].push(v.index());
            }
        }

        Ok(Self {
            points: mesh.vertices().map(|v| Point3::from(points[v])).collect(),
            owned,
        })
    }

    /// Every original vertex the star currently holds.
    fn in_star<'a>(&'a self, star: &'a [StarFace]) -> impl Iterator<Item = u32> + 'a {
        star.iter()
            .flat_map(|f| self.owned[f.f.index() as usize].iter().copied())
    }

    /// Hand the star's absorbed vertices to the faces the collapse left them nearest to.
    ///
    /// `cleared` is every face the star had, including the ones the collapse deleted; `reassigned`
    /// only ever names faces which survived it. Clearing first and filling second is what keeps the
    /// ownership a partition, which is the property the locality argument rests on.
    fn apply(&mut self, cleared: &[FH], reassigned: &[(u32, FH)]) {
        for f in cleared.iter() {
            self.owned[f.index() as usize].clear();
        }
        for (v, f) in reassigned.iter() {
            self.owned[f.index() as usize].push(*v);
        }
    }
}

impl BestEffortDecimator {
    /// Build a best-effort decimator for a half-edge mesh.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the half-edge mesh which will be decimated
    /// * `opts`: how decimation should behave
    ///
    /// returns: `Result<BestEffortDecimator>`
    pub fn new(mesh: &HalfEdgeMesh3, opts: BestEffortOpts) -> Result<Self> {
        opts.validate()?;

        let inner = mesh.as_alum();

        let pinned = {
            let vstatus = inner.vertex_status_prop();
            let vstatus = vstatus
                .try_borrow()
                .map_err(|_| "Failed to borrow vertex status")?;
            inner.vertices().map(|v| vstatus[v].feature()).collect()
        };

        Ok(Self {
            opts,
            quadrics: seed(inner, opts.quadric)?,
            error: vec![0.0; inner.num_vertices()],
            pinned,
            absorbed: Absorbed::new(inner)?,
            structural: Structural::new(opts.max_normal_deviation, opts.min_aspect),
        })
    }

    /// Where every original vertex the star holds would land, and how far away the furthest one
    /// would be, or `None` once that passes `cap`.
    ///
    /// One pass does both jobs, because they are the same computation: the nearest surviving face
    /// of an absorbed point is how far the collapse leaves that point from the surface, and it is
    /// also who should own the point afterwards. See [`Absorbed`].
    ///
    /// `surviving` is the star's faces which are still there after the collapse, paired with the
    /// handles they keep.
    fn coverage(
        &self,
        star: &[StarFace],
        surviving: &[(FH, [Point3; 3])],
        cap: f64,
        mut worst: f64,
    ) -> Option<(f64, Vec<(u32, FH)>)> {
        if surviving.is_empty() {
            return None;
        }

        let mut out = Vec::new();

        for v in self.absorbed.in_star(star) {
            let p = self.absorbed.points[v as usize];

            let mut best = f64::INFINITY;
            let mut owner = surviving[0].0;
            for (f, q) in surviving.iter() {
                let w = closest_barycentric(&q[0], &q[1], &q[2], &p);
                let d = (barycentric_point(&q[0], &q[1], &q[2], w) - p).norm();
                if d < best {
                    best = d;
                    owner = *f;
                    // Nothing nearer can change who owns the point in any way that matters, and
                    // the maximum cannot be raised by it either.
                    if best <= worst {
                        break;
                    }
                }
            }

            if best > worst {
                worst = best;
                if worst > cap {
                    return None;
                }
            }

            out.push((v, owner));
        }

        Some((worst, out))
    }

    /// What this collapse would do, or `None` if no candidate position is acceptable.
    fn evaluate(&self, mesh: &AlumPolyMesh, h: HH) -> Option<Accepted> {
        self.structural.bump_evaluation();
        let v1 = h.tail(mesh);
        let v2 = h.head(mesh);

        let points = mesh.points();
        let points = points.try_borrow().ok()?;

        let p1 = points[v1];
        let p2 = points[v2];

        let q = self.quadrics[v1.index() as usize] + self.quadrics[v2.index() as usize];
        let inherited = self.error[v1.index() as usize].max(self.error[v2.index() as usize]);

        let star = self.structural.stars(mesh, v1, v2, &points)?;
        if star.is_empty() {
            return None;
        }

        // The charge never falls below what the endpoints already carry, so if that alone is over
        // budget no position can rescue the collapse.
        if inherited > self.opts.tol {
            self.structural.bump_error_veto();
            return None;
        }

        let pinned = self
            .pinned
            .get(v2.index() as usize)
            .copied()
            .unwrap_or(false);

        for candidate in candidate_positions(self.opts.placement, pinned, &q, p1, p2) {
            self.structural.bump_acceptance();

            let Some(new) = self.structural.reshape(&star, v1, v2, candidate, &points) else {
                continue;
            };

            // The faces the collapse leaves behind, which are the ones the star's absorbed
            // vertices have to be handed on to.
            let surviving: Vec<(FH, [Point3; 3])> = star
                .iter()
                .zip(new.iter())
                .filter(|(f, _)| !f.vanishing)
                .map(|(f, q)| (f.f, *q))
                .collect();

            let Some((cost, reassigned)) = self.coverage(&star, &surviving, self.opts.tol, 0.0)
            else {
                self.structural.bump_error_veto();
                continue;
            };

            let error = cost.max(inherited);

            return Some(Accepted {
                position: candidate,
                quadric: q,
                cost,
                error,
                reassigned,
                cleared: star.iter().map(|f| f.f).collect(),
            });
        }

        None
    }
}

/// Accumulate the quadric and the area weight of every vertex, in one pass over the mesh.
///
/// The two have to be accumulated together over the same exact faces, or the division which turns
/// a residual back into a length is comparing a sum to the wrong denominator. Faces carrying both
/// endpoints of a collapse are counted twice when the two quadrics are added, which is what
/// Garland-Heckbert does and what `alum` does; the weights double-count the same faces the same
/// way, so the ratio stays right.
fn seed(mesh: &AlumPolyMesh, kind: QuadricKind) -> Result<Vec<Quadric>> {
    let points = mesh.points();
    let points = points
        .try_borrow()
        .map_err(|_| "Failed to borrow points while seeding quadrics")?;

    let n = mesh.num_vertices();
    let mut quadrics = vec![Quadric::default(); n];

    for v in mesh.vertices() {
        let i = v.index() as usize;
        for f in mesh.vf_ccw_iter(v) {
            for vs in mesh.triangulated_face_vertices(f) {
                let (a, b, c) = (points[vs[0]], points[vs[1]], points[vs[2]]);
                quadrics[i] += match kind {
                    QuadricKind::Triangle => Quadric::triangle_quadric(a, b, c),
                    QuadricKind::Probabilistic => {
                        let mean_edge = ((b - a).norm() + (c - b).norm() + (a - c).norm()) / 3.0;
                        Quadric::probabilistic_triangle_quadric(a, b, c, 0.1 * mean_edge)
                    }
                };
            }
        }
    }

    Ok(quadrics)
}

impl CollapseGate for BestEffortDecimator {
    type Accepted = Accepted;

    fn evaluate(&self, mesh: &AlumPolyMesh, h: HH) -> Option<Accepted> {
        BestEffortDecimator::evaluate(self, mesh, h)
    }

    fn cost(accepted: &Accepted) -> f64 {
        accepted.cost
    }

    fn position(accepted: &Accepted) -> Vector3 {
        accepted.position
    }

    fn commit(&mut self, v: VH, accepted: Accepted) -> Result<()> {
        let i = v.index() as usize;
        self.quadrics[i] = accepted.quadric;
        self.error[i] = accepted.error;
        self.absorbed.apply(&accepted.cleared, &accepted.reassigned);
        Ok(())
    }

    fn stats(&self) -> DecimateStats {
        self.structural.stats.get()
    }
}

impl HalfEdgeMesh3 {
    /// Decimate the mesh by estimated deviation rather than by a guaranteed bound.
    ///
    /// Faster and considerably more aggressive than [`HalfEdgeMesh3::decimate_guaranteed`], and
    /// the tolerance it takes is an estimate which the result routinely exceeds. Use it for
    /// display meshes and for geometry which scaffolds a measurement rather than carrying one; use
    /// [`HalfEdgeMesh3::decimate_guaranteed`] when the result is going to be measured.
    ///
    /// Read the `best_effort` module documentation for what the estimate is and the three specific
    /// ways it can be wrong.
    ///
    /// # Arguments
    ///
    /// * `opts`: how decimation should behave
    ///
    /// returns: `Result<DecimateReport>`
    pub fn decimate_best_effort(&mut self, opts: &BestEffortOpts) -> Result<DecimateReport> {
        opts.validate()?;

        set_boundary_locks(self, opts.lock_boundary)?;

        let mut gate = BestEffortDecimator::new(self, *opts)?;

        // Marked before the run rather than after it, because a run which fails partway through has
        // still collapsed edges the error volume knows nothing about. Rejected options never reach
        // here, so a refused call leaves the mesh as it found it.
        self.set_error_volume_stale();
        run_decimation(self, &mut gate, opts.target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::half_edge3::RepairOpts;
    use crate::geom3::half_edge3::difference::ShapeDifference;
    use crate::{Mesh3, Point3};

    /// An `n` by `n` grid of unit squares in the z=0 plane, each split into two triangles.
    fn flat_grid(n: usize) -> Mesh3 {
        let mut points = Vec::new();
        for i in 0..=n {
            for j in 0..=n {
                points.push(Point3::new(i as f64, j as f64, 0.0));
            }
        }

        let mut faces = Vec::new();
        let idx = |i: usize, j: usize| (i * (n + 1) + j) as u32;
        for i in 0..n {
            for j in 0..n {
                faces.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                faces.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }

        Mesh3::new(points, faces, false)
    }

    fn decimate_blade(opts: &BestEffortOpts) -> (Mesh3, Mesh3, DecimateReport) {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
        let report = he.decimate_best_effort(opts).unwrap();
        let out = Mesh3::try_from(&he).unwrap();
        (original, out, report)
    }

    /// A plane is what a quadric represents, so the estimate is zero everywhere and the only
    /// thing left standing between the grid and a handful of triangles is the boundary lock.
    #[test]
    fn a_flat_grid_collapses_to_almost_nothing() {
        let original = flat_grid(20);
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        let report = he
            .decimate_best_effort(&BestEffortOpts::to_tolerance(1.0e-9))
            .unwrap();

        assert!(
            report.faces_after < 100,
            "a flat grid should collapse freely, got {}",
            report
        );

        // Whatever it did, it did in the plane.
        let out = Mesh3::try_from(&he).unwrap();
        for p in out.points().iter() {
            assert!(p.z.abs() < 1.0e-12, "left the plane at {:?}", p);
        }
    }

    /// The one claim this module makes which is not an estimate: every original vertex ends up
    /// within the tolerance of the result.
    ///
    /// This is what the absorbed-vertex tracking buys and it is checked against a measurement which
    /// knows nothing about that tracking. [`ShapeDifference::old_vertices_to_new`] finds each original
    /// vertex's true nearest point on the whole decimated surface; the gate only ever knew the
    /// distance to the one face which owned the point, which is an upper bound on that. So the
    /// measured value has to come in under the tolerance with room to spare, and if it ever does
    /// not, the ownership has stopped being a partition somewhere.
    ///
    /// The reverse direction is deliberately not asserted here. It is the one this path does not
    /// promise, and `the_estimate_overshoots_the_tolerance` below is where that is written down.
    #[test]
    fn every_original_vertex_stays_covered() {
        for tol in [0.005, 0.05, 0.5] {
            let (original, out, report) = decimate_blade(&BestEffortOpts::to_tolerance(tol));
            let diff = ShapeDifference::new(&original, &out).unwrap();

            assert!(
                diff.old_vertices_to_new <= tol,
                "an original vertex ended up {:.6} from the result at a tolerance of {}, \
                 which the absorbed-vertex tracking is supposed to make impossible ({})",
                diff.old_vertices_to_new,
                tol,
                report
            );
        }
    }

    /// At a loose tolerance the exceedance is a thin tail rather than a broadly wrong mesh, which
    /// is the entire value proposition of this path and is not visible in the Hausdorff number.
    ///
    /// The worst point at `tol = 0.1` is 3x the tolerance, which taken alone reads like a failure.
    /// The distribution says otherwise: under a percent of the surface is outside the number at
    /// all and nothing is past twice it. Both facts are needed to describe what this does, so both
    /// are pinned.
    #[test]
    fn a_loose_tolerance_leaves_only_a_thin_tail() {
        let tol = 0.1;
        let (original, out, report) = decimate_blade(&BestEffortOpts::to_tolerance(tol));
        let diff = ShapeDifference::with_spread(&original, &out, tol).unwrap();
        let spread = diff.spread.unwrap();

        assert!(
            spread.p99 <= tol,
            "99% of the surface should be inside the tolerance, p99 was {:.6} ({})",
            spread.p99,
            report
        );
        assert!(
            spread.over[0] < 0.02,
            "under 2% of samples should exceed the tolerance, got {:.3}% ({})",
            spread.over[0] * 100.0,
            report
        );
        assert!(
            spread.over[1] < 0.001,
            "essentially nothing should exceed twice the tolerance, got {:.3}%: {}",
            spread.over[1] * 100.0,
            diff
        );
    }

    /// The tolerance controls the face count. It does not control the deviation, and this is the
    /// test which says so out loud rather than leaving a caller to find out.
    ///
    /// Two separate facts are pinned. The deviation always *exceeds* the tolerance, so nobody can
    /// mistake the number for a bound. And it is bounded on this fixture by a floor plus a loose
    /// multiple, where the floor is the part that matters: below about `0.05` the tolerance stops
    /// moving the deviation at all, because a vertex-only check cannot see what a collapse does
    /// between the vertices and no amount of tightening the number will make it look.
    ///
    /// The numbers are a statement about the engine blade and not about the method. Pinning them
    /// means a change which makes this meaningfully worse has to come and edit them.
    #[test]
    fn the_deviation_is_not_controlled_by_the_tolerance() {
        /// What the deviation does not go below however tight the tolerance gets, on this fixture.
        const FLOOR: f64 = 0.15;

        for tol in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5] {
            let (original, out, report) = decimate_blade(&BestEffortOpts::to_tolerance(tol));
            let worst = ShapeDifference::new(&original, &out).unwrap().hausdorff();

            assert!(
                worst > tol,
                "at tol {} the deviation came in under the tolerance, which would mean this \
                 fixture stopped exercising the difference between the two paths: {:.6}",
                tol,
                worst
            );
            assert!(
                worst < FLOOR + 10.0 * tol,
                "at tol {} the deviation was {:.6}, past the recorded {} + 10x envelope ({})",
                tol,
                worst,
                FLOOR,
                report
            );
        }
    }

    /// More aggressive than the guaranteed path at the same tolerance, which is the entire reason
    /// this module exists.
    #[test]
    fn it_decimates_further_than_the_guaranteed_path() {
        let tol = 0.01;

        let (_, _, best) = decimate_blade(&BestEffortOpts::to_tolerance(tol));

        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
        let guaranteed = he
            .decimate_guaranteed(&crate::geom3::half_edge3::DecimateOpts::to_tolerance(tol))
            .unwrap();

        assert!(
            (best.faces_after as f64) < 0.8 * guaranteed.faces_after as f64,
            "expected the best-effort path to be materially more aggressive; \
             best effort {} vs guaranteed {}",
            best,
            guaranteed
        );
    }

    /// A face count target stops where it was told to.
    #[test]
    fn a_face_count_target_is_respected() {
        let target = 5_000;
        let (_, _, report) = decimate_blade(&BestEffortOpts::to_face_count(1.0, target));

        assert!(
            report.faces_after <= target + 2,
            "expected to stop near {}, got {}",
            target,
            report
        );
    }

    /// The tolerance still binds when a count target asks for more than it allows.
    #[test]
    fn the_gate_wins_over_an_unreachable_count() {
        let (_, _, report) = decimate_blade(&BestEffortOpts::to_face_count(1.0e-5, 100));

        assert!(
            report.faces_after > 100,
            "the tolerance should have stopped this well short of 100 faces, got {}",
            report
        );
    }

    /// Locking the boundary keeps the outline, which matters more here than on the guaranteed path
    /// because a quadric charges nothing for sliding along a plane.
    #[test]
    fn a_locked_boundary_is_not_moved() {
        let original = flat_grid(12);
        let before = original
            .points()
            .iter()
            .filter(|p| p.x <= 0.0 || p.y <= 0.0 || p.x >= 12.0 || p.y >= 12.0)
            .count();

        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
        he.decimate_best_effort(&BestEffortOpts::to_tolerance(1.0e-9))
            .unwrap();
        let out = Mesh3::try_from(&he).unwrap();

        let after = out
            .points()
            .iter()
            .filter(|p| p.x <= 0.0 || p.y <= 0.0 || p.x >= 12.0 || p.y >= 12.0)
            .count();

        assert_eq!(
            before, after,
            "a locked boundary should keep every one of its vertices"
        );
    }

    /// The two paths side by side on the blade, which is where the module doc's table comes from.
    ///
    /// ```text
    /// cargo test -r -p engeom --lib best_effort_profile -- --ignored --nocapture
    /// ```
    /// Where the deviation actually sits, rather than only how bad its worst point is.
    ///
    /// The `max` column is the Hausdorff distance and it is the number a tolerance claim is
    /// judged on, but on a method with no bound it is one sample out of seven million and it
    /// cannot say whether the mesh is broadly fine with a few bad spots or broadly wrong. The
    /// percentiles and the exceedance fractions can.
    ///
    /// ```text
    /// cargo test -r -p engeom --lib best_effort_distribution -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "measurement harness, not an assertion"]
    fn best_effort_distribution() {
        distribution_table(
            "engine blade",
            &crate::tests::engine_blade(),
            &[0.001, 0.005, 0.01, 0.05, 0.1, 0.5],
        );
    }

    /// The same distribution on a much more uniformly triangulated mesh.
    ///
    /// The bunny's edge lengths span a factor of two from p10 to p90 against the blade's ten, and
    /// it is open where the original scans did not reach, so it exercises the boundary handling
    /// the blade cannot.
    ///
    /// ```text
    /// cargo test -r -p engeom --lib bunny_distribution -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "measurement harness, not an assertion"]
    fn bunny_distribution() {
        distribution_table(
            "stanford bunny",
            &crate::tests::stanford_bun_2(),
            &[0.0002, 0.0005, 0.001, 0.005, 0.01, 0.05],
        );
    }

    /// The same distribution on a real industrial scan, which is the case this path is aimed at.
    ///
    /// An ATOS 5 structured light scan of a production part, 734117 faces in millimetres, from the
    /// private multi-align test data. Large, noisy, full of holes, and nothing like a CAD surface.
    /// Slow: the guaranteed path alone is several minutes per tolerance at this size.
    ///
    /// ```text
    /// cargo test -r -p engeom --lib --features ply,private_tests \
    ///     reference_scan_distribution -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "measurement harness, needs the private test data"]
    #[cfg(all(feature = "ply", feature = "private_tests"))]
    fn reference_scan_distribution() {
        let path = std::path::Path::new("../../engeom-test-data")
            .join("private-multi-align")
            .join("sample_00")
            .join("reference.ply");

        let data = crate::io::load_ply_mesh_data(&path).unwrap();
        let mesh = Mesh3::from_data(data, false).unwrap();
        distribution_table(
            "ATOS 5 reference scan",
            &mesh,
            &[0.001, 0.005, 0.01, 0.05, 0.5],
        );
    }

    /// Both paths' deviation distributions on one mesh, over a range of tolerances.
    ///
    /// The mesh is repaired into half-edge form once up front and that is what everything is
    /// measured against, so the numbers are about decimation rather than about repair. Scan data
    /// generally needs it; the blade does not and is unchanged by it.
    fn distribution_table(label: &str, input: &Mesh3, tols: &[f64]) {
        // Flushed after every row. Rust buffers stdout in blocks when it is not a terminal, and
        // these runs are long enough that waiting on the buffer means seeing nothing until the end.
        use std::io::Write;
        let (repaired, _) = HalfEdgeMesh3::from_mesh_repaired(input, &RepairOpts::default())
            .expect("failed to repair the fixture into half-edge form");
        let original = &Mesh3::try_from(&repaired).unwrap();
        drop(repaired);

        let aabb = original.aabb();
        println!(
            "\n=== {} ===\n{} faces, {} points, extents {:.1} x {:.1} x {:.1}\n\n\
             distances pooled over both directions, as multiples of the stated tolerance.\n\n\
             {:>7} {:>10} {:>6} {:>6} {:>6} {:>6} {:>7}   {:>8} {:>8} {:>8} {:>8}",
            label,
            original.faces().len(),
            original.points().len(),
            aabb.extents().x,
            aabb.extents().y,
            aabb.extents().z,
            "tol",
            "faces",
            "p50",
            "p90",
            "p99",
            "p99.9",
            "max",
            "over 1x",
            "over 2x",
            "over 5x",
            "over 10x"
        );

        for &tol in tols {
            for guaranteed in [true, false] {
                let mut he = HalfEdgeMesh3::try_from(original).unwrap();
                let report = if guaranteed {
                    he.decimate_guaranteed(&crate::geom3::half_edge3::DecimateOpts::to_tolerance(
                        tol,
                    ))
                } else {
                    he.decimate_best_effort(&BestEffortOpts::to_tolerance(tol))
                }
                .unwrap();

                let out = Mesh3::try_from(&he).unwrap();
                let g = ShapeDifference::with_spread(original, &out, tol).unwrap();
                let s = g.spread.unwrap();

                println!(
                    "{:>7} {:>10} {:>5.1}% {:>5.2}x {:>5.2}x {:>5.2}x {:>5.2}x {:>6.1}x   \
                     {:>7.3}% {:>7.3}% {:>7.3}% {:>7.3}%",
                    tol,
                    if guaranteed { "guaranteed" } else { "best" },
                    report.ratio() * 100.0,
                    s.p50 / tol,
                    s.p90 / tol,
                    s.p99 / tol,
                    s.p999 / tol,
                    s.max / tol,
                    s.over[0] * 100.0,
                    s.over[1] * 100.0,
                    s.over[2] * 100.0,
                    s.over[3] * 100.0
                );
                std::io::stdout().flush().ok();
            }
        }
    }

    #[test]
    #[ignore = "measurement harness, not an assertion"]
    fn best_effort_profile() {
        let original = crate::tests::engine_blade();
        println!(
            "\n{} faces.  fwd is original -> result, rev is result -> original.\n\n\
             {:>6}  {:>6} {:>6} {:>9} {:>9}  {:>6} {:>6} {:>9} {:>9}",
            original.faces().len(),
            "",
            "faces",
            "time",
            "fwd",
            "rev",
            "faces",
            "time",
            "fwd",
            "rev"
        );
        println!(
            "{:>6}  {:>33}  {:>33}",
            "tol", "----- guaranteed -----", "----- best effort ----"
        );

        for tol in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5] {
            let mut a = HalfEdgeMesh3::try_from(&original).unwrap();
            let start = std::time::Instant::now();
            let ra = a
                .decimate_guaranteed(&crate::geom3::half_edge3::DecimateOpts::to_tolerance(tol))
                .unwrap();
            let ta = start.elapsed().as_secs_f64();
            let ga = ShapeDifference::new(&original, &Mesh3::try_from(&a).unwrap()).unwrap();

            let mut b = HalfEdgeMesh3::try_from(&original).unwrap();
            let start = std::time::Instant::now();
            let rb = b
                .decimate_best_effort(&BestEffortOpts::to_tolerance(tol))
                .unwrap();
            let tb = start.elapsed().as_secs_f64();
            let gb = ShapeDifference::new(&original, &Mesh3::try_from(&b).unwrap()).unwrap();

            println!(
                "{:>6}  {:>5.1}% {:>5.2}s {:>8.6} {:>8.6}  {:>5.1}% {:>5.2}s {:>8.6} {:>8.6}",
                tol,
                ra.ratio() * 100.0,
                ta,
                ga.dense.reference_to_test,
                ga.dense.test_to_reference,
                rb.ratio() * 100.0,
                tb,
                gb.dense.reference_to_test,
                gb.dense.test_to_reference
            );
        }
    }

    #[test]
    fn a_bad_tolerance_is_refused() {
        let mut he = HalfEdgeMesh3::try_from(&flat_grid(4)).unwrap();
        assert!(
            he.decimate_best_effort(&BestEffortOpts::to_tolerance(0.0))
                .is_err()
        );
        assert!(
            he.decimate_best_effort(&BestEffortOpts::to_tolerance(-1.0))
                .is_err()
        );
    }
}
