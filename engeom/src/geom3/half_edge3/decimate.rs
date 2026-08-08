//! Mesh decimation; this module and its submodules provide mesh decimation...perhaps better called
//! "mesh thinning" or "mesh pruning"...which is a set of tools for reducing the complexity of
//! 2-simplices in 3D space while preserving its shape.
//!
//! > This module turned out to be a bad case of "how hard can it be?" syndrome. I had thought,
//! > because mesh decimation is a well studied problem and a commonly implemented tool, that this
//! > wouldn't be all that big of a deal to implement. I figured I'd be wrapping some calls to
//! > `alum`, or at most writing a small trait implementation that plugged into another library's
//! > machinery, and then I'd be off to work on some other problem.  You can guess how that turned
//! > out.
//!
//! # Shape difference, tolerance, and decimation
//!
//! In metrology contexts mesh decimation is supposed to be a measurement-preserving operation. All
//! it should do is reduce the size of the data needed to faithfully represent a measurement, which
//! it does by stripping out geometric redundancy. You provide a tolerance that means _"changes in
//! the measurement below this distance mean nothing in my use case"_, and the decimation prunes
//! faces while adhering to the tolerance.
//!
//! See the "What _does_ difference even mean?" section of the `difference` module to go down the
//! hole of what that means and discuss the two-sided measurement of shape. The quick summary
//! is that while the _idea_ of a tolerance is easy, measuring it from inside an algorithm that's
//! supposed to not violate it is not, and requires a bidirectional measurement of the worst-case
//! closest distance between the entire surface of the original and the entire surface of the
//! simplification, aka the symmetrical Hausdorff distance.
//!
//! Unfortunately, the strong two-sided guarantee isn't what's out there in the world or filling
//! up the literature.  Most decimation implementations (a) use an estimate in place of actual
//! measured distances, and (b) provide at best a one-sided guarantee, usually from the decimated
//! mesh back to the original mesh.
//!
//! > To me, the current benchmark for measurement software aimed at 3D meshes is Zeiss Inspect
//! > (formerly GOM Inspect).
//! >
//! > Out of pure curiosity I took a sample mesh of a camera mounting bracket scanned on an ATOS 5,
//! > with some flat faces and internal threads, and thinned it with a 0.0254mm tolerance.  I used
//! > the bidirectional dense sampling check in `difference` to verify the result. While 99% of
//! > the points were within the tolerance, there were failures in both the original-to-thinned and
//! > thinned-to-original directions, with the former having the worst violations.  The worst
//! > violation not near the mesh boundary was 0.062mm, almost 2.5x the supplied tolerance.
//! >
//! > My take on this is that even metrology-focused software doesn't take the notion of a tolerance
//! > guarantee seriously.
//!
//! # Two implementations
//!
//! In the end, this module came to have two different implementations.  The TL;DR explaination is
//! that, while there is a functioning tolerance-guaranteed version that I'm very proud of
//! implemented here, it is naturally (a) more conservative, and (b) slower.  This produced an
//! internal conflict in the library against the fact that, while a tolerance guarantee is the
//! right tool for reducing any measurement data, that's not the only use of mesh decimation in
//! the metrology world.
//!
//! There are still cases where all you want is to reduce the number of triangles on a mesh,
//! and it's fine so long as it doesn't produce something dramatically different from what you
//! started with.  Visualization is an obvious case, and I've spent at least as much time over the
//! course of my career creating reports and visualizations as I have programming measurements and
//! analysis.  But there are also plenty of uses for reduced meshes when it comes to building
//! fixturing for actual measurements: mesh and patch subsets used for prealignments, constructing
//! sampling or selection regions, partitioning regions of space, and defining rough geometry to
//! pre-seed queries and other operations...to name a few.
//!
//! For that reason there are two submodules of this one:
//!
//! - [`guarantee`] has the metrology path with the tolerance guarantee, reached through
//!   [`HalfEdgeMesh3::decimate_guaranteed`]. Its tolerance is a genuine two-sided bound over the
//!   surfaces, not just at sampled points or the original vertices: no part of the result strays
//!   further than the tolerance from the original, and no part of the original is left further
//!   than the tolerance from the result.
//!
//! - [`best_effort`] is for everything else; it holds every original vertex inside the tolerance
//!   and makes no promises about the surface between them, which is enough to prevent completely
//!   unbounded shape change but is nothing like a two-sided bound. This should be enough for
//!   display meshes and for measurement scaffolding/fixturing applications, and it is two to five
//!   times more aggressive in about two thirds of the time. Set the tolerance loose, around half
//!   the mesh's edge length or more, and 99% of the surface still lands inside it. As of right now,
//!   tight tolerances aren't well respected by the algorithm.  See the module documentation for
//!   the error distributions to get a sense of how the tolerance influences shape change.
//!
//! Both share the guts of this module: the collapse driver, the queue, the structural vetoes, and
//! the quadric bookkeeping. They differ only in what decides a collapse, which is the
//! [`CollapseGate`] trait.
//!
//! ---
//!
//! The rest of this module doc is just here to record why the guaranteed path ended up how it did.
//!
//! # What's in the literature
//!
//! As far as guarantees on distance-based mesh decimation errors go, the bulk of the literature
//! seems to be from the late 90s to early 2010s. Most mesh decimation literature seems to be
//! very old and focused on the development of the data structures and the pruning mechanisms, or
//! more recent and focused on computer graphics and visualization and the use of softer quality
//! metrics and general efficiency.
//!
//! ## Literature on the tolerance guarantee
//!
//! Where there was literature, not all of it was helpful.  For instance, in **Kobbelt, Campagna &
//! Seidel, _A General Framework for Mesh Decimation_, GI '98.**, the authors argue that two-sided
//! measurement isn't necessary because:
//!
//! > Although the two-sided Hausdorff distance fits the intuitive notion of the deviation of one
//! > geometric object from the other very well, it is not appropriate for most of the typical
//! > input data to mesh reduction algorithms.
//! >
//! > The reason for this is that in general only the vertices of the given mesh represent actually
//! > measured points. This is true for laser-range scanned data and data obtained from mechanical
//! > probing. The initial triangulation that recovers the neighborhood relations in the input to
//! > our reduction algorithm has typically been generated by a preprocessing algorithm which
//! > itself is based on heuristic decisions (and not on specific knowledge about the object).
//! > Hence, there is no point in approximating the whole piecewise linear surface but it is enough
//! > to approximate the discrete data points themselves.
//!
//! This is certainly not true anymore, and realistically it may not have been true in 1998 either,
//! when their countrymen at GOM GmbH in Braunschweig were working on structured light scanners
//! that...by the time my career started in the early 00's...had industry trusted polygonization
//! algorithms that used the middle area of mesh faces to carry the measurement evidence of points
//! that were culled during processing.
//!
//! On the other hand, papers like **Klein, Liebich & Strasser, _Mesh Reduction with Error Control_,
//! VIS '96.** recognized the need for the symmetric Hausdorff distance, and their _Figure 1_ in
//! is a perfect demonstration of cases where the one-sided, original-to-simplifed measurement
//! fails to catch a surface that has moved by more than the allowable tolerance.  Their mechanism
//! for dealing with the two-sided measurement was to calculate the upper bound of the Haussdorff
//! distance by dealing with each triangle through the lens of one of three cases, for which the
//! max can be analytically calculated:
//!
//! 1. When all three vertices project onto the same triangle
//! 2. When all three vertices project to two triangles sharing an edge
//! 3. All other cases, resolved by subdividing the triangle until all of its subdivisions fall
//!    into case #1 or #2, and then taking the max value of the set
//!
//! I don't think I got that one implemented correctly, and in any case it seemed to degenerate
//! quickly into recursive subdivisions when the decimation produced a large reduction in faces.
//!
//! ## André Guéziec's method
//!
//! The gem ended up being an IBM research paper from 1997 by André Guéziec, titled  **_Surface
//! Simplification Inside a Tolerance Volume_**.  It's long but a worthwhile read. He explicitly
//! identifies the need for the symmetric Hausdorff distance and builds a method around performing
//! topology operations in a way guaranteed not to violate it.
//!
//! The paper contains a tightly reasoned proof for the soundness of enforcing a tolerance guarantee
//! by tracking error volumes (probably better called 'uncertainty' volumes) as changes are made to
//! a mesh. The volumes are represented by a single scalar value at each vertex of the mesh as it's
//! being simplified, and when an edge collapse occurs, (a) the new volume has to enclose the entire
//! pre-collapse volume, and (b) at the location of the Hausdorff distance the post-collapse sphere
//! must contain the pre-collapse sphere.  If the size of the error volume reaches the tolerance
//! value, the collapse is vetoed.
//!
//! As long as you meet those two critera for each collapse, the entirety of the new and old
//! surfaces are guaranteed to be mutually inside each other's tolerance zones. There are no global
//! queries, no local state being tracked beyond the single scalar, you don't even need to remember
//! the location of the original mesh.
//!
//! The tradeoff for this compactness is that the method achieves its soundness through being
//! conservative: as the algorithm progresses the knowledge of where the original surface was
//! gets fuzzier and fuzzier, and collapses that might have been fine get rejected because the
//! worst-case _possible_ place the Hausdorff distance might be would exceed the tolerance and
//! negate your ability to _prove_ that it was still acceptable. And, because single scalar
//! per-vertex values cannot represent a complicated shape, the error volume is also absorbing space
//! known to _not_ belong to it when it grows.
//!
//! That said, the method does work and provides conservative but guaranteed decimation within the
//! tolerance given.  See the [`guarantee`] module for a full description of the approach.
//!
//! # What a commercial package does with the same question
//!
//! As mentioned above, I tried thinning a scan of a real part in Zeiss Inspect Professional 2026
//! to see what it did with the tolerance value it requires for the operation.
//!
//! The scan was of a camera mounting bracket (the type that gets screwed onto the bottom of a
//! camera to hold it in a tripod), done on an ATOS system. It had been polygonized from the raw
//! shots with some of Zeiss' default postprocessing enabled, so it's important to note that it's
//! already a sparse mesh compared to the original raw point sampling from the scanner, but it's
//! representative of what people using high end structured light systems to measure industrial
//! components are actually getting as the output of 3D scanning. Geometrically, it had a bunch of
//! threaded through-holes, and the scan captured a fair amount of the threads. There were also
//! a large number of topological holes, which is pretty normal for industrial scans where you only
//! worry about capturing the surfaces you need for your original purpose.
//!
//! In the Zeiss Inspect software I thinned it with a 0.0254mm tolerance, then compared it to the
//! original mesh using `difference::ShapeDifference`. I didn't use Zeiss' built-in surface
//! comparison because, if you are familiar with how it works, you know that it only samples
//! difference at mesh vertices, which doesn't tell you anything about the surfaces represented
//! by the middle of the faces.
//!
//! ```text
//! original   82275 faces / 41706 points
//! thinned    31594 faces / 16336 points   (38.4%)
//! stated     0.001 in = 0.025400 mm
//!
//! d(original -> thinned)   0.162967 mm    6.42x    over
//! d(thinned -> original)   0.046689 mm    1.84x    over
//! ```
//!
//! The forward comparison (original to thinned) has a Hausdorff distance more than 6x the
//! tolerance I typed into the user interface.  Thankfully, though, the distribution of differences
//! shows that this is part of a thin tail rather than a wholesale failure to respect the
//! tolerance:
//!
//! ```text
//! p50   0.003135 mm   0.12x
//! p90   0.008920 mm   0.35x
//! p99   0.015751 mm   0.62x
//! max   0.162967 mm   6.42x
//!
//! 1136 of 3519586 forward samples exceed the tolerance   (0.032%)
//! ```
//!
//! The expected usual suspect is boundary simplification, since it's kind of a toss-up whether
//! a software package considers that to be part of a surface tolerance, and I don't actually know
//! how Zeiss treats it.  The four worst violations were all about 0.22mm from a boundary, so there
//! is something there that explains part of it.  That said, the worst interior deviation was
//! 0.062mm, which is pretty large compared to the tolerance.
//!
//! I'm not here to bash Zeiss, to me their software is the current gold standard for metrology on
//! 3D meshes, and the business unit that develops Zeiss Inspect is none other than the original
//! GOM GmbH, purchased by Zeiss in 2019. They've been developing inspection software since
//! the early 1990s.  I did this test to see how this problem has been thought about by people who
//! spent a lot of time thinking about metrology.

mod best_effort;
mod guarantee;

use crate::geom3::half_edge3::{HalfEdgeMesh3, NaAdaptor};
use crate::geom3::mesh::algorithms::normals::compute_face_normal;
use crate::{Point3, Result, Vector3};
use alum::{EditableTopology, FH, HH, Handle, HasIterators, HasTopology, Queue, VH};
use std::cell::Cell;

pub use best_effort::{BestEffortDecimator, BestEffortOpts};
pub use guarantee::{DecimateOpts, ErrorMethod, ToleranceVolumeDecimator};

type Quadric = alum::decimate::quadric::Quadric<NaAdaptor>;
type AlumPolyMesh = alum::PolyMeshT<3, NaAdaptor>;

/// Which quadric formulation to accumulate.
///
/// This mirrors `alum`'s own enum rather than re-exporting it, because that one derives nothing at
/// all and would make every options struct holding it un-`Copy` and un-printable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuadricKind {
    /// The classical Garland-Heckbert triangle quadric.
    #[default]
    Triangle,

    /// The probabilistic triangle quadric of Trettner and Kobbelt, which is more stable on noisy
    /// input because it models the uncertainty of each triangle rather than treating it as exact.
    Probabilistic,
}

/// Where the surviving vertex of a collapse is placed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Placement {
    /// The position which minimizes the combined quadric, falling back to the midpoint and then
    /// the two endpoints if that position is rejected.
    ///
    /// This gives the best surface fit and therefore the most decimation for a given tolerance,
    /// but it introduces positions which were never measured.
    #[default]
    Optimal,

    /// The midpoint of the collapsed edge, falling back to the two endpoints.
    Midpoint,

    /// One of the two original vertices, so every surviving point is one that was measured.
    ///
    /// Choose this when the point cloud must remain a subset of the input, at the cost of a coarser
    /// result for the same tolerance.
    Endpoint,
}

/// When decimation stops.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum DecimateTarget {
    /// Keep collapsing until no collapse is left which the gate permits.
    ///
    /// This is the metrology-shaped target: the caller states how much deviation is acceptable and
    /// the result is as coarse as that allows.
    #[default]
    ToTolerance,

    /// Stop once the mesh is down to this many faces, or sooner if the gate refuses everything.
    FaceCount(usize),

    /// Stop at this fraction of the starting face count, or sooner if the gate refuses everything.
    Ratio(f64),
}

impl DecimateTarget {
    /// The face count this target stops at, given the count the run started with.
    fn face_floor(&self, faces_before: usize) -> usize {
        match self {
            DecimateTarget::ToTolerance => 0,
            DecimateTarget::FaceCount(n) => *n,
            DecimateTarget::Ratio(r) => ((faces_before as f64) * r).round().max(1.0) as usize,
        }
    }

    fn validate(&self) -> Result<()> {
        if let DecimateTarget::Ratio(r) = self
            && (r.is_nan() || *r <= 0.0 || *r > 1.0)
        {
            return Err(format!("A decimation ratio must be in (0, 1], got {}", r).into());
        }

        Ok(())
    }
}

/// Counts of the work a decimation run did and why it refused what it refused.
///
/// These exist to make optimization decisions from measurement rather than from argument. The
/// veto counters in particular say which test is actually binding, which is what decides whether
/// an expensive one is earning its cost.
///
/// A veto is counted at the point it fires, and testing stops there, so the counts are of
/// *decisions* rather than of how many tests each candidate would have failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DecimateStats {
    /// Calls to the full evaluation of a candidate collapse.
    pub evaluations: usize,

    /// Calls to the acceptance test. One evaluation may make several, one per candidate position.
    pub acceptance_tests: usize,

    /// Collapses refused because a face would turn too far.
    pub veto_normal: usize,

    /// Collapses refused because a face would become a sliver.
    pub veto_aspect: usize,

    /// Collapses refused because a face would become degenerate.
    pub veto_degenerate: usize,

    /// Collapses refused by the gate's own error test.
    ///
    /// For the guaranteed path this is the error volume outgrowing the tolerance volume; for the
    /// best-effort path it is the quadric residual exceeding its cap.
    pub veto_error: usize,

    /// Acceptance tests where the chosen error method had no answer at all, so a different one had
    /// to judge the collapse.
    ///
    /// Not a refusal, which is why it is absent from [`vetoes`](DecimateStats::vetoes). The
    /// projected overlay needs both stars to survive one shared projection, and when they do not it
    /// says nothing rather than saying no, leaving the face-to-face map to decide. The count is what
    /// says whether that fallback machinery is servicing a real case or a hypothetical one.
    pub method_not_applicable: usize,
}

impl DecimateStats {
    /// The total number of refusals across every reason.
    pub fn vetoes(&self) -> usize {
        self.veto_normal + self.veto_aspect + self.veto_degenerate + self.veto_error
    }
}

impl std::fmt::Display for DecimateStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} evaluations, {} acceptance tests; vetoes: \
             normal {}, aspect {}, degenerate {}, error {}; not applicable {}",
            self.evaluations,
            self.acceptance_tests,
            self.veto_normal,
            self.veto_aspect,
            self.veto_degenerate,
            self.veto_error,
            self.method_not_applicable
        )
    }
}

/// What a decimation run did.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DecimateReport {
    /// How many edge collapses were performed.
    pub collapses: usize,

    /// The face count before the run.
    pub faces_before: usize,

    /// The face count after the run.
    pub faces_after: usize,

    /// How much work the run did, and which tests refused collapses.
    pub stats: DecimateStats,
}

impl DecimateReport {
    /// The fraction of the original faces which remain, or 1.0 if there were none to begin with.
    pub fn ratio(&self) -> f64 {
        if self.faces_before == 0 {
            return 1.0;
        }
        self.faces_after as f64 / self.faces_before as f64
    }
}

impl std::fmt::Display for DecimateReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} collapses, {} -> {} faces ({:.1}%)",
            self.collapses,
            self.faces_before,
            self.faces_after,
            self.ratio() * 100.0
        )
    }
}

/// The scale-invariant triangle quality used by `min_aspect`, peaking at 1 for an equilateral
/// triangle and falling to 0 for a degenerate one.
fn triangle_quality(p: &[Point3; 3]) -> f64 {
    let e = [p[1] - p[0], p[2] - p[1], p[0] - p[2]];
    let sum_sq: f64 = e.iter().map(|v| v.norm_squared()).sum();
    if sum_sq <= 0.0 {
        return 0.0;
    }

    let area = e[0].cross(&e[1]).norm() * 0.5;
    (4.0 * 3.0_f64.sqrt() * area) / sum_sq
}

/// One face of the edge star, in the form it has now.
///
/// The vertex handles are kept as they stand, with both endpoints of the collapsed edge still
/// distinct. That is what makes this a *correspondence* rather than two unrelated triangles: the
/// face a collapse leaves behind occupies the same slots with the endpoints replaced by the merged
/// vertex, so slot `k` of the old triangle maps to slot `k` of the new one.
///
/// Throughout this module the collapsed edge's endpoints are `v1` and `v2` and the merged position
/// is `p0`, following Guéziec's numbering. `v1` is the halfedge's tail and is deleted, `v2` is its
/// head and survives; `p0` has no handle until the collapse is applied, at which point it lives on
/// `v2`. See the `guarantee` module documentation for the full correspondence.
struct StarFace {
    /// The face itself.
    ///
    /// A collapse deletes the vanishing faces and leaves every other star face with its handle
    /// intact, so this stays valid across the collapse and is how a gate addresses per-face state
    /// it maintains. [`best_effort`] uses it to move the original vertices a face has absorbed onto
    /// whichever face replaces it.
    f: FH,

    /// The three vertices, as handles.
    v: [VH; 3],

    /// Where those three vertices are now.
    p: [Point3; 3],

    /// The face normal as it stands, for the structural vetoes. `None` where the face is already
    /// degenerate and so has nothing worth preserving.
    normal: Option<Vector3>,

    /// Whether this face disappears into the collapsed edge, which is the case when it uses both
    /// endpoints. Such a face has no corresponding new triangle and is handled separately.
    vanishing: bool,
}

/// The tests which are about the mesh remaining a usable mesh, rather than about where the surface
/// is, together with the work counters.
///
/// Both decimators embed one of these. The separation matters because these vetoes are the part
/// which is genuinely common: a flipped triangle or a sliver is equally unacceptable whether or not
/// the caller wanted a bound, and both failures can happen with the surface arbitrarily close to
/// where it belongs.
struct Structural {
    /// The largest angle, in radians, a face may be turned through by a collapse.
    max_normal_deviation: f64,

    /// The lowest triangle quality a collapse may produce.
    min_aspect: f64,

    /// Work counters. `Cell` because the costing path takes `&self`; decimation is single threaded,
    /// and an increment is a few nanoseconds against the cost of the tests around it, so these are
    /// always on rather than feature gated.
    stats: Cell<DecimateStats>,
}

impl Structural {
    fn new(max_normal_deviation: f64, min_aspect: f64) -> Self {
        Self {
            max_normal_deviation,
            min_aspect,
            stats: Cell::new(DecimateStats::default()),
        }
    }

    fn bump_evaluation(&self) {
        let mut s = self.stats.get();
        s.evaluations += 1;
        self.stats.set(s);
    }

    fn bump_acceptance(&self) {
        let mut s = self.stats.get();
        s.acceptance_tests += 1;
        self.stats.set(s);
    }

    fn bump_error_veto(&self) {
        let mut s = self.stats.get();
        s.veto_error += 1;
        self.stats.set(s);
    }

    fn bump_not_applicable(&self) {
        let mut s = self.stats.get();
        s.method_not_applicable += 1;
        self.stats.set(s);
    }

    /// Collect every face of the edge star, as it stands now.
    fn stars(
        &self,
        mesh: &AlumPolyMesh,
        v1: VH,
        v2: VH,
        points: &[Vector3],
    ) -> Option<Vec<StarFace>> {
        let mut out = Vec::with_capacity(16);

        for v in [v1, v2] {
            for f in mesh.vf_ccw_iter(v) {
                let mut vs = [VH::from(0u32); 3];
                let mut n = 0;
                for fv in mesh.fv_ccw_iter(f) {
                    if n >= 3 {
                        return None; // Not a triangle; this decimator only handles triangles.
                    }
                    vs[n] = fv;
                    n += 1;
                }
                if n != 3 {
                    return None;
                }

                // Faces around v2 which also touch v1 were already collected from v1's ring.
                if v == v2 && vs.contains(&v1) {
                    continue;
                }

                let p = [
                    Point3::from(points[vs[0].index() as usize]),
                    Point3::from(points[vs[1].index() as usize]),
                    Point3::from(points[vs[2].index() as usize]),
                ];

                out.push(StarFace {
                    f,
                    v: vs,
                    p,
                    normal: compute_face_normal(&p).map(|n| n.into_inner()),
                    vanishing: vs.contains(&v1) && vs.contains(&v2),
                });
            }
        }

        Some(out)
    }

    /// Build the star the collapse would leave behind, refusing it if any face is unusable.
    ///
    /// Returns the new position triple for every face of the star, in the same order, so slot `k`
    /// of `star[i]` corresponds to slot `k` of `out[i]`. Vanishing faces get an entry too, which is
    /// degenerate by construction and never used as a triangle.
    fn reshape(
        &self,
        star: &[StarFace],
        v1: VH,
        v2: VH,
        p0: Vector3,
        points: &[Vector3],
    ) -> Option<Vec<[Point3; 3]>> {
        let mut out = Vec::with_capacity(star.len());

        for face in star.iter() {
            let q = replacement(face, v1, v2, p0, points);
            out.push(q);

            if face.vanishing {
                continue;
            }
            let Some(normal) = face.normal else {
                continue;
            };

            let mut s = self.stats.get();
            let Some(new_normal) = compute_face_normal(&q) else {
                s.veto_degenerate += 1;
                self.stats.set(s);
                return None;
            };

            if new_normal.angle(&normal) > self.max_normal_deviation {
                s.veto_normal += 1;
                self.stats.set(s);
                return None;
            }

            if triangle_quality(&q) < self.min_aspect {
                s.veto_aspect += 1;
                self.stats.set(s);
                return None;
            }
        }

        Some(out)
    }
}

/// The triangle a face becomes once the collapse is applied, with both endpoints moved to `p0`.
fn replacement(face: &StarFace, v1: VH, v2: VH, p0: Vector3, points: &[Vector3]) -> [Point3; 3] {
    let mut q = [Point3::origin(); 3];
    for (slot, vh) in face.v.iter().enumerate() {
        q[slot] = if *vh == v1 || *vh == v2 {
            Point3::from(p0)
        } else {
            Point3::from(points[vh.index() as usize])
        };
    }
    q
}

/// What a gate needs to answer for the shared collapse driver to run it.
///
/// The two decimators differ only in this. Everything else, the queue, the staleness recheck, the
/// borrow juggling around `alum`'s status properties, and the structural vetoes, is the same either
/// way and lives in [`run_decimation`].
pub(crate) trait CollapseGate {
    /// What this gate records about a collapse it has approved.
    type Accepted;

    /// Whether this collapse is allowed, and if so what it would do.
    fn evaluate(&self, mesh: &AlumPolyMesh, h: HH) -> Option<Self::Accepted>;

    /// The queue key, which orders collapses but does not decide them.
    fn cost(accepted: &Self::Accepted) -> f64;

    /// Where the merged vertex goes.
    fn position(accepted: &Self::Accepted) -> Vector3;

    /// Record an approved collapse against the surviving vertex.
    fn commit(&mut self, v: VH, accepted: Self::Accepted) -> Result<()>;

    /// The work counters accumulated so far.
    fn stats(&self) -> DecimateStats;
}

/// Whether a collapse is structurally allowed, mirroring `alum`'s own private check.
///
/// The feature flag is how `lock_boundary` takes effect, and the two-incident-face requirement is
/// what stops a collapse dismantling an isolated fin.
fn is_collapse_legal(
    mesh: &AlumPolyMesh,
    h: HH,
    estatus: &alum::EPropBuf<alum::Status>,
    vstatus: &mut alum::VPropBuf<alum::Status>,
) -> bool {
    let v = h.tail(mesh);
    !vstatus[v].feature()
        && !estatus[h.edge()].feature()
        && mesh.vf_ccw_iter(v).take(2).count() == 2
        && mesh.check_edge_collapse(h, estatus, vstatus)
}

/// Find the cheapest acceptable collapse around a vertex, queue it, and return the halfedge.
///
/// `locked` is consulted before anything is costed. Every halfedge out of `v` has `v` as its tail,
/// and [`is_collapse_legal`] refuses a collapse whose tail carries the feature flag, so for a locked
/// vertex the whole loop below can only produce candidates that are thrown away later. Skipping is
/// exactly equivalent to running it, not an approximation of it.
///
/// This is worth a check rather than being left to the legality test because the gate is the
/// expensive part: on a flat grid with a locked outline it was 53% of every acceptance test the run
/// performed, and every one of those was on a collapse which could not happen.
fn queue_vertex<G: CollapseGate>(
    mesh: &AlumPolyMesh,
    v: VH,
    gate: &G,
    queue: &mut Queue<VH, f64>,
    locked: &[bool],
) -> Option<HH> {
    if locked.get(v.index() as usize).copied().unwrap_or(false) {
        queue.remove(v);
        return None;
    }

    let mut best: Option<(HH, f64)> = None;

    for h in mesh.voh_ccw_iter(v) {
        let Some(accepted) = gate.evaluate(mesh, h) else {
            continue;
        };
        let cost = G::cost(&accepted);
        if best.is_none_or(|(_, b)| cost < b) {
            best = Some((h, cost));
        }
    }

    match best {
        Some((h, cost)) => {
            queue.insert(v, cost);
            Some(h)
        }
        None => {
            queue.remove(v);
            None
        }
    }
}

/// Set the feature flag on every boundary vertex to `locked`, which is what the collapse legality
/// check reads.
///
/// Every run calls this with its own `lock_boundary`, rather than only the runs that lock. The
/// flag persists on the mesh, so a run that only ever set it would leave a boundary locked by one
/// run locked for every run after it, whatever the later options said. Each run therefore decides
/// the boundary's lock state outright. Feature flags on interior vertices are left alone, so a
/// caller who pinned interior vertices through [`HalfEdgeMesh3::as_alum_mut`] keeps them pinned.
///
/// Doing this once before the run means locking costs nothing per collapse. [`run_decimation`]
/// snapshots the flag afterward, so a vertex locked here is never even costed as a collapse tail.
fn set_boundary_locks(mesh: &mut HalfEdgeMesh3, locked: bool) -> Result<()> {
    let inner = mesh.as_alum_mut();
    let boundary: Vec<VH> = inner.vertices().filter(|v| v.is_boundary(inner)).collect();

    let mut vstatus = inner.vertex_status_prop();
    let mut vstatus = vstatus
        .try_borrow_mut()
        .map_err(|_| "Failed to borrow vertex status")?;

    for v in boundary {
        vstatus[v].set_feature(locked);
    }

    Ok(())
}

/// The collapse loop, shared by both decimators.
///
/// Pops the cheapest queued vertex, re-tests it, applies the collapse, and refreshes its one-ring.
/// The re-test is not redundant: the cost that queued a collapse was computed before some number of
/// neighboring collapses, any of which may have moved geometry or grown a value it depended on.
/// `alum`'s own driver cannot express this, because its `Decimater::before_collapse` hook has no way
/// to refuse, which is why this drives `alum`'s primitives directly.
///
/// # Arguments
///
/// * `mesh`: the half-edge mesh to decimate, in place
/// * `gate`: what decides each collapse
/// * `target`: when to stop
///
/// returns: `Result<DecimateReport>`
fn run_decimation<G: CollapseGate>(
    mesh: &mut HalfEdgeMesh3,
    gate: &mut G,
    target: DecimateTarget,
) -> Result<DecimateReport> {
    let faces_before = mesh.face_count();
    let target_faces = target.face_floor(faces_before);

    let alum_mesh = mesh.as_alum_mut();

    let mut vstatus = alum_mesh.vertex_status_prop();
    let mut hstatus = alum_mesh.halfedge_status_prop();
    let mut estatus = alum_mesh.edge_status_prop();
    let mut fstatus = alum_mesh.face_status_prop();

    let mut queue = Queue::<VH, f64>::new(alum_mesh.num_vertices());
    let mut targets: Vec<Option<HH>> = vec![None; alum_mesh.num_vertices()];
    let mut cache: Vec<HH> = Vec::new();
    let mut one_ring: Vec<VH> = Vec::new();

    // Which vertices can never be the tail of a legal collapse, which is what `lock_boundary`
    // produces. Snapshotted rather than read per call because the feature flag is set once before
    // the run and nothing during it touches it: `alum` never writes the flag, and
    // `lock_boundary_vertices` has already finished.
    let mut locked = vec![false; alum_mesh.num_vertices()];

    {
        let vs = vstatus
            .try_borrow()
            .map_err(|_| "Failed to borrow vertex status")?;
        for v in alum_mesh.vertices() {
            locked[v.index() as usize] = vs[v].feature();
        }
        for v in alum_mesh.vertices() {
            if !vs[v].deleted() {
                targets[v.index() as usize] = queue_vertex(alum_mesh, v, gate, &mut queue, &locked);
            }
        }
    }

    let mut collapses = 0usize;
    let mut faces = alum_mesh.num_faces();

    while let Some((v1, _)) = queue.pop() {
        if faces <= target_faces {
            break;
        }

        let Some(h) = targets[v1.index() as usize] else {
            continue;
        };

        {
            let es = estatus
                .try_borrow()
                .map_err(|_| "Failed to borrow edge status")?;
            let mut vs = vstatus
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow vertex status")?;
            if !is_collapse_legal(alum_mesh, h, &es, &mut vs) {
                continue;
            }
        }

        let Some(accepted) = gate.evaluate(alum_mesh, h) else {
            targets[v1.index() as usize] = queue_vertex(alum_mesh, v1, gate, &mut queue, &locked);
            continue;
        };

        one_ring.clear();
        one_ring.extend(alum_mesh.vv_ccw_iter(v1));

        let v2 = h.head(alum_mesh);
        let on_boundary = h.edge().is_boundary(alum_mesh);
        let position = G::position(&accepted);

        {
            let mut vs = vstatus
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow vertex status")?;
            let mut hs = hstatus
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow halfedge status")?;
            let mut es = estatus
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow edge status")?;
            let mut fs = fstatus
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow face status")?;
            alum_mesh.collapse_edge(h, &mut vs, &mut hs, &mut es, &mut fs, &mut cache);
        }

        {
            let mut points = alum_mesh.points();
            let mut points = points
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow points")?;
            points[v2] = position;
        }

        gate.commit(v2, accepted)?;

        collapses += 1;
        faces -= if on_boundary { 1 } else { 2 };

        for &v in one_ring.iter() {
            targets[v.index() as usize] = queue_vertex(alum_mesh, v, gate, &mut queue, &locked);
        }
    }

    queue.clear();
    let stats = gate.stats();

    // `num_faces` counts elements which are marked deleted but not yet reclaimed, so without this
    // the reported count would be the one we started with.
    mesh.garbage_collect()?;

    Ok(DecimateReport {
        collapses,
        faces_before,
        faces_after: mesh.face_count(),
        stats,
    })
}

/// The candidate positions for the merged vertex, in the order they should be tried.
///
/// Shared by both gates, because this is a rule about the mesh rather than about either error
/// scheme: which positions are legal at all, and which of the legal ones is preferred. A gate then
/// walks the list and takes the first its own test accepts.
///
/// `pinned` says the surviving vertex may not move, which is what `lock_boundary` produces. There
/// is then exactly one legal position for it whatever the placement rule would otherwise prefer.
/// Locking stops a boundary vertex being *deleted*, but on its own it does not stop an interior
/// vertex collapsing *into* it and dragging it off the boundary.
///
/// The quadric minimizer is an unpivoted solve with no rank check, so on a flat or nearly flat
/// region it divides by a value at or near zero and returns a non-finite result. That is filtered
/// here rather than left for the deviation test, which cannot compare a NaN.
fn candidate_positions(
    placement: Placement,
    pinned: bool,
    q: &Quadric,
    p1: Vector3,
    p2: Vector3,
) -> Vec<Vector3> {
    if pinned {
        return vec![p2];
    }

    let mid = (p1 + p2) * 0.5;

    match placement {
        Placement::Optimal => {
            let optimal = q.minimizer();
            let mut out = Vec::with_capacity(4);
            if optimal.iter().all(|c| c.is_finite()) {
                out.push(optimal);
            }
            out.push(mid);
            out.push(p2);
            out.push(p1);
            out
        }
        Placement::Midpoint => vec![mid, p2, p1],
        Placement::Endpoint => vec![p2, p1],
    }
}

/// Accumulate the quadric of a single vertex from the surrounding faces.
///
/// `alum`'s own `QuadricDecimater` keeps its per-vertex quadrics private, so they are accumulated
/// here the same way it does.
fn quadric_of_vertex(mesh: &AlumPolyMesh, v: VH, kind: QuadricKind) -> Result<Quadric> {
    let points = mesh.points();
    let points = points
        .try_borrow()
        .map_err(|_| "Failed to borrow points while seeding quadrics")?;

    let mut total = Quadric::default();

    for f in mesh.vf_ccw_iter(v) {
        for vs in mesh.triangulated_face_vertices(f) {
            let (a, b, c) = (points[vs[0]], points[vs[1]], points[vs[2]]);
            total += match kind {
                QuadricKind::Triangle => Quadric::triangle_quadric(a, b, c),
                QuadricKind::Probabilistic => {
                    // The same edge-length-derived standard deviation `alum` uses, which keeps
                    // the noise model scaled to the local tessellation.
                    let mean_edge = ((b - a).norm() + (c - b).norm() + (a - c).norm()) / 3.0;
                    Quadric::probabilistic_triangle_quadric(a, b, c, 0.1 * mean_edge)
                }
            };
        }
    }

    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Mesh3;
    use crate::geom3::half_edge3::difference::ShapeDifference;

    const TOL: f64 = 0.05;

    /// A face count target stops where it was told to, and the result still holds the bound.
    #[test]
    fn a_face_count_target_is_respected() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        let target = 30_000;
        let report = he
            .decimate_guaranteed(&DecimateOpts::to_face_count(TOL, target))
            .unwrap();

        assert!(
            report.faces_after <= target + 2,
            "expected to stop near {}, got {}",
            target,
            report
        );

        let out = Mesh3::try_from(&he).unwrap();
        ShapeDifference::new(&original, &out)
            .unwrap()
            .assert_within(TOL);
    }

    #[test]
    fn a_ratio_target_is_respected() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        let report = he
            .decimate_guaranteed(&DecimateOpts::to_ratio(TOL * 2.0, 0.5))
            .unwrap();

        assert!(report.ratio() <= 0.51, "got {}", report);
        assert!(report.ratio() > 0.3, "stopped far too early: {}", report);
    }

    /// The gate still binds when a count target asks for more than it allows.
    #[test]
    fn the_gate_wins_over_an_unreachable_count() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        let report = he
            .decimate_guaranteed(&DecimateOpts::to_face_count(0.0002, 100))
            .unwrap();

        assert!(
            report.faces_after > 100,
            "the tolerance should have stopped this well short of 100 faces, got {}",
            report
        );

        let out = Mesh3::try_from(&he).unwrap();
        ShapeDifference::new(&original, &out)
            .unwrap()
            .assert_within(0.0002);
    }
}
