//! A direct check of Guéziec's two invariants, and the hand-built collapses it is run on.
//!
//! # Why this exists
//!
//! Every error method in the parent module is an _argument_ that a radius is large enough: build a
//! correspondence, collect constraints at the corners of its pieces, and appeal to concavity for
//! everything in between. The failure mode of such an argument is not an exception or a wrong
//! answer, it is a region of the surface whose constraints were never collected. Nothing throws.
//! The radius simply comes out too small, and the bound is not a bound.
//!
//! That failure has happened three times in this module's history: twice on Guéziec's Subdivision
//! Method, and once on the projected overlay applied to a star whose outline moves. In all three
//! cases the mesh-level tests kept passing for a while, because the violation has to be both real
//! and large enough to show up in a dense sampling of two whole meshes before anything notices.
//!
//! This module asks the invariants directly instead, on one collapse at a time:
//!
//! - (𝔸) every ball of the old surface lies inside some ball of the new one;
//! - (𝔹) every ball of the new surface contains some ball of the old one.
//!
//! It consults nothing about how the radius was arrived at, which is the whole point. A check which
//! samples the sites a method constrains agrees with that method whether the method is right.
//!
//! # Why the search is exact rather than sampled
//!
//! Both invariants are existential: they ask whether _some_ ball on the other surface works. So
//! each sample needs an optimum over the other surface, and getting that optimum wrong in the safe
//! direction produces false alarms rather than crickets.
//!
//! An earlier attempt sampled the other surface, and separately took the closest point on each
//! triangle, and both were wrong. For a sample `x` and a target triangle, the quantity to maximize
//! is...
//!
//! ```text
//! e(y) - |x - y|
//! ```
//!
//! ...which is affine minus a norm of an affine function, hence concave in `y`. Its maximizer is
//! the closest point only when `x` lies in the triangle's plane. Off the plane it sits away from
//! the closest point, in the direction the radius field grows, by an amount proportional to how
//! far off the plane `x` is. Concavity is also what makes the optimum computable rather than
//! something we'd have to search for. It's either the single interior stationary point, which has
//! a closed form, or it's on an edge, where the problem is one-dimensional and still concave.
//!
//! Getting this exact matters because the binding constraint is tight by construction. A method
//! returns the smallest radius satisfying its constraints, so at the corner that bound it the
//! inequality holds with equality, and a search that lands near the witness rather than on it
//! reports a violation which is not there.

use super::boundary::Outline;
use super::*;
use crate::common::barycentric::{barycentric, barycentric_grid, barycentric_point};
use crate::geom3::half_edge3::decimate::replacement;
use crate::{Point3, Vector3};
use alum::FH;

// The invariant check ───────────────────────────────────────────────────────

/// One triangle of a surface, with the error radius at each corner.
///
/// This is all either invariant needs to know about a surface: the balls are the barycentric blend
/// of these radii swept over the triangle.
#[derive(Debug, Clone, Copy)]
struct Patch {
    p: [Point3; 3],
    e: [f64; 3],
}

impl Patch {
    /// The radius at barycentric weights `w`.
    fn radius(&self, w: [f64; 3]) -> f64 {
        self.e[0] * w[0] + self.e[1] * w[1] + self.e[2] * w[2]
    }

    /// The radius at a point, which for a point off the plane is the radius at its projection.
    fn radius_at(&self, z: &Point3) -> f64 {
        self.radius(barycentric(&self.p[0], &self.p[1], &self.p[2], z))
    }

    /// The same patch with every radius negated, which is how the reaching direction is expressed
    /// as a covering one. See [`best_cover_on`].
    fn negated(&self) -> Self {
        Self {
            p: self.p,
            e: [-self.e[0], -self.e[1], -self.e[2]],
        }
    }
}

/// The largest `patch.radius(y) - |q - y|` over the patch.
///
/// The quantity is exactly "how much room this patch's balls have to spare for a ball centred at
/// `q`": the patch covers a ball of radius `r` at `q` precisely when this is at least `r`.
///
/// Both invariants reduce to it. (𝔸) asks whether some new patch covers the old ball at `x`, which
/// is `best_cover_on(new, x) >= e_old(x)`. (𝔹) asks for the smallest `|y - x| + e_old(x)` over the
/// old surface, and negating the radii turns that minimum into this maximum, so there is one piece
/// of delicate code rather than two.
///
/// Concave in `y`, so the maximum is either at the interior stationary point or on an edge, and both
/// are found rather than sampled.
fn best_cover_on(patch: &Patch, q: &Point3) -> f64 {
    let (a, b, c) = (patch.p[0], patch.p[1], patch.p[2]);
    let Some(n) = (b - a).cross(&(c - a)).try_normalize(1.0e-14) else {
        return f64::NEG_INFINITY; // A degenerate triangle carries no balls worth speaking of.
    };

    let f = |z: &Point3| patch.radius_at(z) - (q - z).norm();

    // The three edges. Each is one dimensional and still concave, so a ternary search is exact up to
    // the interval it closes on, and the corners are picked up as its endpoints.
    let mut best = f64::NEG_INFINITY;
    for (u, v) in [(a, b), (b, c), (c, a)] {
        best = best.max(ternary_max(&u, &v, &f));
    }

    // The interior stationary point. Writing `r = y - q'` for the in-plane offset from `q`'s
    // projection and `g` for the in-plane gradient of the radius field, the derivative vanishes
    // where `r / sqrt(|r|^2 + d^2) = g`, which puts `r` along `g` at distance
    // `|g| d / sqrt(1 - |g|^2)`. A field steeper than one has no stationary point at all and its
    // maximum is on an edge, which the loop above already covered.
    //
    // A flat field is not a special case to skip but the commonest one there is: the step is then
    // zero and the stationary point is the perpendicular foot, which is to say the closest point on
    // the triangle. Treating a zero gradient as "no stationary point" is what a first version of
    // this did, and it made the check disagree with contraction, whose soundness is one line.
    let d = n.dot(&(q - a));
    let qp = q - n * d;

    if let Some(g) = radius_gradient(patch, &n) {
        let gm = g.norm();
        if gm < 1.0 {
            let step = gm * d.abs() / (1.0 - gm * gm).sqrt();
            let z = if gm > 0.0 { qp + g * (step / gm) } else { qp };
            let w = barycentric(&a, &b, &c, &z);
            if w.iter().all(|x| *x >= 0.0) {
                best = best.max(f(&z));
            }
        }
    }

    best
}

/// The in-plane gradient of a patch's radius field, or `None` for a degenerate triangle.
///
/// Solves for the vector whose dot product with each edge is that edge's rise in radius, and which
/// has no component out of the plane.
fn radius_gradient(patch: &Patch, n: &Vector3) -> Option<Vector3> {
    let u = patch.p[1] - patch.p[0];
    let v = patch.p[2] - patch.p[0];

    let m = crate::na::Matrix3::from_rows(&[u.transpose(), v.transpose(), n.transpose()]);
    let rhs = Vector3::new(patch.e[1] - patch.e[0], patch.e[2] - patch.e[0], 0.0);
    m.lu().solve(&rhs)
}

/// The maximum of a concave function along a segment, endpoints included.
fn ternary_max(u: &Point3, v: &Point3, f: &impl Fn(&Point3) -> f64) -> f64 {
    let at = |t: f64| Point3::from(u.coords * (1.0 - t) + v.coords * t);

    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..60 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if f(&at(m1)) < f(&at(m2)) {
            lo = m1;
        } else {
            hi = m2;
        }
    }

    f(&at(0.5 * (lo + hi))).max(f(u)).max(f(v))
}

/// How badly each invariant is violated, in length units, with zero or less meaning it holds.
///
/// The first number is (𝔸), how much larger some new ball would have to be to swallow the worst old
/// ball. The second is (𝔹), how much larger some new ball would have to be to still reach the old
/// surface. Both are reported rather than a pass or fail, because the size of a violation is what
/// says whether it is a real undercharge or an artifact of the sampling.
fn invariant_gaps(old: &[Patch], new: &[Patch], spacing: f64) -> (f64, f64) {
    let mut cover = 0.0f64;
    let mut reach = 0.0f64;

    for patch in old.iter() {
        for w in barycentric_grid(&patch.p[0], &patch.p[1], &patch.p[2], spacing) {
            let x = barycentric_point(&patch.p[0], &patch.p[1], &patch.p[2], w);
            let have = new
                .iter()
                .map(|n| best_cover_on(n, &x))
                .fold(f64::NEG_INFINITY, f64::max);
            cover = cover.max(patch.radius(w) - have);
        }
    }

    for patch in new.iter() {
        for w in barycentric_grid(&patch.p[0], &patch.p[1], &patch.p[2], spacing) {
            let y = barycentric_point(&patch.p[0], &patch.p[1], &patch.p[2], w);
            let need = old
                .iter()
                .map(|o| -best_cover_on(&o.negated(), &y))
                .fold(f64::INFINITY, f64::min);
            reach = reach.max(need - patch.radius(w));
        }
    }

    (cover, reach)
}

// The fixtures ──────────────────────────────────────────────────────────────

/// A 4 by 3 grid of vertices triangulated as two rows of quads, so the middle row is interior and
/// everything around it is on the outline.
///
/// Built by hand rather than taken from a mesh because the point of these tests is which
/// configuration each edge is, and that has to be legible from the fixture to be worth asserting.
/// Every kind the chain extraction distinguishes is present in this one strip.
///
/// ```text
///   8    9   10   11      y=2, all on the outline
///   4    5    6    7      y=1, 5 and 6 interior, 4 and 7 on the outline
///   0    1    2    3      y=0, all on the outline
/// ```
///
/// How the strip is shaped.
#[derive(Debug, Clone, Copy, Default)]
struct Strip {
    /// Lifts the strip into a saddle, curved in both directions at once. A flat strip is the
    /// cleaner fixture to reason about, since every projected distance is the real one, but it
    /// leaves the invariant check's off-plane machinery untested.
    warp: f64,

    /// Pushes the two middle vertices of the bottom row up into the surface, bending that part of
    /// the outline instead of leaving it straight.
    ///
    /// This is not decoration. A **straight** outline cannot exhibit the overlay's undercharge at
    /// all, because the overlay's own location checks catch it first: with the merged vertex above
    /// a straight run, the run's interior vertices necessarily fall below the chain it becomes and
    /// so land outside the new star, and with it below, the merged vertex lands outside the old
    /// one. Either way the overlay declines for a reason that has nothing to do with the premise.
    /// Only a bent run leaves room for every location to succeed while the regions still differ,
    /// which is why the flat grids never showed the bug and the bunny did.
    notch: f64,
}

impl Strip {
    fn flat() -> Self {
        Self::default()
    }

    fn warped() -> Self {
        Self {
            warp: 0.25,
            notch: 0.0,
        }
    }

    /// The shape which can reach the undercharge, with the outline bent by four tenths of an edge.
    fn bent() -> Self {
        Self {
            warp: 0.0,
            notch: 0.4,
        }
    }
}

/// See [`Strip`] for what the two shape parameters do.
fn strip(shape: Strip) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let mut points = Vec::new();
    for j in 0..3 {
        for i in 0..4 {
            let x = i as f64;
            let y = j as f64
                + if j == 0 && (i == 1 || i == 2) {
                    shape.notch
                } else {
                    0.0
                };
            points.push(Point3::new(x, y, shape.warp * (x - 1.5) * (y - 1.0)));
        }
    }

    let mut faces = Vec::new();
    for j in 0..2u32 {
        for i in 0..3u32 {
            let (a, b) = (j * 4 + i, (j + 1) * 4 + i);
            faces.push([a, a + 1, b + 1]);
            faces.push([a, b + 1, b]);
        }
    }

    (points, faces)
}

/// Everything an error method needs for one candidate collapse, assembled by hand.
///
/// This owns its data, where [`Collapse`] borrows, so a fixture outlives the contexts it lends out
/// and one fixture can be judged under more than one outline.
struct Fixture {
    star: Vec<StarFace>,
    new: Vec<[Point3; 3]>,
    points: Vec<Vector3>,
    e: Vec<f64>,
    v1: VH,
    v2: VH,
    p0: Vector3,
}

impl Fixture {
    /// Collapse `v1` into `v2` on the strip, merging them at `p0`, with every vertex starting at the
    /// same error radius.
    ///
    /// A uniform starting radius rather than zero because zero would let a method get the blend
    /// weights wrong and still look right.
    fn new(shape: Strip, v1: u32, v2: u32, p0: Point3, radius: f64) -> Self {
        let (points, faces) = strip(shape);
        let (hv1, hv2) = (VH::from(v1), VH::from(v2));

        // The same set `Structural::stars` collects: every face touching either endpoint.
        let star: Vec<StarFace> = faces
            .iter()
            .filter(|f| f.contains(&v1) || f.contains(&v2))
            .enumerate()
            .map(|(i, f)| StarFace {
                f: FH::from(i as u32),
                v: [VH::from(f[0]), VH::from(f[1]), VH::from(f[2])],
                p: [
                    points[f[0] as usize],
                    points[f[1] as usize],
                    points[f[2] as usize],
                ],
                normal: None,
                vanishing: f.contains(&v1) && f.contains(&v2),
            })
            .collect();

        let coords: Vec<Vector3> = points.iter().map(|q| q.coords).collect();
        let new = star
            .iter()
            .map(|f| replacement(f, hv1, hv2, p0.coords, &coords))
            .collect();

        Self {
            star,
            new,
            e: vec![radius; points.len()],
            points: coords,
            v1: hv1,
            v2: hv2,
            p0: p0.coords,
        }
    }

    fn outline(&self) -> Outline {
        boundary_chain(&self.star, self.v1, self.v2)
    }

    /// The old surface, carrying the radii the collapse starts with.
    fn old_patches(&self) -> Vec<Patch> {
        self.star
            .iter()
            .map(|f| Patch {
                p: f.p,
                e: [0, 1, 2].map(|k| self.e[f.v[k].index() as usize]),
            })
            .collect()
    }

    /// The surface the collapse leaves, with `radius` at the merged vertex.
    ///
    /// Faces which vanish into the collapsed edge are dropped, since they have no successor; the old
    /// side keeps them, because leaving them out would put a hole in the reference surface exactly
    /// where the merged vertex lands.
    fn new_patches(&self, radius: f64) -> Vec<Patch> {
        self.star
            .iter()
            .zip(self.new.iter())
            .filter(|(f, _)| !f.vanishing)
            .map(|(f, q)| Patch {
                p: *q,
                e: [0, 1, 2].map(|k| {
                    let h = f.v[k];
                    if h == self.v1 || h == self.v2 {
                        radius
                    } else {
                        self.e[h.index() as usize]
                    }
                }),
            })
            .collect()
    }

    /// The context an error method is handed, telling it the truth about the outline.
    fn collapse(&self) -> Collapse<'_> {
        self.at(self.outline())
    }

    /// The same, with the outline overridden, which is how a method can be asked to judge a
    /// collapse on a premise that does not hold.
    fn at(&self, outline: Outline) -> Collapse<'_> {
        Collapse {
            star: &self.star,
            new: &self.new,
            points: &self.points,
            e: &self.e,
            v1: self.v1,
            v2: self.v2,
            p0: self.p0,
            outline,
        }
    }

    /// Both invariant gaps for a radius some method produced.
    fn gaps(&self, radius: f64) -> (f64, f64) {
        invariant_gaps(&self.old_patches(), &self.new_patches(radius), 0.15)
    }
}

/// How much slack the invariant check is allowed before it calls a sample violated.
///
/// The binding constraint is tight by construction, so this absorbs nothing but arithmetic. The
/// fixture is unit scale, which is why an absolute figure is enough.
const SLACK: f64 = 1.0e-9;

// The star's boundary chain ─────────────────────────────────────────────────

fn outline_of(v1: u32, v2: u32) -> Outline {
    Fixture::new(Strip::flat(), v1, v2, Point3::origin(), 0.0).outline()
}

/// An edge with both endpoints off the outline has nothing a collapse can move.
///
/// This is every collapse on a watertight mesh, and it is the case the projected overlay was built
/// for.
#[test]
fn an_interior_edge_has_a_closed_outline() {
    assert!(matches!(outline_of(5, 6), Outline::Closed));
}

/// Gueziec's Type II edge: one endpoint on the outline, two incident triangles.
///
/// The run is the two outline edges at that endpoint, so three vertices with the endpoint in the
/// middle. The anchors are the vertices the collapse leaves alone.
#[test]
fn a_type_two_edge_gives_a_three_vertex_chain() {
    let Outline::Open(chain) = outline_of(1, 5) else {
        panic!("expected an open outline for a Type II edge");
    };

    let got: Vec<u32> = chain.iter().map(|(v, _)| v.index()).collect();
    assert_eq!(got, vec![0, 1, 2], "the run should be 0 - 1 - 2");
}

/// Gueziec's Type I edge: the collapsed edge is itself on the outline, carrying one triangle.
///
/// The run picks up that edge as well as the outline edge at each endpoint, so four vertices with
/// both endpoints in the middle. This is the case with no method in the paper.
#[test]
fn a_type_one_edge_gives_a_four_vertex_chain() {
    let Outline::Open(chain) = outline_of(1, 2) else {
        panic!("expected an open outline for a Type I edge");
    };

    let got: Vec<u32> = chain.iter().map(|(v, _)| v.index()).collect();
    assert_eq!(got, vec![0, 1, 2, 3], "the run should be 0 - 1 - 2 - 3");
}

/// Gueziec's Type III edge: both endpoints on the outline, but the edge between them is not.
///
/// He refuses these outright, since collapsing one joins two parts of the outline and leaves a
/// singular vertex on the result. It is refused here too, and for a reason the extraction can see
/// without being told which kind of edge it has: the run through the endpoints is not a single one
/// with only merged vertices inside it, so there is no chain of the required shape for the collapse
/// to produce.
#[test]
fn a_type_three_edge_is_refused() {
    assert!(matches!(outline_of(4, 9), Outline::Unsupported));
}

/// A chain does not move when every vertex it is losing is already where the merged vertex lands,
/// which is what a pinned boundary collapse does.
#[test]
fn a_chain_only_moves_when_a_vertex_it_carries_does() {
    let Outline::Open(chain) = outline_of(1, 5) else {
        panic!("expected an open outline");
    };

    // Vertex 1 sits at (1, 0, 0) and is the only one the collapse relocates, so the run holds still
    // exactly when the merged vertex is pinned there.
    assert!(!chain.moves_under(Vector3::new(1.0, 0.0, 0.0)));
    assert!(chain.moves_under(Vector3::new(1.0, 0.1, 0.0)));

    // Sliding along the run still counts, even though the outline sweeps out no area, because the
    // anchors do not move and the two segments of the run therefore change length.
    assert!(chain.moves_under(Vector3::new(1.5, 0.0, 0.0)));
}

// Calibrating the check ─────────────────────────────────────────────────────

/// The check must agree with the method whose soundness needs no argument.
///
/// Contraction widens the ball at the merged vertex until it swallows the balls at both endpoints,
/// which is the whole of both invariants in one line. If the check disagreed with it, the check
/// would be what is wrong.
#[test]
fn the_invariant_check_agrees_with_contraction() {
    for shape in [Strip::flat(), Strip::warped()] {
        for (v1, v2) in [(5, 6), (1, 5), (1, 2)] {
            let c = Fixture::new(shape, v1, v2, Point3::new(1.4, 0.45, 0.1), 0.03);
            let ErrorBound::Bound(r) = Contraction.bound(&c.collapse()) else {
                panic!("contraction always produces a bound");
            };

            let (cover, reach) = c.gaps(r);
            assert!(
                cover <= SLACK && reach <= SLACK,
                "contraction flagged on shape={shape:?} edge ({v1},{v2}): cover {cover:e}, \
                 reach {reach:e}"
            );
        }
    }
}

/// The face-to-face map holds both invariants, including on a moving outline, which is why it is
/// what the overlay falls back to.
///
/// Its pieces are one triangle onto one triangle, so it never assumed the two stars cover the same
/// region and a boundary collapse costs it nothing in validity. This says so from the outside.
#[test]
fn the_invariant_check_agrees_with_the_face_map() {
    for shape in [Strip::flat(), Strip::warped()] {
        for (v1, v2) in [(5, 6), (1, 5), (1, 2)] {
            let c = Fixture::new(shape, v1, v2, Point3::new(1.4, 0.45, 0.1), 0.03);
            let ErrorBound::Bound(r) = AffineFaceMap.bound(&c.collapse()) else {
                panic!("the face map refused shape={shape:?} edge ({v1},{v2})");
            };

            let (cover, reach) = c.gaps(r);
            assert!(
                cover <= SLACK && reach <= SLACK,
                "the face map flagged on shape={shape:?} edge ({v1},{v2}): cover {cover:e}, \
                 reach {reach:e}"
            );
        }
    }
}

/// The overlay holds both invariants where its premise does, which is an interior collapse.
#[test]
fn the_invariant_check_agrees_with_the_overlay_on_a_closed_outline() {
    for shape in [Strip::flat(), Strip::warped()] {
        let c = Fixture::new(shape, 5, 6, Point3::new(1.4, 1.15, 0.1), 0.03);
        assert!(matches!(c.outline(), Outline::Closed));

        let ErrorBound::Bound(r) = ProjectedOverlay.bound(&c.at(Outline::Closed)) else {
            panic!("the overlay had nothing to say about an interior collapse at shape={shape:?}");
        };

        let (cover, reach) = c.gaps(r);
        assert!(
            cover <= SLACK && reach <= SLACK,
            "the overlay flagged on an interior collapse at shape={shape:?}: cover {cover:e}, \
             reach {reach:e}"
        );
    }
}

/// This runs the overlay on a Type I collapse while telling it the outline holds still, which is
/// what it did before the premise was made explicit. The two stars then cover different regions,
/// the part of each with no counterpart in the other is charged for by no constraint, and the
/// radius comes out far too small.
///
/// The strip is bent rather than straight, and it has to be. See [`Strip::notch`]: a straight run
/// cannot reach this at all, because the overlay's own location checks refuse the collapse first for
/// an unrelated reason. That is why the flat grids never showed the bug and the bunny did.
///
/// It stays flat in `z`, though, so every point of it is coplanar. The overlay's correspondence then
/// moves nothing at all and it returns essentially the radius it started with, while the outline has
/// in fact swung a tenth of an edge length. There is no ambiguity about whether the gap is real.
#[test]
fn the_invariant_check_catches_the_overlay_on_a_moving_outline() {
    let c = Fixture::new(Strip::bent(), 1, 2, Point3::new(1.5, 0.5, 0.0), 0.03);

    assert!(
        matches!(c.outline(), Outline::Open(_)),
        "edge (1,2) should be a Type I boundary edge"
    );

    let ErrorBound::Bound(r) = ProjectedOverlay.bound(&c.at(Outline::Closed)) else {
        panic!("expected the overlay to produce a bound when told the outline holds still");
    };

    // The overlay charges nothing at all here. Every pair it constrains has the same projection and
    // the strip is coplanar, so every distance it looks at is zero and it hands back the radius the
    // collapse started with.
    let (cover, reach) = c.gaps(r);
    assert!(
        cover.max(reach) > r,
        "the check missed an undercharge it should not have: radius {r:e}, cover {cover:e}, \
         reach {reach:e}"
    );

    // What the radius should have been, from the method the check has already certified on this
    // same collapse. The gap is not a rounding matter.
    let ErrorBound::Bound(sound) = AffineFaceMap.bound(&c.collapse()) else {
        panic!("the face map should judge this collapse");
    };
    assert!(
        sound > 3.0 * r,
        "expected the sound radius to dwarf the overlay's, got {sound:e} against {r:e}"
    );

    // Told the truth about the outline, the overlay covers the slivers as well and the check is
    // satisfied. This is the whole of the Type I correspondence in one assertion: the same method on
    // the same collapse, sound once it stops assuming the two stars span the same region.
    let ErrorBound::Bound(fixed) = ProjectedOverlay.bound(&c.collapse()) else {
        panic!("expected the overlay to judge a moving outline rather than decline it");
    };

    let (cover, reach) = c.gaps(fixed);
    assert!(
        cover <= SLACK && reach <= SLACK,
        "the overlay still undercharges a moving outline: radius {fixed:e}, cover {cover:e}, \
         reach {reach:e}"
    );

    // And it is worth having, rather than merely correct: the slivers charge for what actually moved
    // instead of for the whole star, so it stays well inside what the face-to-face map demands.
    assert!(
        fixed < sound,
        "expected the overlay to beat the face map, got {fixed:e} against {sound:e}"
    );
}
