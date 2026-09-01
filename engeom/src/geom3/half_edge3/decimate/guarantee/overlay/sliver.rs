//! Covering the part of a star which has no counterpart on the other one.
//!
//! Private to the overlay, and specific to it. When a collapse moves the mesh outline the two stars
//! stop covering the same region, the bijection the overlay rests on stops existing, and what is
//! left over is a sliver of each star that no constraint would otherwise charge for. A method whose
//! pieces are one triangle onto one triangle never has this problem and would have no use for any
//! of this.
//!
//! The sliver is carried onto the chain the collapse leaves behind by the nearest-point projection,
//! which is piecewise affine over the chain's own Voronoi regions. Those pieces, and the corners
//! which are all that has to be evaluated on them, are what this module produces.

use super::segment_cross;
use crate::common::barycentric::barycentric_within2;
use crate::geom2::Segment2;
use crate::geom2::hull::convex_hull_2d;
use crate::{Point2, Point3, Vector2};

/// Whether a projected point falls anywhere in the region a set of projected triangles covers.
fn inside_any(q: &Point2, faces: &[[Point2; 3]]) -> bool {
    faces
        .iter()
        .any(|t| barycentric_within2(&t[0], &t[1], &t[2], q).is_some())
}

/// Whether a projected point sits on the chain, which is where the two regions part company.
fn on_chain(q: &Point2, chain: &[Point2]) -> bool {
    chain.windows(2).any(|s| {
        let seg = Segment2::new_unchecked(s[0], s[1]);
        (seg.closest_point(q) - q).norm() < 1.0e-9
    })
}

/// The corners of the part of a projected triangle which falls outside the other star's region.
///
/// # Why corners are all that is wanted
///
/// The slack at a point is affine minus convex minus affine under any affine map, hence concave, and
/// **a concave function on a polygon attains its minimum at a vertex even when the polygon is not
/// convex**: the polygon sits inside its own convex hull, whose extreme points are all vertices of
/// the polygon. So the region itself never has to be constructed, only its corner set, which turns
/// what would be polygon clipping into three cheap enumerations.
///
/// # Where the corners are
///
/// The region is bounded by parts of the triangle's own edges and parts of the other star's outline.
/// Only the *chain* part of that outline can bound it: both stars share their link, which lies on
/// the boundary of each, so a triangle of one can touch it but never cross it. That leaves
///
/// - corners of the triangle which are not covered by the other star,
/// - where the triangle's edges meet the chain,
/// - the chain's own turning points, where it bends inside the triangle.
///
/// A corner sitting *on* the chain is kept as well, even though the containment test calls it
/// covered. It is a genuine corner of the region, and it costs nothing: its own foot on the chain is
/// itself, so the constraint it produces is that a radius must be at least itself.
///
/// An empty result means the triangle is covered and has no sliver.
pub(super) fn sliver_corners(
    tri: &[Point2; 3],
    other: &[[Point2; 3]],
    chain: &[Point2],
) -> Vec<Point2> {
    let mut out: Vec<Point2> = Vec::with_capacity(6);

    let push = |out: &mut Vec<Point2>, q: Point2| {
        if !out.iter().any(|z| (z - q).norm() < 1.0e-12) {
            out.push(q);
        }
    };

    for q in tri.iter() {
        if !inside_any(q, other) || on_chain(q, chain) {
            push(&mut out, *q);
        }
    }

    for k in 0..3 {
        let (u, v) = (tri[k], tri[(k + 1) % 3]);
        for s in chain.windows(2) {
            if let Some((t, _)) = segment_cross(&u, &v, &s[0], &s[1]) {
                push(&mut out, Point2::from(u.coords * (1.0 - t) + v.coords * t));
            }
        }
    }

    for q in chain.iter().skip(1).take(chain.len().saturating_sub(2)) {
        if barycentric_within2(&tri[0], &tri[1], &tri[2], q).is_some() {
            push(&mut out, *q);
        }
    }

    out
}

/// Which part of a chain a piece of a sliver is nearest to, and so what carries it there.
#[derive(Debug, Clone, Copy)]
pub(super) enum ChainFeature {
    /// The perpendicular projection onto a segment, which is affine and lands inside it.
    Segment(usize),

    /// The constant map onto a vertex, for the wedge of points whose nearest chain point is that
    /// corner. Constant is as affine as anything, and a single point is inside every triangle having
    /// it as a corner.
    Vertex(usize),
}

/// Split a sliver into the pieces on which "the nearest point of the chain" is one affine map.
///
/// # Why this decomposition and not another
///
/// A sliver has to be carried onto the chain by an affine map, and the concavity argument only holds
/// when *one* map covers the whole piece. The obvious candidates are the projection onto a chain
/// segment, which is only admissible while every foot lands inside that segment, and the constant
/// map onto a chain vertex, which is always admissible but charges the distance to that vertex.
///
/// Letting whole maps compete for a whole sliver does not work: a sliver straddling a bend in the
/// chain has no valid segment projection at all and falls back to a constant map, which then charges
/// it the distance to a corner rather than to the chain. Measured on a bent strip, that came out at
/// two and a half times what the face-to-face map asks for, on a collapse the face map judges
/// perfectly well.
///
/// Splitting by _nearest feature_ instead makes the map the nearest-point projection onto the
/// chain.  Its pieces are the polyline's own Voronoi regions: a slab either side of each segment,
/// and a wedge at each vertex where consecutive slabs leave a gap. Inside a slab a segment
/// projection is admissible by construction, and inside a wedge the vertex *is* the nearest chain
/// point, so the constant map charges the true distance rather than an inflated one. The regions
/// cover the plane, so the pieces cover the sliver.
///
/// Everything here is in the projection. The map is still affine in space, because lifting a chain
/// parameter back onto the chain in space is affine along each segment.
pub(super) fn chain_pieces(
    corners: &[Point2],
    chain: &[Point2],
) -> Vec<(Vec<Point2>, ChainFeature)> {
    // The sliver's own convex hull, which contains it and has a subset of its corners as vertices.
    // Charging on the hull is no looser: a concave function over it is minimised at one of those
    // vertices, which is a corner of the sliver either way. What the hull buys is edges to clip.
    let hull: Vec<Point2> = if corners.len() < 3 {
        corners.to_vec()
    } else {
        convex_hull_2d(corners)
            .into_iter()
            .map(|i| corners[i])
            .collect()
    };

    let dir = |k: usize| chain[k + 1] - chain[k];
    let last = chain.len() - 1;
    let mut out = Vec::with_capacity(2 * chain.len());

    // A half plane as the points where `(q - c) . n` is non-negative.
    let mut cuts: Vec<(Point2, Vector2)> = Vec::with_capacity(2);
    let clip = |cuts: &[(Point2, Vector2)]| -> Vec<Point2> {
        cuts.iter()
            .fold(hull.clone(), |poly, (c, n)| clip_half(&poly, c, n))
    };

    for k in 0..chain.len() {
        // The slab either side of segment `k`, between the perpendiculars at its two ends.
        if k < last {
            cuts.clear();
            cuts.push((chain[k], dir(k)));
            cuts.push((chain[k + 1], -dir(k)));

            let piece = clip(&cuts);
            if !piece.is_empty() {
                out.push((piece, ChainFeature::Segment(k)));
            }
        }

        // The wedge at vertex `k`: past the end of the previous segment and before the start of the
        // next. At a convex bend the two slabs already overlap and this comes back empty.
        cuts.clear();
        if k > 0 {
            cuts.push((chain[k], dir(k - 1)));
        }
        if k < last {
            cuts.push((chain[k], -dir(k)));
        }

        let piece = clip(&cuts);
        if !piece.is_empty() {
            out.push((piece, ChainFeature::Vertex(k)));
        }
    }

    out
}

/// Clip a convex polygon to the half plane of points where `(q - c) . n` is non-negative.
fn clip_half(poly: &[Point2], c: &Point2, n: &Vector2) -> Vec<Point2> {
    let mut out = Vec::with_capacity(poly.len() + 1);
    let side = |q: &Point2| (q - c).dot(n);

    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        let (fa, fb) = (side(&a), side(&b));

        if fa >= 0.0 {
            out.push(a);
        }
        if (fa > 0.0 && fb < 0.0) || (fa < 0.0 && fb > 0.0) {
            let t = fa / (fa - fb);
            out.push(Point2::from(a.coords * (1.0 - t) + b.coords * t));
        }
    }

    out
}

/// Where a point of a piece lands on the chain, and the radius blend waiting for it there.
///
/// `coeff` is the weight the landing puts on the merged vertex, whose radius is the unknown, and
/// `known` is the rest. `radii` carries `None` at the merged vertex and the fixed radius elsewhere.
pub(super) struct Landing {
    pub(super) at: Point3,
    pub(super) coeff: f64,
    pub(super) known: f64,
}

pub(super) fn chain_landing(
    q: &Point2,
    feature: ChainFeature,
    chain2: &[Point2],
    chain3: &[Point3],
    radii: &[Option<f64>],
) -> Landing {
    let blend = |k: usize, w: f64| -> (f64, f64) {
        match radii[k] {
            None => (w, 0.0),
            Some(r) => (0.0, w * r),
        }
    };

    match feature {
        ChainFeature::Vertex(k) => {
            let (coeff, known) = blend(k, 1.0);
            Landing {
                at: chain3[k],
                coeff,
                known,
            }
        }
        ChainFeature::Segment(k) => {
            // The parameter is taken in the projection, where the piece was cut, so it is inside the
            // segment by construction. Lifting it back onto the chain in space is affine along that
            // segment, which is what keeps the whole map affine.
            let t = Segment2::new_unchecked(chain2[k], chain2[k + 1])
                .scalar_projection(q)
                .clamp(0.0, 1.0);

            let (c0, c1) = blend(k, 1.0 - t);
            let (d0, d1) = blend(k + 1, t);

            Landing {
                at: Point3::from(chain3[k].coords * (1.0 - t) + chain3[k + 1].coords * t),
                coeff: c0 + d0,
                known: c1 + d1,
            }
        }
    }
}
