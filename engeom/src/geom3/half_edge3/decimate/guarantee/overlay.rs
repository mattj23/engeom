//! The projected overlay method
//!
//! Both stars are disks with the same link polygon, so projecting both along the old star's average
//! normal gives a bijection between them for "free", with pieces given by the overlay of the two
//! projected triangulations. See the parent module's documentation for why the corners of those
//! pieces are enough.

pub(super) mod sliver;

use super::StarFace;
use super::boundary::Outline;
use super::collapse::{Collapse, ErrorRule};
use super::constraint::{Bound, Constraint, ErrorBound};
use crate::common::barycentric::{barycentric_point, barycentric_within2, signed_area2};
use crate::geom2::Segment2;
use crate::{Point2, Point3, Vector3};
use alum::{Handle, VH};
use sliver::{chain_landing, chain_pieces, sliver_corners};

/// Where two planar segments cross, as a parameter along each, or `None` if the crossing is not
/// strictly between their endpoints.
///
/// [`Segment2::intersection_param`] does the arithmetic and counts a crossing at an endpoint; this
/// adds the strictness. Endpoint crossings are skipped because those points are already vertices of
/// the overlay, and are constrained where the vertices are handled. Admitting them again here would
/// emit a duplicate constraint from a degenerate node.
pub(super) fn segment_cross(a: &Point2, b: &Point2, c: &Point2, d: &Point2) -> Option<(f64, f64)> {
    let first = Segment2::new_unchecked(*a, *b);
    let second = Segment2::new_unchecked(*c, *d);
    let (t, u) = first.intersection_param(&second)?;

    let eps = 1.0e-9;
    let inside = |x: f64| x > eps && x < 1.0 - eps;
    (inside(t) && inside(u)).then_some((t, u))
}

/// Walk the interior edges of a star, visiting each distinct edge exactly once.
///
/// An edge is interior exactly when it touches one of the two endpoints of the collapsed edge; the
/// rest form the link polygon, which both stars share and which therefore cannot be crossed. Two
/// faces carry each interior edge, so the walk dedupes by vertex pair.
///
/// `face_ok` filters at the face level rather than inside `visit`, which matters: skipping a face
/// from within the visitor would still have marked its edges seen, and the other face carrying the
/// same edge would then be dropped as a duplicate.
fn for_each_interior_edge<P, F>(star: &[StarFace], v1: VH, v2: VH, face_ok: P, mut visit: F)
where
    P: Fn(&StarFace) -> bool,
    F: FnMut(usize, usize, usize),
{
    let mut seen: Vec<(u32, u32)> = Vec::with_capacity(24);

    for (i, f) in star.iter().enumerate() {
        if !face_ok(f) {
            continue;
        }

        for k in 0..3 {
            let j = (k + 1) % 3;
            let (hk, hj) = (f.v[k], f.v[j]);
            if !(hk == v1 || hk == v2 || hj == v1 || hj == v2) {
                continue;
            }

            let key = if hk.index() < hj.index() {
                (hk.index(), hj.index())
            } else {
                (hj.index(), hk.index())
            };
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);

            visit(i, k, j);
        }
    }
}

/// The interior edges of the old star, as endpoint positions and radii.
pub(super) fn interior_edges(
    star: &[StarFace],
    v1: VH,
    v2: VH,
    e: &[f64],
) -> Vec<(Point3, f64, Point3, f64)> {
    let mut out = Vec::with_capacity(24);

    for_each_interior_edge(
        star,
        v1,
        v2,
        |_| true,
        |i, k, j| {
            let f = &star[i];
            out.push((
                f.p[k],
                e[f.v[k].index() as usize],
                f.p[j],
                e[f.v[j].index() as usize],
            ));
        },
    );

    out
}

/// The interior edges of the new star, flagging which endpoint is the merged vertex so its unknown
/// radius can be carried symbolically into the constraints.
pub(super) fn interior_edges_new(
    star: &[StarFace],
    new: &[[Point3; 3]],
    v1: VH,
    v2: VH,
    e: &[f64],
) -> Vec<(Point3, bool, f64, Point3, bool, f64)> {
    let mut out = Vec::with_capacity(24);

    for_each_interior_edge(
        star,
        v1,
        v2,
        |f| !f.vanishing,
        |i, k, j| {
            let f = &star[i];
            let (hk, hj) = (f.v[k], f.v[j]);
            let (mk, mj) = (hk == v1 || hk == v2, hj == v1 || hj == v2);

            // Both ends merged means the edge is gone rather than interior.
            if mk && mj {
                return;
            }

            let rk = if mk { 0.0 } else { e[hk.index() as usize] };
            let rj = if mj { 0.0 } else { e[hj.index() as usize] };
            out.push((new[i][k], mk, rk, new[i][j], mj, rj));
        },
    );

    out
}

/// The error radius the merged vertex would take, by overlaying the two stars in a shared
/// projection.
///
/// # The idea
///
/// Both stars are disks with the *same* boundary, the link polygon. Project both along the old
/// star's average normal and, if neither projection folds, each is a triangulation of the same
/// planar polygon. That gives a bijection between them for free: a point of one star
/// corresponds to the point of the other with the same projection. It is piecewise affine, with
/// pieces given by the overlay of the two projected triangulations.
///
/// A bijection is worth more than two separate maps, because a single constraint at each
/// corresponding pair serves both invariants at once. There is no longer a "reach" direction
/// and a "cover" direction to satisfy separately.
///
/// # The premise, and that one weird trick which breaks it
///
/// All of the above rests on the two stars covering the _same_ region. On a closed mesh they
/// always do: the region is bounded by the link polygon, whose vertices the collapse neither
/// moves nor deletes.
///
/// Where the star reaches the edge of the mesh it is bounded partly by mesh boundary edges
/// instead, and those are incident to `v1` or `v2`, so a collapse can move them. When it does,
/// the two regions differ and the correspondence is no longer a bijection. What goes wrong then
/// is silent rather than loud: the part of the new star with no old counterpart is covered by no
/// constraint, so it is never charged for, and the radius comes out too small. Measured on an
/// unlocked bunny outline that was 1.9x over the tolerance, while the face-to-face map, which
/// never assumed a shared region, stayed under it.
///
/// So `outline_moves` is a precondition rather than something detected here, and it is the
/// caller's because only the caller can see which vertices are on the mesh boundary. A locked
/// boundary never sets it: `v1` cannot be a boundary vertex, since the legality check refuses a
/// collapse whose tail is locked, and where `v2` is one the merged position is pinned to it.
/// That is Gueziec's Type II boundary edge, and it is why the overlay handles an open mesh's
/// interior collapses without ever reaching the fallback.
///
/// See [`ErrorBound`] for the three ways this can end.
pub(super) struct ProjectedOverlay;

impl ErrorRule for ProjectedOverlay {
    fn bound(&self, cx: &Collapse<'_>) -> ErrorBound {
        let Collapse {
            star,
            new,
            v1,
            v2,
            p0,
            e,
            outline,
            ..
        } = *cx;

        // Whether this placement actually parts the two regions, and along which run if so.
        let moving = match outline {
            Outline::Unsupported => return ErrorBound::NotApplicable,
            Outline::Open(chain) if chain.moves_under(p0) => Some(chain),
            _ => None,
        };

        let merged = cx.merged();

        // Area weighted normal of the old star, which is the direction least likely to fold either
        // projection.
        let mut normal = Vector3::zeros();
        for f in star.iter() {
            normal += (f.p[1] - f.p[0]).cross(&(f.p[2] - f.p[0]));
        }
        let Some(normal) = normal.try_normalize(1.0e-14) else {
            return ErrorBound::NotApplicable;
        };

        let seed = if normal.x.abs() < 0.9 {
            Vector3::x()
        } else {
            Vector3::y()
        };
        let Some(ex) = normal.cross(&seed).try_normalize(1.0e-14) else {
            return ErrorBound::NotApplicable;
        };
        let ey = normal.cross(&ex);
        let flat = |q: &Point3| Point2::new(q.coords.dot(&ex), q.coords.dot(&ey));

        // Both projections have to be injective for the correspondence to exist at all, which for a
        // triangulated disk means every face keeping the same orientation and none going flat. The
        // projected triangles are kept, since locating points and finding slivers both want them.
        let mut sign = 0.0f64;
        let mut old2: Vec<[Point2; 3]> = Vec::with_capacity(star.len());
        let mut new2: Vec<Option<[Point2; 3]>> = Vec::with_capacity(star.len());

        for (i, f) in star.iter().enumerate() {
            old2.push([flat(&f.p[0]), flat(&f.p[1]), flat(&f.p[2])]);
            new2.push(
                (!f.vanishing).then(|| [flat(&new[i][0]), flat(&new[i][1]), flat(&new[i][2])]),
            );

            for t in [Some(old2[i]), new2[i]].into_iter().flatten() {
                let a = signed_area2(&t[0], &t[1], &t[2]);
                if a.abs() < 1.0e-12 {
                    return ErrorBound::NotApplicable;
                }
                if sign == 0.0 {
                    sign = a.signum();
                } else if a.signum() != sign {
                    return ErrorBound::NotApplicable;
                }
            }
        }

        let new2_present: Vec<[Point2; 3]> = new2.iter().flatten().copied().collect();
        let mut bound = Bound::new();

        // The merged vertex against the old star: same projection, lifted onto whichever old face
        // covers it. The landing is entirely on the merged vertex, so its radius is the only
        // unknown and the constraint is direct.
        let q = flat(&merged);
        let mut located = false;
        for (i, f) in star.iter().enumerate() {
            let t = &old2[i];
            let Some(bary) = barycentric_within2(&t[0], &t[1], &t[2], &q) else {
                continue;
            };
            let x = barycentric_point(&f.p[0], &f.p[1], &f.p[2], bary);
            let mut need = (x - merged).norm();
            for k in 0..3 {
                need += bary[k] * e[f.v[k].index() as usize];
            }
            if !bound.add(Constraint::on_merged(need)) {
                return ErrorBound::Unsatisfiable;
            }
            located = true;
            break;
        }
        if !located && moving.is_none() {
            // With the regions equal every point of one star lies on the other, so a miss here is
            // numerical rather than geometric and there is nothing sound to say. When they differ it
            // means only that the merged vertex is off the old star, and the sliver pass below is
            // what covers it.
            return ErrorBound::NotApplicable;
        }

        // Each endpoint against the new star, the same way round.
        for u in [v1, v2] {
            let Some(idx) = star
                .iter()
                .find_map(|f| f.v.iter().position(|h| *h == u).map(|k| (f, k)))
            else {
                continue;
            };
            let x = idx.0.p[idx.1];
            let q = flat(&x);

            let mut located = false;
            for (i, f) in star.iter().enumerate() {
                let Some(t) = &new2[i] else {
                    continue;
                };
                let Some(bary) = barycentric_within2(&t[0], &t[1], &t[2], &q) else {
                    continue;
                };
                let y = barycentric_point(&new[i][0], &new[i][1], &new[i][2], bary);
                let need = (y - x).norm() + e[u.index() as usize];

                let parts = cx.landing_parts(&f.v, &bary);

                if !bound.add(Constraint::from_landing(need, parts)) {
                    return ErrorBound::Unsatisfiable;
                }
                located = true;
                break;
            }
            if !located && moving.is_none() {
                return ErrorBound::NotApplicable;
            }
        }

        // The crossings. Only the interior edges of each star can cross, since the two share their
        // boundary exactly, and an interior edge cannot leave the polygon.
        let old_edges = interior_edges(star, v1, v2, e);
        let new_edges = interior_edges_new(star, new, v1, v2, e);

        for (a, ea, b, eb) in old_edges.iter() {
            let (a2, b2) = (flat(a), flat(b));
            for (c, cm, ec, d, dm, ed) in new_edges.iter() {
                let (c2, d2) = (flat(c), flat(d));
                let Some((s, t)) = segment_cross(&a2, &b2, &c2, &d2) else {
                    continue;
                };

                let x = Point3::from(a.coords * (1.0 - s) + b.coords * s);
                let e_here = ea * (1.0 - s) + eb * s;
                let y = Point3::from(c.coords * (1.0 - t) + d.coords * t);
                let need = (y - x).norm() + e_here;

                let parts = [(*cm, 1.0 - t, *ec), (*dm, t, *ed)];
                if !bound.add(Constraint::from_landing(need, parts)) {
                    return ErrorBound::Unsatisfiable;
                }
            }
        }

        // Everything above covers the region the two stars share, and on a closed outline that is
        // the whole of both. What is left when the outline moves is a sliver of each star with no
        // counterpart on the other, and the invariants are pointwise existentials, so those slivers
        // only need a cover of their own rather than a place in one global bijection.
        let Some(chain) = moving else {
            return ErrorBound::Bound(bound.finish());
        };

        let last = chain.len - 1;
        let new_chain = [chain.p[0], merged, chain.p[last]];
        let new_radii = [
            Some(e[chain.v[0].index() as usize]),
            None,
            Some(e[chain.v[last].index() as usize]),
        ];
        let new_chain2: Vec<Point2> = new_chain.iter().map(&flat).collect();

        let old_chain: Vec<Point3> = chain.iter().map(|(_, q)| q).collect();
        let old_radii: Vec<Option<f64>> = chain
            .iter()
            .map(|(v, _)| Some(e[v.index() as usize]))
            .collect();
        let old_chain2: Vec<Point2> = old_chain.iter().map(&flat).collect();

        // Invariant (𝔸): the part of each old face the new star does not reach, carried onto the
        // chain the collapse leaves behind.
        for (i, f) in star.iter().enumerate() {
            let corners = sliver_corners(&old2[i], &new2_present, &new_chain2);
            if corners.is_empty() {
                continue;
            }

            for (piece, feature) in chain_pieces(&corners, &new_chain2) {
                for q in piece.iter() {
                    let t = &old2[i];
                    let Some(w) = barycentric_within2(&t[0], &t[1], &t[2], q) else {
                        continue;
                    };

                    let x = barycentric_point(&f.p[0], &f.p[1], &f.p[2], w);
                    let e_x = (0..3)
                        .map(|k| w[k] * e[f.v[k].index() as usize])
                        .sum::<f64>();
                    let site = chain_landing(q, feature, &new_chain2, &new_chain, &new_radii);

                    if !bound.add(Constraint {
                        coeff: site.coeff,
                        known: site.known,
                        need: (x - site.at).norm() + e_x,
                    }) {
                        return ErrorBound::Unsatisfiable;
                    }
                }
            }
        }

        // Invariant (𝔹): the part of each new face the old star does not reach, carried back onto
        // the chain the collapse started from. Here it is the sliver corner which carries the
        // unknown radius and the landing whose radius is already fixed, the mirror of above.
        for (i, f) in star.iter().enumerate() {
            let Some(t) = &new2[i] else {
                continue;
            };

            let corners = sliver_corners(t, &old2, &old_chain2);
            if corners.is_empty() {
                continue;
            }

            for (piece, feature) in chain_pieces(&corners, &old_chain2) {
                for q in piece.iter() {
                    let Some(w) = barycentric_within2(&t[0], &t[1], &t[2], q) else {
                        continue;
                    };

                    let y = barycentric_point(&new[i][0], &new[i][1], &new[i][2], w);
                    let site = chain_landing(q, feature, &old_chain2, &old_chain, &old_radii);
                    let need = (y - site.at).norm() + site.known;

                    let parts = cx.landing_parts(&f.v, &w);

                    if !bound.add(Constraint::from_landing(need, parts)) {
                        return ErrorBound::Unsatisfiable;
                    }
                }
            }
        }

        ErrorBound::Bound(bound.finish())
    }
}
