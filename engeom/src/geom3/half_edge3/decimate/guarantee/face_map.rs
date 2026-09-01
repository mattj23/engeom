//! Face-to-face affine map: carry each face of the star onto the face which replaces it.
//!
//! This is a quick and easy implementation, easy to prove it's sound because its pieces are
//! one triangle onto one triangle. That's why it survives as the projected overlay's fallback even
//! though it decimates pretty badly on its own. See the parent module's documentation.

use super::StarFace;
use super::collapse::{Collapse, ErrorRule};
use super::constraint::{Bound, Constraint, ErrorBound};
use crate::Point3;
use crate::common::barycentric::{barycentric_point, closest_barycentric};
use alum::Handle;

/// The error radius the merged vertex would take, by mapping each face of the star affinely
/// onto the face which replaces it.
///
/// Vertex checks are enough here because of the convexity argument.
pub(super) struct AffineFaceMap;

impl ErrorRule for AffineFaceMap {
    fn bound(&self, cx: &Collapse<'_>) -> ErrorBound {
        match radius(cx) {
            Some(v) => ErrorBound::Bound(v),
            // A constraint no radius can satisfy, which for this map means an endpoint landing that
            // carries no weight on the merged vertex and fixed radii already too small for it.
            None => ErrorBound::Unsatisfiable,
        }
    }
}

fn radius(cx: &Collapse<'_>) -> Option<f64> {
    let Collapse {
        star, new, v1, v2, ..
    } = *cx;

    let mut lower = 0.0f64;

    for (i, face) in star.iter().enumerate() {
        if face.vanishing {
            // No face occupies these slots afterwards. Map onto whichever surviving face shares
            // the third vertex, so that vertex still maps to itself and the map stays affine.
            let third = face.v.iter().position(|v| *v != v1 && *v != v2);
            let Some(third) = third else {
                continue; // Both endpoints and nothing else; nothing to map onto.
            };
            let c = face.v[third];

            let mut best: Option<f64> = None;
            for (j, host) in star.iter().enumerate() {
                if host.vanishing || !host.v.contains(&c) {
                    continue;
                }
                let bound = face_bound(cx, face, &new[j], host)?;
                if best.is_none_or(|b| bound < b) {
                    best = Some(bound);
                }
            }

            lower = lower.max(best?);
        } else {
            lower = lower.max(face_bound(cx, face, &new[i], face)?);
            lower = lower.max(reach_bound(cx, face));
        }
    }

    Some(lower.max(0.0))
}

/// The lower bound on the merged radius imposed by carrying `face` onto the triangle `target`,
/// whose vertex handles are `host`'s with the endpoints merged.
///
/// One constraint per endpoint the face carries. A link vertex maps to itself over zero
/// distance and its radius is unchanged, so it constrains nothing.
///
/// `None` when a constraint cannot be met by any choice of the merged radius, which happens
/// when the landing carries no weight on the merged vertex and the fixed radii there are
/// already too small.
fn face_bound(
    cx: &Collapse<'_>,
    face: &StarFace,
    target: &[Point3; 3],
    host: &StarFace,
) -> Option<f64> {
    let Collapse { v1, v2, e, .. } = *cx;

    let mut bound = Bound::new();

    for k in 0..3 {
        let u = face.v[k];
        if u != v1 && u != v2 {
            continue;
        }

        let bary = closest_barycentric(&target[0], &target[1], &target[2], &face.p[k]);
        let at = barycentric_point(&target[0], &target[1], &target[2], bary);
        let need = (at - face.p[k]).norm() + e[u.index() as usize];

        let parts = cx.landing_parts(&host.v, &bary);

        if !bound.add(Constraint::from_landing(need, parts)) {
            return None;
        }
    }

    Some(bound.finish())
}

/// The lower bound imposed by the merged vertex having to keep reaching back to the old surface.
///
/// This is the other invariant: a ball of the new volume is only a bound on anything if it still
/// contains a ball of the old volume, and so still meets the original surface. Mapping the new
/// face back onto the old one in the same slots is again affine, so the corners suffice, and
/// only the merged vertex's corner says anything.
fn reach_bound(cx: &Collapse<'_>, face: &StarFace) -> f64 {
    let Collapse { p0, e, .. } = *cx;

    let bary = closest_barycentric(&face.p[0], &face.p[1], &face.p[2], &Point3::from(p0));
    let at = barycentric_point(&face.p[0], &face.p[1], &face.p[2], bary);

    let mut reach = (at - Point3::from(p0)).norm();
    for k in 0..3 {
        reach += bary[k] * e[face.v[k].index() as usize];
    }
    reach
}
