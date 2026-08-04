//! Gueziec's Contraction Method: widen the ball at the merged vertex until it swallows both.

use super::collapse::{Collapse, ErrorRule};
use super::constraint::ErrorBound;
use alum::Handle;

/// The error radius the merged vertex would take, by Gueziec's Contraction Method.
///
/// The requirement is that the new error volume contain the old one, and that every ball of the
/// new volume contain a ball of the old. Satisfying both at a single vertex is one line: make
/// the ball at the merged position large enough to swallow the balls at both endpoints, which
/// is `max` over the two of the distance traveled plus the radius already there.
///
/// Sound, and by some distance the loosest of the paper's three schemes: the radius grows by at
/// least the distance the vertex moved on every collapse, whether or not the surface actually
/// went anywhere. Because that distance is at least half the collapsed edge, no edge shorter
/// than about twice the tolerance can ever be collapsed. Kept as an option because it is the
/// one scheme whose validity needs no argument beyond the sentence above, which makes it the
/// reference the others are checked against.
pub(super) struct Contraction;

impl ErrorRule for Contraction {
    fn bound(&self, c: &Collapse<'_>) -> ErrorBound {
        let mut worst = 0.0f64;
        for v in [c.v1, c.v2] {
            let d = (c.p0 - c.points[v.index() as usize]).norm() + c.radius(v);
            worst = worst.max(d);
        }
        ErrorBound::Bound(worst)
    }
}
