//! The shape an error method is handed, and the shape it has to implement.
//!
//! Every method is an *argument* that the radius it returns is large enough: build a
//! correspondence between the two stars, collect constraints at the corners of its pieces, and
//! appeal to concavity for everything in between. They differ only in which constraints their
//! correspondence can find. They do not differ in what they are told about the collapse, and they
//! do not differ in how a set of constraints becomes a radius, which is [`super::constraint`].
//!
//! Keeping those three things apart is what makes a new method cheap: it is a file with one
//! [`ErrorRule`] impl in it, and the enum arm in [`super::ToleranceVolumeDecimator::accept`] that
//! selects it.

use super::StarFace;
use super::boundary::Outline;
use super::constraint::ErrorBound;
use crate::{Point3, Vector3};
use alum::{Handle, VH};

/// Everything an error method is given about one candidate collapse.
///
/// Every field is either borrowed or trivially copied, so this costs nothing to build. It is
/// assembled once per candidate position in
/// [`ToleranceVolumeDecimator::accept`](super::ToleranceVolumeDecimator::accept), on the hot path.
#[derive(Clone, Copy)]
pub(super) struct Collapse<'a> {
    /// Every face touching either endpoint, as it stands before the collapse.
    pub(super) star: &'a [StarFace],

    /// The triangle replacing each face of the star, in the same order. A face which vanishes into
    /// the collapsed edge still has an entry, which is why [`StarFace::vanishing`] has to be
    /// consulted rather than the lengths compared.
    pub(super) new: &'a [[Point3; 3]],

    /// Every vertex position in the mesh, indexed by vertex handle.
    pub(super) points: &'a [Vector3],

    /// The error volume, as one radius per vertex, indexed by vertex handle.
    pub(super) e: &'a [f64],

    /// The endpoints of the collapsed edge, Guéziec's v₁ and v₂.
    ///
    /// `v1` is the halfedge's tail and is deleted; `v2` is its head and survives. No error method
    /// should care which is which: what matters is that both are replaced by one vertex whose
    /// radius is the unknown being solved for, which is what [`Collapse::is_merged`] asks.
    pub(super) v1: VH,
    pub(super) v2: VH,

    /// Where the merged vertex would go, which is Guéziec's v₀.
    ///
    /// A position rather than a handle, because before the collapse is applied there is no handle
    /// for it: a half-edge collapse deletes the tail and moves the head, so v₀ ends up living on
    /// `v2`. This is a *candidate* position, and one collapse is judged at several of them.
    pub(super) p0: Vector3,

    /// What this collapse can do to the star's outline.
    ///
    /// A fact about the collapse rather than about any method, which is why it is carried here.
    /// It cannot be worked out from the star alone by anything downstream, and every method which
    /// builds a correspondence between the two stars needs it for the same reason: the
    /// correspondence rests on the two covering the same region, and this is what says whether
    /// they do.
    pub(super) outline: Outline,
}

impl Collapse<'_> {
    /// Where the merged vertex would go, as a point.
    pub(super) fn merged(&self) -> Point3 {
        Point3::from(self.p0)
    }

    /// Whether a vertex is one the collapse merges away, and so carries the unknown radius.
    pub(super) fn is_merged(&self, v: VH) -> bool {
        v == self.v1 || v == self.v2
    }

    /// The error radius standing at a vertex before the collapse.
    pub(super) fn radius(&self, v: VH) -> f64 {
        self.e[v.index() as usize]
    }

    /// A barycentric landing on a face, split the way [`Constraint::from_landing`] wants it.
    ///
    /// Each corner arrives as `(is_merged, weight, radius)`. The merged parts contribute their
    /// weight to the unknown and their radii are ignored, since what those become is the question.
    ///
    /// [`Constraint::from_landing`]: super::constraint::Constraint::from_landing
    pub(super) fn landing_parts<'b>(
        &'b self,
        v: &'b [VH; 3],
        bary: &'b [f64; 3],
    ) -> impl Iterator<Item = (bool, f64, f64)> + 'b {
        v.iter()
            .zip(bary.iter())
            .map(|(vj, w)| (self.is_merged(*vj), *w, self.radius(*vj)))
    }
}

/// One way of working out the radius the merged vertex has to take.
///
/// Implementations are stateless, so they are unit structs rather than anything constructed. The
/// trait exists to declare the shape, not to be held behind a pointer.
pub(super) trait ErrorRule {
    /// The smallest radius at the merged vertex which this method's correspondence can justify, or
    /// why it has no answer. See [`ErrorBound`] for the three ways this can end.
    fn bound(&self, c: &Collapse<'_>) -> ErrorBound;
}
