//! What a single constraint on the merged vertex's radius looks like, how a set of them folds into
//! a radius, and the three ways a method can end.
//!
//! The methods differ in which constraints their correspondence can find. They do not differ in
//! how a set of constraints becomes a radius, and keeping that difference visible is why this is
//! its own module rather than three private helpers.

/// The three ways any of the error methods can end.
///
/// The distinction between the last two is the whole reason this is an enum rather than an
/// `Option`: one says the collapse is impossible, the other says only that *this method* cannot
/// judge it and something else should.
pub(super) enum ErrorBound {
    /// The smallest radius at the merged vertex which satisfies every constraint.
    Bound(f64),

    /// A constraint no radius can satisfy. The collapse is refused, and no fallback will help.
    Unsatisfiable,

    /// This method has nothing to say about this collapse, so a different one has to judge it.
    ///
    /// The projected overlay reaches this when the two stars have no shared projection, because one
    /// of them folds under it or a point lands outside the link polygon; it falls back to the
    /// face-to-face map, which needs no projection at all.
    NotApplicable,
}

/// The weight below which a landing is treated as carrying none of the merged vertex.
///
/// A constraint with no weight on the unknown is not a constraint on it, so it is either already
/// satisfied by the fixed radii or satisfiable by nothing.
const MERGED_WEIGHT_EPS: f64 = 1.0e-9;

/// One linear constraint on the radius the merged vertex is about to take.
///
/// Under Objective B every radius in the star except the merged one is already fixed, so every
/// constraint any of the methods can produce reduces to this single shape:
///
/// ```text
/// coeff * e_merged + known >= need
/// ```
///
/// `need` is what the landing has to reach, a distance plus whatever radius is already carried on
/// the other side. `known` is the part of the barycentric blend at the landing which sits on link
/// vertices, whose radii this collapse does not change. `coeff` is the weight that same blend puts
/// on the merged vertex, which is the only thing left to solve for.
///
/// Having one shape for all of them is not tidiness. It is what makes two error methods comparable:
/// they differ in which constraints they can find, and not in how a set of constraints is turned
/// into a radius.
#[derive(Debug, Clone, Copy)]
pub(super) struct Constraint {
    pub(super) coeff: f64,
    pub(super) known: f64,
    pub(super) need: f64,
}

impl Constraint {
    /// The constraint `e_merged >= need`, which is what a landing entirely on the merged vertex
    /// leaves.
    pub(super) fn on_merged(need: f64) -> Self {
        Self {
            coeff: 1.0,
            known: 0.0,
            need,
        }
    }

    /// Build from a barycentric landing, splitting the blend into the weight on the merged vertex
    /// and the radius already fixed at the link vertices.
    ///
    /// Parts arrive as `(is_merged, weight, radius)`. A merged part contributes its weight to the
    /// unknown and its radius is ignored, since what that radius becomes is the whole question.
    /// Both endpoints of the collapsed edge count as merged, so their weights combine.
    pub(super) fn from_landing(
        need: f64,
        parts: impl IntoIterator<Item = (bool, f64, f64)>,
    ) -> Self {
        let mut coeff = 0.0;
        let mut known = 0.0;
        for (is_merged, weight, radius) in parts {
            if is_merged {
                coeff += weight;
            } else {
                known += weight * radius;
            }
        }
        Self { coeff, known, need }
    }
}

/// Folds constraints into the smallest merged radius that satisfies all of them, which is Objective
/// B's one-dimensional solve.
///
/// Constraints are folded as they are found rather than collected into a list first. The acceptance
/// test runs several times per collapse on the hot path, and an allocation there is not worth the
/// tidier signature.
#[derive(Debug, Clone, Copy)]
pub(super) struct Bound {
    lower: f64,
}

impl Bound {
    pub(super) fn new() -> Self {
        Self { lower: 0.0 }
    }

    /// Fold one constraint in, returning `false` if no radius can satisfy it.
    #[must_use]
    pub(super) fn add(&mut self, c: Constraint) -> bool {
        if c.coeff > MERGED_WEIGHT_EPS {
            self.lower = self.lower.max((c.need - c.known) / c.coeff);
            true
        } else {
            // Nothing here moves with the merged radius, so the fixed part either already covers
            // what the landing needs or nothing will.
            c.known >= c.need
        }
    }

    /// The smallest radius satisfying everything folded in so far. Radii are never negative.
    pub(super) fn finish(self) -> f64 {
        self.lower.max(0.0)
    }
}
