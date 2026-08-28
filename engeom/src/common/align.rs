pub mod information;
pub mod multi_params;

use crate::common::vec_f64::mean_and_stdev;
use parry3d_f64::na::Isometry;

pub use levenberg_marquardt::TerminationReason;

/// A classification of how a Levenberg-Marquardt solve ended, according to whether the parameters
/// it left behind can be used.
///
/// The `levenberg_marquardt` crate's own [`TerminationReason::was_successful`] draws a single line
/// between "converged" and "everything else", which is too coarse for an alignment.
/// [`TerminationReason::LostPatience`] means the solver hit its evaluation budget, but the
/// parameters it leaves behind are the best ones it found and are perfectly valid geometry. That
/// is a different situation from a solve which encountered a `NaN` and whose parameters are
/// meaningless.
///
/// The distinction matters most for ICP-style alignments, where the correspondences are
/// re-established every time the parameters change. The objective is then only piecewise smooth,
/// with a discontinuity wherever a point's closest match hops from one element of the target to
/// another. Near the solution, a point that sits close to a corner can flip between two matches
/// from one step to the next, and the convergence criteria (which assume a fixed, smooth
/// objective) may never be satisfied no matter how long the solver runs. Exhausting the budget is
/// a routine outcome in that regime rather than a failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveQuality {
    /// The solver met one of its convergence criteria, so the result is trustworthy.
    Converged,

    /// The solver ran out of its evaluation budget without meeting a convergence criterion. The
    /// parameters are the best ones the solver found and are geometrically valid, but convergence
    /// was never demonstrated. Raising the patience may help, though see the note above about why
    /// it often will not.
    Unconverged,

    /// The solve broke down and its parameters cannot be trusted.
    Failed,
}

impl SolveQuality {
    /// Classifies a Levenberg-Marquardt termination reason.
    pub fn from_termination(termination: &TerminationReason) -> Self {
        match termination {
            TerminationReason::Converged { .. }
            | TerminationReason::ResidualsZero
            | TerminationReason::Orthogonal => SolveQuality::Converged,

            TerminationReason::LostPatience => SolveQuality::Unconverged,

            TerminationReason::User(_)
            | TerminationReason::Numerical(_)
            | TerminationReason::NoImprovementPossible(_)
            | TerminationReason::NoParameters
            | TerminationReason::NoResiduals
            | TerminationReason::WrongDimensions(_) => SolveQuality::Failed,
        }
    }

    /// Returns true when the solve left behind parameters that can be used, whether or not
    /// convergence was actually demonstrated.
    pub fn is_usable(&self) -> bool {
        !matches!(self, SolveQuality::Failed)
    }

    /// Returns the worse of two qualities, ordered `Converged` then `Unconverged` then `Failed`.
    pub(crate) fn worse_of(self, other: Self) -> Self {
        match (self, other) {
            (SolveQuality::Failed, _) | (_, SolveQuality::Failed) => SolveQuality::Failed,
            (SolveQuality::Unconverged, _) | (_, SolveQuality::Unconverged) => {
                SolveQuality::Unconverged
            }
            _ => SolveQuality::Converged,
        }
    }
}

/// The reason a robust refinement loop stopped before completing all of the rounds that were
/// requested of it.
///
/// Refinement is an improvement on an alignment which is already usable, so running out of rounds
/// early is not an error. The alignment reported alongside this is the one produced by the last
/// round which did succeed.
#[derive(Debug, PartialEq, Eq)]
pub enum RefinementHalt {
    /// No usable noise bound could be estimated from the residuals, which happens when their
    /// spread has collapsed to zero. There is nothing left to reweight, so the fit is as good as
    /// the data allows.
    NoNoiseEstimate,

    /// Reweighting would have left fewer points carrying nonzero weight than there are free
    /// parameters, making the next solve rank-deficient.
    Underdetermined {
        /// How many points would still have carried nonzero weight.
        weighted: usize,

        /// How many free parameters the problem has.
        params: usize,
    },

    /// A refinement solve broke down. The parameters it produced were discarded and the alignment
    /// rolled back to the result of the previous round.
    SolveFailed(TerminationReason),
}

/// The full outcome of an alignment: the [`Alignment`] itself, plus a record of how each solve
/// that contributed to it terminated.
///
/// An alignment is only reported at all when it is usable, so the presence of this struct already
/// means there is a real answer to work with. What it adds is the ability to tell a proven
/// convergence from a merely plausible one, and to see whether robust refinement ran to
/// completion. Callers who do not care can go straight to [`AlignOutcome::into_alignment`].
#[derive(Debug)]
pub struct AlignOutcome<R, const D: usize> {
    alignment: Alignment<R, D>,
    solves: Vec<TerminationReason>,
    halt: Option<RefinementHalt>,
}

impl<R, const D: usize> AlignOutcome<R, D> {
    pub(crate) fn new(
        alignment: Alignment<R, D>,
        solves: Vec<TerminationReason>,
        halt: Option<RefinementHalt>,
    ) -> Self {
        Self {
            alignment,
            solves,
            halt,
        }
    }

    /// The alignment which was produced.
    pub fn alignment(&self) -> &Alignment<R, D> {
        &self.alignment
    }

    /// Consumes the outcome and returns just the alignment, discarding the diagnostics.
    pub fn into_alignment(self) -> Alignment<R, D> {
        self.alignment
    }

    /// How each solve whose result was kept terminated, beginning with the initial solve and
    /// followed by one entry per completed refinement round.
    ///
    /// A solve which broke down and was rolled back does not appear here; its termination reason
    /// is in [`AlignOutcome::halt`] instead.
    pub fn solves(&self) -> &[TerminationReason] {
        &self.solves
    }

    /// The number of robust refinement rounds which completed and contributed to the result.
    pub fn refinement_rounds(&self) -> usize {
        self.solves.len().saturating_sub(1)
    }

    /// The quality of the weakest solve that contributed to the result.
    ///
    /// This is never [`SolveQuality::Failed`], because a failed solve is never allowed to
    /// contribute. It is [`SolveQuality::Unconverged`] if any contributing solve exhausted its
    /// evaluation budget, which means the result is usable but its convergence was not proven.
    pub fn quality(&self) -> SolveQuality {
        self.solves
            .iter()
            .map(SolveQuality::from_termination)
            .fold(SolveQuality::Converged, SolveQuality::worse_of)
    }

    /// Whether every solve that contributed to the result met a convergence criterion.
    pub fn converged(&self) -> bool {
        self.quality() == SolveQuality::Converged
    }

    /// Why robust refinement stopped early, or `None` if it ran every round it was asked to (or
    /// was not requested at all).
    pub fn halt(&self) -> Option<&RefinementHalt> {
        self.halt.as_ref()
    }
}

/// The full outcome of a simultaneous alignment of several bodies: one [`Alignment`] per body,
/// plus a record of how the solves that produced them terminated.
///
/// Note that the solve record is **shared** rather than per body. A multi-body adjustment is a
/// single least-squares problem over all of the bodies at once, so there is one sequence of solves
/// and one halt reason for the whole thing, not one per body. Only the alignments and their
/// residuals are per body.
#[derive(Debug)]
pub struct MultiAlignOutcome<R, const D: usize> {
    alignments: Vec<Alignment<R, D>>,
    solves: Vec<TerminationReason>,
    halt: Option<RefinementHalt>,
}

impl<R, const D: usize> MultiAlignOutcome<R, D> {
    pub(crate) fn new(
        alignments: Vec<Alignment<R, D>>,
        solves: Vec<TerminationReason>,
        halt: Option<RefinementHalt>,
    ) -> Self {
        Self {
            alignments,
            solves,
            halt,
        }
    }

    /// The alignment produced for each body, in body order.
    pub fn alignments(&self) -> &[Alignment<R, D>] {
        &self.alignments
    }

    /// The alignment produced for one body.
    pub fn alignment(&self, body: usize) -> &Alignment<R, D> {
        &self.alignments[body]
    }

    /// The number of bodies.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    /// Whether there are no bodies at all.
    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }

    /// How each solve whose result was kept terminated, beginning with the initial solve and
    /// followed by one entry per completed refinement round.
    pub fn solves(&self) -> &[TerminationReason] {
        &self.solves
    }

    /// The number of robust refinement rounds which completed and contributed to the result.
    pub fn refinement_rounds(&self) -> usize {
        self.solves.len().saturating_sub(1)
    }

    /// The quality of the weakest solve that contributed to the result. See
    /// [`AlignOutcome::quality`], which this mirrors.
    pub fn quality(&self) -> SolveQuality {
        self.solves
            .iter()
            .map(SolveQuality::from_termination)
            .fold(SolveQuality::Converged, SolveQuality::worse_of)
    }

    /// Whether every solve that contributed to the result met a convergence criterion.
    pub fn converged(&self) -> bool {
        self.quality() == SolveQuality::Converged
    }

    /// Why robust refinement stopped early, or `None` if it ran every round it was asked to.
    pub fn halt(&self) -> Option<&RefinementHalt> {
        self.halt.as_ref()
    }
}

/// A container for the results of an alignment operation, including the full transformation, the
/// various component transformations, and the residuals of the alignment.
#[derive(Debug, Clone)]
pub struct Alignment<R, const D: usize> {
    full: Isometry<f64, R, D>,
    alignment: Isometry<f64, R, D>,
    local_origin: Isometry<f64, R, D>,
    offset: Isometry<f64, R, D>,
    residuals: Vec<f64>,
}

impl<R, const D: usize> Alignment<R, D> {
    pub(crate) fn new(
        full: Isometry<f64, R, D>,
        alignment: Isometry<f64, R, D>,
        local_origin: Isometry<f64, R, D>,
        offset: Isometry<f64, R, D>,
        residuals: Vec<f64>,
    ) -> Self {
        Self {
            full,
            alignment,
            local_origin,
            offset,
            residuals,
        }
    }

    /// Gets the full transformation, which is a single transformation that brings the test entity
    /// geometry directly to the target geometry. This is a composite of the local origin's
    /// inverse, the alignment, and the work offset ($O * A * L^{-1}$).
    ///
    /// This is the one to apply to the test geometry once the alignment is finished. Contrast it
    /// with [`Alignment::local_transform`], which is the same motion stripped of its framing.
    pub fn full_transform(&self) -> &Isometry<f64, R, D> {
        &self.full
    }

    /// Gets the alignment transformation $A$, which is the motion produced by the optimized
    /// parameters (tx, ty, tz, rx, ry, rz) expressed in the frame of the local origin.
    ///
    /// This is *not* the transformation to apply to the test geometry. Reading
    /// $O * A * L^{-1}$ right to left, $L^{-1}$ puts a point into the local origin's frame, $A$
    /// moves it while it is there, and $O$ maps it back out, so $A$ is only meaningful applied to
    /// local-frame coordinates. Use [`Alignment::full_transform`] to move geometry, and this to
    /// inspect what the parameters actually did about the origin you chose.
    pub fn local_transform(&self) -> &Isometry<f64, R, D> {
        &self.alignment
    }

    /// Gets the local origin, which is a transformation from world origin to the local origin
    /// used by the alignment algorithm.
    pub fn local_origin(&self) -> &Isometry<f64, R, D> {
        &self.local_origin
    }

    /// Gets the offset transform, which is a transformation applied to the test entity geometry
    /// after the alignment transformation. It is typically used to counteract the local origin, or
    /// to provide an initial guess.
    pub fn offset(&self) -> &Isometry<f64, R, D> {
        &self.offset
    }

    /// Gets the residuals of the alignment, which are the difference between the target geometry
    /// and the test entity geometry after the alignment transformation. The order is the same as
    /// the order of the target geometry.
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }

    /// Calculates the average residual of the alignment.
    pub fn residual_mean(&self) -> f64 {
        self.residuals.iter().sum::<f64>() / self.residuals.len() as f64
    }

    /// Calculates the mean and standard deviation of the residuals of the alignment.
    pub fn residual_mean_std_dev(&self) -> (f64, f64) {
        mean_and_stdev(&self.residuals).unwrap()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistMode {
    ToPoint,
    ToPlane,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converged_terminations_are_classified() {
        for t in [
            TerminationReason::Converged {
                ftol: true,
                xtol: false,
            },
            TerminationReason::Converged {
                ftol: false,
                xtol: true,
            },
            TerminationReason::ResidualsZero,
            TerminationReason::Orthogonal,
        ] {
            assert_eq!(SolveQuality::from_termination(&t), SolveQuality::Converged);
            assert!(SolveQuality::from_termination(&t).is_usable());
        }
    }

    #[test]
    fn exhausted_budget_is_usable_but_unconverged() {
        // The distinction this whole classification exists for: the solver gave up on proving
        // convergence, but the parameters it left behind are still the best it found.
        let q = SolveQuality::from_termination(&TerminationReason::LostPatience);
        assert_eq!(q, SolveQuality::Unconverged);
        assert!(q.is_usable());
    }

    #[test]
    fn broken_terminations_are_not_usable() {
        for t in [
            TerminationReason::User("residuals"),
            TerminationReason::Numerical("jacobian"),
            TerminationReason::NoImprovementPossible("ftol"),
            TerminationReason::NoParameters,
            TerminationReason::NoResiduals,
            TerminationReason::WrongDimensions("jacobian"),
        ] {
            assert_eq!(SolveQuality::from_termination(&t), SolveQuality::Failed);
            assert!(!SolveQuality::from_termination(&t).is_usable());
        }
    }

    #[test]
    fn worse_of_orders_the_qualities() {
        use SolveQuality::{Converged, Failed, Unconverged};

        assert_eq!(Converged.worse_of(Converged), Converged);
        assert_eq!(Converged.worse_of(Unconverged), Unconverged);
        assert_eq!(Unconverged.worse_of(Converged), Unconverged);
        assert_eq!(Unconverged.worse_of(Failed), Failed);
        assert_eq!(Failed.worse_of(Converged), Failed);
    }
}
