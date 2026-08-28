//! Shared formatting for the solver's own reporting vocabulary, so that the alignment bindings
//! describe a solve the same way in two and in three dimensions.
//!
//! These render the solver crate's termination reasons and this crate's refinement halts as
//! stable snake_case strings rather than leaking `Debug` formatting into the Python API, where it
//! would be neither documented nor safe to depend on.

/// How a single Levenberg-Marquardt solve ended, classified by whether its result can be used.
///
/// `"converged"` means a convergence criterion was met. `"unconverged"` means the solver ran out
/// of its evaluation budget, so the parameters are the best it found but convergence was never
/// demonstrated; the alignment is still valid geometry. The outcome types do not surface a solve
/// that broke down entirely (a failed initial solve is an error and a failed refinement round is
/// rolled back), so the `"failed"` arm below exists for match totality rather than as a value
/// callers should expect to see.
pub(crate) fn quality_str(q: engeom::common::SolveQuality) -> &'static str {
    match q {
        engeom::common::SolveQuality::Converged => "converged",
        engeom::common::SolveQuality::Unconverged => "unconverged",
        engeom::common::SolveQuality::Failed => "failed",
    }
}

/// How a single Levenberg-Marquardt solve terminated, as a stable snake_case string rather than
/// the solver crate's `Debug` formatting.
pub(crate) fn termination_str(t: &engeom::common::TerminationReason) -> String {
    use engeom::common::TerminationReason as T;
    match t {
        T::Converged { ftol, xtol } => match (ftol, xtol) {
            (true, true) => "converged(ftol,xtol)".to_string(),
            (true, false) => "converged(ftol)".to_string(),
            (false, true) => "converged(xtol)".to_string(),
            (false, false) => "converged".to_string(),
        },
        T::ResidualsZero => "residuals_zero".to_string(),
        T::Orthogonal => "orthogonal".to_string(),
        T::LostPatience => "lost_patience".to_string(),
        T::Numerical(s) => format!("numerical({s})"),
        T::User(s) => format!("user({s})"),
        T::NoImprovementPossible(s) => format!("no_improvement_possible({s})"),
        T::NoParameters => "no_parameters".to_string(),
        T::NoResiduals => "no_residuals".to_string(),
        T::WrongDimensions(s) => format!("wrong_dimensions({s})"),
    }
}

/// Why robust refinement stopped before completing every requested round.
pub(crate) fn halt_str(h: &engeom::common::RefinementHalt) -> String {
    match h {
        engeom::common::RefinementHalt::NoNoiseEstimate => "no_noise_estimate".to_string(),
        engeom::common::RefinementHalt::Underdetermined { weighted, params } => {
            format!("underdetermined({weighted} weighted points, {params} parameters)")
        }
        engeom::common::RefinementHalt::SolveFailed(t) => {
            format!("solve_failed({})", termination_str(t))
        }
    }
}
