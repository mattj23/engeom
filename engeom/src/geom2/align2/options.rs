//! Options controlling the behavior of a 2D alignment.

/// Options controlling how a 2D points-to-surface alignment weighs and filters its
/// correspondences.
///
/// Construct with [`AlignOptions2::default`] and override the fields you care about, so that
/// options added later don't break existing call sites.
///
/// # Robust estimation
///
/// By default, the alignment performs [`AlignOptions2::refinement_steps`] rounds of iteratively
/// reweighted least-squares after an initial unweighted solve, using MAGSAC++ noise-marginalized
/// weights. The weights are held fixed within each Levenberg-Marquardt solve so that the analytic
/// jacobian stays consistent with the residual, and are recomputed between solves. This is the
/// same structure the consensus framework in `common::consensus` uses for its refinement step.
///
/// Set `refinement_steps` to `0` to get a plain unweighted least-squares alignment.
#[derive(Clone, Copy, Debug)]
pub struct AlignOptions2<'a> {
    /// If the surface target can tell that a point does not project directly onto the target
    /// (such as when it projects past the end of an open curve or boundary), setting this flag
    /// weights such points at 0.0 to prevent their influence on the alignment.
    pub ignore_off_target: bool,

    /// The number of iteratively reweighted refinement rounds to perform after the initial
    /// unweighted solve. Zero disables robust weighting entirely.
    pub refinement_steps: usize,

    /// The MAGSAC++ upper noise bound.
    ///
    /// When `None`, this is estimated from the residuals of the initial unweighted solve using
    /// the median absolute deviation, scaled by the usual 1.4826 consistency constant. This is
    /// the closest the alignment gets to being free of a hand-tuned threshold, and is a good
    /// default when the noise scale isn't known ahead of time.
    ///
    /// The units depend on whether any uncertainty is in play. When `point_sigma` is supplied, or
    /// the target reports its own uncertainty, the residuals are dimensionless multiples of each
    /// point's combined standard deviation, so this is a count of sigmas and a value around `3.0`
    /// is a sensible explicit choice. Otherwise it is in the same units as the geometry.
    ///
    /// Note that points beyond roughly `3 * sigma_max` receive zero weight, so this is an upper
    /// bound on plausible noise rather than the noise scale itself.
    pub sigma_max: Option<f64>,

    /// Optional per-point measurement uncertainty, as a standard deviation in the same units as
    /// the geometry, with one entry per input point.
    ///
    /// When supplied, each point's residual is divided by its own standard deviation before being
    /// weighted, so that a point the sensor reports as three times noisier contributes one ninth
    /// as much to the squared cost. This is the statistically correct weighting for known
    /// heteroscedastic noise, and it also makes `sigma_max` dimensionless.
    ///
    /// If the target also reports uncertainty at the match position (see
    /// `AlignSurfMatch2::sigma`), the two combine in quadrature as `sqrt(test^2 + target^2)`,
    /// which is the standard deviation of the difference of two independent measurements.
    ///
    /// Every entry must be finite and strictly positive. For a mesh-derived point set this is
    /// what `Mesh3::point_stdev` provides in 3D.
    pub point_sigma: Option<&'a [f64]>,
}

impl Default for AlignOptions2<'_> {
    fn default() -> Self {
        Self {
            ignore_off_target: false,
            refinement_steps: 4,
            sigma_max: None,
            point_sigma: None,
        }
    }
}
