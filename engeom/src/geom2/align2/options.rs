//! Options controlling the behavior of a 2D alignment.

/// Options controlling how a 2D points-to-surface alignment weighs and filters its correspondences.
///
/// Construct with [`AlignOptions2::default`] and override the fields you care about, so that
/// options added later don't break existing call sites.
#[derive(Clone, Copy, Debug, Default)]
pub struct AlignOptions2 {
    /// If the surface target can tell that a point does not project directly onto the target
    /// (such as when it projects past the end of an open curve or boundary), setting this flag
    /// weights such points at 0.0 to prevent their influence on the alignment.
    pub ignore_off: bool,
}
