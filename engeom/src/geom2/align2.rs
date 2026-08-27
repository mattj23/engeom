//! This module contains the tools for aligning 2D geometry using the Levenberg-Marquardt
//! algorithm.
//!
//! # Structure
//!
//! An alignment has three pieces:
//!
//! - [`AlignParams2`] holds the parameters being optimized (`tx`, `ty`, `rz`) and expresses them
//!   as a transformation about an arbitrary local origin, with an optional working offset. This
//!   is what gives the caller control over the center of rotation, the directions the translation
//!   parameters act along, and which degrees of freedom ([`Dof3`]) are free to move at all.
//! - [`SurfaceTarget2`] is the stationary entity being aligned to. It is implemented for `Curve2`
//!   and `Boundary2`, and reports an [`AlignSurfMatch2`] for any query point: the closest
//!   position, its normal, whether the projection actually landed on the target's interior, and
//!   optionally the target's own measurement uncertainty there.
//! - [`points_to_surface2`] runs the solver, with behavior controlled by [`AlignOptions2`].
//!
//! # Reporting
//!
//! The solver returns an [`crate::geom2::AlignOutcome2`], which carries the alignment plus a
//! record of how every solve that contributed to it terminated. An `Err` is reserved for the case
//! where there is no answer at all: rejected arguments, or an initial solve that broke down.
//!
//! Everything short of that is reported rather than raised, because there is still a real
//! alignment to return. In particular, a solve that exhausts its evaluation budget leaves behind
//! the best parameters it found, which is a usable result whose convergence simply was not proven
//! (see [`crate::common::SolveQuality`]). That is a routine outcome here rather than a failure:
//! the correspondences are re-established every time the parameters change, so the objective is
//! only piecewise smooth, and near the solution a point close to a corner can flip between two
//! matches indefinitely without the convergence criteria ever being met. Similarly, a refinement
//! round that breaks down is rolled back to the previous round's result and the reason recorded on
//! the outcome, since refinement is an improvement on an alignment that was already usable.
//!
//! # Robustness
//!
//! The alignment is robust by default. An initial unweighted solve is followed by several rounds
//! of iteratively reweighted least-squares using MAGSAC++ noise-marginalized weights, with the
//! weights held fixed inside each solve so the analytic jacobian stays consistent with the
//! residual it differentiates. The noise bound can be supplied explicitly or estimated from the
//! data via the median absolute deviation.
//!
//! On the single-body path, measurement uncertainty on either the test points or the target
//! combines in quadrature and normalizes the residuals into units of sigma. The multi-curve
//! path carries no uncertainty at all, because nothing in 2D produces it: its residuals stay in
//! the units of the geometry, and so does its noise bound.
//!
//! # Relationship to `align3`
//!
//! This module and the 3D machinery in [`crate::geom3::align3`] are deliberately structural
//! mirrors rather than a shared generic: 2D needs only a single rotation angle, whose partial
//! derivative is a plain 90-degree turn with none of the Euler-angle gimbal correction that 3D
//! requires, and its parameter storage is half the size.
//!
//! Neither module is the authority. This one was written first and the 3D module was brought
//! into line with it, but the multi-body work happened in 3D and came back the other way. A
//! change to either should be considered for the other, and the file layout is kept parallel so
//! that comparison is easy to make.

mod curve;
pub mod jacobian;
mod multi_curve;
mod multi_params;
mod options;
mod params;
mod points_to_surface;
mod target;

use parry2d_f64::na::Vector3;

/// The storage for the three parameters of a 2D alignment problem, in the order tx, ty, rz.
pub type AlignStorage2 = Vector3<f64>;

pub use self::curve::*;
pub use self::multi_curve::{
    MulCurveAlignPoint, MultiAlignOptions2, multi_curve_adjustment,
    multi_curve_adjustment_with_points,
};
pub use self::multi_params::MultiAlignParams2;
pub use self::options::*;
pub use self::params::*;
pub use self::points_to_surface::points_to_surface2;
pub use self::target::*;

/// A struct that handles constraints on degrees of freedom in R^2 space. Each dimension is
/// represented by a bool which specifies if the degree of freedom is _active_.
#[derive(Clone, Copy, Debug)]
pub struct Dof3 {
    pub tx: bool,
    pub ty: bool,
    pub rz: bool,
}

impl Dof3 {
    pub fn new(tx: bool, ty: bool, rz: bool) -> Self {
        Self { tx, ty, rz }
    }

    /// Returns a new Dof3 with all degrees of freedom active.
    pub fn all() -> Self {
        Self {
            tx: true,
            ty: true,
            rz: true,
        }
    }
}

impl Default for Dof3 {
    fn default() -> Self {
        Self::all()
    }
}
