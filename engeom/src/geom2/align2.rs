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
//! # Robustness
//!
//! The alignment is robust by default. An initial unweighted solve is followed by several rounds
//! of iteratively reweighted least-squares using MAGSAC++ noise-marginalized weights, with the
//! weights held fixed inside each solve so the analytic jacobian stays consistent with the
//! residual it differentiates. The noise bound can be supplied explicitly or estimated from the
//! data via the median absolute deviation. Measurement uncertainty on either the test points or
//! the target combines in quadrature and normalizes the residuals into units of sigma.
//!
//! # Relationship to `align3`
//!
//! This module is the reference implementation for `engeom`'s alignment design, and the 3D
//! machinery in [`crate::geom3::align3`] is being brought into line with it. The two are
//! deliberately structural mirrors rather than a shared generic: 2D needs only a single rotation
//! angle, whose partial derivative is a plain 90-degree turn with none of the Euler-angle gimbal
//! correction that 3D requires, and its parameter storage is half the size.

mod jacobian;
mod options;
mod params;
mod points_to_surface;
mod target;

use parry2d_f64::na::Vector3;

/// The storage for the three parameters of a 2D alignment problem, in the order tx, ty, rz.
pub type AlignStorage2 = Vector3<f64>;

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
