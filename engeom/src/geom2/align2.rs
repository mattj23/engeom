//! This module contains the tools for performing geometric align2 on 2D shapes using the
//! Levenberg-Marquardt algorithm.

mod jacobian;
mod params;
mod points_to_curve;
mod rc_params2;

use crate::geom2::Iso2;
use parry2d_f64::na::Vector3;

/// The storage for the three parameters of a 2D alignment problem, in the order tx, ty, rz.
pub type T2Storage = Vector3<f64>;

pub use self::params::*;
pub use points_to_curve::points_to_curve;
pub use rc_params2::RcParams2;

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

/// Produces a 2D transformation from 3 parameters.
pub fn iso2_from_param(p: &T2Storage) -> Iso2 {
    Iso2::translation(p.x, p.y) * Iso2::rotation(p.z)
}

/// Produces 3 parameters from a 2D transformation.
pub fn param_from_iso2(t: &Iso2) -> T2Storage {
    let v = t.translation.vector;
    let z = t.rotation.angle();
    T2Storage::new(v.x, v.y, z)
}
