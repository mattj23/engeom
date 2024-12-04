//! This module contains the tools for performing geometric align2 on 2D shapes using the
//! Levenberg-Marquardt algorithm.

mod jacobian;
mod points_to_curve;
mod rc_params2;

use crate::geom2::Iso2;
use parry2d_f64::na::Vector3;

type T2Storage = Vector3<f64>;

pub use points_to_curve::points_to_curve;
pub use rc_params2::RcParams2;

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
