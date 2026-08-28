//! This module contains the tools for aligning 3D geometry using the Levenberg-Marquardt
//! algorithm.
//!
//! # Structure
//!
//! An alignment has three pieces:
//!
//! - [`AlignParams3`] holds the parameters being optimized (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`)
//!   and expresses them as a transformation about an arbitrary local origin, with an optional
//!   working offset. This is what gives the caller control over the center of rotation, the
//!   directions the translation parameters act along, and which degrees of freedom ([`Dof6`]) are
//!   free to move at all.
//! - [`SurfaceTarget3`] is the stationary entity being aligned to. It is implemented for `Mesh3`,
//!   `ExtrudedBoundary3`, and `RevolvedBoundary3`, and reports an [`AlignSurfMatch3`] for any
//!   query point: the closest position, its normal, whether the projection actually landed on the
//!   target's interior, and optionally the target's own measurement uncertainty there.
//! - [`points_to_surface3`] runs the solver, with behavior controlled by [`AlignOptions3`].
//!
//! # Reporting
//!
//! The solver returns an [`crate::geom3::AlignOutcome3`], which carries the alignment plus a
//! record of how every solve that contributed to it terminated. An `Err` is reserved for the case
//! where there is no answer at all: rejected arguments, or an initial solve that broke down.
//!
//! Everything short of that is reported rather than raised, because there is still a real
//! alignment to return. In particular, a solve that exhausts its evaluation budget leaves behind
//! the best parameters it found, which is a usable result whose convergence simply was not proven
//! (see [`crate::common::SolveQuality`]). That is a routine outcome here rather than a failure:
//! the correspondences are re-established every time the parameters change, so the objective is
//! only piecewise smooth, and near the solution a point close to an edge or a corner can flip
//! between two matches indefinitely without the convergence criteria ever being met. Similarly, a
//! refinement round that breaks down is rolled back to the previous round's result and the reason
//! recorded on the outcome, since refinement is an improvement on an alignment that was already
//! usable.
//!
//! # Robustness
//!
//! The alignment is robust by default. An initial unweighted solve is followed by several rounds
//! of iteratively reweighted least-squares using MAGSAC++ noise-marginalized weights, with the
//! weights held fixed inside each solve so the analytic jacobian stays consistent with the
//! residual it differentiates. The noise bound can be supplied explicitly or estimated from the
//! data via the median absolute deviation. Measurement uncertainty on either the test points or
//! the target combines in quadrature and normalizes the residuals into units of sigma;
//! `Mesh3::point_stdev` is the data source for both sides.
//!
//! # Relationship to `align2`
//!
//! This module and [`crate::geom2::align2`] are deliberately structural mirrors rather than a
//! shared generic: 3D needs three Euler angles whose partial derivatives require a gimbal
//! correction that 2D has no analogue for, and its parameter storage is twice the size.
//!
//! Neither module is the authority. The 2D module was written first and this one was brought
//! into line with it, but the multi-body work happened here and went back the other way. A
//! change to either should be considered for the other, and the file layout is kept parallel so
//! that comparison is easy to make.

mod cloud;
mod information;
pub mod jacobian;
mod mesh;
mod multi_mesh;
mod multi_params;
mod options;
mod params;
mod points_to_surface;

use crate::UnitVec3;
use crate::common::{PCoords, SPCoords};
use crate::geom3::{Iso3, Point3, Vector3};
use crate::na::{SVector, Unit};
use parry3d_f64::na::{Translation3, UnitQuaternion, Vector6};

/// The storage for the six parameters of a 3D alignment problem, in the order tx, ty, tz, rx, ry,
/// rz.
pub type AlignStorage3 = Vector6<f64>;

pub use self::cloud::CloudTarget3;
pub use self::information::*;
pub use self::mesh::*;
pub use self::multi_mesh::{
    MultiAlignOptions3, MultiMeshAlignPoint, multi_mesh_adjustment,
    multi_mesh_adjustment_with_points,
};
pub use self::multi_params::MultiAlignParams3;
pub use self::options::*;
pub use self::params::*;
pub use self::points_to_surface::*;

/// The result of projecting a single point onto a [`SurfaceTarget3`], used as the correspondence
/// for that point during an alignment.
#[derive(Debug, Clone)]
pub struct AlignSurfMatch3 {
    /// The closest point on the target to the query point
    pub point: Point3,

    /// The outward-facing normal of the target at `point`
    pub normal: UnitVec3,

    /// Whether the closest point actually lies on the target's interior, as opposed to having
    /// clamped to an edge or an end of a bounded target. See the target's own `find_align_match`
    /// implementation for the exact clamping rule.
    pub is_on: bool,

    /// A scalar weight for the correspondence, in `[0, 1]`, independent of `is_on`. A target may
    /// use this to de-weight regions of itself that it considers less reliable.
    ///
    /// This is a statement of intent ("care about this correspondence less"), distinct from
    /// [`AlignSurfMatch3::sigma`], which is a statement about measurement noise.
    pub weight: f64,

    /// The measurement uncertainty of the target at `point`, as a standard deviation in the units
    /// of the geometry. Zero means the target is treated as exact, which is the default and is
    /// correct for nominal/theoretical geometry.
    ///
    /// A target built from measured data should report the uncertainty interpolated to `point`,
    /// since the match rarely lands exactly on a vertex. `Mesh3::point_stdev` is the data source
    /// for a measured mesh. The alignment combines this with the test point's own uncertainty in
    /// quadrature, `sqrt(test^2 + target^2)`, which is the variance of the difference of two
    /// independent measurements.
    ///
    /// This is treated as **isotropic**. Real scanner uncertainty is usually one-dimensional
    /// (depth along the sensor axis), and the statistically correct contribution to a residual
    /// measured along direction `d` would be `sigma * |u . d|` for an uncertainty axis `u`.
    /// Nothing currently records `u`, so the isotropic treatment stands in for it. Note that the
    /// approximation only ever under-trusts a point: on a surface at grazing incidence, depth
    /// noise displaces the point along the surface rather than through it, so its true
    /// normal-direction uncertainty is smaller than the scalar suggests.
    pub sigma: f64,
}

impl AlignSurfMatch3 {
    pub fn new(point: Point3, normal: UnitVec3, is_on: bool, weight: f64) -> Self {
        Self {
            point,
            normal,
            is_on,
            weight,
            sigma: 0.0,
        }
    }

    /// Returns a copy of this match carrying the given target-side measurement uncertainty. See
    /// [`AlignSurfMatch3::sigma`] for the semantics and the isotropy caveat.
    pub fn with_sigma(&self, sigma: f64) -> Self {
        Self {
            sigma,
            ..self.clone()
        }
    }
}

impl PCoords<3> for AlignSurfMatch3 {
    fn coords(&self) -> SVector<f64, 3> {
        self.point.coords()
    }
}

impl SPCoords<3> for AlignSurfMatch3 {
    fn normal(&self) -> Unit<SVector<f64, 3>> {
        self.normal
    }
}

impl Default for AlignSurfMatch3 {
    fn default() -> Self {
        Self {
            point: Point3::origin(),
            normal: Vector3::x_axis(),
            is_on: false,
            weight: 0.0,
            sigma: 0.0,
        }
    }
}

/// A stationary 3D entity that a set of points can be aligned to, by projecting each point onto
/// its closest position on the target.
///
/// This takes `&Point3` rather than `&impl PCoords<3>` so that the trait remains object-safe,
/// allowing `Box<dyn SurfaceTarget3>` / `Vec<Box<dyn SurfaceTarget3>>` for cases (such as
/// multi-entity alignment) where the set of targets isn't known at compile time. Every call site
/// already holds a concrete, already-transformed point, so nothing is lost in generality.
///
/// A target derived from measured rather than nominal geometry should populate
/// [`AlignSurfMatch3::sigma`] via [`AlignSurfMatch3::with_sigma`], interpolating its own
/// uncertainty to the match position.
pub trait SurfaceTarget3: Sync + Send {
    fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3;
}

/// A struct that handles constraints on degrees of freedom in R^3 space. Each dimension is
/// represented by a bool which specifies if the degree of freedom is _active_.
#[derive(Clone, Copy, Debug)]
pub struct Dof6 {
    pub tx: bool,
    pub ty: bool,
    pub tz: bool,
    pub rx: bool,
    pub ry: bool,
    pub rz: bool,
}

impl Dof6 {
    pub fn new(tx: bool, ty: bool, tz: bool, rx: bool, ry: bool, rz: bool) -> Self {
        Self {
            tx,
            ty,
            tz,
            rx,
            ry,
            rz,
        }
    }

    /// Returns a new Dof3 with all degrees of freedom active.
    pub fn all() -> Self {
        Self {
            tx: true,
            ty: true,
            tz: true,
            rx: true,
            ry: true,
            rz: true,
        }
    }
}

impl Default for Dof6 {
    fn default() -> Self {
        Self::all()
    }
}

#[derive(Clone, Copy, Debug)]
pub enum SampleMode {
    All,
    Random(usize),
    Poisson(f64),
}

pub fn iso3_from_param(p: &AlignStorage3) -> Iso3 {
    Iso3::from_parts(
        Translation3::new(p.x, p.y, p.z),
        UnitQuaternion::from_euler_angles(p.w, p.a, p.b),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    use approx::assert_relative_eq;

    #[test]
    fn iso3_tx() {
        let storage = AlignStorage3::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(1.0, 0.0, 0.0);
        let p2 = t * p;
        assert_relative_eq!(p2.x, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_ty() {
        let storage = AlignStorage3::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 1.0, 0.0);
        let p2 = t * p;
        assert_relative_eq!(p2.y, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_tz() {
        let storage = AlignStorage3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 0.0, 1.0);
        let p2 = t * p;
        assert_relative_eq!(p2.z, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_rx() {
        let storage = AlignStorage3::new(0.0, 0.0, 0.0, FRAC_PI_2, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 1.0, 0.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::x_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }

    #[test]
    fn iso3_ry() {
        let storage = AlignStorage3::new(0.0, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 0.0, 1.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::y_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }

    #[test]
    fn iso3_rz() {
        let storage = AlignStorage3::new(0.0, 0.0, 0.0, 0.0, 0.0, FRAC_PI_2);
        let t = iso3_from_param(&storage);
        let p = Point3::new(1.0, 0.0, 0.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::z_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }
}
