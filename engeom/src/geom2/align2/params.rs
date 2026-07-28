//! This module contains the parameterization of the 2D alignment problem

use crate::common::PCoords;
use crate::geom2::align2::{AlignStorage2, Dof3};
use crate::na::{Matrix2, Translation2, UnitComplex};
use crate::{Iso2, Point2, Vector2};

// The skew-symmetric generator of a 2D rotation: d/dtheta R(theta) = R(theta) * SK2 = SK2 *
// R(theta), since 2D rotation matrices always commute with each other and with their own
// generator. This is why, unlike the 3D case, there is no gimbal correction to derive: a single
// rotation parameter has nothing "after" it in the rotation order that needs to be undone and
// reapplied.
const SK2: Matrix2<f64> = Matrix2::new(0.0, -1.0, 1.0, 0.0);

/// This struct holds the current values of the alignment problem, including the full
/// transformation from test entity space to target entity space, the vectors of translation
/// associated with the partial derivatives of the translation parameters, and the matrix for
/// calculating the vector of rotation associated with the partial derivative of the rotation
/// parameter. It is created from an [`AlignParams2`] struct based on the current values of the
/// alignment parameters.
///
/// During an alignment, the [`AlignParams2`]'s internal values will be changing as the solver
/// modifies it. This struct is generated to pre-calculate a number of values needed for working
/// with the alignment at that step in the alignment process.
#[derive(Clone, Debug)]
pub struct AlignValues2 {
    /// The transformation of the test entity(s) from their native, local coordinates to the
    /// current position in the target entity's coordinate system. This transformation is a
    /// composite of the active internal parameters and the working transformation.
    pub transform: Iso2,

    /// The transformation created by the tx, ty, and rz parameters.
    pub align: Iso2,

    /// The degrees of freedom that are active
    pub dof: Dof3,

    /// The direction vector for the partial derivative of tx
    pub dtx: Vector2,

    /// The direction vector for the partial derivative of ty
    pub dty: Vector2,

    /// The matrix to extract the partial derivative of rz and rotate it back to the target
    /// entity's coordinate system
    m_drz: Matrix2<f64>,

    /// The transformation from the target entity's coordinate system to the coordinate system
    /// where the rotation angle of the alignment transformation is applied.
    pre_rot: Iso2,
}

impl AlignValues2 {
    /// Given a point in the target entity's coordinate system, return the vector that describes
    /// the instantaneous vector of motion corresponding with the partial derivative of rz
    pub fn drz(&self, point: &impl PCoords<2>) -> Vector2 {
        self.m_drz * (self.pre_rot * Point2::from(point.coords())).coords()
    }
}

/// This is a mechanism of defining a local origin for the test entity.
#[derive(Clone, Copy, Debug)]
pub enum AlignOrigin2 {
    /// The local origin is the origin of the coordinate system that the test entity(s) geometry
    Origin,

    /// The local origin is centered at the given point, but the directions of translation are
    /// the same as the coordinate system of the test entity.
    Center(Point2),

    /// The local origin is defined by full transformation from the origin of the coordinate system
    /// with the test entity(s) geometry, allowing full control over the center of rotation and
    /// the directions of translation.
    Local(Iso2),
}

/// The [`AlignParams2`] struct holds the parameters being optimized in a 2D alignment problem,
/// expressed as a rotation angle and translation problem around an arbitrary local origin.
/// Additionally, a working transformation allows the test entity(s) to be aligned as if they were
/// in a different position without actually requiring them to be transformed ahead of time.
///
/// The parameterization consists of 3 numbers in an owned vector, representing tx, ty, and rz, in
/// that order, as a transformation at the local origin.
///
/// An alignment consists of a target (aka reference) entity and a test entity(s). The alignment
/// problem is to find the transformation that best aligns the test entity(s) to the target entity.
///
/// The entire alignment problem will take place in the target entity's coordinate system, such
/// that any components of the target entity which have positions will retain them wherever they
/// are. For the alignment, the target entity is in the world coordinate system.
///
/// The test entity(s) will be aligned to the geometry of the target entry by applying a 3-DOF
/// transformation ($A$). The transformation consists of two translations and one rotation. By
/// convention rotation is applied first around the origin, and then translation is applied along
/// the cardinal axes.
///
/// To allow for control over the center of rotation and the directions of translation, the
/// alignment can optionally define a local origin $L$ for the test entity(s). The local origin is
/// provided as a transformation from the origin of the coordinate system that the test entity(s)
/// geometry is defined in.
///
/// Lastly, there is an optional offset transformation $O$ that is applied after the alignment
/// transformation. This can be used to start the alignment process at a different position than
/// the current position of the test entity(s), or to counteract the effects of the local origin,
/// or a combination of the two.
///
/// The combined transformation from the test entity geometry to the target entity is:
///
/// $$ O * A * L^{-1} $$
///
#[derive(Clone, Debug)]
pub struct AlignParams2 {
    /// The local origin $L$, defined in the same space as the test entity's geometry. Leave this
    /// at the origin for simplicity, near the center of the test geometry to maximize numerical
    /// stability over rotations, or in a position/orientation to make use of the DOF constraints
    /// for a special case.
    pub local: Iso2,

    /// The degrees of freedom that are active during the alignment. Usually, all degrees of
    /// freedom will be active, but certain cases may require specific ones are locked. The local
    /// origin can be used in conjunction with the DOF constraints to control exactly how the test
    /// entity is allowed to move during alignment.
    pub dof: Dof3,

    /// The current working offset transformation $O$, which is the transformation applied after
    /// the alignment transformation.
    pub offset: Iso2,

    /// The storage for the three parameters
    storage: AlignStorage2,
}

impl AlignParams2 {
    /// Creates an `AlignParams2` with full control over its properties by separately specifying
    /// the local origin, working transformation, and degrees of freedom.
    ///
    /// # Arguments
    ///
    /// * `local`: The local origin $L$, defined in the same space as the test entity's geometry.
    ///   You can leave this as the world origin, pick a rotation center, or specify a full
    ///   transformation with origin and cardinal directions.
    /// * `offset`: An optional working offset transformation $O$, which is the transformation
    ///   applied after the alignment transformation.
    /// * `dof`: Optional constraint on the degrees of freedom. If `None` is provided, all degrees
    ///   of freedom will be active.
    ///
    /// returns: AlignParams2
    pub fn new(local: AlignOrigin2, offset: Option<Iso2>, dof: Option<Dof3>) -> Self {
        let local = match local {
            AlignOrigin2::Origin => Iso2::identity(),
            AlignOrigin2::Center(p) => Iso2::translation(p.x, p.y),
            AlignOrigin2::Local(t) => t,
        };

        let offset = offset.unwrap_or_else(Iso2::identity);
        let dof = dof.unwrap_or_else(Dof3::all);

        Self {
            local,
            dof,
            offset,
            storage: AlignStorage2::zeros(),
        }
    }

    /// Creates an `AlignParams2` which applies its transformation to the test entity(s) at a given
    /// local origin. The local origin and the working offset transformation will be identical.
    ///
    /// The physical interpretation of configuring the parameters this way is that the test
    /// geometry is transformed directly according to the local origin's position and orientation.
    /// That is, tx and ty mean transformation along the local origin's x and y axes. Rotation
    /// from rz means rotation around the local origin's centerpoint. Constraints applied to the
    /// degrees of freedom will refer to the local origin's axes.
    ///
    /// Use this configuration method when the test geometry is already in a good starting location,
    /// and you want to control exactly how the test geometry will move, such as if you want to
    /// apply DOF constraints in some arbitrary direction(s).
    ///
    /// # Arguments
    ///
    /// * `local`: the local origin $L$, defined in the same space as the test entity's geometry.
    /// * `dof`: Optional constraint on the degrees of freedom. If `None` is provided, all degrees
    ///   of freedom will be active.
    pub fn from_local(local: Iso2, dof: Option<Dof3>) -> Self {
        Self::new(AlignOrigin2::Local(local), Some(local), dof)
    }

    /// Creates an `AlignParams2` which applies its rotations to the test entity(s) around a given
    /// rotation center point. In this case the local origin $L$ will be created at the specified
    /// rotation center point but with the same cardinal directions as the world coordinate system.
    /// The working offset transformation will be the same as the local origin.
    ///
    /// Use this configuration method when the test geometry is already in a good starting location,
    /// but you want to provide a rotation center point instead of allowing rotation to happen
    /// around the world origin. This is important for the numerical stability of rotations if the
    /// test geometry is far from the world origin.
    ///
    /// # Arguments
    ///
    /// * `center`: The point around which the test entity(s) will be rotated.
    /// * `dof`: Optional constraint on the degrees of freedom. If `None` is provided, all degrees
    ///   of freedom will be active.
    pub fn from_center(center: Point2, dof: Option<Dof3>) -> Self {
        let local = Iso2::translation(center.x, center.y);
        Self::new(AlignOrigin2::Local(local), Some(local), dof)
    }

    /// Creates an `AlignParams2` with the local and working transformations set to the identity.
    /// Use this configuration method when the test geometry is already in a good starting location
    /// and close enough to the origin that you aren't worried about the numerical stability of
    /// rotations.
    ///
    /// # Arguments
    ///
    /// * `dof`: Optional constraint on the degrees of freedom. If `None` is provided, all degrees
    ///   of freedom will be active.
    ///
    /// returns: AlignParams2
    pub fn from_origin(dof: Option<Dof3>) -> Self {
        Self::new(AlignOrigin2::Origin, None, dof)
    }

    pub fn tx(&self) -> f64 {
        self.storage[0]
    }
    pub fn ty(&self) -> f64 {
        self.storage[1]
    }
    pub fn rz(&self) -> f64 {
        self.storage[2]
    }

    /// Computes the full transformation $O * A * L^{-1}$ that brings the test geometry into the
    /// target's coordinate system, including the local origin and the working offset. If the
    /// problem has converged, this is the final result of the alignment.
    ///
    /// This is a convenience wrapper over [`AlignParams2::compute_values`] for callers that only
    /// need the transform; it does the same work, so prefer `compute_values` if you also need the
    /// partial derivative directions.
    pub fn compute_transform(&self) -> Iso2 {
        self.compute_values().transform
    }

    /// Computes the current alignment values, including the full transform, the translation
    /// directions, and the rotation partial derivative matrix
    pub fn compute_values(&self) -> AlignValues2 {
        let align = align_from_storage(&self.storage);
        let transform = self.offset * align * self.local.inverse();

        let dtx = self.offset * Vector2::new(1.0, 0.0);
        let dty = self.offset * Vector2::new(0.0, 1.0);

        let pre_rot = align.translation.inverse() * self.offset.inverse();
        let post_rot = self.offset.rotation.to_rotation_matrix();

        AlignValues2 {
            transform,
            align,
            dof: self.dof,
            dtx,
            dty,
            m_drz: post_rot * SK2,
            pre_rot,
        }
    }

    pub fn storage(&self) -> AlignStorage2 {
        self.storage
    }

    fn enforce_constraint(&mut self) {
        if !self.dof.tx {
            self.storage[0] = 0.0;
        }
        if !self.dof.ty {
            self.storage[1] = 0.0;
        }
        if !self.dof.rz {
            self.storage[2] = 0.0;
        }
    }

    pub fn set_storage(&mut self, storage: AlignStorage2) {
        self.storage = storage;
        self.enforce_constraint();
    }

    pub fn set_index(&mut self, index: usize, value: f64) {
        self.storage[index] = value;
        self.enforce_constraint();
    }

    pub fn with_tx(&self, tx: f64) -> AlignParams2 {
        self.with_index(0, tx)
    }

    pub fn with_ty(&self, ty: f64) -> AlignParams2 {
        self.with_index(1, ty)
    }

    pub fn with_rz(&self, rz: f64) -> AlignParams2 {
        self.with_index(2, rz)
    }

    fn with_index(&self, index: usize, value: f64) -> AlignParams2 {
        let mut params = self.clone();
        params.set_index(index, value);
        params
    }

    pub fn with_storage(&self, storage: AlignStorage2) -> AlignParams2 {
        let mut params = self.clone();
        params.set_storage(storage);
        params
    }
}

/// Composes the `tx`, `ty`, `rz` parameters into the alignment isometry `A`.
fn align_from_storage(p: &AlignStorage2) -> Iso2 {
    Iso2::from_parts(Translation2::new(p.x, p.y), UnitComplex::new(p.z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry2;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    const ANGLE_EPSILON: f64 = 1e-8;
    const TRANS_EPSILON: f64 = 1e-8;

    // ============================================================================================
    // These first tests are to check some simple and obvious properties of AlignParams2, mirroring
    // the equivalent tests for AlignParams3 in geom3::align3::params.
    // ============================================================================================

    #[test]
    fn rotation_around_origin() {
        let params = AlignParams2::from_origin(None).with_rz(PI / 2.0);
        let test_point = Point2::new(1.0, 2.0);
        let expected = Point2::new(-2.0, 1.0);

        assert_relative_eq!(
            expected,
            params.compute_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    #[test]
    fn rotation_around_center() {
        let rc = Point2::new(1.0, 0.0);
        let params = AlignParams2::from_center(rc, None).with_rz(PI / 2.0);

        let test_point = Point2::new(2.0, 2.0);
        let expected = Point2::new(-1.0, 1.0);

        assert_relative_eq!(
            expected,
            params.compute_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    #[test]
    fn translations_are_in_local_origin_directions() {
        // A local origin whose axes are rotated 90 degrees from world: local +x points in world
        // +y, local +y points in world -x.
        let local = Iso2::rotation(PI / 2.0);
        let params = AlignParams2::from_local(local, None).with_tx(1.0);

        let test_point = Point2::new(1.0, 2.0);
        let expected = Point2::new(1.0, 3.0);

        assert_relative_eq!(
            expected,
            params.compute_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    #[test]
    fn builders_respect_locked_dof() {
        let dof = Dof3::new(false, true, true);
        let params = AlignParams2::from_origin(Some(dof));

        assert_eq!(params.with_tx(1.0).tx(), 0.0);
        assert_eq!(params.with_ty(1.0).ty(), 1.0);
        assert_eq!(params.with_rz(1.0).rz(), 1.0);

        let all = params.with_storage(AlignStorage2::new(1.0, 2.0, 3.0));
        assert_eq!(all.tx(), 0.0);
        assert_eq!(all.ty(), 2.0);
        assert_eq!(all.rz(), 3.0);
    }

    // ============================================================================================
    // These next tests are to check that the current translation directions match the direction
    // that they actually move entities.
    // ============================================================================================

    #[test]
    fn partials_of_translations_at_zero() {
        let params = AlignParams2::from_origin(None);
        let current = params.compute_values();

        assert_relative_eq!(current.dtx, Vector2::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(current.dty, Vector2::y_axis(), epsilon = 1e-12);
    }

    #[test]
    fn partials_of_translations_with_rotation() {
        let params = AlignParams2::from_origin(None).with_rz(0.3);
        let test_point = Point2::new(1.0, 2.0);

        let exp_x = finite_diff(&params, &test_point, 0);
        let exp_y = finite_diff(&params, &test_point, 1);
        let c = params.compute_values();

        // Because translations are applied after rotation, the direction vectors don't change
        // when the rotation parameter is nonzero.
        assert_relative_eq!(c.dtx, Vector2::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(c.dty, Vector2::y_axis(), epsilon = 1e-12);

        assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-8);
        assert_relative_eq!(exp_y, c.dty, epsilon = 1e-8);
    }

    #[test]
    fn stress_partials_translations() {
        let mut rg = RandomGeometry2::new();
        for _ in 0..1000 {
            let local = rg.iso2(10.0);
            let working = rg.iso2(10.0);
            let params = AlignParams2::new(AlignOrigin2::Local(local), Some(working), None)
                .with_storage(rg.vector(PI));

            let test_point = rg.point(10.0);
            let exp_x = finite_diff(&params, &test_point, 0);
            let exp_y = finite_diff(&params, &test_point, 1);

            let c = params.compute_values();
            assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-6);
            assert_relative_eq!(exp_y, c.dty, epsilon = 1e-6);
        }
    }

    // ============================================================================================
    // Partial derivative of the rotation parameter. Unlike the 3D case there's only one rotation
    // axis and no gimbal effect, so there's no analogue of `euler_partials` to verify separately;
    // the stress test below exercises `AlignValues2::drz` directly against a finite-difference
    // estimate over random local origins, offsets, and parameter values.
    // ============================================================================================

    #[test]
    fn partials_of_rotation_at_zero() {
        let params = AlignParams2::from_origin(None);
        let current = params.compute_values();
        let p_local = Point2::new(1.0, 2.0);
        let p = current.transform * p_local;

        let exp_rz = finite_diff(&params, &p_local, 2);
        assert_relative_eq!(exp_rz, current.drz(&p), epsilon = 1e-6);
    }

    #[test]
    fn stress_partials_rotation() {
        let mut rg = RandomGeometry2::new();
        for _ in 0..1000 {
            let local = rg.iso2(10.0);
            let working = rg.iso2(10.0);
            let params = AlignParams2::new(AlignOrigin2::Local(local), Some(working), None)
                .with_storage(rg.vector(PI));
            let c = params.compute_values();

            // The test point in the test entity's space
            let test_point = rg.point(10.0);
            let exp_rz = finite_diff(&params, &test_point, 2);

            // The test point in the target entity's space
            let target_point = c.transform * test_point;

            // As in the 3D stress test, the finite-difference estimate and the analytic value can
            // diverge into the fourth decimal place for points far from the rotation center; this
            // is a plausibility check, not a tight tolerance.
            assert_relative_eq!(exp_rz, c.drz(&target_point), epsilon = 1e-3);
        }
    }

    // ============================================================================================
    // Test support methods
    // ============================================================================================

    /// This function computes the vector of motion of a point in the original coordinate system
    /// of the test entity for a finitely approximated infinitesimal change in one of the
    /// parameters. The parameters by index are: (0: tx, 1: ty, 2: rz)
    fn finite_diff(params: &AlignParams2, point: &Point2, index: usize) -> Vector2 {
        let mut w0 = params.clone();
        let mut w1 = params.clone();
        let stored = params.storage();

        let eps = if index < 2 {
            TRANS_EPSILON
        } else {
            ANGLE_EPSILON
        };

        w0.set_index(index, stored[index] - eps);
        w1.set_index(index, stored[index] + eps);

        let p0 = w0.compute_values().transform * point;
        let p1 = w1.compute_values().transform * point;

        (p1 - p0) / (2.0 * eps)
    }
}
