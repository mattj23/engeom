//! This module contains the parameterization of the alignment problem

use crate::common::PCoords;
use crate::geom3::align3::{Dof6, T3Storage, iso3_from_param};
use crate::na::{Matrix3, UnitQuaternion};
use crate::{Iso3, Point3, Vector3};

// Skew-symmetric matrices
const SK_X: Matrix3<f64> = Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
const SK_Y: Matrix3<f64> = Matrix3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
const SK_Z: Matrix3<f64> = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

/// This struct holds the current values of the alignment problem, including the full
/// transformation from test entity space to target entity space, the vectors of translation
/// associated with the partial derivatives of the translation parameters, and functions for
/// calculating the vectors of rotation associated with the partial derivatives of the rotation
/// parameters. It is created from an [`AlignParams3`] struct based on the current values of
/// the alignment parameters.
///
/// During an alignment, the [`AlignParams3`]'s internal values will be changing as the solver
/// modifies it. This struct is generated to pre-calculate a number of values needed for working
/// with the alignment at that step in the alignment process.
#[derive(Clone, Debug)]
pub struct AlignValues3 {
    /// The transformation of the test entity(s) from their native, local coordinates to the
    /// current position in the target entity's coordinate system. This transformation is a
    /// composite of the active internal parameters and the working transformation.
    pub transform: Iso3,

    /// The transformation created by the tx, ty, tz, rx, ry, and rz parameters.
    pub align: Iso3,

    /// The degrees of freedom that are active
    pub dof: Dof6,

    /// The direction vector for the partial derivative of tx
    pub dtx: Vector3,

    /// The direction vector for the partial derivative of ty
    pub dty: Vector3,

    /// The direction vector for the partial derivative of tz
    pub dtz: Vector3,

    /// The matrix to extract the partial derivative of rx in the euler angle coordinate system
    /// and rotate it back to the target entity's coordinate system
    m_drx: Matrix3<f64>,

    /// The matrix to extract the partial derivative of ry in the euler angle coordinate system
    /// and rotate it back to the target entity's coordinate system
    m_dry: Matrix3<f64>,

    /// The matrix to extract the partial derivative of rz in the euler angle coordinate system
    /// and rotate it back to the target entity's coordinate system
    m_drz: Matrix3<f64>,

    /// The transformation from the target entity's coordinate system to the coordinate system
    /// where the euler angles of the alignment transformation are applied.
    pre_rot: Iso3,
}

impl AlignValues3 {
    /// Given a point in the target entity's coordinate system, return the vector that describes
    /// the instantaneous vector of motion corresponding with the partial derivative of rx
    pub fn drx(&self, point: &impl PCoords<3>) -> Vector3 {
        self.m_drx * (self.pre_rot * Point3::from(point.coords())).coords()
    }

    /// Given a point in the target entity's coordinate system, return the vector that describes
    /// the instantaneous vector of motion corresponding with the partial derivative of ry
    pub fn dry(&self, point: &impl PCoords<3>) -> Vector3 {
        self.m_dry * (self.pre_rot * Point3::from(point.coords())).coords()
    }

    /// Given a point in the target entity's coordinate system, return the vector that describes
    /// the instantaneous vector of motion corresponding with the partial derivative of rz
    pub fn drz(&self, point: &impl PCoords<3>) -> Vector3 {
        self.m_drz * (self.pre_rot * Point3::from(point.coords())).coords()
    }
}

/// This is a mechanism of defining a local origin for the test entity.
pub enum AlignOrigin {
    /// The local origin is the origin of the coordinate system that the test entity(s) geometry
    Origin,

    /// The local origin is centered at the given point, but the directions of translation are
    /// the same as the coordinate system of the test entity.
    Center(Point3),

    /// The local origin is defined by full transformation from the origin of the coordinate system
    /// with the test entity(s) geometry, allowing full control over the center of rotation and
    /// the directions of translation.
    Local(Iso3),
}

/// The [`AlignParams3`] struct holds the parameters being optimized in a 3D alignment problem,
/// expressed as an euler angle rotation and translation problem around an arbitrary local origin.
/// Additionally, a working transformation allows the test entity(s) to be aligned as if they were
/// in a different position without actually requiring them to be transformed ahead of time.
///
/// The parameterization consists of 6 numbers in an owned vector, representing tx, ty, tz, rx, ry,
/// and rz, in that order, as a transformation at the local origin.
///
/// An alignment consists of a target (aka reference) entity and a test entity(s). The alignment
/// problem is to find the transformation that best aligns the test entity(s) to the target entity.
///
/// The entire alignment problem will take place in the target entity's coordinate system, such
/// that any components of the target entity which have positions will retain them wherever they
/// are. For the alignment, the target entity is in the world coordinate system.
///
/// The test entity(s) will be aligned to the geometry of the target entry by applying a 6-DOF
/// transformation ($A$). The transformation consists of three translations and three rotations. By
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
pub struct AlignParams3 {
    /// The local origin $L$, defined in the same space as the test entity's geometry. Leave this
    /// at the origin for simplicity, near the center of the test geometry to maximize numerical
    /// stability over rotations, or in a position/orientation to make use of the DOF constraints
    /// for a special case.
    pub local: Iso3,

    /// The degrees of freedom that are active during the alignment. Usually, all degrees of
    /// freedom will be active, but certain cases may require specific ones are locked. The local
    /// origin can be used in conjunction with the DOF constraints to control exactly how the test
    /// entity is allowed to move during alignment.
    pub dof: Dof6,

    /// The current working offset transformation $O$, which is the transformation applied after
    /// the alignment transformation.
    pub offset: Iso3,

    /// The storage for the six parameters
    storage: T3Storage,
}

impl AlignParams3 {
    /// Creates an `AlignParams3` with full control over its properties by separately specifying
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
    /// returns: AlignParams3
    pub fn new(local: AlignOrigin, offset: Option<Iso3>, dof: Option<Dof6>) -> Self {
        let local = match local {
            AlignOrigin::Origin => Iso3::identity(),
            AlignOrigin::Center(p) => Iso3::translation(p.x, p.y, p.z),
            AlignOrigin::Local(t) => t,
        };

        let working = offset.unwrap_or_else(|| Iso3::identity());
        let dof = dof.unwrap_or_else(|| Dof6::all());

        Self {
            local,
            dof,
            offset: working,
            storage: T3Storage::default(),
        }
    }

    /// Creates an `AlignParams3` which applies its transformation to the test entity(s) at a given
    /// local origin. The local origin and the working offset transformation will be identical.
    ///
    /// The physical interpretation of configuring the parameters this way is that the test
    /// geometry is transformed directly according to the local origin's position and orientation.
    /// That is, tx, ty, and tz mean transformation along the local origin's x, y, and z axes.
    /// Rotation from rx, ry, and rz mean rotation around the local origin's centerpoint by the
    /// local origin's x, y, and z axes. Constraints applied to the degrees of freedom will refer
    /// to the local origin's axes.
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
    pub fn new_at_local(local: Iso3, dof: Option<Dof6>) -> Self {
        let working = local.clone();
        Self::new(AlignOrigin::Local(local), Some(working), dof)
    }

    /// Creates an `AlignParams3` which applies its rotations to the test entity(s) around a given
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
    pub fn new_at_center(center: Point3, dof: Option<Dof6>) -> Self {
        let local = Iso3::translation(center.x, center.y, center.z);
        let working = local.clone();
        Self::new(AlignOrigin::Local(local), Some(working), dof)
    }

    /// Creates an `AlignParams3` with the local and working transformations set to the identity.
    /// Use this configuration method when the test geometry is already in a good starting location
    /// and close enough to the origin that you aren't worried about the numerical stability of
    /// rotations.
    ///
    /// # Arguments
    ///
    /// * `dof`: Optional constraint on the degrees of freedom. If `None` is provided, all degrees
    ///   of freedom will be active.
    ///
    /// returns: AlignParams3
    pub fn new_at_origin(dof: Option<Dof6>) -> Self {
        Self::new(AlignOrigin::Origin, None, dof)
    }

    pub fn tx(&self) -> f64 {
        self.storage[0]
    }
    pub fn ty(&self) -> f64 {
        self.storage[1]
    }
    pub fn tz(&self) -> f64 {
        self.storage[2]
    }
    pub fn rx(&self) -> f64 {
        self.storage[3]
    }
    pub fn ry(&self) -> f64 {
        self.storage[4]
    }
    pub fn rz(&self) -> f64 {
        self.storage[5]
    }

    /// Get the stored isometry without the working transformation. If the problem has converged,
    /// this is the final result.
    pub fn final_result(&self) -> Iso3 {
        self.current_values().transform
    }

    /// Calculate the current alignment values, including the full transform, the translation
    /// directions, and the rotation partial derivative matrices
    pub fn current_values(&self) -> AlignValues3 {
        let align = iso3_from_param(&self.storage);
        let transform = self.offset * align * self.local.inverse();

        let dtx = self.offset * Vector3::new(1.0, 0.0, 0.0);
        let dty = self.offset * Vector3::new(0.0, 1.0, 0.0);
        let dtz = self.offset * Vector3::new(0.0, 0.0, 1.0);

        let (m_drx, m_dry, m_drz) = euler_partials(self.ry(), self.rz());

        let pre_rot = align.translation.inverse() * self.offset.inverse();
        let post_rot = self.offset.rotation.to_rotation_matrix();

        AlignValues3 {
            transform,
            align,
            dof: self.dof,
            dtx,
            dty,
            dtz,
            m_drx: post_rot * m_drx,
            m_dry: post_rot * m_dry,
            m_drz: post_rot * m_drz,
            pre_rot,
        }
    }

    pub fn get_storage(&self) -> T3Storage {
        self.storage
    }

    fn enforce_constraint(&mut self) {
        if !self.dof.tx {
            self.storage[0] = 0.0;
        }
        if !self.dof.ty {
            self.storage[1] = 0.0;
        }
        if !self.dof.tz {
            self.storage[2] = 0.0;
        }
        if !self.dof.rx {
            self.storage[3] = 0.0;
        }
        if !self.dof.ry {
            self.storage[4] = 0.0;
        }
        if !self.dof.rz {
            self.storage[5] = 0.0;
        }
    }

    pub fn set_storage(&mut self, storage: T3Storage) {
        self.storage = storage;
        self.enforce_constraint();
    }

    pub fn set_index(&mut self, index: usize, value: f64) {
        self.storage[index] = value;
        self.enforce_constraint();
    }

    pub(super) fn with_tx(&self, tx: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[0] = tx;
        params
    }

    pub(super) fn with_ty(&self, ty: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[1] = ty;
        params
    }

    pub(super) fn with_tz(&self, tz: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[2] = tz;
        params
    }

    pub(super) fn with_rx(&self, rx: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[3] = rx;
        params
    }

    pub(super) fn with_ry(&self, ry: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[4] = ry;
        params
    }

    pub(super) fn with_rz(&self, rz: f64) -> AlignParams3 {
        let mut params = self.clone();
        params.storage[5] = rz;
        params
    }

    pub fn with_storage(&self, storage: T3Storage) -> AlignParams3 {
        let mut params = self.clone();
        params.storage = storage;
        params
    }
}

/// This function finds the three matrices that, when multiplied by the coordinate vector of a
/// point in space, find the instantaneous velocity vector of that point for an infinitesimal
/// change in the rx, ry, and rz parameters, given the current values of the y and z parameters.
///
/// This assumes that the center of rotation is the origin.
///
/// The skew symmetric matrices are the correct matrices to find the velocity vectors when the
/// existing parameters are all at zero, as they correspond with instantaneous rotations around
/// their respective axes. However, to account for the gimbal effect, we will need to be able to
/// undo the effects of the later rotations for any rotation in the order.
///
/// Since the nalgebra library uses the Euler order X-Y-Z, the Z rotation has nothing following it
/// and can be measured directly with the skew symmetric Z matrix. For the Y rotation, we will need
/// to first undo Z, then apply the skew symmetric Y matrix, and then re-apply the Z rotation. For
/// X we will need to undo then redo both Z and Y. This is why this function requires the ry and
/// rz parameters, but not the rx parameter.
///
/// # Arguments
///
/// * `ry`: the current rotation around the Y axis
/// * `rz`: the current rotation around the Z axis
///
/// returns: (Matrix<f64, Const<3>, Const<3>, ArrayStorage<f64, 3, 3>>, Matrix<f64, Const<3>, Const<3>, ArrayStorage<f64, 3, 3>>, Matrix<f64, Const<3>, Const<3>, ArrayStorage<f64, 3, 3>>)
fn euler_partials(ry: f64, rz: f64) -> (Matrix3<f64>, Matrix3<f64>, Matrix3<f64>) {
    let y = UnitQuaternion::from_euler_angles(0.0, ry, 0.0).to_rotation_matrix();
    let z = UnitQuaternion::from_euler_angles(0.0, 0.0, rz).to_rotation_matrix();

    // Z will always be the skew symmetric matrix, because it is the last rotation applied,
    // so it will always produce a uniform rotation around Z

    // To find the matrix for Y, we need to undo the Z rotation, apply the skew symmetric matrix,
    // then re-rotate back to the original Z
    let dy = z * SK_Y * z.transpose();

    // To find the matrix for X, we need to undo the Z rotation, then the Y rotation, apply the
    // skew symmetric matrix, then re-rotate the Y, and finally re-rotate the Z.
    let dx = z * y * SK_X * y.transpose() * z.transpose();

    (dx, dy, SK_Z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::IsoExtensions3;
    use crate::geom3::align3::params::{AlignOrigin, AlignParams3};
    use crate::geom3::tests::RandomGeometry;
    use crate::na::{Rotation, Rotation3, UnitQuaternion};
    use crate::{Iso3, Point3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    const ANGLE_EPSILON: f64 = 1e-8;
    const TRANS_EPSILON: f64 = 1e-8;

    // ============================================================================================
    // These first tests are to check some simple and obvious properties of AlignParams3
    // ============================================================================================

    #[test]
    fn rotation_around_origin() {
        // Checks that the transformation occurs around the origin when the rotation center is set
        // to the origin.
        let params = AlignParams3::new_at_origin(None).with_rz(PI / 2.0);
        let test_point = Point3::new(1.0, 2.0, 3.0);
        let expected = Point3::new(-2.0, 1.0, 3.0);

        assert_relative_eq!(
            expected,
            params.current_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    #[test]
    fn rotation_around_center() {
        // Checks that the transformation occurs around the rotation center when it is nonzero
        let rc = Point3::new(1.0, 0.0, 0.0);
        let params = AlignParams3::new_at_center(rc, None).with_rz(PI / 2.0);

        let test_point = Point3::new(2.0, 2.0, 3.0);
        let expected = Point3::new(-1.0, 1.0, 3.0);

        assert_relative_eq!(
            expected,
            params.current_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    #[test]
    fn translations_are_in_local_origin_directions() {
        let local = Iso3::try_from_basis_xy(
            &Vector3::new(0.0, -1.0, 0.0),
            &Vector3::new(0.0, 0.0, 1.0),
            None,
        )
        .unwrap();
        let params = AlignParams3::new_at_local(local, None).with_tx(1.0);

        let test_point = Point3::new(1.0, 2.0, 3.0);
        let expected = Point3::new(1.0, 1.0, 3.0);

        assert_relative_eq!(
            expected,
            params.current_values().transform * test_point,
            epsilon = 1e-8
        );
    }

    // ============================================================================================
    // These next tests are to check that the current translation directions match the direction
    // that they actually move entities.
    // ============================================================================================

    #[test]
    fn partials_of_translations_at_zero() {
        let params = AlignParams3::new_at_origin(None);
        let current = params.current_values();

        assert_relative_eq!(current.dtx, Vector3::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(current.dty, Vector3::y_axis(), epsilon = 1e-12);
        assert_relative_eq!(current.dtz, Vector3::z_axis(), epsilon = 1e-12);
    }

    #[test]
    fn partials_of_translations_with_rotations() {
        let params = AlignParams3::new_at_origin(None)
            .with_rx(0.1)
            .with_ry(0.2)
            .with_rz(0.3);
        let test_point = Point3::new(1.0, 2.0, 3.0);

        let exp_x = finite_diff(&params, &test_point, 0);
        let exp_y = finite_diff(&params, &test_point, 1);
        let exp_z = finite_diff(&params, &test_point, 2);
        let c = params.current_values();

        // Because translations are applied after rotations, the direction vectors don't change
        // when rotation parameters are applied
        assert_relative_eq!(c.dtx, Vector3::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(c.dty, Vector3::y_axis(), epsilon = 1e-12);
        assert_relative_eq!(c.dtz, Vector3::z_axis(), epsilon = 1e-12);

        // Sanity check that I've thought through this correctly
        assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-8);
        assert_relative_eq!(exp_y, c.dty, epsilon = 1e-8);
        assert_relative_eq!(exp_z, c.dtz, epsilon = 1e-8);
    }

    #[test]
    fn stress_partials_translations() {
        // This is the main test in this section; it checks that the actual translation directions
        // match the measured directions even with random numbers in every possible parameter.
        let mut rg = RandomGeometry::new();
        for _ in 0..1000 {
            let local = rg.iso3(10.0);
            let working = rg.iso3(10.0);
            let params = AlignParams3::new(AlignOrigin::Local(local), Some(working), None)
                .with_storage(rg.vector(PI));

            let test_point = rg.point3(10.0);
            let exp_x = finite_diff(&params, &test_point, 0);
            let exp_y = finite_diff(&params, &test_point, 1);
            let exp_z = finite_diff(&params, &test_point, 2);

            let c = params.current_values();
            assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-6);
            assert_relative_eq!(exp_y, c.dty, epsilon = 1e-6);
            assert_relative_eq!(exp_z, c.dtz, epsilon = 1e-6);
        }
    }

    // ============================================================================================
    // At this point we have to deal with the partials of the rotation by euler angles. I'm not
    // good enough at linear algebra to derive the direct solution, so I'm hoping to figure it out
    // as we go.
    // ============================================================================================

    #[test]
    fn partials_of_rotation_at_zero() {
        let params = AlignParams3::new_at_origin(None);
        let current = params.current_values();
        let p_local = Point3::new(1.0, 2.0, 3.0);
        let p = current.transform * p_local;

        let exp_rx = finite_diff(&params, &p_local, 3);
        let exp_ry = finite_diff(&params, &p_local, 4);
        let exp_rz = finite_diff(&params, &p_local, 5);

        assert_relative_eq!(exp_rx, current.drx(&p), epsilon = 1e-6);
        assert_relative_eq!(exp_ry, current.dry(&p), epsilon = 1e-6);
        assert_relative_eq!(exp_rz, current.drz(&p), epsilon = 1e-6);
    }

    #[test]
    fn verify_euler_order_is_x_y_z() {
        // This verifies that the order of nalgebra's euler angles is x, then y, then z.
        let rx = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), PI / 4.0).to_rotation_matrix();
        let ry = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), PI / 5.0).to_rotation_matrix();
        let rz = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), PI / 6.0).to_rotation_matrix();

        let expected =
            UnitQuaternion::from_euler_angles(PI / 4.0, PI / 5.0, PI / 6.0).to_rotation_matrix();

        // Remember that matrix multiplication is ordered backwards
        let actual = rz * ry * rx;

        assert_relative_eq!(expected, actual, epsilon = 1e-12);
    }

    #[test]
    fn euler_partial_rx() {
        let f = EulerFixture::new(PI / 4.0, 0.0, 0.0);
        let p_local = Point3::new(1.0, 1.0, 1.0);
        let p_rot = f.rot * p_local;

        let (rdx, rdy, rdz) = euler_partials(f.ry, f.rz);
        assert_relative_eq!(f.exp_z(&p_local), rdz * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_y(&p_local), rdy * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_x(&p_local), rdx * p_rot.coords, epsilon = 1e-6);
    }

    #[test]
    fn euler_partial_ry() {
        let f = EulerFixture::new(0.0, PI / 4.0, 0.0);
        let p_local = Point3::new(1.0, 1.0, 1.0);
        let p_rot = f.rot * p_local;

        let (rdx, rdy, rdz) = euler_partials(f.ry, f.rz);
        assert_relative_eq!(f.exp_z(&p_local), rdz * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_y(&p_local), rdy * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_x(&p_local), rdx * p_rot.coords, epsilon = 1e-6);
    }

    #[test]
    fn euler_partial_rz() {
        let f = EulerFixture::new(0.0, 0.0, PI / 4.0);
        let p_local = Point3::new(1.0, 1.0, 1.0);
        let p_rot = f.rot * p_local;

        let (rdx, rdy, rdz) = euler_partials(f.ry, f.rz);
        assert_relative_eq!(f.exp_z(&p_local), rdz * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_y(&p_local), rdy * p_rot.coords, epsilon = 1e-6);
        assert_relative_eq!(f.exp_x(&p_local), rdx * p_rot.coords, epsilon = 1e-6);
    }

    #[test]
    fn stress_euler_partial_rotated() {
        let mut rg = RandomGeometry::new();
        for _ in 0..1000 {
            let f = EulerFixture::new(rg.f64_sym(PI), rg.f64_sym(PI), rg.f64_sym(PI));
            let p_local = rg.point3(10.0);
            let p_rot = f.rot * p_local;
            let (rdx, rdy, rdz) = euler_partials(f.ry, f.rz);

            assert_relative_eq!(f.exp_z(&p_local), rdz * p_rot.coords, epsilon = 1e-6);
            assert_relative_eq!(f.exp_y(&p_local), rdy * p_rot.coords, epsilon = 1e-6);
            assert_relative_eq!(f.exp_x(&p_local), rdx * p_rot.coords, epsilon = 1e-6);
        }
    }

    // ============================================================================================
    // With the rotation vectors figured out from the euler angles at the origin, we just need to
    // figure out how to make this work in the larger context of the full alignment problem, which
    // may have a local origin _and_ a working offset, and also which needs to receive the point
    // being measured in the target entity's space.
    // ============================================================================================

    #[test]
    fn stress_partials_rotation() {
        let mut rg = RandomGeometry::new();
        for _ in 0..1000 {
            let local = rg.iso3(10.0);
            let working = rg.iso3(10.0);
            let params = AlignParams3::new(AlignOrigin::Local(local), Some(working), None)
                .with_storage(rg.vector(PI));
            let c = params.current_values();

            // The test point in the test entity's space
            let test_point = rg.point3(10.0);
            let exp_x = finite_diff(&params, &test_point, 3);
            let exp_y = finite_diff(&params, &test_point, 4);
            let exp_z = finite_diff(&params, &test_point, 5);

            // The test point in the target entity's space
            let target_point = c.transform * test_point;

            // Because of the possible large distance between the test point and the rotation,
            // the difference between the finite difference estimation and the actual theoretical
            // vectors can be well into the fourth decimal place. Since this is a stress test for
            // plausibility, that's totally fine. The theoretical values will be more accurate
            // than the estimate, but they must be reasonably close.
            assert_relative_eq!(exp_x, c.drx(&target_point), epsilon = 1e-3);
            assert_relative_eq!(exp_y, c.dry(&target_point), epsilon = 1e-3);
            assert_relative_eq!(exp_z, c.drz(&target_point), epsilon = 1e-3);
        }
    }

    // ============================================================================================
    // Test support methods
    // ============================================================================================
    struct EulerFixture {
        rx: f64,
        ry: f64,
        rz: f64,
        rot: Rotation<f64, 3>,
    }

    impl EulerFixture {
        fn new(rx: f64, ry: f64, rz: f64) -> Self {
            let rot = UnitQuaternion::from_euler_angles(rx, ry, rz).to_rotation_matrix();
            Self { rx, ry, rz, rot }
        }

        fn exp_x(&self, point: &Point3) -> Vector3 {
            euler_finite_diff(self.rx, self.ry, self.rz, Ax::X, point)
        }

        fn exp_y(&self, point: &Point3) -> Vector3 {
            euler_finite_diff(self.rx, self.ry, self.rz, Ax::Y, point)
        }

        fn exp_z(&self, point: &Point3) -> Vector3 {
            euler_finite_diff(self.rx, self.ry, self.rz, Ax::Z, point)
        }
    }

    enum Ax {
        X,
        Y,
        Z,
    }

    fn euler_finite_diff(rx: f64, ry: f64, rz: f64, axis: Ax, point: &Point3) -> Vector3 {
        let t0 = match axis {
            Ax::X => unit_quat(rx - ANGLE_EPSILON, ry, rz),
            Ax::Y => unit_quat(rx, ry - ANGLE_EPSILON, rz),
            Ax::Z => unit_quat(rx, ry, rz - ANGLE_EPSILON),
        };
        let t1 = match axis {
            Ax::X => unit_quat(rx + ANGLE_EPSILON, ry, rz),
            Ax::Y => unit_quat(rx, ry + ANGLE_EPSILON, rz),
            Ax::Z => unit_quat(rx, ry, rz + ANGLE_EPSILON),
        };
        let p0 = t0 * point;
        let p1 = t1 * point;
        (p1 - p0) / (2.0 * ANGLE_EPSILON)
    }

    fn unit_quat(rx: f64, ry: f64, rz: f64) -> Rotation3<f64> {
        UnitQuaternion::from_euler_angles(rx, ry, rz).to_rotation_matrix()
    }

    /// This function computes the vector of motion of a point in the original coordinate system
    /// of the test entity for a finitely approximated infinitesimal change in one of the
    /// parameters. The parameters by index are: (0: tx, 1: ty, 2: tz, 3: rx, 4: ry, 5: rz)
    fn finite_diff(params: &AlignParams3, point: &Point3, index: usize) -> Vector3 {
        let mut w0 = params.clone();
        let mut w1 = params.clone();
        let stored = params.get_storage();

        let eps = if index < 3 {
            TRANS_EPSILON
        } else {
            ANGLE_EPSILON
        };

        w0.set_index(index, stored[index] - eps);
        w1.set_index(index, stored[index] + eps);

        let p0 = w0.current_values().transform * point;
        let p1 = w1.current_values().transform * point;

        (p1 - p0) / (2.0 * eps)
    }
}
