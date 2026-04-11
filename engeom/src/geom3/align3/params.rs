//! This module contains the parameterization of the alignment problem

use crate::common::PCoords;
use crate::geom3::align3::{Dof6, T3Storage, iso3_from_param};
use crate::na::{Matrix3, Translation};
use crate::{Iso3, Point3, Vector3};
use parry3d_f64::na::UnitQuaternion;

const EPSILON: f64 = 1e-8;

#[derive(Clone, Debug)]
pub struct AlignValues3 {
    /// The transformation of the test entity(s) from their native, local coordinates to the
    /// current position in the target entity's coordinate system. This transformation is a
    /// composite of the active internal parameters and the working transformation.
    pub transform: Iso3,

    /// The degrees of freedom that are active
    pub dof: Dof6,

    /// The direction vector for the partial derivative of tx
    pub dtx: Vector3,

    /// The direction vector for the partial derivative of ty
    pub dty: Vector3,

    /// The direction vector for the partial derivative of tz
    pub dtz: Vector3,

    /// The location of the rotation center point in the target entity's coordinate system.
    rc: Point3,

    /// Perturbation matrix for the numerical Jacobian of Rx. Transform an already moved point by
    /// this matrix to get the point shifted by dRx.
    t_rx: Iso3,

    /// Perturbation matrix for the numerical Jacobian of Ry. Transform an already moved point by
    /// this matrix to get the point shifted by dRy.
    t_ry: Iso3,

    /// Perturbation matrix for the numerical Jacobian of Rz. Transform an already moved point by
    /// this matrix to get the point shifted by dRz.
    t_rz: Iso3,
}

impl AlignValues3 {
    pub fn drx(&self, point: &impl PCoords<3>) -> Vector3 {
        let m = self.t_rx * Point3::from(point.coords());
        (m.coords() - point.coords()) / EPSILON
    }

    pub fn dry(&self, point: &impl PCoords<3>) -> Vector3 {
        let m = self.t_ry * Point3::from(point.coords());
        (m.coords() - point.coords()) / EPSILON
    }

    pub fn drz(&self, point: &impl PCoords<3>) -> Vector3 {
        let m = self.t_rz * Point3::from(point.coords());
        (m.coords() - point.coords()) / EPSILON
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
    Local(Iso3)
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
/// transformation. The transformation consists of three translations and three rotations. By
/// convention rotation is applied first around the origin, and then translation is applied along
/// the cardinal axes.
///
/// To allow for control over the center of rotation and the directions of translation, the
/// alignment can optionally define a local origin for the test entity(s). The local origin is
/// provided as a transformation from the origin of the coordinate system that the test entity(s)
/// geometry is defined in.
///
/// Lastly, there is an optional working transformation that can be applied to the test entity(s)
/// geometry before the alignment is performed. The working transformation is equivalent to
/// transforming all the test entity(s) geometry before performing the alignment, but it allows
/// the original geometry to remain static. The cost of applying the working transformation is
/// negligible.
#[derive(Clone)]
pub struct AlignParams3 {
    /// The local origin
    pub origin: Iso3,

    /// The degrees of freedom that are active
    pub dof: Dof6,

    /// The current working transformation, which is the transformation of the test entity(s) from
    /// their native, local coordinates to the start position for the alignment process.
    pub working: Iso3,

    /// The storage for the six parameters
    storage: T3Storage,
}

impl AlignParams3 {

    pub fn new(origin: AlignOrigin, working: Option<Iso3>, dof: Option<Dof6>) -> Self {
        let origin = match origin {
            AlignOrigin::Origin => Iso3::identity(),
            AlignOrigin::Center(p) => Iso3::translation(p.x, p.y, p.z),
            AlignOrigin::Local(t) => t,
        };

        let working = working.unwrap_or_else(|| Iso3::identity());
        let dof = dof.unwrap_or_else(|| Dof6::all());

        Self {
            origin,
            dof,
            working,
            storage: T3Storage::default(),
        }
    }

    pub fn new_local(origin: AlignOrigin, dof: Option<Dof6>) -> Self {
        Self::new(origin, None, dof)
    }

    /// Get the stored isometry without the working transformation. If the problem has converged,
    /// this is the final result.
    pub fn final_result(&self) -> Iso3 {
        let local = iso3_from_param(&self.storage);
        // self.origin_to_center * local * self.center_to_origin
        todo!()
    }

    pub fn current_values(&self) -> AlignValues3 {
        let local = iso3_from_param(&self.storage);
        let local_drx = numeric_perterb(&self.storage, 3);
        let local_dry = numeric_perterb(&self.storage, 4);
        let local_drz = numeric_perterb(&self.storage, 5);

        let transform = self.origin_to_center * local * self.center_to_origin * self.working;
        let t_drx = self.origin_to_center * local_drx * self.center_to_origin * self.working;
        let t_dry = self.origin_to_center * local_dry * self.center_to_origin * self.working;
        let t_drz = self.origin_to_center * local_drz * self.center_to_origin * self.working;
        let t_rx = transform.inverse() * t_drx;
        let t_ry = transform.inverse() * t_dry;
        let t_rz = transform.inverse() * t_drz;

        let rc = transform * self.rc;

        let dtx = Vector3::x_axis().into_inner();
        let dty = Vector3::y_axis().into_inner();
        let dtz = Vector3::z_axis().into_inner();

        AlignValues3 {
            transform,
            dof: self.dof,
            rc,
            dtx,
            dty,
            dtz,
            t_rx,
            t_ry,
            t_rz,
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
}

fn numeric_perterb(x: &T3Storage, index: usize) -> Iso3 {
    let mut x = x.clone();
    x[index] += EPSILON;
    iso3_from_param(&x)
}

#[cfg(test)]
mod tests {
    use crate::geom3::align3::params::{AlignParams3, EPSILON};
    use crate::geom3::align3::{Dof6, T3Storage};
    use crate::{Iso3, Point3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    use crate::na::{Translation3, UnitQuaternion};

    

    // #[test]
    // fn working_transform_is_applied() {
    //     // The idea of the working transformation is that it has the same effect as moving the
    //     // entire test entity(s) to a new position before
    //
    //     let working = Iso3::from_parts(
    //         Translation3::new(2.0, 3.0, 4.0),
    //         UnitQuaternion::from_euler_angles(PI / 3.0, PI / 4.0, PI / 5.0),
    //     );
    //
    // }
    //
    // #[test]
    // fn rotation_around_origin() {
    //     // Checks that the transformation occurs around the origin when the rotation center is set
    //     // to the origin.
    //     let params = AlignParams3::new_local(Point3::origin(), Dof6::all()).with_rz(PI / 2.0);
    //
    //     let test_point = Point3::new(1.0, 2.0, 3.0);
    //     let expected = Point3::new(-2.0, 1.0, 3.0);
    //
    //     assert_relative_eq!(
    //         expected,
    //         params.current_values().transform * test_point,
    //         epsilon = 1e-8
    //     );
    // }
    //
    // #[test]
    // fn rotation_around_center() {
    //     // Checks that the transformation occurs around the rotation center when it is nonzero
    //     let rc = Point3::new(1.0, 0.0, 0.0);
    //     let params = AlignParams3::new_local(rc, Dof6::all()).with_rz(PI / 2.0);
    //
    //     let test_point = Point3::new(2.0, 2.0, 3.0);
    //     let expected = Point3::new(-1.0, 1.0, 3.0);
    //
    //     assert_relative_eq!(
    //         expected,
    //         params.current_values().transform * test_point,
    //         epsilon = 1e-8
    //     );
    // }
    //
    // #[test]
    // fn numeric_matrices() {
    //     let params = AlignParams3::new_with_values(
    //         Point3::new(1.0, 2.0, 3.0),
    //         Iso3::identity(),
    //         Dof6::default(),
    //         T3Storage::new(0.1, 0.2, 0.3, 0.1, 0.2, 0.3),
    //     );
    //
    //     let p0 = Point3::new(3.3, 2.2, 1.1);
    //
    //     let current = params.current_values();
    //
    //     let moved = current.transform * p0;
    //     let moved_dx = current.t_rx * moved;
    //     let moved_dy = current.t_ry * moved;
    //     let moved_dz = current.t_rz * moved;
    //
    //     let mut p_x = params.clone();
    //     p_x.set_index(3, params.storage[3] + EPSILON);
    //     let expected_x = p_x.current_values().transform * p0;
    //
    //     let mut p_y = params.clone();
    //     p_y.set_index(4, params.storage[4] + EPSILON);
    //     let expected_y = p_y.current_values().transform * p0;
    //
    //     let mut p_z = params.clone();
    //     p_z.set_index(5, params.storage[5] + EPSILON);
    //     let expected_z = p_z.current_values().transform * p0;
    //
    //     assert_relative_eq!(expected_x, moved_dx, epsilon = 1e-8);
    //     assert_relative_eq!(expected_y, moved_dy, epsilon = 1e-8);
    //     assert_relative_eq!(expected_z, moved_dz, epsilon = 1e-8);
    // }
}
