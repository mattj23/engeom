//! This module contains the parameterization of the alignment problem

use crate::geom3::align3::rotations::Euler;
use crate::geom3::align3::{Dof6, T3Storage, iso3_from_param};
use crate::na::{Matrix3, Translation};
use crate::{Iso3, Point3};
use parry3d_f64::na::UnitQuaternion;

const EPSILON: f64 = 1e-8;
const SKS_X: Matrix3<f64> = Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
const SKS_Y: Matrix3<f64> = Matrix3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
const SKS_Z: Matrix3<f64> = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

#[derive(Clone, Debug)]
pub struct AlignValues3 {
    /// The transformation of the test entity(s) from their native, local coordinates to the
    /// current position in the target entity's coordinate system. This transformation is a
    /// composite of the active internal parameters and the working transformation.
    pub transform: Iso3,

    /// The location of the rotation center point in the target entity's coordinate system.
    pub rc: Point3,

    /// the numerical jacobian matrix for X
    pub t_rx: Iso3,

    /// the numerical jacobian matrix for Y
    pub t_ry: Iso3,

    /// the numerical jacobian matrix for Z
    pub t_rz: Iso3,

    /// The current working transformation, which is the transformation of the test entity(s) from
    /// their native, local coordinates to the start position in the target entity's coordinate
    /// system.
    pub working: Iso3,
}

#[derive(Clone)]
pub struct AlignParams3 {
    /// The rotation center point in the same coordinate system as the test entity(s)
    pub rc: Point3,

    /// The degrees of freedom that are active
    pub dof: Dof6,

    /// The current working transformation, which is the transformation of the test entity(s) from
    /// their native, local coordinates to the start position for the alignment process.
    pub working: Iso3,

    /// The storage for the six parameters
    storage: T3Storage,

    center_to_origin: Iso3,

    origin_to_center: Iso3,
}

impl AlignParams3 {
    pub fn new(rc: Point3, working: Iso3, dof: Dof6) -> Self {
        let shift = Iso3::translation(rc.x, rc.y, rc.z);
        AlignParams3 {
            rc,
            dof,
            working,
            storage: T3Storage::default(),
            center_to_origin: shift.inverse(),
            origin_to_center: shift,
        }
    }

    pub fn new_with_values(rc: Point3, working: Iso3, dof: Dof6, storage: T3Storage) -> Self {
        let shift = Iso3::translation(rc.x, rc.y, rc.z);
        AlignParams3 {
            rc,
            dof,
            working,
            storage,
            center_to_origin: shift.inverse(),
            origin_to_center: shift,
        }
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

        AlignValues3 {
            transform,
            rc,
            t_rx,
            t_ry,
            t_rz,
            working: self.working,
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

    // pub fn set_tx(&mut self, value: f64) {
    //     self.storage[0] = value;
    //     self.enforce_constraint();
    // }
    //
    // pub fn set_ty(&mut self, value: f64) {
    //     self.storage[1] = value;
    //     self.enforce_constraint();
    // }
    //
    // pub fn set_tz(&mut self, value: f64) {
    //     self.storage[2] = value;
    //     self.enforce_constraint();
    // }
    //
    // pub fn set_rx(&mut self, value: f64) {
    //     self.storage[3] = value;
    //     self.enforce_constraint();
    // }
    //
    // pub fn set_ry(&mut self, value: f64) {
    //     self.storage[4] = value;
    //     self.enforce_constraint();
    // }
    //
    // pub fn set_rz(&mut self, value: f64) {
    //     self.storage[5] = value;
    //     self.enforce_constraint();
    // }
}

fn numeric_perterb(x: &T3Storage, index: usize) -> Iso3 {
    let mut x = x.clone();
    x[index] += EPSILON;
    iso3_from_param(&x)
}

// fn local_jacobians(rx: f64, ry: f64, rz: f64) -> (Iso3, Iso3, Iso3) {
fn local_jacobians(rx: f64, ry: f64, rz: f64) -> (Matrix3<f64>, Matrix3<f64>, Matrix3<f64>) {
    let x = UnitQuaternion::from_euler_angles(rx, 0.0, 0.0);
    let y = UnitQuaternion::from_euler_angles(0.0, ry, 0.0);
    let z = UnitQuaternion::from_euler_angles(0.0, 0.0, rz);

    let q = x * y * z;

    let m = to_matrix(&q);
    let ck = to_matrix(&z);

    let dx = SKS_X * m;
    let dy = m * ck.transpose() * SKS_Y * ck;
    let dz = m * SKS_Z;

    let q_inv = q.inverse().to_rotation_matrix();

    // (wrap_iso3(dx), wrap_iso3(dy), wrap_iso3(dz))
    // (dx * q_inv, dy * q_inv, dz * q_inv)
    (dx, dy, dz)
}

fn wrap_iso3(m: Matrix3<f64>) -> Iso3 {
    Iso3::from_parts(Translation::identity(), UnitQuaternion::from_matrix(&m))
}

fn to_matrix(q: &UnitQuaternion<f64>) -> Matrix3<f64> {
    let m = q.to_rotation_matrix();
    *m.matrix()
}

#[cfg(test)]
mod tests {
    use crate::geom3::align3::params::{AlignParams3, EPSILON};
    use crate::geom3::align3::{Dof6, T3Storage};
    use crate::{Iso3, Point3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    const NUMERIC_EPSILON: f64 = 1e-8;

    fn finite_diff(params: &AlignParams3, point: &Point3, index: usize) -> Vector3 {
        let mut w0 = params.clone();
        let mut w1 = params.clone();

        w0.set_index(index, params.storage[index] - NUMERIC_EPSILON);
        w1.set_index(index, params.storage[index] + NUMERIC_EPSILON);

        let p0 = w0.current_values().transform * point;
        let p1 = w1.current_values().transform * point;

        (p1 - p0) / (2.0 * NUMERIC_EPSILON)
    }

    fn numeric_jacobian(params: &AlignParams3, point: &Point3) -> (Vector3, Vector3, Vector3) {
        (
            finite_diff(params, point, 3),
            finite_diff(params, point, 4),
            finite_diff(params, point, 5),
        )
    }

    #[test]
    fn numeric_matrices() {
        let params = AlignParams3::new_with_values(
            Point3::new(1.0, 2.0, 3.0),
            Iso3::identity(),
            Dof6::default(),
            T3Storage::new(0.1, 0.2, 0.3, 0.1, 0.2, 0.3),
        );

        let p0 = Point3::new(3.3, 2.2, 1.1);

        let current = params.current_values();

        let moved = current.transform * p0;
        let moved_dx = current.t_rx * moved;
        let moved_dy = current.t_ry * moved;
        let moved_dz = current.t_rz * moved;

        let mut p_x = params.clone();
        p_x.set_index(3, params.storage[3] + EPSILON);
        let expected_x = p_x.current_values().transform * p0;

        let mut p_y = params.clone();
        p_y.set_index(4, params.storage[4] + EPSILON);
        let expected_y = p_y.current_values().transform * p0;

        let mut p_z = params.clone();
        p_z.set_index(5, params.storage[5] + EPSILON);
        let expected_z = p_z.current_values().transform * p0;

        assert_relative_eq!(expected_x, moved_dx, epsilon = 1e-8);
        assert_relative_eq!(expected_y, moved_dy, epsilon = 1e-8);
        assert_relative_eq!(expected_z, moved_dz, epsilon = 1e-8);
    }

    /*
    #[test]
    fn jacobian_simple_at_origin() {
        let params = AlignParams3::new(Point3::origin(), Iso3::identity(), Dof6::default());
        let p = Point3::new(1.0, 0.0, 0.0);

        let (ex, ey, ez) = numeric_jacobian(&params, &p);

        let current = params.current_values();

        let tx = current.d_rx * p.coords;
        let ty = current.d_ry * p.coords;
        let tz = current.d_rz * p.coords;

        assert_relative_eq!(tx, ex, epsilon = 1e-12);
        assert_relative_eq!(ty, ey, epsilon = 1e-12);
        assert_relative_eq!(tz, ez, epsilon = 1e-12);
    }

    #[test]
    fn jacobian_simple_at_rc() {
        let params = AlignParams3::new(
            Point3::new(1.0, 0.0, 0.0),
            Iso3::identity(),
            Dof6::default(),
        );
        // The point is in original test entity space
        let p = Point3::new(2.0, 0.0, 0.0);

        let (ex, ey, ez) = numeric_jacobian(&params, &p);
        let current = params.current_values();

        let moved = current.transform * p;
        assert_relative_eq!(
            moved - current.rc,
            Vector3::new(1.0, 0.0, 0.0),
            epsilon = 1e-12
        );

        let m = current.dx_i * moved;
        let tx = current.d_rx * m.coords;
        let ty = current.d_ry * m.coords;
        let tz = current.d_rz * m.coords;

        assert_relative_eq!(tx, ex, epsilon = 1e-12);
        assert_relative_eq!(ty, ey, epsilon = 1e-12);
        assert_relative_eq!(tz, ez, epsilon = 1e-12);
    }

    #[test]
    fn jacobian_simple_shift_at_rc() {
        let params = AlignParams3::new_with_values(
            Point3::new(1.0, 0.0, 0.0),
            Iso3::identity(),
            Dof6::default(),
            T3Storage::new(2.0, 2.0, 2.0, 0.0, 0.0, 0.0),
        );
        let p = Point3::new(2.0, 0.0, 0.0);

        let (ex, ey, ez) = numeric_jacobian(&params, &p);
        let current = params.current_values();

        let moved = current.transform * p;
        assert_relative_eq!(
            moved - current.rc,
            Vector3::new(1.0, 0.0, 0.0),
            epsilon = 1e-12
        );

        let m = current.dx_i * moved;
        let tx = current.d_rx * m.coords;
        let ty = current.d_ry * m.coords;
        let tz = current.d_rz * m.coords;

        assert_relative_eq!(tx, ex, epsilon = 1e-8);
        assert_relative_eq!(ty, ey, epsilon = 1e-8);
        assert_relative_eq!(tz, ez, epsilon = 1e-8);
    }

    #[test]
    fn jacobian_simple_rot_at_rc() {
        let params = AlignParams3::new_with_values(
            Point3::new(1.0, 0.0, 0.0),
            Iso3::identity(),
            Dof6::default(),
            T3Storage::new(0.0, 0.0, 0.0, 0.1, 0.2, 0.3),
        );
        let p = Point3::new(2.0, 0.0, 0.0);

        let (ex, ey, ez) = numeric_jacobian(&params, &p);
        let current = params.current_values();

        let moved = current.transform * p;

        let m = current.dx_i * moved;
        assert_relative_eq!(m, Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
        let check_x = current.d_rx * Vector3::new(1.0, 0.0, 0.0);
        assert_relative_eq!(check_x, Vector3::new(0.0, 0.0, 0.0), epsilon = 1e-12);


        let tx = current.d_rx * m.coords;
        let ty = current.d_ry * m.coords;
        let tz = current.d_rz * m.coords;

        assert_relative_eq!(tx, ex, epsilon = 1e-8);
        assert_relative_eq!(ty, ey, epsilon = 1e-8);
        assert_relative_eq!(tz, ez, epsilon = 1e-8);
    }

     */

    // #[test]
    // fn rot_at_origin() {
    //     let mut params = AlignParams3::new(&Point3::origin(), Dof6::default());
    //     params.set_rz(PI / 2.0);
    //     let p0 = Point3::new(1.0, 0.0, 0.0);
    //
    //     let current = params.current_values(&Iso3::identity());
    //     let test = current.transform * p0;
    //     assert_relative_eq!(test, Point3::new(0.0, 1.0, 0.0), epsilon = 1e-12);
    // }
    //
    // #[test]
    // fn rot_at_rc() {
    //     let mut params = AlignParams3::new(&Point3::new(1.0, 0.0, 0.0), Dof6::default());
    //     params.set_rz(PI / 2.0);
    //     let p0 = Point3::new(2.0, 0.0, 0.0);
    //
    //     let current = params.current_values(&Iso3::identity());
    //     let test = current.transform * p0;
    //     assert_relative_eq!(test, Point3::new(1.0, 1.0, 0.0), epsilon = 1e-12);
    // }
}
