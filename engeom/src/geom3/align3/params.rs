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

    /// Get the stored isometry without the working transformation. If the problem has converged,
    /// this is the final result.
    pub fn final_result(&self) -> Iso3 {
        let local = iso3_from_param(&self.storage);
        self.origin_to_center * local * self.center_to_origin
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
            dof: self.dof,
            rc,
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
}
