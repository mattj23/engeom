//! This module contains the parameterization of the alignment problem

use crate::geom3::align3::{Dof6, RotationMatrices, T3Storage, iso3_from_param};
use crate::na::Matrix3;
use crate::{Iso3, Point3};

#[derive(Clone, Debug)]
pub struct AlignValues3 {
    /// The transformation of the test entity(s) from their native, local coordinates to the
    /// current position in the target entity's coordinate system. This transformation is a
    /// composite of the active internal parameters and the working transformation.
    pub transform: Iso3,

    /// The location of the rotation center point in the target entity's coordinate system.
    pub rc: Point3,

    /// The jacobian of RX
    pub d_rx: Matrix3<f64>,

    /// The jacobian of RY
    pub d_ry: Matrix3<f64>,

    /// The jacobian of RZ
    pub d_rz: Matrix3<f64>,

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

    pub fn current_values(&self) -> AlignValues3 {
        let transform = self.origin_to_center
            * iso3_from_param(&self.storage)
            * self.center_to_origin
            * self.working;

        let rc = transform * self.rc;

        AlignValues3 { transform, rc }
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

    pub fn set_tx(&mut self, value: f64) {
        self.storage[0] = value;
        self.enforce_constraint();
    }

    pub fn set_ty(&mut self, value: f64) {
        self.storage[1] = value;
        self.enforce_constraint();
    }

    pub fn set_tz(&mut self, value: f64) {
        self.storage[2] = value;
        self.enforce_constraint();
    }

    pub fn set_rx(&mut self, value: f64) {
        self.storage[3] = value;
        self.enforce_constraint();
    }

    pub fn set_ry(&mut self, value: f64) {
        self.storage[4] = value;
        self.enforce_constraint();
    }

    pub fn set_rz(&mut self, value: f64) {
        self.storage[5] = value;
        self.enforce_constraint();
    }
}

#[cfg(test)]
mod tests {
    use crate::geom3::align3::Dof6;
    use crate::geom3::align3::params::AlignParams3;
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
