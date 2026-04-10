//! This module contains the parameterization of the alignment problem

use crate::common::PCoords;
use crate::geom3::align3::{Dof6, T3Storage, iso3_from_param};
use crate::{Iso3, Point3};

#[derive(Clone, Debug)]
pub struct AlignValues3 {
    pub transform: Iso3,
    pub inverse: Iso3,
    pub rc: Point3,
}

#[derive(Clone)]
pub struct AlignParams3 {
    /// The rotation center point in the same coordinate system as the test entity(s)
    pub rc: Point3,

    pub dof: Dof6,

    /// The storage for the six parameters
    storage: T3Storage,

    center_to_origin: Iso3,

    origin_to_center: Iso3,
}

impl AlignParams3 {
    pub fn new(rc: &impl PCoords<3>, dof: Dof6) -> Self {
        let rc: Point3 = rc.coords().into();
        let shift = Iso3::translation(rc.x, rc.y, rc.z);
        AlignParams3 {
            rc,
            dof,
            storage: T3Storage::default(),
            center_to_origin: shift.inverse(),
            origin_to_center: shift,
        }
    }

    pub fn current_values(&self, working_iso: &Iso3) -> AlignValues3 {
        let transform = self.origin_to_center
            * iso3_from_param(&self.storage)
            * self.center_to_origin
            * working_iso;
        let inverse = transform.inverse();
        let rc = transform * self.rc;

        AlignValues3 {
            transform,
            inverse,
            rc,
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
    use crate::{Iso3, Point3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn rot_at_origin() {
        let mut params = AlignParams3::new(&Point3::origin(), Dof6::default());
        params.set_rz(PI / 2.0);
        let p0 = Point3::new(1.0, 0.0, 0.0);

        let current = params.current_values(&Iso3::identity());
        let test = current.transform * p0;
        assert_relative_eq!(test, Point3::new(0.0, 1.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn rot_at_rc() {
        let mut params = AlignParams3::new(&Point3::new(1.0, 0.0, 0.0), Dof6::default());
        params.set_rz(PI / 2.0);
        let p0 = Point3::new(2.0, 0.0, 0.0);

        let current = params.current_values(&Iso3::identity());
        let test = current.transform * p0;
        assert_relative_eq!(test, Point3::new(1.0, 1.0, 0.0), epsilon = 1e-12);
    }
}
