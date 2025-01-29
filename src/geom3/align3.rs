mod jacobian;
mod multi_param;
mod points_to_mesh;
mod rotations;

use crate::geom3::{Iso3, Point3, Vector3};
use parry3d_f64::na::{Translation3, UnitQuaternion, Vector6};

type T3Storage = Vector6<f64>;

pub use self::points_to_mesh::points_to_mesh;
pub use self::rotations::RotationMatrices;

#[derive(Clone, Copy, Debug)]
pub enum SampleMode {
    All,
    Random(usize),
    Poisson(f64),
}

pub fn iso3_from_param(p: &T3Storage) -> Iso3 {
    Iso3::from_parts(
        Translation3::new(p.x, p.y, p.z),
        UnitQuaternion::from_euler_angles(p.w, p.a, p.b),
    )
}

pub fn param_from_iso3(t: &Iso3) -> T3Storage {
    let v = t.translation.vector;
    let e = t.rotation.euler_angles();
    T3Storage::new(v.x, v.y, v.z, e.0, e.1, e.2)
}

pub fn distance_weight(d: f64, threshold: f64) -> f64 {
    (threshold - d).ceil().clamp(0.0, 1.0)
}

pub fn normal_weight(n: Vector3, n_ref: Vector3) -> f64 {
    // If the normals are pointing in opposite directions, the dot product will be negative,
    // so we clamp it to 0.0, otherwise we want to return 1
    n.dot(&n_ref).ceil().max(0.0)
}

/// This struct manages the parameters for a transformation which is expressed as rotations around
/// a rotation center point that is not at the origin, but with the cardinal axes pointing in the
/// same directions as the global coordinate system.  This lowers the scalar values of parameters
/// on alignments happening far from the origins by largely decoupling the translation and rotation
/// parameters.
///
/// To work, the RcParams struct must be initialized with the rotation center point and it will
/// manage the storage of the parameters and the conversion to and from the Iso3 transformation.
///
/// However, this is complicated when an initial transformation is provided.
#[derive(Clone)]
pub struct RcParams3 {
    /// The rotation center point in the same coordinate system as the test entity(s)
    pub rc: Point3,

    /// The shift from the rotation center point to the origin
    shift0: Iso3,

    /// The shift from the origin to the initial transformed rotation center point
    shift1: Iso3,

    /// The storage for the 6 parameters
    x: T3Storage,

    /// The currently active transformation computed from the parameters `x`
    transform: Iso3,

    /// The currently active inverse transformation computed from the parameters `x`
    inverse: Iso3,

    /// The currently active rotation matrices computed from the parameters `x`
    rotations: RotationMatrices,

    /// The currently active center of rotation, computed by transforming the rotation center point
    /// `rc` by the current transformation `transform`
    current_rc: Point3,
}

impl RcParams3 {
    pub fn from_initial(initial: &Iso3, rc: &Point3) -> Self {
        let rc_d = initial * rc;
        let rotations = RotationMatrices::from_rotation(&initial.rotation);
        let x = T3Storage::new(0.0, 0.0, 0.0, rotations.r.x, rotations.r.y, rotations.r.z);

        let mut item = Self {
            rc: *rc,
            shift0: Iso3::translation(-rc.x, -rc.y, -rc.z),
            shift1: Iso3::translation(rc_d.x, rc_d.y, rc_d.z),
            x,
            transform: Iso3::identity(),
            inverse: Iso3::identity(),
            rotations,
            current_rc: rc_d,
        };

        item.compute();
        item
    }

    pub fn rotations(&self) -> &RotationMatrices {
        &self.rotations
    }

    pub fn current_rc(&self) -> &Point3 {
        &self.current_rc
    }

    pub fn set(&mut self, x: &T3Storage) {
        self.x = *x;
        self.compute();
    }

    pub fn x(&self) -> &T3Storage {
        &self.x
    }

    fn compute(&mut self) {
        self.rotations = RotationMatrices::from_euler(self.x[3], self.x[4], self.x[5]);

        // 1. A translation from the source rotation center point to the origin
        // 2. The transformation encoded by the parameters
        // 3. A translation from the origin to the destination rotation center point
        let p = Iso3::from_parts(
            Translation3::new(self.x.x, self.x.y, self.x.z),
            self.rotations.q,
        );

        self.transform = self.shift1 * p * self.shift0;
        self.inverse = self.transform.inverse();
        self.current_rc = self.transform * self.rc;
    }

    pub fn transform(&self) -> &Iso3 {
        &self.transform
    }

    pub fn inverse(&self) -> &Iso3 {
        &self.inverse
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Vector3;
    use approx::assert_relative_eq;
    use rand::distributions::Uniform;
    use rand::prelude::*;
    use std::f64::consts::PI;

    fn random_iso3() -> Iso3 {
        let mut rng = thread_rng();
        let v = Vector3::new(
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
        );
        let e = Vector3::new(
            Uniform::from(-PI..PI).sample(&mut rng),
            Uniform::from(-PI..PI).sample(&mut rng),
            Uniform::from(-PI..PI).sample(&mut rng),
        );
        Iso3::from_parts(
            Translation3::from(v),
            UnitQuaternion::from_euler_angles(e.x, e.y, e.z),
        )
    }

    fn random_point3() -> Point3 {
        let mut rng = thread_rng();
        Point3::new(
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
        )
    }

    #[test]
    fn test_iso3_param_round_trips_stress_test() {
        for _ in 0..10000 {
            let t = random_iso3();
            let p = param_from_iso3(&t);
            let t2 = iso3_from_param(&p);

            assert_relative_eq!(t.to_matrix(), t2.to_matrix(), epsilon = 1e-10);
        }
    }

    #[test]
    fn test_iso3_param_round_trips_stress_test_rc() {
        for _ in 0..10000 {
            let t = random_iso3();
            let rc = random_point3();
            let p = RcParams3::from_initial(&t, &rc);
            let t2 = p.transform();

            assert_relative_eq!(t.to_matrix(), t2.to_matrix(), epsilon = 1e-10);
        }
    }
}
