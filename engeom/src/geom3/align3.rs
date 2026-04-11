pub mod jacobian;
mod mesh;
mod mesh_overlap;
mod mesh_to_mesh;
mod multi_mesh;
pub mod multi_param;
mod params;
mod partials;
mod point_stability;
mod points_to_cloud;
mod points_to_mesh;
mod rotations;

use crate::geom3::{Iso3, Point3, Vector3};
use parry3d_f64::na::{Translation3, UnitQuaternion, Vector6};

type T3Storage = Vector6<f64>;

pub use self::mesh::*;
pub use self::mesh_to_mesh::mesh_to_mesh_iterative;
pub use self::multi_mesh::{
    MMOpts, MulMeshAlignPoint, multi_mesh_adjustment, multi_mesh_adjustment_with_points,
};
pub use self::point_stability::{StabilityResult, point_stability, point_stability_reduce};
pub use self::points_to_cloud::points_to_cloud;
pub use self::points_to_mesh::{points_to_mesh, ransac_points_to_mesh};
pub use self::rotations::RotationMatrices;

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

/// This function returns 0.0 if the distance `d` is greater than the `threshold`, otherwise it
/// returns 1.0. It is used for turning off the residuals of sample points that are beyond a
/// distance threshold.
pub fn distance_weight(d: f64, threshold: f64) -> f64 {
    // Branchless version of returning 0.0 if d > threshold, otherwise returning (threshold - d)
    (threshold - d).ceil().clamp(0.0, 1.0)
}

/// This function returns 0.0 if the normals `n` and `n_ref` are pointing in opposite directions,
/// otherwise it returns 1.0. It is used for turning off the residuals of sample points that have
/// normals pointing into different half-spaces.
pub fn normal_weight(n: &Vector3, n_ref: &Vector3) -> f64 {
    // If the normals are pointing in opposite directions, the dot product will be negative,
    // so we clamp it to 0.0, otherwise we want to return 1
    n.dot(n_ref).ceil().max(0.0)
}

/// This struct manages the parameters for a transformation which is expressed as rotations around
/// a rotation center point that is not at the origin, but with the cardinal axes pointing in the
/// same directions as the global coordinate system.  This lowers the scalar values of parameters
/// on alignments happening far from the origins by largely decoupling the translation and rotation
/// parameters.
///
/// To work, the RcParams struct must be initialized with the rotation center point, and it will
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

    pub fn set_index(&mut self, index: usize, value: f64) {
        self.x[index] = value;
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
    use crate::common::linear_space;
    use std::f64::consts::FRAC_PI_2;

    use crate::geom3::tests::RandomGeometry;
    use approx::assert_relative_eq;

    #[test]
    fn iso3_tx() {
        let storage = T3Storage::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(1.0, 0.0, 0.0);
        let p2 = t * p;
        assert_relative_eq!(p2.x, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_ty() {
        let storage = T3Storage::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 1.0, 0.0);
        let p2 = t * p;
        assert_relative_eq!(p2.y, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_tz() {
        let storage = T3Storage::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 0.0, 1.0);
        let p2 = t * p;
        assert_relative_eq!(p2.z, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn iso3_rx() {
        let storage = T3Storage::new(0.0, 0.0, 0.0, FRAC_PI_2, 0.0, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 1.0, 0.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::x_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }

    #[test]
    fn iso3_ry() {
        let storage = T3Storage::new(0.0, 0.0, 0.0, 0.0, FRAC_PI_2, 0.0);
        let t = iso3_from_param(&storage);
        let p = Point3::new(0.0, 0.0, 1.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::y_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }

    #[test]
    fn iso3_rz() {
        let storage = T3Storage::new(0.0, 0.0, 0.0, 0.0, 0.0, FRAC_PI_2);
        let t = iso3_from_param(&storage);
        let p = Point3::new(1.0, 0.0, 0.0);
        let test = t * p;
        let expected = Iso3::rotation(Vector3::z_axis().into_inner() * FRAC_PI_2) * p;
        assert_relative_eq!(test, expected, epsilon = 1e-10);
    }

    #[test]
    fn check_distance_weight() {
        let threshold = 30.0;
        for x in linear_space(0.0, 50.0, 1000).iter() {
            let ex = if *x > threshold { 0.0 } else { 1.0 };
            let w = distance_weight(*x, threshold);
            assert_relative_eq!(w, ex, epsilon = 1e-10);
        }

        let threshold = 0.5;
        for x in linear_space(0.0, 1.0, 1000).iter() {
            let ex = if *x > threshold { 0.0 } else { 1.0 };
            let w = distance_weight(*x, threshold);
            assert_relative_eq!(w, ex, epsilon = 1e-10);
        }
    }

    #[test]
    fn iso3_param_round_trips_stress_test() {
        let mut rg = RandomGeometry::new();
        for _ in 0..10000 {
            let t = rg.iso3(10.0);
            let p = param_from_iso3(&t);
            let t2 = iso3_from_param(&p);

            assert_relative_eq!(t.to_matrix(), t2.to_matrix(), epsilon = 1e-10);
        }
    }

    #[test]
    fn iso3_param_round_trips_stress_test_rc() {
        let mut rg = RandomGeometry::new();
        for _ in 0..10000 {
            let t = rg.iso3(10.0);
            let rc = rg.point3(10.0);
            let p = RcParams3::from_initial(&t, &rc);
            let t2 = p.transform();

            assert_relative_eq!(t.to_matrix(), t2.to_matrix(), epsilon = 1e-10);
        }
    }
}
