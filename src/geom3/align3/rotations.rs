use parry3d_f64::na::{Matrix3, UnitQuaternion};
use std::f64::consts::PI;

// These are the skew symmetric matrices
const P_X: Matrix3<f64> = Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
const P_Y: Matrix3<f64> = Matrix3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
const P_Z: Matrix3<f64> = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
const EPSILON: f64 = 1e-8;

#[derive(Clone)]
pub struct Euler<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Euler<T> {
    fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }
}

#[derive(Clone)]
pub struct RotationMatrices {
    pub r: Euler<f64>,
    pub q: UnitQuaternion<f64>,
    pub d: Euler<Matrix3<f64>>,
    pub rd: Euler<Matrix3<f64>>,
}

impl RotationMatrices {
    pub fn from_rotation(q: &UnitQuaternion<f64>) -> Self {
        let m = to_matrix(q);
        let (w, p, r) = to_wpr(&m);
        Self::from_euler(w, p, r)
    }

    pub fn from_euler(rx: f64, ry: f64, rz: f64) -> Self {
        let x = UnitQuaternion::from_euler_angles(rx, 0.0, 0.0);
        let y = UnitQuaternion::from_euler_angles(0.0, ry, 0.0);
        let z = UnitQuaternion::from_euler_angles(0.0, 0.0, rz);

        let q = x * y * z;

        let m = to_matrix(&q);
        let ck = to_matrix(&z);

        let r = Euler::new(rx, ry, rz);
        let d = Euler::new(P_X * m, m * ck.transpose() * P_Y * ck, m * P_Z);

        let qi = q.inverse().to_rotation_matrix();
        let rd = Euler::new(d.x * qi, d.y * qi, d.z * qi);

        Self { r, q, d, rd }
    }
}

fn to_matrix(q: &UnitQuaternion<f64>) -> Matrix3<f64> {
    let m = q.to_rotation_matrix();
    *m.matrix()
}

fn to_wpr(m: &Matrix3<f64>) -> (f64, f64, f64) {
    // https://www.geometrictools.com/Documentation/EulerAngles.pdf
    let sin_y = m[(0, 2)];

    if sin_y > 1.0 - EPSILON {
        let ry = PI / 2.0;
        let rx = m[(1, 0)].atan2(m[(1, 1)]);
        let rz = 0.0;
        (rx, ry, rz)
    } else if sin_y < EPSILON - 1.0 {
        let ry = -PI / 2.0;
        let rx = -(m[(1, 0)].atan2(m[(1, 1)]));
        let rz = 0.0;
        (rx, ry, rz)
    } else {
        let ry = sin_y.asin();
        let rx = (-m[(1, 2)]).atan2(m[(2, 2)]);
        let rz = (-m[(0, 1)]).atan2(m[(0, 0)]);
        (rx, ry, rz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::{Point3, Vector3};
    use approx::assert_relative_eq;
    use rand::{thread_rng, Rng};
    use std::f64::consts::PI;

    const NUMERIC_EPSILON: f64 = 1e-8;
    const TEST_ANGLE: f64 = PI / 4.0;

    /// Creates a rotation matrix from euler angles, but the angles are specified in the order of
    /// `[roll, pitch, yaw]`, while they will be combined in the order of yaw, pitch, roll
    fn euler(x: &[f64; 3]) -> Matrix3<f64> {
        // Internally, the `from_euler_angles` function orders the rotation as roll, pitch, yaw,
        // so we need to decompose them and assemble them in ypr so that they match the order of
        // the jacobian derivations
        let rx = UnitQuaternion::from_euler_angles(x[0], 0.0, 0.0).to_rotation_matrix();
        let ry = UnitQuaternion::from_euler_angles(0.0, x[1], 0.0).to_rotation_matrix();
        let rz = UnitQuaternion::from_euler_angles(0.0, 0.0, x[2]).to_rotation_matrix();
        let q = rx * ry * rz;
        *q.matrix()
    }

    fn finite_diff_motion(p: &Point3, x0: &[f64; 3], index: usize) -> Vector3 {
        // const EPSILON: f64 = 1e-4;
        let e0 = euler(x0);
        let mut x1 = *x0;
        x1[index] += NUMERIC_EPSILON;
        let e1 = euler(&x1);

        let p0 = e0 * *p;
        let p1 = e1 * *p;
        (p1 - p0) / NUMERIC_EPSILON
    }

    #[test]
    fn test_wpr_round_trip_only_rx() {
        let rot = RotationMatrices::from_euler(TEST_ANGLE, 0.0, 0.0);
        let (rx, _ry, _rz) = to_wpr(&to_matrix(&rot.q));
        assert_relative_eq!(rx, TEST_ANGLE, epsilon = 1e-4);
    }

    #[test]
    fn test_wpr_round_trip_only_ry() {
        let rot = RotationMatrices::from_euler(0.0, TEST_ANGLE, 0.0);
        let (_rx, ry, _rz) = to_wpr(&to_matrix(&rot.q));
        assert_relative_eq!(ry, TEST_ANGLE, epsilon = 1e-4);
    }

    #[test]
    fn test_wpr_round_trip_only_rz() {
        let rot = RotationMatrices::from_euler(0.0, 0.0, TEST_ANGLE);
        let (_rx, _ry, rz) = to_wpr(&to_matrix(&rot.q));
        assert_relative_eq!(rz, TEST_ANGLE, epsilon = 1e-4);
    }

    #[test]
    fn test_wpr_round_trip_stress() {
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let rx = rnd.gen_range(-PI..PI);
            let ry = rnd.gen_range(-PI..PI);
            let rz = rnd.gen_range(-PI..PI);

            let m0 = RotationMatrices::from_euler(rx, ry, rz);
            let t0 = to_matrix(&m0.q);

            let (rx1, ry1, rz1) = to_wpr(&t0);
            let m1 = RotationMatrices::from_euler(rx1, ry1, rz1);
            let t1 = to_matrix(&m1.q);

            assert_relative_eq!(t0, t1, epsilon = 1e-4);
        }
    }

    #[test]
    fn test_wpr_rot_mat_round_trip_stress() {
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let m0 = RotationMatrices::from_euler(
                rnd.gen_range(-PI..PI),
                rnd.gen_range(-PI..PI),
                rnd.gen_range(-PI..PI),
            );
            let t0 = to_matrix(&m0.q);

            let m1 = RotationMatrices::from_rotation(&m0.q);
            let t1 = to_matrix(&m1.q);

            assert_relative_eq!(t0, t1, epsilon = 1e-4);
        }
    }

    /// Demonstrates and verifies the method for numerically computing the motion vector of a
    /// point under the influence of a rotation matrix.
    #[test]
    fn test_point_motion_finite_difference() {
        let p = Point3::new(1.0, 1.0, 1.0);
        let theta = 0.0;
        let v = finite_diff_motion(&p, &[theta, 0.0, 0.0], 0);
        assert_relative_eq!(Vector3::new(0.0, -1.0, 1.0), v, epsilon = 1e-4);

        let theta = PI / 2.0;
        let v = finite_diff_motion(&p, &[theta, 0.0, 0.0], 0);
        assert_relative_eq!(Vector3::new(0.0, -1.0, -1.0), v, epsilon = 1e-4);
    }

    #[test]
    fn test_point_motion_analytical() {
        let r = RotationMatrices::from_euler(0.0, 0.0, 0.0);
        let p = Point3::new(1.0, 1.0, 1.0);
        let v = r.d.x * p;
        assert_relative_eq!(Vector3::new(0.0, -1.0, 1.0), v.coords);
    }

    #[test]
    fn test_point_motion_rx() {
        let p = Point3::new(1.0, 1.0, 1.0);
        let e = &[0.7, 0.5, 0.3];

        let v = finite_diff_motion(&p, e, 0);
        let r = RotationMatrices::from_euler(e[0], e[1], e[2]);
        let test = r.d.x * p;

        assert_relative_eq!(v, test.coords, epsilon = 1e-4);
    }

    #[test]
    fn test_point_motion_rz() {
        let p = Point3::new(1.0, 1.0, 1.0);
        let e = &[0.7, 0.5, 0.3];

        let v = finite_diff_motion(&p, e, 2);
        let r = RotationMatrices::from_euler(e[0], e[1], e[2]);
        let test = r.d.z * p;

        assert_relative_eq!(v, test.coords, epsilon = 1e-4);
    }

    #[test]
    fn test_point_motion_rz_neg_90() {
        let p = Point3::new(0.0, 1.0, 0.0);
        let e = &[0.0, 0.0, -PI / 2.0];
        let v = finite_diff_motion(&p, e, 2);
        let r = RotationMatrices::from_euler(e[0], e[1], e[2]);
        let test = r.d.z * p;

        assert_relative_eq!(v, test.coords, epsilon = 1e-4);
    }

    #[test]
    fn test_point_motion_ry() {
        let p = Point3::new(1.0, 1.0, 1.0);
        let e = &[0.7, 0.5, 0.3];

        let v = finite_diff_motion(&p, e, 1);
        let r = RotationMatrices::from_euler(e[0], e[1], e[2]);
        let test = r.d.y * p;

        assert_relative_eq!(v, test.coords, epsilon = 1e-4);
    }
}
