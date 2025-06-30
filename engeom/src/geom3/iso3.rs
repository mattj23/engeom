//! This module has additional tools and functions for working with 3D isometries

use crate::{Iso3, Point3, Result, UnitVec3, Vector3};
use parry3d_f64::na::{Matrix3, Translation3};
use parry3d_f64::na::{Matrix4, UnitQuaternion, try_convert};

pub trait IsoExtensions3 {
    fn flip_around_x(&self) -> Iso3;
    fn flip_around_y(&self) -> Iso3;
    fn flip_around_z(&self) -> Iso3;
    fn try_from_array(array: &[f64; 16]) -> Result<Iso3>;

    fn try_from_basis_xy(e0: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3>;
    fn try_from_basis_xz(e0: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3>;
    fn try_from_basis_yz(e1: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3>;
    fn try_from_basis_yx(e1: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3>;
    fn try_from_basis_zx(e2: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3>;
    fn try_from_basis_zy(e2: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

    fn from_rx(angle: f64) -> Iso3;
    fn from_ry(angle: f64) -> Iso3;
    fn from_rz(angle: f64) -> Iso3;
}

impl IsoExtensions3 for Iso3 {
    /// Rotate the isometry in place by 180 degrees around the x-axis. The location of the origin
    /// is not changed, but the y and z directions are reversed.
    fn flip_around_x(&self) -> Self {
        let r = Iso3::rotation(Vector3::x() * std::f64::consts::PI);
        self.translation * r * self.rotation
    }

    /// Rotate the isometry in place by 180 degrees around the y-axis. The location of the origin
    /// is not changed, but the x and z directions are reversed.
    fn flip_around_y(&self) -> Self {
        let r = Iso3::rotation(Vector3::y() * std::f64::consts::PI);
        self.translation * r * self.rotation
    }

    /// Rotate the isometry in place by 180 degrees around the z-axis. The location of the origin
    /// is not changed, but the x and y directions are reversed.
    fn flip_around_z(&self) -> Self {
        let r = Iso3::rotation(Vector3::z() * std::f64::consts::PI);
        self.translation * r * self.rotation
    }

    /// Try to convert a 16 element array into an Iso3. The array is expected to be in row-major
    /// order.
    fn try_from_array(array: &[f64; 16]) -> Result<Self> {
        try_convert(Matrix4::from_row_slice(array)).ok_or("Could not convert to Iso3".into())
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the x-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the y-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e0`: A vector in the world coordinate system that will become the x-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e1`: A vector in the world coordinate system whose component linearly independent of `e0`
    ///   will become the y-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_xy(e0: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e0 = e0.try_normalize(1e-10).ok_or("Could not normalize e0")?;
        let e2 = e0
            .cross(e1)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        let e1 = e2
            .cross(&e0)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e1")?;

        from_bases(e0, e1, e2, origin)
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the x-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the z-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e0`: A vector in the world coordinate system that will become the x-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e2`: A vector in the world coordinate system whose component linearly independent of `e0`
    ///   will become the z-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_xz(e0: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e0 = e0.try_normalize(1e-10).ok_or("Could not normalize e0")?;
        let e1 = e2
            .cross(&e0)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e1")?;
        let e2 = e0
            .cross(&e1)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        from_bases(e0, e1, e2, origin)
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the y-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the z-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e1`: A vector in the world coordinate system that will become the y-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e2`: A vector in the world coordinate system whose component linearly independent of `e1`
    ///   will become the z-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_yz(e1: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e1 = e1.try_normalize(1e-10).ok_or("Could not normalize e1")?;
        let e0 = e1
            .cross(&e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        let e2 = e0
            .cross(&e1)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        from_bases(e0, e1, e2, origin)
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the y-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the x-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e1`: A vector in the world coordinate system that will become the y-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e0`: A vector in the world coordinate system whose component linearly independent of `e1`
    ///   will become the x-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_yx(e1: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e1 = e1.try_normalize(1e-10).ok_or("Could not normalize e1")?;
        let e2 = e0
            .cross(&e1)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        let e0 = e1
            .cross(&e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        from_bases(e0, e1, e2, origin)
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the z-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the x-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e2`: A vector in the world coordinate system that will become the z-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e0`: A vector in the world coordinate system whose component linearly independent of `e2`
    ///   will become the x-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_zx(e2: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e2 = e2.try_normalize(1e-10).ok_or("Could not normalize e2")?;
        let e1 = e2
            .cross(&e0)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        let e0 = e1
            .cross(&e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        from_bases(e0, e1, e2, origin)
    }

    /// Try to create an isometry from two basis vectors and an optional origin. The primary basis
    /// vector will become the z-axis in the isometry, the secondary basis vector will be projected
    /// onto the primary and the remaining component will be the y-axis. The final axis will be
    /// computed by cross product for a right-handed coordinate system.
    ///
    /// The isometry produced by this method will move a point in the basis coordinate system to
    /// where it would be located in the world coordinate system.
    ///
    /// If you want to take features in the world coordinate system and move them into the basis
    /// coordinate system, you need to use the inverse of the isometry.
    ///
    /// # Arguments
    ///
    /// * `e2`: A vector in the world coordinate system that will become the z-axis in the basis
    ///   coordinate system, will be normalized to unit length automatically.
    /// * `e1`: A vector in the world coordinate system whose component linearly independent of `e2`
    ///   will become the y-axis in the basis coordinate system, will be normalized to unit length.
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    fn try_from_basis_zy(e2: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e2 = e2.try_normalize(1e-10).ok_or("Could not normalize e2")?;
        let e0 = e1
            .cross(&e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        let e1 = e2
            .cross(&e0)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        from_bases(e0, e1, e2, origin)
    }

    fn from_rx(angle: f64) -> Iso3 {
        Iso3::rotation(Vector3::x() * angle)
    }

    fn from_ry(angle: f64) -> Iso3 {
        Iso3::rotation(Vector3::y() * angle)
    }

    fn from_rz(angle: f64) -> Iso3 {
        Iso3::rotation(Vector3::z() * angle)
    }
}

fn from_bases(e0: Vector3, e1: Vector3, e2: Vector3, origin: Option<Point3>) -> Result<Iso3> {
    let rot_m = Matrix3::from_columns(&[e0, e1, e2]);
    let r = UnitQuaternion::from_matrix(&rot_m);
    let t = if let Some(o) = origin {
        Translation3::from(o.coords)
    } else {
        Translation3::identity()
    };

    Ok(Iso3::from_parts(t, r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point3;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    struct BasisCheck {
        o: Point3,
        e0: Vector3,
        e1: Vector3,
        e2: Vector3,
        fwd: Iso3,
    }

    impl BasisCheck {
        fn new() -> Self {
            let o = Point3::new(1.0, 2.0, 3.0);
            let angle = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0));
            let fwd = Iso3::from_parts(
                Translation3::from(o.coords),
                UnitQuaternion::new(PI / 4.0 * angle.into_inner()),
            );
            let e0 = fwd * Vector3::x();
            let e1 = fwd * Vector3::y();
            let e2 = fwd * Vector3::z();

            Self { o, e0, e1, e2, fwd }
        }
    }

    #[test]
    fn iso3_try_from_basis_xy() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_xy(&check.e0, &check.e1, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_xz() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_xz(&check.e0, &check.e2, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_yx() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_yx(&check.e1, &check.e0, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_yz() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_yz(&check.e1, &check.e2, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_zx() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_zx(&check.e2, &check.e0, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_zy() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::try_from_basis_zy(&check.e2, &check.e1, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn iso3_try_from_basis_xy_manual() {
        let o = Point3::new(1.0, 2.0, 3.0);
        let e0 = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 0.0));
        let e1 = Vector3::new(0.0, 1.0, 1.0);

        let iso = Iso3::try_from_basis_xy(&e0.into_inner(), &e1, Some(o)).unwrap();

        assert_relative_eq!(iso * Point3::origin(), o, epsilon = 1e-6);
        assert_relative_eq!(
            iso * Point3::new(1.0, 0.0, 0.0),
            o + e0.into_inner() * 1.0,
            epsilon = 1e-6
        );
    }

    #[test]
    fn iso3_try_from_array_simple() {
        let array = [
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let iso = Iso3::try_from_array(&array).unwrap();
        let m = iso.to_matrix();
        let expected = Matrix4::new(
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        );
        assert_relative_eq!(m, expected);
    }

    #[test]
    fn iso3_flip_x() {
        let iso = Iso3::new(Vector3::new(1.0, 2.0, 3.0), Vector3::new(0.0, 0.0, 0.0));
        let flipped = iso.flip_around_x();

        let p = Point3::new(0.0, 0.0, 0.0);
        assert_relative_eq!(flipped * p, Point3::new(1.0, 2.0, 3.0));

        let p1 = Point3::new(1.0, 0.0, 0.0);
        assert_relative_eq!(flipped * p1, Point3::new(2.0, 2.0, 3.0));

        let p2 = Point3::new(0.0, 1.0, 0.0);
        assert_relative_eq!(flipped * p2, Point3::new(1.0, 1.0, 3.0));
    }
}
