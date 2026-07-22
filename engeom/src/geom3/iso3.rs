//! This module has additional tools and functions for working with 3D isometries

use crate::{Iso3, Line3, Point3, Result, UnitVec3, Vector3};
use consts::PI;
use parry3d_f64::na::{Matrix3, Translation3};
use parry3d_f64::na::{Matrix4, UnitQuaternion, try_convert};
use std::f64::consts;

pub trait IsoExtensions3 {
    /// Tries to create an isometry from a 16-element array, which should contain the values of
    /// the transformation matrix in row-major order (the first four values are the first row,
    /// the next four are the second row, etc.).
    ///
    /// This will return an `Err` if the matrix values don't form a legitimate isometry.
    ///
    /// # Arguments
    ///
    /// * `array`: a 16-element floating point array containing the matrix values
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let array = [
    ///     1.0, 0.0, 0.0, 1.0,
    ///     0.0, 1.0, 0.0, 2.0,
    ///     0.0, 0.0, 1.0, 3.0,
    ///     0.0, 0.0, 0.0, 1.0,
    /// ];
    /// let iso = Iso3::from_array(&array).unwrap();
    /// let moved = iso * Point3::origin();
    /// assert_relative_eq!(moved, Point3::new(1.0, 2.0, 3.0), epsilon = 1e-10);
    /// ```
    fn from_array(array: &[f64; 16]) -> Result<Iso3>;

    /// Create an isometry from an arbitrary orthonormal basis spanning the plane perpendicular to
    /// `z`. There is no guarantee about the orientation of the `x` and `y` axes around the given
    /// `z` beyond their mutual perpendicularity.
    ///
    /// If no origin is specified, the global origin will be used.
    ///
    /// Use this only when the orientation of `x` and `y` really don't matter.
    ///
    /// # Arguments
    ///
    /// * `z`: a unit vector pointing in the direction of the result's Z axis
    /// * `origin`: an optional origin, if not provided the global origin will be used
    ///
    /// returns: Isometry<f64, Unit<Quaternion<f64>>, 3>
    fn from_z_arbitrary_xy(z: &UnitVec3, origin: Option<Point3>) -> Iso3;

    /// Try to create an isometry that rotates around an axis represented by a `Line3` entity by
    /// a given angle. This will only return an `Err` if the length of the line is zero.
    ///
    /// # Arguments
    ///
    /// * `axis`: the axis of rotation for the isometry. When the line direction is pointed at the
    ///   viewer, positive angles will produce rotations viewed as counter-clockwise
    /// * `angle`: the angle of rotation, specified in radians
    ///
    /// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Line3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let axis = Line3::new(Point3::new(1.0, 1.0, 0.0), Vector3::z());
    /// let iso = Iso3::from_rot_axis(&axis, PI / 2.0).unwrap();
    ///
    /// let result = iso * Point3::new(2.0, 1.0, 0.0);
    /// assert_relative_eq!(result, Point3::new(1.0, 2.0, 0.0), epsilon = 1e-6);
    /// ```
    fn from_rot_axis(axis: &Line3, angle: f64) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e0` is tilted 45 degrees between the global x and y axes, while `e1` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e0 = Vector3::new(1.0, 1.0, 0.0);
    /// let e1 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_xy(&e0, &e1, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's x-axis (1, 0, 0) lands on the normalized `e0`
    /// assert_relative_eq!(iso * Point3::new(1.0, 0.0, 0.0), origin + e0.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's y-axis (0, 1, 0) lands on `e1`, unchanged since it was
    /// // already orthogonal to `e0`
    /// assert_relative_eq!(iso * Point3::new(0.0, 1.0, 0.0), origin + e1, epsilon = 1e-10);
    /// ```
    fn from_basis_xy(e0: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e0` is tilted 45 degrees between the global x and y axes, while `e2` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e0 = Vector3::new(1.0, 1.0, 0.0);
    /// let e2 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_xz(&e0, &e2, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's x-axis (1, 0, 0) lands on the normalized `e0`
    /// assert_relative_eq!(iso * Point3::new(1.0, 0.0, 0.0), origin + e0.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's z-axis (0, 0, 1) lands on `e2`, unchanged since it was
    /// // already orthogonal to `e0`
    /// assert_relative_eq!(iso * Point3::new(0.0, 0.0, 1.0), origin + e2, epsilon = 1e-10);
    /// ```
    fn from_basis_xz(e0: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e1` is tilted 45 degrees between the global x and y axes, while `e2` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e1 = Vector3::new(1.0, 1.0, 0.0);
    /// let e2 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_yz(&e1, &e2, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's y-axis (0, 1, 0) lands on the normalized `e1`
    /// assert_relative_eq!(iso * Point3::new(0.0, 1.0, 0.0), origin + e1.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's z-axis (0, 0, 1) lands on `e2`, unchanged since it was
    /// // already orthogonal to `e1`
    /// assert_relative_eq!(iso * Point3::new(0.0, 0.0, 1.0), origin + e2, epsilon = 1e-10);
    /// ```
    fn from_basis_yz(e1: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e1` is tilted 45 degrees between the global x and y axes, while `e0` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e1 = Vector3::new(1.0, 1.0, 0.0);
    /// let e0 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_yx(&e1, &e0, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's y-axis (0, 1, 0) lands on the normalized `e1`
    /// assert_relative_eq!(iso * Point3::new(0.0, 1.0, 0.0), origin + e1.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's x-axis (1, 0, 0) lands on `e0`, unchanged since it was
    /// // already orthogonal to `e1`
    /// assert_relative_eq!(iso * Point3::new(1.0, 0.0, 0.0), origin + e0, epsilon = 1e-10);
    /// ```
    fn from_basis_yx(e1: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e2` is tilted 45 degrees between the global x and y axes, while `e0` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e2 = Vector3::new(1.0, 1.0, 0.0);
    /// let e0 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_zx(&e2, &e0, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's z-axis (0, 0, 1) lands on the normalized `e2`
    /// assert_relative_eq!(iso * Point3::new(0.0, 0.0, 1.0), origin + e2.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's x-axis (1, 0, 0) lands on `e0`, unchanged since it was
    /// // already orthogonal to `e2`
    /// assert_relative_eq!(iso * Point3::new(1.0, 0.0, 0.0), origin + e0, epsilon = 1e-10);
    /// ```
    fn from_basis_zx(e2: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

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
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // `e2` is tilted 45 degrees between the global x and y axes, while `e1` is already
    /// // orthogonal to it, so the two are easy to tell apart in the result.
    /// let e2 = Vector3::new(1.0, 1.0, 0.0);
    /// let e1 = Vector3::new(0.0, 0.0, 1.0);
    /// let origin = Point3::new(1.0, 2.0, 3.0);
    /// let iso = Iso3::from_basis_zy(&e2, &e1, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point3::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's z-axis (0, 0, 1) lands on the normalized `e2`
    /// assert_relative_eq!(iso * Point3::new(0.0, 0.0, 1.0), origin + e2.normalize(), epsilon = 1e-10);
    /// // The basis coordinate system's y-axis (0, 1, 0) lands on `e1`, unchanged since it was
    /// // already orthogonal to `e2`
    /// assert_relative_eq!(iso * Point3::new(0.0, 1.0, 0.0), origin + e1, epsilon = 1e-10);
    /// ```
    fn from_basis_zy(e2: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3>;

    /// Create an isometry that rotates around the x-axis by a given angle. Positive angles rotate
    /// in the clockwise direction.
    ///
    /// # Arguments
    ///
    /// * `angle`: angle of rotation, specified in radians
    ///
    /// returns: Isometry<f64, Unit<Quaternion<f64>>, 3>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_rx(PI / 2.0);
    /// let rotated = iso * Vector3::y();
    /// assert_relative_eq!(rotated, Vector3::z(), epsilon = 1e-10);
    /// ```
    fn from_rx(angle: f64) -> Iso3;

    /// Create an isometry that rotates around the y-axis by a given angle. Positive angles rotate
    /// in the clockwise direction.
    ///
    /// # Arguments
    ///
    /// * `angle`: angle of rotation, specified in radians
    ///
    /// returns: Isometry<f64, Unit<Quaternion<f64>>, 3>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_ry(PI / 2.0);
    /// let rotated = iso * Vector3::z();
    /// assert_relative_eq!(rotated, Vector3::x(), epsilon = 1e-10);
    /// ```
    fn from_ry(angle: f64) -> Iso3;

    /// Create an isometry that rotates around the z-axis by a given angle. Positive angles rotate
    /// in the clockwise direction.
    ///
    /// # Arguments
    ///
    /// * `angle`: angle of rotation, specified in radians
    ///
    /// returns: Isometry<f64, Unit<Quaternion<f64>>, 3>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_rz(PI / 2.0);
    /// let rotated = iso * Vector3::x();
    /// assert_relative_eq!(rotated, Vector3::y(), epsilon = 1e-10);
    /// ```
    fn from_rz(angle: f64) -> Iso3;

    /// Return a copy of the isometry rotated by 180 degrees around its own x-axis. The location of
    /// the origin and the direction of its x-axis with respect to the global coordinate system
    /// are unchanged, but its y and z directions are reversed.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // A tilted, off-origin isometry, so its origin and axis directions don't coincide with
    /// // the global coordinate system and the effect of the flip is easy to see.
    /// let iso = Iso3::from_basis_xy(
    ///     &Vector3::new(1.0, 1.0, 1.0),
    ///     &Vector3::new(-1.0, 1.0, 1.0),
    ///     Some(Point3::new(1.0, 2.0, 3.0)),
    /// )
    /// .unwrap();
    /// let flipped = iso.flipped_around_x();
    ///
    /// assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.x(), flipped.x(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.z(), -flipped.z(), epsilon = 1e-6);
    /// ```
    fn flipped_around_x(&self) -> Iso3;

    /// Return a copy of the isometry rotated by 180 degrees around its own y-axis. The location of
    /// the origin and the direction of its y-axis with respect to the global coordinate system are
    /// unchanged, but its x and z directions are reversed.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // A tilted, off-origin isometry, so its origin and axis directions don't coincide with
    /// // the global coordinate system and the effect of the flip is easy to see.
    /// let iso = Iso3::from_basis_xy(
    ///     &Vector3::new(1.0, 1.0, 1.0),
    ///     &Vector3::new(-1.0, 1.0, 1.0),
    ///     Some(Point3::new(1.0, 2.0, 3.0)),
    /// )
    /// .unwrap();
    /// let flipped = iso.flipped_around_y();
    ///
    /// assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.y(), flipped.y(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.z(), -flipped.z(), epsilon = 1e-6);
    /// ```
    fn flipped_around_y(&self) -> Iso3;

    /// Return a copy of the isometry rotated by 180 degrees around its own z-axis. The location of
    /// the origin and the direction of its z-axis with respect to the global coordinate system
    /// are unchanged, but its x and y directions are reversed.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// // A tilted, off-origin isometry, so its origin and axis directions don't coincide with
    /// // the global coordinate system and the effect of the flip is easy to see.
    /// let iso = Iso3::from_basis_xy(
    ///     &Vector3::new(1.0, 1.0, 1.0),
    ///     &Vector3::new(-1.0, 1.0, 1.0),
    ///     Some(Point3::new(1.0, 2.0, 3.0)),
    /// )
    /// .unwrap();
    /// let flipped = iso.flipped_around_z();
    ///
    /// assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-6);
    /// assert_relative_eq!(iso.z(), flipped.z(), epsilon = 1e-6);
    /// ```
    fn flipped_around_z(&self) -> Iso3;

    /// Return the location of the isometry's origin in the global coordinate system. This is what
    /// you get when you transform the point (0, 0, 0) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Point3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::translation(1.0, 2.0, 3.0);
    /// assert_relative_eq!(iso.origin(), Point3::new(1.0, 2.0, 3.0), epsilon = 1e-10);
    /// ```
    fn origin(&self) -> Point3;

    /// Return the direction of the isometry's x-axis in the global coordinate system. This is what
    /// you get when you transform the unit vector (1, 0, 0) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_rz(PI / 2.0);
    /// assert_relative_eq!(iso.x(), Vector3::y_axis(), epsilon = 1e-10);
    /// ```
    fn x(&self) -> UnitVec3;

    /// Return the direction of the isometry's y-axis in the global coordinate system. This is what
    /// you get when you transform the unit vector (0, 1, 0) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_rx(PI / 2.0);
    /// assert_relative_eq!(iso.y(), Vector3::z_axis(), epsilon = 1e-10);
    /// ```
    fn y(&self) -> UnitVec3;

    /// Return the direction of the isometry's z-axis in the global coordinate system. This is what
    /// you get when you transform the unit vector (0, 0, 1) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso3, Vector3};
    /// use engeom::geom3::IsoExtensions3;
    ///
    /// let iso = Iso3::from_ry(PI / 2.0);
    /// assert_relative_eq!(iso.z(), Vector3::x_axis(), epsilon = 1e-10);
    /// ```
    fn z(&self) -> UnitVec3;
}

impl IsoExtensions3 for Iso3 {
    fn from_array(array: &[f64; 16]) -> Result<Self> {
        try_convert(Matrix4::from_row_slice(array)).ok_or("Could not convert to Iso3".into())
    }

    fn from_z_arbitrary_xy(z: &UnitVec3, origin: Option<Point3>) -> Iso3 {
        let n = z.into_inner();
        let reference = if n.z.abs() < 0.9 {
            Vector3::z()
        } else {
            Vector3::x()
        };
        let x_axis = reference.cross(&n).normalize();
        let y_axis = n.cross(&x_axis);

        let rot_m = Matrix3::from_columns(&[x_axis, y_axis, n]);
        let r = UnitQuaternion::from_matrix(&rot_m);
        let t = if let Some(o) = origin {
            Translation3::from(o.coords)
        } else {
            Translation3::identity()
        };

        Iso3::from_parts(t, r)
    }

    fn from_rot_axis(axis: &Line3, angle: f64) -> Result<Iso3> {
        let dir =
            UnitVec3::try_new(axis.direction, 1e-10).ok_or("Could not normalize axis direction")?;
        let r = UnitQuaternion::from_axis_angle(&dir, angle);
        let t = axis.origin.coords - r * axis.origin.coords;
        Ok(Iso3::from_parts(Translation3::from(t), r))
    }

    fn from_basis_xy(e0: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
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

    fn from_basis_xz(e0: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
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

    fn from_basis_yz(e1: &Vector3, e2: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e1 = e1.try_normalize(1e-10).ok_or("Could not normalize e1")?;
        let e0 = e1
            .cross(e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        let e2 = e0
            .cross(&e1)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        from_bases(e0, e1, e2, origin)
    }

    fn from_basis_yx(e1: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
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

    fn from_basis_zx(e2: &Vector3, e0: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
        let e2 = e2.try_normalize(1e-10).ok_or("Could not normalize e2")?;
        let e1 = e2
            .cross(e0)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e2")?;
        let e0 = e1
            .cross(&e2)
            .try_normalize(1e-10)
            .ok_or("Could not normalize e0")?;
        from_bases(e0, e1, e2, origin)
    }

    fn from_basis_zy(e2: &Vector3, e1: &Vector3, origin: Option<Point3>) -> Result<Iso3> {
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

    fn flipped_around_x(&self) -> Self {
        self * Iso3::rotation(Vector3::x() * PI)
    }

    fn flipped_around_y(&self) -> Self {
        self * Iso3::rotation(Vector3::y() * PI)
    }

    fn flipped_around_z(&self) -> Self {
        self * Iso3::rotation(Vector3::z() * PI)
    }

    fn origin(&self) -> Point3 {
        self * Point3::origin()
    }

    fn x(&self) -> UnitVec3 {
        self * Vector3::x_axis()
    }

    fn y(&self) -> UnitVec3 {
        self * Vector3::y_axis()
    }

    fn z(&self) -> UnitVec3 {
        self * Vector3::z_axis()
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
    use crate::{Point3, UnitVec3};
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
    fn from_basis_xy() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_xy(&check.e0, &check.e1, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_xz() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_xz(&check.e0, &check.e2, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_yx() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_yx(&check.e1, &check.e0, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_yz() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_yz(&check.e1, &check.e2, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_zx() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_zx(&check.e2, &check.e0, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_zy() -> Result<()> {
        let check = BasisCheck::new();
        let iso = Iso3::from_basis_zy(&check.e2, &check.e1, Some(check.o))?;
        assert_relative_eq!(iso, check.fwd, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn from_basis_xy_manual() {
        let o = Point3::new(1.0, 2.0, 3.0);
        let e0 = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 0.0));
        let e1 = Vector3::new(0.0, 1.0, 1.0);

        let iso = Iso3::from_basis_xy(&e0.into_inner(), &e1, Some(o)).unwrap();

        assert_relative_eq!(iso * Point3::origin(), o, epsilon = 1e-6);
        assert_relative_eq!(
            iso * Point3::new(1.0, 0.0, 0.0),
            o + e0.into_inner() * 1.0,
            epsilon = 1e-6
        );
    }

    #[test]
    fn from_array_simple() {
        let array = [
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let iso = Iso3::from_array(&array).unwrap();
        let m = iso.to_matrix();
        let expected = Matrix4::new(
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        );
        assert_relative_eq!(m, expected);
    }

    #[test]
    fn flip_x_semantics() {
        let iso = flip_cs();
        let flipped = iso.flipped_around_x();

        assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
        assert_relative_eq!(iso.x(), flipped.x(), epsilon = 1e-6);
        assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-6);
        assert_relative_eq!(iso.z(), -flipped.z(), epsilon = 1e-6);
    }

    #[test]
    fn flip_y_semantics() {
        let iso = flip_cs();
        let flipped = iso.flipped_around_y();

        assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
        assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-6);
        assert_relative_eq!(iso.y(), flipped.y(), epsilon = 1e-6);
        assert_relative_eq!(iso.z(), -flipped.z(), epsilon = 1e-6);
    }

    #[test]
    fn flip_z_semantics() {
        let iso = flip_cs();
        let flipped = iso.flipped_around_z();

        assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
        assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-6);
        assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-6);
        assert_relative_eq!(iso.z(), flipped.z(), epsilon = 1e-6);
    }

    fn flip_cs() -> Iso3 {
        Iso3::from_basis_xy(
            &Vector3::new(1.0, 1.0, 1.0),
            &Vector3::new(-1.0, 1.0, 1.0),
            Some(Point3::new(1.0, 2.0, 3.0)),
        )
        .unwrap()
    }

    #[test]
    fn rotation_axis() {
        let axis = Line3::new(Point3::new(1.0, 1.0, 0.0), Vector3::z());
        let test = Point3::new(2.0, 1.0, 0.0);
        let expected = Point3::new(1.0, 2.0, 0.0);

        let iso = Iso3::from_rot_axis(&axis, PI / 2.0).unwrap();
        let result = iso * test;
        assert_relative_eq!(result, expected, epsilon = 1e-6);
    }
}
