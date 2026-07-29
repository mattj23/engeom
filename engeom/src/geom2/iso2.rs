//! This module has additional tools and functions for working with 2D isometries

use crate::{Iso2, Point2, Result, UnitVec2, Vector2};
use parry2d_f64::na::{Matrix3, Translation2, UnitComplex, try_convert};
use std::f64::consts::PI;

pub trait IsoExtensions2 {
    /// Create an isometry that translates by the given x and y distances, with no rotation.
    ///
    /// This is a pass-through to nalgebra's `Iso2::translation`, provided so that the
    /// constructor family on `Iso2` follows this project's `from_<description>` naming
    /// convention. It also avoids the ambiguity of the underlying name, which nalgebra uses for
    /// both this constructor and the isometry's public `translation` field.
    ///
    /// # Arguments
    ///
    /// * `x`: the distance to translate along the x-axis
    /// * `y`: the distance to translate along the y-axis
    ///
    /// returns: Isometry<f64, Unit<Complex<f64>>, 2>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let iso = Iso2::from_translation(1.0, 2.0);
    /// assert_relative_eq!(iso * Point2::origin(), Point2::new(1.0, 2.0), epsilon = 1e-10);
    /// ```
    fn from_translation(x: f64, y: f64) -> Iso2;

    /// Create an isometry that rotates around the global origin by a given angle, with no
    /// translation. To rotate around a different point, use `from_rotation_about`.
    ///
    /// This is a pass-through to nalgebra's `Iso2::rotation`, provided so that the constructor
    /// family on `Iso2` follows this project's `from_<description>` naming convention.
    ///
    /// # Arguments
    ///
    /// * `angle`: the angle of rotation, specified in radians. Positive angles rotate
    ///   counter-clockwise.
    ///
    /// returns: Isometry<f64, Unit<Complex<f64>>, 2>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let iso = Iso2::from_rotation(FRAC_PI_2);
    /// assert_relative_eq!(iso * Vector2::x(), Vector2::y(), epsilon = 1e-10);
    /// ```
    fn from_rotation(angle: f64) -> Iso2;

    /// Tries to create an isometry from a 9-element array, which should contain the values of
    /// the transformation matrix in row-major order (the first three values are the first row,
    /// the next three are the second row, etc.).
    ///
    /// This will return an `Err` if the matrix values don't form a legitimate isometry.
    ///
    /// # Arguments
    ///
    /// * `array`: a 9-element floating point array containing the matrix values
    ///
    /// returns: Result<Isometry<f64, Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let array = [
    ///     1.0, 0.0, 1.0,
    ///     0.0, 1.0, 2.0,
    ///     0.0, 0.0, 1.0,
    /// ];
    /// let iso = Iso2::from_array(&array).unwrap();
    /// let moved = iso * Point2::origin();
    /// assert_relative_eq!(moved, Point2::new(1.0, 2.0), epsilon = 1e-10);
    /// ```
    fn from_array(array: &[f64; 9]) -> Result<Iso2>;

    /// Create an isometry that rotates around an arbitrary point by a given angle. Unlike
    /// `Iso2::rotation`, which always rotates around the global origin, this allows the center of
    /// rotation to be any point in the plane.
    ///
    /// In 3D, the equivalent operation (`IsoExtensions3::from_rot_axis`) takes a `Line3` because
    /// the axis of rotation can point in any direction. In 2D the axis of rotation is always
    /// perpendicular to the plane, so it is fully specified by the single point it passes through.
    ///
    /// # Arguments
    ///
    /// * `point`: the point in the plane that the isometry will rotate around
    /// * `angle`: the angle of rotation, specified in radians. Positive angles rotate
    ///   counter-clockwise.
    ///
    /// returns: Isometry<f64, Unit<Complex<f64>>, 2>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::PI;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let center = Point2::new(1.0, 1.0);
    /// let iso = Iso2::from_rotation_about(&center, PI / 2.0);
    ///
    /// let result = iso * Point2::new(2.0, 1.0);
    /// assert_relative_eq!(result, Point2::new(1.0, 2.0), epsilon = 1e-6);
    /// ```
    fn from_rotation_about(point: &Point2, angle: f64) -> Iso2;

    /// Try to create an isometry from a single basis vector and an optional origin. The vector
    /// will become the x-axis of the isometry (after being normalized to unit length), and the
    /// y-axis will be the x-axis rotated 90 degrees counter-clockwise.
    ///
    /// Unlike the 3D `from_basis_xy`/`from_basis_xz`/etc. family, a single vector is enough to
    /// fully determine a 2D isometry's rotation, since there is only one rotational degree of
    /// freedom and the second axis is always a fixed 90 degree turn away from the first.
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
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let origin = Point2::new(1.0, 2.0);
    /// let e0 = Vector2::new(1.0, 1.0);
    /// let iso = Iso2::from_basis_x(&e0, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point2::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's x-axis (1, 0) lands on the normalized `e0`
    /// assert_relative_eq!(iso * Point2::new(1.0, 0.0), origin + e0.normalize(), epsilon = 1e-10);
    /// ```
    fn from_basis_x(e0: &Vector2, origin: Option<Point2>) -> Result<Iso2>;

    /// Try to create an isometry from a single basis vector and an optional origin. The vector
    /// will become the y-axis of the isometry (after being normalized to unit length), and the
    /// x-axis will be the y-axis rotated 90 degrees clockwise.
    ///
    /// Unlike the 3D `from_basis_xy`/`from_basis_xz`/etc. family, a single vector is enough to
    /// fully determine a 2D isometry's rotation, since there is only one rotational degree of
    /// freedom and the second axis is always a fixed 90-degree turn away from the first.
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
    /// * `origin`: An optional point in the world coordinate system that will be the origin of the
    ///   basis coordinate system. If not provided, the origin of the basis coordinate system will
    ///   be coincident with the origin of the world coordinate system.
    ///
    /// returns: Result<Isometry<f64, Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let origin = Point2::new(1.0, 2.0);
    /// let e1 = Vector2::new(-1.0, 1.0);
    /// let iso = Iso2::from_basis_y(&e1, Some(origin)).unwrap();
    ///
    /// assert_relative_eq!(iso * Point2::origin(), origin, epsilon = 1e-10);
    /// // The basis coordinate system's y-axis (0, 1) lands on the normalized `e1`
    /// assert_relative_eq!(iso * Point2::new(0.0, 1.0), origin + e1.normalize(), epsilon = 1e-10);
    /// ```
    fn from_basis_y(e1: &Vector2, origin: Option<Point2>) -> Result<Iso2>;

    /// Return a copy of the isometry rotated 180 degrees in-plane around its own origin.
    ///
    /// This is a **rotation, not a mirror/reflection**. The location of the origin is unchanged,
    /// but both the x-axis and y-axis directions are reversed together, since in 2D the only axis
    /// a proper (determinant +1) 180-degree rotation can turn around is the implicit axis
    /// perpendicular to the plane. There is no way to reverse only one of the two in-plane axes
    /// and remain a rigid-body transformation: doing so would be a reflection (determinant -1),
    /// which `Iso2` cannot represent at all, since its rotation component is stored as a unit
    /// complex number.
    ///
    /// This is the 2D analog of `IsoExtensions3::flipped_around_z`; there is no 2D equivalent of
    /// `flipped_around_x`/`flipped_around_y`, because those reverse only one in-plane axis and are
    /// therefore reflections rather than rotations.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// // A tilted, off-origin isometry, so its origin and axis directions don't coincide with
    /// // the global coordinate system and the effect of the flip is easy to see.
    /// let iso = Iso2::from_basis_x(&Vector2::new(1.0, 1.0), Some(Point2::new(1.0, 2.0))).unwrap();
    /// let flipped = iso.flipped();
    ///
    /// assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-10);
    /// assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-10);
    /// assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-10);
    /// ```
    fn flipped(&self) -> Iso2;

    /// Return the location of the isometry's origin in the global coordinate system. This is what
    /// you get when you transform the point (0, 0) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Point2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let iso = Iso2::from_translation(1.0, 2.0);
    /// assert_relative_eq!(iso.origin(), Point2::new(1.0, 2.0), epsilon = 1e-10);
    /// ```
    fn origin(&self) -> Point2;

    /// Return the direction of the isometry's x-axis in the global coordinate system. This is what
    /// you get when you transform the unit vector (1, 0) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let iso = Iso2::from_rotation(FRAC_PI_2);
    /// assert_relative_eq!(iso.x(), Vector2::y_axis(), epsilon = 1e-10);
    /// ```
    fn x(&self) -> UnitVec2;

    /// Return the direction of the isometry's y-axis in the global coordinate system. This is what
    /// you get when you transform the unit vector (0, 1) by this isometry.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    /// use approx::assert_relative_eq;
    /// use engeom::{Iso2, Vector2};
    /// use engeom::geom2::IsoExtensions2;
    ///
    /// let iso = Iso2::from_rotation(FRAC_PI_2);
    /// assert_relative_eq!(iso.y(), -Vector2::x_axis(), epsilon = 1e-10);
    /// ```
    fn y(&self) -> UnitVec2;
}

impl IsoExtensions2 for Iso2 {
    fn from_translation(x: f64, y: f64) -> Iso2 {
        Iso2::translation(x, y)
    }

    fn from_rotation(angle: f64) -> Iso2 {
        Iso2::rotation(angle)
    }

    fn from_array(array: &[f64; 9]) -> Result<Self> {
        try_convert(Matrix3::from_row_slice(array)).ok_or("Could not convert to Iso2".into())
    }

    fn from_rotation_about(point: &Point2, angle: f64) -> Iso2 {
        let r = UnitComplex::new(angle);
        let t = point.coords - r * point.coords;
        Iso2::from_parts(Translation2::from(t), r)
    }

    fn from_basis_x(e0: &Vector2, origin: Option<Point2>) -> Result<Iso2> {
        let e0 = e0.try_normalize(1e-10).ok_or("Could not normalize e0")?;
        let r = UnitComplex::from_cos_sin_unchecked(e0.x, e0.y);
        Ok(from_parts(r, origin))
    }

    fn from_basis_y(e1: &Vector2, origin: Option<Point2>) -> Result<Iso2> {
        let e1 = e1.try_normalize(1e-10).ok_or("Could not normalize e1")?;
        let r = UnitComplex::from_cos_sin_unchecked(e1.y, -e1.x);
        Ok(from_parts(r, origin))
    }

    fn flipped(&self) -> Self {
        self * Iso2::rotation(PI)
    }

    fn origin(&self) -> Point2 {
        self * Point2::origin()
    }

    fn x(&self) -> UnitVec2 {
        self * Vector2::x_axis()
    }

    fn y(&self) -> UnitVec2 {
        self * Vector2::y_axis()
    }
}

fn from_parts(r: UnitComplex<f64>, origin: Option<Point2>) -> Iso2 {
    let t = if let Some(o) = origin {
        Translation2::from(o.coords)
    } else {
        Translation2::identity()
    };

    Iso2::from_parts(t, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn from_translation_matches_nalgebra() {
        let iso = Iso2::from_translation(1.0, 2.0);

        assert_relative_eq!(
            iso * Point2::origin(),
            Point2::new(1.0, 2.0),
            epsilon = 1e-10
        );
        assert_relative_eq!(iso, Iso2::translation(1.0, 2.0), epsilon = 1e-10);
    }

    #[test]
    fn from_rotation_matches_nalgebra() {
        let iso = Iso2::from_rotation(PI / 2.0);

        assert_relative_eq!(iso * Vector2::x(), Vector2::y(), epsilon = 1e-10);
        assert_relative_eq!(iso, Iso2::rotation(PI / 2.0), epsilon = 1e-10);
    }

    #[test]
    fn from_array_simple() {
        let array = [1.0, 0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 1.0];
        let iso = Iso2::from_array(&array).unwrap();
        assert_relative_eq!(
            iso * Point2::origin(),
            Point2::new(1.0, 2.0),
            epsilon = 1e-10
        );
    }

    #[test]
    fn from_rotation_about_point() {
        let center = Point2::new(1.0, 1.0);
        let test = Point2::new(2.0, 1.0);
        let expected = Point2::new(1.0, 2.0);

        let iso = Iso2::from_rotation_about(&center, PI / 2.0);
        let result = iso * test;
        assert_relative_eq!(result, expected, epsilon = 1e-6);
    }

    #[test]
    fn from_basis_x_matches_angle() {
        let origin = Point2::new(1.0, 2.0);
        let e0 = Vector2::new(1.0, 1.0);
        let iso = Iso2::from_basis_x(&e0, Some(origin)).unwrap();

        assert_relative_eq!(iso * Point2::origin(), origin, epsilon = 1e-10);
        assert_relative_eq!(
            iso * Point2::new(1.0, 0.0),
            origin + e0.normalize(),
            epsilon = 1e-10
        );
        assert_relative_eq!(
            iso * Point2::new(0.0, 1.0),
            origin + Vector2::new(-e0.normalize().y, e0.normalize().x),
            epsilon = 1e-10
        );
    }

    #[test]
    fn from_basis_y_matches_angle() {
        let origin = Point2::new(1.0, 2.0);
        let e1 = Vector2::new(-1.0, 1.0);
        let iso = Iso2::from_basis_y(&e1, Some(origin)).unwrap();

        assert_relative_eq!(iso * Point2::origin(), origin, epsilon = 1e-10);
        assert_relative_eq!(
            iso * Point2::new(0.0, 1.0),
            origin + e1.normalize(),
            epsilon = 1e-10
        );
        assert_relative_eq!(
            iso * Point2::new(1.0, 0.0),
            origin + Vector2::new(e1.normalize().y, -e1.normalize().x),
            epsilon = 1e-10
        );
    }

    #[test]
    fn flip_semantics() {
        let iso = Iso2::from_basis_x(&Vector2::new(1.0, 1.0), Some(Point2::new(1.0, 2.0))).unwrap();
        let flipped = iso.flipped();

        assert_relative_eq!(iso.origin(), flipped.origin(), epsilon = 1e-6);
        assert_relative_eq!(iso.x(), -flipped.x(), epsilon = 1e-6);
        assert_relative_eq!(iso.y(), -flipped.y(), epsilon = 1e-6);
    }
}
