use crate::common::PCoords;
use crate::geom3::UnitVec3;
use crate::geom3::line3::Line3;
use crate::{Iso3, Point3, SurfacePoint3, SvdBasis3, Vector3};
use std::ops;

#[derive(Debug, Clone)]
pub struct Plane3 {
    pub normal: UnitVec3,
    pub d: f64,
}

impl Plane3 {
    /// Creates a plane with normal along the x-axis and offset 0.0
    pub fn yz() -> Self {
        Self::new(Vector3::x_axis(), 0.0)
    }

    /// Creates a plane with normal along the y-axis and offset 0.0
    pub fn xz() -> Self {
        Self::new(Vector3::y_axis(), 0.0)
    }

    /// Creates a plane with normal along the z-axis and offset 0.0
    pub fn xy() -> Self {
        Self::new(Vector3::z_axis(), 0.0)
    }

    pub fn new(normal: UnitVec3, d: f64) -> Self {
        Self { normal, d }
    }

    /// Create a new plane which is in the same position as the input plane, but with the normal
    /// direction inverted.
    pub fn inverted_normal(&self) -> Self {
        Self::new(-self.normal, -self.d)
    }

    /// Measure and return the signed distance from the plane to a point in 3D space. The sign of
    /// the distance indicates whether the point is above or below the plane according to the
    /// plane's normal vector.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn signed_distance_to_point(&self, point: &impl PCoords<3>) -> f64 {
        self.normal.dot(&point.coords()) - self.d
    }

    /// Returns true if the point lies in the positive half-space defined by the plane's normal and
    /// offset from the origin. This will return true if the signed distance to the point is >= 0.0
    ///
    /// # Arguments
    ///
    /// * `point`: the point to check against the plane
    ///
    /// returns: bool
    pub fn point_is_positive(&self, point: &impl PCoords<3>) -> bool {
        self.signed_distance_to_point(point) >= 0.0
    }

    /// Measure and return the distance from the plane to a point in 3D space. The distance is
    /// always positive and indicates the shortest distance from the point to the plane. If you
    /// need to know whether the point is above or below the plane, use `signed_distance_to_point`.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn distance_to_point(&self, point: &Point3) -> f64 {
        self.signed_distance_to_point(point).abs()
    }

    /// Project a point onto the plane, returning a point in 3D space which lies on the plane. This
    /// is also the closest point on the plane to the input point.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: OPoint<f64, Const<3>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn project_point(&self, point: &Point3) -> Point3 {
        point - self.normal.into_inner() * self.signed_distance_to_point(point)
    }

    /// Projects a vector onto the plane by removing the component along the plane's normal.
    pub fn project_vector(&self, v: &Vector3) -> Vector3 {
        v - self.normal.into_inner() * self.normal.dot(v)
    }

    /// Intersects this plane with another plane, returning the line of intersection, or `None` if
    /// the planes are parallel (or coincident).
    ///
    /// The line's direction is `self.normal × other.normal` (not normalized). The origin is the
    /// point on the intersection line closest to the world origin.
    pub fn intersect_plane(&self, other: &Plane3) -> Option<Line3> {
        let direction = self.normal.cross(&other.normal);
        let denom = direction.norm_squared();
        if denom < 1e-20 {
            return None;
        }
        // Point closest to origin on the intersection line via Lagrange multipliers:
        //   p = λ1*n1 + λ2*n2, solving n1·p = d1, n2·p = d2
        let k = self.normal.dot(&other.normal);
        let denom_lm = 1.0 - k * k;
        let l1 = (self.d - k * other.d) / denom_lm;
        let l2 = (other.d - k * self.d) / denom_lm;
        let origin = Point3::from(self.normal.into_inner() * l1 + other.normal.into_inner() * l2);
        Some(Line3::new(origin, direction))
    }

    pub fn intersection_distance(&self, sp: &SurfacePoint3) -> Option<f64> {
        let p0 = Point3::from(self.normal.into_inner() * self.d);

        let denom = self.normal.dot(&sp.normal);
        if denom <= 1e-6 {
            None
        } else {
            Some((p0 - sp.point).dot(&self.normal) / denom)
        }
    }

    /// Create a new plane parallel to this one, displaced along its normal direction by the given
    /// distance. A positive value moves in the normal direction; negative moves opposite.
    ///
    /// # Arguments
    ///
    /// * `shift`: The distance to shift the plane along its normal vector.
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::geom3::{Plane3, Point3, Vector3};
    /// use approx::assert_relative_eq;
    /// let plane = Plane3::new(Vector3::x_axis(), -5.0);
    /// let moved = plane.new_parallel(2.0);
    ///
    /// assert_relative_eq!(moved.signed_distance_to_point(&Point3::origin()), 3.0, epsilon = 1e-6);
    /// ```
    pub fn new_parallel(&self, shift: f64) -> Self {
        Self::new(self.normal, self.d + shift)
    }

    /// Returns a new plane transformed by the given isometry.
    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        let pos = self.normal.into_inner() * self.d;
        let repr = SurfacePoint3::new(pos.into(), self.normal);
        let new_repr = repr.transformed(iso);
        Self::from(&new_repr)
    }

    /// Transforms this plane in place by the given isometry.
    pub fn transform_by(&mut self, iso: &Iso3) {
        *self = self.new_transformed_by(iso);
    }
}

impl From<&SvdBasis3> for Plane3 {
    /// Create a Plane3 from a SvdBasis3 using the third basis vector as the normal and the mean
    /// point to calculate d. If a `SvdBasis3` has been constructed from a set of planar points,
    /// this will create a plane that best fits those points.
    ///
    /// # Arguments
    ///
    /// * `svd`: The SvdBasis3 to create the plane from
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::geom3::{Plane3, SvdBasis3, Point3};
    /// let points = vec![
    ///    Point3::new(5.0, 10.0, 15.0),
    ///    Point3::new(5.0, 11.0, 16.0),
    ///    Point3::new(5.0, 10.0, 16.0),
    ///    Point3::new(5.0, 11.0, 15.0),
    /// ];
    /// let svd = SvdBasis3::from_points(&points, None).unwrap();
    /// let plane = Plane3::from(&svd);
    /// assert_relative_eq!(plane.normal.x, 1.0, epsilon = 1e-6);
    /// assert_relative_eq!(plane.d, 5.0, epsilon = 1e-6);
    /// ```
    fn from(svd: &SvdBasis3) -> Self {
        let normal = UnitVec3::new_normalize(svd.basis[2]);
        let d = normal.dot(&svd.center.coords);
        Self::new(normal, d)
    }
}

// TODO: should this be a Result?
impl From<(&Point3, &Point3, &Point3)> for Plane3 {
    /// Create a Plane3 from three points
    ///
    /// # Arguments
    ///
    /// * `(p1, p2, p3)`:
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn from((p1, p2, p3): (&Point3, &Point3, &Point3)) -> Self {
        let normal = UnitVec3::new_normalize((p2 - p1).cross(&(p3 - p1)));
        Self::from((&normal, p1))
    }
}

impl From<(&UnitVec3, &Point3)> for Plane3 {
    fn from((normal, point): (&UnitVec3, &Point3)) -> Self {
        let d = normal.dot(&point.coords);
        Self::new(*normal, d)
    }
}

impl From<&SurfacePoint3> for Plane3 {
    fn from(surface_point: &SurfacePoint3) -> Self {
        Self::from((&surface_point.normal, &surface_point.point))
    }
}

impl ops::Mul<Plane3> for Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: Plane3) -> Plane3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<&Plane3> for Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: &Plane3) -> Plane3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<Plane3> for &Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: Plane3) -> Plane3 {
        rhs.new_transformed_by(self)
    }
}

impl ops::Mul<&Plane3> for &Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: &Plane3) -> Plane3 {
        rhs.new_transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry;
    use approx::assert_relative_eq;

    #[test]
    fn intersect_xy_xz_gives_x_axis() {
        // xy-plane (z=0) ∩ xz-plane (y=0) should be the X axis
        let line = Plane3::xy().intersect_plane(&Plane3::xz()).unwrap();
        // direction must be parallel to X
        let dir = line.direction().normalize();
        assert_relative_eq!(dir.x.abs(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(dir.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(dir.z, 0.0, epsilon = 1e-12);
        // origin must lie on both planes
        assert_relative_eq!(line.origin().y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(line.origin().z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn intersect_parallel_planes_returns_none() {
        let p1 = Plane3::new(Vector3::z_axis(), 1.0);
        let p2 = Plane3::new(Vector3::z_axis(), 3.0);
        assert!(p1.intersect_plane(&p2).is_none());
    }

    #[test]
    fn intersect_same_plane_returns_none() {
        let p = Plane3::xy();
        assert!(p.intersect_plane(&p).is_none());
    }

    #[test]
    fn stress_intersection_line_lies_on_both_planes() {
        let mut rg = RandomGeometry::new();
        for _ in 0..500 {
            let iso1 = rg.iso3(10.0);
            let iso2 = rg.iso3(10.0);
            let p1 = Plane3::xy().new_transformed_by(&iso1);
            let p2 = Plane3::xy().new_transformed_by(&iso2);

            if let Some(line) = p1.intersect_plane(&p2) {
                for t in [-5.0, -1.0, 0.0, 1.0, 5.0] {
                    let pt = line.at(t);
                    assert_relative_eq!(
                        p1.signed_distance_to_point(&pt).abs(),
                        0.0,
                        epsilon = 1e-8
                    );
                    assert_relative_eq!(
                        p2.signed_distance_to_point(&pt).abs(),
                        0.0,
                        epsilon = 1e-8
                    );
                }
            }
        }
    }
}
