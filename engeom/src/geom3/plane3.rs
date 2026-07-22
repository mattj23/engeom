use crate::common::PCoords;
use crate::common::svd_basis::SvdBasis;
use crate::geom3::UnitVec3;
use crate::geom3::line3::Line3;
use crate::{Iso3, Point3, Result, SurfacePoint3, Vector3};
use std::ops;

/// A plane in 3D space, defined by a unit normal and a signed offset from the origin.
///
/// This is one of `engeom`'s 3D geometric primitives.
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

    /// Create a new plane from a unit vector and an offset from the origin. The unit vector
    /// components ux, uy, and uz are the equivalent to a, b, and c in the traditional a, b, c, d
    /// representation of a plane.
    ///
    /// # Arguments
    ///
    /// * `normal`: The plane normal
    /// * `d`: The distance between the plane and the origin in the direction of the plane normal
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new(normal: UnitVec3, d: f64) -> Self {
        Self { normal, d }
    }

    /// Fit a plane to a set of points using singular value decomposition, resulting in a
    /// least-squares fitting. Optional weights may be provided in a slice of `f64` with the same
    /// number of elements as `points`, where the weight `i` corresponds with the point `i`.
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of coordinates to fit the plane to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Plane3, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_fit(points: &[impl PCoords<3>], weights: Option<&[f64]>) -> Result<Self> {
        let basis = SvdBasis::from_points(points, weights)
            .ok_or("Failed to fit plane with singular value decomposition")?;
        Ok(Plane3::from_point_normal(&basis.center, &basis.smallest()))
    }

    /// Create a Plane3 from three points, with the normal following the right-hand rule from
    /// `p1` to `p2` to `p3`. Returns an error if the points are collinear (or coincident).
    ///
    /// # Arguments
    ///
    /// * `p1`: the first point
    /// * `p2`: the second point
    /// * `p3`: the third point
    ///
    /// returns: Result<Plane3>
    pub fn from_3_points(p1: &Point3, p2: &Point3, p3: &Point3) -> Result<Self> {
        let cross = (p2 - p1).cross(&(p3 - p1));
        let normal = UnitVec3::try_new(cross, 1e-10).ok_or("Points are collinear")?;
        Ok(Self::from_point_normal(p1, &normal))
    }

    /// Create a Plane3 from a point on the plane and a unit normal direction.
    ///
    /// # Arguments
    ///
    /// * `point`: a point lying on the plane
    /// * `normal`: the unit normal of the plane
    ///
    /// returns: Plane3
    pub fn from_point_normal(point: &Point3, normal: &UnitVec3) -> Self {
        let d = normal.dot(&point.coords);
        Self::new(*normal, d)
    }

    /// Create a Plane3 from a `SurfacePoint3`, using its point and normal directly.
    ///
    /// # Arguments
    ///
    /// * `surface_point`: the surface point to create the plane from
    ///
    /// returns: Plane3
    pub fn from_surface_point(surface_point: &SurfacePoint3) -> Self {
        Self::from_point_normal(&surface_point.point, &surface_point.normal)
    }

    /// Returns a new plane in the same position as this one, but with the normal direction
    /// reversed, without modifying the original.
    pub fn normal_reversed(&self) -> Self {
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

    pub fn intersect_distance(&self, sp: &SurfacePoint3) -> Option<f64> {
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
    /// let moved = plane.offset_by(2.0);
    ///
    /// assert_relative_eq!(moved.signed_distance_to_point(&Point3::origin()), 3.0, epsilon = 1e-6);
    /// ```
    pub fn offset_by(&self, shift: f64) -> Self {
        Self::new(self.normal, self.d + shift)
    }

    /// Returns a new plane transformed by the given isometry, without modifying the original.
    pub fn transformed_by(&self, iso: &Iso3) -> Self {
        let pos = self.normal.into_inner() * self.d;
        let repr = SurfacePoint3::new(pos.into(), self.normal);
        let new_repr = repr.transformed(iso);
        Self::from_surface_point(&new_repr)
    }
}

impl ops::Mul<Plane3> for Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: Plane3) -> Plane3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<&Plane3> for Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: &Plane3) -> Plane3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<Plane3> for &Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: Plane3) -> Plane3 {
        rhs.transformed_by(self)
    }
}

impl ops::Mul<&Plane3> for &Iso3 {
    type Output = Plane3;
    fn mul(self, rhs: &Plane3) -> Plane3 {
        rhs.transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::tests::RandomGeometry3;
    use approx::assert_relative_eq;

    #[test]
    fn normal_reversed_negates_normal_and_distance() {
        let plane = Plane3::new(Vector3::z_axis(), 2.0);
        let reversed = plane.normal_reversed();
        assert_relative_eq!(reversed.normal.into_inner(), -plane.normal.into_inner());
        assert_relative_eq!(reversed.d, -plane.d);

        // The plane occupies the same position in space: signed distance is negated everywhere.
        let mut rg = RandomGeometry3::new();
        for _ in 0..100 {
            let p = rg.point3(10.0);
            assert_relative_eq!(
                reversed.signed_distance_to_point(&p),
                -plane.signed_distance_to_point(&p),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn intersect_xy_xz_gives_x_axis() {
        // xy-plane (z=0) ∩ xz-plane (y=0) should be the X axis
        let line = Plane3::xy().intersect_plane(&Plane3::xz()).unwrap();
        // direction must be parallel to X
        let dir = line.direction.normalize();
        assert_relative_eq!(dir.x.abs(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(dir.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(dir.z, 0.0, epsilon = 1e-12);
        // origin must lie on both planes
        assert_relative_eq!(line.origin.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(line.origin.z, 0.0, epsilon = 1e-12);
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
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let iso1 = rg.iso3(10.0);
            let iso2 = rg.iso3(10.0);
            let p1 = Plane3::xy().transformed_by(&iso1);
            let p2 = Plane3::xy().transformed_by(&iso2);

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

    // -------------------------------------------------------------------------
    // from_fit tests
    // -------------------------------------------------------------------------

    #[test]
    fn from_fit_recovers_flat_plane() {
        let points = [
            Point3::new(0.0, 0.0, 5.0),
            Point3::new(1.0, 0.0, 5.0),
            Point3::new(0.0, 1.0, 5.0),
            Point3::new(1.0, 1.0, 5.0),
        ];
        let plane = Plane3::from_fit(&points, None).unwrap();
        assert_relative_eq!(plane.normal.into_inner().z.abs(), 1.0, epsilon = 1e-10);
        for p in &points {
            assert_relative_eq!(plane.distance_to_point(p), 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() {
        let points = [
            Point3::new(0.0, 0.0, 5.0),
            Point3::new(1.0, 0.3, 5.2),
            Point3::new(0.2, 1.0, 4.8),
            Point3::new(1.0, 1.0, 5.1),
        ];
        let unweighted = Plane3::from_fit(&points, None).unwrap();
        let weights = [1.0; 4];
        let weighted = Plane3::from_fit(&points, Some(&weights)).unwrap();
        assert_relative_eq!(
            weighted.normal.into_inner(),
            unweighted.normal.into_inner(),
            epsilon = 1e-10
        );
        assert_relative_eq!(weighted.d, unweighted.d, epsilon = 1e-10);
    }

    #[test]
    fn from_fit_heavily_weighted_point_pulls_plane_toward_it() {
        // A mostly-flat cluster at z=0 plus one outlier point; in the limit of a very large
        // weight on the outlier, the least-squares fit is forced to pass through it almost
        // exactly.
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.3, 0.7, 0.0),
            Point3::new(0.8, 0.2, 0.0),
            Point3::new(3.0, 3.0, 2.0),
        ];
        let outlier = points[6];

        let unweighted = Plane3::from_fit(&points, None).unwrap();
        let mut weights = [1.0; 7];
        weights[6] = 1000.0;
        let weighted = Plane3::from_fit(&points, Some(&weights)).unwrap();

        assert!(weighted.distance_to_point(&outlier) < unweighted.distance_to_point(&outlier));
        assert_relative_eq!(weighted.distance_to_point(&outlier), 0.0, epsilon = 1e-2);
    }

    #[test]
    fn from_fit_empty_points_is_error() {
        let points: [Point3; 0] = [];
        assert!(Plane3::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_insufficient_points_is_error() {
        // Two points can't uniquely determine a plane orientation; the SVD degenerates.
        let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        assert!(Plane3::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_coincident_points_is_error() {
        let points = [Point3::new(1.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)];
        assert!(Plane3::from_fit(&points, None).is_err());
    }

    #[test]
    fn stress_from_fit_recovers_known_plane() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..200 {
            let normal = rg.unit_vec3();
            let point_on_plane = rg.point3(10.0);
            let true_plane = Plane3::from_point_normal(&point_on_plane, &normal);

            // Build two orthonormal in-plane basis vectors to generate points that lie exactly
            // on the plane.
            let n = normal.into_inner();
            let reference = if n.z.abs() < 0.9 {
                Vector3::z()
            } else {
                Vector3::x()
            };
            let u = reference.cross(&n).normalize();
            let v = n.cross(&u);

            let points: Vec<Point3> = (0..12)
                .map(|_| point_on_plane + u * rg.f64(-5.0, 5.0) + v * rg.f64(-5.0, 5.0))
                .collect();

            let fit = Plane3::from_fit(&points, None).unwrap();
            for p in &points {
                assert_relative_eq!(fit.distance_to_point(p), 0.0, epsilon = 1e-8);
            }
            assert_relative_eq!(
                fit.normal.dot(&true_plane.normal).abs(),
                1.0,
                epsilon = 1e-8
            );
        }
    }
}
