use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::common::svd_basis::SvdBasis;
use crate::geom3::line3::Line3;
use crate::geom3::{IsoExtensions3, UnitVec3};
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

    /// Fit a plane to a set of points using MAGSAC++ robust consensus estimation.
    ///
    /// Unlike an ordinary least-squares fit ([`Plane3::from_fit`]), this rejects gross outliers by
    /// taking an upper bound on the inlier noise (`sigma_max`) rather than a hard inlier/outlier
    /// threshold, and refines each candidate with noise-marginalized iteratively reweighted least
    /// squares. It is substantially less sensitive to `sigma_max` than RANSAC is to its threshold,
    /// as long as `sigma_max` is not chosen smaller than the actual noise.
    ///
    /// The resulting plane passes through the centroid of the inlier set with a unit normal. The
    /// direction of the normal is not meaningful (it may point to either side of the plane).
    ///
    /// # Arguments
    ///
    /// * `points`: the points to fit the plane to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `options`: an optional [`Magsac`] configuration to override the iteration count, refinement
    ///   steps, confidence, or RNG seed. Its `sigma_max` field is overridden by the `sigma_max`
    ///   argument.
    ///
    /// returns: Result<Plane3, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point3],
        sigma_max: f64,
        options: Option<Magsac>,
    ) -> Result<Self> {
        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit = magsac.fit::<3, Plane3>(points)?;
        Ok(fit.model)
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

    /// Computes a deterministic orthonormal coordinate frame lying in this plane, as the isometry
    /// which takes a point expressed in that frame to where it lies in the world.
    ///
    /// This is what gives the plane a two-dimensional coordinate system, so that geometry can be
    /// brought into it and back out again reproducibly. Use the inverse to go from the world into
    /// the plane, where a point on the plane has a zero z coordinate.
    ///
    /// The frame is built so that:
    ///
    /// - its origin is the point of the plane closest to the world origin,
    /// - its z axis is the plane normal, so `x`, `y`, `z` stay right-handed, and
    /// - its x axis is the world x axis projected onto the plane, falling back to the world y axis
    ///   when the plane normal is too close to x for that to be meaningful.
    ///
    /// The consequence worth knowing is that **`Plane3::xy().compute_frame()` is the identity**, so
    /// bringing geometry into the x-y plane is exactly dropping the z coordinate, and doing it to
    /// any plane `z = k` with a `+z` normal leaves x and y untouched. A plane with a `-z` normal is
    /// instead seen from its own normal side, which negates y; that is not an accident, it is what
    /// keeps the frame right-handed and the inside/outside sense of a section intact.
    ///
    /// Unlike [`crate::geom3::IsoExtensions3::from_z_arbitrary_xy`], which promises nothing about
    /// where x and y land, this is stable for a given plane and is safe to depend on.
    ///
    /// The frame is derived from the plane's fields rather than stored, so a caller using it in a
    /// loop should compute it once and keep it.
    ///
    /// returns: Isometry<f64, Unit<Quaternion<f64>>, 3>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Plane3, Point3};
    ///
    /// // Into the plane `z = 5`, whose frame differs from the world only by that offset.
    /// let plane = Plane3::xy().offset_by(5.0);
    /// let to_plane = plane.compute_frame().inverse();
    ///
    /// let in_plane = to_plane * Point3::new(2.0, 3.0, 5.0);
    /// assert_relative_eq!(in_plane, Point3::new(2.0, 3.0, 0.0), epsilon = 1e-12);
    /// ```
    pub fn compute_frame(&self) -> Iso3 {
        let n = self.normal.into_inner();

        // Chosen by the normal rather than by testing the projection afterwards, so that the
        // vector being normalized is never a near-zero remainder: this keeps `|u|` at or above
        // 0.43 in every case.
        let reference = if n.x.abs() < 0.9 {
            Vector3::x()
        } else {
            Vector3::y()
        };

        let u = self.project_vector(&reference);
        let v = n.cross(&u);
        let origin = Point3::from(n * self.d);

        // `u` and `v` are perpendicular, both well away from zero, and neither is parallel to the
        // other, so the only way this fails is if the invariant above is broken.
        Iso3::from_basis_xy(&u, &v, Some(origin))
            .expect("a plane frame is built from two non-degenerate perpendicular vectors")
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
        let new_repr = repr.transformed_by(iso);
        Self::from_surface_point(&new_repr)
    }
}

impl ConsensusModel<3> for Plane3 {
    type Point = Point3;
    const SAMPLE_SIZE: usize = 3;

    fn from_sample(sample: &[Point3]) -> Option<Self> {
        // `from_3_points` already rejects collinear (and coincident) samples.
        Plane3::from_3_points(&sample[0], &sample[1], &sample[2]).ok()
    }

    fn residual(&self, point: &Point3) -> f64 {
        // Signed distance keeps the residual smooth through the plane for the least-squares
        // refinement; only its magnitude is used for scoring.
        self.signed_distance_to_point(point)
    }

    fn refine_weighted(points: &[Point3], weights: &[f64], _initial: &Self) -> Option<Self> {
        // The MAGSAC++ refinement step is a single weighted least-squares fit, which for a plane is
        // exactly the weighted SVD fit provided by `from_fit`.
        Plane3::from_fit(points, Some(weights)).ok()
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
    use crate::common::consensus::Magsac;
    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;

    // ============================================================================================
    // The in-plane coordinate frame
    // ============================================================================================

    /// The whole point of the frame: projecting onto the x-y plane has to be the same operation as
    /// dropping the z coordinate, or nothing built on top of it lines up with `To2D`.
    #[test]
    fn the_xy_plane_frame_is_the_identity() {
        assert_relative_eq!(
            Plane3::xy().compute_frame().to_matrix(),
            Iso3::identity().to_matrix(),
            epsilon = 1e-12
        );
    }

    /// A plane parallel to x-y differs from the world only by the offset, so x and y have to come
    /// through a round trip untouched.
    #[test]
    fn an_offset_xy_plane_leaves_x_and_y_alone() {
        let plane = Plane3::xy().offset_by(5.0);
        let to_plane = plane.compute_frame().inverse();

        let p = to_plane * Point3::new(2.0, 3.0, 5.0);
        assert_relative_eq!(p, Point3::new(2.0, 3.0, 0.0), epsilon = 1e-12);
    }

    /// Reversing the normal looks at the plane from the other side, which has to negate exactly one
    /// in-plane axis. If it negated both, the frame would be a rotation rather than a reflection and
    /// a section's inside/outside sense would silently flip.
    #[test]
    fn a_reversed_normal_negates_one_axis() {
        let frame = Plane3::xy().normal_reversed().compute_frame();

        assert_relative_eq!(
            frame * Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            frame * Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            epsilon = 1e-12
        );
    }

    /// The frame has to be a frame *of the plane*: its origin on the plane and its z axis the
    /// normal, whatever the plane is.
    #[test]
    fn stress_a_frame_sits_on_its_plane() {
        let mut rg = RandomGeometry3::from_seed(31);

        for _ in 0..500 {
            let plane = Plane3::from_point_normal(
                &Point3::new(rg.f64_sym(50.0), rg.f64_sym(50.0), rg.f64_sym(50.0)),
                &rg.unit_vec(),
            );
            let frame = plane.compute_frame();

            let origin = frame * Point3::origin();
            assert_relative_eq!(plane.signed_distance_to_point(&origin), 0.0, epsilon = 1e-9);

            // The z axis is the normal, and the other two therefore lie in the plane.
            let z = frame * Vector3::z();
            assert_relative_eq!(z, plane.normal.into_inner(), epsilon = 1e-9);

            let x = frame * Vector3::x();
            let y = frame * Vector3::y();
            assert_relative_eq!(x.dot(&plane.normal), 0.0, epsilon = 1e-9);
            assert_relative_eq!(y.dot(&plane.normal), 0.0, epsilon = 1e-9);

            // ...and it is orthonormal and right-handed.
            assert_relative_eq!(x.dot(&y), 0.0, epsilon = 1e-9);
            assert_relative_eq!(x.cross(&y), z, epsilon = 1e-9);
        }
    }

    /// A point on the plane has a zero z coordinate in the frame, which is what makes dropping it
    /// a projection rather than a truncation.
    #[test]
    fn stress_points_on_a_plane_have_no_local_z() {
        let mut rg = RandomGeometry3::from_seed(77);

        for _ in 0..300 {
            let plane = Plane3::from_point_normal(
                &Point3::new(rg.f64_sym(20.0), rg.f64_sym(20.0), rg.f64_sym(20.0)),
                &rg.unit_vec(),
            );
            let to_plane = plane.compute_frame().inverse();

            let on_plane = plane.project_point(&Point3::new(
                rg.f64_sym(20.0),
                rg.f64_sym(20.0),
                rg.f64_sym(20.0),
            ));
            assert_relative_eq!((to_plane * on_plane).z, 0.0, epsilon = 1e-9);
        }
    }

    /// The frame has to survive a normal pointing along the axis the x fallback exists for.
    #[test]
    fn a_plane_normal_to_x_falls_back_to_the_y_axis() {
        let frame = Plane3::yz().compute_frame();

        assert_relative_eq!(frame * Vector3::x(), Vector3::y(), epsilon = 1e-12);
        assert_relative_eq!(frame * Vector3::y(), Vector3::z(), epsilon = 1e-12);
        assert_relative_eq!(frame * Vector3::z(), Vector3::x(), epsilon = 1e-12);
    }

    /// Nothing about the frame may depend on how the plane was built or on call order.
    #[test]
    fn a_frame_is_reproducible() {
        let plane = Plane3::from_point_normal(
            &Point3::new(1.0, -2.0, 3.0),
            &UnitVec3::new_normalize(Vector3::new(0.3, -0.7, 0.5)),
        );

        assert_relative_eq!(
            plane.compute_frame().to_matrix(),
            plane.compute_frame().to_matrix(),
            epsilon = 1e-15
        );
    }

    /// Build two orthonormal in-plane basis vectors spanning `plane`.
    fn in_plane_basis(plane: &Plane3) -> (Vector3, Vector3) {
        let n = plane.normal.into_inner();
        let reference = if n.z.abs() < 0.9 {
            Vector3::z()
        } else {
            Vector3::x()
        };
        let u = reference.cross(&n).normalize();
        let v = n.cross(&u);
        (u, v)
    }

    /// Generate `n` points scattered across `plane` with isotropic Gaussian noise `sigma`.
    fn plane_noise(
        rg: &mut RandomGeometry3,
        plane: &Plane3,
        n: usize,
        span: f64,
        sigma: f64,
    ) -> Vec<Point3> {
        let origin = Point3::from(plane.normal.into_inner() * plane.d);
        let (u, v) = in_plane_basis(plane);
        (0..n)
            .map(|_| {
                origin
                    + u * rg.f64_sym(span)
                    + v * rg.f64_sym(span)
                    + rg.gaussian_vector::<3>(sigma)
            })
            .collect()
    }

    #[test]
    fn normal_reversed_negates_normal_and_distance() {
        let plane = Plane3::new(Vector3::z_axis(), 2.0);
        let reversed = plane.normal_reversed();
        assert_relative_eq!(reversed.normal.into_inner(), -plane.normal.into_inner());
        assert_relative_eq!(reversed.d, -plane.d);

        // The plane occupies the same position in space: signed distance is negated everywhere.
        let mut rg = RandomGeometry3::new();
        for _ in 0..100 {
            let p = rg.point(10.0);
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
            let normal = rg.unit_vec();
            let point_on_plane = rg.point(10.0);
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

    // -------------------------------------------------------------------------
    // from_consensus tests
    // -------------------------------------------------------------------------

    #[test]
    fn from_consensus_convenience_recovers_plane() {
        // A clean set of coplanar points; the convenience method should recover the plane.
        let mut rg = RandomGeometry3::from_seed(7);
        let true_plane = Plane3::from_point_normal(&Point3::new(1.0, -2.0, 3.0), &rg.unit_vec());
        let points = plane_noise(&mut rg, &true_plane, 80, 5.0, 0.0);

        let plane = Plane3::from_consensus(&points, 0.01, None).unwrap();

        for p in &points {
            assert!(plane.distance_to_point(p) < 1e-6);
        }
        assert_relative_eq!(
            plane.normal.dot(&true_plane.normal).abs(),
            1.0,
            epsilon = 1e-6
        );
    }

    #[test]
    fn from_consensus_rejects_outliers() {
        let mut rg = RandomGeometry3::from_seed(101);

        // Inliers lie on the z = 2 plane with a small amount of noise.
        let true_plane = Plane3::new(Vector3::z_axis(), 2.0);
        let inliers = plane_noise(&mut rg, &true_plane, 200, 10.0, 0.01);
        let mut points = inliers.clone();

        // A dense cluster of gross outliers well off the plane.
        let center = Point3::new(-4.0, 5.0, 9.0);
        for _ in 0..60 {
            points.push(center + rg.gaussian_vector::<3>(1.0));
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(400),
            refinement_steps: 4,
            confidence: 0.99,
            seed: Some(42),
        };
        let fit = magsac.fit::<3, Plane3>(&points).unwrap();

        // Every inlier should lie very close to the recovered plane.
        for i in &inliers {
            assert!(fit.model.distance_to_point(i) < 0.01 * 6.0);
        }

        // The recovered normal should be parallel (or anti-parallel) to the true plane's normal.
        assert_relative_eq!(
            fit.model.normal.dot(&true_plane.normal).abs(),
            1.0,
            epsilon = 1e-2
        );

        // No outlier should be classified as an inlier.
        assert!(fit.inliers.iter().all(|&i| i < inliers.len()));
    }
}
