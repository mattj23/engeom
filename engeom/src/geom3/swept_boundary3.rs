use std::f64::consts::TAU;

use crate::common::triangulation::ParallelBuilder;
use crate::common::{PCoords, transform_points};
use crate::geom2::Boundary2;
use crate::geom3::IsoExtensions3;
use crate::geom3::align3::{AlignSurfMatch3, SurfaceTarget3};
use crate::{Iso3, Mesh, Point2, Point3, Result, To2D, To3D, UnitVec3, Vector3};

/// An `ExtrudedBoundary3` is a means of representing a surface in 3D space using a 2D [`Boundary2`]
/// entity and an arbitrary 3D position, direction, and length. The surface is the set of all points
/// which the boundary passes through in an imaginary sweep of the given length and in the
/// specified direction.
pub struct ExtrudedBoundary3 {
    /// A boundary in the X-Y plane which is the basis of the sweep. It does not matter if the
    /// boundary is closed or open.
    shape: Boundary2,

    /// A transformation that brings the X-Y boundary to an arbitrary position in 3D space.
    start: Iso3,

    /// The inverse of the start isometry
    start_inv: Iso3,

    /// The length that the boundary is extruded along the `start` isometry's Z axis in the 3D
    /// world coordinates.
    length: f64,
}

impl ExtrudedBoundary3 {
    /// Create a surface representation consisting of a boundary (defined in the 2D, X-Y plane)
    /// transformed to a given position and orientation in 3D space, and then swept along by a
    /// given length.
    ///
    /// # Arguments
    ///
    /// * `shape`: the 2d boundary, defined in the X-Y plane
    /// * `start`: an isometry to transform the 2D boundary to the starting position and orientation
    ///   in 3D space. The boundary will be swept along the positive Z direction of this starting
    ///   isometry.
    /// * `length`: the distance along the `start` argument's positive Z direction to sweep the
    ///   boundary. The boundary will exist between 0.0 and `length`.
    ///
    /// returns: ExtrudedBoundary3
    pub fn new(shape: Boundary2, start: Iso3, length: f64) -> Self {
        let start_inv = start.inverse();
        Self {
            shape,
            start,
            start_inv,
            length,
        }
    }

    pub fn to_mesh(&self, tol: f64) -> Result<Mesh> {
        let points = self.shape.to_points(tol)?.to_3d();
        let mut builder = ParallelBuilder::new(points.len(), false);

        let p0 = transform_points(&points, &self.start);
        builder.push(&p0)?;

        let iso2 = self.start * Iso3::translation(0.0, 0.0, self.length);
        let p1 = transform_points(&points, &iso2);
        builder.push(&p1)?;

        let (points, faces) = builder.take();

        Ok(Mesh::new(points, faces, false))
    }

    pub fn transform_by(&mut self, iso: &Iso3) {
        self.start = iso * self.start;
        self.start_inv = self.start.inverse();
    }
}

impl SurfaceTarget3 for ExtrudedBoundary3 {
    fn align_surf_closest_to(&self, p: &impl PCoords<3>) -> AlignSurfMatch3 {
        let p = Point3::from(p.coords());

        // First we want to bring the test point into the local coordinates of the start isometry
        // so we can work with the 2D boundary.
        let lp3 = self.start_inv * p;
        let lp2 = lp3.to_2d();

        // Now we get the closest point on the boundary and create the local 3d point and normal
        let (_, m) = self.shape.at_closest_to_point(&lp2);
        let lc3 = Point3::new(m.point.x, m.point.y, lp3.z);
        let ln3 = UnitVec3::new_normalize(m.normal.into_inner().to_3d());

        // Now we figure out if we're off the boundary by either missing the z interval or by
        // projecting onto a boundary end.
        let is_off = lp3.z < 0.0
            || lp3.z > self.length
            || (!self.shape.is_closed()
                && (m.l <= f64::EPSILON || m.l >= self.shape.length() - f64::EPSILON));

        AlignSurfMatch3::new(self.start * lc3, self.start * ln3, !is_off, 1.0)
    }
}

/// A `RevolvedBoundary3` is a means of representing a surface in 3D space using a 2D [`Boundary2`]
/// entity and an arbitrary 3D position, orientation, and sweep angle.  The boundary is defined in
/// the 2D X-Y plane, and transformed to 3D space by the `start` isometry. There it is rotated
/// around the `start` isometry's Y axis, which is the equivalent of the Y axis in the original
/// 2D boundary definition.  It is rotated by angle `theta`.
pub struct RevolvedBoundary3 {
    shape: Boundary2,
    start: Iso3,
    start_inv: Iso3,
    theta: f64,
}

impl RevolvedBoundary3 {
    /// Create a surface representation consisting of a boundary (defined in the 2D, X-Y plane)
    /// transformed to a given position and orientation in 3D space, and then revolved around the
    /// Y axis of the `start` isometry by a given angle.
    ///
    /// In the 2D profile, the X coordinate is the radial distance from the Y axis, and the Y
    /// coordinate is the position along the rotation axis.
    ///
    /// # Arguments
    ///
    /// * `shape`: the 2D boundary profile, defined in the X-Y plane
    /// * `start`: an isometry to transform the 2D boundary to a starting position and orientation
    ///   in 3D space. The profile will be revolved around the positive Y axis of this isometry.
    /// * `theta`: the angle in radians to sweep the boundary around the Y axis, starting from
    ///   angle 0. The boundary will exist between 0.0 and `theta`.
    ///
    /// returns: RevolvedBoundary3
    pub fn new(shape: Boundary2, start: Iso3, theta: f64) -> Self {
        let theta = if theta.abs() > TAU {
            TAU * theta.signum()
        } else {
            theta
        };

        let start_inv = start.inverse();
        Self {
            shape,
            start,
            start_inv,
            theta,
        }
    }

    pub fn to_mesh(&self, tol: f64) -> Result<Mesh> {
        let points = self.shape.to_points(tol)?.to_3d();

        // Find the largest radius
        let r_max = points
            .iter()
            .map(|p| p.x.abs())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .ok_or("Boundary has no points".to_string())?;

        let max_theta = 2.0 * (1.0 - tol / r_max).acos();
        let n_segments = (self.theta / max_theta).ceil().max(1.0) as usize;
        let angle_step = self.theta / n_segments as f64;

        let mut builder = ParallelBuilder::new(points.len(), false);

        let p0 = transform_points(&points, &self.start);
        builder.push(&p0)?;

        for i in 0..n_segments {
            let t = (i + 1) as f64 * angle_step;
            let iso = self.start * Iso3::from_ry(t);
            let pn = transform_points(&points, &iso);
            builder.push(&pn)?;
        }

        let (points, faces) = builder.take();

        Ok(Mesh::new(points, faces, false))
    }

    pub fn transform_by(&mut self, iso: &Iso3) {
        self.start = iso * self.start;
        self.start_inv = self.start.inverse();
    }
}

impl SurfaceTarget3 for RevolvedBoundary3 {
    fn align_surf_closest_to(&self, p: &impl PCoords<3>) -> AlignSurfMatch3 {
        let p = Point3::from(p.coords());

        // Transform to local coordinates where the rotation axis is Y and the profile starts
        // in the X-Y half-plane (positive X)
        let lp3 = self.start_inv * p;

        // Cylindrical coordinates: radial distance from Y axis and azimuthal angle
        let r = lp3.x.hypot(lp3.z);
        let phi = f64::atan2(lp3.z, lp3.x).rem_euclid(TAU);

        // Query the 2D profile boundary with (radius, axial) coordinates
        let lp2 = Point2::new(r, lp3.y);
        let (_, m) = self.shape.at_closest_to_point(&lp2);

        // Reconstruct the 3D closest point and surface normal at the same azimuthal angle
        let (cos_phi, sin_phi) = (phi.cos(), phi.sin());
        let r_c = m.point.x;
        let y_c = m.point.y;
        let n = m.normal.into_inner();
        let lc3 = Point3::new(r_c * cos_phi, y_c, r_c * sin_phi);
        let ln3 = UnitVec3::new_normalize(Vector3::new(n.x * cos_phi, n.y, n.x * sin_phi));

        // Off the surface if outside the angular sweep, or at ends of an open profile
        let is_off = phi > self.theta
            || (!self.shape.is_closed()
                && (m.l <= f64::EPSILON || m.l >= self.shape.length() - f64::EPSILON));

        AlignSurfMatch3::new(self.start * lc3, self.start * ln3, !is_off, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::geom2::{BoundaryData2, BoundaryEditor};
    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;

    fn open_extruded() -> ExtrudedBoundary3 {
        let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
        data.add_seg_xy(0.0, 0.0);
        data.add_seg_xy(0.0, 1.0);
        let b = data.try_to_boundary().unwrap();
        ExtrudedBoundary3::new(b, Iso3::identity(), 1.0)
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn extruded_points_off() {
        let target = open_extruded();

        assert_eq!(false, target.align_surf_closest_to(&p(2.0, 0.0, 0.5)).is_on);
        assert_eq!(false, target.align_surf_closest_to(&p(0.0, 2.0, 0.5)).is_on);
        assert_eq!(
            false,
            target.align_surf_closest_to(&p(0.0, 0.0, -1.0)).is_on
        );
        assert_eq!(false, target.align_surf_closest_to(&p(0.0, 0.0, 2.0)).is_on);
    }

    #[test]
    fn extruded_points_on() {
        let target = open_extruded();
        assert_eq!(true, target.align_surf_closest_to(&p(0.5, 0.1, 0.5)).is_on);
        assert_eq!(true, target.align_surf_closest_to(&p(0.1, 0.5, 0.5)).is_on);
    }

    #[test]
    fn extruded_known_points() {
        let pairs = vec![
            ([0.0, 2.0, 0.5], [0.0, 1.0, 0.5], None),
            ([2.0, 0.0, 0.5], [1.0, 0.0, 0.5], None),
            ([0.5, 0.1, 0.3], [0.5, 0.0, 0.3], Some([0.0, 1.0, 0.0])),
            ([0.1, 0.5, 0.3], [0.0, 0.5, 0.3], Some([1.0, 0.0, 0.0])),
            ([0.0, 2.0, -1.0], [0.0, 1.0, -1.0], None),
            ([2.0, 0.0, -1.0], [1.0, 0.0, -1.0], None),
            ([0.0, 2.0, 2.0], [0.0, 1.0, 2.0], None),
            ([2.0, 0.0, 2.0], [1.0, 0.0, 2.0], None),
            ([-1.0, -1.0, 0.5], [0.0, 0.0, 0.5], None),
            ([-1.0, -1.0, -1.0], [0.0, 0.0, -1.0], None),
            ([-1.0, -1.0, 2.0], [0.0, 0.0, 2.0], None),
        ];

        let mut rg = RandomGeometry3::new();

        for _ in 0..100 {
            let iso = rg.iso3(10.0);

            let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
            data.add_seg_xy(0.0, 0.0);
            data.add_seg_xy(0.0, 1.0);
            let b = data.try_to_boundary().unwrap();
            let target = ExtrudedBoundary3::new(b, iso.clone(), 1.0);

            for (t, e, n) in pairs.iter() {
                let t = iso * Point3::from(*t);
                let e = iso * Point3::from(*e);

                let c = target.align_surf_closest_to(&t);
                assert_relative_eq!(c.point, e, epsilon = 1.0e-6);

                if let Some(n) = n {
                    let n = iso * UnitVec3::new_normalize(Vector3::from_column_slice(n));
                    assert_relative_eq!(c.normal, n, epsilon = 1.0e-6);
                }
            }
        }
    }

    // Profile: open segment from (1,0) to (1,1) in 2D, a cylindrical surface at radius 1,
    // height 0..1, revolved by the given theta around the Y axis.
    fn open_revolved(theta: f64) -> RevolvedBoundary3 {
        let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
        data.add_seg_xy(1.0, 1.0);
        let b = data.try_to_boundary().unwrap();
        RevolvedBoundary3::new(b, Iso3::identity(), theta)
    }

    #[test]
    fn revolved_points_on() {
        let full = open_revolved(TAU);
        // Front, side, and back of the cylinder mid-height
        assert_eq!(true, full.align_surf_closest_to(&p(2.0, 0.5, 0.0)).is_on);
        assert_eq!(true, full.align_surf_closest_to(&p(0.0, 0.5, 2.0)).is_on);
        assert_eq!(true, full.align_surf_closest_to(&p(-2.0, 0.5, 0.0)).is_on);

        let half = open_revolved(std::f64::consts::PI);
        // phi=0 and phi=PI/2 are within [0, PI]
        assert_eq!(true, half.align_surf_closest_to(&p(2.0, 0.5, 0.0)).is_on);
        assert_eq!(true, half.align_surf_closest_to(&p(0.0, 0.5, 2.0)).is_on);
    }

    #[test]
    fn revolved_points_off() {
        let full = open_revolved(TAU);
        // Off at the open profile ends (above and below the height range)
        assert_eq!(false, full.align_surf_closest_to(&p(2.0, -1.0, 0.0)).is_on);
        assert_eq!(false, full.align_surf_closest_to(&p(2.0, 2.0, 0.0)).is_on);

        let half = open_revolved(std::f64::consts::PI);
        // phi = atan2(-2, 0) = -PI/2, rem_euclid = 3*PI/2 > PI: off angular sweep
        assert_eq!(false, half.align_surf_closest_to(&p(0.0, 0.5, -2.0)).is_on);
        // phi slightly greater than PI is off the half sweep
        assert_eq!(false, half.align_surf_closest_to(&p(-1.0, 0.5, -0.1)).is_on);
    }

    #[test]
    fn revolved_known_points() {
        use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

        // (test_point, expected_closest, Option<expected_normal>)
        // Profile: open segment (1,0)→(1,1) in 2D, revolved fully (TAU) around the Y axis.
        // Rotation axis is Y; X is radial; Z is the other radial axis.
        // Note: off-surface points project onto the nearest profile point (not the query y).
        let pairs: Vec<([f64; 3], [f64; 3], Option<[f64; 3]>)> = vec![
            // +X side, mid-height
            ([2.0, 0.5, 0.0], [1.0, 0.5, 0.0], Some([1.0, 0.0, 0.0])),
            // +Z side (phi = PI/2)
            ([0.0, 0.5, 2.0], [0.0, 0.5, 1.0], Some([0.0, 0.0, 1.0])),
            // -X side (phi = PI)
            ([-2.0, 0.5, 0.0], [-1.0, 0.5, 0.0], Some([-1.0, 0.0, 0.0])),
            // Diagonal at phi = PI/4; test point r=2, closest r=1
            (
                [SQRT_2, 0.5, SQRT_2],
                [FRAC_1_SQRT_2, 0.5, FRAC_1_SQRT_2],
                Some([FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2]),
            ),
            // Below profile (off, open end): closest projects to profile start at y=0
            ([2.0, -1.0, 0.0], [1.0, 0.0, 0.0], Some([1.0, 0.0, 0.0])),
            // Above profile (off, open end): closest projects to profile end at y=1
            ([2.0, 2.0, 0.0], [1.0, 1.0, 0.0], Some([1.0, 0.0, 0.0])),
        ];

        let mut rg = RandomGeometry3::new();

        for _ in 0..100 {
            let iso = rg.iso3(10.0);

            let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
            data.add_seg_xy(1.0, 1.0);
            let b = data.try_to_boundary().unwrap();
            let target = RevolvedBoundary3::new(b, iso.clone(), TAU);

            for (t, e, n) in pairs.iter() {
                let t = iso * Point3::from(*t);
                let e = iso * Point3::from(*e);

                let c = target.align_surf_closest_to(&t);
                assert_relative_eq!(c.point, e, epsilon = 1.0e-6);

                if let Some(n) = n {
                    let n = iso * UnitVec3::new_normalize(Vector3::from_column_slice(n));
                    assert_relative_eq!(c.normal, n, epsilon = 1.0e-6);
                }
            }
        }
    }
}
