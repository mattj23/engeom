//! This module contains two traits, `To2D` and `To3D`, which are used to convert between 2D and
//! 3D constructs by dropping or adding a Z component.

use crate::common::SurfacePoint;
use crate::geom2::{Point2, UnitVec2, Vector2};
use crate::geom3::{Point3, UnitVec3, Vector3};

/// A trait for converting a 3D construct to a 2D construct by dropping the Z component.
pub trait To2D {
    type T2D;

    fn to_2d(&self) -> Self::T2D;
}

impl To2D for &[Point3] {
    type T2D = Vec<Point2>;

    /// Converts a slice of 3D points to a slice of 2D points by dropping the Z component of each
    /// point.
    fn to_2d(&self) -> Self::T2D {
        self.iter().map(|p| p.to_2d()).collect()
    }
}

impl To2D for Vec<Point3> {
    type T2D = Vec<Point2>;

    /// Converts a vector of 3D points to a vector of 2D points by dropping the Z component of each
    /// point.
    fn to_2d(&self) -> Self::T2D {
        self.iter().map(|p| p.to_2d()).collect()
    }
}

impl To2D for UnitVec3 {
    type T2D = UnitVec2;

    /// Converts a 3D unit vector to a 2D unit vector by dropping the Z component and re-normalizing
    /// the vector so that it has a magnitude of 1 in the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        UnitVec2::new_normalize(Vector2::new(self.x, self.y))
    }
}

impl To2D for Point3 {
    type T2D = Point2;

    /// Converts a 3D point to a 2D point by dropping the Z component, effectively projecting it
    /// to the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        Point2::new(self.x, self.y)
    }
}

impl To2D for Vector3 {
    type T2D = Vector2;

    /// Converts a 3D vector to a 2D vector by dropping the Z component, effectively projecting it
    /// into the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        Vector2::new(self.x, self.y)
    }
}

impl To2D for SurfacePoint<3> {
    type T2D = SurfacePoint<2>;

    /// Converts a 3D surface point to a 2D surface point by dropping the Z component of both
    /// the point and the normal, and re-normalizing the normal vector so that it has a magnitude
    /// of 1 in the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        let p0 = Point2::new(self.point.x, self.point.y);
        let n0 = UnitVec2::new_normalize(Vector2::new(self.normal.x, self.normal.y));
        Self::T2D::new(p0, n0)
    }
}

impl To2D for &[SurfacePoint<3>] {
    type T2D = Vec<SurfacePoint<2>>;

    /// Converts a slice of 3D surface points to a slice of 2D surface points by dropping the Z
    /// component of both the point and the normal, and re-normalizing the normal vector so that
    /// it has a magnitude of 1 in the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        self.iter().map(|p| p.to_2d()).collect()
    }
}

impl To2D for &Vec<SurfacePoint<3>> {
    type T2D = Vec<SurfacePoint<2>>;

    /// Converts a vector of 3D surface points to a vector of 2D surface points by dropping the Z
    /// component of both the point and the normal, and re-normalizing the normal vector so that
    /// it has a magnitude of 1 in the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        self.iter().map(|p| p.to_2d()).collect()
    }
}

// ================================================================================================

/// A trait for converting a 2D construct to a 3D construct by adding a zero-valued Z component.
pub trait To3D {
    type T3D;

    fn to_3d(&self) -> Self::T3D;
}

impl To3D for &[Point2] {
    type T3D = Vec<Point3>;

    /// Converts a slice of 2D points to a slice of 3D points by adding a zero-valued Z component
    /// to each point.
    fn to_3d(&self) -> Self::T3D {
        self.iter().map(|p| p.to_3d()).collect()
    }
}

impl To3D for Vec<Point2> {
    type T3D = Vec<Point3>;

    /// Converts a slice of 2D points to a slice of 3D points by adding a zero-valued Z component
    /// to each point.
    fn to_3d(&self) -> Self::T3D {
        self.iter().map(|p| p.to_3d()).collect()
    }
}

impl To3D for Point2 {
    type T3D = Point3;

    /// Converts a 2D point to a 3D point by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        Point3::new(self.x, self.y, 0.0)
    }
}

impl To3D for Vector2 {
    type T3D = Vector3;

    /// Converts a 2D vector to a 3D vector by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        Vector3::new(self.x, self.y, 0.0)
    }
}

impl To3D for UnitVec2 {
    type T3D = UnitVec3;

    /// Converts a 2D unit vector to a 3D unit vector by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        UnitVec3::new_normalize(self.into_inner().to_3d())
    }
}

impl To3D for SurfacePoint<2> {
    type T3D = SurfacePoint<3>;

    /// Converts a 2D surface point to a 3D surface point by adding a zero-valued Z component to
    /// both the point and the normal.
    fn to_3d(&self) -> Self::T3D {
        let p0 = self.point.to_3d();
        let n0 = self.normal.to_3d();
        Self::T3D::new(p0, n0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn point3_to_2d_drops_z() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let p2 = p.to_2d();
        assert_relative_eq!(p2.x, 1.0);
        assert_relative_eq!(p2.y, 2.0);
    }

    #[test]
    fn vector3_to_2d_drops_z() {
        let v = Vector3::new(3.0, 4.0, 5.0);
        let v2 = v.to_2d();
        assert_relative_eq!(v2.x, 3.0);
        assert_relative_eq!(v2.y, 4.0);
    }

    #[test]
    fn unit_vec3_to_2d_renormalizes() {
        // A vector pointing at 45° in XY with a Z component — result should be unit length
        let uv3 = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 10.0));
        let uv2 = uv3.to_2d();
        assert_relative_eq!(uv2.norm(), 1.0, epsilon = 1e-10);
        // X and Y components should be equal (both came from 1.0)
        assert_relative_eq!(uv2.x, uv2.y, epsilon = 1e-10);
    }

    #[test]
    fn slice_of_point3_to_2d() {
        let pts: Vec<Point3> = vec![Point3::new(1.0, 2.0, 3.0), Point3::new(4.0, 5.0, 6.0)];
        let pts2 = pts.as_slice().to_2d();
        assert_eq!(pts2.len(), 2);
        assert_relative_eq!(pts2[0].x, 1.0);
        assert_relative_eq!(pts2[0].y, 2.0);
        assert_relative_eq!(pts2[1].x, 4.0);
        assert_relative_eq!(pts2[1].y, 5.0);
    }

    #[test]
    fn vec_of_point3_to_2d() {
        let pts: Vec<Point3> = vec![Point3::new(7.0, 8.0, 9.0)];
        let pts2 = pts.to_2d();
        assert_eq!(pts2.len(), 1);
        assert_relative_eq!(pts2[0].x, 7.0);
        assert_relative_eq!(pts2[0].y, 8.0);
    }

    #[test]
    fn surface_point3_to_2d_drops_z_and_renormalizes() {
        let sp = SurfacePoint::<3>::new_normalize(
            Point3::new(1.0, 2.0, 3.0),
            Vector3::new(0.0, 0.0, 1.0),
        );
        // Normal is (0,0,1) — projecting to XY gives (0,0), which is degenerate but the
        // function still calls new_normalize; test a non-degenerate case instead.
        let sp2 = SurfacePoint::<3>::new_normalize(
            Point3::new(1.0, 2.0, 3.0),
            Vector3::new(1.0, 0.0, 0.5),
        )
        .to_2d();
        assert_relative_eq!(sp2.point.x, 1.0);
        assert_relative_eq!(sp2.point.y, 2.0);
        assert_relative_eq!(sp2.normal.norm(), 1.0, epsilon = 1e-10);
        // Normal should point in +X direction after dropping Z and renormalizing
        assert_relative_eq!(sp2.normal.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(sp2.normal.y, 0.0, epsilon = 1e-10);
        let _ = sp; // suppress unused warning
    }

    #[test]
    fn slice_of_surface_point3_to_2d() {
        let pts = vec![
            SurfacePoint::<3>::new_normalize(
                Point3::new(0.0, 1.0, 2.0),
                Vector3::new(1.0, 0.0, 0.0),
            ),
            SurfacePoint::<3>::new_normalize(
                Point3::new(3.0, 4.0, 5.0),
                Vector3::new(0.0, 1.0, 0.0),
            ),
        ];
        let pts2 = pts.as_slice().to_2d();
        assert_eq!(pts2.len(), 2);
        assert_relative_eq!(pts2[0].point.x, 0.0);
        assert_relative_eq!(pts2[0].point.y, 1.0);
        assert_relative_eq!(pts2[1].point.x, 3.0);
        assert_relative_eq!(pts2[1].point.y, 4.0);
    }

    #[test]
    fn ref_vec_of_surface_point3_to_2d() {
        let pts = vec![SurfacePoint::<3>::new_normalize(
            Point3::new(1.0, 2.0, 3.0),
            Vector3::new(1.0, 0.0, 0.0),
        )];
        let pts2 = (&pts).to_2d();
        assert_eq!(pts2.len(), 1);
        assert_relative_eq!(pts2[0].normal.norm(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn point2_to_3d_adds_zero_z() {
        let p = Point2::new(3.0, 4.0);
        let p3 = p.to_3d();
        assert_relative_eq!(p3.x, 3.0);
        assert_relative_eq!(p3.y, 4.0);
        assert_relative_eq!(p3.z, 0.0);
    }

    #[test]
    fn vector2_to_3d_adds_zero_z() {
        let v = Vector2::new(5.0, 6.0);
        let v3 = v.to_3d();
        assert_relative_eq!(v3.x, 5.0);
        assert_relative_eq!(v3.y, 6.0);
        assert_relative_eq!(v3.z, 0.0);
    }

    #[test]
    fn unit_vec2_to_3d_stays_unit_length() {
        let uv2 = UnitVec2::new_normalize(Vector2::new(1.0, 0.0));
        let uv3 = uv2.to_3d();
        assert_relative_eq!(uv3.norm(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(uv3.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(uv3.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(uv3.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn slice_of_point2_to_3d() {
        let pts: Vec<Point2> = vec![Point2::new(1.0, 2.0), Point2::new(3.0, 4.0)];
        let pts3 = pts.as_slice().to_3d();
        assert_eq!(pts3.len(), 2);
        assert_relative_eq!(pts3[0].z, 0.0);
        assert_relative_eq!(pts3[1].x, 3.0);
        assert_relative_eq!(pts3[1].z, 0.0);
    }

    #[test]
    fn vec_of_point2_to_3d() {
        let pts: Vec<Point2> = vec![Point2::new(9.0, 8.0)];
        let pts3 = pts.to_3d();
        assert_eq!(pts3.len(), 1);
        assert_relative_eq!(pts3[0].x, 9.0);
        assert_relative_eq!(pts3[0].y, 8.0);
        assert_relative_eq!(pts3[0].z, 0.0);
    }

    #[test]
    fn surface_point2_to_3d_adds_zero_z() {
        let sp = SurfacePoint::<2>::new_normalize(Point2::new(1.0, 2.0), Vector2::new(0.0, 1.0));
        let sp3 = sp.to_3d();
        assert_relative_eq!(sp3.point.x, 1.0);
        assert_relative_eq!(sp3.point.y, 2.0);
        assert_relative_eq!(sp3.point.z, 0.0);
        assert_relative_eq!(sp3.normal.norm(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(sp3.normal.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(sp3.normal.y, 1.0, epsilon = 1e-10);
        assert_relative_eq!(sp3.normal.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn point3_round_trip_xy_preserved() {
        let p = Point3::new(1.5, -2.5, 99.0);
        let rt = p.to_2d().to_3d();
        assert_relative_eq!(rt.x, p.x);
        assert_relative_eq!(rt.y, p.y);
        assert_relative_eq!(rt.z, 0.0);
    }

    #[test]
    fn surface_point_round_trip_xy_preserved() {
        let sp = SurfacePoint::<3>::new_normalize(
            Point3::new(2.0, 3.0, 4.0),
            Vector3::new(1.0, 0.0, 0.0),
        );
        let rt = sp.to_2d().to_3d();
        assert_relative_eq!(rt.point.x, 2.0);
        assert_relative_eq!(rt.point.y, 3.0);
        assert_relative_eq!(rt.point.z, 0.0);
        assert_relative_eq!(rt.normal.norm(), 1.0, epsilon = 1e-10);
    }
}
