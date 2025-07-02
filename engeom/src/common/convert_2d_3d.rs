//! This module contains two traits, `To2D` and `To3D`, which are used to convert between 2D and
//! 3D constructs by dropping or adding a Z component.

use crate::common::SurfacePoint;
use crate::geom2::{Point2, UnitVec2, Vector2};
use crate::geom3::Point3;

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

impl To2D for crate::geom3::UnitVec3 {
    type T2D = UnitVec2;

    /// Converts a 3D unit vector to a 2D unit vector by dropping the Z component and re-normalizing
    /// the vector so that it has a magnitude of 1 in the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        UnitVec2::new_normalize(Vector2::new(self.x, self.y))
    }
}

impl To2D for crate::geom3::Point3 {
    type T2D = Point2;

    /// Converts a 3D point to a 2D point by dropping the Z component, effectively projecting it
    /// to the X-Y plane.
    fn to_2d(&self) -> Self::T2D {
        Point2::new(self.x, self.y)
    }
}

impl To2D for crate::geom3::Vector3 {
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

impl To3D for Point2 {
    type T3D = crate::geom3::Point3;

    /// Converts a 2D point to a 3D point by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        crate::geom3::Point3::new(self.x, self.y, 0.0)
    }
}

impl To3D for Vector2 {
    type T3D = crate::geom3::Vector3;

    /// Converts a 2D vector to a 3D vector by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        crate::geom3::Vector3::new(self.x, self.y, 0.0)
    }
}

impl To3D for UnitVec2 {
    type T3D = crate::geom3::UnitVec3;

    /// Converts a 2D unit vector to a 3D unit vector by adding a zero-valued Z component.
    fn to_3d(&self) -> Self::T3D {
        crate::geom3::UnitVec3::new_normalize(self.into_inner().to_3d())
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
