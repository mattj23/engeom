pub mod align3;
mod circle3;
mod curve3;
mod iso3;
mod line3;
mod manifold;
pub mod mesh;
mod plane3;
pub mod point_cloud;
mod sphere3;
mod unroll_transform;
mod xyzwpr;

use parry3d_f64::na::UnitQuaternion;

use crate::TransformBy;
use crate::common::surface_point::{SurfacePoint, SurfacePointCollection};
use crate::common::svd_basis::SvdBasis;
pub use circle3::Circle3;
pub use curve3::{Curve3, CurveStation3};
pub use iso3::IsoExtensions3;
pub use line3::Line3;
pub use manifold::Manifold1Pos3;
pub use mesh::{Mesh, MeshCollisionSet, UvMapping};
use parry3d_f64::query::Ray;
pub use plane3::Plane3;
pub use point_cloud::{PointCloud, PointCloudFeatures, PointCloudKdTree, PointCloudOverlap};
pub use sphere3::Sphere3;
use std::ops;
pub use unroll_transform::UnrollTransform;
pub use xyzwpr::XyzWpr;

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type SurfacePoint3 = SurfacePoint<3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;
pub type KdTree3 = crate::common::kd_tree::KdTree<3>;

pub type SvdBasis3 = SvdBasis<3>;
pub type Align3 = crate::common::align::Alignment<UnitQuaternion<f64>, 3>;

pub type Aabb3 = parry3d_f64::bounding_volume::Aabb;

impl ops::Mul<SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: SurfacePoint3) -> Self::Output {
        rhs.transformed(self)
    }
}

impl ops::Mul<&SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: &SurfacePoint3) -> Self::Output {
        rhs.transformed(self)
    }
}

impl SurfacePointCollection<3> for &Vec<SurfacePoint3> {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl SurfacePointCollection<3> for Vec<SurfacePoint3> {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl SurfacePointCollection<3> for &[SurfacePoint3] {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl TransformBy<Iso3, Vec<Point3>> for &[Point3] {
    fn transform_by(&self, transform: &Iso3) -> Vec<Point3> {
        self.iter().map(|p| transform * p).collect()
    }
}

impl TransformBy<Iso3, Vec<Point3>> for &Vec<Point3> {
    fn transform_by(&self, transform: &Iso3) -> Vec<Point3> {
        self.iter().map(|p| transform * p).collect()
    }
}

impl From<&SurfacePoint3> for Ray {
    fn from(value: &SurfacePoint3) -> Self {
        Ray::new(value.point, value.normal.into_inner())
    }
}

impl Default for SurfacePoint3 {
    fn default() -> Self {
        SurfacePoint3::new(Point3::origin(), Vector3::x_axis())
    }
}

#[cfg(test)]
pub mod tests {
    use crate::na::{SVector, Translation3, UnitQuaternion, Vector};
    use crate::{Iso3, Point3, UnitVec3, Vector3};
    use rand::distr::{Distribution, Uniform};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::f64::consts::PI;

    enum RngSource {
        Seeded(StdRng),
        Thread,
    }

    /// A helper for generating random geometric entities in tests. Can be constructed with either
    /// a fixed seed (for reproducible tests) or using the default thread RNG.
    pub struct RandomGeometry {
        rng: RngSource,
    }

    impl RandomGeometry {
        /// Create a `RandomGeometry` that uses the default thread RNG.
        pub fn new() -> Self {
            Self {
                rng: RngSource::Thread,
            }
        }

        /// Create a `RandomGeometry` seeded with a fixed value for reproducible tests.
        pub fn from_seed(seed: u64) -> Self {
            Self {
                rng: RngSource::Seeded(StdRng::seed_from_u64(seed)),
            }
        }

        pub fn f64(&mut self, lo: f64, hi: f64) -> f64 {
            let u = Uniform::new(lo, hi).unwrap();
            match &mut self.rng {
                RngSource::Seeded(r) => u.sample(r),
                RngSource::Thread => u.sample(&mut rand::rng()),
            }
        }

        pub fn f64_sym(&mut self, hi: f64) -> f64 {
            let u = Uniform::new(-hi, hi).unwrap();
            match &mut self.rng {
                RngSource::Seeded(r) => u.sample(r),
                RngSource::Thread => u.sample(&mut rand::rng()),
            }
        }

        pub fn vector<const D: usize>(&mut self, limit: f64) -> SVector<f64, D> {
            let mut v = SVector::zeros();
            for i in 0..D {
                v[i] = self.f64(-limit, limit);
            }
            v
        }

        /// Returns a random `Iso3` with translation components in `[-10, 10]` and arbitrary
        /// rotation.
        pub fn iso3(&mut self, t: f64) -> Iso3 {
            let tx = self.f64(-t, t);
            let ty = self.f64(-t, t);
            let tz = self.f64(-t, t);
            let rx = self.f64(-PI, PI);
            let ry = self.f64(-PI, PI);
            let rz = self.f64(-PI, PI);
            Iso3::from_parts(
                Translation3::from(Vector3::new(tx, ty, tz)),
                UnitQuaternion::from_euler_angles(rx, ry, rz),
            )
        }

        /// Returns a random `Point3` with each component in `[-limit, limit]`.
        pub fn point3(&mut self, limit: f64) -> Point3 {
            Point3::new(
                self.f64(-limit, limit),
                self.f64(-limit, limit),
                self.f64(-limit, limit),
            )
        }

        /// Returns a random `Vector3` with each component in `[-limit, limit]`.
        pub fn vector3(&mut self, limit: f64) -> Vector3 {
            Vector3::new(
                self.f64(-limit, limit),
                self.f64(-limit, limit),
                self.f64(-limit, limit),
            )
        }

        /// Returns a random unit vector.
        pub fn unit_vec3(&mut self) -> UnitVec3 {
            self.iso3(1.0).rotation * Vector3::x_axis()
        }
    }
}
