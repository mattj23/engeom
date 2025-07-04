use three_d::Vec3;
use crate::{Point3, Vector3};

pub trait ToCgVec3 {
    fn to_cg(&self) -> Vec3;
}

impl ToCgVec3 for Vector3 {
    fn to_cg(&self) -> Vec3 {
        Vec3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

impl ToCgVec3 for Point3 {
    fn to_cg(&self) -> Vec3 {
        Vec3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

pub trait ToEngeom3 {
    fn to_engeom(&self) -> Vector3;
}

impl ToEngeom3 for Vec3 {
    fn to_engeom(&self) -> Vector3 {
        Vector3::new(self.x as f64, self.y as f64, self.z as f64)
    }
}
