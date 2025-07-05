use crate::{Point3, Vector3};
use three_d::{CpuMaterial, Srgba, Vec3};

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

pub fn cpu_mat(
    red: u8,
    green: u8,
    blue: u8,
    alpha: u8,
    roughness: f32,
    metallic: f32,
) -> CpuMaterial {
    let albedo = Srgba::new(red, green, blue, alpha);
    CpuMaterial {
        albedo,
        roughness,
        metallic,
        ..Default::default()
    }
}

pub enum ModState {
    None,
    ShiftOnly,
    CtrlOnly,
    AltOnly,
    ShiftCtrl,
    ShiftAlt,
    CtrlAlt,
    ShiftCtrlAlt,
}

pub fn mod_state(modifiers: &three_d::Modifiers) -> ModState {
    match (modifiers.shift, modifiers.ctrl, modifiers.alt) {
        (false, false, false) => ModState::None,
        (true, false, false) => ModState::ShiftOnly,
        (false, true, false) => ModState::CtrlOnly,
        (false, false, true) => ModState::AltOnly,
        (true, true, false) => ModState::ShiftCtrl,
        (true, false, true) => ModState::ShiftAlt,
        (false, true, true) => ModState::CtrlAlt,
        (true, true, true) => ModState::ShiftCtrlAlt,
    }
}
