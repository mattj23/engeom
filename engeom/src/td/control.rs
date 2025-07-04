use crate::na::UnitQuaternion;
use crate::td::{ToCgVec3, ToEngeom3};
use crate::{Iso3, Point3, UnitVec3, Vector3};
use three_d::{MouseButton, radians};

enum ModifierState {
    None,
    ShiftOnly,
    CtrlOnly,
    AltOnly,
    ShiftCtrl,
    ShiftAlt,
    CtrlAlt,
    ShiftCtrlAlt,
}

fn mod_state(modifiers: &three_d::Modifiers) -> ModifierState {
    match (modifiers.shift, modifiers.ctrl, modifiers.alt) {
        (false, false, false) => ModifierState::None,
        (true, false, false) => ModifierState::ShiftOnly,
        (false, true, false) => ModifierState::CtrlOnly,
        (false, false, true) => ModifierState::AltOnly,
        (true, true, false) => ModifierState::ShiftCtrl,
        (true, false, true) => ModifierState::ShiftAlt,
        (false, true, true) => ModifierState::CtrlAlt,
        (true, true, true) => ModifierState::ShiftCtrlAlt,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CameraControl {
    move_speed: f64,
    rot_speed: f64,
    view: Iso3,

    change_flag: bool,
}

impl CameraControl {
    /// Creates a new camera control.
    pub fn new(move_speed: f64, rot_speed: f64) -> Self {
        Self {
            move_speed,
            rot_speed,
            view: Iso3::translation(0.0, 0.0, 100.0),
            change_flag: false,
        }
    }

    pub fn set_view(&self, camera: &mut three_d::Camera) {
        let origin = self.view * Point3::origin();
        let look = self.view * -Vector3::z_axis();
        let up = self.view * Vector3::y_axis();

        // The view is based on the internal isometry, in which -Z is the view direction, +X is the
        // right direction, and +Y is the up direction.
        camera.set_view(origin.to_cg(), (origin + look.into_inner()).to_cg(), up.to_cg());
    }

    fn roll(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 1800.0;
        let roll = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), rot);
        self.view = self.view * roll;
        self.change_flag = true;
    }

    fn pitch(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 1800.0;
        let pitch = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), rot);
        self.view = self.view * pitch;
        self.change_flag = true;
    }

    fn yaw(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 1800.0;
        let yaw = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), rot);
        self.view = self.view * yaw;
        self.change_flag = true;
    }

    fn move_up(&mut self, input: f32) {
        let move_dist = input as f64 * self.move_speed / 10.0;
        let translation = Iso3::translation(0.0, move_dist, 0.0);
        self.view = self.view * translation;
        self.change_flag = true;
    }

    fn move_right(&mut self, input: f32) {
        let move_dist = input as f64 * self.move_speed / 10.0;
        let translation = Iso3::translation(-move_dist, 0.0, 0.0);
        self.view = self.view * translation;
        self.change_flag = true;
    }

    fn move_forward(&mut self, input: f32) {
        let move_dist = input as f64 * self.move_speed / 10.0;
        let translation = Iso3::translation(0.0, 0.0, -move_dist);
        self.view = self.view * translation;
        self.change_flag = true;
    }

    /// Handles the events. Must be called each frame.
    pub fn handle_events(
        &mut self,
        camera: &mut three_d::Camera,
        events: &mut [three_d::Event],
    ) -> bool {
        self.change_flag = false;

        for event in events.iter_mut() {
            match event {
                three_d::Event::MouseMotion {
                    delta,
                    button,
                    modifiers,
                    handled,
                    ..
                } => {
                    if *handled {
                        continue;
                    }

                    match (button, mod_state(modifiers)) {
                        // With no modifiers, dragging the left mouse button yaws and pitches
                        (Some(MouseButton::Left), ModifierState::None) => {
                            self.yaw(delta.0);
                            self.pitch(delta.1);
                            *handled = true;
                        },
                        // With Shift, dragging the left mouse button rolls
                        (Some(MouseButton::Left), ModifierState::ShiftOnly) => {
                            self.roll(delta.0);
                            *handled = true;
                        },

                        // With no modifiers, dragging the right mouse button moves the camera
                        // left/right/up/down
                        (Some(MouseButton::Right), ModifierState::None) => {
                            self.move_right(delta.0);
                            self.move_up(delta.1);
                            *handled = true;
                        },

                        // With Shift, dragging the right mouse button moves the camera
                        // forward/backward/left/right
                        (Some(MouseButton::Right), ModifierState::ShiftOnly) => {
                            self.move_right(delta.0);
                            self.move_forward(delta.1);
                            *handled = true;
                        },

                        _ => {}
                    }


                }
                _ => {}
            }
        }

        if self.change_flag {
            self.set_view(camera);
        }

        self.change_flag
    }
}
