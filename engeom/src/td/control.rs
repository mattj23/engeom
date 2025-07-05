use crate::na::{Translation, UnitQuaternion};
use crate::td::{ModState, ToCgVec3, mod_state};
use crate::{Iso3, Point3, UnitVec3, Vector3};
use three_d::{Event, MouseButton};

#[derive(Clone, Copy, Debug)]
/// A controller that performs orbiting and zooming around a movable point in 3D space while
/// avoiding gimbal lock by continuously moving the "up" direction.
///
/// The view center point can be updated externally (usually by the user selecting a point in the
/// scene), which will cause the camera to stay in its current position but update the direction
/// that it's looking.  The camera can be moved around the view center point by dragging the mouse
/// to perform yaw and pitch, by shift+dragging to roll, and by using the mouse wheel to zoom.
///
/// The camera's position in-out from the view center is controlled by the mouse wheel. The wheel
/// does not change the zoom of the camera, but instead translates the camera towards or away from
/// the view center point. Movement is exponential, so camera motion will slow as it approaches the
/// view center and accelerate as it moves away.
///
/// Additionally, there is a light vector that can be rotated around the camera using ALT+mouse
/// dragging. This light vector can be used by a scene to update a directional light source so
/// that it stays in relation to the camera, which is useful for scientific applications.
///
/// Internally, the view is represented as an isometry _at the view center point_ with the Y axis
/// pointing up, the Z axis pointing towards the camera, and the X axis pointing to the right. The
/// camera position is `distance` units away from the view center point in +Z.
pub struct CameraControl {
    /// A scaling factor for moving in and out
    move_speed: f64,
    rot_speed: f64,
    view: Iso3,
    distance: f64,
    light_vector: Vector3,
    change_flag: bool,
}

impl CameraControl {
    /// Creates a new camera control.
    pub fn new(move_speed: f64, rot_speed: f64, view: Iso3, distance: f64) -> Self {
        Self {
            move_speed,
            rot_speed,
            view,
            distance,
            light_vector: Vector3::new(-1.0, -0.25, -0.5).normalize(),
            change_flag: false,
        }
    }

    pub fn world_light_vector(&self) -> Vector3 {
        // Return the light vector in world space.
        self.view * self.light_vector
    }

    /// Return the current position of the camera in world space
    pub fn position(&self) -> Point3 {
        self.view * Point3::new(0.0, 0.0, self.distance)
    }

    pub fn changed(&self) -> bool {
        // Return true if the view has changed since the last call to reset_change_flag.
        self.change_flag
    }

    /// Return the view center of the camera in world space
    pub fn view_center(&self) -> Point3 {
        self.view * Point3::origin()
    }

    /// Return the up direction of the camera in world space
    pub fn up(&self) -> UnitVec3 {
        self.view * Vector3::y_axis()
    }

    pub fn right(&self) -> UnitVec3 {
        self.view * Vector3::x_axis()
    }

    /// Return the look direction of the camera in world space
    pub fn look(&self) -> UnitVec3 {
        self.view * -Vector3::z_axis()
    }

    pub fn reset_change_flag(&mut self) {
        // Reset the change flag to false, indicating that the view has not yet changed.
        self.change_flag = false;
    }

    pub fn set_center(&mut self, center: Point3) {
        // Set the view center to the given point, keeping the current position and orientation.
        let translation = center - self.view_center();
        self.view = Translation::from(translation) * self.view;
        self.change_flag = true;
    }

    pub fn set_view(&self, camera: &mut three_d::Camera) {
        // The view is based on the internal isometry, in which -Z is the view direction, +X is the
        // right direction, and +Y is the up direction.
        camera.set_view(
            self.position().to_cg(),
            self.view_center().to_cg(),
            self.up().to_cg(),
        );
    }

    fn roll(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 180.0;
        let roll = UnitQuaternion::from_axis_angle(&Vector3::z_axis(), rot);
        self.view = self.view * roll;
        self.change_flag = true;
    }

    fn pitch(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 180.0;
        let pitch = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), rot);
        self.view = self.view * pitch;
        self.change_flag = true;
    }

    fn yaw(&mut self, input: f32) {
        let rot = input as f64 * self.rot_speed / 180.0;
        let yaw = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), rot);
        self.view = self.view * yaw;
        self.change_flag = true;
    }

    fn pitch_light(&mut self, input: f32) {
        // Pitch the light vector
        let rot = input as f64 * self.rot_speed / 180.0;
        let pitch = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), rot);
        let updated = pitch * self.light_vector;
        if Vector3::z_axis().dot(&updated) < 0.0 {
            self.light_vector = pitch * self.light_vector;
            self.change_flag = true;
        }
    }

    fn yaw_light(&mut self, input: f32) {
        // Yaw the light vector
        let rot = input as f64 * self.rot_speed / 180.0;
        let yaw = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), rot);
        let updated = yaw * self.light_vector;
        if Vector3::z_axis().dot(&updated) < 0.0 {
            self.light_vector = yaw * self.light_vector;
            self.change_flag = true;
        }
    }

    /// Handles the events. Must be called each frame.
    pub fn handle_events(&mut self, camera: &mut three_d::Camera, events: &mut [Event]) -> bool {
        for event in events.iter_mut() {
            match event {
                Event::MouseMotion {
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
                        (Some(MouseButton::Left), ModState::None) => {
                            self.yaw(-delta.0);
                            self.pitch(-delta.1);
                            *handled = true;
                        }
                        // With Shift, dragging the left mouse button rolls
                        (Some(MouseButton::Left), ModState::ShiftOnly) => {
                            self.roll(delta.0);
                            *handled = true;
                        }
                        // With alt, dragging the left mouse button yaws and pitches the light
                        (Some(MouseButton::Left), ModState::AltOnly) => {
                            self.yaw_light(delta.0);
                            self.pitch_light(delta.1);
                            *handled = true;
                        }
                        _ => {}
                    }
                }
                Event::MouseWheel {
                    delta,
                    modifiers,
                    handled,
                    ..
                } => {
                    if *handled {
                        continue;
                    }

                    match mod_state(modifiers) {
                        // With no modifiers, the mouse wheel zooms in and out
                        ModState::None => {
                            self.distance *= 1.01_f64.powf(delta.1 as f64 * self.move_speed);
                            self.change_flag = true;
                            *handled = true;
                        }
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
