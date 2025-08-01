use crate::td::{CameraControl, ModState, ToCgVec3, ToEngeom3, mod_state, ToCpuMesh, cpu_mat};
use crate::{Iso3, Point3, Result};
use itertools::Itertools;
use std::collections::HashMap;
use three_d::{
    AmbientLight, Angle, AxisAlignedBoundingBox, Camera, ClearState, ColorMaterial, Context,
    CoreError, CpuMaterial, CpuMesh, DirectionalLight, Event, FrameOutput, Gm, Material, Mesh,
    MouseButton, Object, PhysicalMaterial, Srgba, Window, WindowSettings, degrees, pick, vec3,
};

pub struct SimpleViewer {
    items: HashMap<usize, Box<dyn Object>>,
    next_id: usize,
    window: Window,
    context: Context,
    pub shadows_on: bool,
}

impl SimpleViewer {
    /// Creates a new simple viewing window for inspecting 3D objects.  View controls are a simple
    /// free orbiting camera with no gimbal lock.
    ///
    /// Has very simple view controls:
    ///   - Right-click on an object to center the camera on the clicked point
    ///   - Left mouse drag to orbit around the center point
    ///   - Mouse wheel to translate the camera forward/backwards
    ///   - Shift + left mouse drag to roll the camera left/right
    ///   - Press 'S' to toggle shadows on/off
    ///   - Alt + left mouse drag to move the light direction
    ///
    /// # Arguments
    ///
    /// * `title`: A title for the window
    /// * `max_size`: An optional maximum size for the window. If `None`, the window will be
    ///   created with the default size.
    ///
    /// returns: Result<SimpleViewer, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use engeom::{Point3, Mesh};
    /// use engeom::td::{SimpleViewer, ToCpuMesh, cpu_mat};
    ///
    /// let p0 = Point3::new(5.0, 0.0, 0.0);
    /// let m0 = Mesh::create_capsule(&(-p0), &p0, 1.0, 100, 100);
    /// let mut view = SimpleViewer::new("Demo", None).unwrap();
    /// view.add_mesh(m0.to_cpu_mesh(), cpu_mat(150, 150, 150, 255, 0.7, 0.8));
    /// view.display().unwrap();
    /// ```
    pub fn new(title: &str, max_size: Option<(u32, u32)>) -> Result<SimpleViewer> {
        let window = Window::new(WindowSettings {
            title: title.to_string(),
            max_size,
            ..Default::default()
        })?;

        let context = window.gl();

        Ok(SimpleViewer {
            items: HashMap::new(),
            next_id: 0,
            window,
            context,
            shadows_on: true,
        })
    }

    pub fn context(&self) -> &Context {
        &self.context
    }

    pub fn window(&self) -> &Window {
        &self.window
    }

    pub fn add_object(&mut self, object: Box<dyn Object>) -> usize {
        let id = self.next_id;
        self.items.insert(id, object);
        self.next_id += 1;
        id
    }

    pub fn add_mesh(&mut self, mesh: CpuMesh, cpu_material: CpuMaterial) -> usize {
        let mat = if cpu_material.albedo.a < 255 {
            PhysicalMaterial::new_transparent(&self.context, &cpu_material)
        } else {
            PhysicalMaterial::new_opaque(&self.context, &cpu_material)
        };

        let item = Box::new(Gm::new(Mesh::new(&self.context, &mesh), mat)) as Box<dyn Object>;
        let id = self.next_id;
        self.items.insert(id, item);
        self.next_id += 1;
        id
    }

    pub fn add_polyline(
        &mut self,
        points: &[Point3],
        color: (u8, u8, u8),
        thickness: f64,
    ) -> Result<()> {
        if points.len() < 2 {
            return Ok(()); // Not enough points to create a polyline
        }

        let mut mesh = crate::Mesh::create_cylinder_between(&points[0], &points[1], thickness, 6);
        for i in 2..points.len() {
            let next_mesh =
                crate::Mesh::create_cylinder_between(&points[i - 1], &points[i], thickness, 6);
            mesh.append(&next_mesh)?
        }

        let cpu_mesh = mesh.to_cpu_mesh();
        self.add_mesh(cpu_mesh, cpu_mat(color.0, color.1, color.2, 255, 1.0, 0.0));


        Ok(())
    }

    pub fn remove_item(&mut self, id: usize) {
        self.items.remove(&id);
    }

    pub fn get_item(&self, id: usize) -> Option<&Box<dyn Object>> {
        self.items.get(&id)
    }

    pub fn items(&self) -> &HashMap<usize, Box<dyn Object>> {
        &self.items
    }

    pub fn display(self) -> Result<()> {
        let camera_angle: f64 = 45.0;

        // Find the combined bounding box of all items to center the camera
        let mut bounds = AxisAlignedBoundingBox::EMPTY;
        for item in self.items.values() {
            bounds.expand_with_aabb(item.aabb());
        }

        let e = bounds.size().to_engeom().norm();
        let t = (camera_angle / 2.0).to_radians().tan();
        let d = e / (2.0 * t);
        let center = bounds.center().to_engeom();

        let mut camera = Camera::new_perspective(
            self.window.viewport(),
            vec3(-500.0, 250.0, 200.0), // Position of the camera
            vec3(0.0, 0.0, 0.0),        // Target point the camera is looking at
            vec3(0.0, 0.0, 1.0),        // Up vector of the camera
            degrees(camera_angle as f32),
            0.1,
            10000.0,
        );

        let mut shadows_on = self.shadows_on;
        let mut control = CameraControl::new(1.0, 1.0, Iso3::from(center), d);
        control.set_view(&mut camera);

        let ambient = AmbientLight::new(&self.context, 0.7, Srgba::WHITE);
        let mut light0 =
            DirectionalLight::new(&self.context, 2.0, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
        let mut light1 =
            DirectionalLight::new(&self.context, 2.0, Srgba::WHITE, vec3(1.0, 1.0, 1.0));

        self.window.render_loop(move |mut frame_input| {
            let mut redraw = frame_input.first_frame;
            control.reset_change_flag();

            redraw |= camera.set_viewport(frame_input.viewport);

            // Handle general display viewer events
            for event in frame_input.events.iter() {
                match event {
                    Event::KeyPress {
                        kind,
                        modifiers,
                        handled,
                        ..
                    } => {
                        if *handled {
                            continue;
                        }
                        match (kind, mod_state(modifiers)) {
                            (three_d::Key::S, ModState::None) => {
                                // Toggle shadows on/off
                                shadows_on = !shadows_on;
                                if !shadows_on {
                                    light0.clear_shadow_map()
                                }

                                redraw = true;
                            }
                            _ => {}
                        }
                    }
                    Event::MousePress {
                        button,
                        position,
                        modifiers,
                        handled,
                        ..
                    } => {
                        if *handled {
                            continue;
                        }
                        match (button, mod_state(modifiers)) {
                            (MouseButton::Right, ModState::None) => {
                                for item in self.items.values() {
                                    if let Some(pick) = pick(
                                        &self.context,
                                        &camera,
                                        *position,
                                        std::iter::once(item.as_ref()),
                                    ) {
                                        control.set_center(Point3::from(pick.position.to_engeom()));
                                    }
                                }
                            }

                            _ => {}
                        }
                    }
                    _ => {}
                }
            }

            // Handle camera control events
            redraw |= control.handle_events(&mut camera, &mut frame_input.events);

            if redraw {
                light0.direction = control.world_light_vector().to_cg();

                if shadows_on {
                    for item in self.items.values() {
                        light0.generate_shadow_map(2048, std::iter::once(item.as_ref()));
                    }
                }

                frame_input
                    .screen()
                    .clear(ClearState::color_and_depth(0.1, 0.1, 0.1, 1.0, 1.0))
                    .write(|| {
                        for (_, item) in self.items.iter() {
                            item.render(&camera, &[&ambient, &light0]);
                        }

                        Ok::<(), CoreError>(())
                    })
                    .unwrap();
            }
            FrameOutput {
                swap_buffers: redraw,
                ..Default::default()
            }
        });

        Ok(())
    }
}
