use crate::td::{CameraControl, ModState, ToEngeom3, mod_state, ToCgVec3};
use crate::{Iso3, Point3, Result};
use std::collections::HashMap;
use itertools::Itertools;
use three_d::{
    AmbientLight, Camera, ClearState, Context, CoreError, CpuMaterial, CpuMesh, DirectionalLight,
    Event, FrameOutput, Geometry, Gm, Mesh, MouseButton, Object, PhysicalMaterial, Srgba, Window,
    WindowSettings, degrees, pick, vec3,
};

pub struct SimpleViewer {
    items: HashMap<usize, Box<dyn Object>>,
    next_id: usize,
    window: Window,
    context: Context,
}

impl SimpleViewer {
    pub fn new(max_size: Option<(u32, u32)>, title: String) -> Result<SimpleViewer> {
        let window = Window::new(WindowSettings {
            title,
            max_size,
            ..Default::default()
        })?;

        let context = window.gl();

        Ok(SimpleViewer {
            items: HashMap::new(),
            next_id: 0,
            window,
            context,
        })
    }

    pub fn context(&self) -> &Context {
        &self.context
    }

    pub fn window(&self) -> &Window {
        &self.window
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
        let mut camera = Camera::new_perspective(
            self.window.viewport(),
            vec3(-500.0, 250.0, 200.0), // Position of the camera
            vec3(0.0, 0.0, 0.0),        // Target point the camera is looking at
            vec3(0.0, 0.0, 1.0),        // Up vector of the camera
            degrees(45.0),
            0.1,
            10000.0,
        );

        let mut shadows_on = true;

        let mut control = CameraControl::new(1.0, 1.0, Iso3::identity(), 500.0);
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
                        kind, modifiers, handled, ..
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
