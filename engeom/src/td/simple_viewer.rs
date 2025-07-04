use crate::Result;
use crate::td::CameraControl;
use std::collections::HashMap;
use three_d::{AmbientLight, Camera, ClearState, Context, CoreError, CpuMaterial, CpuMesh, DirectionalLight, Event, FrameOutput, Geometry, MouseButton, Object, Srgba, Window, WindowSettings, degrees, pick, vec3, PhysicalMaterial, Gm, Mesh};

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

    pub fn add_mesh(
        &mut self,
        mesh: CpuMesh,
        cpu_material: CpuMaterial,
    ) -> usize {
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

        let mut control = CameraControl::new(1.0, 2.0);
        control.set_view(&mut camera);

        let ambient = AmbientLight::new(&self.context, 0.7, Srgba::WHITE);
        let mut directional0 =
            DirectionalLight::new(&self.context, 2.0, Srgba::WHITE, vec3(-1.0, -1.0, -1.0));
        let mut directional1 =
            DirectionalLight::new(&self.context, 2.0, Srgba::WHITE, vec3(1.0, 1.0, 1.0));

        self.window.render_loop(move |mut frame_input| {
            let mut redraw = frame_input.first_frame;
            redraw |= camera.set_viewport(frame_input.viewport);
            redraw |= control.handle_events(&mut camera, &mut frame_input.events);

            for event in frame_input.events.iter() {
                if let Event::MousePress {
                    button, position, ..
                } = *event
                {
                    if button == MouseButton::Left {
                        // Pick
                        for item in self.items.values() {
                            if let Some(pick) =
                                pick(&self.context, &camera, position, std::iter::once(item.as_ref()))
                            {
                                println!("{:?}", pick);
                            }
                        }
                    }
                }
            }

            if redraw {
                // directional0.generate_shadow_map(1024, &view_mesh);
                // directional1.generate_shadow_map(1024, &view_mesh);

                frame_input
                    .screen()
                    .clear(ClearState::color_and_depth(0.1, 0.1, 0.1, 1.0, 1.0))
                    .write(|| {
                        for (_, item) in self.items.iter() {
                            item.render(&camera, &[&ambient, &directional0, &directional1]);
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
