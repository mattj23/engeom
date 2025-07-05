use crate::Mesh;
use crate::geom3::mesh::HalfEdgeMesh;
use alum::{Handle, HasIterators, HasTopology};
use three_d::CpuMesh;

pub trait ToCpuMesh {
    fn to_cpu_mesh(&self) -> CpuMesh;
}

impl ToCpuMesh for Mesh {
    fn to_cpu_mesh(&self) -> CpuMesh {
        let points = self
            .vertices()
            .iter()
            .map(|v| three_d::vec3(v.x, v.y, v.z))
            .collect::<Vec<_>>();

        let indices = self
            .faces()
            .iter()
            .flat_map(|x| x.iter().map(|v| *v))
            .collect();

        let mut cm = CpuMesh {
            positions: three_d::Positions::F64(points),
            indices: three_d::Indices::U32(indices),
            ..Default::default()
        };

        cm.compute_normals();
        cm
    }
}

impl ToCpuMesh for HalfEdgeMesh {
    /// Generates a `CpuMesh` from the `HalfEdgeMesh`.
    fn to_cpu_mesh(&self) -> CpuMesh {
        let point_prop = self.points();
        let points = point_prop.try_borrow().expect("Cannot borrow points");
        let f_status_prop = self.face_status_prop();
        let f_status = f_status_prop.try_borrow().unwrap();

        let mut cm = CpuMesh {
            positions: three_d::Positions::F64(
                points
                    .iter()
                    .map(|p| three_d::vec3(p.x, p.y, p.z))
                    .collect::<Vec<_>>(),
            ),
            indices: three_d::Indices::U32(
                self.triangulated_vertices(&f_status)
                    .flatten()
                    .map(|v| v.index())
                    .collect(),
            ),
            ..Default::default()
        };

        // I've found that `alum`'s normal computation sometimes produces strange results that
        // render with artifacts. The `compute_normals` method of `CpuMesh` seems to produce
        // better results, so we use that instead.
        cm.compute_normals();
        cm
    }
}
