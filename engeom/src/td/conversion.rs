use crate::geom3::mesh::HalfEdgeMesh;
use alum::{Handle, HasIterators, HasTopology};
use three_d::CpuMesh;
use crate::Mesh;


pub trait ToCpuMesh {
    fn to_cpu_mesh(&self) -> CpuMesh;
}

impl ToCpuMesh for Mesh {
    fn to_cpu_mesh(&self) -> CpuMesh {
        let points = self.vertices()
            .iter()
            .map(|v| three_d::vec3(v.x, v.y, v.z))
            .collect::<Vec<_>>();

        let vtx_normals = self.get_vertex_normals()
            .iter()
            .map(|v| three_d::vec3(v.x as f32, v.y as f32, v.z as f32))
            .collect::<Vec<_>>();

        let indices = self.faces().iter()
            .flat_map(|x| {
                x.iter().map(|v| *v)
            })
            .collect();

        CpuMesh {
            positions: three_d::Positions::F64(points),
            indices: three_d::Indices::U32(indices),
            normals: Some(vtx_normals),
            ..Default::default()
        }
    }
}

impl ToCpuMesh for HalfEdgeMesh {
    /// Generates a `CpuMesh` from the `HalfEdgeMesh`.  Be sure to call one of the
    /// `update_vertex_normals_*` methods before calling this method to ensure the
    /// vertex normals are computed correctly.
    fn to_cpu_mesh(&self) -> CpuMesh {
        let point_prop = self.points();
        let points = point_prop.try_borrow().expect("Cannot borrow points");

        let vtx_normal_prop = self.vertex_normals().unwrap();
        let vtx_normals = vtx_normal_prop.try_borrow().unwrap();
        let vtx_normals = vtx_normals
            .iter()
            .map(|n| three_d::vec3(n.x as f32, n.y as f32, n.z as f32))
            .collect::<Vec<_>>();

        let f_status_prop = self.face_status_prop();
        let f_status = f_status_prop.try_borrow().unwrap();

        CpuMesh {
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
            normals: Some(vtx_normals),
            ..Default::default()
        }
    }
}
