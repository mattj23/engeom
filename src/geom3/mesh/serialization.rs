use crate::geom2::Point2;
use crate::geom3::mesh::{Mesh, UvMapping};
use crate::geom3::Point3;
use crate::utility::unflatten_points;
use crate::Result;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::error::Error;

/// A serializable, editable representation of a mesh.
#[derive(Deserialize, Serialize, Clone)]
pub struct MeshFlatData {
    pub vertices: Vec<f64>,
    pub triangles: Vec<u32>,
}

impl MeshFlatData {
    pub fn new(vertices: Vec<f64>, triangles: Vec<u32>) -> Self {
        Self {
            vertices,
            triangles,
        }
    }
}

impl TryFrom<MeshFlatData> for Mesh {
    type Error = Box<dyn Error>;

    fn try_from(value: MeshFlatData) -> Result<Self> {
        let vertices = unflatten_points::<3>(&value.vertices)?;
        let triangles = value
            .triangles
            .iter()
            .tuples()
            .map(|(a, b, c)| [*a, *b, *c])
            .collect();
        Ok(Mesh::new(vertices, triangles, false))
    }
}

/// A serializable, editable representation of a mesh.
#[derive(Deserialize, Serialize, Clone)]
pub struct MeshData {
    pub vertices: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,

    #[serde(default)]
    pub uv: Option<Vec<Point2>>,

    #[serde(default)]
    pub is_solid: bool,
}

impl MeshData {
    pub fn new(
        vertices: Vec<Point3>,
        triangles: Vec<[u32; 3]>,
        uv: Option<Vec<Point2>>,
        is_solid: bool,
    ) -> Self {
        Self {
            vertices,
            triangles,
            uv,
            is_solid,
        }
    }
}

impl From<MeshData> for Mesh {
    fn from(value: MeshData) -> Self {
        let uv_map = value
            .uv
            .map(|uv| UvMapping::new_from_vertices(&uv).expect("Failed to create UV mapping"));
        Mesh::new_with_uv(value.vertices, value.triangles, value.is_solid, uv_map)
    }
}

impl From<&MeshData> for Mesh {
    fn from(value: &MeshData) -> Self {
        let uv_map = value
            .uv
            .as_ref()
            .map(|uv| UvMapping::new_from_vertices(uv).expect("Failed to create UV mapping"));

        Mesh::new_with_uv(
            value.vertices.clone(),
            value.triangles.clone(),
            value.is_solid,
            uv_map,
        )
    }
}

impl From<&Mesh> for MeshData {
    fn from(value: &Mesh) -> Self {
        let uv = value.uv().map(|uv| uv.to_vertices());
        let tri_mesh = value.tri_mesh();
        Self::new(
            tri_mesh.vertices().iter().map(|v| *v).collect_vec(),
            tri_mesh.indices().iter().map(|t| *t).collect_vec(),
            uv,
            value.is_solid(),
        )
    }
}
