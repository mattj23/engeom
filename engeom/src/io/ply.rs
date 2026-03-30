#![cfg(feature = "ply")]

use std::fs::File;
use std::io::BufReader;
use crate::{Mesh, Point3, Result};
use std::path::Path;
use serde::{Deserialize, Serialize};
use serde_ply::from_reader;

/// Load a triangle mesh from a ply file. The ply file must have x, y, and z coordinates for each
/// vertex and a list of triangles.
///
/// # Arguments
///
/// * `path`: Path to the ply file.
///
/// returns: Result<Mesh, Box<dyn Error, Global>>
pub fn load_ply_mesh(path: &Path) -> Result<Mesh> {

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let ply_data: PlyData = from_reader(&mut reader)?;
    let vertices = ply_data.vertex.iter()
        .map(|v| Point3::new(v.x, v.y, v.z))
    .collect::<Vec<_>>();
    let faces = ply_data.face.iter()
        .map(|f| [f.vertex_indices[0], f.vertex_indices[1], f.vertex_indices[2]])
    .collect::<Vec<_>>();

    Ok(Mesh::new(vertices, faces, false))
}

#[derive(Serialize, Deserialize, Debug)]
struct PlyVertex {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Serialize, Deserialize, Debug)]
struct PlyFace {
    vertex_indices: Vec<u32>,
}

#[derive(Serialize, Deserialize, Debug)]
struct PlyData {
    vertex: Vec<PlyVertex>,
    face: Vec<PlyFace>,
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::get_test_file_path;
    use approx::assert_relative_eq;

    #[test]
    fn load_ply_mesh_test() -> Result<()> {
        let path = get_test_file_path("bun_zipper_res4.ply");
        let mesh = load_ply_mesh(&path)?;

        assert_eq!(mesh.faces().len(), 948);

        assert_relative_eq!(mesh.vertices()[0], Point3::new(-0.0312216,0.126304, 0.00514924), epsilon = 1e-6);
        assert_relative_eq!(mesh.vertices()[300], Point3::new(0.0414633,0.0456301,-0.00354513), epsilon = 1e-6);
        assert_relative_eq!(mesh.vertices()[400], Point3::new(-0.0422796,0.0699944,-0.0161463), epsilon = 1e-6);
        Ok(())
    }
}