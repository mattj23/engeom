#![cfg(feature = "ply")]

use crate::{Mesh, Point3, Result};
use serde::{Deserialize, Serialize};
use serde_ply::from_reader;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

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
    let vertices = ply_data
        .vertex
        .iter()
        .map(|v| Point3::new(v.x, v.y, v.z))
        .collect::<Vec<_>>();
    let faces = ply_data
        .face
        .iter()
        .map(|f| {
            [
                f.vertex_indices[0],
                f.vertex_indices[1],
                f.vertex_indices[2],
            ]
        })
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
    use crate::tests::{get_test_file_path, stanford_bun_4};
    use approx::assert_relative_eq;

    #[test]
    fn load_ply_mesh_test() -> Result<()> {
        let path = get_test_file_path("bun_zipper_res4.ply");
        let mesh = load_ply_mesh(&path)?;
        let expected = stanford_bun_4();

        assert_eq!(mesh.faces().len(), expected.faces().len());
        assert_eq!(mesh.vertices().len(), expected.vertices().len());
        for (a, b) in mesh.faces().iter().zip(expected.faces().iter()) {
            assert_eq!(a, b);
        }

        for (a, b) in mesh.vertices().iter().zip(expected.vertices().iter()) {
            assert_relative_eq!(a, b, epsilon = 0.000002);
        }

        Ok(())
    }
}
