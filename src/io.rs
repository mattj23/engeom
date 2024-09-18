use crate::Result;
use parry3d_f64::na::{Point3, Vector3};
use serde::Serialize;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

#[cfg(feature = "stl")]
use stl_io;

#[cfg(feature = "stl")]
use crate::geom3::Mesh;

#[cfg(feature = "stl")]
pub fn read_mesh_stl(path: &Path) -> Result<Mesh> {
    let mut file = OpenOptions::new().read(true).open(path)?;
    let mesh = stl_io::read_stl(&mut file)?;

    let vertices = mesh
        .vertices
        .iter()
        .map(|v| Point3::new(v[0] as f64, v[1] as f64, v[2] as f64))
        .collect::<Vec<_>>();

    let triangles = mesh
        .faces
        .iter()
        .map(|f| {
            [
                f.vertices[0] as u32,
                f.vertices[1] as u32,
                f.vertices[2] as u32,
            ]
        })
        .collect::<Vec<_>>();

    Ok(Mesh::new(vertices, triangles, false))
}

#[cfg(feature = "stl")]
pub fn write_mesh_stl(path: &Path, mesh: &Mesh) -> Result<()> {
    let mut faces = Vec::new();
    for triangle in mesh.tri_mesh().triangles() {
        if let Some(normal) = triangle.normal() {
            let vertices = [
                stl_io::Vertex::new([
                    triangle.a.x as f32,
                    triangle.a.y as f32,
                    triangle.a.z as f32,
                ]),
                stl_io::Vertex::new([
                    triangle.b.x as f32,
                    triangle.b.y as f32,
                    triangle.b.z as f32,
                ]),
                stl_io::Vertex::new([
                    triangle.c.x as f32,
                    triangle.c.y as f32,
                    triangle.c.z as f32,
                ]),
            ];
            faces.push(stl_io::Triangle {
                normal: stl_io::Normal::new([normal.x as f32, normal.y as f32, normal.z as f32]),
                vertices,
            });
        }
    }

    if path.exists() {
        std::fs::remove_file(path)?;
    }

    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;

    stl_io::write_stl(&mut file, faces.iter())?;
    Ok(())
}

pub fn write_xyz(path: &Path, points: &[Point3<f64>]) -> Result<()> {
    if path.exists() {
        std::fs::remove_file(path)?;
    }

    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;

    for point in points {
        writeln!(&mut file, "{} {} {}", point.x, point.y, point.z)?;
    }

    Ok(())
}

pub fn write_xyzn(path: &Path, points: &[Point3<f64>], normals: &[Vector3<f64>]) -> Result<()> {
    if path.exists() {
        std::fs::remove_file(path)?;
    }

    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut buffered = BufWriter::new(file);

    for (point, normal) in points.iter().zip(normals.iter()) {
        writeln!(
            &mut buffered,
            "{} {} {} {} {} {}",
            point.x, point.y, point.z, normal.x, normal.y, normal.z
        )?;
    }

    Ok(())
}

pub fn json_elements_save<T>(path: &Path, item: &T) -> Result<()>
where
    T: Serialize + ?Sized,
{
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let bytes = serde_json::to_vec_pretty(item)?;
    writer.write_all(&bytes)?;
    Ok(())
}
