//! This module has I/O functionality for points clouds

use crate::{PointCloud, PointCloudFeatures, Result, UnitVec3, Vector3};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;
use std::ptr::write;

fn write_to_f32(writer: &mut impl Write, value: f64) -> Result<()> {
    let bytes = (value as f32).to_le_bytes();
    writer.write_all(&bytes)?;
    Ok(())
}

fn read_from_f32(reader: &mut impl Read) -> Result<f64> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    let value = f32::from_le_bytes(bytes);
    Ok(value as f64)
}

pub fn write_bxyz(path: &Path, cloud: &PointCloud) -> Result<()> {
    // The fields are X, Y, Z, NX, NY, NZ, R, G, B, STDEV
    let mut header_bytes = vec![
        cloud.normals().is_some(),
        cloud.colors().is_some(),
        cloud.std_devs().is_some(),
    ];

    let point_count = cloud.points().len() as u32;

    let file = File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);

    for byte in header_bytes {
        writer.write_all(&[byte as u8])?;
    }
    writer.write_all(&point_count.to_le_bytes())?;

    for (i, point) in cloud.points().iter().enumerate() {
        write_to_f32(&mut writer, point.x)?;
        write_to_f32(&mut writer, point.y)?;
        write_to_f32(&mut writer, point.z)?;
        // writer.write_all(&point.x.to_le_bytes())?;
        // writer.write_all(&point.y.to_le_bytes())?;
        // writer.write_all(&point.z.to_le_bytes())?;

        if let Some(normals) = cloud.normals() {
            let normal = &normals[i];
            write_to_f32(&mut writer, normal.x)?;
            write_to_f32(&mut writer, normal.y)?;
            write_to_f32(&mut writer, normal.z)?;
            // writer.write_all(&(normal.x as f32).to_le_bytes())?;
            // writer.write_all(&(normal.y as f32).to_le_bytes())?;
            // writer.write_all(&(normal.z as f32).to_le_bytes())?;
        }

        if let Some(colors) = cloud.colors() {
            let color = &colors[i];
            writer.write_all(&color[0].to_le_bytes())?;
            writer.write_all(&color[1].to_le_bytes())?;
            writer.write_all(&color[2].to_le_bytes())?;
        }

        if let Some(std_devs) = cloud.std_devs() {
            write_to_f32(&mut writer, std_devs[i])?;
        }
    }

    Ok(())
}

pub fn load_bxyz(path: &Path) -> Result<PointCloud> {
    let file = File::open(path)?;
    let mut reader = std::io::BufReader::new(file);
    let mut header = [0u8; 3];
    reader.read_exact(&mut header)?;
    let has_normals = header[0] != 0;
    let has_colors = header[1] != 0;
    let has_std_devs = header[2] != 0;

    let mut count_bytes = [0u8; 4];
    reader.read_exact(&mut count_bytes)?;
    let point_count = u32::from_le_bytes(count_bytes) as usize;

    let mut points = Vec::with_capacity(point_count);
    let mut normals = if has_normals {
        Some(Vec::with_capacity(point_count))
    } else {
        None
    };
    let mut colors = if has_colors {
        Some(Vec::with_capacity(point_count))
    } else {
        None
    };
    let mut std_devs = if has_std_devs {
        Some(Vec::with_capacity(point_count))
    } else {
        None
    };

    for _ in 0..point_count {
        let x = read_from_f32(&mut reader)?;
        let y = read_from_f32(&mut reader)?;
        let z = read_from_f32(&mut reader)?;
        // let mut point_bytes = [0u8; 24];
        // reader.read_exact(&mut point_bytes)?;
        // let x = f64::from_le_bytes(point_bytes[0..8].try_into().unwrap());
        // let y = f64::from_le_bytes(point_bytes[8..16].try_into().unwrap());
        // let z = f64::from_le_bytes(point_bytes[16..24].try_into().unwrap());
        points.push(crate::Point3::new(x, y, z));

        if has_normals {
            let nx = read_from_f32(&mut reader)?;
            let ny = read_from_f32(&mut reader)?;
            let nz = read_from_f32(&mut reader)?;

            // let mut normal_bytes = [0u8; 12];
            // reader.read_exact(&mut normal_bytes)?;
            // let nx = f32::from_le_bytes(normal_bytes[0..4].try_into().unwrap());
            // let ny = f32::from_le_bytes(normal_bytes[4..8].try_into().unwrap());
            // let nz = f32::from_le_bytes(normal_bytes[8..12].try_into().unwrap());
            if let Some(normals_vec) = &mut normals {
                normals_vec.push(UnitVec3::new_normalize(Vector3::new(nx, ny, nz)));
            }
        }

        if has_colors {
            let mut color_bytes = [0u8; 3];
            reader.read_exact(&mut color_bytes)?;
            if let Some(colors_vec) = &mut colors {
                colors_vec.push([color_bytes[0], color_bytes[1], color_bytes[2]]);
            }
        }

        if has_std_devs {
            let std_dev = read_from_f32(&mut reader)?;
            // let mut std_dev_bytes = [0u8; 4];
            // reader.read_exact(&mut std_dev_bytes)?;
            // let std_dev = f32::from_le_bytes(std_dev_bytes) as f64;
            if let Some(std_devs_vec) = &mut std_devs {
                std_devs_vec.push(std_dev);
            }
        }
    }

    PointCloud::try_new(points, normals, colors, std_devs)
}
