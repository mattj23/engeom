//! This module contains functionality for working with a compact file format for storing 3D point
//! data taken from a laser profile triangulation scanner. The format is simple and designed to be
//! similar to the intermediate representation produced by sensors in their raw output.
//!
//! The format is identified by the extension `.lptf3` and the header structured as follows:
//!
//! - bytes 0-5: magic number b"LPTF3" to identify the file type
//! - bytes 6-7: version number (currently 1)
//! - byte 8-9: data flags
//!   - bit 0: Bytes per point coordinate (0=16 bit, 1=32 bit)
//!   - bit 1: Color data present (0=none, 1=single uint8)
//! - byte 10: motion type
//!   - 0: fixed y translation
//!   - 1-255: not implemented
//!
//! After the header, the next set of bytes will depend on the motion type:
//!
//! - If the motion type is 0, the next four bytes will be a 32 bit uint representing the y
//!   translation of the scanner per frame in nanometers.
//!
//! Following the motion type values, the file will contain a repeating sequence consisting of
//! a frame header and a variable number of point entries. The frame header consists of the
//! following 24 bytes:
//!
//! - bytes 0-3: frame number (uint32)
//! - bytes 4-7: number of points in the frame (uint32)
//! - bytes 8-11: x offset for all frame points in micrometers (int32)
//! - bytes 12-15: z offset for all frame points in micrometers (int32)
//! - bytes 16-19: x resolution for all frame points in nanometers (uint32)
//! - bytes 20-23: z resolution for all frame points in nanometers (uint32)
//!
//! Following the frame header, there will be the number of individual point entries specified
//! by the frame header. Each point entry consists of the following:
//!
//! - x coordinate (16 or 32-bit signed integer, depending on the data flags)
//! - z coordinate (16 or 32-bit signed integer, depending on the data flags)
//! - color (optional, 8-bit unsigned integer if color data is present)
//!
//! At the end of the point entries, there will be either another frame header or the end of the
//! file.

mod downsample;
mod loader;
mod mesh;
mod uncertainty;

use self::downsample::load_lptf3_downfilter;
use crate::geom3::mesh::HalfEdgeMesh;
use crate::io::lptf3::downsample::load_downsample_filter_lptf3;
use crate::io::lptf3::mesh::load_lptf3_mesh_core;
use crate::{Point3, PointCloud, Result};
use std::path::Path;

pub use self::downsample::Lptf3DsParams;
pub use self::loader::Lptf3Loader;
pub use self::uncertainty::*;

/// The Lptf3Load enum defines the different ways to load data from a LPTF3 file, and is used to
/// pass these options to loading functions.
#[derive(Debug, Clone, Copy)]
pub enum Lptf3Load {
    All,
    TakeEveryN(u32),
    SmoothSample(Lptf3DsParams),
}

/// This trait defines the interface for uncertainty models that can be used when loading LPTF3
/// data from a file.
pub trait Lptf3UncertaintyModel {
    fn value(&self, x: f64, z: f64) -> f64;
}

/// Read a lptf3 (Laser Profile Triangulation Format 3D) file and return a `PointCloud`.
///
/// This function reads a LPTF3 file, which is a compact file format for storing 3D point data
/// taken from a laser profile triangulation scanner. The format is simple and compact, capable
/// of practically storing about 200k points (with an 8-bit color value each) per MB when using a
/// 16-bit coordinate format, or half that when using a 32-bit coordinate format.
///
/// There are a few different ways to load the data, controlled by the `Lptf3Load` enum:
///   - `Lptf3Load::All`: Load all points from the file.
///   - `Lptf3Load::TakeEveryN(n)`: Load every Nth row from the file. The loader will attempt to
///     roughly match the x spacing of the points to the gap distance between rows, resulting in a
///     grid-like point cloud with an approximately uniform point spacing when viewed from the
///     X-Y plane.  This is a very fast method of retrieving a downsampled point cloud.
///   - `Lptf3Load::SmoothSample(params)`: Load the points using a downsampling filter, which
///     downsamples the point cloud similar to the `TakeEveryN` method, but also performs a gaussian
///     smoothing step using the full original cloud.  This takes the longest time, but can remove
///     a significant amount of noise from the data by making use of an adjacency structure that
///     will be lost once the points are turned into a cloud.
///
/// # Arguments
///
/// * `file_path`: A path to the LPTF3 file to load.
/// * `load`: An enum specifying how to load the data from the file.
///
/// returns: Result<PointCloud, Box<dyn Error, Global>>
pub fn load_lptf3(file_path: &Path, load: Lptf3Load) -> Result<PointCloud> {
    match load {
        Lptf3Load::All => load_take_every(file_path, None),
        Lptf3Load::TakeEveryN(n) => load_take_every(file_path, Some(n)),
        Lptf3Load::SmoothSample(params) => load_lptf3_downfilter(file_path, params),
    }
}

pub fn load_lptf3_mesh_uncertainty(
    file_path: &Path,
    load: Lptf3Load,
    uncertainty_model: &dyn Lptf3UncertaintyModel,
) -> Result<(HalfEdgeMesh, Vec<f64>)> {
    let (mesh, uncert) = load_lptf3_mesh_core(file_path, load, Some(uncertainty_model))?;
    let u = uncert.ok_or("Mesh loader did not return uncertainty values")?;
    Ok((mesh, u))
}

/// Read a lptf3 (Laser Profile Triangulation Format 3D) file and return a `HalfEdgeMesh`.
///
/// This function reads a LPTF3 file, which is a compact file format for storing 3D point data
/// taken from a laser profile triangulation scanner. The format is simple and compact, capable
/// of practically storing about 200k points (with an 8-bit color value each) per MB when using a
/// 16-bit coordinate format, or half that when using a 32-bit coordinate format.
///
/// There are a few different ways to load the data, controlled by the `Lptf3Load` enum:
///   - `Lptf3Load::All`: Load all points from the file.
///   - `Lptf3Load::TakeEveryN(n)`: Load every Nth row from the file. The loader will attempt to
///     roughly match the x spacing of the points to the gap distance between rows, resulting in a
///     grid-like point cloud with an approximately uniform point spacing when viewed from the
///     X-Y plane.  This is a very fast method of retrieving a downsampled point cloud.
///   - `Lptf3Load::SmoothSample(params)`: Load the points using a downsampling filter, which
///     downsamples the point cloud similar to the `TakeEveryN` method, but also performs a gaussian
///     smoothing step using the full original cloud.  This takes the longest time, but can remove
///     a significant amount of noise from the data by making use of an adjacency structure that
///     will be lost once the points are turned into a cloud.
///
/// Once the points are loaded, they will be converted into a triangle mesh by connecting points in
/// adjacent rows with triangles that meet certain edge length criterial. The result is a fast mesh
/// that can be built using knowledge of the LPTF3's internal structure rather than having to rely
/// on more general techniques that can build meshes from arbitrary point clouds.
///
/// # Arguments
///
/// * `file_path`: A path to the LPTF3 file to load.
/// * `load`: An enum specifying how to load the data from the file.
///
/// returns: Result<PolyMeshT<3, NaAdaptor>, Box<dyn Error, Global>>
pub fn load_lptf3_mesh(file_path: &Path, load: Lptf3Load) -> Result<HalfEdgeMesh> {
    let (result, _) = load_lptf3_mesh_core(file_path, load, None)?;
    Ok(result)
}

fn load_take_every(file_path: &Path, take_every: Option<u32>) -> Result<PointCloud> {
    let mut loader = Lptf3Loader::new(file_path, take_every, false)?;
    let mut points = Vec::new();
    let mut colors = Vec::new();

    while let Some(full) = loader.get_next_frame_points()? {
        if full.header.skip {
            // If this frame is skipped, we don't add any points to the point cloud.
            continue;
        }

        for i in full.to_take.iter() {
            points.push(full.points[*i].at_y(full.y_pos));

            if let Some(color) = full.points[*i].color {
                colors.push([color; 3]);
            }
        }
    }

    let c = if loader.has_color { Some(colors) } else { None };

    PointCloud::try_new(points, None, c)
}

fn get_loader_point_rows(
    file_path: &Path,
    take_every: Option<u32>,
) -> Result<(u32, f64, Vec<Vec<Point3>>)> {
    let mut loader = Lptf3Loader::new(file_path, take_every, true)?;
    let mut point_rows = Vec::new();

    while let Some(full) = loader.get_next_frame_points()? {
        let mut row = Vec::new();

        for i in full.to_take.iter() {
            let p = full.points[*i].at_y(full.y_pos);
            row.push(p);
        }
        if !row.is_empty() {
            point_rows.push(row);
        }
    }

    Ok((take_every.unwrap_or(1), loader.y_translation, point_rows))
}

fn get_downfilter_point_rows(
    file_path: &Path,
    params: Lptf3DsParams,
) -> Result<(u32, f64, Vec<Vec<Point3>>)> {
    let result = load_downsample_filter_lptf3(file_path, params)?;
    Ok((params.take_every, result.y_translation, result.rows))
}

fn expand_colors(colors: &[u8]) -> Vec<[u8; 3]> {
    colors.iter().map(|&c| [c, c, c]).collect()
}
