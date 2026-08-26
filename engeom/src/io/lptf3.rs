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

mod comprehensive;
mod downsample;
mod loader;
mod mesh;
mod uncertainty;

use crate::geom3::mesh::MeshData3;
use crate::io::lptf3::downsample::load_downsample_filter_lptf3;
use crate::io::lptf3::mesh::load_lptf3_mesh_data_core;
use crate::{Point3, PointCloud3, Result};
use std::path::Path;

pub use self::comprehensive::*;
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

/// Quickly get the point distribution for the LPTF3 file, which returns a vector of point counts
/// for bins of the given size. For certain applications this is a quick way to heuristically check
/// if a scan of a known scene roughly matches the expected data.
///
/// # Arguments
///
/// * `file_path`:
/// * `bin_size`:
///
/// returns: Result<Vec<usize, Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn lptf3_point_distribution(file_path: &Path, bin_size: f64) -> Result<Vec<usize>> {
    let mut distribution = Vec::new();
    let mut loader = Lptf3Loader::new(file_path, None, true)?;
    while let Some(full) = loader.get_next_frame_points()? {
        if full.y_pos < 0.0 {
            // This shouldn't happen
            continue;
        }

        let index = (full.y_pos / bin_size).floor() as usize;
        if index >= distribution.len() {
            distribution.resize(index + 1, 0);
        }

        distribution[index] += full.points.len();
    }

    Ok(distribution)
}

/// Read a lptf3 (Laser Profile Triangulation Format 3D) file and return a `PointCloud3`.
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
/// # Attributes
///
/// If the file has a color channel it is expanded to grayscale and stored as `point_colors`. The
/// channel is a single 8-bit laser return intensity, not a true color.
///
/// # Arguments
///
/// * `file_path`: A path to the LPTF3 file to load.
/// * `load`: An enum specifying how to load the data from the file.
///
/// returns: `Result<PointCloud3>`
pub fn load_lptf3(file_path: &Path, load: Lptf3Load) -> Result<PointCloud3> {
    let (points, colors) = match load {
        Lptf3Load::All => load_take_every(file_path, None)?,
        Lptf3Load::TakeEveryN(n) => load_take_every(file_path, Some(n))?,
        Lptf3Load::SmoothSample(params) => load_downfilter_points(file_path, params)?,
    };

    let mut data = PointCloud3::new(points);
    if colors.is_some() {
        data.set_point_colors(colors)?;
    }

    Ok(data)
}

/// Read a lptf3 (Laser Profile Triangulation Format 3D) file and return a `MeshData3`.
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
/// adjacent rows with triangles that meet certain edge length criteria. The result is a fast mesh
/// that can be built using knowledge of the LPTF3's internal structure rather than having to rely
/// on more general techniques that can build meshes from arbitrary point clouds.
///
/// # Orphan Points
///
/// **Points which end up in no face are discarded.** The edge criteria reject the triangles around
/// a dropout or a depth discontinuity, and the first and last rows of a scan may have nothing to
/// join to, so a real scan always leaves some points unconnected. They are real measurements, but
/// they are not part of a surface, and what is being asked for here is a mesh.
///
/// The consequence is that the point buffer of the returned mesh is a *subset* of the points that
/// `load_lptf3` returns for the same file and load mode, in the same order but with gaps. Use
/// `load_lptf3` instead when every measured point matters. A scan which produces no faces at all,
/// such as a single row, returns an empty mesh rather than a bag of loose points.
///
/// # Attributes
///
/// Data which the sensor recorded alongside the geometry is carried onto the mesh rather than
/// discarded:
///
/// - If the file has a color channel, it is expanded to grayscale and stored as `point_colors`.
///   The channel is a single 8-bit laser return intensity, not a true color.
/// - If an uncertainty model is supplied, it is evaluated at every point and the results are
///   stored as `point_stdev`. A model returning a negative or non-finite value is an error, since
///   that is not a standard deviation.
///
/// Note that the uncertainty is evaluated at the *final* position of each point, so under
/// `Lptf3Load::SmoothSample` the model sees the smoothed coordinates. The model is also a function
/// of the sensor geometry alone, so it describes the inherent triangulation uncertainty of the
/// measurement and knows nothing about how much noise the smoothing removed.
///
/// # Arguments
///
/// * `file_path`: A path to the LPTF3 file to load.
/// * `load`: An enum specifying how to load the data from the file.
/// * `uncertainty`: An optional model of the sensor's measurement uncertainty. When supplied, the
///   standard deviation it predicts at each point is stored on the mesh as `point_stdev`.
///
/// returns: Result<MeshData3, Box<dyn Error, Global>>
pub fn load_lptf3_mesh_data(
    file_path: &Path,
    load: Lptf3Load,
    uncertainty: Option<&dyn Lptf3UncertaintyModel>,
) -> Result<MeshData3> {
    load_lptf3_mesh_data_core(file_path, load, uncertainty)
}

/// Load every `take_every`th row of a file, returning the points and, if the file has a color
/// channel, one grayscale triple per point.
type LoadedPoints = (Vec<Point3>, Option<Vec<[u8; 3]>>);

fn load_take_every(file_path: &Path, take_every: Option<u32>) -> Result<LoadedPoints> {
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

    Ok((points, c))
}

/// The same as `load_take_every` but running the gaussian downsampling filter, which flattens the
/// row structure the filter works over once it is done with it.
fn load_downfilter_points(file_path: &Path, params: Lptf3DsParams) -> Result<LoadedPoints> {
    let downsampled = load_downsample_filter_lptf3(file_path, params)?;
    let points = downsampled.rows.into_iter().flatten().collect::<Vec<_>>();

    let colors = downsampled.colors.map(|c| {
        let flat = c.into_iter().flatten().collect::<Vec<_>>();
        expand_colors(&flat)
    });

    Ok((points, colors))
}

/// The rows of points recovered from a LPTF3 file, in the order they were scanned. Each inner
/// vector is one profile, sorted in ascending order by X coordinate.
///
/// When `colors` is `Some`, it has exactly the same row and column structure as `points`, so a
/// flattening of the two stays aligned element for element.
struct Lptf3Rows {
    take_every: u32,
    y_translation: f64,
    points: Vec<Vec<Point3>>,
    colors: Option<Vec<Vec<u8>>>,
}

fn get_loader_point_rows(file_path: &Path, take_every: Option<u32>) -> Result<Lptf3Rows> {
    let mut loader = Lptf3Loader::new(file_path, take_every, true)?;
    let mut point_rows = Vec::new();
    let mut color_rows = Vec::new();

    while let Some(full) = loader.get_next_frame_points()? {
        let mut row = Vec::new();
        let mut c_row = Vec::new();

        for i in full.to_take.iter() {
            let p = full.points[*i].at_y(full.y_pos);
            row.push(p);
            c_row.push(full.points[*i].color.unwrap_or(0));
        }
        if !row.is_empty() {
            point_rows.push(row);
            color_rows.push(c_row);
        }
    }

    let colors = if loader.has_color {
        Some(color_rows)
    } else {
        None
    };

    Ok(Lptf3Rows {
        take_every: take_every.unwrap_or(1),
        y_translation: loader.y_translation,
        points: point_rows,
        colors,
    })
}

fn get_downfilter_point_rows(file_path: &Path, params: Lptf3DsParams) -> Result<Lptf3Rows> {
    let result = load_downsample_filter_lptf3(file_path, params)?;
    Ok(Lptf3Rows {
        take_every: params.take_every,
        y_translation: result.y_translation,
        points: result.rows,
        colors: result.colors,
    })
}

fn expand_colors(colors: &[u8]) -> Vec<[u8; 3]> {
    colors.iter().map(|&c| [c, c, c]).collect()
}
