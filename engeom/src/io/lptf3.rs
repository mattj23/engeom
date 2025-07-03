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

use crate::common::triangulation::parallel_row2::{build_parallel_row_strip, StripRowPoint};
use crate::geom3::mesh::HalfEdgeMesh;
use crate::{Point3, PointCloud, Result};
use alum::Handle;
use rayon::prelude::*;
use std::io::{Read, Seek};
use std::path::Path;

use self::downsample::load_lptf3_downfilter;
pub use self::downsample::Lptf3DsParams;
pub use self::loader::Lptf3Loader;

#[derive(Debug, Clone, Copy)]
pub enum Lptf3Load {
    All,
    TakeEveryN(u32),
    SmoothSample(Lptf3DsParams),
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

pub fn load_lptf3_mesh_original(file_path: &Path, take_every: Option<u32>) -> Result<HalfEdgeMesh> {
    let strip_r = 3.0; // The maximum edge ratio for the strip triangulation.
    let world_r = 8.0; // The maximum edge ratio for world triangulation.

    let mut loader = Lptf3Loader::new(file_path, take_every, true)?;

    let mut last_delaunay_row: Option<(Vec<StripRowPoint>, f64)> = None;

    let mut mesh = HalfEdgeMesh::new();
    let max_spacing = take_every.unwrap_or(1) as f64 * loader.y_translation * 2.0;

    while let Some(full) = loader.get_next_frame_points()? {
        let mut row = Vec::new();

        for i in full.to_take.iter() {
            let ih = mesh
                .add_vertex(full.points[*i].at_y(full.y_pos).coords)
                .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
            row.push(StripRowPoint::new(full.points[*i].x, ih));
        }

        if let Some((last_row, last_y_pos)) = &last_delaunay_row {
            if (full.y_pos - *last_y_pos).abs() < max_spacing {
                let r = build_parallel_row_strip(last_row, *last_y_pos, &row, full.y_pos, strip_r)?;
                for (i0, i1, i2) in r {
                    // Check the edge ratio on actual points
                    let pa: Point3 = mesh
                        .point(i0)
                        .map_err(|e| format!("Failed to get point {}: {:?}", i0, e))?
                        .into();
                    let pb: Point3 = mesh
                        .point(i1)
                        .map_err(|e| format!("Failed to get point {}: {:?}", i1, e))?
                        .into();
                    let pc: Point3 = mesh
                        .point(i2)
                        .map_err(|e| format!("Failed to get point {}: {:?}", i2, e))?
                        .into();
                    let ea = (pa - pb).norm();
                    let eb = (pb - pc).norm();
                    let ec = (pc - pa).norm();

                    let edge_ratio = ea.max(eb).max(ec) / max_spacing;
                    if edge_ratio < world_r {
                        mesh.add_tri_face(i0, i1, i2)
                            .map_err(|e| format!("Failed to add face: {:?}", e))?;
                    }
                }
            }
        }

        last_delaunay_row = Some((row, full.y_pos));
    }

    Ok(mesh)
}

pub fn load_lptf3_mesh(file_path: &Path, params: Lptf3DsParams) -> Result<HalfEdgeMesh> {
    // let strip_r = 3.0; // The maximum edge ratio for the strip triangulation.
    // let world_r = 8.0; // The maximum edge ratio for world triangulation.
    //
    // let result = load_downsample_filter_lptf3(file_path, params)?;
    // let max_spacing = take_every as f64 * result.y_translation * 2.0;
    //
    // // let mut last_delaunay_row: Option<(Vec<StripRowPoint>, f64)> = None;
    //
    // // First build the mesh vertices and the corresponding rows of strip row points
    // let mut mesh = HalfEdgeMesh::new();
    // let mut strip_rows = Vec::new();
    // for row in result.rows.iter() {
    //     let mut strip_row = Vec::new();
    //     for p in row.iter() {
    //         let ih = mesh
    //             .add_vertex(p.coords)
    //             .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
    //         strip_row.push(StripRowPoint::new(p.x, ih));
    //     }
    //     strip_rows.push(strip_row);
    // }
    //
    // // Now iterate through the strip rows and build the mesh
    // for row_i in 0..strip_rows.len() - 1 {
    //     if result.rows[row_i].is_empty() || result.rows[row_i + 1].is_empty() {
    //         continue; // Skip empty rows
    //     }
    //
    //     let y0 = result.rows[row_i][0].y;
    //     let y1 = result.rows[row_i + 1][0].y;
    //
    //     // If the rows are too far apart, skip the triangulation
    //     if (y1 - y0).abs() > max_spacing {
    //         continue;
    //     }
    //
    //     let row0 = &strip_rows[row_i];
    //     let row1 = &strip_rows[row_i + 1];
    //
    //     // Build the strip triangulation between the two rows
    //     let r = build_parallel_row_strip(row0, y0, row1, y1, strip_r)?;
    //     for (i0, i1, i2) in r {
    //         // Check the edge ratio on actual points
    //         let pa: Point3 = mesh
    //             .point(i0)
    //             .map_err(|e| format!("Failed to get point {}: {:?}", i0, e))?
    //             .into();
    //         let pb: Point3 = mesh
    //             .point(i1)
    //             .map_err(|e| format!("Failed to get point {}: {:?}", i1, e))?
    //             .into();
    //         let pc: Point3 = mesh
    //             .point(i2)
    //             .map_err(|e| format!("Failed to get point {}: {:?}", i2, e))?
    //             .into();
    //         let ea = (pa - pb).norm();
    //         let eb = (pb - pc).norm();
    //         let ec = (pc - pa).norm();
    //
    //         let edge_ratio = ea.max(eb).max(ec) / max_spacing;
    //         if edge_ratio < world_r {
    //             mesh.add_tri_face(i0, i1, i2)
    //                 .map_err(|e| format!("Failed to add face: {:?}", e))?;
    //         }
    //     }
    // }
    //
    // Ok(mesh)
    todo!()
}

fn expand_colors(colors: &[u8]) -> Vec<[u8; 3]> {
    colors.iter().map(|&c| [c, c, c]).collect()
}
