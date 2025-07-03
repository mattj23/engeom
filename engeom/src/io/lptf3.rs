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

use crate::common::triangulation::VertexBuilder;
use crate::common::triangulation::parallel_row2::{StripRowPoint, build_parallel_row_strip};
use crate::{Mesh, Point3, PointCloud, Result};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::Path;

/// Read a lptf3 (Laser Profile Triangulation Format 3D) file and return a `PointCloud`.
///
/// This function reads a LPTF3 file, which is a compact file format for storing 3D point data
/// taken from a laser profile triangulation scanner. The format is simple and compact, capable
/// of practically storing about 200k points (with an 8-bit color value each) per MB when using a
/// 16-bit coordinate format, or half that when using a 32-bit coordinate format.
///
/// # Arguments
///
/// * `file_path`: the path to the LPTF3 file to read.
/// * `take_every`: skips frames that are not multiples of this value, while also skipping points
///   that are not roughly spaced by the same gap between frames. To take only every nth frame,
///   set this to `Some(n)`. If `None`, or `Some(1)` all frames are taken.
///
/// returns: Result<PointCloud, Box<dyn Error, Global>>
pub fn load_lptf3(file_path: &Path, take_every: Option<u32>) -> Result<PointCloud> {
    let mut loader = Lptf3Loader::new(file_path, take_every)?;
    let mut points = Vec::new();
    let mut colors = Vec::new();

    while let Some((y_pos, frame_points)) = loader.get_next_frame_points()? {
        for p in frame_points {
            points.push(p.at_y(y_pos));

            if let Some(color) = p.color {
                colors.push([color; 3]);
            }
        }
    }

    let c = if loader.has_color { Some(colors) } else { None };

    PointCloud::try_new(points, None, c)
}

pub fn load_lptf3_mesh(file_path: &Path, take_every: Option<u32>) -> Result<Mesh> {
    let strip_r = 3.0; // The maximum edge ratio for the strip triangulation.
    let world_r = 8.0; // The maximum edge ratio for world triangulation.

    let mut loader = Lptf3Loader::new(file_path, take_every)?;

    let mut last_delaunay_row: Option<(Vec<StripRowPoint>, f64)> = None;

    let mut vertices = VertexBuilder::new();
    let max_spacing = take_every.unwrap_or(1) as f64 * loader.y_translation * 2.0;

    let mut faces = Vec::new();
    while let Some((y_pos, frame_points)) = loader.get_next_frame_points()? {
        let mut row = Vec::new();
        for p in frame_points {
            let i = vertices.push(p.at_y(y_pos));
            row.push(StripRowPoint::new(p.x, i));
        }

        if let Some((last_row, last_y_pos)) = &last_delaunay_row {
            if (y_pos - *last_y_pos).abs() < max_spacing {
                let r = build_parallel_row_strip(last_row, *last_y_pos, &row, y_pos, strip_r)?;
                for f in r {
                    // Check the edge ratio on actual points
                    let pa = &vertices.points()[f[0]];
                    let pb = &vertices.points()[f[1]];
                    let pc = &vertices.points()[f[2]];
                    let ea = (pa - pb).norm();
                    let eb = (pb - pc).norm();
                    let ec = (pc - pa).norm();

                    let edge_ratio = ea.max(eb).max(ec) / max_spacing;
                    if edge_ratio < world_r {
                        faces.push([f[0] as u32, f[1] as u32, f[2] as u32]);
                    }
                }
            }
        }

        last_delaunay_row = Some((row, y_pos));
    }

    if faces.is_empty() {
        return Err("No valid faces found in the LPTF3 file".into());
    }

    let mesh = Mesh::new(vertices.take_points(), faces, false);

    Ok(mesh)
}

/// This struct offers frame-by-frame loading of LPTF3 files, which allows the generalized loading
/// process to be shared between different loading strategies. For instance, a LPTF3 file can be
/// loaded naively, returning a point cloud.  However, loading it while knowing the sensor geometry
/// can allow for normal direction and point quality estimation, or even color adjustment scalars
/// based on reflection.
///
/// This struct contains the core logic for reading the files and is intended to be used by
/// different loading mechanisms.
pub struct Lptf3Loader {
    file: BufReader<File>,
    take_every: Option<u32>,
    bytes_per_point: u32,
    y_translation: f64,
    skip_spacing: Option<f64>,
    has_color: bool,
    is_32_bit: bool,
}

impl Lptf3Loader {
    /// Creates a new instance of the Lptf3Loader.
    pub fn new(file_path: &Path, take_every: Option<u32>) -> Result<Self> {
        let path_str = file_path
            .to_str()
            .ok_or_else(|| format!("Invalid path: {}", file_path.display()))?;

        let raw_file = File::open(file_path)
            .map_err(|e| format!("Failed to open file '{}': {}", path_str, e))?;
        let mut f = BufReader::new(raw_file);

        // Read the magic number
        let mut magic = [0; 5];
        f.read_exact(&mut magic)?;
        if &magic != b"LPTF3" {
            return Err(format!("Invalid magic number in file '{}'", path_str).into());
        }

        // Read the version number
        let version = read_u16(&mut f)?;
        if version != 1 {
            return Err(format!("Unsupported version {} in file '{}'", version, path_str).into());
        }

        // Read the data flags
        let data_flags = read_u16(&mut f)?;
        let is_32_bit = (data_flags & 0x0001) != 0;
        let has_color = (data_flags & 0x0002) != 0;

        // Read the motion type
        let motion_type = read_u8(&mut f)?;
        if motion_type != 0 {
            return Err(format!(
                "Unsupported motion type {} in file '{}'",
                motion_type, path_str
            )
            .into());
        }

        // Read the y translation and skip distance for motion type 0
        let y_translation = (read_u32(&mut f)? as f64) / 1_000_000.0; // Convert from nanometers to mm
        let skip_spacing = take_every.map(|t| t as f64 * y_translation);

        // Prepare the point and color vectors
        // Calculate the number of bytes per point
        let bytes_per_point = if is_32_bit { 8 } else { 4 } + if has_color { 1 } else { 0 };

        let take = match take_every {
            Some(1) => None,
            Some(n) => Some(n),
            _ => None,
        };

        Ok(Self {
            file: f,
            take_every: take,
            bytes_per_point: bytes_per_point as u32,
            y_translation,
            skip_spacing,
            has_color,
            is_32_bit,
        })
    }

    /// This function reads the frame header at the current file position and does one of the
    /// following actions:
    ///
    /// - If it can't read the frame header, it returns `HdrRd::EndOfFile`.
    /// - If there's a reason to not take the data from this frame, it returns `HdrRd::Skip` and
    ///   seeks forward to the position of the next frame header.
    /// - If the frame header is valid and the number of points is greater than zero, it returns
    ///   `HdrRd::Valid(header)` with the parsed frame header. It leaves the file cursor at the
    ///   position of the first point in the frame.
    /// - If it encounters an error (other than the end of the file), it returns the error.
    fn read_next_frame_header(&mut self) -> Result<HdrRd> {
        let mut buffer = [0; 24];
        let read_result = self.file.read_exact(&mut buffer);
        if read_result.is_err() {
            return Ok(HdrRd::EndOfFile);
        }

        // If we can read the frame header, parse it
        let header = FrameHeader::from_buffer(&buffer)?;

        // If the number of points is zero, we return Skip (no need to seek, there are no points)
        if header.num_points == 0 {
            // If the number of points is zero, we skip this frame
            return Ok(HdrRd::Skip);
        }

        // If this frame is being skipped, we seek to the next frame and return Skip
        if let Some(take_n) = self.take_every {
            if header.frame_index % take_n != 0 {
                let skip_bytes = self.bytes_per_point * header.num_points;
                self.file.seek_relative(skip_bytes as i64)?;
                return Ok(HdrRd::Skip);
            }
        }

        Ok(HdrRd::Valid(header))
    }

    /// Seek the file cursor forward to the next valid frame header. This function will do one of
    /// three things:
    /// - If it encounters an error, it returns `Err`.
    /// - If it reaches the end of the file, it returns `Ok(None)`.
    /// - If it finds a valid frame header, it returns `Ok(Some(header))` with the parsed frame
    ///   header.
    fn seek_next_valid_frame(&mut self) -> Result<Option<FrameHeader>> {
        loop {
            match self.read_next_frame_header()? {
                HdrRd::Valid(header) => return Ok(Some(header)),
                HdrRd::Skip => continue, // Skip this frame and read the next one
                HdrRd::EndOfFile => return Ok(None), // Reached the end of the file
            }
        }
    }

    /// This function reads forward in the file looking for a frame to read.  It will skip over
    /// frames that have zero points, or would be skipped by the `take_every` parameter.  Under
    /// normal circumstances, it will either return an Ok(None) if it has reached the end of the
    /// file, or an Ok(Some(f64, Vec<FramePoint>)) if it has found a valid frame to read.
    ///
    /// In the return from a valid frame, the f64 value is the y position of the frame, and the
    /// vector of `FramePoint` structs contains the x/z/color values for each point in the frame.
    /// The points are sorted by their x coordinate.
    ///
    /// If an error occurs during the operation, it will return an `Err` with the error message.
    pub fn get_next_frame_points(&mut self) -> Result<Option<(f64, Vec<FramePoint>)>> {
        let mut points = Vec::new();

        let header = match self.seek_next_valid_frame()? {
            Some(h) => h,
            None => return Ok(None), // No more frames to read
        };

        let y_pos = header.frame_index as f64 * self.y_translation;
        let skip_int = self.skip_spacing.map(|s| (s / header.x_res) as i32);
        let mut last_skip_index = i32::MIN;
        let mut skip_offset = i32::MIN;

        for _ in 0..header.num_points {
            let (x_raw, z_raw) = read_raw_point(&mut self.file, self.is_32_bit)?;

            if let Some(skip_i) = skip_int {
                // We have to calculate the skip offset based on the first point in order to
                // pick a value large enough to ensure that the skip index will never be less than
                // zero, otherwise it will produce a missing row when it crosses the zero boundary.
                if skip_offset == i32::MIN {
                    skip_offset = skip_i * ((-x_raw / skip_i) + 1);
                }
            }

            let c = if self.has_color {
                Some(read_u8(&mut self.file)?)
            } else {
                None
            };

            if let Some(skip_i) = skip_int {
                let skip_index = (x_raw + skip_offset) / skip_i;
                if skip_index <= last_skip_index {
                    continue;
                }
                last_skip_index = skip_index;
            }

            let p = FramePoint {
                x: (x_raw as f64) * header.x_res + header.x_offset,
                z: (z_raw as f64) * header.z_res + header.z_offset,
                color: c,
            };

            points.push(p);
        }

        // Sort points by x coordinate
        points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal));

        Ok(Some((y_pos, points)))
    }
}

enum HdrRd {
    Valid(FrameHeader),
    Skip,
    EndOfFile,
}

pub struct FramePoint {
    pub x: f64,
    pub z: f64,
    pub color: Option<u8>,
}

impl FramePoint {
    pub fn new(x: f64, z: f64, color: Option<u8>) -> Self {
        Self { x, z, color }
    }

    pub fn at_zero(&self) -> Point3 {
        self.at_y(0.0)
    }

    pub fn at_y(&self, y: f64) -> Point3 {
        Point3::new(self.x, y, self.z)
    }
}

struct FrameHeader {
    frame_index: u32,
    num_points: u32,
    x_offset: f64,
    z_offset: f64,
    x_res: f64,
    z_res: f64,
}

impl FrameHeader {
    fn from_buffer(buffer: &[u8; 24]) -> Result<Self> {
        if buffer.len() != 24 {
            return Err("Invalid frame header size".into());
        }

        let frame_index = u32::from_le_bytes(buffer[0..4].try_into()?);
        let num_points = u32::from_le_bytes(buffer[4..8].try_into()?);
        let x_offset = read_offset(&buffer[8..12])?;
        let z_offset = read_offset(&buffer[12..16])?;
        let x_res = read_res(&buffer[16..20])?;
        let z_res = read_res(&buffer[20..24])?;

        Ok(Self {
            frame_index,
            num_points,
            x_offset,
            z_offset,
            x_res,
            z_res,
        })
    }
}

fn read_raw_point<R: Read>(reader: &mut R, is_32bit: bool) -> Result<(i32, i32)> {
    let (x, z) = if is_32bit {
        (read_i32(reader)?, read_i32(reader)?)
    } else {
        (read_i16(reader)? as i32, read_i16(reader)? as i32)
    };
    Ok((x, z))
}

fn read_res(buffer: &[u8]) -> Result<f64> {
    // Convert from nanometers to millimeters
    Ok(u32::from_le_bytes(buffer[0..4].try_into()?) as f64 / 1_000_000.0)
}

fn read_offset(buffer: &[u8]) -> Result<f64> {
    // Convert from micrometers to millimeters
    Ok(i32::from_le_bytes(buffer[0..4].try_into()?) as f64 / 1_000.0)
}

fn read_u16<R: Read>(reader: &mut R) -> Result<u16> {
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_i32<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_i16<R: Read>(reader: &mut R) -> Result<i16> {
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(i16::from_le_bytes(buf))
}
