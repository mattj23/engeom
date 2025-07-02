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

use crate::{Point3, PointCloud, Result};
use std::io::Read;
use std::path::Path;

pub fn load_lptf3(file_path: &Path, take_every: Option<u32>) -> Result<PointCloud> {
    let path_str = file_path
        .to_str()
        .ok_or_else(|| format!("Invalid path: {}", file_path.display()))?;

    let raw_file = std::fs::File::open(file_path)
        .map_err(|e| format!("Failed to open file '{}': {}", path_str, e))?;
    let mut f = std::io::BufReader::new(raw_file);

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
    let mut points = Vec::new();
    let mut colors = Vec::new();

    // Calculate the number of bytes per point
    let bytes_per_point = if is_32_bit { 8 } else { 4 } + if has_color { 1 } else { 0 };

    // Read the frames
    // Frame header size is 4 (frame index) + 4 (number of points) + 4 (x offset) +
    let mut frame_header = [0; 24];
    while let Ok(()) = f.read_exact(&mut frame_header) {
        let frame_index = u32::from_le_bytes(frame_header[0..4].try_into().unwrap());
        let num_points = u32::from_le_bytes(frame_header[4..8].try_into().unwrap());

        if let Some(take_n) = take_every {
            if frame_index % take_n != 0 {
                f.seek_relative(bytes_per_point * num_points as i64)?;
                continue;
            }
        }

        let x_offset = read_offset(&frame_header[8..12])?;
        let z_offset = read_offset(&frame_header[12..16])?;
        let x_res = read_res(&frame_header[16..20])?;
        let z_res = read_res(&frame_header[20..24])?;
        let y_pos = y_translation * (frame_index as f64);

        let skip_int = skip_spacing.map(|s| (s / x_res) as i32);
        let mut last_skip_index = i32::MIN;
        let mut skip_offset = i32::MIN;

        for _ in 0..num_points {
            let (x_raw, z_raw) = read_raw_point(&mut f, is_32_bit)?;

            if let Some(skip_i) = skip_int {
                // We have to calculate the skip offset based on the first point in order to
                // pick a value large enough to ensure that the skip index will never be less than
                // zero, otherwise it will produce a missing row when it crosses the zero boundary.
                if skip_offset == i32::MIN {
                    skip_offset = skip_i * ((-x_raw / skip_i) + 1);
                }
            }

            let c = if has_color {
                Some(read_u8(&mut f)?)
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

            points.push(Point3::new(
                (x_raw as f64) * x_res + x_offset,
                y_pos,
                (z_raw as f64) * z_res + z_offset,
            ));

            if let Some(color) = c {
                colors.push([color; 3]);
            }
        }
    }

    PointCloud::try_new(points, None, Some(colors))
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
