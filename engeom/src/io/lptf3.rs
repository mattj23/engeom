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

use crate::airfoil::AirfoilGeometry;
use crate::common::points::{dist, mean_point};
use crate::common::triangulation::VertexBuilder;
use crate::common::triangulation::parallel_row2::{StripRowPoint, build_parallel_row_strip};
use crate::geom3::mesh::HalfEdgeMesh;
use crate::{
    Iso3, Mesh, Plane3, Point3, PointCloud, Result, SurfacePoint3, SvdBasis3, UnitVec3, Vector3,
};
use alum::Handle;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::Path;

pub fn load_lptf3_downfilter(file_path: &Path, take_every: u32) -> Result<PointCloud> {
    if take_every < 2 {
        return Err("take_every must be at least 2".into());
    }

    let mut loader = Lptf3Loader::new(file_path, Some(take_every), true)?;

    // Prepare the full point cloud into a set of rows
    // =========================================================================================
    // At the end of this step, `all_points` will contain a vector of vectors, where each inner
    // vector is the 3d point data for a single row of points, sorted in ascending order by x
    // coordinate.  The `row_data` will contain the indices of the points that are destined for the
    // final point cloud. The color vector will have been filled with the color values in the order
    // of the points matching with a flattening of `all_points`.
    let mut all_points = Vec::new();
    let mut all_colors = Vec::new();
    let mut row_data = Vec::new();

    let mut total_points = 0;

    while let Some(full) = loader.get_next_frame_points()? {
        let mut row = Vec::new();
        let mut c_row = Vec::new();
        for p in full.points.iter() {
            row.push(p.at_y(full.y_pos));
            c_row.push(p.color.unwrap_or(0));
        }
        all_points.push(row);
        all_colors.push(c_row);
        total_points += full.to_take.len();
        row_data.push(full.to_take)
    }

    // Sample the final point cloud
    // =========================================================================================
    // We will iterate through each row of points and for each index in the `row_data` element
    // we will find all points within a sampling distance of the point at that index.  We will
    // perform a gaussian weighted SVD on those points and then correct the working point's z value
    // to intersect with the plane.

    // The number of rows to look forward and backwards when sampling the point cloud.
    let look_rows = if take_every % 2 == 0 {
        take_every / 2
    } else {
        (take_every + 1) / 2
    } as i32;

    let look_dist = look_rows as f64 * loader.y_translation * 1.5;
    let mut final_points = Vec::new();
    let mut final_colors = Vec::new();

    for (row_i, to_take) in row_data.iter().enumerate() {
        for col_i in to_take.iter() {
            // The point that we're working on
            let p = all_points[row_i][*col_i];

            let mut samples = Vec::new();
            for check_i in (row_i as i32 - look_rows)..=(row_i as i32 + look_rows) {
                if check_i < 0 || check_i >= all_points.len() as i32 {
                    continue; // Skip rows that are out of bounds
                }
                let check_row = &all_points[check_i as usize];

                // Binary search for the first point that is p.x - look_dist
                let target = p.x - look_dist;
                let start = check_row
                    .binary_search_by(|a| a.x.total_cmp(&target))
                    .unwrap_or_else(|i| i);

                for col in start..check_row.len() {
                    let check_p = &check_row[col];
                    if (check_p.x - p.x).abs() <= look_dist {
                        samples.push(*check_p);
                    }
                    if check_p.x > p.x + look_dist {
                        // If the point is beyond the look distance, we can stop checking this row
                        break;
                    }
                }
            }

            final_points.push(adjust_point(&p, &samples, look_dist));
            final_colors.push(all_colors[row_i][*col_i]);
        }
    }

    let c = if loader.has_color {
        Some(expand_colors(&final_colors))
    } else {
        None
    };
    PointCloud::try_new(final_points, None, c)
}

fn gaussian_weight(x: f64, sigma: f64) -> f64 {
    (-0.5 * (x.powi(2) / sigma.powi(2))).exp()
}

fn adjust_point(p: &Point3, samples: &[Point3], look_dist: f64) -> Point3 {
    if samples.len() < 3 {
        return *p;
    }

    let mut weights = Vec::with_capacity(samples.len());
    for cp in samples.iter() {
        let d = dist(p, cp);
        weights.push(gaussian_weight(d, look_dist));
    }

    // Calculate the SVD basis from the samples
    let sp = SurfacePoint3::new(*p, UnitVec3::new_unchecked(Vector3::z()));
    let basis = SvdBasis3::from_points(&samples, Some(&weights));

    // Check that the standard deviation of the second basis is at least 1/3rd of the first so that
    // we didn't just best-fit a plane to a line
    let stdev = basis.basis_stdevs();
    if stdev[0] > stdev[1] * 3.0 {
        return *p;
    }

    let plane = Plane3::from(&basis);
    if let Some(t) = plane.intersection_distance(&sp) {
        sp.at_distance(t)
    } else {
        *p
    }
}

fn expand_colors(colors: &[u8]) -> Vec<[u8; 3]> {
    colors.iter().map(|&c| [c, c, c]).collect()
}

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
    let mut loader = Lptf3Loader::new(file_path, take_every, true)?;
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

pub fn load_lptf3_mesh(file_path: &Path, take_every: Option<u32>) -> Result<HalfEdgeMesh> {
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

    /// If set to true, the loader will return all frames and all points in those frames, regardless
    /// of what `take_every` is set to, however, the `to_take` field in the `FullFrame` will have
    /// no indices in it indicating that this is a skipped frame. This data is used for neighbor
    /// based smoothing during downsampling.
    return_all: bool,
}

impl Lptf3Loader {
    /// Creates a new instance of the Lptf3Loader.
    pub fn new(file_path: &Path, take_every: Option<u32>, return_all: bool) -> Result<Self> {
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
            Some(n) if n > 1 => Some(n),
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
            return_all,
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
        let mut header = FrameHeader::from_buffer(&buffer)?;

        // If the number of points is zero, we return Skip (no need to seek, there are no points)
        if header.num_points == 0 {
            // If the number of points is zero, we skip this frame
            return Ok(HdrRd::Skip);
        }

        // If this frame is being skipped, we seek to the next frame and return Skip
        if let Some(take_n) = self.take_every {
            if header.frame_index % take_n != 0 {
                if self.return_all {
                    header.skip = true;
                } else {
                    let skip_bytes = self.bytes_per_point * header.num_points;
                    self.file.seek_relative(skip_bytes as i64)?;
                    return Ok(HdrRd::Skip);
                }
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
    pub fn get_next_frame_points(&mut self) -> Result<Option<FullFrame>> {
        let mut points = Vec::new();

        let header = match self.seek_next_valid_frame()? {
            Some(h) => h,
            None => return Ok(None), // No more frames to read
        };

        let y_pos = header.frame_index as f64 * self.y_translation;
        let skip_int = self.skip_spacing.map(|s| (s / header.x_res) as i32);

        // We're going to start by reading all the points and then sorting them by x coordinate.
        // ========================================================================================
        let mut raw_points = Vec::with_capacity(header.num_points as usize);
        for _ in 0..header.num_points {
            let (x_raw, z_raw) = read_raw_point(&mut self.file, self.is_32_bit)?;
            let c = if self.has_color {
                read_u8(&mut self.file)?
            } else {
                0
            };
            raw_points.push(RawPoint {
                x: x_raw,
                z: z_raw,
                c,
            });
        }
        raw_points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());

        // Now we'll iterate through the raw points to mark the indices to keep
        // ========================================================================================
        // We have to calculate the skip offset based on the first point in order to
        // pick a value large enough to ensure that the skip index will never be less than
        // zero, otherwise it will produce a missing row when it crosses the zero boundary.
        let skip_offset = skip_int
            .map(|s| s * (-raw_points[0].x / s + 1))
            .unwrap_or(i32::MIN);
        let mut last_skip_index = i32::MIN;
        let mut to_take = Vec::new();

        for (i, raw) in raw_points.into_iter().enumerate() {
            if let Some(skip_i) = skip_int {
                let skip_index = (raw.x + skip_offset) / skip_i;
                if skip_index > last_skip_index {
                    last_skip_index = skip_index;
                    if !header.skip {
                        to_take.push(i);
                    }
                }
            } else {
                // If we're not skipping, we take every point
                to_take.push(i);
            }

            let p = FramePoint {
                x: (raw.x as f64) * header.x_res + header.x_offset,
                z: (raw.z as f64) * header.z_res + header.z_offset,
                color: if self.has_color { Some(raw.c) } else { None },
            };

            points.push(p);
        }

        // Sort points by x coordinate
        let result = FullFrame::new(header, points, y_pos, to_take);

        Ok(Some(result))
    }
}

enum HdrRd {
    Valid(FrameHeader),
    Skip,
    EndOfFile,
}

pub struct FullFrame {
    pub header: FrameHeader,
    pub points: Vec<FramePoint>,
    pub y_pos: f64,
    pub to_take: Vec<usize>,
}

impl FullFrame {
    pub fn new(
        header: FrameHeader,
        points: Vec<FramePoint>,
        y_pos: f64,
        take_indices: Vec<usize>,
    ) -> Self {
        Self {
            header,
            points,
            y_pos,
            to_take: take_indices,
        }
    }
}

struct RawPoint {
    x: i32,
    z: i32,
    c: u8,
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

#[derive(Clone)]
struct FrameHeader {
    frame_index: u32,
    num_points: u32,
    x_offset: f64,
    z_offset: f64,
    x_res: f64,
    z_res: f64,

    /// This flag is set during the loading process to indicate that this frame consists entirely
    /// of points that would be skipped based on the `take_every` parameter.  If the loader is set
    /// to return all frames, this flag is what distinguishes a skipped frame from a valid frame
    pub skip: bool,
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
            skip: false,
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
