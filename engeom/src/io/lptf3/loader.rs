use crate::Point3;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

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
    pub take_every: Option<u32>,
    pub bytes_per_point: u32,
    pub y_translation: f64,
    pub skip_spacing: Option<f64>,
    pub has_color: bool,
    pub is_32_bit: bool,

    /// If set to true, the loader will return all frames and all points in those frames, regardless
    /// of what `take_every` is set to, however, the `to_take` field in the `FullFrame` will have
    /// no indices in it indicating that this is a skipped frame. This data is used for neighbor
    /// based smoothing during downsampling.
    pub return_all: bool,
}

impl Lptf3Loader {
    /// Creates a new instance of the Lptf3Loader.
    pub fn new(file_path: &Path, take_every: Option<u32>, return_all: bool) -> crate::Result<Self> {
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
    fn read_next_frame_header(&mut self) -> crate::Result<HdrRd> {
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
    fn seek_next_valid_frame(&mut self) -> crate::Result<Option<FrameHeader>> {
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
    pub fn get_next_frame_points(&mut self) -> crate::Result<Option<FullFrame>> {
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
pub struct FrameHeader {
    pub frame_index: u32,
    pub num_points: u32,
    pub x_offset: f64,
    pub z_offset: f64,
    pub x_res: f64,
    pub z_res: f64,

    /// This flag is set during the loading process to indicate that this frame consists entirely
    /// of points that would be skipped based on the `take_every` parameter.  If the loader is set
    /// to return all frames, this flag is what distinguishes a skipped frame from a valid frame
    pub skip: bool,
}

impl FrameHeader {
    fn from_buffer(buffer: &[u8; 24]) -> crate::Result<Self> {
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

fn read_raw_point<R: Read>(reader: &mut R, is_32bit: bool) -> crate::Result<(i32, i32)> {
    let (x, z) = if is_32bit {
        (read_i32(reader)?, read_i32(reader)?)
    } else {
        (read_i16(reader)? as i32, read_i16(reader)? as i32)
    };
    Ok((x, z))
}

fn read_res(buffer: &[u8]) -> crate::Result<f64> {
    // Convert from nanometers to millimeters
    Ok(u32::from_le_bytes(buffer[0..4].try_into()?) as f64 / 1_000_000.0)
}

fn read_offset(buffer: &[u8]) -> crate::Result<f64> {
    // Convert from micrometers to millimeters
    Ok(i32::from_le_bytes(buffer[0..4].try_into()?) as f64 / 1_000.0)
}

fn read_u16<R: Read>(reader: &mut R) -> crate::Result<u16> {
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u32<R: Read>(reader: &mut R) -> crate::Result<u32> {
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_i32<R: Read>(reader: &mut R) -> crate::Result<i32> {
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u8<R: Read>(reader: &mut R) -> crate::Result<u8> {
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_i16<R: Read>(reader: &mut R) -> crate::Result<i16> {
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(i16::from_le_bytes(buf))
}
