//! Reader for the Point Cloud Library (PCL) Point Cloud Data (`.pcd`) format. Open3D, ROS tooling,
//! and many scanner and lidar vendors also write this format.

use crate::geom3::attributes3::{Attr3, PointAttrSet3};
use crate::{Point3, Result, UnitVec3, Vector3};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Load the points and per-point attributes from a PCD file.
///
/// See `read_pcd_points` for the supported encodings, how fields map onto attributes, and which
/// points and fields are dropped or refused.
///
/// # Arguments
///
/// * `path`: the path to the `.pcd` file
///
/// returns: `Result<(Vec<Point3>, PointAttrSet3)>`
pub fn load_pcd_points(path: &Path) -> Result<(Vec<Point3>, PointAttrSet3)> {
    let file = BufReader::new(File::open(path)?);
    read_pcd_points(file)
}

/// Read the points and per-point attributes of a PCD file from a buffered source.
///
/// # Format
///
/// A PCD file consists of an ASCII header followed by a payload. The header has one keyword per
/// line and ends at the `DATA` line. Lines starting with `#` are comments.
///
/// - `VERSION`: ignored because the payload layout does not depend on it
/// - `FIELDS`: the name of each field
/// - `SIZE`: the size in bytes of one value of each field
/// - `TYPE`: `I` (signed integer), `U` (unsigned integer), or `F` (floating point) for each field
/// - `COUNT`: the number of values each field holds per point. Each field defaults to 1 value if
///   `COUNT` is absent.
/// - `WIDTH`, `HEIGHT`: the grid dimensions of an organized cloud, or the point count and 1 for an
///   unorganized one. `HEIGHT` defaults to 1 if absent.
/// - `VIEWPOINT`: the sensor pose as a translation and a quaternion. It is validated but not
///   applied, because the points are stored in their own frame (PCL does not apply it on load
///   either).
/// - `POINTS`: the point count, defaulting to `WIDTH * HEIGHT` if absent
/// - `DATA`: one of the three payload encodings described below
///
/// The payload encodings are:
///
/// - `ascii`: one point per line. Values are separated by whitespace and appear in field order,
///   with each field contributing its `COUNT` values.
/// - `binary`: one little-endian record per point, with the fields packed in order and no padding.
/// - `binary_compressed`: a `u32` compressed size, a `u32` uncompressed size, and then LZF
///   compressed data. The decompressed data is field-major: every point's values for the first
///   field, then every point's values for the second, and so on.
///
/// # Mapping to engeom
///
/// - `x`, `y`, `z` are required and become the point positions. A point whose position is not
///   finite is **dropped**. PCL marks the invalid cells of an organized cloud with NaN positions,
///   and the unordered point set this reader produces has no grid to hold them in.
/// - `normal_x`, `normal_y`, `normal_z` become the point normals. A file must have all three or
///   none, and a kept point with a non-finite or zero-length normal is an error.
/// - `rgb` or `rgba` becomes the point colors. Both are packed into four bytes as `0xAARRGGBB`
///   whether the field is declared `F` or `U`, and the alpha byte is discarded. In an ASCII
///   payload, an all-digit token that fits in a `u32` is the packed integer (which is how PCL writes
///   it). Any other token is a float whose bit pattern provides the packing.
/// - `_` is PCL's padding field and is skipped.
/// - Every other field goes into the open attribute map. Floating point fields become
///   `Attr3::Scalar`. Integer fields become `Attr3::Label` if every kept value fits in a `u32`, or
///   `Attr3::Scalar` if every kept value is exactly representable as an `f64`, and are refused
///   otherwise rather than rounded.
/// - A field with a `COUNT` above 1 is split into one attribute per value, named `name_0`,
///   `name_1`, and so on. Multi-value PCD fields are usually histograms or descriptors, and storing
///   a triple as an `Attr3::Vector` would cause it to rotate as a spatial direction when the cloud
///   is transformed.
/// - A field whose name is reserved by the attribute set (PCL's `label`, for example) is refused.
///
/// # Arguments
///
/// * `source`: a buffered reader positioned at the start of the PCD file
///
/// returns: `Result<(Vec<Point3>, PointAttrSet3)>`
pub fn read_pcd_points<R: BufRead>(mut source: R) -> Result<(Vec<Point3>, PointAttrSet3)> {
    let header = read_header(&mut source)?;

    let raw = match header.data {
        DataKind::Ascii => read_ascii(&mut source, &header)?,
        DataKind::Binary => read_binary(&mut source, &header)?,
        DataKind::BinaryCompressed => read_binary_compressed(&mut source, &header)?,
    };

    route(&header, &raw)
}

// ===============================================================================================
// Header
// ===============================================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FieldKind {
    Signed,
    Unsigned,
    Float,
}

impl FieldKind {
    fn parse(token: &str) -> Result<Self> {
        match token {
            "I" | "i" => Ok(Self::Signed),
            "U" | "u" => Ok(Self::Unsigned),
            "F" | "f" => Ok(Self::Float),
            _ => Err(format!("PCD header TYPE '{token}' is not one of I, U, or F").into()),
        }
    }

    fn letter(self) -> char {
        match self {
            Self::Signed => 'I',
            Self::Unsigned => 'U',
            Self::Float => 'F',
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DataKind {
    Ascii,
    Binary,
    BinaryCompressed,
}

impl DataKind {
    fn parse(values: &[&str]) -> Result<Self> {
        let [value] = values else {
            return Err("PCD header DATA line must name a single encoding".into());
        };

        match value.to_ascii_lowercase().as_str() {
            "ascii" => Ok(Self::Ascii),
            "binary" => Ok(Self::Binary),
            "binary_compressed" => Ok(Self::BinaryCompressed),
            other => Err(format!(
                "PCD DATA encoding '{other}' is not one of ascii, binary, or binary_compressed"
            )
            .into()),
        }
    }
}

/// One declared field: its name, the byte size and kind of each value, and how many values each
/// point holds.
#[derive(Debug, Clone)]
struct FieldDef {
    name: String,
    size: usize,
    kind: FieldKind,
    count: usize,
}

impl FieldDef {
    /// The number of bytes this field occupies in one point's record.
    fn stride(&self) -> usize {
        self.size * self.count
    }

    /// The bytes of one value in this field's field-major buffer.
    fn value_bytes<'a>(&self, raw: &'a [u8], point: usize, slot: usize) -> &'a [u8] {
        let start = (point * self.count + slot) * self.size;
        &raw[start..start + self.size]
    }

    /// Decode one value as a float, whatever the field's declared kind.
    fn float(&self, bytes: &[u8]) -> f64 {
        match (self.kind, self.size) {
            (FieldKind::Float, 4) => f32::from_le_bytes(bytes.try_into().unwrap()) as f64,
            (FieldKind::Float, _) => f64::from_le_bytes(bytes.try_into().unwrap()),
            _ => self.integer(bytes) as f64,
        }
    }

    /// Decode one value of an integer field. An `i128` holds every signed and unsigned integer
    /// the format can declare without loss.
    fn integer(&self, bytes: &[u8]) -> i128 {
        let mut wide = [0u8; 8];
        wide[..self.size].copy_from_slice(bytes);
        let unsigned = u64::from_le_bytes(wide);

        match self.kind {
            FieldKind::Unsigned => unsigned as i128,
            FieldKind::Signed => {
                let shift = 64 - 8 * self.size as u32;
                (((unsigned << shift) as i64) >> shift) as i128
            }
            FieldKind::Float => unreachable!("integer() is only called on integer fields"),
        }
    }
}

#[derive(Debug)]
struct PcdHeader {
    fields: Vec<FieldDef>,
    points: usize,
    data: DataKind,

    /// The number of bytes in one point's record, summed over the fields.
    point_stride: usize,

    /// The number of bytes in the whole uncompressed payload.
    payload_len: usize,
}

fn read_header<R: BufRead>(source: &mut R) -> Result<PcdHeader> {
    let mut names: Option<Vec<String>> = None;
    let mut sizes: Option<Vec<usize>> = None;
    let mut kinds: Option<Vec<FieldKind>> = None;
    let mut counts: Option<Vec<usize>> = None;
    let mut width: Option<usize> = None;
    let mut height: Option<usize> = None;
    let mut points: Option<usize> = None;
    let mut buffer = Vec::new();

    // Read the header byte by byte up to each newline because a binary payload begins immediately
    // after the DATA line.
    let data = loop {
        buffer.clear();
        if source.read_until(b'\n', &mut buffer)? == 0 {
            return Err("PCD file ended before its header reached a DATA line".into());
        }

        let line = String::from_utf8_lossy(&buffer);
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut tokens = line.split_whitespace();
        let Some(keyword) = tokens.next() else {
            continue;
        };
        let values: Vec<&str> = tokens.collect();

        match keyword.to_ascii_uppercase().as_str() {
            "FIELDS" => names = Some(values.iter().map(|v| v.to_string()).collect()),
            "SIZE" => sizes = Some(parse_list("SIZE", &values)?),
            "TYPE" => {
                kinds = Some(
                    values
                        .iter()
                        .map(|v| FieldKind::parse(v))
                        .collect::<Result<_>>()?,
                )
            }
            "COUNT" => counts = Some(parse_list("COUNT", &values)?),
            "WIDTH" => width = Some(parse_single("WIDTH", &values)?),
            "HEIGHT" => height = Some(parse_single("HEIGHT", &values)?),
            "POINTS" => points = Some(parse_single("POINTS", &values)?),
            "VIEWPOINT" => check_viewpoint(&values)?,
            "DATA" => break DataKind::parse(&values)?,
            // VERSION, and any keyword this reader does not know, cannot change how the payload is
            // laid out, so they are skipped rather than refused.
            _ => {}
        }
    };

    let names = names.ok_or("PCD header has no FIELDS line")?;
    let sizes = sizes.ok_or("PCD header has no SIZE line")?;
    let kinds = kinds.ok_or("PCD header has no TYPE line")?;
    let counts = counts.unwrap_or_else(|| vec![1; names.len()]);

    if names.is_empty() {
        return Err("PCD header FIELDS line declares no fields".into());
    }

    for (keyword, len) in [
        ("SIZE", sizes.len()),
        ("TYPE", kinds.len()),
        ("COUNT", counts.len()),
    ] {
        if len != names.len() {
            return Err(format!(
                "PCD header lists {} FIELDS but {len} {keyword} values",
                names.len()
            )
            .into());
        }
    }

    let mut fields: Vec<FieldDef> = Vec::with_capacity(names.len());
    let mut point_stride = 0usize;
    for (((name, size), kind), count) in names.into_iter().zip(sizes).zip(kinds).zip(counts) {
        let valid = match kind {
            FieldKind::Float => matches!(size, 4 | 8),
            _ => matches!(size, 1 | 2 | 4 | 8),
        };
        if !valid {
            return Err(format!(
                "PCD field '{name}' is declared as TYPE {} with SIZE {size}, which is not a valid \
                 combination",
                kind.letter()
            )
            .into());
        }

        if count == 0 {
            return Err(format!("PCD field '{name}' declares COUNT 0").into());
        }

        if name != "_" && fields.iter().any(|f| f.name == name) {
            return Err(format!("PCD header declares the field '{name}' more than once").into());
        }

        point_stride = size
            .checked_mul(count)
            .and_then(|s| point_stride.checked_add(s))
            .ok_or_else(|| format!("PCD field '{name}' declares a COUNT too large to address"))?;

        fields.push(FieldDef {
            name,
            size,
            kind,
            count,
        });
    }

    let grid = width.map(|w| w.checked_mul(height.unwrap_or(1)));
    let points = match (points, grid) {
        (Some(p), Some(g)) => {
            if g != Some(p) {
                return Err(format!(
                    "PCD header declares POINTS {p}, but WIDTH {} x HEIGHT {} does not match it",
                    width.unwrap_or(0),
                    height.unwrap_or(1)
                )
                .into());
            }
            p
        }
        (Some(p), None) => p,
        (None, Some(g)) => g.ok_or("PCD header WIDTH x HEIGHT is too large to address")?,
        (None, None) => {
            return Err(
                "PCD header has neither POINTS nor WIDTH, so the point count is unknown".into(),
            );
        }
    };

    let payload_len = points
        .checked_mul(point_stride)
        .ok_or("PCD header declares a payload too large to address")?;

    Ok(PcdHeader {
        fields,
        points,
        data,
        point_stride,
        payload_len,
    })
}

fn parse_list(keyword: &str, values: &[&str]) -> Result<Vec<usize>> {
    let mut parsed = Vec::with_capacity(values.len());
    for v in values {
        let n = v.parse::<usize>().map_err(|_| {
            format!("PCD header {keyword} value '{v}' is not a non-negative integer")
        })?;
        parsed.push(n);
    }
    Ok(parsed)
}

fn parse_single(keyword: &str, values: &[&str]) -> Result<usize> {
    let [value] = values else {
        return Err(format!("PCD header {keyword} line must hold a single value").into());
    };

    Ok(value.parse::<usize>().map_err(|_| {
        format!("PCD header {keyword} value '{value}' is not a non-negative integer")
    })?)
}

fn check_viewpoint(values: &[&str]) -> Result<()> {
    if values.len() != 7 || values.iter().any(|v| v.parse::<f64>().is_err()) {
        return Err(format!(
            "PCD header VIEWPOINT must be seven numbers (a translation and a quaternion), but is \
             '{}'",
            values.join(" ")
        )
        .into());
    }
    Ok(())
}

// ===============================================================================================
// Payloads
// ===============================================================================================
//
// Every payload reader produces the same thing: one field-major byte buffer per declared field,
// holding that field's values for every point in the little-endian layout of its declared type.
// That keeps the field routing below independent of the encoding.

/// Read an ASCII payload, encoding each token into its field's binary layout.
fn read_ascii<R: BufRead>(source: &mut R, header: &PcdHeader) -> Result<Vec<Vec<u8>>> {
    let n = header.points;
    let mut fields: Vec<Vec<u8>> = vec![Vec::new(); header.fields.len()];
    let mut line = String::new();
    let mut read = 0;

    while read < n {
        line.clear();
        if source.read_line(&mut line)? == 0 {
            return Err(format!(
                "PCD ascii payload holds only {read} of the {n} points the header declares"
            )
            .into());
        }

        let text = line.trim();
        if text.is_empty() || text.starts_with('#') {
            continue;
        }

        let mut tokens = text.split_whitespace();
        for (def, out) in header.fields.iter().zip(fields.iter_mut()) {
            for _ in 0..def.count {
                let token = tokens.next().ok_or_else(|| {
                    format!("PCD ascii point {read} has fewer values than its fields declare")
                })?;
                push_ascii_value(def, token, out)
                    .map_err(|e| format!("PCD ascii point {read}, field '{}': {e}", def.name))?;
            }
        }

        if tokens.next().is_some() {
            return Err(
                format!("PCD ascii point {read} has more values than its fields declare").into(),
            );
        }

        read += 1;
    }

    Ok(fields)
}

/// Parse one ASCII token and append it to a field buffer in the field's binary layout.
fn push_ascii_value(
    def: &FieldDef,
    token: &str,
    out: &mut Vec<u8>,
) -> std::result::Result<(), String> {
    let invalid = || {
        format!(
            "'{token}' is not a valid {} value of size {}",
            def.kind.letter(),
            def.size
        )
    };

    match def.kind {
        FieldKind::Float if def.size == 4 => {
            // A packed color written as its integer bit pattern, which is how PCL writes it.
            if is_packed_color(&def.name)
                && token.bytes().all(|b| b.is_ascii_digit())
                && let Ok(bits) = token.parse::<u32>()
            {
                out.extend_from_slice(&bits.to_le_bytes());
            } else {
                let v: f32 = token.parse().map_err(|_| invalid())?;
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        FieldKind::Float => {
            let v: f64 = token.parse().map_err(|_| invalid())?;
            out.extend_from_slice(&v.to_le_bytes());
        }
        FieldKind::Signed => {
            let v: i64 = token.parse().map_err(|_| invalid())?;
            let bits = 8 * def.size as u32;
            if bits < 64 && (v < -(1i64 << (bits - 1)) || v >= (1i64 << (bits - 1))) {
                return Err(invalid());
            }
            out.extend_from_slice(&v.to_le_bytes()[..def.size]);
        }
        FieldKind::Unsigned => {
            let v: u64 = token.parse().map_err(|_| invalid())?;
            let bits = 8 * def.size as u32;
            if bits < 64 && v >= (1u64 << bits) {
                return Err(invalid());
            }
            out.extend_from_slice(&v.to_le_bytes()[..def.size]);
        }
    }

    Ok(())
}

/// Read a binary payload and split its point-major records into field-major buffers.
fn read_binary<R: Read>(source: &mut R, header: &PcdHeader) -> Result<Vec<Vec<u8>>> {
    let data = read_len(source, header.payload_len, "binary")?;

    let mut fields: Vec<Vec<u8>> = header
        .fields
        .iter()
        .map(|f| Vec::with_capacity(f.stride() * header.points))
        .collect();

    for record in data.chunks_exact(header.point_stride) {
        let mut offset = 0;
        for (def, out) in header.fields.iter().zip(fields.iter_mut()) {
            out.extend_from_slice(&record[offset..offset + def.stride()]);
            offset += def.stride();
        }
    }

    Ok(fields)
}

/// Read a binary_compressed payload, whose decompressed data is already field-major.
fn read_binary_compressed<R: Read>(source: &mut R, header: &PcdHeader) -> Result<Vec<Vec<u8>>> {
    // PCL cannot write a compressed payload for an empty cloud, so there may be nothing to read.
    if header.points == 0 {
        return Ok(vec![Vec::new(); header.fields.len()]);
    }

    let sizes = read_len(source, 8, "binary_compressed")?;
    let compressed_len = u32::from_le_bytes([sizes[0], sizes[1], sizes[2], sizes[3]]) as usize;
    let uncompressed_len = u32::from_le_bytes([sizes[4], sizes[5], sizes[6], sizes[7]]) as usize;

    if uncompressed_len != header.payload_len {
        return Err(format!(
            "PCD binary_compressed payload declares {uncompressed_len} uncompressed bytes, but {} \
             points of {} bytes each need {}",
            header.points, header.point_stride, header.payload_len
        )
        .into());
    }

    let compressed = read_len(source, compressed_len, "binary_compressed")?;
    let data = lzf_decompress(&compressed, uncompressed_len)?;

    let mut fields = Vec::with_capacity(header.fields.len());
    let mut offset = 0;
    for def in header.fields.iter() {
        let len = def.stride() * header.points;
        fields.push(data[offset..offset + len].to_vec());
        offset += len;
    }

    Ok(fields)
}

/// Read a known number of bytes, failing if the source ends first.
fn read_len<R: Read>(source: &mut R, len: usize, encoding: &str) -> Result<Vec<u8>> {
    // Reading through `take` rather than preallocating means a header which declares far more
    // data than the file holds fails on the missing bytes instead of on a huge allocation.
    let mut buffer = Vec::new();
    source.take(len as u64).read_to_end(&mut buffer)?;

    if buffer.len() != len {
        return Err(format!(
            "PCD {encoding} payload ended after {} bytes, but {len} were expected",
            buffer.len()
        )
        .into());
    }

    Ok(buffer)
}

/// Decompress an LZF stream, as produced by liblzf's `lzf_compress`, which PCL uses for its
/// binary_compressed payload.
///
/// The stream is a sequence of control bytes. A control byte below 32 starts a literal run of
/// `ctrl + 1` bytes copied from the input. Any other control byte starts a back reference: its top
/// three bits are the length minus 2 (with 7 meaning an extra length byte follows), and its low
/// five bits and the next byte form the distance minus 1 back into the output. A reference may
/// overlap the bytes it is producing.
fn lzf_decompress(input: &[u8], out_len: usize) -> Result<Vec<u8>> {
    const TRUNCATED: &str = "PCD binary_compressed payload ends in the middle of an LZF sequence";

    // The expected length comes from the file, so it only guides the initial allocation up to a
    // bound set by the data actually present.
    let mut out: Vec<u8> = Vec::with_capacity(out_len.min(input.len().saturating_mul(4)));
    let mut ip = 0;

    while ip < input.len() {
        let ctrl = input[ip] as usize;
        ip += 1;

        if ctrl < 32 {
            let len = ctrl + 1;
            let literal = input.get(ip..ip + len).ok_or(TRUNCATED)?;
            if out.len() + len > out_len {
                return Err(lzf_overrun(out_len));
            }
            out.extend_from_slice(literal);
            ip += len;
        } else {
            let mut len = ctrl >> 5;
            if len == 7 {
                len += *input.get(ip).ok_or(TRUNCATED)? as usize;
                ip += 1;
            }
            let low = *input.get(ip).ok_or(TRUNCATED)? as usize;
            ip += 1;

            let len = len + 2;
            let distance = ((ctrl & 0x1f) << 8) + low + 1;
            if distance > out.len() {
                return Err(format!(
                    "PCD binary_compressed payload has an LZF reference {distance} bytes back \
                     from byte {}, before the start of the data",
                    out.len()
                )
                .into());
            }
            if out.len() + len > out_len {
                return Err(lzf_overrun(out_len));
            }

            // Byte by byte, since the reference may overlap the bytes it is producing.
            let start = out.len() - distance;
            for k in 0..len {
                let b = out[start + k];
                out.push(b);
            }
        }
    }

    if out.len() != out_len {
        return Err(format!(
            "PCD binary_compressed payload decompressed to {} bytes, but {out_len} were expected",
            out.len()
        )
        .into());
    }

    Ok(out)
}

fn lzf_overrun(out_len: usize) -> Box<dyn std::error::Error> {
    format!(
        "PCD binary_compressed payload decompresses to more than the {out_len} bytes it declares"
    )
    .into()
}

// ===============================================================================================
// Field routing
// ===============================================================================================

/// Field names consumed by positions, normals, or colors rather than carried into the open
/// attribute map.
const CONSUMED: &[&str] = &[
    "x", "y", "z", "normal_x", "normal_y", "normal_z", "rgb", "rgba",
];

fn is_packed_color(name: &str) -> bool {
    name == "rgb" || name == "rgba"
}

fn route(header: &PcdHeader, raw: &[Vec<u8>]) -> Result<(Vec<Point3>, PointAttrSet3)> {
    let n = header.points;
    let find = |name: &str| header.fields.iter().position(|f| f.name == name);

    let coordinate = |name: &str| -> Result<Vec<f64>> {
        let i = find(name).ok_or_else(|| format!("PCD file has no '{name}' field"))?;
        single_values(&header.fields[i], &raw[i], n)
    };
    let x = coordinate("x")?;
    let y = coordinate("y")?;
    let z = coordinate("z")?;

    let keep: Vec<usize> = (0..n)
        .filter(|&i| x[i].is_finite() && y[i].is_finite() && z[i].is_finite())
        .collect();
    let m = keep.len();

    let points = keep
        .iter()
        .map(|&i| Point3::new(x[i], y[i], z[i]))
        .collect();
    let mut attrs = PointAttrSet3::empty();

    match (find("normal_x"), find("normal_y"), find("normal_z")) {
        (Some(ix), Some(iy), Some(iz)) => {
            let nx = single_values(&header.fields[ix], &raw[ix], n)?;
            let ny = single_values(&header.fields[iy], &raw[iy], n)?;
            let nz = single_values(&header.fields[iz], &raw[iz], n)?;

            let mut normals = Vec::with_capacity(m);
            for &p in &keep {
                let v = Vector3::new(nx[p], ny[p], nz[p]);
                if !v.iter().all(|c| c.is_finite()) {
                    return Err(format!(
                        "PCD point {p} has a finite position but a non-finite normal ({}, {}, {})",
                        v.x, v.y, v.z
                    )
                    .into());
                }
                let unit = UnitVec3::try_new(v, 1.0e-12).ok_or_else(|| {
                    format!(
                        "PCD point {p} has a normal of zero length, which cannot be a direction"
                    )
                })?;
                normals.push(unit);
            }
            attrs.set_normals(Some(normals), m)?;
        }
        (None, None, None) => {}
        _ => {
            return Err(
                "PCD file has only some of the 'normal_x', 'normal_y', and 'normal_z' \
                        fields, but normals need all three"
                    .into(),
            );
        }
    }

    let color_field = match (find("rgb"), find("rgba")) {
        (Some(_), Some(_)) => {
            return Err(
                "PCD file has both 'rgb' and 'rgba' fields, so its point colors are ambiguous"
                    .into(),
            );
        }
        (i, None) | (None, i) => i,
    };
    if let Some(i) = color_field {
        let def = &header.fields[i];
        if def.size != 4 || def.count != 1 {
            return Err(format!(
                "PCD field '{}' must be a single 4-byte packed color, but declares SIZE {} and \
                 COUNT {}",
                def.name, def.size, def.count
            )
            .into());
        }

        let colors = keep
            .iter()
            .map(|&p| {
                let bits = u32::from_le_bytes(def.value_bytes(&raw[i], p, 0).try_into().unwrap());
                [(bits >> 16) as u8, (bits >> 8) as u8, bits as u8]
            })
            .collect();
        attrs.set_colors(Some(colors), m)?;
    }

    for (def, bytes) in header.fields.iter().zip(raw.iter()) {
        if def.name == "_" || CONSUMED.contains(&def.name.as_str()) {
            continue;
        }

        for slot in 0..def.count {
            let name = if def.count == 1 {
                def.name.clone()
            } else {
                format!("{}_{slot}", def.name)
            };

            if attrs.attr(&name).is_some() {
                return Err(format!(
                    "PCD field '{}' would be stored as the attribute '{name}', which another \
                     field already uses",
                    def.name
                )
                .into());
            }

            let attr = open_attr(def, bytes, slot, &keep)
                .map_err(|e| format!("PCD field '{}' {e}", def.name))?;
            attrs.insert_attr(&name, attr, m).map_err(|e| {
                format!(
                    "PCD field '{}' cannot be stored as a point attribute: {e}",
                    def.name
                )
            })?;
        }
    }

    Ok((points, attrs))
}

/// Decode every point's value of a field which must hold a single value per point.
fn single_values(def: &FieldDef, raw: &[u8], n: usize) -> Result<Vec<f64>> {
    if def.count != 1 {
        return Err(format!(
            "PCD field '{}' declares COUNT {}, but it must hold a single value per point",
            def.name, def.count
        )
        .into());
    }

    Ok((0..n)
        .map(|p| def.float(def.value_bytes(raw, p, 0)))
        .collect())
}

/// Build the open attribute for one value slot of a field, over the kept points only.
fn open_attr(def: &FieldDef, raw: &[u8], slot: usize, keep: &[usize]) -> Result<Attr3> {
    if def.kind == FieldKind::Float {
        let values = keep
            .iter()
            .map(|&p| def.float(def.value_bytes(raw, p, slot)))
            .collect();
        return Ok(Attr3::Scalar(values));
    }

    let values: Vec<i128> = keep
        .iter()
        .map(|&p| def.integer(def.value_bytes(raw, p, slot)))
        .collect();

    if values.iter().all(|v| (0..=u32::MAX as i128).contains(v)) {
        return Ok(Attr3::Label(values.iter().map(|v| *v as u32).collect()));
    }

    if let Some(v) = values.iter().find(|v| (**v as f64) as i128 != **v) {
        return Err(format!(
            "holds the integer {v}, which is outside the range of a label and cannot be stored \
             as a scalar without rounding"
        )
        .into());
    }

    Ok(Attr3::Scalar(values.iter().map(|v| *v as f64).collect()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::get_test_file_path;
    use approx::assert_relative_eq;
    use std::io::Cursor;

    fn error_of(bytes: impl AsRef<[u8]>) -> String {
        match read_pcd_points(Cursor::new(bytes.as_ref())) {
            Ok(_) => panic!("expected the PCD data to be refused"),
            Err(e) => e.to_string(),
        }
    }

    // --- The same cloud in every encoding ------------------------------------------------------

    /// One point of the sample cloud.
    struct SampleRow {
        xyz: [f32; 3],
        normal: [f32; 3],
        rgb: u32,
        intensity: f32,
        ring: u16,
        offset: i16,
        hist: [f32; 2],
    }

    /// Covers every kind of field the router distinguishes: positions, normals, a float-declared
    /// packed color, a float scalar, an unsigned label, a signed field with a negative value, PCL
    /// padding, and a multi-value field.
    const SAMPLE_HEADER: &str = "# .PCD v0.7 - Point Cloud Data file format\n\
        VERSION 0.7\n\
        FIELDS x y z normal_x normal_y normal_z rgb intensity ring offset _ hist\n\
        SIZE 4 4 4 4 4 4 4 4 2 2 1 4\n\
        TYPE F F F F F F F F U I U F\n\
        COUNT 1 1 1 1 1 1 1 1 1 1 3 2\n\
        WIDTH 3\n\
        HEIGHT 1\n\
        VIEWPOINT 0 0 0 1 0 0 0\n\
        POINTS 3\n";

    fn sample_rows() -> Vec<SampleRow> {
        vec![
            SampleRow {
                xyz: [1.0, 2.0, 3.0],
                normal: [0.0, 0.0, 1.0],
                rgb: 0xFF10_2030,
                intensity: 0.25,
                ring: 7,
                offset: -3,
                hist: [0.5, 1.5],
            },
            // An invalid point, whose values must not leak into any attribute.
            SampleRow {
                xyz: [f32::NAN; 3],
                normal: [f32::NAN; 3],
                rgb: 0,
                intensity: 9.0,
                ring: 60000,
                offset: 0,
                hist: [0.0, 0.0],
            },
            SampleRow {
                xyz: [-4.5, 0.125, 1.0e3],
                normal: [0.6, 0.8, 0.0],
                rgb: 0xFFFF_8000,
                intensity: 1.0,
                ring: 2,
                offset: 12,
                hist: [2.0, -1.0],
            },
        ]
    }

    /// The bytes of each field of one row, in declaration order.
    fn sample_fields(r: &SampleRow) -> Vec<Vec<u8>> {
        let mut fields: Vec<Vec<u8>> = Vec::new();
        for v in r.xyz.iter().chain(r.normal.iter()) {
            fields.push(v.to_le_bytes().to_vec());
        }
        fields.push(r.rgb.to_le_bytes().to_vec());
        fields.push(r.intensity.to_le_bytes().to_vec());
        fields.push(r.ring.to_le_bytes().to_vec());
        fields.push(r.offset.to_le_bytes().to_vec());
        fields.push(vec![0u8; 3]);
        fields.push([r.hist[0].to_le_bytes(), r.hist[1].to_le_bytes()].concat());
        fields
    }

    fn sample_ascii() -> Vec<u8> {
        let mut text = format!("{SAMPLE_HEADER}DATA ascii\n");
        for r in sample_rows() {
            text += &format!(
                "{} {} {} {} {} {} {} {} {} {} 0 0 0 {} {}\n",
                r.xyz[0],
                r.xyz[1],
                r.xyz[2],
                r.normal[0],
                r.normal[1],
                r.normal[2],
                r.rgb,
                r.intensity,
                r.ring,
                r.offset,
                r.hist[0],
                r.hist[1]
            );
        }
        text.into_bytes()
    }

    fn sample_binary() -> Vec<u8> {
        let mut bytes = format!("{SAMPLE_HEADER}DATA binary\n").into_bytes();
        for r in sample_rows() {
            bytes.extend(sample_fields(&r).concat());
        }
        bytes
    }

    /// The compressed payload here uses only LZF literal runs, which is a valid stream. Back
    /// references are covered by the LZF tests and by the pypcd4 fixture.
    fn sample_compressed() -> Vec<u8> {
        let rows: Vec<Vec<Vec<u8>>> = sample_rows().iter().map(sample_fields).collect();
        let mut field_major = Vec::new();
        for f in 0..rows[0].len() {
            for row in rows.iter() {
                field_major.extend_from_slice(&row[f]);
            }
        }

        let mut bytes = format!("{SAMPLE_HEADER}DATA binary_compressed\n").into_bytes();
        let stream = lzf_literals(&field_major);
        bytes.extend((stream.len() as u32).to_le_bytes());
        bytes.extend((field_major.len() as u32).to_le_bytes());
        bytes.extend(stream);
        bytes
    }

    fn lzf_literals(data: &[u8]) -> Vec<u8> {
        let mut stream = Vec::new();
        for chunk in data.chunks(32) {
            stream.push((chunk.len() - 1) as u8);
            stream.extend_from_slice(chunk);
        }
        stream
    }

    fn check_sample(bytes: Vec<u8>) -> Result<()> {
        let (points, attrs) = read_pcd_points(Cursor::new(bytes))?;

        assert_eq!(
            points,
            vec![Point3::new(1.0, 2.0, 3.0), Point3::new(-4.5, 0.125, 1.0e3)]
        );

        let normals = attrs.normals().expect("the sample declares normals");
        assert_relative_eq!(normals[0].into_inner(), Vector3::new(0.0, 0.0, 1.0));
        assert_relative_eq!(
            normals[1].into_inner(),
            Vector3::new(0.6, 0.8, 0.0),
            epsilon = 1.0e-6
        );

        assert_eq!(
            attrs.colors(),
            Some(&[[0x10, 0x20, 0x30], [0xFF, 0x80, 0x00]][..])
        );

        let scalar = |name: &str| {
            attrs
                .attr(name)
                .and_then(|a| a.as_scalar())
                .map(|v| v.to_vec())
        };
        assert_eq!(scalar("intensity"), Some(vec![0.25, 1.0]));
        assert_eq!(scalar("offset"), Some(vec![-3.0, 12.0]));
        assert_eq!(scalar("hist_0"), Some(vec![0.5, 2.0]));
        assert_eq!(scalar("hist_1"), Some(vec![1.5, -1.0]));
        assert_eq!(
            attrs.attr("ring").and_then(|a| a.as_label()),
            Some(&[7, 2][..])
        );

        // The padding field, and the fields consumed by typed attributes, leave nothing behind.
        let mut names: Vec<&str> = attrs.attr_names().collect();
        names.sort();
        assert_eq!(names, ["hist_0", "hist_1", "intensity", "offset", "ring"]);

        Ok(())
    }

    #[test]
    fn reads_the_sample_as_ascii() -> Result<()> {
        check_sample(sample_ascii())
    }

    #[test]
    fn reads_the_sample_as_binary() -> Result<()> {
        check_sample(sample_binary())
    }

    #[test]
    fn reads_the_sample_as_binary_compressed() -> Result<()> {
        check_sample(sample_compressed())
    }

    // --- A file from a real encoder ------------------------------------------------------------

    /// A `binary_compressed` file written by pypcd4. Its LZF stream comes from liblzf and contains
    /// back references. The file holds the 453 vertices of `bun_zipper_res4.ply` with two NaN points
    /// spliced in at file indices 100 and 301. Its colors and ring numbers are derived from each
    /// vertex's index in the PLY, and its normals point away from the vertex centroid.
    #[test]
    fn reads_a_binary_compressed_file_from_pypcd4() -> Result<()> {
        let (points, attrs) = load_pcd_points(&get_test_file_path("bun_zipper_res4.pcd"))?;

        assert_eq!(points.len(), 453);
        assert_relative_eq!(points[0].x, -0.0312216, epsilon = 1.0e-6);
        assert_relative_eq!(points[0].y, 0.126304, epsilon = 1.0e-6);
        assert_relative_eq!(points[0].z, 0.00514924, epsilon = 1.0e-6);

        // Any misalignment from dropping the NaN points would break these relationships.
        let colors = attrs.colors().expect("the file declares an 'rgb' field");
        let ring = attrs
            .attr("ring")
            .and_then(|a| a.as_label())
            .expect("a u16 field should become labels");
        for i in 0..points.len() {
            let expected = [(i % 256) as u8, (7 * i % 256) as u8, (13 * i % 256) as u8];
            assert_eq!(colors[i], expected, "color of vertex {i}");
            assert_eq!(ring[i], (i % 16) as u32, "ring of vertex {i}");
        }

        assert_eq!(attrs.normals().map(|n| n.len()), Some(453));
        let intensity = attrs.attr("intensity").and_then(|a| a.as_scalar());
        assert_eq!(intensity.map(|v| v[0]), Some(0.5));

        Ok(())
    }

    #[cfg(feature = "ply")]
    #[test]
    fn pypcd4_file_agrees_with_the_ply_it_was_made_from() -> Result<()> {
        let (expected, ply_attrs, _) =
            crate::io::load_ply_points(&get_test_file_path("bun_zipper_res4.ply"))?;
        let (points, attrs) = load_pcd_points(&get_test_file_path("bun_zipper_res4.pcd"))?;

        // Both files store float32 positions and intensities, so they agree without tolerance.
        assert_eq!(points, expected);
        assert_eq!(
            attrs.attr("intensity").and_then(|a| a.as_scalar()),
            ply_attrs.attr("intensity").and_then(|a| a.as_scalar())
        );

        let centroid = expected
            .iter()
            .fold(Vector3::zeros(), |acc, p| acc + p.coords)
            / expected.len() as f64;
        let normals = attrs.normals().expect("the file declares normals");
        for (p, n) in expected.iter().zip(normals) {
            let direction = (p.coords - centroid).normalize();
            assert_relative_eq!(n.into_inner(), direction, epsilon = 1.0e-6);
        }

        Ok(())
    }

    // --- Header and payload details ------------------------------------------------------------

    /// Covers omitted VERSION, COUNT, and POINTS fields, comments, blank lines, CRLF endings, and
    /// an organized cloud with an invalid cell.
    #[test]
    fn header_defaults_and_layout_tolerance() -> Result<()> {
        let text = "# written by hand\r\nFIELDS x y z\r\nSIZE 8 8 8\r\nTYPE F F F\r\n\r\n\
                    WIDTH 2\r\nHEIGHT 2\r\nDATA ascii\r\n1 2 3\r\nnan nan nan\r\n\r\n4 5 6\r\n\
                    7 8 9\r\n";
        let (points, attrs) = read_pcd_points(Cursor::new(text))?;

        assert_eq!(
            points,
            vec![
                Point3::new(1.0, 2.0, 3.0),
                Point3::new(4.0, 5.0, 6.0),
                Point3::new(7.0, 8.0, 9.0)
            ]
        );
        assert!(attrs.is_empty());

        Ok(())
    }

    #[test]
    fn reads_an_empty_cloud() -> Result<()> {
        let text = "FIELDS x y z\nSIZE 4 4 4\nTYPE F F F\nWIDTH 0\nHEIGHT 1\nPOINTS 0\n\
                    DATA binary_compressed\n";
        let (points, attrs) = read_pcd_points(Cursor::new(text))?;

        assert!(points.is_empty());
        assert!(attrs.is_empty());

        Ok(())
    }

    /// Confirms that a packed color has the same value when written as an unsigned integer, as
    /// PCL's integer representation of a float-declared color, or as float text.
    #[test]
    fn packed_colors_read_the_same_however_they_are_written() -> Result<()> {
        let header = |kind: char| {
            format!("FIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F {kind}\nWIDTH 1\nDATA ascii\n")
        };

        let variants = [
            format!("{}0 0 0 {}\n", header('U'), 0x7F10_2030u32),
            format!("{}0 0 0 {}\n", header('F'), 0x7F10_2030u32),
            // An all-digit token too large for a u32, which can only be a float.
            format!("{}0 0 0 {}\n", header('F'), f32::from_bits(0x7F10_2030)),
            format!("{}0 0 0 {}\n", header('F'), f32::from_bits(0xFF10_2030)),
        ];

        for text in variants {
            let (_, attrs) = read_pcd_points(Cursor::new(text.as_str()))?;
            assert_eq!(attrs.colors(), Some(&[[0x10, 0x20, 0x30]][..]), "{text}");
        }

        Ok(())
    }

    #[test]
    fn wide_integers_become_scalars_only_when_exact() -> Result<()> {
        let header = "FIELDS x y z stamp\nSIZE 4 4 4 8\nTYPE F F F U\nWIDTH 2\nDATA ascii\n";

        // This value is outside the label range but is exactly representable as a double. The
        // dropped point's value is not exactly representable and must not affect the decision.
        let text = format!("{header}1 2 3 5000000000\nnan nan nan 9007199254740993\n");
        let (points, attrs) = read_pcd_points(Cursor::new(text))?;
        assert_eq!(points.len(), 1);
        assert_eq!(
            attrs.attr("stamp").and_then(|a| a.as_scalar()),
            Some(&[5.0e9][..])
        );

        let text = format!("{header}1 2 3 9007199254740993\n4 5 6 1\n");
        let err = error_of(text);
        assert!(err.contains("9007199254740993"), "{err}");

        Ok(())
    }

    #[test]
    fn refuses_malformed_or_unrepresentable_files() {
        let xyz = "FIELDS x y z\nSIZE 4 4 4\nTYPE F F F\n";
        let normals =
            "FIELDS x y z normal_x normal_y normal_z\nSIZE 4 4 4 4 4 4\nTYPE F F F F F F\n";

        let cases = [
            (
                "FIELDS x y z\nSIZE 4 4\nTYPE F F F\nWIDTH 1\nDATA ascii\n1 2 3\n".to_string(),
                "2 SIZE values",
            ),
            (
                "FIELDS x y z\nSIZE 4 4 2\nTYPE F F F\nWIDTH 1\nDATA ascii\n1 2 3\n".to_string(),
                "not a valid combination",
            ),
            (
                format!("{xyz}WIDTH 2\nHEIGHT 1\nPOINTS 3\nDATA ascii\n"),
                "POINTS 3",
            ),
            (format!("{xyz}WIDTH 1\n"), "DATA line"),
            (format!("{xyz}WIDTH 1\nDATA binary_lz4\n"), "binary_lz4"),
            (
                format!("{xyz}VIEWPOINT 0 0 0\nWIDTH 1\nDATA ascii\n1 2 3\n"),
                "VIEWPOINT",
            ),
            (
                format!("{xyz}DATA ascii\n1 2 3\n"),
                "neither POINTS nor WIDTH",
            ),
            (
                "FIELDS x y\nSIZE 4 4\nTYPE F F\nWIDTH 1\nDATA ascii\n1 2\n".to_string(),
                "no 'z' field",
            ),
            (
                "FIELDS x y z x\nSIZE 4 4 4 4\nTYPE F F F F\nWIDTH 1\nDATA ascii\n1 2 3 4\n"
                    .to_string(),
                "more than once",
            ),
            (
                "FIELDS x y z normal_x\nSIZE 4 4 4 4\nTYPE F F F F\nWIDTH 1\nDATA ascii\n1 2 3 4\n"
                    .to_string(),
                "need all three",
            ),
            (
                "FIELDS x y z rgb rgba\nSIZE 4 4 4 4 4\nTYPE F F F U U\nWIDTH 1\nDATA ascii\n\
                 1 2 3 4 5\n"
                    .to_string(),
                "both 'rgb' and 'rgba'",
            ),
            (
                format!("{xyz}WIDTH 2\nDATA ascii\n1 2 3\n"),
                "only 1 of the 2",
            ),
            (format!("{xyz}WIDTH 1\nDATA ascii\n1 2\n"), "fewer values"),
            (
                format!("{xyz}WIDTH 1\nDATA ascii\n1 2 3 4\n"),
                "more values",
            ),
            (format!("{xyz}WIDTH 1\nDATA ascii\n1 2 three\n"), "'three'"),
            (
                "FIELDS x y z ring\nSIZE 4 4 4 1\nTYPE F F F U\nWIDTH 1\nDATA ascii\n1 2 3 256\n"
                    .to_string(),
                "'256'",
            ),
            (
                format!("{xyz}WIDTH 1\nDATA binary\n\0\0"),
                "binary payload ended",
            ),
            (
                format!("{xyz}WIDTH 1\nDATA binary_compressed\n\x08\0\0\0\x0c\0\0\0\0"),
                "binary_compressed payload ended",
            ),
            (
                format!("{xyz}WIDTH 1\nDATA binary_compressed\n\x01\0\0\0\x0b\0\0\0\0"),
                "declares 11 uncompressed bytes",
            ),
            (
                "FIELDS x y z label\nSIZE 4 4 4 4\nTYPE F F F U\nWIDTH 1\nDATA ascii\n1 2 3 4\n"
                    .to_string(),
                "'label' is a reserved attribute name",
            ),
            (
                format!("{normals}WIDTH 1\nDATA ascii\n1 2 3 nan 0 0\n"),
                "non-finite normal",
            ),
            (
                format!("{normals}WIDTH 1\nDATA ascii\n1 2 3 0 0 0\n"),
                "zero length",
            ),
        ];

        for (text, expected) in cases {
            let err = error_of(&text);
            assert!(
                err.contains(expected),
                "expected '{expected}' in: {err}\nfor:\n{text}"
            );
        }
    }

    // --- LZF -----------------------------------------------------------------------------------

    #[test]
    fn lzf_back_references_copy_from_earlier_output() -> Result<()> {
        // "abc" as a literal run, then a reference of length 5 at distance 3, which overlaps the
        // bytes it produces.
        assert_eq!(
            lzf_decompress(&[0x02, b'a', b'b', b'c', 0x60, 0x02], 8)?,
            b"abcabcab"
        );

        // One literal, then a reference whose length needs the extra byte (7 + 11 + 2).
        assert_eq!(
            lzf_decompress(&[0x00, b'z', 0xE0, 11, 0x00], 21)?,
            vec![b'z'; 21]
        );

        // A distance which needs the high bits in the control byte ((1 << 8) + 0x21 + 1 = 290).
        let literal: Vec<u8> = (0..300).map(|i| (i % 251) as u8).collect();
        let mut stream = lzf_literals(&literal);
        stream.extend([0x41, 0x21]);
        let out = lzf_decompress(&stream, 304)?;
        assert_eq!(&out[..300], &literal[..]);
        assert_eq!(&out[300..], &literal[10..14]);

        Ok(())
    }

    #[test]
    fn lzf_refuses_malformed_streams() {
        let cases: [(&[u8], usize, &str); 4] = [
            (&[0x20, 0x00], 3, "before the start"),
            (&[0x05, 1, 2], 6, "middle of an LZF sequence"),
            (&[0x01, 1, 2], 3, "decompressed to 2 bytes"),
            (&[0x02, 1, 2, 3], 2, "more than the 2 bytes"),
        ];

        for (stream, len, expected) in cases {
            let err = lzf_decompress(stream, len).unwrap_err().to_string();
            assert!(err.contains(expected), "expected '{expected}' in: {err}");
        }
    }
}
