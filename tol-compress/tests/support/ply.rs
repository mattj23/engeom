//! A minimal reader for the binary PLY fixtures in `tests/data`.
//!
//! Hand-rolled rather than pulling a PLY crate in as a dev-dependency, so that `tol-compress` has
//! no dependencies at all, not even for development. It handles exactly the subset the fixtures
//! use, which is the subset they were normalized into, and rejects anything else rather than
//! guessing.
//!
//! Shared by `tests/real_fixtures.rs` and `examples/size_report.rs` through `#[path]`.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// A triangle mesh loaded from a fixture.
pub struct Mesh {
    pub points: Vec<[f64; 3]>,
    pub faces: Vec<[u32; 3]>,
}

impl Mesh {
    /// The span of the mesh along each axis.
    pub fn extents(&self) -> [f64; 3] {
        let mut mins = [f64::INFINITY; 3];
        let mut maxs = [f64::NEG_INFINITY; 3];
        for p in &self.points {
            for i in 0..3 {
                mins[i] = mins[i].min(p[i]);
                maxs[i] = maxs[i].max(p[i]);
            }
        }
        std::array::from_fn(|i| maxs[i] - mins[i])
    }
}

/// The directory holding the committed fixtures.
pub fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Load a fixture by file name.
pub fn load(name: &str) -> Mesh {
    let path = data_dir().join(name);
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    parse(&bytes).unwrap_or_else(|e| panic!("parsing {}: {e}", path.display()))
}

fn parse(bytes: &[u8]) -> Result<Mesh, String> {
    let marker = b"end_header\n";
    let header_end = bytes
        .windows(marker.len())
        .position(|w| w == marker)
        .ok_or("no end_header")?
        + marker.len();

    let header = std::str::from_utf8(&bytes[..header_end]).map_err(|e| e.to_string())?;

    if !header.contains("format binary_little_endian 1.0") {
        return Err("only binary_little_endian 1.0 is supported".into());
    }

    let mut vertex_count = 0usize;
    let mut face_count = 0usize;
    let mut vertex_props: Vec<&str> = Vec::new();
    let mut element = "";

    for line in header.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        match parts.as_slice() {
            ["element", "vertex", n] => {
                element = "vertex";
                vertex_count = n.parse().map_err(|_| "bad vertex count")?;
            }
            ["element", "face", n] => {
                element = "face";
                face_count = n.parse().map_err(|_| "bad face count")?;
            }
            ["property", ty, name] if element == "vertex" => {
                if *ty != "float" {
                    return Err(format!("vertex property {name} is {ty}, expected float"));
                }
                vertex_props.push(name);
            }
            ["property", "list", "uchar", "int", _] if element == "face" => {}
            ["property", ..] if element == "face" => {
                return Err("unexpected face property".into());
            }
            _ => {}
        }
    }

    if vertex_props != ["x", "y", "z"] {
        return Err(format!("expected exactly x, y, z; found {vertex_props:?}"));
    }

    let mut cursor = header_end;
    let mut take = |n: usize| -> Result<&[u8], String> {
        let end = cursor + n;
        if end > bytes.len() {
            return Err("truncated".into());
        }
        let slice = &bytes[cursor..end];
        cursor = end;
        Ok(slice)
    };

    let mut points = Vec::with_capacity(vertex_count);
    for _ in 0..vertex_count {
        let raw = take(12)?;
        points.push(std::array::from_fn(|i| {
            let mut b = [0u8; 4];
            b.copy_from_slice(&raw[i * 4..i * 4 + 4]);
            f32::from_le_bytes(b) as f64
        }));
    }

    let mut faces = Vec::with_capacity(face_count);
    for _ in 0..face_count {
        let n = take(1)?[0];
        if n != 3 {
            return Err(format!("only triangles are supported, found a {n}-gon"));
        }
        let raw = take(12)?;
        faces.push(std::array::from_fn(|i| {
            let mut b = [0u8; 4];
            b.copy_from_slice(&raw[i * 4..i * 4 + 4]);
            u32::from_le_bytes(b)
        }));
    }

    for f in &faces {
        for &i in f {
            if i as usize >= points.len() {
                return Err("face references a vertex that does not exist".into());
            }
        }
    }

    Ok(Mesh { points, faces })
}
