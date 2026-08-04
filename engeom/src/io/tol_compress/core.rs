use crate::Result;
use std::io::{Read, Write};

// The 2d half of this module is gone: curves now go through the `tol-compress` crate, and nothing
// else ever used it. The 3d half survives only until meshes follow.
mod indices;
mod points3;

pub use indices::*;
pub use points3::*;

pub(crate) fn read_u8<R: Read>(reader: &mut R) -> crate::Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

pub(crate) fn read_u32<R: Read>(reader: &mut R) -> crate::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn read_f64<R: Read>(reader: &mut R) -> crate::Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

pub(crate) fn bytes_for_tol(range: f64, tol: f64) -> u8 {
    let mut bytes = 2;
    while range / 2f64.powi((bytes as i32) * 8) > tol {
        bytes += 1;
    }

    bytes
}

pub(crate) fn write_var_width<W: Write>(writer: &mut W, fraction: f64, bytes: u8) -> Result<()> {
    // Using max_int = 2^(bytes*8) - 1 ensures u fits in `bytes` bytes (avoids overflow at
    // fraction = 1.0 which would produce 2^(bytes*8) and silently truncate to 0).
    let max_int = 2u64.pow((bytes as u32) * 8) - 1;
    let u = (fraction * max_int as f64).round() as u64;

    for i in 0..bytes {
        let byte = ((u >> (i * 8)) & 0xFF) as u8;
        writer.write_all(&[byte])?;
    }

    Ok(())
}

pub(crate) fn read_var_width<R: Read>(reader: &mut R, bytes: u8) -> Result<f64> {
    let mut u: u64 = 0;
    for i in 0..bytes {
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf)?;
        u |= (buf[0] as u64) << (i * 8);
    }
    let max_int = 2u64.pow((bytes as u32) * 8) - 1;
    Ok(u as f64 / max_int as f64)
}
