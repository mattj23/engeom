//! Fixed-width little-endian scalars for the byte-oriented parts of a stream.
//!
//! Headers stay byte-aligned and plainly readable; only block payloads are bit-packed.

use crate::error::Result;
use std::io::{Read, Write};

/// Ceiling on how much a decoder will pre-allocate from a length field before seeing any data.
///
/// A corrupt or hostile stream can claim four billion elements. Growing the vector as elements
/// actually arrive costs a few reallocations on genuinely large blocks and bounds the damage on
/// bad ones.
pub(crate) const MAX_PREALLOC: usize = 1 << 16;

pub(crate) fn write_u8<W: Write>(writer: &mut W, value: u8) -> Result<()> {
    writer.write_all(&[value])?;
    Ok(())
}

pub(crate) fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

pub(crate) fn write_u32<W: Write>(writer: &mut W, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

pub(crate) fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn write_f64<W: Write>(writer: &mut W, value: f64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

pub(crate) fn read_f64<R: Read>(reader: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn scalars_round_trip() {
        let mut buf = Vec::new();
        write_u8(&mut buf, 0xAB).unwrap();
        write_u32(&mut buf, 0xDEAD_BEEF).unwrap();
        write_f64(&mut buf, -12.375).unwrap();

        assert_eq!(buf.len(), 13);

        let mut c = Cursor::new(&buf);
        assert_eq!(read_u8(&mut c).unwrap(), 0xAB);
        assert_eq!(read_u32(&mut c).unwrap(), 0xDEAD_BEEF);
        assert_eq!(read_f64(&mut c).unwrap(), -12.375);
    }

    #[test]
    fn byte_order_is_little_endian() {
        let mut buf = Vec::new();
        write_u32(&mut buf, 1).unwrap();
        assert_eq!(buf, vec![1, 0, 0, 0]);
    }
}
