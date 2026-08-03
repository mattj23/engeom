//! Optional key-value metadata, carried by files and by individual items.
//!
//! The format stores geometry. Everything a particular caller wants to say _about_ that geometry,
//! this crate has no real business worrying about. To that purpose, metadata is an opt-in map,
//! and any item that doesn't have one doesn't have the extra space.
//!
//! Values are straight pass-throughs with no interpretation on the library side. There are only a
//! few basic binary types supported, because this was meant to be a smal, simple format and I'm
//! worried about feature creep.
//!
//! Keys are arbitrary UTF-8, but independent users of the same file will collide unless they
//! namespace. The convention is a dotted prefix owned by the writer, `engeom.chord_tol` or
//! `acme.scanner.serial`, but that only matters if files of these types end up being shared, if
//! you're just managing your own files feel free to do as you please.
//!
//! [`Metadata`] is a [`BTreeMap`], so entries serialize in sorted key order. Writing the same map
//! twice produces the same bytes, which might matter if you're saving files in a repository of some
//! type and they get rewritten with identical data.

use crate::error::{Error, Result};
use crate::raw::{read_u8, read_u32, write_u8, write_u32};
use std::collections::BTreeMap;
use std::io::{Read, Write};

/// A map of metadata keys to values.
///
/// A [`BTreeMap`] rather than a `HashMap` so serialization is deterministic, and rather than a
/// `Vec` of pairs so keys are unique by construction.
pub type Metadata = BTreeMap<String, Value>;

/// Largest key this decoder will allocate for.
const MAX_KEY_LEN: u32 = 64 * 1024;

/// Largest [`Value::Text`] or [`Value::Bytes`] payload this decoder will allocate for.
///
/// Metadata annotates geometry; bulk arrays belong in attribute blocks. This is generous for the
/// former and bounds what a corrupt length field can demand.
const MAX_PAYLOAD_LEN: u32 = 16 * 1024 * 1024;

const TAG_BOOL: u8 = 1;
const TAG_I64: u8 = 2;
const TAG_F64: u8 = 3;
const TAG_TEXT: u8 = 4;
const TAG_BYTES: u8 = 5;

/// A metadata value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A flag.
    Bool(bool),
    /// A signed integer: a count, an identifier, a timestamp.
    I64(i64),
    /// A real number: a tolerance, a scale, a physical quantity.
    F64(f64),
    /// UTF-8 text: units, a label, provenance.
    Text(String),
    /// Arbitrary bytes, for anything the tags above do not cover. The layout is the writer's to
    /// define and the reader's to understand.
    Bytes(Vec<u8>),
}

impl Value {
    /// The value if it is a [`Value::Bool`].
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Bool(v) => Some(*v),
            _ => None,
        }
    }

    /// The value if it is a [`Value::I64`].
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Value::I64(v) => Some(*v),
            _ => None,
        }
    }

    /// The value if it is a [`Value::F64`].
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::F64(v) => Some(*v),
            _ => None,
        }
    }

    /// The value if it is [`Value::Text`].
    pub fn as_text(&self) -> Option<&str> {
        match self {
            Value::Text(v) => Some(v),
            _ => None,
        }
    }

    /// The value if it is [`Value::Bytes`].
    pub fn as_bytes(&self) -> Option<&[u8]> {
        match self {
            Value::Bytes(v) => Some(v),
            _ => None,
        }
    }

    fn tag(&self) -> u8 {
        match self {
            Value::Bool(_) => TAG_BOOL,
            Value::I64(_) => TAG_I64,
            Value::F64(_) => TAG_F64,
            Value::Text(_) => TAG_TEXT,
            Value::Bytes(_) => TAG_BYTES,
        }
    }
}

impl From<bool> for Value {
    fn from(v: bool) -> Self {
        Value::Bool(v)
    }
}

impl From<i64> for Value {
    fn from(v: i64) -> Self {
        Value::I64(v)
    }
}

impl From<f64> for Value {
    fn from(v: f64) -> Self {
        Value::F64(v)
    }
}

impl From<String> for Value {
    fn from(v: String) -> Self {
        Value::Text(v)
    }
}

impl From<&str> for Value {
    fn from(v: &str) -> Self {
        Value::Text(v.to_string())
    }
}

impl From<Vec<u8>> for Value {
    fn from(v: Vec<u8>) -> Self {
        Value::Bytes(v)
    }
}

/// Write a metadata map. Callers set the presence flag; an empty map should not reach here.
pub(crate) fn write_metadata<W: Write>(writer: &mut W, metadata: &Metadata) -> Result<()> {
    let count = u16::try_from(metadata.len())
        .map_err(|_| Error::Malformed("metadata holds more entries than a u16 can count"))?;
    writer.write_all(&count.to_le_bytes())?;

    // BTreeMap iterates in sorted key order, which is what makes the output deterministic.
    for (key, value) in metadata {
        write_len_prefixed(writer, key.as_bytes(), MAX_KEY_LEN, "metadata key")?;
        write_u8(writer, value.tag())?;

        match value {
            Value::Bool(v) => write_u8(writer, u8::from(*v))?,
            Value::I64(v) => writer.write_all(&v.to_le_bytes())?,
            Value::F64(v) => writer.write_all(&v.to_le_bytes())?,
            Value::Text(v) => {
                write_len_prefixed(writer, v.as_bytes(), MAX_PAYLOAD_LEN, "metadata text")?
            }
            Value::Bytes(v) => write_len_prefixed(writer, v, MAX_PAYLOAD_LEN, "metadata bytes")?,
        }
    }

    Ok(())
}

/// Read a metadata map.
pub(crate) fn read_metadata<R: Read>(reader: &mut R) -> Result<Metadata> {
    let mut count_bytes = [0u8; 2];
    reader.read_exact(&mut count_bytes)?;
    let count = u16::from_le_bytes(count_bytes);

    let mut out = Metadata::new();
    for _ in 0..count {
        let key_bytes = read_len_prefixed(reader, MAX_KEY_LEN, "metadata key")?;
        let key = String::from_utf8(key_bytes)
            .map_err(|_| Error::Malformed("metadata key is not valid UTF-8"))?;

        let value = match read_u8(reader)? {
            TAG_BOOL => Value::Bool(read_u8(reader)? != 0),
            TAG_I64 => {
                let mut b = [0u8; 8];
                reader.read_exact(&mut b)?;
                Value::I64(i64::from_le_bytes(b))
            }
            TAG_F64 => {
                let mut b = [0u8; 8];
                reader.read_exact(&mut b)?;
                Value::F64(f64::from_le_bytes(b))
            }
            TAG_TEXT => {
                let bytes = read_len_prefixed(reader, MAX_PAYLOAD_LEN, "metadata text")?;
                Value::Text(
                    String::from_utf8(bytes)
                        .map_err(|_| Error::Malformed("metadata text is not valid UTF-8"))?,
                )
            }
            TAG_BYTES => Value::Bytes(read_len_prefixed(
                reader,
                MAX_PAYLOAD_LEN,
                "metadata bytes",
            )?),
            _ => return Err(Error::Malformed("unknown metadata value type")),
        };

        // A duplicate would silently overwrite, losing whichever entry lost the race.
        if out.insert(key, value).is_some() {
            return Err(Error::Malformed("metadata contains a duplicate key"));
        }
    }

    Ok(out)
}

fn write_len_prefixed<W: Write>(
    writer: &mut W,
    bytes: &[u8],
    max: u32,
    what: &'static str,
) -> Result<()> {
    let len = u32::try_from(bytes.len())
        .ok()
        .filter(|&n| n <= max)
        .ok_or(match what {
            "metadata key" => Error::Malformed("metadata key is too long to store"),
            "metadata text" => Error::Malformed("metadata text is too long to store"),
            _ => Error::Malformed("metadata byte string is too long to store"),
        })?;
    write_u32(writer, len)?;
    writer.write_all(bytes)?;
    Ok(())
}

fn read_len_prefixed<R: Read>(reader: &mut R, max: u32, what: &'static str) -> Result<Vec<u8>> {
    let len = read_u32(reader)?;
    if len > max {
        return Err(match what {
            "metadata key" => Error::Malformed("metadata key length is implausible"),
            "metadata text" => Error::Malformed("metadata text length is implausible"),
            _ => Error::Malformed("metadata byte string length is implausible"),
        });
    }
    let mut bytes = vec![0u8; len as usize];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn round_trip(metadata: &Metadata) -> Metadata {
        let mut buf = Vec::new();
        write_metadata(&mut buf, metadata).unwrap();

        let mut cursor = Cursor::new(&buf);
        let back = read_metadata(&mut cursor).unwrap();
        assert_eq!(
            cursor.position() as usize,
            buf.len(),
            "decoder left bytes unread"
        );
        back
    }

    #[test]
    fn every_value_type_round_trips() {
        let mut m = Metadata::new();
        m.insert("flag".into(), Value::Bool(true));
        m.insert("off".into(), Value::Bool(false));
        m.insert("count".into(), Value::I64(-9_223_372_036_854_775_808));
        m.insert("tol".into(), Value::F64(1e-4));
        m.insert("units".into(), Value::Text("mm".into()));
        m.insert("blob".into(), Value::Bytes(vec![0, 1, 2, 255]));

        assert_eq!(round_trip(&m), m);
    }

    #[test]
    fn accessors_return_only_their_own_type() {
        assert_eq!(Value::F64(1.5).as_f64(), Some(1.5));
        assert_eq!(Value::F64(1.5).as_i64(), None);
        assert_eq!(Value::Text("mm".into()).as_text(), Some("mm"));
        assert_eq!(Value::Text("mm".into()).as_bool(), None);
        assert_eq!(Value::Bool(true).as_bool(), Some(true));
        assert_eq!(Value::I64(7).as_i64(), Some(7));
        assert_eq!(Value::Bytes(vec![1, 2]).as_bytes(), Some(&[1u8, 2][..]));
    }

    #[test]
    fn conversions_are_ergonomic() {
        let mut m = Metadata::new();
        m.insert("a".into(), 1e-4.into());
        m.insert("b".into(), "mm".into());
        m.insert("c".into(), 42i64.into());
        m.insert("d".into(), true.into());

        assert_eq!(m["a"].as_f64(), Some(1e-4));
        assert_eq!(m["b"].as_text(), Some("mm"));
        assert_eq!(m["c"].as_i64(), Some(42));
        assert_eq!(m["d"].as_bool(), Some(true));
    }

    /// Fixtures get committed and rewritten. Two writes of the same map must produce the same
    /// bytes, or every rewrite shows up as a spurious diff.
    #[test]
    fn output_is_deterministic() {
        let mut m = Metadata::new();
        for key in ["zulu", "alpha", "mike", "bravo", "yankee"] {
            m.insert(key.into(), Value::I64(key.len() as i64));
        }

        let mut first = Vec::new();
        let mut second = Vec::new();
        write_metadata(&mut first, &m).unwrap();
        write_metadata(&mut second, &m).unwrap();
        assert_eq!(first, second);

        // And insertion order must not matter, only content.
        let mut reordered = Metadata::new();
        for key in ["mike", "yankee", "alpha", "zulu", "bravo"] {
            reordered.insert(key.into(), Value::I64(key.len() as i64));
        }
        let mut third = Vec::new();
        write_metadata(&mut third, &reordered).unwrap();
        assert_eq!(first, third);
    }

    /// The crate stores what it is given. A caller's `NaN` is a caller's business.
    #[test]
    fn values_are_stored_not_interpreted() {
        let mut m = Metadata::new();
        m.insert("nan".into(), Value::F64(f64::NAN));
        m.insert("inf".into(), Value::F64(f64::NEG_INFINITY));
        m.insert("empty".into(), Value::Text(String::new()));
        m.insert("nothing".into(), Value::Bytes(Vec::new()));
        m.insert(String::new(), Value::Bool(false));

        let back = round_trip(&m);
        assert!(back["nan"].as_f64().unwrap().is_nan());
        assert_eq!(back["inf"].as_f64(), Some(f64::NEG_INFINITY));
        assert_eq!(back["empty"].as_text(), Some(""));
        assert_eq!(back["nothing"].as_bytes(), Some(&[][..]));
        assert_eq!(back[""].as_bool(), Some(false));
    }

    #[test]
    fn an_empty_map_is_two_bytes() {
        let mut buf = Vec::new();
        write_metadata(&mut buf, &Metadata::new()).unwrap();
        assert_eq!(buf, vec![0, 0]);
        assert!(read_metadata(&mut Cursor::new(&buf)).unwrap().is_empty());
    }

    #[test]
    fn non_ascii_keys_and_values_survive() {
        let mut m = Metadata::new();
        m.insert(
            "secci\u{f3}n".into(),
            Value::Text("\u{65ad}\u{9762}".into()),
        );

        assert_eq!(round_trip(&m), m);
    }

    #[test]
    fn an_unknown_value_type_is_rejected() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u16.to_le_bytes());
        write_u32(&mut buf, 1).unwrap();
        buf.push(b'k');
        buf.push(99);

        assert!(matches!(
            read_metadata(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    /// A duplicate key would silently overwrite on insert, losing one of the two entries.
    #[test]
    fn a_duplicate_key_is_rejected() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&2u16.to_le_bytes());
        for _ in 0..2 {
            write_u32(&mut buf, 1).unwrap();
            buf.push(b'k');
            buf.push(TAG_BOOL);
            buf.push(1);
        }

        assert!(matches!(
            read_metadata(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn an_absurd_length_does_not_allocate() {
        for tag in [TAG_TEXT, TAG_BYTES] {
            let mut buf = Vec::new();
            buf.extend_from_slice(&1u16.to_le_bytes());
            write_u32(&mut buf, 1).unwrap();
            buf.push(b'k');
            buf.push(tag);
            write_u32(&mut buf, u32::MAX).unwrap();

            assert!(
                matches!(
                    read_metadata(&mut Cursor::new(&buf)),
                    Err(Error::Malformed(_))
                ),
                "tag {tag} should run out of input, not memory"
            );
        }

        // And the key length itself.
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u16.to_le_bytes());
        write_u32(&mut buf, u32::MAX).unwrap();
        assert!(matches!(
            read_metadata(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn invalid_utf8_is_rejected() {
        // An invalid text value.
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u16.to_le_bytes());
        write_u32(&mut buf, 1).unwrap();
        buf.push(b'k');
        buf.push(TAG_TEXT);
        write_u32(&mut buf, 2).unwrap();
        buf.extend_from_slice(&[0xFF, 0xFE]);

        assert!(matches!(
            read_metadata(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));

        // And an invalid key.
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u16.to_le_bytes());
        write_u32(&mut buf, 2).unwrap();
        buf.extend_from_slice(&[0xFF, 0xFE]);
        assert!(matches!(
            read_metadata(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn truncation_is_an_error() {
        let mut m = Metadata::new();
        m.insert("tol".into(), Value::F64(1e-4));
        m.insert("units".into(), Value::Text("mm".into()));

        let mut buf = Vec::new();
        write_metadata(&mut buf, &m).unwrap();

        for cut in 0..buf.len() {
            assert!(
                read_metadata(&mut Cursor::new(&buf[..cut])).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
    }
}
