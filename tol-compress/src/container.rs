//! File framing: the header that turns blocks into a file, and the per-item metadata shared by
//! every kind of container.
//!
//! # Every container is a collection
//!
//! A container declares one kind and an item count, and holds that many items. A mesh file
//! normally holds exactly one, but the layout does not distinguish that case, which permanently
//! removes the question of what to do when a file needs a second item. Reading a single item is an
//! API convenience rather than a separate format, and it **errors** rather than silently taking the
//! first when a file holds several.
//!
//! # Layout
//!
//! ```text
//! container header (12 bytes, plus metadata if present)
//!   4 bytes    magic "TOLC"
//!   u8         format version
//!   u8         kind
//!   u8         compression, always 0, reserved
//!   u8         file flags: bit0 has metadata
//!   u32        item count
//!   [metadata] if the flag is set
//!
//! per item:
//!   u8         item flags: bit0 closed (polylines), bit1 has name, bit2 has metadata
//!   [metadata] if the flag is set
//!   [name]     u32 byte length then UTF-8, if the flag is set
//!   u8         attribute block count, always 0, reserved
//!   ...        the kind's geometry blocks
//! ```
//!
//! Two fields are reserved and rejected if nonzero: `compression` and the per-item attribute
//! count. Both cost one byte and let a later version add the feature without a version bump, which
//! matters because files written once should never need rewriting.
//!
//! # Metadata rather than named fields
//!
//! Anything a caller wants to record *about* the geometry goes in a [`crate::metadata::Metadata`]
//! map, at file level or per item. Nothing domain-specific gets a flag bit of its own, so a caller
//! who does not need chord tolerances or scanner serials pays nothing for them. An absent map
//! costs zero bytes.
//!
//! Kind records the dimension even though the points block records it too. The kind byte is what
//! lets [`probe`] identify a file without parsing any geometry, and the readers check that the two
//! agree.

use crate::error::{Error, Result};
use crate::metadata::{Metadata, read_metadata, write_metadata};
use crate::raw::{read_u8, read_u32, write_u8, write_u32};
use std::io::{Read, Write};

/// Leading bytes of every tol-compress file.
pub const MAGIC: [u8; 4] = *b"TOLC";

/// The format version this build writes.
pub const VERSION: u8 = 2;

/// What a container holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    /// Triangle meshes in 3D.
    Mesh3,
    /// Ordered polylines in 2D.
    Polyline2,
    /// Ordered polylines in 3D.
    Polyline3,
    /// Unordered point sets in 2D.
    Cloud2,
    /// Unordered point sets in 3D.
    Cloud3,
}

impl Kind {
    /// The byte written to the header.
    pub fn as_byte(self) -> u8 {
        match self {
            Kind::Mesh3 => 1,
            Kind::Polyline2 => 2,
            Kind::Polyline3 => 3,
            Kind::Cloud2 => 4,
            Kind::Cloud3 => 5,
        }
    }

    /// Recover a kind from its header byte.
    ///
    /// # Errors
    ///
    /// [`Error::Malformed`] for any byte this version does not define.
    pub fn from_byte(byte: u8) -> Result<Self> {
        match byte {
            1 => Ok(Kind::Mesh3),
            2 => Ok(Kind::Polyline2),
            3 => Ok(Kind::Polyline3),
            4 => Ok(Kind::Cloud2),
            5 => Ok(Kind::Cloud3),
            _ => Err(Error::Malformed("unknown container kind")),
        }
    }

    /// The dimension of the points this kind stores.
    pub fn dimension(self) -> usize {
        match self {
            Kind::Polyline2 | Kind::Cloud2 => 2,
            Kind::Mesh3 | Kind::Polyline3 | Kind::Cloud3 => 3,
        }
    }

    /// The conventional file extension, without a leading dot.
    pub fn extension(self) -> &'static str {
        match self {
            Kind::Mesh3 => "tcmesh",
            Kind::Polyline2 => "tccurve2",
            Kind::Polyline3 => "tccurve3",
            Kind::Cloud2 => "tccloud2",
            Kind::Cloud3 => "tccloud3",
        }
    }
}

/// What a container header says about the file.
#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    /// Format version the file was written at.
    pub version: u8,
    /// What the file holds.
    pub kind: Kind,
    /// How many items follow.
    pub count: u32,
    /// File-level metadata, empty when the file carries none.
    pub metadata: Metadata,
}

/// File-level flag bits.
mod file_flags {
    /// A metadata map follows the item count.
    pub const HAS_METADATA: u8 = 1 << 0;
    /// Every bit this version understands.
    pub const KNOWN: u8 = HAS_METADATA;
}

/// Write a container header.
pub fn write_header<W: Write>(
    writer: &mut W,
    kind: Kind,
    count: u32,
    metadata: &Metadata,
) -> Result<()> {
    writer.write_all(&MAGIC)?;
    write_u8(writer, VERSION)?;
    write_u8(writer, kind.as_byte())?;
    write_u8(writer, 0)?;

    let flags = if metadata.is_empty() {
        0
    } else {
        file_flags::HAS_METADATA
    };
    write_u8(writer, flags)?;
    write_u32(writer, count)?;

    if !metadata.is_empty() {
        write_metadata(writer, metadata)?;
    }

    Ok(())
}

/// Read a container header without caring what kind it turns out to be.
///
/// This is the cheap identification path: eleven bytes, no geometry touched. Use it to decide which
/// reader a file wants, rather than guessing from its extension.
///
/// ```
/// use tol_compress::{Kind, Polyline3, polyline, probe};
///
/// let mut buf = Vec::new();
/// let lines = [Polyline3::new(vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], false)];
/// polyline::write_to(&mut buf, &lines, 1e-4)?;
///
/// let header = probe(&mut buf.as_slice())?;
/// assert_eq!(header.kind, Kind::Polyline3);
/// assert_eq!(header.kind.dimension(), 3);
/// assert_eq!(header.kind.extension(), "tccurve3");
/// assert_eq!(header.count, 1);
/// # Ok::<(), tol_compress::Error>(())
/// ```
///
/// # Errors
///
/// [`Error::BadMagic`], [`Error::UnsupportedVersion`], [`Error::UnsupportedCompression`], or
/// [`Error::Malformed`] for an unknown kind.
pub fn probe<R: Read>(reader: &mut R) -> Result<Header> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if magic != MAGIC {
        return Err(Error::BadMagic);
    }

    let version = read_u8(reader)?;
    if version != VERSION {
        return Err(Error::UnsupportedVersion(version));
    }

    let kind = Kind::from_byte(read_u8(reader)?)?;

    let compression = read_u8(reader)?;
    if compression != 0 {
        return Err(Error::UnsupportedCompression(compression));
    }

    let flags = read_u8(reader)?;
    if flags & !file_flags::KNOWN != 0 {
        return Err(Error::Malformed("file sets unknown flag bits"));
    }

    let count = read_u32(reader)?;

    let metadata = if flags & file_flags::HAS_METADATA != 0 {
        read_metadata(reader)?
    } else {
        Metadata::new()
    };

    Ok(Header {
        version,
        kind,
        count,
        metadata,
    })
}

/// Read a container header and require it to be the kind the caller expects.
///
/// # Errors
///
/// Everything [`probe`] can return, plus [`Error::Malformed`] when the file holds a different kind.
pub fn read_header<R: Read>(reader: &mut R, expected: Kind) -> Result<Header> {
    let header = probe(reader)?;
    if header.kind != expected {
        return Err(Error::Malformed(
            "container holds a different kind than requested",
        ));
    }
    Ok(header)
}

/// An item that can be identified by name within a collection.
pub trait Named {
    /// The item's name, if it has one.
    fn name(&self) -> Option<&str>;
}

/// The first item in a collection with the given name.
///
/// The format permits duplicate names, because enforcing uniqueness is a caller's policy rather
/// than a container's. This returns the first match; iterate directly when every match matters.
pub fn find_by_name<'a, T: Named>(items: &'a [T], name: &str) -> Option<&'a T> {
    items.iter().find(|item| item.name() == Some(name))
}

/// The framing every item carries, whatever kind it is.
///
/// Grouped rather than scattered because the three kind modules that consume it all write the same
/// preamble in the same order, and an inconsistency between them would be a format bug rather than
/// a local one.
pub(crate) mod item {
    // Consumed by the kind modules, which arrive in the next stage.
    #![allow(dead_code)]

    use super::*;

    /// Largest name this decoder will allocate for.
    ///
    /// Names identify items; anything approaching this is corruption or hostility rather than a
    /// real name, and a `u32` length field would otherwise let a malformed file ask for four
    /// gigabytes.
    pub const MAX_NAME_LEN: u32 = 64 * 1024;

    /// The polyline is closed. Meaningless for other kinds.
    pub const CLOSED: u8 = 1 << 0;
    /// A name follows.
    pub const HAS_NAME: u8 = 1 << 1;
    /// A metadata map follows.
    pub const HAS_METADATA: u8 = 1 << 2;
    /// Every flag bit this version understands.
    pub const KNOWN_FLAGS: u8 = CLOSED | HAS_NAME | HAS_METADATA;

    /// Everything an item carries before its geometry.
    #[derive(Debug, Default)]
    pub struct Preamble {
        /// Optional identifier.
        pub name: Option<String>,
        /// Optional metadata, empty when absent.
        pub metadata: Metadata,
        /// Whether the polyline closes. Always false for kinds that have no such notion.
        pub closed: bool,
    }

    /// Write an item's preamble: flags, metadata, name, and the reserved attribute count.
    ///
    /// Every kind writes this, in this order. Centralized because a divergence between kinds would
    /// be a format bug rather than a local one.
    pub fn write_preamble<W: Write>(
        writer: &mut W,
        name: Option<&str>,
        metadata: &Metadata,
        closed: bool,
    ) -> Result<()> {
        let mut flags = 0u8;
        if closed {
            flags |= CLOSED;
        }
        if name.is_some() {
            flags |= HAS_NAME;
        }
        if !metadata.is_empty() {
            flags |= HAS_METADATA;
        }

        write_flags(writer, flags)?;
        if !metadata.is_empty() {
            write_metadata(writer, metadata)?;
        }
        write_name(writer, name)?;
        write_attribute_count(writer)?;

        Ok(())
    }

    /// Read an item's preamble.
    ///
    /// `closable` is false for kinds with no notion of closure, which then reject the flag rather
    /// than ignoring it.
    pub fn read_preamble<R: Read>(reader: &mut R, closable: bool) -> Result<Preamble> {
        let flags = read_flags(reader)?;
        if !closable && flags & CLOSED != 0 {
            return Err(Error::Malformed(
                "item sets the closed flag, which only applies to polylines",
            ));
        }

        let metadata = if flags & HAS_METADATA != 0 {
            read_metadata(reader)?
        } else {
            Metadata::new()
        };
        let name = read_name(reader, flags)?;
        read_attribute_count(reader)?;

        Ok(Preamble {
            name,
            metadata,
            closed: flags & CLOSED != 0,
        })
    }

    /// Write an item's flag byte.
    pub fn write_flags<W: Write>(writer: &mut W, flags: u8) -> Result<()> {
        write_u8(writer, flags)
    }

    /// Read an item's flag byte, rejecting bits this version does not define.
    ///
    /// Rejecting unknown bits rather than masking them off means a file using a future feature
    /// fails loudly instead of being silently misread as one that does not.
    pub fn read_flags<R: Read>(reader: &mut R) -> Result<u8> {
        let flags = read_u8(reader)?;
        if flags & !KNOWN_FLAGS != 0 {
            return Err(Error::Malformed("item sets unknown flag bits"));
        }
        Ok(flags)
    }

    /// Write an item name, if there is one. The caller must have set [`HAS_NAME`] to match.
    pub fn write_name<W: Write>(writer: &mut W, name: Option<&str>) -> Result<()> {
        let Some(name) = name else {
            return Ok(());
        };
        let bytes = name.as_bytes();
        let len = u32::try_from(bytes.len())
            .ok()
            .filter(|&n| n <= MAX_NAME_LEN)
            .ok_or(Error::Malformed("item name is too long to store"))?;
        write_u32(writer, len)?;
        writer.write_all(bytes)?;
        Ok(())
    }

    /// Read an item name if `flags` says one is present.
    pub fn read_name<R: Read>(reader: &mut R, flags: u8) -> Result<Option<String>> {
        if flags & HAS_NAME == 0 {
            return Ok(None);
        }

        let len = read_u32(reader)?;
        if len > MAX_NAME_LEN {
            return Err(Error::Malformed("item name length is implausible"));
        }

        let mut bytes = vec![0u8; len as usize];
        reader.read_exact(&mut bytes)?;
        String::from_utf8(bytes)
            .map(Some)
            .map_err(|_| Error::Malformed("item name is not valid UTF-8"))
    }

    /// Write the reserved attribute-block count, which is always zero for now.
    pub fn write_attribute_count<W: Write>(writer: &mut W) -> Result<()> {
        write_u8(writer, 0)
    }

    /// Read and validate the reserved attribute-block count.
    ///
    /// # Errors
    ///
    /// [`Error::UnsupportedAttributes`] if the item declares any, since this version has no
    /// definition of what an attribute block contains and skipping past one is impossible.
    pub fn read_attribute_count<R: Read>(reader: &mut R) -> Result<()> {
        let count = read_u8(reader)?;
        if count != 0 {
            return Err(Error::UnsupportedAttributes(count));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::item::*;
    use super::*;
    use crate::metadata::Value;
    use std::io::Cursor;

    const ALL_KINDS: [Kind; 5] = [
        Kind::Mesh3,
        Kind::Polyline2,
        Kind::Polyline3,
        Kind::Cloud2,
        Kind::Cloud3,
    ];

    #[test]
    fn header_round_trips_for_every_kind() {
        for kind in ALL_KINDS {
            for count in [0u32, 1, 12, u32::MAX] {
                let mut buf = Vec::new();
                write_header(&mut buf, kind, count, &Metadata::new()).unwrap();
                assert_eq!(buf.len(), 12, "{kind:?}");

                let header = read_header(&mut Cursor::new(&buf), kind).unwrap();
                assert_eq!(header.version, VERSION);
                assert_eq!(header.kind, kind);
                assert_eq!(header.count, count);
            }
        }
    }

    #[test]
    fn kind_bytes_are_stable() {
        // These are on-disk values, not an implementation detail. Renumbering them silently
        // reinterprets every existing file.
        assert_eq!(Kind::Mesh3.as_byte(), 1);
        assert_eq!(Kind::Polyline2.as_byte(), 2);
        assert_eq!(Kind::Polyline3.as_byte(), 3);
        assert_eq!(Kind::Cloud2.as_byte(), 4);
        assert_eq!(Kind::Cloud3.as_byte(), 5);

        for kind in ALL_KINDS {
            assert_eq!(Kind::from_byte(kind.as_byte()).unwrap(), kind);
        }
    }

    #[test]
    fn kinds_report_their_dimension_and_extension() {
        assert_eq!(Kind::Mesh3.dimension(), 3);
        assert_eq!(Kind::Polyline2.dimension(), 2);
        assert_eq!(Kind::Cloud2.dimension(), 2);
        assert_eq!(Kind::Polyline3.dimension(), 3);
        assert_eq!(Kind::Cloud3.dimension(), 3);

        assert_eq!(Kind::Mesh3.extension(), "tcmesh");
        assert_eq!(Kind::Polyline2.extension(), "tccurve2");
        assert_eq!(Kind::Polyline3.extension(), "tccurve3");
    }

    #[test]
    fn probe_identifies_a_file_without_being_told_the_kind() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Polyline2, 7, &Metadata::new()).unwrap();

        let header = probe(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(header.kind, Kind::Polyline2);
        assert_eq!(header.count, 7);
    }

    #[test]
    fn wrong_magic_is_rejected() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        buf[0] = b'X';

        assert!(matches!(
            probe(&mut Cursor::new(&buf)),
            Err(Error::BadMagic)
        ));
    }

    #[test]
    fn a_future_version_is_rejected() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        buf[4] = VERSION + 1;

        assert!(matches!(
            probe(&mut Cursor::new(&buf)),
            Err(Error::UnsupportedVersion(v)) if v == VERSION + 1
        ));
    }

    #[test]
    fn an_unknown_kind_is_rejected() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        buf[5] = 99;

        assert!(matches!(
            probe(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    /// The reserved byte must fail loudly rather than being ignored, otherwise a future compressed
    /// file would be read as an uncompressed one and produce garbage geometry.
    #[test]
    fn a_declared_compression_is_refused_rather_than_ignored() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        buf[6] = 1;

        assert!(matches!(
            probe(&mut Cursor::new(&buf)),
            Err(Error::UnsupportedCompression(1))
        ));
    }

    #[test]
    fn asking_for_the_wrong_kind_is_rejected() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Cloud3, 1, &Metadata::new()).unwrap();

        assert!(matches!(
            read_header(&mut Cursor::new(&buf), Kind::Mesh3),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn truncation_at_every_offset_is_an_error() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 3, &Metadata::new()).unwrap();

        for cut in 0..buf.len() {
            assert!(
                probe(&mut Cursor::new(&buf[..cut])).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
        assert!(probe(&mut Cursor::new(&buf)).is_ok());
    }

    #[test]
    fn names_round_trip() {
        for name in [
            None,
            Some(""),
            Some("section"),
            Some("coupe transversale n\u{b0}3"),
            Some("\u{65ad}\u{9762}"),
        ] {
            let mut buf = Vec::new();
            let f = if name.is_some() { HAS_NAME } else { 0 };
            write_flags(&mut buf, f).unwrap();
            write_name(&mut buf, name).unwrap();

            let mut cursor = Cursor::new(&buf);
            let read_f = read_flags(&mut cursor).unwrap();
            let got = read_name(&mut cursor, read_f).unwrap();

            assert_eq!(got.as_deref(), name, "name {name:?}");
            assert_eq!(cursor.position() as usize, buf.len());
        }
    }

    /// An empty name is a name, distinct from having none. Conflating them would make the
    /// `Option` in the API a lie.
    #[test]
    fn an_empty_name_is_not_the_same_as_no_name() {
        let mut with = Vec::new();
        write_flags(&mut with, HAS_NAME).unwrap();
        write_name(&mut with, Some("")).unwrap();

        let mut without = Vec::new();
        write_flags(&mut without, 0).unwrap();
        write_name(&mut without, None).unwrap();

        assert_ne!(with, without);
        assert_eq!(with.len(), 5);
        assert_eq!(without.len(), 1);
    }

    #[test]
    fn unknown_flag_bits_are_rejected() {
        let mut buf = Vec::new();
        write_flags(&mut buf, 1 << 7).unwrap();

        assert!(matches!(
            read_flags(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn an_invalid_utf8_name_is_rejected() {
        let mut buf = Vec::new();
        write_flags(&mut buf, HAS_NAME).unwrap();
        write_u32(&mut buf, 2).unwrap();
        buf.extend_from_slice(&[0xFF, 0xFE]);

        let mut cursor = Cursor::new(&buf);
        let f = read_flags(&mut cursor).unwrap();
        assert!(matches!(
            read_name(&mut cursor, f),
            Err(Error::Malformed(_))
        ));
    }

    /// A `u32` length field is attacker controlled. Claiming four gigabytes must not allocate.
    #[test]
    fn an_absurd_name_length_does_not_allocate() {
        let mut buf = Vec::new();
        write_flags(&mut buf, HAS_NAME).unwrap();
        write_u32(&mut buf, u32::MAX).unwrap();

        let mut cursor = Cursor::new(&buf);
        let f = read_flags(&mut cursor).unwrap();
        assert!(matches!(
            read_name(&mut cursor, f),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn file_level_metadata_round_trips() {
        let mut meta = Metadata::new();
        meta.insert("units".into(), "mm".into());
        meta.insert("machine.id".into(), Value::I64(42));

        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 3, &meta).unwrap();

        let header = probe(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(header.metadata, meta);
        assert_eq!(header.count, 3);
    }

    /// Absent metadata must cost nothing, so a caller who does not use it pays no bytes.
    #[test]
    fn absent_file_metadata_costs_nothing() {
        let mut bare = Vec::new();
        write_header(&mut bare, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        assert_eq!(bare.len(), 12);

        let mut meta = Metadata::new();
        meta.insert("k".into(), Value::Bool(true));
        let mut tagged = Vec::new();
        write_header(&mut tagged, Kind::Mesh3, 1, &meta).unwrap();
        assert!(tagged.len() > bare.len());
    }

    #[test]
    fn unknown_file_flag_bits_are_rejected() {
        let mut buf = Vec::new();
        write_header(&mut buf, Kind::Mesh3, 1, &Metadata::new()).unwrap();
        buf[7] = 1 << 6;

        assert!(matches!(
            probe(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn an_item_preamble_round_trips() {
        let mut meta = Metadata::new();
        meta.insert("engeom.chord_tol".into(), Value::F64(1e-4));

        for (name, closed) in [
            (None, false),
            (Some("outer"), true),
            (Some(""), false),
            (None, true),
        ] {
            let mut buf = Vec::new();
            write_preamble(&mut buf, name, &meta, closed).unwrap();

            let mut cursor = Cursor::new(&buf);
            let back = read_preamble(&mut cursor, true).unwrap();

            assert_eq!(back.name.as_deref(), name);
            assert_eq!(back.closed, closed);
            assert_eq!(back.metadata, meta);
            assert_eq!(cursor.position() as usize, buf.len());
        }
    }

    /// Kinds with no notion of closure reject the flag rather than ignoring it, so a mesh cannot
    /// quietly carry a bit that means nothing for it.
    #[test]
    fn a_closed_flag_is_refused_where_it_makes_no_sense() {
        let mut buf = Vec::new();
        write_preamble(&mut buf, None, &Metadata::new(), true).unwrap();

        assert!(read_preamble(&mut Cursor::new(&buf), true).is_ok());
        assert!(matches!(
            read_preamble(&mut Cursor::new(&buf), false),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn the_attribute_count_round_trips_as_zero() {
        let mut buf = Vec::new();
        write_attribute_count(&mut buf).unwrap();

        assert_eq!(buf, vec![0]);
        assert!(read_attribute_count(&mut Cursor::new(&buf)).is_ok());
    }

    /// Attribute blocks have no defined layout yet, so a decoder cannot skip past them. Declaring
    /// any must fail rather than leaving the stream misaligned for whatever follows.
    #[test]
    fn declared_attributes_are_refused() {
        assert!(matches!(
            read_attribute_count(&mut Cursor::new([2u8])),
            Err(Error::UnsupportedAttributes(2))
        ));
    }
}
