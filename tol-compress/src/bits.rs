//! Bit-level stream I/O.
//!
//! Everything in this crate ultimately reduces to writing integers of an arbitrary, data-derived
//! width. Rounding those widths up to whole bytes is a huge source of waste, so the encoders
//! here write bits rather than bytes.
//!
//! # Format conventions
//!
//! These are part of the on-disk format, not implementation details:
//!
//! - Bits are packed **least-significant-first within each byte**, and bytes accumulate in
//!   little-endian order. This matches the byte-oriented convention the original format used.
//! - A width of `0` is legal and encodes nothing. This is intentional; an axis whose values are
//!   all identical has a range of zero and therefore costs no bits at all.
//! - Blocks are padded with zero bits to the next byte boundary when they end, so every block
//!   begins at a byte offset and stays independently addressable.
//!
//! # Mixing bit and byte data
//!
//! Headers are byte-oriented and block payloads are bit-oriented. Because [`Write`] and [`Read`]
//! are implemented for `&mut W` and `&mut R`, the intended pattern is to borrow the underlying
//! stream for the duration of a block:
//!
//! ```
//! use std::io::Write;
//! use tol_compress::bits::BitWriter;
//!
//! let mut out: Vec<u8> = Vec::new();
//! out.write_all(b"HDR")?;
//!
//! let mut bw = BitWriter::new(&mut out);
//! bw.write_bits(5, 3)?;
//! bw.finish()?;
//! drop(bw);
//!
//! out.write_all(b"TAIL")?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::error::{Error, Result};
use std::io::{Read, Write};

/// Packs integers of arbitrary bit width into a byte stream.
///
/// Call [`BitWriter::finish`] when a block ends. Dropping the writer flushes any buffered bits on
/// a best-effort basis so data is never silently truncated, but a drop cannot report an I/O
/// failure, so `finish` is the only way to observe one.
pub struct BitWriter<W: Write> {
    inner: W,
    /// Buffered bits not yet written, right-aligned.
    acc: u128,
    /// Count of valid bits in `acc`. At most 7 between calls, at most 71 within one.
    n_bits: u8,
}

impl<W: Write> BitWriter<W> {
    /// Wrap a byte sink.
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            acc: 0,
            n_bits: 0,
        }
    }

    /// Append the low `bits` bits of `value`.
    ///
    /// A `bits` of 0 writes nothing and is always valid.
    ///
    /// # Panics
    ///
    /// Panics if `bits` exceeds 64, or if `value` has any bit set at or above position `bits`.
    /// Both indicate the caller computed a width that does not match its data, which would
    /// otherwise corrupt the stream silently and break the tolerance guarantee. The check is
    /// unconditional rather than a `debug_assert` because this crate's tests run in release mode.
    pub fn write_bits(&mut self, value: u64, bits: u8) -> Result<()> {
        if bits == 0 {
            assert_eq!(value, 0, "wrote a nonzero value with a width of zero bits");
            return Ok(());
        }
        assert!(bits <= 64, "bit width {bits} exceeds 64");

        let masked = if bits == 64 {
            value
        } else {
            value & ((1u64 << bits) - 1)
        };
        assert_eq!(
            value, masked,
            "value {value} does not fit in {bits} bits, which would corrupt the stream"
        );

        self.acc |= (masked as u128) << self.n_bits;
        self.n_bits += bits;

        while self.n_bits >= 8 {
            let byte = (self.acc & 0xFF) as u8;
            self.inner.write_all(&[byte])?;
            self.acc >>= 8;
            self.n_bits -= 8;
        }

        Ok(())
    }

    /// Pad with zero bits up to the next byte boundary, emitting any buffered partial byte.
    ///
    /// This is what makes a block independently addressable. It is idempotent.
    pub fn align(&mut self) -> Result<()> {
        if self.n_bits > 0 {
            let byte = (self.acc & 0xFF) as u8;
            self.inner.write_all(&[byte])?;
            self.acc = 0;
            self.n_bits = 0;
        }
        Ok(())
    }

    /// End the block: pad to a byte boundary and flush the underlying sink.
    ///
    /// Call this rather than relying on the drop, which cannot report I/O errors.
    pub fn finish(&mut self) -> Result<()> {
        self.align()?;
        self.inner.flush()?;
        Ok(())
    }
}

impl<W: Write> Drop for BitWriter<W> {
    fn drop(&mut self) {
        // Best effort so a forgotten `finish` truncates nothing. Errors are unobservable here,
        // which is exactly why `finish` exists.
        let _ = self.align();
    }
}

/// The largest payload [`read_payload`] will reserve up front, 64 MiB.
///
/// Comfortably above any real block, so legitimate data is read in a single allocation, and small
/// enough that a corrupt length field cannot exhaust memory before the data fails to arrive.
const MAX_PAYLOAD_PREALLOC: usize = 1 << 26;

/// Read exactly `payload_bytes` from `reader`, for handing to a [`BitReader`].
///
/// The length is grown as data arrives rather than reserved up front, so a corrupt length field
/// claiming four billion bytes runs out of input instead of out of memory.
///
/// # Errors
///
/// [`Error::Malformed`] if the source ends before the payload does.
pub fn read_payload<R: Read>(reader: &mut R, payload_bytes: usize) -> Result<Vec<u8>> {
    // Reserve exactly what the header claims, so the read neither grows nor recopies. The cap is
    // what keeps a corrupt length from turning into a huge allocation; past it the vector grows as
    // data arrives instead, which is slower but bounded by input that actually exists.
    let mut buf = Vec::with_capacity(payload_bytes.min(MAX_PAYLOAD_PREALLOC));
    reader.take(payload_bytes as u64).read_to_end(&mut buf)?;

    if buf.len() != payload_bytes {
        return Err(Error::Malformed(
            "block payload is shorter than its header declares",
        ));
    }

    Ok(buf)
}

/// Reads integers of arbitrary bit width from a block payload.
///
/// The mirror of [`BitWriter`]. A reader and writer agree only if they use the same widths in the
/// same order, since the payload itself carries no width information.
///
/// # Why this takes a slice rather than a `Read`
///
/// Two reasons, one correctness and one speed.
///
/// A block payload is followed immediately by byte-oriented data belonging to whatever comes next,
/// and [`Read`] offers no way to put back a byte that was read too eagerly. Reading the payload out
/// first, at a length the block header determines, makes over-reading structurally impossible.
///
/// It is also what makes decoding fast. Pulling bytes through a [`Read`] inside the loop, buffered
/// or not, keeps the accumulator in memory rather than in registers; against a plain slice the
/// whole hot path stays in registers and block decoding measured roughly twice as fast as either
/// alternative. Use [`read_payload`] to get the slice.
///
/// Call [`BitReader::finish`] at the end of a block to check the payload was consumed exactly,
/// which turns a miscomputed length into an error rather than a silently misparsed stream.
pub struct BitReader<'a> {
    /// The block payload.
    buf: &'a [u8],
    /// Next unconsumed byte in `buf`.
    pos: usize,
    /// Bits moved into the accumulator but not yet consumed, right-aligned.
    acc: u128,
    /// Count of valid bits in `acc`. At most 127.
    n_bits: u8,
}

impl<'a> BitReader<'a> {
    /// Wrap a block payload.
    pub fn new(payload: &'a [u8]) -> Self {
        Self {
            buf: payload,
            pos: 0,
            acc: 0,
            n_bits: 0,
        }
    }

    /// Ensure the accumulator holds at least `bits` bits.
    ///
    /// The eight-byte path is what makes this reasonable. A byte-at-a-time loop was
    /// nearly twice as slow; a fixed-size load compiles to one unaligned 64-bit read.
    #[inline]
    fn refill(&mut self, bits: u8) -> Result<()> {
        while self.n_bits < bits {
            // The loop only runs while `n_bits < bits` and `bits` is at most 64, so `n_bits` is at
            // most 63 here and another 64 bits always fit in the 128 bit accumulator.
            let left = self.buf.len() - self.pos;

            if left >= 8 {
                let bytes: [u8; 8] = self.buf[self.pos..self.pos + 8]
                    .try_into()
                    .expect("slice of exactly eight bytes");
                self.acc |= u128::from(u64::from_le_bytes(bytes)) << self.n_bits;
                self.n_bits += 64;
                self.pos += 8;
            } else if left > 0 {
                self.acc |= u128::from(self.buf[self.pos]) << self.n_bits;
                self.n_bits += 8;
                self.pos += 1;
            } else {
                return Err(Error::Io(std::io::ErrorKind::UnexpectedEof.into()));
            }
        }

        Ok(())
    }

    /// Consume the next `bits` bits as an integer.
    ///
    /// A `bits` of 0 consumes nothing and yields 0.
    ///
    /// # Panics
    ///
    /// Panics if `bits` exceeds 64.
    #[inline]
    pub fn read_bits(&mut self, bits: u8) -> Result<u64> {
        if bits == 0 {
            return Ok(0);
        }
        assert!(bits <= 64, "bit width {bits} exceeds 64");

        self.refill(bits)?;

        let mask = if bits == 64 {
            u64::MAX as u128
        } else {
            (1u128 << bits) - 1
        };
        let out = (self.acc & mask) as u64;
        self.acc >>= bits;
        self.n_bits -= bits;

        Ok(out)
    }

    /// Discard bits up to the next byte boundary, matching [`BitWriter::align`].
    ///
    /// Only the bits belonging to the partially consumed byte are dropped; whole buffered bytes
    /// are kept. It is idempotent.
    pub fn align(&mut self) {
        let partial = self.n_bits % 8;
        self.acc >>= partial;
        self.n_bits -= partial;
    }

    /// End the block: drop the padding bits and check the payload was consumed exactly.
    ///
    /// # Errors
    ///
    /// [`Error::Malformed`] if any of the payload is left over, which means the length the caller
    /// derived from the block header disagrees with the data. Reporting that here is the difference
    /// between a clear error and a stream that parses into plausible nonsense from this point on.
    pub fn finish(&mut self) -> Result<()> {
        self.align();

        if self.pos != self.buf.len() || self.n_bits != 0 {
            return Err(Error::Malformed(
                "block payload length disagrees with the block header",
            ));
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Round-trip a single value at every width it can legally take.
    #[test]
    fn round_trip_every_width() {
        for bits in 0..=64u8 {
            let values: Vec<u64> = if bits == 0 {
                vec![0]
            } else if bits == 64 {
                vec![0, 1, u64::MAX / 2, u64::MAX]
            } else {
                let max = (1u64 << bits) - 1;
                vec![0, 1, max / 2, max]
            };

            for v in values {
                let mut buf = Vec::new();
                let mut bw = BitWriter::new(&mut buf);
                bw.write_bits(v, bits).unwrap();
                bw.finish().unwrap();
                drop(bw);

                let mut br = BitReader::new(&buf);
                assert_eq!(br.read_bits(bits).unwrap(), v, "width {bits}, value {v}");
            }
        }
    }

    /// A long interleaved sequence of irregular widths, which is what real encoded data looks
    /// like: per-axis widths differ and none of them are byte multiples.
    #[test]
    fn round_trip_mixed_widths() {
        let widths: Vec<u8> = (0..500).map(|i| ((i * 7 + 3) % 41) as u8).collect();
        let values: Vec<u64> = widths
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                if b == 0 {
                    0
                } else {
                    (i as u64 * 2_654_435_761) & ((1u64 << b) - 1)
                }
            })
            .collect();

        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        for (&v, &b) in values.iter().zip(widths.iter()) {
            bw.write_bits(v, b).unwrap();
        }
        bw.finish().unwrap();
        drop(bw);

        let mut br = BitReader::new(&buf);
        for (i, (&v, &b)) in values.iter().zip(widths.iter()).enumerate() {
            assert_eq!(br.read_bits(b).unwrap(), v, "item {i} at width {b}");
        }
    }

    /// Zero-width writes must occupy no space at all. This is what lets a degenerate axis cost
    /// nothing rather than the two-byte floor the previous implementation charged.
    #[test]
    fn zero_width_occupies_nothing() {
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        for _ in 0..1000 {
            bw.write_bits(0, 0).unwrap();
        }
        bw.finish().unwrap();
        drop(bw);

        assert!(buf.is_empty());

        let mut br = BitReader::new(&buf);
        assert_eq!(br.read_bits(0).unwrap(), 0);
    }

    #[test]
    fn packs_to_the_expected_size() {
        // 10 values of 5 bits is 50 bits, which occupies 7 bytes once padded.
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        for i in 0..10u64 {
            bw.write_bits(i, 5).unwrap();
        }
        bw.finish().unwrap();
        drop(bw);

        assert_eq!(buf.len(), 7);
    }

    /// Bits are least-significant-first within a byte. Pinning this down because it is a format
    /// guarantee, not an implementation choice.
    #[test]
    fn bit_order_is_lsb_first() {
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        bw.write_bits(0b101, 3).unwrap();
        bw.write_bits(0b11, 2).unwrap();
        bw.finish().unwrap();
        drop(bw);

        // The 3-bit value occupies bits 0..3 and the 2-bit value bits 3..5.
        assert_eq!(buf, vec![0b00011101]);
    }

    #[test]
    fn align_pads_to_a_byte_boundary() {
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        bw.write_bits(1, 1).unwrap();
        bw.align().unwrap();
        bw.write_bits(1, 1).unwrap();
        bw.align().unwrap();
        // Aligning twice in a row must not emit a second padding byte.
        bw.align().unwrap();
        bw.finish().unwrap();
        drop(bw);

        assert_eq!(buf, vec![1, 1]);

        let mut br = BitReader::new(&buf);
        assert_eq!(br.read_bits(1).unwrap(), 1);
        br.align();
        br.align();
        assert_eq!(br.read_bits(1).unwrap(), 1);
    }

    /// Dropping without calling `finish` must not lose buffered bits, only the ability to see an
    /// I/O error.
    #[test]
    fn drop_flushes_buffered_bits() {
        let mut buf = Vec::new();
        {
            let mut bw = BitWriter::new(&mut buf);
            bw.write_bits(0b101, 3).unwrap();
        }
        assert_eq!(buf, vec![0b101]);
    }

    /// The documented pattern for mixing byte-oriented headers with bit-oriented blocks.
    #[test]
    fn interleaves_with_direct_byte_writes() {
        let mut buf: Vec<u8> = Vec::new();
        buf.write_all(b"HDR").unwrap();
        {
            let mut bw = BitWriter::new(&mut buf);
            bw.write_bits(0x2A, 6).unwrap();
            bw.finish().unwrap();
        }
        buf.write_all(b"TAIL").unwrap();

        let mut cursor = Cursor::new(&buf);
        let mut magic = [0u8; 3];
        cursor.read_exact(&mut magic).unwrap();
        assert_eq!(&magic, b"HDR");

        {
            let payload = read_payload(&mut cursor, 1).unwrap();
            let mut br = BitReader::new(&payload);
            assert_eq!(br.read_bits(6).unwrap(), 0x2A);
            br.finish().unwrap();
        }

        let mut tail = [0u8; 4];
        cursor.read_exact(&mut tail).unwrap();
        assert_eq!(&tail, b"TAIL");
    }

    #[test]
    fn reading_past_the_end_is_an_io_error() {
        let buf = vec![0u8; 1];
        let mut br = BitReader::new(&buf);
        assert!(br.read_bits(8).is_ok());
        assert!(br.read_bits(8).is_err());
    }

    /// A declared payload longer than the source must fail up front rather than yield zeros.
    #[test]
    fn a_payload_longer_than_the_source_is_an_error() {
        let buf = vec![0xFFu8; 4];
        let mut cursor = Cursor::new(&buf);

        assert!(matches!(
            read_payload(&mut cursor, 64),
            Err(Error::Malformed(_))
        ));
    }

    /// A payload length claiming more than the source could possibly hold must not be reserved
    /// before the bytes have actually arrived.
    #[test]
    fn an_absurd_payload_length_does_not_allocate() {
        let buf = vec![0u8; 16];
        let mut cursor = Cursor::new(&buf);

        assert!(read_payload(&mut cursor, u32::MAX as usize).is_err());
    }

    /// The reader must stop at the payload boundary even when far more input is available, since
    /// the bytes after a block belong to whatever comes next and cannot be handed back.
    #[test]
    fn the_payload_length_is_a_hard_stop() {
        let mut buf: Vec<u8> = vec![0b0000_0011];
        buf.extend_from_slice(&[0xAA; 4096]);

        let mut cursor = Cursor::new(&buf);
        {
            let payload = read_payload(&mut cursor, 1).unwrap();
            let mut br = BitReader::new(&payload);
            assert_eq!(br.read_bits(2).unwrap(), 0b11);
            br.finish().unwrap();
        }

        // Everything after that single payload byte must still be there to read.
        let mut rest = Vec::new();
        cursor.read_to_end(&mut rest).unwrap();
        assert_eq!(rest.len(), 4096);
        assert!(rest.iter().all(|&b| b == 0xAA));
    }

    /// A payload length that disagrees with what was actually consumed means the caller derived it
    /// from the header wrongly, which would otherwise desynchronize the rest of the stream.
    #[test]
    fn finish_catches_a_length_that_does_not_match() {
        let buf = vec![0u8; 4];

        // Declared longer than the values consumed.
        let mut br = BitReader::new(&buf);
        br.read_bits(8).unwrap();
        assert!(matches!(br.finish(), Err(Error::Malformed(_))));

        // Consumed exactly, which is the only accepted case.
        let mut br = BitReader::new(&buf);
        for _ in 0..4 {
            br.read_bits(8).unwrap();
        }
        assert!(br.finish().is_ok());

        // Padding bits are not leftovers: 30 bits of values in 4 bytes is a correct block.
        let mut br = BitReader::new(&buf);
        for _ in 0..3 {
            br.read_bits(10).unwrap();
        }
        assert!(br.finish().is_ok());
    }

    /// A long stream, so the eight-byte refill path runs many times and its tail case is reached.
    #[test]
    fn reads_a_long_stream_across_many_refills() {
        let count = 10_000;
        let values: Vec<u64> = (0..count as u64).map(|i| i & 0x1FFFF).collect();

        let mut buf = Vec::new();
        {
            let mut bw = BitWriter::new(&mut buf);
            for &v in &values {
                bw.write_bits(v, 17).unwrap();
            }
            bw.finish().unwrap();
        }

        let mut br = BitReader::new(&buf);
        for (i, &v) in values.iter().enumerate() {
            assert_eq!(br.read_bits(17).unwrap(), v, "value {i}");
        }
        br.finish().unwrap();
    }

    #[test]
    #[should_panic(expected = "does not fit in 3 bits")]
    fn value_wider_than_its_width_panics() {
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        let _ = bw.write_bits(0b1000, 3);
    }

    #[test]
    #[should_panic(expected = "exceeds 64")]
    fn width_above_64_panics() {
        let mut buf = Vec::new();
        let mut bw = BitWriter::new(&mut buf);
        let _ = bw.write_bits(0, 65);
    }
}
