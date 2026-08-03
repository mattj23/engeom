//! The crate's error type.
//!
//! Deliberately hand-written rather than derived, because a derive macro would mean a dependency
//! and this crate carries none. The type implements [`std::error::Error`], so callers using a
//! boxed-error alias absorb it through `?` with no glue of their own.

use std::fmt;

/// Result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// An error produced while encoding or decoding tolerance-compressed data.
#[non_exhaustive]
#[derive(Debug)]
pub enum Error {
    /// An underlying read or write failed.
    Io(std::io::Error),

    /// The leading bytes of a stream did not identify the expected format.
    BadMagic,

    /// The stream identifies a format version this build does not know how to decode.
    UnsupportedVersion(u8),

    /// The stream declares a compression method this build does not implement.
    ///
    /// The header byte is reserved but always written as 0; nothing compresses tol-compress output
    /// usefully, because quantizing to a tolerance is itself the compression.
    UnsupportedCompression(u8),

    /// An item declares attribute blocks, which this version cannot decode.
    ///
    /// The count is reserved in the layout so attributes can be added without a version bump.
    UnsupportedAttributes(u8),

    /// The stream was structurally invalid. The payload names what was wrong.
    Malformed(&'static str),

    /// The requested tolerance cannot be guaranteed over the requested range.
    ///
    /// Representing `range` to within `tol` would need more bits than an `f64` mantissa can hold,
    /// so encoding would silently return values outside the tolerance. Since that guarantee is the
    /// entire point of the format, this is an error rather than a clamp.
    ToleranceNotRepresentable {
        /// The span of the data along the offending axis.
        range: f64,
        /// The per-axis tolerance that could not be met.
        tol: f64,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "io error: {e}"),
            Error::BadMagic => write!(f, "not a tol-compress stream: invalid magic bytes"),
            Error::UnsupportedVersion(v) => {
                write!(f, "unsupported tol-compress format version {v}")
            }
            Error::UnsupportedCompression(c) => write!(
                f,
                "stream declares compression method {c}, which this build does not support"
            ),
            Error::UnsupportedAttributes(n) => write!(
                f,
                "item declares {n} attribute blocks, which this format version cannot decode"
            ),
            Error::Malformed(what) => write!(f, "malformed tol-compress stream: {what}"),
            Error::ToleranceNotRepresentable { range, tol } => write!(
                f,
                "cannot guarantee a tolerance of {tol} over a range of {range}: \
                 the ratio exceeds what an f64 can represent exactly"
            ),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(value: std::io::Error) -> Self {
        Error::Io(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn io_error_converts_and_chains() {
        let io = std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "ran out");
        let err: Error = io.into();

        assert!(matches!(err, Error::Io(_)));
        assert!(std::error::Error::source(&err).is_some());
        assert!(err.to_string().contains("ran out"));
    }

    #[test]
    fn messages_name_the_problem() {
        assert!(Error::BadMagic.to_string().contains("magic"));
        assert!(
            Error::UnsupportedVersion(7)
                .to_string()
                .contains("version 7")
        );
        assert!(
            Error::Malformed("index out of range")
                .to_string()
                .contains("index out of range")
        );

        // Both magnitudes have to appear so the message says which request was impossible, not
        // merely that something was.
        let (range, tol) = (250.0f64, 1e-20f64);
        let msg = Error::ToleranceNotRepresentable { range, tol }.to_string();
        assert!(msg.contains(&range.to_string()), "{msg}");
        assert!(msg.contains(&tol.to_string()), "{msg}");
    }

    /// The point of implementing `std::error::Error` is that a caller with a boxed-error alias,
    /// which is what engeom uses, can absorb ours through `?` without any conversion glue.
    #[test]
    fn absorbs_into_a_boxed_error_alias() {
        fn boxed() -> std::result::Result<(), Box<dyn std::error::Error>> {
            Err(Error::BadMagic)?;
            Ok(())
        }

        assert!(boxed().is_err());
    }
}
