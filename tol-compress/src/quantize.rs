//! Mapping real values onto integer lattices narrow enough to meet a stated tolerance.
//!
//! This is where the format's guarantee is actually decided. Everything else is packing.

use crate::error::{Error, Result};

/// The largest lattice width the crate will use.
///
/// An `f64` has a 53-bit significand, so `max_int` above this stops being exactly representable
/// and the encode/decode arithmetic starts introducing error of its own.
const MAX_BITS: u8 = 53;

/// How much finer than the requested tolerance the underlying `f64` resolution must be.
///
/// Quantization error is analyzed assuming exact arithmetic; representation error is what is left
/// over. Requiring the latter to be two orders of magnitude below the tolerance keeps it out of
/// the error budget entirely rather than leaving the guarantee resting on it. In practice this is
/// never the binding constraint: it permits a range-to-tolerance ratio around 4.5e13, which is a
/// kilometre-long part measured to a fraction of an angstrom.
const REPRESENTATION_MARGIN: f64 = 100.0;

/// The number of bits needed to represent `range` such that any value within it round-trips to
/// within `tol`.
///
/// Returns 0 when `range` is 0: a dimension whose values are all identical costs nothing to store.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if `range` and `tol` are not both finite and sensibly
/// signed, or if meeting `tol` would need a lattice finer than an `f64` can carry.
pub fn bits_for_tol(range: f64, tol: f64) -> Result<u8> {
    let fail = || Error::ToleranceNotRepresentable { range, tol };

    if !range.is_finite() || range < 0.0 || !tol.is_finite() || tol <= 0.0 {
        return Err(fail());
    }
    if range == 0.0 {
        return Ok(0);
    }

    // Representation error near the top of the range is about `range * 2^-53`. Reject up front if
    // that is not comfortably below the tolerance, rather than issuing a width that cannot deliver.
    if range * (f64::EPSILON / 2.0) * REPRESENTATION_MARGIN > tol {
        return Err(fail());
    }

    // Rounding to the nearest lattice point puts the worst case at a midpoint, half a spacing
    // away, so the requirement is `range / max_int <= 2 * tol`. Stepping the width one bit at a
    // time with an exact comparison avoids the off-by-one a `log2` would risk, and a width that
    // came out one bit low would silently break the guarantee.
    for bits in 1..=MAX_BITS {
        let max_int = ((1u64 << bits) - 1) as f64;
        if range <= 2.0 * tol * max_int {
            return Ok(bits);
        }
    }

    Err(fail())
}

/// Maps values in a fixed interval to and from integer codes of a fixed width.
///
/// Holding the interval and the width together is deliberate: an encoder and decoder that disagree
/// on either would produce silently wrong positions rather than an error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quantizer {
    min: f64,
    range: f64,
    bits: u8,
    /// The largest code, `2^bits - 1`, or 0 when `bits` is 0.
    max_int: u64,
}

impl Quantizer {
    /// Build a quantizer over `[min, min + range]` at an explicit width.
    ///
    /// This is the escape hatch for re-encoding data that is already on a known lattice, where
    /// deriving the width from a tolerance would re-quantize it. Use [`Quantizer::from_tol`] when
    /// encoding fresh data.
    ///
    /// # Panics
    ///
    /// Panics if `bits` exceeds 53.
    pub fn new(min: f64, range: f64, bits: u8) -> Self {
        assert!(bits <= MAX_BITS, "bit width {bits} exceeds {MAX_BITS}");
        Self {
            min,
            range,
            bits,
            max_int: if bits == 0 { 0 } else { (1u64 << bits) - 1 },
        }
    }

    /// Build the narrowest quantizer over `[min, min + range]` that meets `tol`.
    ///
    /// # Errors
    ///
    /// Propagates [`bits_for_tol`].
    pub fn from_tol(min: f64, range: f64, tol: f64) -> Result<Self> {
        Ok(Self::new(min, range, bits_for_tol(range, tol)?))
    }

    /// The lattice width in bits.
    pub fn bits(&self) -> u8 {
        self.bits
    }

    /// The largest code this quantizer emits.
    pub fn max_int(&self) -> u64 {
        self.max_int
    }

    /// Map a value to its nearest lattice code.
    ///
    /// Values outside the interval clamp to its ends rather than producing a code too wide for the
    /// declared width, which would corrupt the stream.
    pub fn encode(&self, value: f64) -> u64 {
        if self.bits == 0 || self.range == 0.0 {
            return 0;
        }
        let fraction = ((value - self.min) / self.range).clamp(0.0, 1.0);
        (fraction * self.max_int as f64).round() as u64
    }

    /// Map a lattice code back to a value.
    pub fn decode(&self, code: u64) -> f64 {
        if self.bits == 0 || self.range == 0.0 {
            return self.min;
        }
        self.min + (code as f64 / self.max_int as f64) * self.range
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_range_costs_nothing() {
        assert_eq!(bits_for_tol(0.0, 1e-6).unwrap(), 0);

        // And it still round-trips, which is what makes a degenerate axis free rather than broken.
        let q = Quantizer::new(7.5, 0.0, 0);
        assert_eq!(q.encode(7.5), 0);
        assert_eq!(q.decode(0), 7.5);
    }

    #[test]
    fn widths_match_a_hand_computed_table() {
        // With round-to-nearest the requirement is range <= 2 * tol * (2^bits - 1).
        // tol = 0.5 makes the arithmetic legible: range <= 2^bits - 1.
        assert_eq!(bits_for_tol(1.0, 0.5).unwrap(), 1); // 1 <= 1
        assert_eq!(bits_for_tol(2.0, 0.5).unwrap(), 2); // 2 <= 3
        assert_eq!(bits_for_tol(3.0, 0.5).unwrap(), 2); // 3 <= 3
        assert_eq!(bits_for_tol(4.0, 0.5).unwrap(), 3); // 4 <= 7
        assert_eq!(bits_for_tol(255.0, 0.5).unwrap(), 8);
        assert_eq!(bits_for_tol(256.0, 0.5).unwrap(), 9);
    }

    /// The case the whole format exists for: a 100 mm part at 1 micron.
    #[test]
    fn a_realistic_metrology_case_is_narrow() {
        let bits = bits_for_tol(100.0, 1e-3).unwrap();
        assert_eq!(bits, 16);
        // Which is to say a 3D point costs 6 bytes, where an f64 triple would cost 24.
    }

    #[test]
    fn rejects_nonsense_inputs() {
        assert!(bits_for_tol(-1.0, 1e-6).is_err());
        assert!(bits_for_tol(f64::NAN, 1e-6).is_err());
        assert!(bits_for_tol(f64::INFINITY, 1e-6).is_err());
        assert!(bits_for_tol(1.0, 0.0).is_err());
        assert!(bits_for_tol(1.0, -1e-6).is_err());
        assert!(bits_for_tol(1.0, f64::NAN).is_err());
    }

    #[test]
    fn rejects_a_tolerance_finer_than_f64_can_carry() {
        // A 1 metre range asked to hold 1e-18 m needs more resolution than an f64 has.
        assert!(bits_for_tol(1.0, 1e-18).is_err());

        // The boundary sits around a ratio of 4.5e13, so either side of it should behave.
        assert!(bits_for_tol(1.0, 1e-12).is_ok());
        assert!(bits_for_tol(1.0, 1e-15).is_err());
    }

    /// The load-bearing test. Sweeps ranges and tolerances across many orders of magnitude and
    /// checks the worst case over values placed where quantization error peaks: exactly halfway
    /// between lattice points. If the halved bound were wrong, this is what would catch it.
    #[test]
    fn round_trip_error_never_exceeds_tolerance() {
        let ranges = [
            1e-6, 1e-3, 0.5, 1.0, 2.5, 100.0, 1234.567, 1e6, 1e9, 6.02e12,
        ];
        let tols = [1e-9, 1e-6, 1e-4, 1e-3, 0.01, 0.5, 1.0];
        let mins = [0.0, -1.0, 1e5, -3.7e8];

        let mut worst_ratio = 0.0f64;
        let mut checked = 0usize;

        for &range in &ranges {
            for &tol in &tols {
                let Ok(bits) = bits_for_tol(range, tol) else {
                    continue;
                };
                for &min in &mins {
                    let q = Quantizer::new(min, range, bits);
                    let max_int = q.max_int().max(1);
                    let spacing = range / max_int as f64;

                    // Lattice midpoints are the worst case, plus both endpoints, plus a spread of
                    // interior samples that are not aligned to the lattice at all.
                    let mut samples: Vec<f64> = vec![min, min + range];
                    for k in 0..64u64 {
                        let idx = (k * max_int / 64).min(max_int.saturating_sub(1));
                        samples.push(min + (idx as f64 + 0.5) * spacing);
                    }
                    for k in 0..64 {
                        samples.push(min + range * (k as f64 * 0.987_654_321 % 1.0));
                    }

                    for v in samples {
                        let back = q.decode(q.encode(v));
                        let err = (back - v).abs();
                        assert!(
                            err <= tol,
                            "range {range}, tol {tol}, bits {bits}, min {min}: \
                             value {v} recovered as {back}, error {err}"
                        );
                        worst_ratio = worst_ratio.max(err / tol);
                        checked += 1;
                    }
                }
            }
        }

        assert!(checked > 10_000, "sweep did not actually cover much");
        // Should approach but not reach 1.0. Far below would mean widths are wastefully wide.
        assert!(
            worst_ratio > 0.25,
            "worst observed error was only {worst_ratio} of tolerance, \
             which suggests the widths are wider than necessary"
        );
    }

    /// One bit narrower must actually fail, which is what proves the chosen width is minimal
    /// rather than merely sufficient.
    #[test]
    fn the_chosen_width_is_minimal() {
        let cases = [(100.0, 1e-3), (1.0, 1e-6), (2.5, 0.01), (1e6, 1.0)];

        for (range, tol) in cases {
            let bits = bits_for_tol(range, tol).unwrap();
            assert!(bits > 0);

            let narrower = Quantizer::new(0.0, range, bits - 1);
            let max_int = narrower.max_int().max(1);
            let spacing = range / max_int as f64;

            let midpoint = spacing * 0.5;
            let err = (narrower.decode(narrower.encode(midpoint)) - midpoint).abs();
            assert!(
                err > tol,
                "range {range} at tol {tol} claims to need {bits} bits, \
                 but {} bits achieved an error of only {err}",
                bits - 1
            );
        }
    }

    #[test]
    fn codes_stay_inside_the_declared_width() {
        let q = Quantizer::from_tol(10.0, 4.0, 0.01).unwrap();
        let limit = q.max_int();

        // Including values well outside the interval, which must clamp rather than overflow the
        // width and corrupt the bitstream.
        for v in [-1e9, 9.9, 10.0, 12.0, 14.0, 14.1, 1e9] {
            assert!(
                q.encode(v) <= limit,
                "value {v} produced an out-of-range code"
            );
        }
    }

    #[test]
    fn endpoints_are_exact() {
        let q = Quantizer::from_tol(-5.0, 12.0, 1e-4).unwrap();
        assert_eq!(q.encode(-5.0), 0);
        assert_eq!(q.decode(0), -5.0);
        assert_eq!(q.encode(7.0), q.max_int());
        assert!((q.decode(q.max_int()) - 7.0).abs() < 1e-12);
    }

    /// `Quantizer::new` has to reproduce a lattice exactly, since increment 4a depends on it to
    /// re-encode existing files without shifting any value.
    #[test]
    fn explicit_widths_reproduce_a_lattice_exactly() {
        let q = Quantizer::new(2.0, 8.0, 20);
        for code in [0u64, 1, 512, 65_535, q.max_int() - 1, q.max_int()] {
            assert_eq!(
                q.encode(q.decode(code)),
                code,
                "code {code} did not survive"
            );
        }
    }

    #[test]
    #[should_panic(expected = "exceeds 53")]
    fn width_above_the_cap_panics() {
        Quantizer::new(0.0, 1.0, 54);
    }
}
