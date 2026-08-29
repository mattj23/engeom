//! This module contains common constructs for working with angles

use crate::Result;
use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};

pub const ANGLE_TOL: f64 = 1.0e-12;

/// The minimum number of segments used to approximate a full circle when working with chordal
/// tolerances.
pub const MIN_FULL_CIRCLE_SEGMENTS: usize = 8;

/// Calculate the number of equal segments needed so that a chord inscribed on an arc of `radius`,
/// swept through `sweep` radians, deviates from the true arc by no more than `tol` (the maximum
/// sagitta). The segment endpoints lie on the arc, so the chords sag inward.
///
/// A full circle always uses at least [`MIN_FULL_CIRCLE_SEGMENTS`] segments. Partial sweeps receive
/// a proportional share of that minimum, so even a tolerance loose enough to permit a single
/// chord still produces a recognizable arc.
///
/// # Arguments
///
/// * `radius`: radius of the arc, which must be positive and finite
/// * `sweep`: the swept angle in radians. The sign is ignored and the magnitude is clamped to a
///   full circle.
/// * `tol`: the maximum allowed chordal deviation, which must be positive and finite. A tolerance
///   so tight that the count would exceed a `u32` index range is an error.
///
/// returns: `Result<usize>`
///
/// # Examples
///
/// ```
/// use engeom::common::arc_segments_for_tol;
/// use std::f64::consts::TAU;
///
/// // A chord across 1/16th of a unit circle sags by 1 - cos(pi/16), so using that
/// // tolerance produces 16 segments.
/// let tol = 1.0 - (TAU / 32.0).cos();
/// assert_eq!(arc_segments_for_tol(1.0, TAU, tol).unwrap(), 16);
///
/// // A tolerance looser than the radius falls back to the minimum guard.
/// assert_eq!(arc_segments_for_tol(1.0, TAU, 10.0).unwrap(), 8);
/// ```
pub fn arc_segments_for_tol(radius: f64, sweep: f64, tol: f64) -> Result<usize> {
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(format!("Arc radius must be positive and finite, got {radius}").into());
    }
    if !(tol > 0.0 && tol.is_finite()) {
        return Err(format!("Chordal tolerance must be positive and finite, got {tol}").into());
    }

    let sweep = sweep.abs().min(TAU);

    // The largest arc a single chord can span while its sagitta stays within tol. This is
    // 2 * acos(1 - tol / radius) (the form in `Arc2::to_points`) rewritten through the half-angle
    // identity 1 - cos(a) = 2 * sin^2(a / 2), which stays well conditioned when the tolerance is
    // many orders of magnitude below the radius. The min handles tol >= 2 * radius, where a
    // single chord satisfies the tolerance and the minimum guard decides the count.
    let max_theta = 4.0 * (tol / (2.0 * radius)).sqrt().min(1.0).asin();

    // Segment counts beyond a u32 index range have no physical meaning and exceed what the
    // u32-indexed meshes downstream could hold, so treat them as an unresolvable tolerance. The
    // finiteness check catches a ratio that overflowed to infinity or became NaN.
    let ratio = sweep / max_theta;
    if !ratio.is_finite() || ratio > u32::MAX as f64 {
        return Err(format!(
            "Tolerance {tol} is too small relative to radius {radius} to tessellate"
        )
        .into());
    }

    // The slack prevents a ratio within floating error of an integer, such as a tolerance computed
    // from an exact segment count, from being increased by a whole segment by the ceil. It is orders
    // of magnitude below the significance of any chordal tolerance.
    let n = (ratio * (1.0 - 1.0e-12)).ceil() as usize;
    let floor = ((MIN_FULL_CIRCLE_SEGMENTS as f64) * sweep / TAU).ceil() as usize;
    Ok(n.max(floor).max(1))
}

/// Lists the two possible directions of rotation, clockwise and counter-clockwise.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum AngleDir {
    /// Clockwise rotation (negative angular direction).
    Cw,
    /// Counter-clockwise rotation (positive angular direction).
    Ccw,
}

impl AngleDir {
    /// Returns `-1.0` for `Cw` and `1.0` for `Ccw`, suitable for use as a sign multiplier.
    pub fn to_sign(self) -> f64 {
        match self {
            AngleDir::Cw => -1.0,
            AngleDir::Ccw => 1.0,
        }
    }

    /// Returns `Cw` if `sign < 0.0`, otherwise `Ccw`.
    pub fn from_sign(sign: f64) -> Self {
        if sign < 0.0 {
            AngleDir::Cw
        } else {
            AngleDir::Ccw
        }
    }

    /// Returns the opposite direction: `Cw` <-> `Ccw`.
    pub fn opposite(self) -> Self {
        match self {
            AngleDir::Cw => AngleDir::Ccw,
            AngleDir::Ccw => AngleDir::Cw,
        }
    }

    pub fn is_cw(self) -> bool {
        matches!(self, AngleDir::Cw)
    }

    pub fn is_ccw(self) -> bool {
        matches!(self, AngleDir::Ccw)
    }
}

/// Get the counterclockwise angle from `radians0` to `radians1`
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the destination angle, in radians
///
/// returns: f64
pub fn angle_ccw_to(radians0: f64, radians1: f64) -> f64 {
    let t0 = angle_signed_pi(radians0);
    let t1 = angle_signed_pi(radians1);
    let t1 = if t1 < t0 { t1 + 2.0 * PI } else { t1 };
    t1 - t0
}

/// Get the clockwise angle from `radians0` to `radians1`
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the destination angle, in radians
///
/// returns: f64
pub fn angle_cw_to(radians0: f64, radians1: f64) -> f64 {
    let t0 = angle_signed_pi(radians0);
    let t1 = angle_signed_pi(radians1);
    let t1 = if t1 > t0 { t1 - 2.0 * PI } else { t1 };
    t0 - t1
}

/// Calculates the angle between two angles in a given rotational direction. The angle returned
/// is the angle in the given rotational direction (clockwise or counter-clockwise) which `radians0`
/// would need to be rotated to align with `radians1`. The result will always be positive, in the
/// range [0, 2π].
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the ending angle, in radians
/// * `angle_dir`: the rotational direction to consider
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::{angle_in_direction, AngleDir};
/// use std::f64::consts::{PI, FRAC_PI_2};
/// use approx::assert_relative_eq;
///
/// // Going CCW from 0 to π/2 is a quarter turn
/// let ccw = angle_in_direction(0.0, FRAC_PI_2, AngleDir::Ccw);
/// assert_relative_eq!(ccw, FRAC_PI_2, epsilon = 1.0e-10);
///
/// // Going CW from 0 to π/2 is three-quarter turns (the long way around)
/// let cw = angle_in_direction(0.0, FRAC_PI_2, AngleDir::Cw);
/// assert_relative_eq!(cw, 3.0 * FRAC_PI_2, epsilon = 1.0e-10);
/// ```
pub fn angle_in_direction(radians0: f64, radians1: f64, angle_dir: AngleDir) -> f64 {
    match angle_dir {
        AngleDir::Cw => angle_cw_to(radians0, radians1),
        AngleDir::Ccw => angle_ccw_to(radians0, radians1),
    }
}

/// Returns the signed shortest angular distance from `radians0` to `radians1`.
///
/// A positive result means the shortest path is counter-clockwise; a negative result means the
/// shortest path is clockwise. The magnitude is always in the range [0, π].
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the ending angle, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::shortest_angle_between;
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
/// use approx::assert_relative_eq;
///
/// // CCW from 0 to π/2 is the short way, so result is positive
/// let a = shortest_angle_between(0.0, FRAC_PI_2);
/// assert_relative_eq!(a, FRAC_PI_2, epsilon = 1.0e-10);
///
/// // CW from 0 to -π/4 is the short way, so result is negative
/// let b = shortest_angle_between(0.0, -FRAC_PI_4);
/// assert_relative_eq!(b, -FRAC_PI_4, epsilon = 1.0e-10);
/// ```
pub fn shortest_angle_between(radians0: f64, radians1: f64) -> f64 {
    let cw = angle_in_direction(radians0, radians1, AngleDir::Cw);
    let ccw = angle_in_direction(radians0, radians1, AngleDir::Ccw);
    if cw < ccw { -cw } else { ccw }
}

/// Re-expresses an angle, specified in radians, in the range [-π, π).  If the angle was already
/// in the range [-π, π), it is returned unchanged.
///
/// > **Note:** π (and any angle equivalent to π) maps to -π, not π, since the range is
/// > half-open. If your code depends on π being returned for inputs at exactly ±π, use a
/// > different normalization.
///
/// # Arguments
///
/// * `radians`: the angle to re-express, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::angle_signed_pi;
/// use std::f64::consts::PI;
/// use approx::assert_relative_eq;
/// let new_angle = angle_signed_pi(2.5 * PI);
/// assert_relative_eq!(new_angle, PI / 2.0, epsilon = 1.0e-10);
/// ```
pub fn angle_signed_pi(radians: f64) -> f64 {
    (radians + PI).rem_euclid(2.0 * PI) - PI
}

/// Re-expresses an angle, specified in radians, in the range [0, 2pi].  If the angle was already
/// in the range [0, 2pi], it is returned unchanged.
///
/// # Arguments
///
/// * `radians`: The angle to re-express, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::angle_to_2pi;
/// use std::f64::consts::PI;
/// use approx::assert_relative_eq;
/// let new_angle = angle_to_2pi(-PI);
/// assert_relative_eq!(new_angle, PI, epsilon = 1.0e-10);
/// ```
pub fn angle_to_2pi(radians: f64) -> f64 {
    let mut angle = radians % (2.0 * PI);
    if angle < 0.0 {
        angle += 2.0 * PI;
    }
    angle
}

/// Returns the signed complement of an angle with respect to a full rotation (2π).
///
/// Given an angle, this returns the remaining arc needed to complete a full circle, with the same
/// sign convention: a positive input returns a negative complement (since the remaining arc goes
/// the other way), and vice versa. The result is in the range (-2π, 2π].
///
/// This is useful when you have an arc angle and need the arc that covers the rest of the circle.
///
/// # Arguments
///
/// * `radians`: the angle to complement, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::signed_compliment_2pi;
/// use std::f64::consts::{PI, FRAC_PI_2};
/// use approx::assert_relative_eq;
///
/// // A quarter turn CCW leaves three-quarter turns in the CW direction
/// let a = signed_compliment_2pi(FRAC_PI_2);
/// assert_relative_eq!(a, -3.0 * FRAC_PI_2, epsilon = 1.0e-10);
///
/// // A quarter turn CW leaves three-quarter turns in the CCW direction
/// let b = signed_compliment_2pi(-FRAC_PI_2);
/// assert_relative_eq!(b, 3.0 * FRAC_PI_2, epsilon = 1.0e-10);
/// ```
pub fn signed_compliment_2pi(radians: f64) -> f64 {
    if radians >= 0.0 {
        (-2.0 * PI) + radians
    } else {
        2.0 * PI + radians
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::{Circle2, Iso2};
    use approx::assert_relative_eq;
    use rand::RngExt;
    use test_case::test_case;

    /// When the tolerance corresponds to the sagitta of a chord spanning exactly 1/n of the
    /// circle, the answer must be n.
    #[test_case(8)]
    #[test_case(16)]
    #[test_case(100)]
    fn arc_segments_exact_division(n: usize) {
        let tol = 1.0 - (PI / n as f64).cos();
        assert_eq!(arc_segments_for_tol(1.0, TAU, tol).unwrap(), n);
    }

    /// At the returned value of n, the sagitta of a chord spanning 1/n of the circle must be within
    /// tol. At n-1 it must exceed tol; otherwise, the count is either insufficient or wasteful.
    #[test]
    fn arc_segments_are_sufficient_and_not_excessive() {
        let radius = 2.5;
        for &tol in &[0.5, 0.1, 0.01, 1.0e-4, 1.0e-8] {
            let n = arc_segments_for_tol(radius, TAU, tol).unwrap();
            let sag = |k: usize| radius * (1.0 - (PI / k as f64).cos());
            assert!(sag(n) <= tol, "n={n} insufficient for tol={tol}");
            if n > MIN_FULL_CIRCLE_SEGMENTS {
                assert!(sag(n - 1) > tol, "n={n} excessive for tol={tol}");
            }
        }
    }

    #[test]
    fn arc_segments_minimum_guard() {
        // Loose tolerances fall back to the floor, proportional to the sweep
        assert_eq!(arc_segments_for_tol(1.0, TAU, 10.0).unwrap(), 8);
        assert_eq!(arc_segments_for_tol(1.0, PI, 10.0).unwrap(), 4);
        assert_eq!(arc_segments_for_tol(1.0, PI / 2.0, 10.0).unwrap(), 2);
        // Even a tiny sweep still gets one segment
        assert_eq!(arc_segments_for_tol(1.0, 1.0e-6, 10.0).unwrap(), 1);
    }

    #[test]
    fn arc_segments_sweep_is_unsigned_and_clamped() {
        let tol = 1.0 - (PI / 16.0).cos();
        assert_eq!(
            arc_segments_for_tol(1.0, -TAU, tol).unwrap(),
            arc_segments_for_tol(1.0, TAU, tol).unwrap()
        );
        assert_eq!(
            arc_segments_for_tol(1.0, 5.0 * TAU, tol).unwrap(),
            arc_segments_for_tol(1.0, TAU, tol).unwrap()
        );
    }

    #[test]
    fn arc_segments_invalid_arguments() {
        assert!(arc_segments_for_tol(1.0, TAU, 0.0).is_err());
        assert!(arc_segments_for_tol(1.0, TAU, -1.0).is_err());
        assert!(arc_segments_for_tol(1.0, TAU, f64::NAN).is_err());
        assert!(arc_segments_for_tol(1.0, TAU, f64::INFINITY).is_err());
        assert!(arc_segments_for_tol(0.0, TAU, 0.1).is_err());
        assert!(arc_segments_for_tol(-1.0, TAU, 0.1).is_err());
        assert!(arc_segments_for_tol(f64::NAN, TAU, 0.1).is_err());
    }

    /// A tolerance so tight that the segment count would exceed a u32 index range must be an
    /// error rather than an absurd allocation request downstream.
    #[test]
    fn arc_segments_unresolvable_tolerance() {
        assert!(arc_segments_for_tol(1.0, TAU, 1.0e-20).is_err());
    }

    #[test]
    fn angle_dir_to_sign() {
        assert_eq!(AngleDir::Cw.to_sign(), -1.0);
        assert_eq!(AngleDir::Ccw.to_sign(), 1.0);
    }

    #[test]
    fn angle_dir_from_sign() {
        assert_eq!(AngleDir::from_sign(-1.0), AngleDir::Cw);
        assert_eq!(AngleDir::from_sign(1.0), AngleDir::Ccw);
        assert_eq!(AngleDir::from_sign(0.0), AngleDir::Ccw);
    }

    #[test]
    fn angle_dir_opposite() {
        assert_eq!(AngleDir::Cw.opposite(), AngleDir::Ccw);
        assert_eq!(AngleDir::Ccw.opposite(), AngleDir::Cw);
    }

    #[test]
    fn angle_dir_roundtrip() {
        assert_eq!(AngleDir::from_sign(AngleDir::Cw.to_sign()), AngleDir::Cw);
        assert_eq!(AngleDir::from_sign(AngleDir::Ccw.to_sign()), AngleDir::Ccw);
    }

    #[test_case(90.0, -270.0)]
    #[test_case(180.0, -180.0)]
    #[test_case(270.0, -90.0)]
    #[test_case(-91.0, 269.0)]
    #[test_case(-181.0, 179.0)]
    #[test_case(-271.0, 89.0)]
    fn test_signed_compliment_0(angle: f64, compliment: f64) {
        let test = signed_compliment_2pi(angle.to_radians());
        assert_relative_eq!(test, compliment.to_radians(), epsilon = 1.0e-10);
    }

    #[test]
    fn stress_test_angle_to_2pi() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let angle = rnd.random_range(-8.0 * PI..8.0 * PI);
            let test = angle_to_2pi(angle);
            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle to 2pi: angle={}, test={}",
                angle,
                test
            );

            assert_relative_eq!(f64::sin(angle), f64::sin(test), epsilon = 1.0e-10);
            assert_relative_eq!(f64::cos(angle), f64::cos(test), epsilon = 1.0e-10);
        }
    }

    #[test]
    fn stress_test_angle_signed_pi() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let angle = rnd.random_range(-8.0 * PI..8.0 * PI);
            let test = angle_signed_pi(angle);
            assert!(
                (-PI..PI).contains(&test),
                "Failed Angle Signed Pi: angle={}, test={}",
                angle,
                test
            );

            assert_relative_eq!(f64::sin(angle), f64::sin(test), epsilon = 1.0e-10);
            assert_relative_eq!(f64::cos(angle), f64::cos(test), epsilon = 1.0e-10);
        }
    }

    #[test]
    fn stress_test_angle_in_direction_counterclockwise() {
        let mut rnd = rand::rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let end = rnd.random_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Ccw);

            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle in Direction CW: start={}, end={}, test={} not in [0, 2pi]",
                start,
                end,
                test
            );

            // By rotating the start vector by the test angle, we should get the end vector
            let rot = Iso2::rotation(test);
            let test0 = rot * v0;

            assert_relative_eq!(test0.x, v1.x, epsilon = 1.0e-10);
            assert_relative_eq!(test0.y, v1.y, epsilon = 1.0e-10);
        }
    }
    #[test]
    fn stress_test_angle_in_direction_clockwise() {
        let mut rnd = rand::rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let end = rnd.random_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Cw);

            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle in Direction CW: start={}, end={}, test={} not in [0, 2pi]",
                start,
                end,
                test
            );

            // By rotating the start vector by the test angle, we should get the end vector
            let rot = Iso2::rotation(-test);
            let test0 = rot * v0;

            assert_relative_eq!(test0.x, v1.x, epsilon = 1.0e-10);
            assert_relative_eq!(test0.y, v1.y, epsilon = 1.0e-10);
        }
    }
}
