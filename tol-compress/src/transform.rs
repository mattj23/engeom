//! This module is for reading, writing, and applying transformations to point partitions so that
//! the points inside it can be quantized over a smaller space, rather than whatever they happened
//! to be in world axes.
//!
//! # The idea in theory
//!
//! In theory, a scanned surface is a typically a hollow shell, so a partition ends up being
//! something akin to a patch of it (depending on point order, naturally). A patch that lies
//! obliquely to the world axes will fill it's axis-aligned box almost completely, while the same
//! patch looked at in a frame tangent to the surface will usually end up with at least one thin axis.
//!
//! That's the goal of this module; to use rotation to produce bounding schemes that are less sparse
//! so the encoding isn't paying for a lot of empty space.  It hasn't quite worked out so far; see
//! the top-level section below.
//!
//! ## The rotation is the only thing stored
//!
//! A rotated partition needs a frame, a box in that frame, and widths. Everything but the rotation
//! is already derivable from what the block records:
//!
//! - The frame's origin is the anchor's center, so no partition needs to store it seperately. Local
//!   coordinates are `R^T (p - c)`, which keeps the f64 arithmetic proportional to the block's own
//!   extent instead of how far the part sits from the world origin.
//! - The lattice its corners are placed on is the rotated anchor: the axis-aligned box that encloses
//!   the anchor once turned into the local frame. Both sides can compute it, so
//!   [`crate::segment::corner_codes`] and [`crate::segment::decode_corners`] are used unchanged with
//!   that box in place of the anchor. Every point is inside the anchor, so its local box is inside
//!   the rotated anchor and nothing needs to clamp.
//!
//! The rotated lattice is at most `sqrt(N)` coarser than the world one, which can widen a local box
//! by a step or two. That costs a fraction of a bit and doesn't reduce accuracy.
//!
//! ## Quantization doesn't degrade recovery accuracy
//!
//! The encoder quantizes the rotation first and then works in the frame the quantized rotation
//! defines. The decoder rebuilds that same rotation from the same bits, and a rotation is an
//! isometry, so error measured in the local frame is error in the world.
//!
//! So a coarse rotation costs a slightly looser box but not less accurate recovery. That's why the
//! twelve bit encoding is reasonable, it quantizes to about 3e-4 radius, which over a patch a few
//! millimeters wide would grow a box extent by about a micron.
//!
//! [`Rotation`] can only be built from its packed code, by [`Rotation::from_code`], which every
//! other constructor goes through. There is no way to hold a rotation the writer could use and the
//! reader could not reproduce.
//!
//! ## Dimension
//!
//! Two and three dimensions only. A plane rotation is one angle and a spatial one is a quaternion;
//! above that there is no compact general form worth defining for a format. [`supported`] is the
//! predicate, [`Rotation::from_axes`] and [`Rotation::from_direction`] are the builders, and
//! [`crate::points`] refuses to read a transform flag at any other dimension.
//!
//! # The idea in practice, unused so far
//!
//! Both the read and write ends of this are implemented and neither is used. The format carries the
//! flag, the writer emites a rotation when handed one, and the reader does indeed apply it.
//! However, no estimator currently runs, so as of right now every partition is written on world
//! axes and no file produced by this crate has the flag set yet.
//!
//! The problem wasn't that it didn't shrink the files, it definitely did. The size reduction wasn't
//! huge, but on dense, realistic scan I did see the point block shrink by about 14%, which was
//! still worth keeping.
//!
//! The problem turned out to be that it was caught by the test suite for breaking the stability
//! guarantee I've been operating under so far: that repeated read/write loops of the same data
//! would eventually settle to a fixed value roughly close to the original tolerance.  For some
//! reason rotation broke that and I haven't yet decided how to proceed.
//!
//! ## What was measured
//!
//! I ended up trying with two different estimators.  The first was PCA based, fitting the
//! eigenvectors of covariance, which is known to produce oriented bounding boxes that are at best
//! mid, but also happens to be quick and easy to implement.
//!
//! The second method I tried was DiTO-14, from Larsson and Kallberg's "Fast Computation of
//! Tight-Fitting Oriented Bounding Boxes", which builds candidate frames from extremal points and
//! is both cheaper and better regarded. On coordinate blocks the two came out level: the 370k point
//! turbine stator at 1e-3 gave 1,787,809 bytes under the PCA fit and 1,789,821 under DiTO, against
//! 1,879,074 on world axes. So DiTO reaches the same size without an eigensolver, and if this is
//! picked up again it is the one to pick up.
//!
//! The win was real but not large once it reached a whole file, because a mesh's index block dwarfs
//! its coordinates: **-4.9% of coordinates and -1.9% of the file** on that stator, and -14.3% and
//! -3.8% on an unthinned phone scan at 1e-2. The scan does better because it is mostly large smooth
//! panels, which is the shape a tangent frame flatters.
//!
//! ## Why it is not switched on
//!
//! Writing a file, reading it back and writing it again stopped settling. Without frames, that
//! cycle stops in two to four passes and every point stays inside about 1.5x the tolerance of where
//! it started. With frames it did not close at all on some inputs, and points wandered off at
//! `sqrt(passes)`, reaching thirty times the tolerance and still climbing. No single write ever
//! broke its promise, but a pipeline that stores and restores repeatedly could walk a part out
//! indefinitely, and my early intentions for the format was to make that impossible.
//!
//! There were two causes that I _was_ able to find and fix:
//!
//! - **A corner sitting on a lattice point.** [`crate::segment::corner_codes`] steps a box outward
//!   until it encloses, and the comparison it steps on was exact, so which side of a lattice point
//!   the arithmetic happened to land decided whether the box grew by a whole step. It now compares
//!   with a margin. That fix is in and stands on its own.
//! - **Sliver triangles.** A candidate frame takes its second axis from a triangle's normal, and a
//!   nearly-collinear triangle names that direction barely at all, so a nudge smaller than the
//!   tolerance swings it a long way. This is the same ill-conditioning that makes an eigenvector
//!   unusable when two eigenvalues meet, and it is why swapping PCA for DiTO did not help: the two
//!   estimators fail on the same geometry for the same reason, through different doors. Rejecting
//!   thin triangles fixed the synthetic cases and cost about half the size win.
//!
//! But, even with that, the Stanford bunny at 3e-6 still fails to settle, and drifts about 1/3rd of
//! its tolerance each path...but steadily not randomly, so it's acting like a ratchet and not a
//! random walk. I'm not really sure at this point what's causing it, and the space savings from the
//! transforms didn't feel big enough to be worth dumping a ton of time into it right now.
//!
//! ## What would have to change to finish it
//!
//! An estimator alone is encoder-side work and slots into [`crate::segment`]'s `finish_run` without
//! touching anything here. But if the suspicion above is right, the fix is to give an oriented
//! partition a box of its own rather than corners on a shared derived lattice, and that changes what
//! the reader parses. Since no file has the flag set, such a change costs an edit and nothing else:
//! there is nothing in the wild to stay compatible with.
//!
//! The other thing that could change is that I could decide that the rewrite stabilization
//! guarantee isn't realistic and/or useful.  I need to check what happens if you transform a mesh
//! around before rewriting it, or remove points/faces, because if that drifts in the same way and
//! isn't fixable than the current guarantee probably isn't realistic...there isn't really a lot of
//! loading/saving the same exact data, since in most real-world programs you're just going to copy
//! the file whole.

use crate::bits::{BitReader, BitWriter};
use crate::bounds::Bounds;
use crate::error::Result;
use crate::quantize::Quantizer;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};
use std::io::Write;

/// Bits per stored quaternion component, and bits for the plane angle.
///
/// Both are the same number by coincidence rather than for a shared reason.
const COMPONENT_BITS: u8 = 12;
const ANGLE_BITS: u8 = 12;

/// Whether a rotation can be encoded at this dimension.
pub const fn supported(dimension: usize) -> bool {
    dimension == 2 || dimension == 3
}

/// The bits a rotation occupies at this dimension, or zero where none can be written.
///
/// Two bits name the largest quaternion component in three dimensions and the other three follow at
/// twelve bits each. A plane rotation is a single angle, also at twelve bits.
pub const fn bits(dimension: usize) -> u8 {
    match dimension {
        2 => ANGLE_BITS,
        3 => 2 + 3 * COMPONENT_BITS,
        _ => 0,
    }
}

/// The frame a partition's points are written in.
///
/// The matrix takes local coordinates to world ones: its columns are the local axes. It is always
/// derived from `code`, never from whatever was fitted, so the writer and the reader cannot end up
/// holding different rotations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rotation<const N: usize> {
    /// The packed form, exactly as it goes into the stream.
    code: u64,
    /// `m[i][j]` is row `i`, column `j`. World is `centre + m * local`.
    m: [[f64; N]; N],
}

impl<const N: usize> Rotation<N> {
    /// Rebuild a rotation from its packed code.
    ///
    /// Total, by construction. Every code is a rotation: the plane angle is read straight out, and
    /// a quaternion whose components decode to something longer than a unit is normalized rather
    /// than rejected. This is deliberate, since it means a corrupt transform produces geometry in
    /// the wrong place rather than an error path the reader has to carry.
    ///
    /// # Panics
    ///
    /// Panics if `N` is not [`supported`]. Callers reading from a stream must check the dimension
    /// first, which [`crate::points`] does before it reads the flag's payload.
    pub fn from_code(code: u64) -> Self {
        assert!(supported(N), "no rotation is defined in {N} dimensions");

        // The two builders work at a fixed dimension, so each copies into the generic array a cell
        // at a time. `N` is 2 or 3 here, so the loops stay inside both.
        let mut m = [[0.0f64; N]; N];
        match N {
            2 => {
                let fixed = plane_matrix(code);
                for i in 0..2 {
                    for j in 0..2 {
                        m[i][j] = fixed[i][j];
                    }
                }
            }
            3 => {
                let fixed = spatial_matrix(code);
                for i in 0..3 {
                    for j in 0..3 {
                        m[i][j] = fixed[i][j];
                    }
                }
            }
            _ => unreachable!("guarded by the assertion above"),
        }

        Self { code, m }
    }

    /// A world point in this frame, relative to `centre`.
    pub fn to_local(&self, p: &[f64; N], centre: &[f64; N]) -> [f64; N] {
        let d: [f64; N] = std::array::from_fn(|i| p[i] - centre[i]);
        // The transpose, since the matrix takes local to world and this goes the other way.
        std::array::from_fn(|i| (0..N).map(|j| self.m[j][i] * d[j]).sum())
    }

    /// The inverse of [`Rotation::to_local`].
    pub fn to_world(&self, q: &[f64; N], centre: &[f64; N]) -> [f64; N] {
        std::array::from_fn(|i| centre[i] + (0..N).map(|j| self.m[i][j] * q[j]).sum::<f64>())
    }

    /// Append the rotation to a bitstream.
    pub fn write<W: Write>(&self, bw: &mut BitWriter<W>) -> Result<()> {
        bw.write_bits(self.code, bits(N))
    }

    /// Read a rotation from a bitstream.
    ///
    /// # Panics
    ///
    /// Panics if `N` is not [`supported`], as [`Rotation::from_code`] does.
    pub fn read(br: &mut BitReader) -> Result<Self> {
        Ok(Self::from_code(br.read_bits(bits(N))?))
    }
}

impl Rotation<2> {
    /// The frame whose first axis points along `direction`.
    ///
    /// `direction` need not be a unit vector and its sign does not matter: a box is symmetric under
    /// a quarter turn, so only the direction modulo ninety degrees is kept.
    pub fn from_direction(direction: [f64; 2]) -> Self {
        Self::from_code(plane_code(direction))
    }
}

impl Rotation<3> {
    /// The frame whose first two axes are `along` and `up`.
    ///
    /// Both are normalized here, and the third axis is their cross product, so what comes back is a
    /// rotation rather than a reflection. `up` is squared off against `along` first, which makes a
    /// pair that is only approximately perpendicular usable rather than an error: a caller deriving
    /// axes from geometry has floating point in the way.
    ///
    /// Returns `None` if either vector has no length or the two are parallel, which name no frame.
    pub fn from_axes(along: [f64; 3], up: [f64; 3]) -> Option<Self> {
        let unit = |v: [f64; 3]| {
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            (n > 0.0).then(|| [v[0] / n, v[1] / n, v[2] / n])
        };
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };

        let u = unit(along)?;
        let dot = u[0] * up[0] + u[1] * up[1] + u[2] * up[2];
        let v = unit([up[0] - dot * u[0], up[1] - dot * u[1], up[2] - dot * u[2]])?;
        let w = cross(u, v);

        // Columns are the axes, so the matrix takes a local coordinate back to the world.
        let m = [[u[0], v[0], w[0]], [u[1], v[1], w[1]], [u[2], v[2], w[2]]];
        Some(Self::from_code(pack_quaternion(quaternion_of(&m))))
    }
}

/// The box that encloses `anchor` once turned into `rotation`'s frame, centred on the origin.
///
/// This is the lattice a rotated partition's corners are placed on, and it is what makes the corner
/// coding identical to the axis-aligned case. Rotating a box gives half-extents
/// `h'_i = sum_j |m[j][i]| h_j`, which needs no corner enumeration.
pub fn rotated_anchor<const N: usize>(anchor: &Bounds<N>, rotation: &Rotation<N>) -> Bounds<N> {
    let half: [f64; N] = std::array::from_fn(|j| (anchor.maxs[j] - anchor.mins[j]) / 2.0);

    let out: [f64; N] =
        std::array::from_fn(|i| (0..N).map(|j| rotation.m[j][i].abs() * half[j]).sum());

    Bounds {
        mins: std::array::from_fn(|i| -out[i]),
        maxs: out,
    }
}

/// The centre of a box, which is where a rotated partition's frame sits.
pub fn centre_of<const N: usize>(anchor: &Bounds<N>) -> [f64; N] {
    std::array::from_fn(|i| (anchor.mins[i] + anchor.maxs[i]) / 2.0)
}

/// The quantizer the plane angle is stored through, over a quarter turn.
///
/// A quarter turn is the whole domain because a box is symmetric under one: turning the frame by
/// ninety degrees swaps two axes, which trades one width for another and costs nothing.
fn angle_quantizer() -> Quantizer {
    Quantizer::new(0.0, FRAC_PI_2, ANGLE_BITS)
}

/// How far either side of zero a stored quaternion component's codes run.
///
/// Three of the four components are stored and the largest is left out, so each stored one is at
/// most `1/sqrt(2)` in magnitude. The codes are laid out symmetrically about this offset rather
/// than through a [`Quantizer`], which spreads `2^bits` codes across a range and so never places
/// one on the midpoint. A component of zero has to be storable as zero, or the identity rotation
/// would come back with a tilt of about 2e-4 radians and an axis-aligned frame would never be
/// exactly axis aligned. One code out of the 4096 goes unused, which buys that.
const COMPONENT_HALF: i64 = (1 << (COMPONENT_BITS - 1)) - 1;

fn encode_component(value: f64) -> u64 {
    let scaled = (value / FRAC_1_SQRT_2).clamp(-1.0, 1.0) * COMPONENT_HALF as f64;
    (scaled.round() as i64 + COMPONENT_HALF) as u64
}

fn decode_component(code: u64) -> f64 {
    (code as i64 - COMPONENT_HALF) as f64 / COMPONENT_HALF as f64 * FRAC_1_SQRT_2
}

/// Fold a direction in the plane into the quarter turn the angle coding spans, and pack it.
fn plane_code(direction: [f64; 2]) -> u64 {
    let angle = direction[1].atan2(direction[0]).rem_euclid(FRAC_PI_2);
    angle_quantizer().encode(angle)
}

fn plane_matrix(code: u64) -> [[f64; 2]; 2] {
    let angle = angle_quantizer().decode(code);
    let (s, c) = angle.sin_cos();
    [[c, -s], [s, c]]
}

fn spatial_matrix(code: u64) -> [[f64; 3]; 3] {
    matrix_of(unpack_quaternion(code))
}

/// The unit quaternion of a rotation matrix, as `[x, y, z, w]`.
///
/// Shepperd's method: build from whichever component the trace says is largest, so the square root
/// is never taken of something near zero.
fn quaternion_of(m: &[[f64; 3]; 3]) -> [f64; 4] {
    let trace = m[0][0] + m[1][1] + m[2][2];

    let q = if trace > 0.0 {
        let s = (trace + 1.0).sqrt() * 2.0;
        [
            (m[2][1] - m[1][2]) / s,
            (m[0][2] - m[2][0]) / s,
            (m[1][0] - m[0][1]) / s,
            0.25 * s,
        ]
    } else if m[0][0] > m[1][1] && m[0][0] > m[2][2] {
        let s = (1.0 + m[0][0] - m[1][1] - m[2][2]).sqrt() * 2.0;
        [
            0.25 * s,
            (m[0][1] + m[1][0]) / s,
            (m[0][2] + m[2][0]) / s,
            (m[2][1] - m[1][2]) / s,
        ]
    } else if m[1][1] > m[2][2] {
        let s = (1.0 + m[1][1] - m[0][0] - m[2][2]).sqrt() * 2.0;
        [
            (m[0][1] + m[1][0]) / s,
            0.25 * s,
            (m[1][2] + m[2][1]) / s,
            (m[0][2] - m[2][0]) / s,
        ]
    } else {
        let s = (1.0 + m[2][2] - m[0][0] - m[1][1]).sqrt() * 2.0;
        [
            (m[0][2] + m[2][0]) / s,
            (m[1][2] + m[2][1]) / s,
            0.25 * s,
            (m[1][0] - m[0][1]) / s,
        ]
    };

    normalized(q)
}

fn matrix_of(q: [f64; 4]) -> [[f64; 3]; 3] {
    let [x, y, z, w] = q;

    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

fn normalized(q: [f64; 4]) -> [f64; 4] {
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if norm == 0.0 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm]
}

/// Two bits naming the largest component, then the other three at [`COMPONENT_BITS`].
///
/// The largest is left out because it can be rebuilt from the other three, and it is the largest
/// that is left out so the rebuilt one is the least sensitive to their rounding. A quaternion and
/// its negation are the same rotation, so the sign is normalized away rather than stored.
fn pack_quaternion(q: [f64; 4]) -> u64 {
    let mut largest = 0usize;
    for i in 1..4 {
        if q[i].abs() > q[largest].abs() {
            largest = i;
        }
    }

    let sign = if q[largest] < 0.0 { -1.0 } else { 1.0 };

    let mut code = largest as u64;
    let mut slot = 0;
    for (i, &value) in q.iter().enumerate() {
        if i == largest {
            continue;
        }
        code |= encode_component(sign * value) << (2 + slot * u32::from(COMPONENT_BITS));
        slot += 1;
    }

    code
}

fn unpack_quaternion(code: u64) -> [f64; 4] {
    let largest = (code & 0b11) as usize;
    let mask = (1u64 << COMPONENT_BITS) - 1;

    let mut q = [0.0f64; 4];
    let mut slot = 0;
    for (i, value) in q.iter_mut().enumerate() {
        if i == largest {
            continue;
        }
        let bits = (code >> (2 + slot * u32::from(COMPONENT_BITS))) & mask;
        *value = decode_component(bits);
        slot += 1;
    }

    // The one left out, taken positive. Rounding can push the other three past a unit between them,
    // which leaves nothing under the root; the normalization below rescues that case rather than
    // treating it as corrupt.
    let rest = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    q[largest] = (1.0 - rest).max(0.0).sqrt();

    normalized(q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testgen::Rng;

    fn determinant(m: &[[f64; 3]; 3]) -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// A rotation from an axis and an angle, for building test cases by hand.
    fn about(axis: [f64; 3], angle: f64) -> Rotation<3> {
        let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let (s, c) = (angle / 2.0).sin_cos();
        let q = [axis[0] / n * s, axis[1] / n * s, axis[2] / n * s, c];
        Rotation::<3>::from_code(pack_quaternion(q))
    }

    #[test]
    fn the_widths_are_what_the_format_says() {
        assert_eq!(bits(2), 12);
        assert_eq!(bits(3), 38);
        assert_eq!(bits(1), 0);
        assert_eq!(bits(4), 0);
        assert!(supported(2) && supported(3));
        assert!(!supported(1) && !supported(4));
    }

    /// The claim the whole design rests on: a rotation survives its packed form closely enough that
    /// the box it describes is barely looser, and nothing about accuracy depends on how closely.
    #[test]
    fn a_packed_rotation_comes_back_within_a_thousandth_of_a_degree() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0001);
        let mut worst: f64 = 0.0;

        for _ in 0..20_000 {
            // A uniformly distributed rotation, by normalizing a gaussian quaternion.
            let q = normalized([
                rng.gaussian(1.0),
                rng.gaussian(1.0),
                rng.gaussian(1.0),
                rng.gaussian(1.0),
            ]);

            let packed = pack_quaternion(q);
            let back = unpack_quaternion(packed);

            // A quaternion and its negation are the same rotation, so compare the rotations.
            let a = matrix_of(q);
            let b = matrix_of(back);
            let mut err: f64 = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    err = err.max((a[i][j] - b[i][j]).abs());
                }
            }
            worst = worst.max(err);
        }

        // 12 bits over a component range of sqrt(2) is a step of 3.5e-4, and the matrix entries
        // move by about that much.
        assert!(
            worst < 2e-3,
            "a packed rotation moved a matrix entry by {worst}"
        );
    }

    /// Every code has to name a rotation, because the reader has no way to reject one and should
    /// not need one. This walks a wide spread of codes including deliberately impossible ones.
    #[test]
    fn every_code_is_a_rotation() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0002);

        let mut codes: Vec<u64> = vec![0, 1, 2, 3, (1 << 38) - 1];
        // Three components at their extremes, which claims a quaternion longer than a unit.
        let full = (1u64 << COMPONENT_BITS) - 1;
        for largest in 0..4u64 {
            codes.push(largest | (full << 2) | (full << 14) | (full << 26));
            codes.push(largest);
        }
        for _ in 0..5_000 {
            codes.push(rng.next_u64() & ((1 << 38) - 1));
        }

        for code in codes {
            let r = Rotation::<3>::from_code(code);

            // Orthonormal: the columns are unit length and mutually perpendicular.
            for i in 0..3 {
                let len: f64 = (0..3).map(|k| r.m[k][i] * r.m[k][i]).sum::<f64>().sqrt();
                assert!(
                    (len - 1.0).abs() < 1e-12,
                    "code {code}: column {i} has length {len}"
                );
                for j in (i + 1)..3 {
                    let dot: f64 = (0..3).map(|k| r.m[k][i] * r.m[k][j]).sum();
                    assert!(
                        dot.abs() < 1e-12,
                        "code {code}: columns {i} and {j} are not perpendicular, dot {dot}"
                    );
                }
            }

            assert!(
                (determinant(&std::array::from_fn(|i| std::array::from_fn(|j| r.m[i][j]))) - 1.0)
                    .abs()
                    < 1e-12,
                "code {code} is a reflection rather than a rotation"
            );
        }
    }

    /// Local and world are inverses of each other, which is what makes the round trip exact apart
    /// from the lattice.
    #[test]
    fn local_and_world_are_inverses() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0003);
        let r = about([1.0, 2.0, -0.5], 0.9);
        let centre = [12.0, -3.0, 400.0];

        for _ in 0..10_000 {
            let p: [f64; 3] = rng.point(-50.0, 50.0);
            let back = r.to_world(&r.to_local(&p, &centre), &centre);
            for i in 0..3 {
                assert!(
                    (back[i] - p[i]).abs() < 1e-9,
                    "axis {i}: {} came back as {}",
                    p[i],
                    back[i]
                );
            }
        }
    }

    /// A rotation is an isometry, which is the reason none of this touches the tolerance budget.
    /// If distance were not preserved, error measured in the local frame would not be error in the
    /// world and the whole guarantee would move.
    #[test]
    fn the_frame_preserves_distance() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0004);
        let centre = [0.0, 0.0, 0.0];
        let mut worst: f64 = 0.0;

        for _ in 0..2_000 {
            let r = Rotation::<3>::from_code(rng.next_u64() & ((1 << 38) - 1));
            let a: [f64; 3] = rng.point(-100.0, 100.0);
            let b: [f64; 3] = rng.point(-100.0, 100.0);

            let world = (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt();
            let (la, lb) = (r.to_local(&a, &centre), r.to_local(&b, &centre));
            let local = (0..3).map(|i| (la[i] - lb[i]).powi(2)).sum::<f64>().sqrt();

            worst = worst.max((world - local).abs() / world);
        }

        assert!(worst < 1e-14, "distance moved by {worst} of itself");
    }

    /// The frame builders are what an estimator would produce and what the writer would then use,
    /// so they have to hold up on their own now that nothing calls them in anger.
    #[test]
    fn a_frame_built_from_axes_has_those_axes() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_000A);

        for _ in 0..5_000 {
            let along: [f64; 3] = rng.point(-1.0, 1.0);
            let up: [f64; 3] = rng.point(-1.0, 1.0);

            let Some(frame) = Rotation::<3>::from_axes(along, up) else {
                continue;
            };

            // The first axis is the first column, up to the rotation lattice and to `along` being
            // free to point either way along itself once packed.
            let length = (along[0] * along[0] + along[1] * along[1] + along[2] * along[2]).sqrt();
            let unit = [along[0] / length, along[1] / length, along[2] / length];
            let dot: f64 = (0..3).map(|i| frame.m[i][0] * unit[i]).sum();
            assert!(
                dot.abs() > 0.999,
                "the first axis came back at {dot} of where it was put"
            );

            for i in 0..3 {
                let len: f64 = (0..3)
                    .map(|k| frame.m[k][i] * frame.m[k][i])
                    .sum::<f64>()
                    .sqrt();
                assert!((len - 1.0).abs() < 1e-9, "column {i} has length {len}");
            }
        }
    }

    /// Degenerate axes name no frame, and saying so beats returning one built from a zero vector.
    #[test]
    fn a_frame_needs_two_directions_that_differ() {
        assert!(Rotation::<3>::from_axes([0.0; 3], [0.0, 0.0, 1.0]).is_none());
        assert!(Rotation::<3>::from_axes([1.0, 0.0, 0.0], [0.0; 3]).is_none());
        assert!(Rotation::<3>::from_axes([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]).is_none());
        assert!(Rotation::<3>::from_axes([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]).is_none());
    }

    /// A pair that is only approximately perpendicular is what geometry actually hands you, so the
    /// second axis is squared off rather than trusted.
    #[test]
    fn an_approximate_pair_is_squared_off() {
        let frame = Rotation::<3>::from_axes([1.0, 0.0, 0.0], [0.02, 1.0, 0.0]).unwrap();
        let dot: f64 = (0..3).map(|i| frame.m[i][0] * frame.m[i][1]).sum();
        assert!(dot.abs() < 1e-9, "axes are not perpendicular, dot {dot}");
    }

    /// The plane builder, whose quarter-turn symmetry means opposite directions are the same frame.
    #[test]
    fn a_plane_frame_points_along_its_direction() {
        let frame = Rotation::<2>::from_direction([1.0, 1.0]);
        let angle = frame.m[1][0].atan2(frame.m[0][0]);
        assert!(
            (angle - std::f64::consts::FRAC_PI_4).abs() < 1e-3,
            "a diagonal came back at {angle} radians"
        );

        // A box does not care which way along an axis it is measured, so these must agree.
        assert_eq!(
            Rotation::<2>::from_direction([1.0, 1.0]),
            Rotation::<2>::from_direction([-1.0, -1.0])
        );
    }

    /// The rotated anchor has to contain every local point, or a partition's corners would be
    /// placed on a lattice too small to hold them and the coding would clamp.
    #[test]
    fn the_rotated_anchor_contains_what_the_anchor_did() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0007);

        for _ in 0..500 {
            let anchor = Bounds {
                mins: [-3.0, 10.0, 0.5],
                maxs: [7.0, 14.0, 22.0],
            };
            let centre = centre_of(&anchor);
            let r = Rotation::<3>::from_code(rng.next_u64() & ((1 << 38) - 1));
            let local_anchor = rotated_anchor(&anchor, &r);

            for _ in 0..50 {
                let p = [
                    rng.in_range(anchor.mins[0], anchor.maxs[0]),
                    rng.in_range(anchor.mins[1], anchor.maxs[1]),
                    rng.in_range(anchor.mins[2], anchor.maxs[2]),
                ];
                let q = r.to_local(&p, &centre);
                for (i, &v) in q.iter().enumerate() {
                    assert!(
                        v >= local_anchor.mins[i] && v <= local_anchor.maxs[i],
                        "axis {i}: {v} outside [{}, {}]",
                        local_anchor.mins[i],
                        local_anchor.maxs[i]
                    );
                }
            }
        }
    }

    /// A rotation that does nothing must leave the anchor alone, so an oriented partition that
    /// happens to be axis aligned pays nothing for the frame beyond the bits.
    #[test]
    fn an_identity_rotation_leaves_the_anchor_where_it_was() {
        let anchor = Bounds {
            mins: [-3.0, 10.0, 0.5],
            maxs: [7.0, 14.0, 22.0],
        };
        let r = Rotation::<3>::from_code(pack_quaternion([0.0, 0.0, 0.0, 1.0]));
        let rotated = rotated_anchor(&anchor, &r);

        for i in 0..3 {
            let half = (anchor.maxs[i] - anchor.mins[i]) / 2.0;
            assert!((rotated.maxs[i] - half).abs() < 1e-12);
            assert!((rotated.mins[i] + half).abs() < 1e-12);
        }
    }

    /// A round trip through the stream, which is the only way a rotation ever reaches a reader.
    #[test]
    fn a_rotation_survives_the_bitstream() {
        let mut rng = Rng::new(0x5EED_1234_ABCD_0008);

        for _ in 0..1_000 {
            let original = Rotation::<3>::from_code(rng.next_u64() & ((1 << 38) - 1));

            let mut buf = Vec::new();
            {
                let mut bw = BitWriter::new(&mut buf);
                original.write(&mut bw).unwrap();
                bw.finish().unwrap();
            }
            assert_eq!(buf.len(), 5, "38 bits occupy five bytes once padded");

            let mut br = BitReader::new(&buf);
            assert_eq!(Rotation::<3>::read(&mut br).unwrap(), original);
        }
    }

    /// A frame that is already axis aligned has to survive as one.
    ///
    /// This is what the offset component coding buys. Spread across a [`Quantizer`]'s range instead,
    /// zero falls between two codes and the identity comes back tilted by about 2e-4 radians, which
    /// costs nothing in accuracy but is a strange thing for a format to be unable to say.
    #[test]
    fn an_axis_aligned_frame_is_exact() {
        let identity = Rotation::<3>::from_code(pack_quaternion([0.0, 0.0, 0.0, 1.0]));

        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert_eq!(identity.m[i][j], want, "entry ({i}, {j})");
            }
        }

        // A quarter turn about each axis is the other frame a box can be written in without
        // changing shape. These are exact in the coding but not in the matrix: the quaternion of a
        // quarter turn has components of `1/sqrt(2)`, and squaring that in binary gives
        // `0.5000000000000001`, so an entry that should be zero lands two ulps away instead.
        for axis in 0..3 {
            let mut e = [0.0; 3];
            e[axis] = 1.0;
            let turned = about(e, FRAC_PI_2);
            for i in 0..3 {
                for j in 0..3 {
                    let v = turned.m[i][j];
                    let nearest = v.round();
                    assert!(
                        (v - nearest).abs() < 1e-15 && nearest.abs() <= 1.0,
                        "a quarter turn about axis {axis} produced {v} at ({i}, {j})"
                    );
                }
            }
        }
    }

    #[test]
    fn a_plane_rotation_survives_the_bitstream() {
        for code in 0..(1u64 << ANGLE_BITS) {
            let original = Rotation::<2>::from_code(code);

            let mut buf = Vec::new();
            {
                let mut bw = BitWriter::new(&mut buf);
                original.write(&mut bw).unwrap();
                bw.finish().unwrap();
            }
            assert_eq!(buf.len(), 2);

            let mut br = BitReader::new(&buf);
            assert_eq!(Rotation::<2>::read(&mut br).unwrap(), original);
        }
    }

    #[test]
    #[should_panic(expected = "no rotation is defined in 4 dimensions")]
    fn a_dimension_without_an_encoding_panics_rather_than_inventing_one() {
        let _ = Rotation::<4>::from_code(0);
    }
}
