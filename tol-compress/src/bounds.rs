//! Axis-aligned bounding intervals in any number of dimensions.
//!
//! The encoders need somewhere to anchor each dimension's lattice and something to divide by. That
//! is the whole job, which is why this is a plain min/max pair rather than a dependency on a
//! geometry library's bounding volume type.

use crate::error::{Error, Result};

/// An axis-aligned box in `N` dimensions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds<const N: usize> {
    /// The lower corner.
    pub mins: [f64; N],
    /// The upper corner.
    pub maxs: [f64; N],
}

impl<const N: usize> Bounds<N> {
    /// Build from an explicit pair of corners.
    ///
    /// # Panics
    ///
    /// Panics if any `maxs` component is below its `mins` counterpart, or if any component is not
    /// finite. An inverted or non-finite box would produce a negative or `NaN` range and silently
    /// poison every value quantized against it.
    pub fn new(mins: [f64; N], maxs: [f64; N]) -> Self {
        for i in 0..N {
            assert!(
                mins[i].is_finite() && maxs[i].is_finite(),
                "non-finite bound on axis {i}: [{}, {}]",
                mins[i],
                maxs[i]
            );
            assert!(
                maxs[i] >= mins[i],
                "inverted bound on axis {i}: [{}, {}]",
                mins[i],
                maxs[i]
            );
        }
        Self { mins, maxs }
    }

    /// Compute the bounds enclosing a set of points.
    ///
    /// # Errors
    ///
    /// [`Error::Malformed`] if `points` is empty, since an empty set has no bounding box, or if
    /// any coordinate is not finite. Callers that legitimately handle empty input, such as the
    /// point encoder writing a count of zero, should check for it before calling.
    pub fn from_points(points: &[[f64; N]]) -> Result<Self> {
        let Some(first) = points.first() else {
            return Err(Error::Malformed("cannot bound an empty point set"));
        };

        let mut mins = *first;
        let mut maxs = *first;

        for p in points {
            for i in 0..N {
                if !p[i].is_finite() {
                    return Err(Error::Malformed(
                        "point set contains a non-finite coordinate",
                    ));
                }
                if p[i] < mins[i] {
                    mins[i] = p[i];
                }
                if p[i] > maxs[i] {
                    maxs[i] = p[i];
                }
            }
        }

        Ok(Self { mins, maxs })
    }

    /// The span along each axis. Components may be zero, which is the degenerate case the encoders
    /// store for free.
    pub fn extents(&self) -> [f64; N] {
        let mut out = [0.0; N];
        for (i, e) in out.iter_mut().enumerate() {
            *e = self.maxs[i] - self.mins[i];
        }
        out
    }

    /// Whether a point lies within the box, inclusive of its faces.
    pub fn contains(&self, point: &[f64; N]) -> bool {
        (0..N).all(|i| point[i] >= self.mins[i] && point[i] <= self.maxs[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_a_point_set() {
        let points = [[1.0, -2.0, 3.0], [-4.0, 5.0, 0.5], [0.0, 0.0, 9.0]];
        let b = Bounds::from_points(&points).unwrap();

        assert_eq!(b.mins, [-4.0, -2.0, 0.5]);
        assert_eq!(b.maxs, [1.0, 5.0, 9.0]);
        assert_eq!(b.extents(), [5.0, 7.0, 8.5]);
    }

    #[test]
    fn a_single_point_gives_a_degenerate_box() {
        let b = Bounds::from_points(&[[2.0, 3.0]]).unwrap();

        assert_eq!(b.mins, b.maxs);
        assert_eq!(b.extents(), [0.0, 0.0]);
        assert!(b.contains(&[2.0, 3.0]));
    }

    /// Planar data, which the previous implementation charged two bytes per point to store on the
    /// axis with no information in it.
    #[test]
    fn a_flat_axis_has_zero_extent() {
        let points = [[0.0, 0.0, 5.0], [1.0, 2.0, 5.0], [-3.0, 1.0, 5.0]];
        let b = Bounds::from_points(&points).unwrap();

        assert_eq!(b.extents()[2], 0.0);
    }

    #[test]
    fn empty_input_is_an_error_not_a_box() {
        let empty: [[f64; 3]; 0] = [];
        assert!(matches!(
            Bounds::from_points(&empty),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn non_finite_coordinates_are_rejected() {
        assert!(Bounds::from_points(&[[0.0, 1.0], [f64::NAN, 0.0]]).is_err());
        assert!(Bounds::from_points(&[[0.0, 1.0], [f64::INFINITY, 0.0]]).is_err());

        // A NaN in the first point must be caught too, not adopted as the seed and compared away:
        // every comparison against NaN is false, so a careless implementation keeps it.
        assert!(Bounds::from_points(&[[f64::NAN, 0.0], [1.0, 1.0]]).is_err());
    }

    #[test]
    fn works_beyond_three_dimensions() {
        let points = [[1.0, 2.0, 3.0, 4.0, 5.0], [-1.0, 0.0, 8.0, 4.0, 2.0]];
        let b = Bounds::from_points(&points).unwrap();

        assert_eq!(b.extents(), [2.0, 2.0, 5.0, 0.0, 3.0]);
    }

    #[test]
    fn contains_is_inclusive() {
        let b = Bounds::new([0.0, 0.0], [1.0, 1.0]);

        assert!(b.contains(&[0.0, 0.0]));
        assert!(b.contains(&[1.0, 1.0]));
        assert!(b.contains(&[0.5, 0.5]));
        assert!(!b.contains(&[1.000_001, 0.5]));
        assert!(!b.contains(&[-0.5, 0.5]));
    }

    #[test]
    #[should_panic(expected = "inverted bound on axis 1")]
    fn inverted_bounds_panic() {
        Bounds::new([0.0, 5.0], [1.0, 2.0]);
    }

    #[test]
    #[should_panic(expected = "non-finite bound on axis 0")]
    fn non_finite_bounds_panic() {
        Bounds::new([f64::NAN, 0.0], [1.0, 1.0]);
    }
}
