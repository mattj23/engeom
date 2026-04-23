//! This module has tools for partitioning curves into sub-curves.

use crate::common::points::dist;
use crate::geom2::LineOps2;
use crate::na::Unit;
use crate::{Curve2, Point2};
use parry3d_f64::query::SplitResult;

impl Curve2 {
    /// Returns a curve portion between the station at lengths `a` and `b` *which includes* the
    /// part of the curve passing through the control length.  This is useful for partitioning a
    /// closed curve into a part when you know the endpoints and any point in the middle, but you
    /// don't necessarily know the order of the points.
    ///
    /// # Arguments
    ///
    /// * `a`: a length along the curve to be the start or end of the new curve
    /// * `b`: a length along the curve to be the start or end of the new curve
    /// * `control`: a length along the curve which must be included in the new curve
    ///
    /// returns: Option<Curve2>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn between_lengths_by_control(&self, a: f64, b: f64, control: f64) -> Option<Self> {
        if control > self.length() {
            return None;
        }

        let lower = a.min(b);
        let upper = a.max(b);
        if lower < control && control < upper {
            self.between_lengths(lower, upper)
        } else if control < lower || control > upper && self.is_closed {
            self.between_lengths(upper, lower)
        } else {
            None
        }
    }

    /// Returns a curve portion between the section at length l0 and l1. If the curve is not closed,
    /// the case where l1 < l0 will return None. If the curve is closed, the portion of the curve
    /// which is returned will depend on whether l0 is larger or smaller than l1.
    ///
    /// The new curve will begin at the point corresponding with l0.
    ///
    /// # Arguments
    ///
    /// * `l0`:
    /// * `l1`:
    ///
    /// returns: Option<Curve2>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn between_lengths(&self, l0: f64, l1: f64) -> Option<Curve2> {
        // If either the distance between l1 and l0 are less than the curve tolerance or the orders
        // are inverted when the curve isn't closed, we have a poorly conditioned request and we
        // can return None

        let start = self.at_length(l0)?;
        let end = self.at_length(l1)?;
        let mut wrap = end.length_along() < start.length_along();

        let last_index = if self.is_closed {
            self.count() - 2
        } else {
            self.count() - 1
        };

        if (l1 - l0).abs() < self.tol || (!self.is_closed && wrap) {
            None
        } else {
            let mut points = Vec::new();
            let mut working = start;

            loop {
                points.push(working.point);

                // Advance to the next index
                let next_index = working.index + 1;
                if next_index > last_index {
                    // Terminal condition if we're not wrapping, otherwise we go to the beginning
                    if !wrap {
                        break;
                    } else {
                        wrap = false;
                        working = self.at_front();
                    }
                } else if working.length_along() <= end.length_along() && next_index > end.index {
                    break;
                } else {
                    working = self.at_vertex(next_index);
                }
            }

            if dist(&end.point, points.last().unwrap()) > self.tol {
                points.push(end.point);
            }

            Curve2::from_points(&points, self.tol, false).ok()
        }
    }

    /// If the curve is open, this will attempt to split it into two parts with the division at the
    /// specified length from the curve start.  If the curve is closed, or if this results in a
    /// segment which has fewer than 2 points, this will return an error.
    ///
    /// # Arguments
    ///
    /// * `length`: The length along the curve at which to split
    ///
    /// returns: Result<(Curve2, Curve2), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::geom2::{Curve2, Point2};
    /// let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0), Point2::new(1.0, 2.0)];
    /// let curve = Curve2::from_points(&points, 1e-6, false).unwrap();
    /// let (a, b) = curve.split_open_at_length(2.0).unwrap();
    ///
    /// assert_relative_eq!(a.length(), 2.0, epsilon = 1e-6);
    /// assert_relative_eq!(b.length(), 1.0, epsilon = 1e-6);
    /// ```
    pub fn split_open_at_length(&self, length: f64) -> crate::Result<(Self, Self)> {
        if self.is_closed {
            Err("Cannot split_open_at_length a closed curve"
                .to_string()
                .into())
        } else {
            let a = self.between_lengths(0.0, length).ok_or(format!(
                "Failed to extract curve start 0.0->{} of {}",
                length,
                self.length()
            ))?;
            let b = self.between_lengths(length, self.length()).ok_or(format!(
                "Failed to extract curve end {}->{}",
                length,
                self.length()
            ))?;

            Ok((a, b))
        }
    }

    /// If the curve is closed, this will attempt to split it into two parts with the divisions at
    /// the specified lengths from the curve start.  If the curve is open, or if this results in a
    /// segment which has fewer than 2 points, this will return an error.
    ///
    /// # Arguments
    ///
    /// * `length0`: the first length along the curve at which to split
    /// * `length1`: the second length along the curve at which to split
    ///
    /// returns: Result<(<unknown>, <unknown>), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::geom2::{Curve2, Point2};
    /// let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0), Point2::new(1.0, 2.0),
    ///                   Point2::new(0.0, 2.0)];
    /// let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
    /// let (a, b) = curve.split_closed_at_lengths(2.0, 4.0).unwrap();
    ///
    /// assert_relative_eq!(a.length(), 2.0, epsilon = 1e-6);
    /// assert_relative_eq!(b.length(), 4.0, epsilon = 1e-6);
    /// ```
    pub fn split_closed_at_lengths(
        &self,
        length0: f64,
        length1: f64,
    ) -> crate::Result<(Self, Self)> {
        if !self.is_closed {
            Err("Cannot split_closed_at_lengths an open curve"
                .to_string()
                .into())
        } else {
            let a = self.between_lengths(length0, length1).ok_or(format!(
                "Failed to extract curve between {}->{}",
                length0, length1
            ))?;
            let b = self.between_lengths(length1, length0).ok_or(format!(
                "Failed to extract curve between {}->{}",
                length1, length0
            ))?;

            Ok((a, b))
        }
    }

    pub fn split_across_line(&self, line: &impl LineOps2) -> SplitResult<Vec<Self>> {
        todo!()
    }

    /// Trim a specified amount of length off of the curve's front, returning a new curve if the
    /// operation is successful.
    ///
    /// # Arguments
    ///
    /// * `length`: The amount of length to remove from the front of the curve
    ///
    /// returns: Option<Curve2>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn trim_front(&self, length: f64) -> Option<Curve2> {
        self.between_lengths(length, self.length())
    }

    /// Trim a specified amount of length off of the curve's back, returning a new curve if the
    /// operation is successful.
    ///
    /// # Arguments
    ///
    /// * `length`: the amount of length to remove from the back of the curve
    ///
    /// returns: Option<Curve2>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn trim_back(&self, length: f64) -> Option<Curve2> {
        self.between_lengths(0.0, self.length() - length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Vector2;
    use approx::assert_relative_eq;

    use test_case::test_case;

    use rand::distr::Uniform;
    use rand::prelude::Distribution;
    use rand::rng;

    fn sample1() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    }

    fn sample2() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    }

    fn sample_points(p: &[(f64, f64)]) -> Vec<Point2> {
        p.iter().map(|(a, b)| Point2::new(*a, *b)).collect()
    }

    fn sample_points_scaled(p: &[(f64, f64)], f: f64) -> Vec<Point2> {
        p.iter().map(|(a, b)| Point2::new(*a * f, *b * f)).collect()
    }

    #[test]
    fn stress_between_lengths_by_control() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let mut rn = rng();
        let dist = Uniform::new(0.0, curve.length()).unwrap();

        for _ in 0..5000 {
            let a = dist.sample(&mut rn);
            let mut b = dist.sample(&mut rn);
            while (a - b).abs() < 1e-6 {
                b = dist.sample(&mut rn);
            }

            let mut c = dist.sample(&mut rn);
            while (a - c).abs() < 1e-6 || (b - c).abs() < 1e-6 {
                c = dist.sample(&mut rn);
            }

            let _p_a = curve.at_length(a).unwrap().point;
            let _p_b = curve.at_length(b).unwrap().point;
            let p_c = curve.at_length(c).unwrap().point;

            let segment = curve.between_lengths_by_control(a, b, c);
            assert!(segment.is_some());
            let segment = segment.unwrap();

            let cp = segment.at_closest_to_point(&p_c).point();
            assert_relative_eq!((cp - p_c).norm(), 0.0, epsilon = 1e-6);
        }
    }

    #[test_case((0.1, 1.2), false, vec![1])] //             (0) |->  (1)  ->| (2)      (3)      O/C
    #[test_case((0.1, 2.2), false, vec![1, 2])] //          (0) |->  (1)  ->  (2)  ->| (3)      O/C
    #[test_case((0.7, 0.2), true, vec![1, 2, 3, 0])] //     (0)->||->(1)  ->  (2)  ->  (3)      C
    #[test_case((1.7, 1.2), true, vec![2, 3, 0, 1])] //     (0)  ->  (1)->||->(2)  ->  (3)      C
    #[test_case((2.7, 2.2), true, vec![3, 0, 1, 2])] //     (0)  ->  (1)  ->  (2)->||->(3) ->   C
    #[test_case((3.7, 3.2), true, vec![0, 1, 2, 3])] //     (0)  ->  (1)  ->  (2)  ->  (3)->||->C
    #[test_case((1.2, 0.7), true, vec![2, 3, 0])] //        (0)  ->| (1) |->  (2)  ->  (3) ->   C
    #[test_case((3.2, 0.7), true, vec![0])] //              (0)  ->| (1)      (2)      (3) ->|  C
    #[test_case((0.2, 3.7), true, vec![1, 2, 3])] //        (0) |->  (1)  ->  (2)  ->  (3) ->|  C
    #[test_case((0.1, 0.2), false, Vec::<usize>::new())] // (0) |->| (1)      (2)      (3)     O/C
    #[test_case((0.1, 0.2), true, Vec::<usize>::new())] //  (0) |->| (1)      (2)      (3)     O/C
    #[test_case((1.1, 1.8), false, Vec::<usize>::new())] // (0)      (1) |->| (2)      (3)     O/C
    #[test_case((1.1, 1.8), true, Vec::<usize>::new())] //  (0)      (1) |->| (2)      (3)     O/C
    #[test_case((3.1, 3.8), true, Vec::<usize>::new())] //  (0)      (1)      (2)      (3)|->| C
    fn portioning(l: (f64, f64), c: bool, i: Vec<usize>) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, c).unwrap();
        let p0 = curve.at_length(l.0).unwrap().point;
        let p1 = curve.at_length(l.1).unwrap().point;
        let result = curve.between_lengths(l.0, l.1).unwrap();

        let e_l = if l.1 > l.0 {
            l.1 - l.0
        } else {
            curve.length() - (l.0 - l.1)
        };

        assert_relative_eq!(e_l, result.length(), epsilon = result.tol);

        let first = result.at_front();
        let last = result.at_back();
        assert_relative_eq!(p0.x, first.point.x, epsilon = result.tol);
        assert_relative_eq!(p0.y, first.point.y, epsilon = result.tol);
        assert_relative_eq!(p1.x, last.point.x, epsilon = result.tol);
        assert_relative_eq!(p1.y, last.point.y, epsilon = result.tol);

        for index in i {
            assert!(has_vertex(&points[index], result.shape.vertices()));
        }
    }

    fn has_vertex(v: &Point2, c: &[Point2]) -> bool {
        for t in c.iter() {
            if dist(t, v) < 1e-6 {
                return true;
            }
        }
        false
    }
}
