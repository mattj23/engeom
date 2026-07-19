//! This module has tools for partitioning curves into sub-curves.

use crate::common::points::dist;
use crate::common::{Intersection, PCoords};
use crate::geom2::{Aabb2, LineOps2};
use crate::{Circle2, Curve2, Line2};
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

    /// Partitions this curve into sub-curves by a [`CurvePartitioner2`] boundary, returning a
    /// [`SplitResult`] that classifies the resulting pieces.
    ///
    /// The boundary divides 2D space into a *positive* half and a *negative* half via
    /// [`CurvePartitioner2::is_pos`].  Each segment of the curve that lies entirely on one side
    /// becomes one sub-curve.  Wherever the curve crosses the boundary, a precise intersection
    /// point is computed and inserted as the shared endpoint of the two adjacent sub-curves.
    ///
    /// If the curve is closed and the first and last sub-curves end up on the same side of
    /// the boundary, they are spliced together into a single sub-curve so that the topology of
    /// the closed curve is preserved.
    ///
    /// The return value is a [`SplitResult`]:
    /// - `SplitResult::Positive`: the entire curve lies on the positive side (no crossings).
    /// - `SplitResult::Negative`: the entire curve lies on the negative side (no crossings).
    /// - `SplitResult::Pair(negatives, positives)`: the curve crosses the boundary at least
    ///   once; `negatives` holds the sub-curves on the negative side and `positives` holds those
    ///   on the positive side.
    ///
    /// # Built-in partitioners
    ///
    /// [`Line2`], [`Aabb2`], and [`Circle2`] all implement [`CurvePartitioner2`] out of the box.
    /// For `Line2` the positive side is to the right of the line's direction of travel; for
    /// `Aabb2` and `Circle2` the positive side is the exterior of the shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::geom2::{Curve2, Point2};
    /// use engeom::Line2;
    /// use parry3d_f64::query::SplitResult;
    ///
    /// let pts = vec![
    ///     Point2::new(0.0, 0.0), Point2::new(0.3, 0.5),
    ///     Point2::new(0.7, 0.5), Point2::new(1.0, 0.0),
    /// ];
    /// let curve = Curve2::from_points(&pts, 1e-6, false).unwrap();
    /// let divider = Line2::new([0.5, 0.0].into(), [0.0, 1.0].into());
    ///
    /// match curve.partition_by(&divider) {
    ///     SplitResult::Pair(neg, pos) => {
    ///         assert_eq!(neg.len(), 1); // x < 0.5 portion
    ///         assert_eq!(pos.len(), 1); // x > 0.5 portion
    ///     }
    ///     _ => panic!("expected a split"),
    /// }
    /// ```
    pub fn partition_by(&self, boundary: &impl CurvePartitioner2) -> SplitResult<Vec<Self>> {
        // We will work our way from start to end, building the curves based on when and where they
        // cross the partitioning boundary.

        let mut groups = Vec::new();
        let mut working = Vec::new();
        let mut is_pos = boundary.is_pos(&self.points()[0]);

        for s in self.iter() {
            // Add the current point to the working list
            working.push(s.point);

            // Check if we're at the end of the curve.
            if s.fraction > 1.0 - f64::EPSILON {
                // We're at the last point, so we'll end the working group.
                groups.push(working);
                break;
            }

            // Now we'll check if we're about to pass through the partitioning boundary.
            let ns = s.at_next_index();
            if dist(&s, &ns) < f64::EPSILON {
                continue;
            }

            let mut current = Line2::new(s.point, (ns.point - s.point).normalize());
            let mut to_next = current.scalar_project(&ns);

            // Sanity check
            // let s_pos = boundary.is_pos(&s);
            // let ns_pos = boundary.is_pos(&ns);
            // if s_pos != ns_pos {
            //     if !boundary.next_intersection(&current, to_next).is_some() {
            //         println!("Inconsistency: Boundary sign change without intersection")
            //     }
            // }

            while let Some(t) = find_next_crossing(boundary, &current, to_next, is_pos) {
                let end = current.at(t);
                working.push(end);
                groups.push(working);

                is_pos = !is_pos;

                working = vec![end];
                current = current.shifted_origin(t);
                to_next = current.scalar_project(&ns);
            }
        }

        // Build them into curves
        let mut curves = groups
            .into_iter()
            .filter_map(|group| Curve2::from_points(&group, self.tol, false).ok())
            .collect::<Vec<_>>();

        // If we have a closed curve, we may need to splice together the last and first groups
        let finalized = if curves.len() > 1 && self.is_closed {
            let last_pos = boundary.is_pos(&curves.last().unwrap().at_fraction(0.5).unwrap());
            let first_pos = boundary.is_pos(&curves.first().unwrap().at_fraction(0.5).unwrap());

            if last_pos == first_pos {
                let last = curves.pop().unwrap();
                let first = curves.remove(0);
                let combined_points = [last.points(), first.points()].concat();
                let combined = Curve2::from_points(&combined_points, self.tol, false).unwrap();
                curves.push(combined);
            }

            curves
        } else {
            curves
        };

        let mut positives = Vec::new();
        let mut negatives = Vec::new();
        for curve in finalized {
            let is_pos = boundary.is_pos(&curve.at_fraction(0.5).unwrap());
            if is_pos {
                positives.push(curve);
            } else {
                negatives.push(curve);
            }
        }

        match (negatives.len(), positives.len()) {
            (0, 0) => panic!("No curves found, should not be possible!"),
            (0, _) => SplitResult::Positive,
            (_, 0) => SplitResult::Negative,
            (_, _) => SplitResult::Pair(negatives, positives),
        }
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

/// Given a curve and a line, find the next crossing, where the is_pos value changes. The value
/// will lie between 0.0 and 1.0. The provided line should not be normalized.
fn find_next_crossing(
    bnd: &impl CurvePartitioner2,
    line: &Line2,
    max_dist: f64,
    current_pos: bool,
) -> Option<f64> {
    // The reason that this function exists is that the concept of an intersection does not
    // necessarily correspond with a crossing of the boundary. Rather, crossings occur next to the
    // intersection, but the intersection itself is usually contained within the boundary. This
    // becomes a problem when endpoints are at the boundary.
    let mut candidates = bnd
        .all_intersections(line)
        .into_iter()
        .filter(|t| *t >= 0.0 && *t <= max_dist)
        .collect::<Vec<_>>();

    if candidates.is_empty() {
        return None;
    }

    candidates.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Now we want to find the first intersection that ends on the opposite side of the boundary.
    for t in candidates.iter() {
        let (_, end_pos) = find_crossing_info(bnd, line, *t);
        if end_pos != current_pos {
            return Some(*t);
        }
    }

    // If we didn't find one, there is no crossing
    None
}

fn find_crossing_info(bnd: &impl CurvePartitioner2, line: &Line2, t: f64) -> (bool, bool) {
    let at_t = bnd.is_pos(&line.at(t));
    let mut delta = 1e-10;

    let mut at_t_plus = bnd.is_pos(&line.at(t + delta));
    let mut at_t_minus = bnd.is_pos(&line.at(t - delta));

    // One of these should be different from at_t
    while at_t == at_t_plus && at_t == at_t_minus {
        delta *= 2.0;
        at_t_plus = bnd.is_pos(&line.at(t + delta));
        at_t_minus = bnd.is_pos(&line.at(t - delta));

        if delta > 1e-6 {
            // This intersection isn't going to result in a crossing.
            return (at_t, at_t);
        }
    }

    if at_t != at_t_plus {
        (at_t, at_t_plus)
    } else {
        (at_t_minus, at_t_plus)
    }
}

pub trait CurvePartitioner2 {
    fn is_pos(&self, point: &impl PCoords<2>) -> bool;

    /// Provide all intersections between the entity and the line. The results should be presented
    /// as a vector of t-values for the intersection points along `line`.
    ///
    /// # Arguments
    ///
    /// * `line`: The line, provided by the caller.
    ///
    /// returns: Vec<f64, Global>
    fn all_intersections(&self, line: &Line2) -> Vec<f64>;
}

impl<T: CurvePartitioner2> CurvePartitioner2 for &T {
    fn is_pos(&self, point: &impl PCoords<2>) -> bool {
        (**self).is_pos(point)
    }

    fn all_intersections(&self, line: &Line2) -> Vec<f64> {
        (**self).all_intersections(line)
    }
}

impl CurvePartitioner2 for Aabb2 {
    fn is_pos(&self, point: &impl PCoords<2>) -> bool {
        !self.contains_local_point(&point.coords().into())
    }

    fn all_intersections(&self, line: &Line2) -> Vec<f64> {
        line.intersection(*self)
    }
}

impl CurvePartitioner2 for Line2 {
    fn is_pos(&self, point: &impl PCoords<2>) -> bool {
        self.signed_projection_dist(point).is_sign_positive()
    }

    fn all_intersections(&self, line: &Line2) -> Vec<f64> {
        if let Some((t0, _t1)) = line.intersection_params(self) {
            vec![t0]
        } else {
            Vec::new()
        }
    }
}

impl CurvePartitioner2 for Circle2 {
    fn is_pos(&self, point: &impl PCoords<2>) -> bool {
        !self.contains_point(point)
    }

    fn all_intersections(&self, line: &Line2) -> Vec<f64> {
        self.intersection(line)
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::*;
    use super::*;
    use approx::assert_relative_eq;

    use test_case::test_case;

    use crate::{Line2, Point2};
    use rand::distr::Uniform;
    use rand::prelude::Distribution;
    use rand::rng;

    #[test]
    fn split_on_circle() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, false).unwrap();
        let circle = Circle2::new(0.5, 0.5, 0.6);

        match curve.partition_by(&circle) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 3);
                assert_eq!(positives.len(), 4);
            }
        }
    }

    #[test]
    fn split_on_line_open() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, false).unwrap();
        let line = Line2::new([0.5, 0.0].into(), [0.0, 1.0].into());

        match curve.partition_by(&line) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 2);
                assert_eq!(positives.len(), 1);
            }
        }
    }

    #[test]
    fn split_on_line_closed() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, true).unwrap();
        let line = Line2::new([0.5, 0.0].into(), [0.0, 1.0].into());

        match curve.partition_by(&line) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 1);
                assert_eq!(positives.len(), 1);
            }
        }
    }

    #[test]
    fn split_on_box_open_end() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, false).unwrap();
        let aabb = Aabb2::new([-0.5, 0.5].into(), [0.5, 1.5].into());

        match curve.partition_by(&aabb) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 1);
                assert_eq!(positives.len(), 1);
            }
        }
    }

    #[test]
    fn split_on_box_open_center() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, false).unwrap();
        let aabb = Aabb2::new([0.25, -1.0].into(), [0.75, 1.5].into());

        match curve.partition_by(&aabb) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 2);
                assert_eq!(positives.len(), 3);
            }
        }
    }

    #[test]
    fn split_on_box_closed_center() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, true).unwrap();
        let aabb = Aabb2::new([0.25, -1.0].into(), [0.75, 1.5].into());

        match curve.partition_by(&aabb) {
            SplitResult::Positive => assert!(false, "Should not be positive"),
            SplitResult::Negative => assert!(false, "Should not be negative"),
            SplitResult::Pair(negatives, positives) => {
                assert_eq!(negatives.len(), 2);
                assert_eq!(positives.len(), 2);
            }
        }
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
