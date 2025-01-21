use super::polyline2::{polyline_intersections, spanning_ray, SpanningRay};
use crate::common::points::{
    dist, max_point_in_direction, ramer_douglas_peucker, transform_points,
};
use crate::common::{Intersection, Resample};
use crate::errors::InvalidGeometry;
use crate::geom2::hull::convex_hull_2d;
use crate::geom2::line2::Segment2;
use crate::geom2::{Iso2, Point2, SurfacePoint2, UnitVec2};
use crate::{Arc2, Circle2, Result, Series1, Vector2};
use parry2d_f64::na::Unit;
use parry2d_f64::query::{PointQueryWithLocation, Ray};
use parry2d_f64::shape::{ConvexPolygon, Polyline};

/// A `CurveStation2` is a convenience struct which represents a location on the manifold defined
/// by the curve. It has a point, a direction, and a normal. It has an index and a fraction which
/// represent where the location is in the underlying data structure, and it maintains a
/// reference back to the original curve for computing other properties efficiently.
#[derive(Copy, Clone)]
pub struct CurveStation2<'a> {
    /// The point in space where this station is located
    point: Point2,

    /// The direction of the curve at this station. If the station is on an edge of the polyline,
    /// it points in the direction of the edge. If the station is at a vertex, it points in the
    /// averaged direction of the two edges that meet at the vertex.
    direction: UnitVec2,

    /// The index of the vertex in the underlying polyline which directly precedes this station
    index: usize,

    /// The fraction of the distance between the vertex at index and the next vertex at which this
    /// station is located.  This is an alternate way of representing the location of the station
    /// in the manifold's space which is more convenient for some operations.
    fraction: f64,

    /// A reference back to the original curve for computing other properties efficiently.
    curve: &'a Curve2,
}

impl<'a> CurveStation2<'a> {
    fn new(
        point: Point2,
        direction: UnitVec2,
        index: usize,
        fraction: f64,
        curve: &'a Curve2,
    ) -> Self {
        Self {
            point,
            direction,
            index,
            fraction,
            curve,
        }
    }

    /// Get the curve station at the vertex which is or directly precedes this station. If this
    /// station is already at the vertex (fraction == 0.0), this will return an identical station.
    pub fn at_index(&self) -> Self {
        self.curve.at_vertex(self.index)
    }

    /// Get the curve station at the vertex which directly follows this station. If this station is
    /// at the last vertex, this will return an identical station.
    pub fn at_next_index(&self) -> Self {
        if self.index == self.curve.count() - 1 {
            self.curve.at_back()
        } else {
            self.curve.at_vertex(self.index + 1)
        }
    }

    /// Get the curve station at the vertex which directly precedes this station. If this station
    /// is beyond its vertex (fraction > 0.0), this will return the same station as `at_index()`,
    /// but if the station is directly on the vertex (fraction == 0.0), this will return the
    /// station at the vertex before this one.  If it is at the first vertex, this will return
    /// `None`.
    pub fn previous(&self) -> Option<Self> {
        // TODO: needs tests
        if self.fraction > 0.0 {
            Some(self.curve.at_vertex(self.index))
        } else if self.index > 0 {
            Some(self.curve.at_vertex(self.index - 1))
        } else {
            None
        }
    }

    /// Get the curve station at the vertex which directly follows this station. If this station is
    /// beyond its vertex (fraction < 1.0), this will return the same station as `at_next_index()`,
    /// but if the station is already at the end of the curve it will return `None`.
    pub fn next(&self) -> Option<Self> {
        if self.fraction < 1.0 {
            if self.index < self.curve.count() - 1 {
                Some(self.curve.at_vertex(self.index + 1))
            } else {
                Some(self.curve.at_back())
            }
        } else if self.index == self.curve.count() - 2 {
            Some(self.curve.at_back())
        } else {
            None
        }
    }

    /// Returns the point in space where this station is located
    pub fn point(&self) -> Point2 {
        self.point
    }

    /// Returns the direction of the curve at this station. If the station is on an edge of the
    /// polyline, it points in the direction of the edge. If the station is at a vertex, it points
    /// in the averaged direction of the two edges that meet at the vertex.
    pub fn direction(&self) -> UnitVec2 {
        self.direction
    }

    /// Returns the normal of the curve at this station. This is the direction of the curve rotated
    /// -90 degrees
    pub fn normal(&self) -> UnitVec2 {
        let t = Iso2::rotation(-std::f64::consts::FRAC_PI_2);
        t * self.direction
    }

    /// Returns the index of the vertex in the underlying polyline which directly precedes this
    /// station
    pub fn index(&self) -> usize {
        self.index
    }

    /// Returns the fraction of the distance between the vertex at index and the next vertex at
    /// which this station is located.  This is an alternate way of representing the location of
    /// the station in the manifold's space which is more convenient for some operations.
    pub fn fraction(&self) -> f64 {
        self.fraction
    }

    /// Returns the total length of the curve (in world units) up to this station. This would be
    /// the length of the curve if it were cut at this station and then straightened out.
    pub fn length_along(&self) -> f64 {
        let l = &self.curve.lengths;
        l[self.index] + (l[self.index + 1] - l[self.index]) * self.fraction
    }

    /// Create a SurfacePoint2 from this station, where the point is the same as the station's
    /// point, and the normal is the same as the station's normal.
    pub fn surface_point(&self) -> SurfacePoint2 {
        SurfacePoint2::new(self.point, self.normal())
    }

    /// Create a SurfacePoint2 from this station, where the point is the same as the station's
    /// point and the direction is the same as the station's direction
    pub fn direction_point(&self) -> SurfacePoint2 {
        SurfacePoint2::new(self.point, self.direction())
    }

    /// Creates a `SurfacePoint2` from this station similar to `surface_point()`, but where the
    /// normal is the linearly interpolated normal between the previous vertex and the next
    /// vertex based on the fraction of the distance this station is between the two vertices.
    /// This can be used to estimate smoothed curvature along the parent curve.
    pub fn interpolated_surface_point(&self) -> SurfacePoint2 {
        let n0_o = self.previous();
        let n1_o = self.next();

        if let (Some(n0), Some(n1)) = (n0_o, n1_o) {
            let n0 = n0.normal();
            let n1 = n1.normal();
            let n = n0.slerp(&n1, self.fraction);
            SurfacePoint2::new(self.point, n)
        } else {
            self.surface_point()
        }
    }

    /// Creates a `SurfacePoint2` from this station similar to `surface_point()`, but where the
    /// normal is the linearly interpolated direction between the previous vertex and the next
    /// vertex based on the fraction of the distance this station is between the two vertices.
    /// This can be used to estimate smoothed curvature along the parent curve.
    pub fn interpolated_direction_point(&self) -> SurfacePoint2 {
        let n0_o = self.previous();
        let n1_o = self.next();

        if let (Some(n0), Some(n1)) = (n0_o, n1_o) {
            let n0 = n0.direction();
            let n1 = n1.direction();
            let n = n0.slerp(&n1, self.fraction);
            SurfacePoint2::new(self.point, n)
        } else {
            self.surface_point()
        }
    }
}

/// A Curve2 is a 2-dimensional polygonal chain in which its points are connected. It optionally
/// may include normals. This struct and its methods allow for convenient handling of distance
/// searches, transformations, resampling, and splitting.
#[derive(Clone)]
pub struct Curve2 {
    line: Polyline,

    /// These will be the lengths of each vertex along the curve, where the first length is 0.0,
    /// and the last is the total length of the curve.  This is used for efficient length-based
    /// searches.
    lengths: Vec<f64>,
    is_closed: bool,
    tol: f64,
}

impl Curve2 {
    pub fn points(&self) -> &[Point2] {
        self.line.vertices()
    }

    pub fn vtx(&self, i: usize) -> Point2 {
        self.line.vertices()[i].clone()
    }

    /// Builds a Curve2 from a sequence of points. The points will be de-duplicated within the
    /// tolerance.  If the first and last points are within the tolerance *or* the `force_closed`
    /// flag is set the curve will be considered closed.
    pub fn from_points(points: &[Point2], tol: f64, force_closed: bool) -> Result<Self> {
        let mut pts = points.to_vec();
        pts.dedup_by(|a, b| dist(a, b) <= tol);

        if pts.len() < 2 {
            return Err(Box::try_from(InvalidGeometry::NotEnoughPoints).unwrap());
        }

        // Check if the curve is supposed to be closed
        if let (true, Some(start), Some(end)) = (force_closed, pts.first(), pts.last()) {
            if dist(start, end) > tol {
                pts.push(*start);
            }
        }

        let is_closed = pts.len() >= 2 && dist(&pts[0], pts.last().unwrap()) <= tol;

        // Because we will not actually pass indices into the polyline creation method we can
        // trust that the edges will match the vertex indices.  There will be one less edge than
        // there is vertices, and each edge i will join vertex i with vertex i+1
        let line = Polyline::new(pts, None);
        let v = line.vertices();

        let mut lengths: Vec<f64> = vec![0.0];
        for i in 0..v.len() - 1 {
            let d = dist(&v[i + 1], &v[i]);
            lengths.push(d + lengths.last().unwrap_or(&0.0));
        }

        Ok(Curve2 {
            line,
            lengths,
            is_closed,
            tol,
        })
    }

    /// Create a `Curve2` from a sequence of points, but ensure that the curve is oriented such
    /// that the points are in counter-clockwise order compared to their convex hull.  This is
    /// useful for enforcing that a closed curve is oriented so that the normals are pointing
    /// outwards.
    ///
    /// Internally, this works by computing the convex hull and then checking if more points on the
    /// hull have an index greater than or less than the index of their immediate neighbor.
    /// Because the convex hull is always oriented counter-clockwise, ascending indices indicate
    /// that the curve is oriented counter-clockwise, and descending indices indicate that the
    /// curve is oriented clockwise.  If the curve is oriented clockwise, it will be reversed.
    ///
    /// # Arguments
    ///
    /// * `points`: The points to build the curve from, must be in sequence
    /// * `tol`: The general tolerance to use for the curve, used for de-duplication of points
    /// * `force_closed`: If true, the curve will be closed even if the first and last points are
    /// not within the tolerance of each other
    ///
    /// returns: Result<Curve2, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_points_ccw(points: &[Point2], tol: f64, force_closed: bool) -> Result<Self> {
        let mut d_sum = 0;
        let hull = convex_hull_2d(points);
        for i in 0..hull.len() {
            let j = (i + 1) % hull.len();
            let d = hull[j] as i32 - hull[i] as i32;
            d_sum += d.signum();
        }

        if d_sum > 0 {
            Curve2::from_points(points, tol, force_closed)
        } else {
            let points2 = points.iter().rev().copied().collect::<Vec<_>>();
            Curve2::from_points(&points2, tol, force_closed)
        }
    }

    /// Builds a Curve2 from a sequence of SurfacePoints. The points will be de-duplicated within
    /// the tolerance.  The direction of the curve will be such that the majority of the normals
    /// are pointing in the same half-space as the corresponding normals of the curve.  This is not
    /// a guarantee that the normals of the surface points will all match the normals of the curve.
    pub fn from_surf_points(
        points: &[SurfacePoint2],
        tol: f64,
        force_closed: bool,
    ) -> Result<Self> {
        let pts: Vec<Point2> = points.iter().map(|p| p.point).collect();
        let c = Self::from_points(&pts, tol, force_closed)?;

        let mut votes = 0.0;
        for p in points {
            let s = c.at_closest_to_point(&p.point);
            if s.normal().dot(&p.normal) > 0.0 {
                votes += 1.0;
            } else {
                votes -= 1.0;
            }
        }

        if votes < 0.0 {
            Ok(c.reversed())
        } else {
            Ok(c)
        }
    }

    pub fn count(&self) -> usize {
        self.line.vertices().len()
    }

    pub fn lengths(&self) -> &Vec<f64> {
        &self.lengths
    }

    pub fn tol(&self) -> f64 {
        self.tol
    }

    fn dir_of_edge(&self, edge_index: usize) -> UnitVec2 {
        let v0 = self.vtx(edge_index);
        let v1 = self.vtx(edge_index + 1);
        Unit::new_normalize(v1 - v0)
    }

    fn dir_of_vertex(&self, index: usize) -> UnitVec2 {
        let v = self.line.vertices();
        let is_first = index == 0;
        let is_last = index == v.len() - 1;

        if self.is_closed && (is_first || is_last) {
            let d0 = self.dir_of_edge(0).into_inner();
            let d1 = self.dir_of_edge(v.len() - 2).into_inner();
            // TODO: this will fail on a curve that doubles back, use angles?
            Unit::new_normalize(d0 + d1)
        } else if is_first {
            self.dir_of_edge(0)
        } else if is_last {
            self.dir_of_edge(v.len() - 2)
        } else {
            let d0 = self.dir_of_edge(index - 1).into_inner();
            let d1 = self.dir_of_edge(index).into_inner();
            // TODO: this will fail on a curve that doubles back, use angles?
            Unit::new_normalize(d0 + d1)
        }
    }

    fn at_vertex(&self, index: usize) -> CurveStation2 {
        let v = self.line.vertices();
        let (i, f) = if index == v.len() - 1 {
            (index - 1, 1.0)
        } else {
            (index, 0.0)
        };

        CurveStation2::new(self.vtx(index), self.dir_of_vertex(index), i, f, self)
    }

    pub fn at_front(&self) -> CurveStation2 {
        self.at_vertex(0)
    }

    pub fn at_back(&self) -> CurveStation2 {
        self.at_vertex(self.line.vertices().len() - 1)
    }

    pub fn at_length(&self, length: f64) -> Option<CurveStation2> {
        if length < 0.0 || length > self.length() {
            None
        } else {
            let search = self
                .lengths
                .binary_search_by(|a| a.partial_cmp(&length).unwrap());
            match search {
                Ok(index) => Some(self.at_vertex(index)),
                Err(next_index) => {
                    // next_index will be the index of the first element greater than the value we
                    // were searching for, and the first length in the lengths vector will be 0.0,
                    // so we are guaranteed that next_index is greater than 1
                    let index = next_index - 1;
                    let dir = self.dir_of_edge(index);
                    let remaining_len = length - self.lengths[index];
                    let f = remaining_len / (self.lengths[index + 1] - self.lengths[index]);
                    let point = self.vtx(index) + dir.into_inner() * remaining_len;
                    Some(CurveStation2::new(point, dir, index, f, self))
                }
            }
        }
    }

    pub fn at_fraction(&self, fraction: f64) -> Option<CurveStation2> {
        self.at_length(fraction * self.length())
    }

    pub fn at_closest_to_point(&self, test_point: &Point2) -> CurveStation2 {
        let (prj, loc) = self
            .line
            .project_local_point_and_get_location(test_point, false);
        let (edge_index, sp) = loc;
        let dir = self.dir_of_edge(edge_index as usize);

        CurveStation2::new(
            prj.point,
            dir,
            edge_index as usize,
            sp.barycentric_coordinates()[1],
            self,
        )
    }

    pub fn dist_to_point(&self, test_point: &Point2) -> f64 {
        let (prj, _) = self
            .line
            .project_local_point_and_get_location(test_point, false);
        dist(&prj.point, test_point)
    }

    pub fn is_closed(&self) -> bool {
        self.is_closed
    }

    pub fn length(&self) -> f64 {
        *self.lengths.last().unwrap_or(&0.0)
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

            if let Ok(c) = Curve2::from_points(&points, self.tol, false) {
                Some(c)
            } else {
                None
            }
        }
    }

    pub fn reversed(&self) -> Self {
        let mut points = self.clone_points();
        points.reverse();
        Curve2::from_points(&points, self.tol, false).unwrap()
    }

    pub fn make_hull(&self) -> Option<ConvexPolygon> {
        ConvexPolygon::from_convex_hull(self.line.vertices())
    }

    /// Find the point on the curve which is the maximum in the direction of the given vector.
    ///
    /// # Arguments
    ///
    /// * `vector`: The direction vector to search for the maximum point in
    ///
    /// returns: Option<(usize, OPoint<f64, Const<{ D }>>)>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn max_point_in_direction(&self, vector: &Vector2) -> Option<(usize, Point2)> {
        // TODO: there is probably a much more efficient way to do this using the bvh
        max_point_in_direction(&self.points(), &vector)
    }

    /// Find the maximum distance of any point on the curve in the direction of the given surface
    /// point normal, and return that distance. The maximum point is found identically to
    /// `max_point_in_direction()`, and then the distance is computed as the scalar projection of
    /// that maximum point onto the surface point. The result is the component of distance only in
    /// the direction of the normal.
    ///
    /// # Arguments
    ///
    /// * `surf_point`: the point and normal to use for the measurement
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::geom2::{Curve2, Point2, SurfacePoint2, Vector2};
    /// let p = vec![Point2::new(-10.0, 0.0), Point2::new(-5.0, 11.0), Point2::new(-10.0, 12.0)];
    /// let curve = Curve2::from_points(&p, 1e-6, false).unwrap();
    ///
    /// let test = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(1.0, 0.0));
    ///
    /// let d = curve.max_dist_in_direction(&test);
    /// assert_relative_eq!(d, -5.0, epsilon = 1e-6);
    /// ```
    pub fn max_dist_in_direction(&self, surf_point: &SurfacePoint2) -> f64 {
        if let Some((_, p)) = self.max_point_in_direction(&surf_point.normal) {
            surf_point.scalar_projection(&p)
        } else {
            f64::NEG_INFINITY
        }
    }

    /// Clone the points of the curve into a new vector
    pub fn clone_points(&self) -> Vec<Point2> {
        self.line.vertices().to_vec()
    }

    /// Simplify a curve by removing points such that the largest difference between the new curve
    /// and the old curve is less than or equal to the tolerance specified.  This uses the
    /// Ramer-Douglas-Peucker algorithm.
    ///
    /// # Arguments
    ///
    /// * `tol`: The maximum allowable distance between any point on the old curve and its closest
    /// projection onto the new curve
    ///
    /// returns: Curve2
    pub fn simplify(&self, tol: f64) -> Self {
        let new_points = ramer_douglas_peucker(self.line.vertices(), tol);
        Curve2::from_points(&new_points, self.tol, self.is_closed).unwrap()
    }

    pub fn resample(&self, mode: Resample) -> Result<Self> {
        match mode {
            Resample::ByCount(n) => resample_by_count(self, n),
            Resample::BySpacing(l) => resample_by_spacing(self, l),
            Resample::ByMaxSpacing(lm) => resample_by_max_spacing(self, lm),
        }
    }

    pub fn iter(&self) -> Curve2Iterator {
        Curve2Iterator {
            curve: self,
            index: 0,
        }
    }

    /// Create a new curve which is the result of transforming this curve by the given isometry.
    ///
    /// # Arguments
    ///
    /// * `transform`:
    ///
    /// returns: Curve2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn transformed(&self, transform: &Iso2) -> Self {
        let points = transform_points(self.line.vertices(), transform);
        Curve2::from_points(&points, self.tol, self.is_closed).unwrap()
    }

    /// Perform a ray cast against the curve, returning a list of intersections, each as a pair of
    /// values representing the distance along the ray and the index of the edge where the
    /// intersection occurs.
    ///
    /// # Arguments
    ///
    /// * `ray`: The 2D ray to cast against the curve.
    ///
    /// returns: Vec<(f64, usize), Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn ray_intersections(&self, ray: &Ray) -> Vec<(f64, usize)> {
        polyline_intersections(&self.line, ray)
    }

    /// Given a full ray, this method will attempt to create a `SpanningRay` which crosses the
    /// curve.  This will require that the ray intersects the curve at exactly two points,
    /// otherwise `None` will be returned.  The resulting ray, if successful, will have the same
    /// direction as the full ray, but the origin will be the earliest point of intersection with
    /// the curve.
    ///
    /// # Arguments
    ///
    /// * `full_ray`: the full test ray to attempt to span the curve
    ///
    /// returns: Option<SpanningRay>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn try_create_spanning_ray(&self, full_ray: &Ray) -> Option<SpanningRay> {
        spanning_ray(&self.line, full_ray)
    }

    /// Gets the curvature of three points as turning angle per unit length, found by the reciprocal
    /// of the radius of the circle.
    pub fn get_curvature(&self, i0: usize, i1: usize, i2: usize) -> f64 {
        if let Ok(circle) = Circle2::from_3_points(self.vtx(i0), self.vtx(i1), self.vtx(i2)) {
            1.0 / circle.r()
        } else {
            0.0
        }
    }

    /// Get a `Series1` representing the curvature of the curve.  The x values of the series will
    /// be the length along the curve, and the y values will be the turning angle per unit length
    /// at that location. Larger values indicate higher curvature, while a 0 value indicates a
    /// straight section.
    ///
    /// This works by first attempting to calculate the radius of curvature at each vertex, by
    /// taking a three point circle of each vertex and its two neighbors.  If successful, the
    /// curvature is 1/r, otherwise it is 0 (as it would be if r was infinite)
    ///
    /// If the curve is closed, the end points will be treated as
    /// contiguous, otherwise the end points will have the same value as their neighbors.
    pub fn get_curvature_series(&self) -> Series1 {
        let xs = self.lengths.clone();
        let mut ys = vec![0.0; self.count()];

        for i in 1..self.count() - 1 {
            if i == 0 {
                if self.is_closed {
                    ys[i] = self.get_curvature(self.count() - 2, i, i + 1);
                }
            }

            if i == self.count() - 1 {
                if self.is_closed {
                    ys[i] = self.get_curvature(i - 1, i, 1);
                }
            }

            ys[i] = self.get_curvature(i - 1, i, i + 1);
        }

        if !self.is_closed {
            ys[0] = ys[1];
            ys[self.count() - 1] = ys[self.count() - 2];
        }

        Series1::try_new(xs, ys).unwrap()
    }

    /// Searches the curve for regions which have equivalent arcs within the given tolerance and
    /// with at least the given number of points.  The result is a list of tuples, each containing
    /// the start and end indices of the equivalent arc and the arc itself.
    ///
    /// # Arguments
    ///
    /// * `tol`:
    ///
    /// returns: Vec<(usize, usize, Arc2), Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn equivalent_arcs(&self, tol: f64, min_points: usize) -> Vec<(usize, usize, Arc2)> {
        // TODO: handle the wrapping case
        let mut arcs = Vec::new();
        let mut i = 0;

        while i < self.count() {
            if let Some((i0, i1, arc)) = self.equivalent_arc_at(i, tol) {
                if i1 - i0 + 1 >= min_points {
                    arcs.push((i0, i1, arc));
                    i = i1 + 1;
                }
            }

            i += 1;
        }

        arcs
    }

    pub fn equivalent_arc_at(&self, seed_index: usize, tol: f64) -> Option<(usize, usize, Arc2)> {
        // TODO: this is hacked together, fix it

        // Find the initial valid ranges
        let mut pos = 2;
        let mut neg = 0;
        let mut best_arc = None;

        while let Some(arc) = self.in_tol_arc(seed_index, pos, neg, tol) {
            pos += 1;
            best_arc = Some(arc);
        }

        while let Some(arc) = self.in_tol_arc(seed_index, pos, neg, tol) {
            neg += 1;
            best_arc = Some(arc);
        }

        best_arc.map(|arc| (seed_index - neg, seed_index + pos, arc))
    }

    fn in_tol_arc(&self, start: usize, pos: usize, neg: usize, tol: f64) -> Option<Arc2> {
        if neg > start || start + pos >= self.count() {
            return None;
        }

        let i0 = start - neg;
        let i2 = start + pos;
        let i1 = (i0 + i2) / 2;

        let c = Circle2::from_3_points(self.vtx(i0), self.vtx(i1), self.vtx(i2));

        if c.is_err() {
            return None;
        }

        let c = c.unwrap();

        for i in i0..=i2 {
            if c.distance_to(&self.vtx(i)) > tol {
                return None;
            }
        }

        Some(Arc2::three_points(self.vtx(i0), self.vtx(i1), self.vtx(i2)))
    }
}

impl Intersection<&Circle2, Vec<Point2>> for Curve2 {
    /// Computes a set of intersections between this curve and a circle.  The result is a list of
    /// points which are on both the curve and the circle.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: Vec<OPoint<f64, Const<2>>, Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn intersection(&self, other: &Circle2) -> Vec<Point2> {
        let mut points = Vec::new();

        // TODO: We should be able to use the bvh to prune the segments we need to check
        for i in 0..self.count() - 1 {
            if let Ok(seg) = Segment2::try_new(self.vtx(i), self.vtx(i + 1)) {
                for p in other.intersection(&seg) {
                    points.push(p);
                }
            }
        }

        points
    }
}

impl Intersection<&SurfacePoint2, Vec<f64>> for Curve2 {
    /// Generates all intersections between this curve and a surface point, where an intersection
    /// is represented as a distance from the origin of the surface point along the normal of the
    /// surface point where the intersection occurs.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: Vec<f64, Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn intersection(&self, other: &SurfacePoint2) -> Vec<f64> {
        let ray = Ray::new(other.point, other.normal.into_inner());
        self.ray_intersections(&ray)
            .iter()
            .map(|(t, _)| *t)
            .collect()
    }
}

pub struct Curve2Iterator<'a> {
    curve: &'a Curve2,
    index: usize,
}

impl<'a> Iterator for Curve2Iterator<'a> {
    type Item = CurveStation2<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.curve.count() {
            let r = self.curve.at_vertex(self.index);
            self.index += 1;
            Some(r)
        } else {
            None
        }
    }
}

fn resample_by_max_spacing(curve: &Curve2, max_spacing: f64) -> Result<Curve2> {
    let n = (curve.length() / max_spacing).ceil() as usize;
    resample_by_count(curve, n)
}

fn resample_by_spacing(curve: &Curve2, spacing: f64) -> Result<Curve2> {
    let mut positions = Vec::new();
    let mut length = 0.0;
    while length < curve.length() {
        positions.push(length);
        length += spacing;
    }

    let padding = (curve.length() - positions.last().unwrap()) / 2.0;
    for p in &mut positions {
        *p += padding;
    }

    resample_at_positions(curve, &positions)
}

fn resample_by_count(curve: &Curve2, count: usize) -> Result<Curve2> {
    let mut positions = Vec::new();
    for i in 0..count {
        positions.push(i as f64 / (count - 1) as f64);
    }
    resample_at_positions(curve, &positions)
}

fn resample_at_positions(curve: &Curve2, positions: &[f64]) -> Result<Curve2> {
    let mut points = Vec::new();
    for p in positions {
        points.push(curve.at_length(*p).unwrap().point);
    }
    Curve2::from_points(&points, curve.tol, curve.is_closed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Vector2;
    use approx::assert_relative_eq;
    use num_traits::Signed;
    use test_case::test_case;

    use rand::distributions::Uniform;
    use rand::prelude::Distribution;
    use rand::thread_rng;

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
        let mut rng = thread_rng();
        let dist = Uniform::new(0.0, curve.length());

        for _ in 0..5000 {
            let a = dist.sample(&mut rng);
            let mut b = dist.sample(&mut rng);
            while (a - b).abs() < 1e-6 {
                b = dist.sample(&mut rng);
            }

            let mut c = dist.sample(&mut rng);
            while (a - c).abs() < 1e-6 || (b - c).abs() < 1e-6 {
                c = dist.sample(&mut rng);
            }

            let p_a = curve.at_length(a).unwrap().point;
            let p_b = curve.at_length(b).unwrap().point;
            let p_c = curve.at_length(c).unwrap().point;

            let segment = curve.between_lengths_by_control(a, b, c);
            assert!(segment.is_some());
            let segment = segment.unwrap();

            let cp = segment.at_closest_to_point(&p_c).point();
            assert_relative_eq!((cp - p_c).norm(), 0.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_closest_point() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        let p = curve.at_closest_to_point(&Point2::new(2.0, 0.5));
        assert_eq!(p.index, 1);
        assert_relative_eq!(1.0, p.point.x, epsilon = 1e-8);
        assert_relative_eq!(0.5, p.point.y, epsilon = 1e-8);
    }

    #[test]
    fn test_create_open() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();

        assert!(!curve.is_closed());
        assert_relative_eq!(3.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_force_closed() {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }

    #[test]
    fn test_create_naturally_closed() {
        let points = sample_points(&sample2());
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();

        assert!(curve.is_closed());
        assert_relative_eq!(4.0, curve.length(), epsilon = 1e-10);
    }

    #[test_case(0.5, 0, 0.5)]
    #[test_case(0.0, 0, 0.0)]
    #[test_case(2.0, 2, 0.0)]
    #[test_case(2.25, 2, 0.25)]
    fn test_lengths(l: f64, ei: usize, ef: f64) {
        let points = sample_points_scaled(&sample1(), 0.5);
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();

        let r = curve.at_length(l * 0.5).unwrap();
        assert_eq!(ei, r.index);
        assert_relative_eq!(ef, r.fraction, epsilon = 1e-8);
    }

    #[test_case(0.5, (0.5, 0.0))]
    #[test_case(0.0, (0.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.25, (0.75, 1.0))]
    fn test_points_at_length(l: f64, e: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let result = curve.at_length(l).unwrap();

        assert_relative_eq!(e.0, result.point.x, epsilon = 1e-8);
        assert_relative_eq!(e.1, result.point.y, epsilon = 1e-8);
    }

    #[test_case(0.0, (-1.0, -1.0))]
    #[test_case(0.5, (0.0, -1.0))]
    #[test_case(1.0, (1.0, -1.0))]
    #[test_case(1.5, (1.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.5, (0.0, 1.0))]
    #[test_case(3.0, (-1.0, 1.0))]
    #[test_case(3.5, (-1.0, 0.0))]
    #[test_case(4.0, (-1.0, -1.0))]
    fn test_normals_at_length_closed(l: f64, ec: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let e = Unit::new_normalize(Vector2::new(ec.0, ec.1));
        let n = curve.at_length(l).unwrap().normal();

        assert_relative_eq!(e.x, n.x, epsilon = 1e-8);
        assert_relative_eq!(e.y, n.y, epsilon = 1e-8);
    }

    #[test_case(0.0, (0.0, -1.0))]
    #[test_case(0.5, (0.0, -1.0))]
    #[test_case(1.0, (1.0, -1.0))]
    #[test_case(1.5, (1.0, 0.0))]
    #[test_case(2.0, (1.0, 1.0))]
    #[test_case(2.5, (0.0, 1.0))]
    #[test_case(3.0, (0.0, 1.0))]
    fn test_normals_at_length_open(l: f64, ec: (f64, f64)) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, false).unwrap();
        let e = Unit::new_normalize(Vector2::new(ec.0, ec.1));
        let n = curve.at_length(l).unwrap().normal();

        assert_relative_eq!(e.x, n.x, epsilon = 1e-8);
        assert_relative_eq!(e.y, n.y, epsilon = 1e-8);
    }

    fn has_vertex(v: &Point2, c: &[Point2]) -> bool {
        for t in c.iter() {
            if dist(t, v) < 1e-6 {
                return true;
            }
        }
        false
    }

    #[test_case(0.0)]
    #[test_case(0.5)]
    #[test_case(0.75)]
    #[test_case(2.0)]
    #[test_case(2.1)]
    #[test_case(3.9)]
    fn test_distance_along(l: f64) {
        let points = sample_points(&sample1());
        let curve = Curve2::from_points(&points, 1e-6, true).unwrap();
        let p = curve.at_length(l).unwrap().point;
        let d = curve.at_closest_to_point(&p);

        assert_relative_eq!(l, d.length_along(), epsilon = 1e-6);
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
    fn test_portioning(l: (f64, f64), c: bool, i: Vec<usize>) {
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
            assert!(has_vertex(&points[index], result.line.vertices()));
        }
    }
}
