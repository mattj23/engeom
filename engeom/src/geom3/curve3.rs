use crate::common::To2D;
use crate::common::domain_window::DomainWindowIter;
use crate::common::points::{dist, ramer_douglas_peucker};
use crate::errors::InvalidGeometry;
use crate::geom2::Curve2;
use crate::geom3::{Aabb3, Iso3, Plane3, Point3, UnitVec3};
use crate::{Func1, Polynomial, Resample, Result, Smoothing, SurfacePoint3, SvdBasis3};
use parry3d_f64::na::Unit;
use parry3d_f64::query::PointQueryWithLocation;
use parry3d_f64::shape::Polyline;

#[derive(Copy, Clone)]
pub struct CurveStation3<'a> {
    point: Point3,

    direction: UnitVec3,

    index: usize,

    fraction: f64,

    curve: &'a Curve3,
}

impl<'a> CurveStation3<'a> {
    fn new(
        point: Point3,
        direction: UnitVec3,
        index: usize,
        fraction: f64,
        curve: &'a Curve3,
    ) -> Self {
        Self {
            point,
            direction,
            index,
            fraction,
            curve,
        }
    }

    pub fn point(&self) -> Point3 {
        self.point
    }

    pub fn direction(&self) -> UnitVec3 {
        self.direction
    }

    /// Returns a SurfacePoint3 at the same position in 3d space as the station, but with a normal
    /// pointing in the direction of the next vertex on the curve.
    pub fn direction_point(&self) -> SurfacePoint3 {
        SurfacePoint3::new(self.point, self.direction)
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn fraction(&self) -> f64 {
        self.fraction
    }

    pub fn curve(&self) -> &'a Curve3 {
        self.curve
    }

    pub fn at_index(&self) -> Self {
        self.curve.at_vertex(self.index)
    }

    pub fn length_along(&self) -> f64 {
        let l = self.curve.lengths();
        l[self.index] + (l[self.index + 1] - l[self.index]) * self.fraction
    }

    pub fn plane(&self) -> Plane3 {
        Plane3::from_point_normal(&self.point, &self.direction)
    }

    pub fn is_front(&self) -> bool {
        self.length_along() < self.curve.tol
    }

    pub fn is_back(&self) -> bool {
        self.length_along() > self.curve.length() - self.curve.tol
    }
}

#[derive(Clone)]
pub struct Curve3 {
    line: Polyline,
    lengths: Vec<f64>,
    tol: f64,
}

impl Curve3 {
    pub fn line(&self) -> &Polyline {
        &self.line
    }

    pub fn vtx(&self, i: usize) -> Point3 {
        self.line.vertices()[i]
    }

    /// Returns the axis-aligned bounding box enclosing all curve vertices.
    pub fn aabb(&self) -> Aabb3 {
        self.line.local_aabb()
    }

    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        let points = self
            .line
            .vertices()
            .iter()
            .map(|p| iso * p)
            .collect::<Vec<_>>();
        Self::from_points(&points, self.tol).unwrap()
    }

    /// Projects this curve onto a plane and returns the result as a two-dimensional curve,
    /// expressed in that plane's own coordinate frame.
    ///
    /// This is the general form of [`crate::common::To2D::to_2d`], which only ever drops the z
    /// coordinate. Because [`Plane3::compute_frame`] is built so that the x-y plane's frame is the
    /// identity, projecting onto `Plane3::xy()` here is exactly that same operation, and the two
    /// agree vertex for vertex.
    ///
    /// The curve's chordal tolerance carries over, and the result is closed when the projected
    /// first and last vertices land within it of each other. A section of a mesh therefore comes
    /// back closed if it was a loop, since such a curve repeats its first vertex as its last.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to project onto, whose frame gives the result its coordinates
    ///
    /// returns: Result<Curve2, Box<dyn Error, Global>>
    ///
    /// # Failure
    ///
    /// The projection can bring vertices closer together than the tolerance, and the
    /// de-duplication that follows may leave fewer than two distinct points. A curve running along
    /// the plane normal collapses to a point this way and is an error rather than a degenerate
    /// curve.
    pub fn to_2d_in_plane(&self, plane: &Plane3) -> Result<Curve2> {
        let to_plane = plane.compute_frame().inverse();
        let points = self
            .vertices()
            .iter()
            .map(|p| (to_plane * p).to_2d())
            .collect::<Vec<_>>();

        Curve2::from_points(&points, self.tol, false)
    }

    pub fn from_points(points: &[Point3], tol: f64) -> Result<Self> {
        let mut points = points.to_vec();
        points.dedup_by(|a, b| dist(a, b) <= tol);
        if points.len() < 2 {
            return Err(Box::from(InvalidGeometry::NotEnoughPoints));
        }

        let line = Polyline::new(points, None);
        let v = line.vertices();

        let mut lengths = vec![0.0];
        for i in 0..v.len() - 1 {
            let d = dist(&v[i + 1], &v[i]);
            lengths.push(lengths[i] + d);
        }

        Ok(Self { line, lengths, tol })
    }

    pub fn count(&self) -> usize {
        self.line.vertices().len()
    }

    pub fn lengths(&self) -> &[f64] {
        &self.lengths
    }

    pub fn length(&self) -> f64 {
        *self.lengths.last().unwrap()
    }

    fn dir_of_edge(&self, edge_index: usize) -> UnitVec3 {
        Unit::new_normalize(self.vtx(edge_index + 1) - self.vtx(edge_index))
    }

    fn dir_of_vertex(&self, index: usize) -> UnitVec3 {
        if index == self.line.vertices().len() - 1 {
            self.dir_of_edge(index - 1)
        } else {
            self.dir_of_edge(index)
        }
    }

    fn at_vertex(&self, index: usize) -> CurveStation3<'_> {
        let (i, f) = if index == self.line.vertices().len() - 1 {
            (index - 1, 1.0)
        } else {
            (index, 0.0)
        };

        CurveStation3::new(
            self.line.vertices()[index],
            self.dir_of_vertex(index),
            i,
            f,
            self,
        )
    }

    pub fn at_length(&self, length: f64) -> Option<CurveStation3<'_>> {
        if length < 0.0 || length > self.length() {
            None
        } else {
            let search = self
                .lengths
                .binary_search_by(|l| l.partial_cmp(&length).unwrap());
            match search {
                Ok(index) => Some(self.at_vertex(index)),
                Err(next_index) => {
                    let index = next_index - 1;
                    let dir = self.dir_of_edge(index);
                    let remaining = length - self.lengths[index];
                    let f = remaining / (self.lengths[index + 1] - self.lengths[index]);
                    let point = self.vtx(index) + dir.into_inner() * remaining;
                    Some(CurveStation3::new(point, dir, index, f, self))
                }
            }
        }
    }

    pub fn at_fraction(&self, fraction: f64) -> Option<CurveStation3<'_>> {
        self.at_length(fraction * self.length())
    }

    pub fn at_closest_to_point(&self, test_point: &Point3) -> CurveStation3<'_> {
        let (prj, loc) = self
            .line
            .project_local_point_and_get_location(test_point, false);
        let (edge_index, sp) = loc;
        CurveStation3::new(
            prj.point,
            self.dir_of_edge(edge_index as usize),
            edge_index as usize,
            sp.barycentric_coordinates()[1],
            self,
        )
    }

    pub fn resample(&self, mode: Resample) -> Self {
        match mode {
            Resample::ByCount(n) => resample_by_count(self, n),
            Resample::BySpacing(l) => resample_by_spacing(self, l),
            Resample::ByMaxSpacing(lm) => resample_by_max_spacing(self, lm),
        }
    }

    pub fn dist_to_point(&self, test_point: &Point3) -> f64 {
        let (prj, _) = self
            .line
            .project_local_point_and_get_location(test_point, false);
        dist(test_point, &prj.point)
    }

    pub fn tol(&self) -> f64 {
        self.tol
    }

    pub fn vertices(&self) -> &[Point3] {
        self.line.vertices()
    }

    pub fn clone_points(&self) -> Vec<Point3> {
        self.line.vertices().to_vec()
    }

    pub fn at_front(&self) -> CurveStation3<'_> {
        self.at_vertex(0)
    }

    pub fn at_back(&self) -> CurveStation3<'_> {
        self.at_vertex(self.line.vertices().len() - 1)
    }

    pub fn iter(&self) -> Curve3Iterator<'_> {
        Curve3Iterator {
            curve: self,
            index: 0,
        }
    }

    pub fn simplify(&self, tol: f64) -> Self {
        let new_points = ramer_douglas_peucker(self.line.vertices(), tol);
        Self::from_points(&new_points, tol).unwrap()
    }

    pub fn smoothed(&self, mode: Smoothing) -> Result<Self> {
        match mode {
            Smoothing::Gaussian(_sigma) => {
                todo!()
            }
            Smoothing::Quadratic(window) => smooth_by_polynomial::<3>(self, window),
            Smoothing::Cubic(window) => smooth_by_polynomial::<4>(self, window),
        }
    }

    pub fn window_iter(&self, window_size: f64) -> DomainWindowIter<'_> {
        DomainWindowIter::new(self.lengths(), window_size)
    }
}

fn smooth_by_polynomial<const D: usize>(curve: &Curve3, window_size: f64) -> Result<Curve3> {
    let mut new_points = Vec::new();
    for window in curve.window_iter(window_size) {
        let points = window
            .iter()
            .map(|i| curve.line.vertices()[i])
            .collect::<Vec<_>>();

        if points.len() < 3 {
            new_points.push(curve.line.vertices()[window.index]);
            continue;
        }

        let svd = SvdBasis3::from_points(&points, None).ok_or("Failed to compute SVD basis")?;
        let svd_t = Iso3::from(&svd);

        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for p in &points {
            let pt = svd_t * p;
            xs.push(pt.x);
            ys.push(pt.y);
        }

        let Some(poly) = Polynomial::<D>::least_squares(&xs, &ys, None) else {
            return Err("Failed to fit polynomial".into());
        };

        let x = (svd_t * curve.line.vertices()[window.index]).x;
        let y = poly.f(x);
        new_points.push(svd_t.inverse() * Point3::new(x, y, 0.0));
    }

    Curve3::from_points(&new_points, curve.tol())
}

fn resample_by_max_spacing(curve: &Curve3, max_spacing: f64) -> Curve3 {
    // Find the number of points it will take to ensure that the spacing is less than the max
    // spacing
    let n = (curve.length() / max_spacing).ceil() as usize;
    resample_by_count(curve, n)
}

fn resample_by_spacing(curve: &Curve3, spacing: f64) -> Curve3 {
    let mut positions = Vec::new();
    let mut length = 0.0;
    while length < curve.length() {
        positions.push(length);
        length += spacing;
    }

    // Center
    let padding = (curve.length() - positions.last().unwrap()) / 2.0;
    for p in &mut positions {
        *p += padding;
    }

    resample_at_positions(curve, &positions)
}

fn resample_by_count(curve: &Curve3, n: usize) -> Curve3 {
    let mut positions = Vec::new();
    for i in 0..n {
        let f = i as f64 / (n - 1) as f64;
        positions.push(f * curve.length());
    }

    resample_at_positions(curve, &positions)
}

fn resample_at_positions(curve: &Curve3, positions: &[f64]) -> Curve3 {
    let mut points = Vec::new();
    for p in positions {
        let station = curve.at_length(*p).unwrap();
        points.push(station.point());
    }
    Curve3::from_points(&points, curve.tol()).unwrap()
}

pub struct Curve3Iterator<'a> {
    curve: &'a Curve3,
    index: usize,
}

impl<'a> Iterator for Curve3Iterator<'a> {
    type Item = CurveStation3<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.curve.line.vertices().len() {
            let result = self.curve.at_vertex(self.index);
            self.index += 1;
            Some(result)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Vector3;
    use approx::assert_relative_eq;

    /// A closed unit square in the plane `z = z0`, with its first vertex repeated as its last,
    /// which is how a section of a mesh represents a loop.
    fn square_at(z0: f64) -> Curve3 {
        let points = [
            Point3::new(0.0, 0.0, z0),
            Point3::new(1.0, 0.0, z0),
            Point3::new(1.0, 1.0, z0),
            Point3::new(0.0, 1.0, z0),
            Point3::new(0.0, 0.0, z0),
        ];
        Curve3::from_points(&points, 1e-8).unwrap()
    }

    /// The x-y plane's frame is the identity, so this projection has to agree with `to_2d` exactly
    /// rather than merely closely. If these ever diverge, one of the two is wrong.
    #[test]
    fn projecting_onto_the_xy_plane_is_the_same_as_dropping_z() {
        let curve = square_at(3.0);

        let dropped = curve.to_2d().unwrap();
        let projected = curve.to_2d_in_plane(&Plane3::xy()).unwrap();

        assert_eq!(dropped.points().len(), projected.points().len());
        for (a, b) in dropped.points().iter().zip(projected.points().iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
    }

    /// A curve lying in the plane it is projected onto keeps its shape, so its length is unchanged.
    #[test]
    fn projecting_onto_its_own_plane_preserves_length() {
        let curve = square_at(0.0);
        let plane = Plane3::from_point_normal(
            &Point3::new(0.0, 0.0, 0.0),
            &UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
        );

        let projected = curve.to_2d_in_plane(&plane).unwrap();
        assert_relative_eq!(projected.length(), curve.length(), epsilon = 1e-12);
    }

    /// A loop has to come back closed, which is the property the multi-curve alignment depends on
    /// and the reason sections repeat their first vertex.
    #[test]
    fn a_loop_projects_to_a_closed_curve() {
        let projected = square_at(0.0).to_2d_in_plane(&Plane3::xy()).unwrap();
        assert!(projected.is_closed());
    }

    /// ...and an open curve must not be closed by the projection.
    #[test]
    fn an_open_curve_stays_open() {
        let curve = Curve3::from_points(
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            1e-8,
        )
        .unwrap();

        assert!(!curve.to_2d_in_plane(&Plane3::xy()).unwrap().is_closed());
    }

    /// A curve built in an arbitrary plane, projected back onto it, has to reproduce the 2D curve
    /// it was lifted from. This is the round trip the section pipeline actually performs.
    #[test]
    fn a_curve_lifted_into_a_plane_projects_back_onto_itself() {
        let plane = Plane3::from_point_normal(
            &Point3::new(3.0, -1.0, 7.0),
            &UnitVec3::new_normalize(Vector3::new(0.4, -0.6, 0.7)),
        );
        let frame = plane.compute_frame();

        // Four corners of a square, expressed in the plane's own frame and lifted into the world.
        let flat = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.5, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ];
        let lifted = flat.iter().map(|p| frame * p).collect::<Vec<_>>();
        let curve = Curve3::from_points(&lifted, 1e-9).unwrap();

        let back = curve.to_2d_in_plane(&plane).unwrap();
        assert!(back.is_closed());
        assert_eq!(back.points().len(), flat.len());
        for (want, got) in flat.iter().zip(back.points().iter()) {
            assert_relative_eq!(want.to_2d(), *got, epsilon = 1e-9);
        }
    }

    /// A curve running along the plane normal collapses to a point, which is an error rather than
    /// a degenerate curve.
    #[test]
    fn a_curve_along_the_normal_collapses_and_errors() {
        let curve = Curve3::from_points(
            &[Point3::new(1.0, 2.0, 0.0), Point3::new(1.0, 2.0, 10.0)],
            1e-8,
        )
        .unwrap();

        assert!(curve.to_2d_in_plane(&Plane3::xy()).is_err());
    }

    /// Two planes differing only in offset give the same in-plane coordinates, since the frame
    /// origin moves with the plane.
    #[test]
    fn a_parallel_plane_offset_does_not_move_the_projection() {
        let curve = square_at(0.0);

        let near = curve.to_2d_in_plane(&Plane3::xy()).unwrap();
        let far = curve
            .to_2d_in_plane(&Plane3::xy().offset_by(100.0))
            .unwrap();

        for (a, b) in near.points().iter().zip(far.points().iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-12);
        }
    }
}
