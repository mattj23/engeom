//! A mapping between the world and the two-dimensional coordinate system of a plane, and the
//! choices a caller has for how that coordinate system is laid out.
//!
//! A planar section of a mesh is nearly always followed by the same three steps: bring the
//! section into the x-y plane, work on it with the tools in [`crate::geom2`], and bring the 2D
//! results back into the original 3D space. [`PlanarMap`] carries both directions of that trip, so
//! the frame is chosen once and every conversion uses it consistently. [`PlaneFrame`] describes
//! the available ways to choose that frame when taking a section.
//!
//! The planar map is a thin wrapper around an isometry. Sections against other entities (a
//! cylinder, for instance) will need a mapping that is not an isometry and whose vector conversions
//! depend on where the vector sits. The per-entity API here is designed so such mappings can offer
//! the same operations later.

use crate::common::PCoords;
use crate::geom2::{
    Curve2, CurveGroup2, Line2, Point2, Segment2, SurfacePoint2, UnitVec2, Vector2,
};
use crate::geom3::{
    Curve3, CurveGroup3, Iso3, IsoExtensions3, Line3, Plane3, Point3, Segment3, SurfacePoint3,
    SvdBasis3, UnitVec3, Vector3,
};
use crate::{Result, To2D, To3D};

/// The smallest fraction of a direction vector's length that may survive projection into the
/// plane before the direction is treated as parallel to the normal.
const MIN_PROJECTED_FRACTION: f64 = 1e-9;

/// How to lay out the two-dimensional coordinate system of a planar section.
///
/// A plane on its own fixes only the z axis of its frame. Where the origin sits and which way x
/// points within the plane remain free choices. These choices matter whenever the caller needs to
/// inspect the 2D result or relate it to something else, so the section accepts this option and
/// returns the resulting [`PlanarMap`].
#[derive(Debug, Clone, Default)]
pub enum PlaneFrame {
    /// The plane's own deterministic frame, [`Plane3::compute_frame`]. The origin is the point
    /// of the plane nearest the world origin and x is the world x axis projected into the plane.
    /// Use this when the layout does not matter; it is stable for a given plane and costs
    /// nothing to compute.
    #[default]
    Auto,

    /// The principal axes of the section itself: the origin at its centroid and x along its
    /// direction of greatest extent, from a length-weighted SVD of the section's vertices. The
    /// weighting makes the result independent of the curves' tessellation density. The x
    /// direction is ill-conditioned for a section with no dominant direction, such as a circle,
    /// where it is still deterministic but not stable against small changes in the data.
    Svd,

    /// An origin and an x direction supplied by the caller. The origin is projected onto the
    /// plane and the direction into it, so neither needs to lie there already, but the direction
    /// must not be parallel to the plane normal.
    Oriented { origin: Point3, x: Vector3 },
}

/// Maps geometry between the world and the two-dimensional coordinate system of a plane.
///
/// The map stores the plane's frame—the isometry that takes a point expressed in the plane's
/// coordinates to its position in the world—together with the inverse frame. Mapping to 2D applies
/// the inverse and then drops z; mapping to 3D adds a zero z and then applies the frame. Storing both
/// directions lets a map be reused without repeatedly computing the inverse.
///
/// Going to 3D never fails. Going to 2D loses the component along the plane normal, so any entity
/// with a direction can collapse. Examples include a unit vector or normal parallel to the plane
/// normal, a segment running along it, or a curve whose vertices fall together. Those conversions
/// return an error rather than a degenerate result, matching [`crate::common::To2D`] for curves.
///
/// Every conversion here takes and returns `engeom` types directly. A point off the plane is
/// projected onto it when mapped to 2D; use [`inverse`](Self::inverse) directly if the distance
/// from the plane is needed as well.
#[derive(Debug, Clone)]
pub struct PlanarMap {
    frame: Iso3,
    inverse: Iso3,
}

impl PlanarMap {
    fn from_frame(frame: Iso3) -> Self {
        Self {
            inverse: frame.inverse(),
            frame,
        }
    }

    /// Creates the map for a plane's own deterministic frame, [`Plane3::compute_frame`]. This
    /// is the [`PlaneFrame::Auto`] choice.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane whose frame the map uses
    ///
    /// returns: PlanarMap
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Plane3, PlanarMap, Point2, Point3};
    ///
    /// let map = PlanarMap::from_plane(&Plane3::xy().offset_by(5.0));
    /// let p2 = map.point_to_2d(&Point3::new(2.0, 3.0, 5.0));
    /// assert_relative_eq!(p2, Point2::new(2.0, 3.0), epsilon = 1e-12);
    /// assert_relative_eq!(map.point_to_3d(&p2), Point3::new(2.0, 3.0, 5.0), epsilon = 1e-12);
    /// ```
    pub fn from_plane(plane: &Plane3) -> Self {
        Self::from_frame(plane.compute_frame())
    }

    /// Creates a map whose origin and x direction are chosen by the caller. This is the
    /// [`PlaneFrame::Oriented`] choice.
    ///
    /// The origin is projected onto the plane and `x` is projected into it, so neither needs to
    /// lie there already. The frame's z axis is the plane normal and y completes a right-handed
    /// set.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane the map lies in
    /// * `origin`: the point that becomes `(0, 0)` in the plane after being projected onto it
    /// * `x`: the direction that becomes the plane's x axis after being projected into it
    ///
    /// returns: Result<PlanarMap, Box<dyn Error, Global>>
    ///
    /// # Failure
    ///
    /// Fails if `x` is zero or parallel to the plane normal, since nothing is left of it in the
    /// plane to point the x axis along.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Plane3, PlanarMap, Point2, Point3, Vector3};
    ///
    /// // The origin is off the plane and the direction is tilted out of it; both are projected.
    /// let map = PlanarMap::from_plane_oriented(
    ///     &Plane3::xy(),
    ///     &Point3::new(1.0, 1.0, 7.0),
    ///     &Vector3::new(0.0, 1.0, 3.0),
    /// ).unwrap();
    ///
    /// assert_relative_eq!(map.point_to_2d(&Point3::new(1.0, 1.0, 0.0)), Point2::origin(), epsilon = 1e-12);
    /// assert_relative_eq!(map.point_to_2d(&Point3::new(1.0, 3.0, 0.0)), Point2::new(2.0, 0.0), epsilon = 1e-12);
    /// ```
    pub fn from_plane_oriented(plane: &Plane3, origin: &Point3, x: &Vector3) -> Result<Self> {
        let u = plane.project_vector(x);
        if u.norm() < MIN_PROJECTED_FRACTION * x.norm() || x.norm() == 0.0 {
            return Err("the x direction is parallel to the plane normal".into());
        }
        let v = plane.normal.cross(&u);
        let origin = plane.project_point(origin);

        Ok(Self::from_frame(Iso3::from_basis_xy(&u, &v, Some(origin))?))
    }

    /// Creates a map from the principal axes of a set of points lying in (or near) the plane.
    /// This is the [`PlaneFrame::Svd`] choice.
    ///
    /// The origin is the (weighted) centroid of the points projected onto the plane, and x is
    /// the direction of greatest extent projected into it, with its sign chosen so that it points
    /// the same way as the x axis of [`from_plane`](Self::from_plane) would. This keeps the result
    /// deterministic for a given plane and set of points. The z axis is the plane normal rather
    /// than the smallest SVD axis, which is only approximately the normal.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane the map lies in
    /// * `points`: the points whose principal axes decide the frame
    /// * `weights`: optional per-point weights, one per point. A section assigns each vertex half
    ///   the length of every segment it touches so that the frame does not depend on tessellation
    ///   density.
    ///
    /// returns: Result<PlanarMap, Box<dyn Error, Global>>
    ///
    /// # Failure
    ///
    /// Fails if there are no points, if the weights do not match the points in number, or if the
    /// direction of greatest extent is parallel to the plane normal, meaning the points were not
    /// in the plane at all.
    ///
    /// The x direction is ill-conditioned when the two largest singular values are close, as they
    /// are for a circle. The result is deterministic, but a small change in the data can alter the
    /// direction substantially. That is not an error, since the frame is still valid, but a caller
    /// who needs a stable direction in that case should use
    /// [`from_plane_oriented`](Self::from_plane_oriented).
    pub fn from_plane_svd(
        plane: &Plane3,
        points: &[impl PCoords<3>],
        weights: Option<&[f64]>,
    ) -> Result<Self> {
        if let Some(w) = weights
            && w.len() != points.len()
        {
            return Err(
                format!("{} weights were given for {} points", w.len(), points.len()).into(),
            );
        }

        let basis = SvdBasis3::from_points(points, weights)
            .ok_or("the principal axes could not be computed from the points")?;

        let mut x = basis.basis[0];
        if x.dot(&plane.compute_frame().x()) < 0.0 {
            x = -x;
        }

        Self::from_plane_oriented(plane, &basis.center, &x).map_err(|e| {
            format!("the direction of greatest extent is not in the plane: {e}").into()
        })
    }

    /// The frame isometry, which takes a point expressed in the plane's coordinates (with z zero
    /// for a point on the plane) to where it lies in the world.
    pub fn frame(&self) -> &Iso3 {
        &self.frame
    }

    /// The inverse frame, which takes a world point into the plane's coordinates. A point on the
    /// plane has a zero z coordinate, while the z coordinate of any other point is its signed
    /// distance from the plane.
    pub fn inverse(&self) -> &Iso3 {
        &self.inverse
    }

    /// The plane this map lies in, whose normal is the frame's z axis.
    pub fn plane(&self) -> Plane3 {
        Plane3::from_point_normal(&self.frame.origin(), &self.frame.z())
    }

    // --------------------------------------------------------------------------------------------
    // Points and vectors
    // --------------------------------------------------------------------------------------------

    /// Maps a point into the plane's coordinates, dropping its distance from the plane.
    pub fn point_to_2d(&self, point: &impl PCoords<3>) -> Point2 {
        (self.inverse * Point3::from(point.coords())).to_2d()
    }

    /// Maps a point in the plane's coordinates to where it lies on the plane in the world.
    pub fn point_to_3d(&self, point: &impl PCoords<2>) -> Point3 {
        self.frame * Point2::from(point.coords()).to_3d()
    }

    /// Maps each point into the plane's coordinates. See [`point_to_2d`](Self::point_to_2d).
    pub fn points_to_2d(&self, points: &[impl PCoords<3>]) -> Vec<Point2> {
        points.iter().map(|p| self.point_to_2d(p)).collect()
    }

    /// Maps each point onto the plane in the world. See [`point_to_3d`](Self::point_to_3d).
    pub fn points_to_3d(&self, points: &[impl PCoords<2>]) -> Vec<Point3> {
        points.iter().map(|p| self.point_to_3d(p)).collect()
    }

    /// Maps a vector into the plane's coordinates, dropping its component along the normal.
    pub fn vector_to_2d(&self, vector: &Vector3) -> Vector2 {
        (self.inverse * vector).to_2d()
    }

    /// Maps a vector in the plane's coordinates to the world.
    pub fn vector_to_3d(&self, vector: &Vector2) -> Vector3 {
        self.frame * vector.to_3d()
    }

    /// Maps a unit vector into the plane's coordinates, dropping its component along the normal
    /// and renormalizing what remains.
    ///
    /// # Failure
    ///
    /// Fails if the vector is parallel to the plane normal, since nothing of it remains in the
    /// plane.
    pub fn unit_vector_to_2d(&self, vector: &UnitVec3) -> Result<UnitVec2> {
        UnitVec2::try_new((self.inverse * vector).into_inner().to_2d(), 1e-12)
            .ok_or_else(|| "the unit vector is parallel to the plane normal".into())
    }

    /// Maps a unit vector in the plane's coordinates to the world.
    pub fn unit_vector_to_3d(&self, vector: &UnitVec2) -> UnitVec3 {
        self.frame * vector.to_3d()
    }

    // --------------------------------------------------------------------------------------------
    // Surface points, lines and segments
    // --------------------------------------------------------------------------------------------

    /// Maps a surface point into the plane's coordinates. The point is projected onto the plane
    /// and the normal into it.
    ///
    /// # Failure
    ///
    /// Fails if the normal is parallel to the plane normal.
    pub fn surface_point_to_2d(&self, sp: &SurfacePoint3) -> Result<SurfacePoint2> {
        let normal = self
            .unit_vector_to_2d(&sp.normal)
            .map_err(|_| "the surface point's normal is parallel to the plane normal")?;
        Ok(SurfacePoint2::new(self.point_to_2d(&sp.point), normal))
    }

    /// Maps a surface point in the plane's coordinates to the world.
    pub fn surface_point_to_3d(&self, sp: &SurfacePoint2) -> SurfacePoint3 {
        SurfacePoint3::new(
            self.point_to_3d(&sp.point),
            self.unit_vector_to_3d(&sp.normal),
        )
    }

    /// Maps a line into the plane's coordinates. The origin is projected onto the plane and the
    /// direction into it. The direction's length is whatever survives the projection, so a line
    /// with a unit direction is not unit-speed afterwards unless it lay in the plane.
    ///
    /// # Failure
    ///
    /// Fails if the direction is parallel to the plane normal.
    pub fn line_to_2d(&self, line: &Line3) -> Result<Line2> {
        let direction = self.vector_to_2d(&line.direction);
        if direction.norm() < MIN_PROJECTED_FRACTION * line.direction.norm()
            || line.direction.norm() == 0.0
        {
            return Err("the line direction is parallel to the plane normal".into());
        }
        Ok(Line2::new(self.point_to_2d(&line.origin), direction))
    }

    /// Maps a line in the plane's coordinates to the world.
    pub fn line_to_3d(&self, line: &Line2) -> Line3 {
        Line3::new(
            self.point_to_3d(&line.origin),
            self.vector_to_3d(&line.direction),
        )
    }

    /// Maps a segment into the plane's coordinates by projecting both endpoints.
    ///
    /// # Failure
    ///
    /// Fails if the endpoints fall together, which happens when the segment runs along the plane
    /// normal.
    pub fn segment_to_2d(&self, segment: &Segment3) -> Result<Segment2> {
        Segment2::new(&self.point_to_2d(&segment.a), &self.point_to_2d(&segment.b))
            .map_err(|_| "the segment runs along the plane normal".into())
    }

    /// Maps a segment in the plane's coordinates to the world.
    pub fn segment_to_3d(&self, segment: &Segment2) -> Segment3 {
        Segment3::new_unchecked(self.point_to_3d(&segment.a), self.point_to_3d(&segment.b))
    }

    // --------------------------------------------------------------------------------------------
    // Curves and curve groups
    // --------------------------------------------------------------------------------------------

    /// Projects a curve onto the plane and returns it as a two-dimensional curve in the plane's
    /// coordinates.
    ///
    /// For a map built on `Plane3::xy()` this agrees vertex for vertex with
    /// [`crate::common::To2D::to_2d`], which only ever drops the z coordinate, because that
    /// plane's frame is the identity.
    ///
    /// The curve's chordal tolerance carries over, and the result is closed when the projected
    /// first and last vertices are within that tolerance of each other. A section of a mesh
    /// therefore comes back closed if it was a loop, since such a curve repeats its first vertex
    /// as its last.
    ///
    /// # Failure
    ///
    /// The projection can bring vertices closer together than the tolerance, and the
    /// de-duplication that follows may leave fewer than two distinct points. A curve running
    /// along the plane normal collapses to a point this way and is an error rather than a
    /// degenerate curve.
    pub fn curve_to_2d(&self, curve: &Curve3) -> Result<Curve2> {
        let points = self.points_to_2d(curve.vertices());
        Curve2::from_points(&points, curve.tol(), false)
    }

    /// Lifts a two-dimensional curve in the plane's coordinates onto the plane in the world.
    ///
    /// The tolerance carries over. A closed `Curve2` stores a repeated first vertex as its last,
    /// so the result repeats it too. This lets [`curve_to_2d`](Self::curve_to_2d) recover the
    /// closure on the way back.
    ///
    /// # Failure
    ///
    /// Fails only if the curve cannot be rebuilt from its vertices, which a valid `Curve2` does
    /// not cause. The result type is retained so callers do not have to unwrap this conversion.
    pub fn curve_to_3d(&self, curve: &Curve2) -> Result<Curve3> {
        let points = self.points_to_3d(curve.points());
        Curve3::from_points(&points, curve.tol())
    }

    /// Projects every member of a group onto the plane and returns the result as a
    /// two-dimensional group in the plane's coordinates. Member order is preserved.
    ///
    /// This is the standard way to bring a planar section of a mesh into two dimensions for the
    /// tools in [`crate::geom2`], including the multi-curve alignment, which works on
    /// [`CurveGroup2`] bodies.
    ///
    /// # Failure
    ///
    /// A member that collapses under the projection is an error for the whole group rather than
    /// being dropped from it. Silently returning a smaller group would renumber the members, and
    /// a member that collapses indicates that the group was not in this plane to begin with. The
    /// error is therefore reported instead of hidden.
    pub fn curve_group_to_2d(&self, group: &CurveGroup3) -> Result<CurveGroup2> {
        let mut curves = Vec::with_capacity(group.len());
        for (i, curve) in group.curves().iter().enumerate() {
            let projected = self.curve_to_2d(curve).map_err(|e| {
                format!("member curve {i} could not be projected onto the plane: {e}")
            })?;
            curves.push(projected);
        }

        CurveGroup2::new(curves)
    }

    /// Lifts every member of a two-dimensional group onto the plane in the world. Member order
    /// is preserved.
    ///
    /// # Failure
    ///
    /// Fails if any member cannot be rebuilt from its vertices; see
    /// [`curve_to_3d`](Self::curve_to_3d).
    pub fn curve_group_to_3d(&self, group: &CurveGroup2) -> Result<CurveGroup3> {
        let mut curves = Vec::with_capacity(group.len());
        for (i, curve) in group.curves().iter().enumerate() {
            let lifted = self
                .curve_to_3d(curve)
                .map_err(|e| format!("member curve {i} could not be lifted to the plane: {e}"))?;
            curves.push(lifted);
        }

        CurveGroup3::new(curves)
    }
}

/// Per-vertex weights that make a weighted fit independent of how densely a curve is tessellated.
/// Each vertex receives half the length of every segment it touches, so every segment's length is
/// split between its endpoints and the weights sum to the curve's length. A loop that repeats its
/// first vertex as its last splits that corner's weight across the two copies, keeping the weighted
/// centroid at the loop's centroid.
pub(crate) fn vertex_length_weights(curve: &Curve3) -> Vec<f64> {
    let v = curve.vertices();
    let mut weights = vec![0.0; v.len()];
    for i in 1..v.len() {
        let half = (v[i] - v[i - 1]).norm() / 2.0;
        weights[i - 1] += half;
        weights[i] += half;
    }
    weights
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn tilted_plane() -> Plane3 {
        Plane3::from_point_normal(
            &Point3::new(1.0, -2.0, 3.0),
            &UnitVec3::new_normalize(Vector3::new(1.0, 2.0, 3.0)),
        )
    }

    /// An axis-aligned rectangle `w` by `h` centered on `c`, in the plane's coordinates, lifted
    /// into the world by the map. Its first vertex is repeated as its last, as a section's is.
    fn rectangle_in(map: &PlanarMap, c: Point2, w: f64, h: f64) -> Curve3 {
        let corners = [
            Point2::new(c.x - w / 2.0, c.y - h / 2.0),
            Point2::new(c.x + w / 2.0, c.y - h / 2.0),
            Point2::new(c.x + w / 2.0, c.y + h / 2.0),
            Point2::new(c.x - w / 2.0, c.y + h / 2.0),
            Point2::new(c.x - w / 2.0, c.y - h / 2.0),
        ];
        Curve3::from_points(&map.points_to_3d(&corners), 1e-9).unwrap()
    }

    // --- Auto -------------------------------------------------------------------------------

    #[test]
    fn xy_plane_map_is_the_identity() {
        let map = PlanarMap::from_plane(&Plane3::xy());
        let p = Point3::new(1.5, -2.5, 0.0);
        assert_relative_eq!(map.point_to_2d(&p), Point2::new(1.5, -2.5), epsilon = 1e-15);
        assert_relative_eq!(map.point_to_3d(&Point2::new(1.5, -2.5)), p, epsilon = 1e-15);
        assert_relative_eq!(
            map.vector_to_3d(&Vector2::new(3.0, 4.0)),
            Vector3::new(3.0, 4.0, 0.0),
            epsilon = 1e-15
        );
    }

    /// The x-y plane's frame is the identity, so this projection must agree exactly with `to_2d`.
    /// If they ever diverge, one of the two is wrong.
    #[test]
    fn projecting_onto_the_xy_plane_is_the_same_as_dropping_z() {
        let map = PlanarMap::from_plane(&Plane3::xy());
        let curve = rectangle_in(&map, Point2::new(0.5, 0.5), 1.0, 1.0);
        let raised = Curve3::from_points(
            &curve
                .vertices()
                .iter()
                .map(|p| Point3::new(p.x, p.y, 3.0))
                .collect::<Vec<_>>(),
            1e-9,
        )
        .unwrap();

        let dropped = raised.to_2d().unwrap();
        let projected = map.curve_to_2d(&raised).unwrap();

        assert_eq!(dropped.points().len(), projected.points().len());
        for (a, b) in dropped.points().iter().zip(projected.points().iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
    }

    /// Two planes differing only in offset give the same in-plane coordinates, since the frame
    /// origin moves with the plane.
    #[test]
    fn a_parallel_plane_offset_does_not_move_the_projection() {
        let near = PlanarMap::from_plane(&Plane3::xy());
        let far = PlanarMap::from_plane(&Plane3::xy().offset_by(100.0));
        let curve = rectangle_in(&near, Point2::new(0.5, 0.5), 1.0, 1.0);

        let a = near.curve_to_2d(&curve).unwrap();
        let b = far.curve_to_2d(&curve).unwrap();
        for (p, q) in a.points().iter().zip(b.points().iter()) {
            assert_relative_eq!(p, q, epsilon = 1e-12);
        }
    }

    #[test]
    fn from_plane_agrees_with_compute_frame() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let frame = plane.compute_frame();
        assert_relative_eq!(map.frame(), &frame, epsilon = 1e-15);
        assert_relative_eq!(map.inverse(), &frame.inverse(), epsilon = 1e-15);
    }

    #[test]
    fn plane_accessor_recovers_the_plane() {
        let plane = tilted_plane();
        let recovered = PlanarMap::from_plane(&plane).plane();
        assert_relative_eq!(recovered.normal, plane.normal, epsilon = 1e-12);
        assert_relative_eq!(recovered.d, plane.d, epsilon = 1e-12);
    }

    #[test]
    fn point_to_2d_projects_onto_the_plane() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let off = Point3::new(4.0, 5.0, 6.0);
        let back = map.point_to_3d(&map.point_to_2d(&off));
        assert_relative_eq!(back, plane.project_point(&off), epsilon = 1e-12);
        assert_relative_eq!(
            (map.inverse() * off).z,
            plane.signed_distance_to_point(&off),
            epsilon = 1e-12
        );
    }

    // --- Oriented ---------------------------------------------------------------------------

    #[test]
    fn oriented_places_origin_and_x_axis() {
        let plane = tilted_plane();
        let origin = Point3::new(0.5, 0.5, 0.5);
        let x = Vector3::new(0.0, 0.0, 1.0);
        let map = PlanarMap::from_plane_oriented(&plane, &origin, &x).unwrap();

        assert_relative_eq!(
            map.frame().origin(),
            plane.project_point(&origin),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            map.frame().x().into_inner(),
            plane.project_vector(&x).normalize(),
            epsilon = 1e-12
        );
        assert_relative_eq!(map.frame().z(), plane.normal, epsilon = 1e-12);
        // The frame is right-handed, with y in the plane.
        assert_relative_eq!(
            map.frame().x().cross(&map.frame().y()),
            map.frame().z().into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn oriented_rejects_x_along_the_normal() {
        let plane = tilted_plane();
        let x = plane.normal.into_inner() * 2.0;
        assert!(PlanarMap::from_plane_oriented(&plane, &Point3::origin(), &x).is_err());
        assert!(
            PlanarMap::from_plane_oriented(&plane, &Point3::origin(), &Vector3::zeros()).is_err()
        );
    }

    // --- Svd --------------------------------------------------------------------------------

    #[test]
    fn svd_puts_x_along_the_long_side_and_origin_at_the_center() {
        let plane = tilted_plane();
        let auto = PlanarMap::from_plane(&plane);
        let center = Point2::new(3.0, -1.0);

        // A rectangle whose long side runs along the auto frame's y axis.
        let rect = rectangle_in(&auto, center, 2.0, 10.0);
        let weights = vertex_length_weights(&rect);
        let map = PlanarMap::from_plane_svd(&plane, rect.vertices(), Some(&weights)).unwrap();

        assert_relative_eq!(
            map.frame().origin(),
            auto.point_to_3d(&center),
            epsilon = 1e-9
        );
        assert_relative_eq!(
            map.frame().x().dot(&auto.frame().y()).abs(),
            1.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(map.frame().z(), plane.normal, epsilon = 1e-12);

        // The rectangle is now axis-aligned with its long side along x.
        let flat = map.curve_to_2d(&rect).unwrap();
        let aabb = flat.aabb();
        assert_relative_eq!(aabb.maxs.x - aabb.mins.x, 10.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.y - aabb.mins.y, 2.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.center(), Point2::origin(), epsilon = 1e-9);
    }

    #[test]
    fn svd_sign_follows_the_auto_frame() {
        let plane = tilted_plane();
        let auto = PlanarMap::from_plane(&plane);
        // Long side along the auto x axis: the SVD x must point the same way, not opposite.
        let rect = rectangle_in(&auto, Point2::origin(), 10.0, 2.0);
        let weights = vertex_length_weights(&rect);
        let map = PlanarMap::from_plane_svd(&plane, rect.vertices(), Some(&weights)).unwrap();
        assert!(map.frame().x().dot(&auto.frame().x()) > 0.99);
    }

    #[test]
    fn svd_rejects_mismatched_weights_and_no_points() {
        let plane = Plane3::xy();
        let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        assert!(PlanarMap::from_plane_svd(&plane, &points, Some(&[1.0])).is_err());
        let none: [Point3; 0] = [];
        assert!(PlanarMap::from_plane_svd(&plane, &none, None).is_err());
    }

    #[test]
    fn svd_rejects_points_along_the_normal() {
        let plane = Plane3::xy();
        let points = [Point3::new(0.0, 0.0, -1.0), Point3::new(0.0, 0.0, 1.0)];
        assert!(PlanarMap::from_plane_svd(&plane, &points, None).is_err());
    }

    #[test]
    fn length_weights_split_each_segment_between_its_ends() {
        let map = PlanarMap::from_plane(&Plane3::xy());
        let rect = rectangle_in(&map, Point2::origin(), 2.0, 10.0);
        let w = vertex_length_weights(&rect);
        assert_eq!(w.len(), 5);
        assert_relative_eq!(w.iter().sum::<f64>(), rect.length(), epsilon = 1e-12);
        // The repeated closing vertex splits its corner's weight across the two copies.
        assert_relative_eq!(w[0] + w[4], 6.0, epsilon = 1e-12);
        for &wi in &w[1..4] {
            assert_relative_eq!(wi, 6.0, epsilon = 1e-12);
        }
    }

    /// Length weighting prevents refining one side of a loop from dragging the frame toward it.
    #[test]
    fn svd_frame_is_independent_of_tessellation_density() {
        let plane = tilted_plane();
        let auto = PlanarMap::from_plane(&plane);
        let coarse = rectangle_in(&auto, Point2::new(2.0, 3.0), 6.0, 3.0);

        // The same rectangle with one long side subdivided into many vertices.
        let mut pts = Vec::new();
        pts.push(Point2::new(-1.0, 1.5));
        pts.push(Point2::new(-1.0, 4.5));
        for i in 0..=40 {
            pts.push(Point2::new(-1.0 + 6.0 * i as f64 / 40.0, 4.5));
        }
        pts.push(Point2::new(5.0, 1.5));
        pts.push(Point2::new(-1.0, 1.5));
        let fine = Curve3::from_points(&auto.points_to_3d(&pts), 1e-9).unwrap();

        let wc = vertex_length_weights(&coarse);
        let wf = vertex_length_weights(&fine);
        let mc = PlanarMap::from_plane_svd(&plane, coarse.vertices(), Some(&wc)).unwrap();
        let mf = PlanarMap::from_plane_svd(&plane, fine.vertices(), Some(&wf)).unwrap();

        assert_relative_eq!(mc.frame().origin(), mf.frame().origin(), epsilon = 1e-9);
        assert_relative_eq!(mc.frame().x(), mf.frame().x(), epsilon = 1e-9);
    }

    // --- Vectors, surface points, lines, segments -------------------------------------------

    #[test]
    fn in_plane_vectors_round_trip() {
        let map = PlanarMap::from_plane(&tilted_plane());
        let v2 = Vector2::new(0.3, -0.7);
        assert_relative_eq!(
            map.vector_to_2d(&map.vector_to_3d(&v2)),
            v2,
            epsilon = 1e-12
        );

        let u2 = UnitVec2::new_normalize(Vector2::new(2.0, 1.0));
        let u3 = map.unit_vector_to_3d(&u2);
        assert_relative_eq!(u3.norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(map.unit_vector_to_2d(&u3).unwrap(), u2, epsilon = 1e-12);
    }

    #[test]
    fn vector_to_2d_drops_the_normal_component() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let v = Vector3::new(1.0, 2.0, 3.0);
        let back = map.vector_to_3d(&map.vector_to_2d(&v));
        assert_relative_eq!(back, plane.project_vector(&v), epsilon = 1e-12);
        assert!(map.unit_vector_to_2d(&plane.normal).is_err());
    }

    #[test]
    fn surface_points_round_trip_and_reject_a_normal_along_the_plane_normal() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let sp2 = SurfacePoint2::new_normalize(Point2::new(1.0, 2.0), Vector2::new(-1.0, 1.0));
        let sp3 = map.surface_point_to_3d(&sp2);
        assert_relative_eq!(plane.distance_to_point(&sp3.point), 0.0, epsilon = 1e-12);
        assert_relative_eq!(sp3.normal.dot(&plane.normal), 0.0, epsilon = 1e-12);

        let back = map.surface_point_to_2d(&sp3).unwrap();
        assert_relative_eq!(back.point, sp2.point, epsilon = 1e-12);
        assert_relative_eq!(back.normal, sp2.normal, epsilon = 1e-12);

        let bad = SurfacePoint3::new(sp3.point, plane.normal);
        assert!(map.surface_point_to_2d(&bad).is_err());
    }

    #[test]
    fn lines_round_trip_and_reject_a_direction_along_the_normal() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let l2 = Line2::new(Point2::new(1.0, 1.0), Vector2::new(2.0, -1.0));
        let l3 = map.line_to_3d(&l2);
        assert_relative_eq!(plane.distance_to_point(&l3.origin), 0.0, epsilon = 1e-12);
        assert_relative_eq!(l3.direction.dot(&plane.normal), 0.0, epsilon = 1e-12);

        let back = map.line_to_2d(&l3).unwrap();
        assert_relative_eq!(back.origin, l2.origin, epsilon = 1e-12);
        assert_relative_eq!(back.direction, l2.direction, epsilon = 1e-12);

        let bad = Line3::new(l3.origin, plane.normal.into_inner());
        assert!(map.line_to_2d(&bad).is_err());
    }

    #[test]
    fn segments_round_trip_and_reject_one_along_the_normal() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let s2 = Segment2::new(&Point2::new(0.0, 0.0), &Point2::new(3.0, 4.0)).unwrap();
        let s3 = map.segment_to_3d(&s2);
        assert_relative_eq!(s3.length(), 5.0, epsilon = 1e-12);

        let back = map.segment_to_2d(&s3).unwrap();
        assert_relative_eq!(back.a, s2.a, epsilon = 1e-12);
        assert_relative_eq!(back.b, s2.b, epsilon = 1e-12);

        let bad = Segment3::new_unchecked(s3.a, s3.a + plane.normal.into_inner());
        assert!(map.segment_to_2d(&bad).is_err());
    }

    // --- Curves and groups ------------------------------------------------------------------

    #[test]
    fn curve_round_trip_keeps_closure_and_length() {
        let map = PlanarMap::from_plane(&tilted_plane());
        let rect = rectangle_in(&map, Point2::new(1.0, 1.0), 4.0, 2.0);

        let flat = map.curve_to_2d(&rect).unwrap();
        assert!(flat.is_closed());
        assert_relative_eq!(flat.length(), rect.length(), epsilon = 1e-9);

        let back = map.curve_to_3d(&flat).unwrap();
        assert_eq!(back.count(), rect.count());
        for (a, b) in back.vertices().iter().zip(rect.vertices()) {
            assert_relative_eq!(a, b, epsilon = 1e-9);
        }
        assert!(map.curve_to_2d(&back).unwrap().is_closed());
    }

    #[test]
    fn open_curve_stays_open() {
        let map = PlanarMap::from_plane(&tilted_plane());
        let pts = map.points_to_3d(&[
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
        ]);
        let curve = Curve3::from_points(&pts, 1e-9).unwrap();
        assert!(!map.curve_to_2d(&curve).unwrap().is_closed());
    }

    #[test]
    fn curve_along_the_normal_is_an_error() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let n = plane.normal.into_inner();
        let base = Point3::new(1.0, 2.0, 3.0);
        let curve = Curve3::from_points(&[base, base + n, base + n * 2.0], 1e-6).unwrap();
        assert!(map.curve_to_2d(&curve).is_err());
    }

    #[test]
    fn group_round_trip_preserves_member_order_and_names_a_failing_member() {
        let plane = tilted_plane();
        let map = PlanarMap::from_plane(&plane);
        let big = rectangle_in(&map, Point2::origin(), 10.0, 10.0);
        let small = rectangle_in(&map, Point2::origin(), 2.0, 2.0);
        let group = CurveGroup3::new(vec![big.clone(), small.clone()]).unwrap();

        let flat = map.curve_group_to_2d(&group).unwrap();
        assert_eq!(flat.len(), 2);
        assert_relative_eq!(flat.curves()[0].length(), big.length(), epsilon = 1e-9);
        assert_relative_eq!(flat.curves()[1].length(), small.length(), epsilon = 1e-9);

        let back = map.curve_group_to_3d(&flat).unwrap();
        assert_eq!(back.len(), 2);
        assert_relative_eq!(back.length(), group.length(), epsilon = 1e-9);

        let n = plane.normal.into_inner();
        let base = Point3::new(1.0, 2.0, 3.0);
        let bad = Curve3::from_points(&[base, base + n], 1e-6).unwrap();
        let group = CurveGroup3::new(vec![big, bad]).unwrap();
        let err = match map.curve_group_to_2d(&group) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("a member along the normal should fail the group"),
        };
        assert!(err.contains("member curve 1"), "{err}");
    }
}
