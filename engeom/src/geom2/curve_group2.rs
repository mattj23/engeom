//! A collection of disjoint curves treated as a single rigid entity.

use crate::Result;
use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::curve2::CurveStation2;
use crate::geom2::{Aabb2, Curve2, Iso2};
use parry2d_f64::bounding_volume::BoundingVolume;

/// A `CurveGroup2` is a collection of disjoint `Curve2` polylines treated as a single rigid
/// entity, such as the loops and open segments produced by a planar section of a `Mesh3`.
///
/// Where a `Curve2` is a single connected polyline, a group holds any number of them, closed or
/// open in any mixture, and answers the whole-entity questions that treating them one at a time
/// would get wrong: the bounding box of everything, the closest position on *any* member, and a
/// rigid transform applied to all of them together. It plays the same role for curves that a
/// multi-component `Mesh3` plays for surfaces, where disconnected patches are still one body.
///
/// Members keep their identity: queries that land on the group report which member they landed
/// on, and per-vertex data (such as measurement uncertainty) is addressed by concatenating the
/// members' vertices in member order, with [`CurveGroup2::vertex_offset`] locating each member's
/// span.
///
/// A group is never empty; construction rejects an empty collection because a group with no
/// members has no bounding box or closest point, and every consumer would otherwise need to
/// handle a state that means nothing.
#[derive(Clone)]
pub struct CurveGroup2 {
    curves: Vec<Curve2>,

    /// The start index of each member's vertices in a concatenated per-vertex array, in member
    /// order. Computed once at construction.
    vertex_offsets: Vec<usize>,
}

impl CurveGroup2 {
    /// Creates a new `CurveGroup2` from a collection of curves.
    ///
    /// The members are kept in the order given, and that order defines both the member indices
    /// reported by queries and the concatenation order of per-vertex data.
    ///
    /// # Arguments
    ///
    /// * `curves`: the member curves. At least one is required.
    ///
    /// returns: Result<CurveGroup2, Box<dyn Error, Global>>
    pub fn new(curves: Vec<Curve2>) -> Result<Self> {
        if curves.is_empty() {
            return Err("a curve group must have at least one member curve".into());
        }

        let mut vertex_offsets = Vec::with_capacity(curves.len());
        let mut total = 0;
        for c in &curves {
            vertex_offsets.push(total);
            total += c.count();
        }

        Ok(Self {
            curves,
            vertex_offsets,
        })
    }

    /// The member curves, in member order.
    pub fn curves(&self) -> &[Curve2] {
        &self.curves
    }

    /// The number of member curves.
    pub fn len(&self) -> usize {
        self.curves.len()
    }

    /// Whether the group has no members, which construction guarantees it never does.
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty()
    }

    /// The axis-aligned bounding box enclosing every member curve.
    pub fn aabb(&self) -> Aabb2 {
        let mut aabb = self.curves[0].aabb();
        for c in self.curves.iter().skip(1) {
            aabb.merge(&c.aabb());
        }
        aabb
    }

    /// The total arc length of all member curves.
    pub fn length(&self) -> f64 {
        self.curves.iter().map(|c| c.length()).sum()
    }

    /// The total number of vertices across all member curves.
    pub fn vertex_count(&self) -> usize {
        self.vertex_offsets.last().unwrap() + self.curves.last().unwrap().count()
    }

    /// The index at which a member's vertices begin in an array of per-vertex values covering the
    /// whole group, concatenated in member order.
    pub fn vertex_offset(&self, member: usize) -> usize {
        self.vertex_offsets[member]
    }

    /// Finds the closest position on any member curve to the test point, returning the index of
    /// the owning member along with the station on it.
    ///
    /// Ties between members are broken in favor of the lower member index. The scan is linear
    /// over the members, each of which answers through its own spatial index; a section produces
    /// a handful of members, so a group-level index has nothing to earn.
    pub fn at_closest_to_point(&self, test_point: &impl PCoords<2>) -> (usize, CurveStation2<'_>) {
        self.curves
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let station = c.at_closest_to_point(test_point);
                let d = dist(test_point, &station.point());
                (i, d, station)
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _, station)| (i, station))
            .expect("a curve group is never empty")
    }

    /// Returns a new group with every member transformed by the isometry. Member order, and with
    /// it vertex addressing, is preserved.
    pub fn new_transformed_by(&self, iso: &Iso2) -> Self {
        Self {
            curves: self
                .curves
                .iter()
                .map(|c| c.new_transformed_by(iso))
                .collect(),
            vertex_offsets: self.vertex_offsets.clone(),
        }
    }
}

impl From<Curve2> for CurveGroup2 {
    /// Wraps a single curve as a group of one, which is what alignment against a lone curve
    /// uses. Infallible because the one-member group is always valid.
    fn from(curve: Curve2) -> Self {
        Self {
            vertex_offsets: vec![0],
            curves: vec![curve],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Point2, Vector2};
    use approx::assert_relative_eq;

    /// A closed unit square with its lower-left corner at (x0, y0), wound counter-clockwise. It
    /// has a length of 4 and, like every closed `Curve2`, a closing vertex duplicating the first,
    /// so it counts 5 vertices.
    fn square_at(x0: f64, y0: f64) -> Curve2 {
        let points = [
            Point2::new(x0, y0),
            Point2::new(x0 + 1.0, y0),
            Point2::new(x0 + 1.0, y0 + 1.0),
            Point2::new(x0, y0 + 1.0),
        ];
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    /// An open two-vertex segment from (0, 5) to (1, 5).
    fn segment() -> Curve2 {
        Curve2::from_points(&[Point2::new(0.0, 5.0), Point2::new(1.0, 5.0)], 1e-8, false).unwrap()
    }

    #[test]
    fn an_empty_group_is_rejected() {
        assert!(CurveGroup2::new(Vec::new()).is_err());
    }

    #[test]
    fn a_single_curve_becomes_a_group_of_one() {
        let group = CurveGroup2::from(square_at(0.0, 0.0));

        assert_eq!(group.len(), 1);
        assert!(!group.is_empty());
        assert_eq!(group.vertex_offset(0), 0);
        assert_eq!(group.vertex_count(), 5);
    }

    #[test]
    fn the_aabb_spans_every_member() {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), square_at(10.0, 0.0)]).unwrap();
        let aabb = group.aabb();

        assert_relative_eq!(aabb.mins.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.mins.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs.x, 11.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs.y, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn the_length_is_the_sum_of_the_members() {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), segment()]).unwrap();
        assert_relative_eq!(group.length(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn vertex_offsets_concatenate_the_members_in_order() {
        // A closed square carries 5 vertices (its seam vertex is stored twice), the open
        // segment 2.
        let group =
            CurveGroup2::new(vec![square_at(0.0, 0.0), segment(), square_at(3.0, 0.0)]).unwrap();

        assert_eq!(group.vertex_offset(0), 0);
        assert_eq!(group.vertex_offset(1), 5);
        assert_eq!(group.vertex_offset(2), 7);
        assert_eq!(group.vertex_count(), 12);
    }

    #[test]
    fn the_closest_point_reports_its_owning_member() {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), square_at(10.0, 0.0)]).unwrap();

        // Just left of the first square.
        let (m, station) = group.at_closest_to_point(&Point2::new(-1.0, 0.5));
        assert_eq!(m, 0);
        assert_relative_eq!(station.point(), Point2::new(0.0, 0.5), epsilon = 1e-12);

        // Just right of the second square.
        let (m, station) = group.at_closest_to_point(&Point2::new(12.0, 0.5));
        assert_eq!(m, 1);
        assert_relative_eq!(station.point(), Point2::new(11.0, 0.5), epsilon = 1e-12);

        // Between the two, nearer the first.
        let (m, _) = group.at_closest_to_point(&Point2::new(4.0, 0.5));
        assert_eq!(m, 0);
    }

    #[test]
    fn a_transformed_group_moves_every_member_together() {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), segment()]).unwrap();
        let iso = Iso2::new(Vector2::new(100.0, -50.0), std::f64::consts::FRAC_PI_2);
        let moved = group.new_transformed_by(&iso);

        assert_eq!(moved.len(), 2);
        assert_relative_eq!(moved.length(), group.length(), epsilon = 1e-10);
        assert_eq!(moved.vertex_offset(1), group.vertex_offset(1));
        assert_eq!(moved.vertex_count(), group.vertex_count());

        // Spot-check a vertex of each member against the isometry applied directly.
        for (a, b) in group.curves().iter().zip(moved.curves().iter()) {
            for (p, q) in a.points().iter().zip(b.points().iter()) {
                assert_relative_eq!(iso * p, *q, epsilon = 1e-12);
            }
        }
    }
}
