//! A collection of disjoint 3D curves treated as a single rigid entity.

use crate::Result;
use crate::common::points::dist;
use crate::geom3::curve3::CurveStation3;
use crate::geom3::{Aabb3, Curve3, Iso3, Point3};
use parry3d_f64::bounding_volume::BoundingVolume;

/// A `CurveGroup3` is a collection of disjoint `Curve3` polylines treated as a single rigid
/// entity, such as the loops and open segments produced by a planar section of a `Mesh3` before
/// they are brought into two dimensions.
///
/// Where a `Curve3` is a single connected polyline, a group holds any number of them, and answers
/// the whole-entity questions that treating them one at a time would get wrong: the bounding box
/// of everything, the closest position on *any* member, and a rigid transform applied to all of
/// them together. It plays the same role for curves that a multi-component `Mesh3` plays for
/// surfaces, where disconnected patches are still one body.
///
/// Members keep their identity: queries that land on the group report which member they landed
/// on.
///
/// A group is never empty; construction rejects an empty collection because a group with no
/// members has no bounding box or closest point, and every consumer would otherwise need to
/// handle a state that means nothing.
///
/// This is the 3D counterpart of [`crate::geom2::CurveGroup2`] and is deliberately structured
/// identically to it. Unlike the 2D type, a `Curve3` carries no notion of being closed, so a loop
/// is represented by its first and last vertices coinciding.
#[derive(Clone)]
pub struct CurveGroup3 {
    curves: Vec<Curve3>,
}

impl CurveGroup3 {
    /// Creates a new `CurveGroup3` from a collection of curves.
    ///
    /// The members are kept in the order given, and that order defines both the member indices
    /// reported by queries and the concatenation order of per-vertex data.
    ///
    /// # Arguments
    ///
    /// * `curves`: the member curves. At least one is required.
    ///
    /// returns: Result<CurveGroup3, Box<dyn Error, Global>>
    pub fn new(curves: Vec<Curve3>) -> Result<Self> {
        if curves.is_empty() {
            return Err("a curve group must have at least one member curve".into());
        }

        Ok(Self { curves })
    }

    /// The member curves, in member order.
    pub fn curves(&self) -> &[Curve3] {
        &self.curves
    }

    /// Consumes the group and returns its member curves, in member order.
    pub fn into_curves(self) -> Vec<Curve3> {
        self.curves
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
    pub fn aabb(&self) -> Aabb3 {
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

    /// Finds the closest position on any member curve to the test point, returning the index of
    /// the owning member along with the station on it.
    ///
    /// Ties between members are broken in favor of the lower member index. The scan is linear
    /// over the members, each of which answers through its own spatial index; a section produces
    /// a handful of members, so a group-level index has nothing to earn.
    pub fn at_closest_to_point(&self, test_point: &Point3) -> (usize, CurveStation3<'_>) {
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

    /// Returns a new group with every member transformed by the isometry. Member order is
    /// preserved.
    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            curves: self
                .curves
                .iter()
                .map(|c| c.new_transformed_by(iso))
                .collect(),
        }
    }
}

impl From<Curve3> for CurveGroup3 {
    /// Wraps a single curve as a group of one. Infallible because the one-member group is always
    /// valid.
    fn from(curve: Curve3) -> Self {
        Self {
            curves: vec![curve],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Vector3;
    use approx::assert_relative_eq;

    /// A closed unit square in the plane `z = z0`, with its lower-left corner at (x0, y0), wound
    /// counter-clockwise, with a length of 4. The first vertex is repeated as the last, which is
    /// how a `Curve3` represents a loop.
    fn square_at(x0: f64, y0: f64, z0: f64) -> Curve3 {
        let points = [
            Point3::new(x0, y0, z0),
            Point3::new(x0 + 1.0, y0, z0),
            Point3::new(x0 + 1.0, y0 + 1.0, z0),
            Point3::new(x0, y0 + 1.0, z0),
            Point3::new(x0, y0, z0),
        ];
        Curve3::from_points(&points, 1e-8).unwrap()
    }

    /// An open two-vertex segment from (0, 5, 0) to (1, 5, 0).
    fn segment() -> Curve3 {
        Curve3::from_points(
            &[Point3::new(0.0, 5.0, 0.0), Point3::new(1.0, 5.0, 0.0)],
            1e-8,
        )
        .unwrap()
    }

    #[test]
    fn an_empty_group_is_rejected() {
        assert!(CurveGroup3::new(Vec::new()).is_err());
    }

    #[test]
    fn a_single_curve_becomes_a_group_of_one() {
        let group = CurveGroup3::from(square_at(0.0, 0.0, 0.0));

        assert_eq!(group.len(), 1);
        assert!(!group.is_empty());
    }

    #[test]
    fn the_aabb_spans_every_member() {
        let group =
            CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), square_at(10.0, 0.0, 3.0)]).unwrap();
        let aabb = group.aabb();

        assert_relative_eq!(aabb.mins.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.mins.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.mins.z, 0.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs.x, 11.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs.y, 1.0, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs.z, 3.0, epsilon = 1e-12);
    }

    #[test]
    fn the_length_is_the_sum_of_the_members() {
        let group = CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), segment()]).unwrap();
        assert_relative_eq!(group.length(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn the_closest_point_reports_its_owning_member() {
        let group =
            CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), square_at(10.0, 0.0, 0.0)]).unwrap();

        // Just left of the first square.
        let (m, station) = group.at_closest_to_point(&Point3::new(-1.0, 0.5, 0.0));
        assert_eq!(m, 0);
        assert_relative_eq!(station.point(), Point3::new(0.0, 0.5, 0.0), epsilon = 1e-12);

        // Just right of the second square.
        let (m, station) = group.at_closest_to_point(&Point3::new(12.0, 0.5, 0.0));
        assert_eq!(m, 1);
        assert_relative_eq!(
            station.point(),
            Point3::new(11.0, 0.5, 0.0),
            epsilon = 1e-12
        );

        // Between the two, nearer the first.
        let (m, _) = group.at_closest_to_point(&Point3::new(4.0, 0.5, 0.0));
        assert_eq!(m, 0);
    }

    #[test]
    fn a_transformed_group_moves_every_member_together() {
        let group = CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), segment()]).unwrap();
        let iso = Iso3::new(
            Vector3::new(100.0, -50.0, 7.0),
            Vector3::z() * std::f64::consts::FRAC_PI_2,
        );
        let moved = group.new_transformed_by(&iso);

        assert_eq!(moved.len(), 2);
        assert_relative_eq!(moved.length(), group.length(), epsilon = 1e-10);

        // Spot-check every vertex of each member against the isometry applied directly.
        for (a, b) in group.curves().iter().zip(moved.curves().iter()) {
            for (p, q) in a.vertices().iter().zip(b.vertices().iter()) {
                assert_relative_eq!(iso * p, *q, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn a_group_keeps_its_members_in_order() {
        let group = CurveGroup3::new(vec![
            square_at(0.0, 0.0, 0.0),
            segment(),
            square_at(10.0, 0.0, 0.0),
        ])
        .unwrap();

        let lengths: Vec<f64> = group.curves().iter().map(|c| c.length()).collect();
        assert_relative_eq!(lengths[0], 4.0, epsilon = 1e-12);
        assert_relative_eq!(lengths[1], 1.0, epsilon = 1e-12);
        assert_relative_eq!(lengths[2], 4.0, epsilon = 1e-12);

        let back = group.clone().into_curves();
        assert_eq!(back.len(), 3);
        assert_relative_eq!(back[1].length(), 1.0, epsilon = 1e-12);
    }
}
