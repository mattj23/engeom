//! A collection of disjoint curves treated as a single rigid entity.

use crate::Result;
use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::curve2::CurveStation2;
use crate::geom2::{Aabb2, Curve2, Iso2};
use crate::io::{read_tc_curves2_file, write_tc_curves2_file};
use parry2d_f64::bounding_volume::BoundingVolume;
use std::path::Path;

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
/// on.
///
/// A group is never empty; construction rejects an empty collection because a group with no
/// members has no bounding box or closest point, and every consumer would otherwise need to
/// handle a state that means nothing.
#[derive(Clone)]
pub struct CurveGroup2 {
    curves: Vec<Curve2>,
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

        Ok(Self { curves })
    }

    /// The member curves, in member order.
    pub fn curves(&self) -> &[Curve2] {
        &self.curves
    }

    /// Consumes the group and returns its member curves, in member order.
    pub fn into_curves(self) -> Vec<Curve2> {
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

    /// Read a group from a tolerance-compressed `.tccurve2` file, taking every curve in it as a
    /// member in the order the file stores them.
    ///
    /// This is the counterpart of [`CurveGroup2::save_tccurve2`], and it also reads a file written
    /// by [`Curve2::save_tccurve2`], which arrives as a group of one.
    ///
    /// # Failure
    ///
    /// A file holding no curves at all is refused, since a group is never empty.
    pub fn load_tccurve2(path: &Path) -> Result<Self> {
        Self::new(read_tc_curves2_file(path)?)
    }

    /// Write this group to a single tolerance-compressed `.tccurve2` file, one item per member.
    ///
    /// Member order is the file item order, so a group read back has the same member indices it
    /// was saved with. Each member keeps its own chord tolerance and closed state, so a group
    /// mixing closed loops with open strands round-trips as one.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `tol`: the largest acceptable round-trip position error for any vertex of any member, in
    ///   the same units as the coordinates
    ///
    /// returns: `Result<()>`
    pub fn save_tccurve2(&self, path: &Path, tol: f64) -> Result<()> {
        write_tc_curves2_file(path, &self.curves, tol)
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

    /// Returns a new group in which open members whose ends meet have been joined into single
    /// curves. This reassembles a section through separately meshed but touching patches into the
    /// loops and strands it describes.
    ///
    /// This is [`Curve2::chain_merge`] applied to the members: repeatedly, the pair of open
    /// members with the smallest distance from the end of one to the start of the other is
    /// joined end-to-start, until no pair is within `max_dist`. Only end-to-start joins are made,
    /// so every member keeps its direction and no member is reversed. Closed members take no
    /// part and pass through unchanged. A chain whose last vertex comes within the first member's
    /// tolerance of its first vertex becomes a closed curve; a chain that stops
    /// short of that stays open, and [`Curve2::closed_within`] is the way to close it across a
    /// larger gap.
    ///
    /// The result has one member for each remaining chain, so a section through a part with a
    /// hole reduces to two loops, not one. A result containing a single curve has `len() == 1`.
    ///
    /// # Arguments
    ///
    /// * `max_dist`: the largest end-to-start distance that may be bridged by a join, or `None`
    ///   to keep joining the closest pair until nothing open is left to join
    ///
    /// returns: `Result<CurveGroup2>`
    pub fn chain_merged(&self, max_dist: Option<f64>) -> Result<Self> {
        Self::new(Curve2::chain_merge(self.curves.clone(), max_dist)?)
    }

    /// Returns a new group with every member transformed by the isometry. Member order is
    /// preserved.
    pub fn new_transformed_by(&self, iso: &Iso2) -> Self {
        Self {
            curves: self
                .curves
                .iter()
                .map(|c| c.new_transformed_by(iso))
                .collect(),
        }
    }
}

impl From<Curve2> for CurveGroup2 {
    /// Wraps a single curve as a group of one, which is what alignment against a lone curve
    /// uses. Infallible because the one-member group is always valid.
    fn from(curve: Curve2) -> Self {
        Self {
            curves: vec![curve],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Point2, Vector2};
    use approx::assert_relative_eq;

    /// A closed unit square with its lower-left corner at (x0, y0), wound counter-clockwise,
    /// with a length of 4.
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

    /// A path in the temp directory, tagged so concurrent tests cannot collide on it.
    fn temp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("engeom_curve_group2_{tag}.tccurve2"))
    }

    /// The four sides of a unit square, listed out of order as separate open members. Each side is
    /// shortened by `gap` at its end so that consecutive sides do not quite touch.
    fn square_sides_with_gap(gap: f64) -> Vec<Curve2> {
        let corners = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        let mut sides = (0..4)
            .map(|i| {
                let a = corners[i];
                let b = corners[(i + 1) % 4];
                let end = a + (b - a) * (1.0 - gap);
                Curve2::from_points(&[a, end], 1e-8, false).unwrap()
            })
            .collect::<Vec<_>>();
        sides.swap(0, 2);
        sides.swap(1, 3);
        sides
    }

    /// Sides that meet exactly form one chain. Because the chain returns to its own start, the
    /// result is closed without any further step.
    #[test]
    fn touching_sides_merge_into_one_closed_loop() -> Result<()> {
        let group = CurveGroup2::new(square_sides_with_gap(0.0))?;
        let merged = group.chain_merged(Some(1e-6))?;

        assert_eq!(merged.len(), 1);
        let loop_ = &merged.curves()[0];
        assert!(loop_.is_closed());
        assert_relative_eq!(loop_.area().unwrap(), 1.0, epsilon = 1e-12);
        Ok(())
    }

    /// Sides that stop short of each other form one open curve when the gap is within the limit.
    /// Calling `closed_within` on that curve then turns it into a loop.
    #[test]
    fn gapped_sides_merge_into_one_open_chain() -> Result<()> {
        let group = CurveGroup2::new(square_sides_with_gap(0.05))?;
        let merged = group.chain_merged(Some(0.1))?;

        assert_eq!(merged.len(), 1);
        let chain = &merged.curves()[0];
        assert!(!chain.is_closed());
        assert_eq!(chain.count(), 8);

        let closed = chain.closed_within(0.1)?;
        assert!(closed.is_closed());
        assert_relative_eq!(closed.area().unwrap(), 1.0, epsilon = 0.1);
        Ok(())
    }

    #[test]
    fn a_gap_beyond_the_limit_is_not_bridged() -> Result<()> {
        let group = CurveGroup2::new(square_sides_with_gap(0.05))?;
        let merged = group.chain_merged(Some(0.01))?;
        assert_eq!(merged.len(), 4);
        Ok(())
    }

    /// Closed members take no part in merging and come through unchanged alongside the chains.
    #[test]
    fn closed_members_pass_through_a_merge() -> Result<()> {
        let mut members = square_sides_with_gap(0.0);
        members.push(square_at(10.0, 0.0));
        let merged = CurveGroup2::new(members)?.chain_merged(None)?;

        assert_eq!(merged.len(), 2);
        assert!(merged.curves().iter().all(|c| c.is_closed()));
        assert_eq!(merged.length(), 8.0);
        Ok(())
    }

    /// The whole point of putting a group in one file: order, count, closedness and each member's
    /// own chord tolerance all have to survive, or the group that comes back is a different body
    /// from the one that went in.
    #[test]
    fn a_group_round_trips_through_a_tccurve2_file() -> Result<()> {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), segment(), square_at(10.0, 0.0)])?;
        let path = temp_path("round_trip");
        let tol = 1e-6;

        group.save_tccurve2(&path, tol)?;
        let back = CurveGroup2::load_tccurve2(&path)?;

        assert_eq!(back.len(), group.len());
        for (i, (a, b)) in group.curves().iter().zip(back.curves().iter()).enumerate() {
            assert_eq!(a.is_closed(), b.is_closed(), "member {i} closed state");
            assert_eq!(a.points().len(), b.points().len(), "member {i} point count");
            assert_relative_eq!(a.tol(), b.tol(), epsilon = 1e-15);
            for (p, q) in a.points().iter().zip(b.points().iter()) {
                assert_relative_eq!(p, q, epsilon = tol);
            }
        }

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// Members are told apart by their index, so a group whose members are the same shape in a
    /// different order must not come back reordered.
    #[test]
    fn member_order_is_the_file_order() -> Result<()> {
        let a = CurveGroup2::new(vec![segment(), square_at(0.0, 0.0)])?;
        let path = temp_path("order");

        a.save_tccurve2(&path, 1e-6)?;
        let back = CurveGroup2::load_tccurve2(&path)?;

        // The segment is open and one unit long; the square is closed and four.
        assert!(!back.curves()[0].is_closed());
        assert_relative_eq!(back.curves()[0].length(), 1.0, epsilon = 1e-9);
        assert!(back.curves()[1].is_closed());
        assert_relative_eq!(back.curves()[1].length(), 4.0, epsilon = 1e-9);

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// A single-curve file is just a collection of one, so loading it as a group has to work
    /// rather than being a different format the caller has to know about in advance.
    #[test]
    fn a_single_curve_file_loads_as_a_group_of_one() -> Result<()> {
        let curve = square_at(0.0, 0.0);
        let path = temp_path("single");

        curve.save_tccurve2(&path, 1e-6)?;
        let group = CurveGroup2::load_tccurve2(&path)?;

        assert_eq!(group.len(), 1);
        assert!(group.curves()[0].is_closed());

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// The reverse is not true, and must not silently become true: reading one curve out of a
    /// multi-curve file would discard the rest with nothing to show for it.
    #[test]
    fn a_group_file_is_refused_by_the_single_curve_reader() -> Result<()> {
        let group = CurveGroup2::new(vec![square_at(0.0, 0.0), segment()])?;
        let path = temp_path("refuse_single");
        group.save_tccurve2(&path, 1e-6)?;

        // Matched rather than unwrapped: `Curve2` is not `Debug`, so `unwrap_err` will not compile.
        let refused = Curve2::load_tccurve2(&path).is_err();
        assert!(refused, "a two-curve file must not read back as one curve");

        let _ = std::fs::remove_file(&path);
        Ok(())
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

        // Spot-check a vertex of each member against the isometry applied directly.
        for (a, b) in group.curves().iter().zip(moved.curves().iter()) {
            for (p, q) in a.points().iter().zip(b.points().iter()) {
                assert_relative_eq!(iso * p, *q, epsilon = 1e-12);
            }
        }
    }
}
