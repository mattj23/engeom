//! A collection of disjoint 3D curves treated as a single rigid entity.

use crate::Result;
use crate::common::points::dist;
use crate::geom2::CurveGroup2;
use crate::geom3::curve3::CurveStation3;
use crate::geom3::{Aabb3, Curve3, Iso3, Plane3, Point3};
use crate::io::{read_tc_curves3_file, write_tc_curves3_file};
use parry3d_f64::bounding_volume::BoundingVolume;
use std::path::Path;

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

    /// Read a group from a tolerance-compressed `.tccurve3` file, taking every curve in it as a
    /// member in the order the file stores them.
    ///
    /// This is the counterpart of [`CurveGroup3::save_tccurve3`], and it also reads a file written
    /// by [`Curve3::save_tccurve3`], which arrives as a group of one.
    ///
    /// # Failure
    ///
    /// A file holding no curves at all is refused, since a group is never empty.
    pub fn load_tccurve3(path: &Path) -> Result<Self> {
        Self::new(read_tc_curves3_file(path)?)
    }

    /// Write this group to a single tolerance-compressed `.tccurve3` file, one item per member.
    ///
    /// Member order is the file item order, so a group read back has the same member indices it
    /// was saved with. Each member keeps its own chord tolerance.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to write to, which is overwritten if it already exists
    /// * `tol`: the largest acceptable round-trip position error for any vertex of any member, in
    ///   the same units as the coordinates
    ///
    /// returns: `Result<()>`
    pub fn save_tccurve3(&self, path: &Path, tol: f64) -> Result<()> {
        write_tc_curves3_file(path, &self.curves, tol)
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

    /// Projects every member onto a plane and returns the result as a two-dimensional group,
    /// expressed in that plane's own coordinate frame. Member order is preserved.
    ///
    /// This is the ordinary way to bring a planar section of a mesh into two dimensions for the
    /// tools in [`crate::geom2`], including the multi-curve alignment, which works on
    /// [`CurveGroup2`] bodies. Passing the same plane that produced the section is the usual case,
    /// and it recovers a faithful 2D copy of it.
    ///
    /// See [`Curve3::to_2d_in_plane`] for how each member is converted, and
    /// [`Plane3::compute_frame`] for the frame the result is expressed in.
    ///
    /// # Failure
    ///
    /// A member which collapses under the projection is an error for the whole group rather than
    /// being dropped from it. Silently returning a smaller group would renumber the members, and a
    /// member that collapses is a sign the group was not in this plane to begin with, so it is
    /// better raised than hidden.
    pub fn to_2d_in_plane(&self, plane: &Plane3) -> Result<CurveGroup2> {
        let mut curves = Vec::with_capacity(self.curves.len());
        for (i, curve) in self.curves.iter().enumerate() {
            let projected = curve.to_2d_in_plane(plane).map_err(|e| {
                format!("member curve {i} could not be projected onto the plane: {e}")
            })?;
            curves.push(projected);
        }

        CurveGroup2::new(curves)
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

    /// A path in the temp directory, tagged so concurrent tests cannot collide on it.
    fn temp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("engeom_curve_group3_{tag}.tccurve3"))
    }

    /// Order, count and each member's own chord tolerance have to survive the file, the same way
    /// they do in 2D.
    #[test]
    fn a_group_round_trips_through_a_tccurve3_file() -> Result<()> {
        let group = CurveGroup3::new(vec![
            square_at(0.0, 0.0, 0.0),
            segment(),
            square_at(10.0, 0.0, 3.0),
        ])?;
        let path = temp_path("round_trip");
        let tol = 1e-6;

        group.save_tccurve3(&path, tol)?;
        let back = CurveGroup3::load_tccurve3(&path)?;

        assert_eq!(back.len(), group.len());
        for (i, (a, b)) in group.curves().iter().zip(back.curves().iter()).enumerate() {
            assert_eq!(a.count(), b.count(), "member {i} vertex count");
            assert_relative_eq!(a.tol(), b.tol(), epsilon = 1e-15);
            for (p, q) in a.vertices().iter().zip(b.vertices().iter()) {
                assert_relative_eq!(p, q, epsilon = tol);
            }
        }

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// A loop is a loop because its first and last vertices coincide, and that is the only record
    /// of its closure a `Curve3` has. If the file did not preserve it, projecting a restored
    /// section would silently produce open curves.
    #[test]
    fn a_restored_loop_still_projects_to_a_closed_curve() -> Result<()> {
        let group = CurveGroup3::from(square_at(0.0, 0.0, 0.0));
        let path = temp_path("loop");

        group.save_tccurve3(&path, 1e-9)?;
        let back = CurveGroup3::load_tccurve3(&path)?;

        let flat = back.to_2d_in_plane(&Plane3::xy())?;
        assert!(flat.curves()[0].is_closed());

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// A single-curve file is a collection of one, and loads as a group of one.
    #[test]
    fn a_single_curve_file_loads_as_a_group_of_one() -> Result<()> {
        let path = temp_path("single");
        square_at(0.0, 0.0, 0.0).save_tccurve3(&path, 1e-6)?;

        assert_eq!(CurveGroup3::load_tccurve3(&path)?.len(), 1);

        let _ = std::fs::remove_file(&path);
        Ok(())
    }

    /// ...and a group file is refused by the single-curve reader rather than truncated to its
    /// first member.
    #[test]
    fn a_group_file_is_refused_by_the_single_curve_reader() -> Result<()> {
        let group = CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), segment()])?;
        let path = temp_path("refuse_single");
        group.save_tccurve3(&path, 1e-6)?;

        assert!(Curve3::load_tccurve3(&path).is_err());

        let _ = std::fs::remove_file(&path);
        Ok(())
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

    // ============================================================================================
    // Projection into a plane
    // ============================================================================================

    /// Every member has to come through, in order, and the group's own invariants with it.
    #[test]
    fn a_projected_group_keeps_its_members_and_their_order() {
        let group = CurveGroup3::new(vec![
            square_at(0.0, 0.0, 0.0),
            segment(),
            square_at(10.0, 0.0, 0.0),
        ])
        .unwrap();

        let flat = group.to_2d_in_plane(&Plane3::xy()).unwrap();

        assert_eq!(flat.len(), 3);
        assert_relative_eq!(flat.curves()[0].length(), 4.0, epsilon = 1e-12);
        assert_relative_eq!(flat.curves()[1].length(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(flat.curves()[2].length(), 4.0, epsilon = 1e-12);

        // The closed members stay closed and the open one stays open.
        assert!(flat.curves()[0].is_closed());
        assert!(!flat.curves()[1].is_closed());
        assert!(flat.curves()[2].is_closed());
    }

    /// A member that collapses fails the whole group rather than quietly renumbering the rest.
    #[test]
    fn a_collapsing_member_fails_the_group() {
        let along_normal = Curve3::from_points(
            &[Point3::new(5.0, 5.0, 0.0), Point3::new(5.0, 5.0, 9.0)],
            1e-8,
        )
        .unwrap();
        let group = CurveGroup3::new(vec![square_at(0.0, 0.0, 0.0), along_normal]).unwrap();

        let err = match group.to_2d_in_plane(&Plane3::xy()) {
            Err(e) => e.to_string(),
            Ok(g) => panic!("a collapsing member yielded a group of {}", g.len()),
        };
        assert!(err.contains("member curve 1"), "unexpected message: {err}");
    }

    /// The end-to-end path this exists for: cut a mesh, bring the section into the plane it came
    /// from, and get back a faithful 2D copy of it.
    #[test]
    fn a_mesh_section_projects_back_into_its_own_plane() {
        use crate::Mesh3;

        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let plane = Plane3::xy();

        let section = mesh.section_with_plane(&plane, Some(1e-9), None).unwrap();
        let flat = section.to_2d_in_plane(&plane).unwrap();

        assert_eq!(flat.len(), 1);

        let loop_2d = &flat.curves()[0];
        assert!(loop_2d.is_closed(), "a box sections into a closed loop");

        // The box is 2 by 4 in the cutting plane, so the loop is its perimeter.
        assert_relative_eq!(loop_2d.length(), 12.0, epsilon = 1e-6);

        let aabb = loop_2d.aabb();
        assert_relative_eq!(aabb.mins.x, -1.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.x, 1.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.mins.y, -2.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.y, 2.0, epsilon = 1e-9);
    }

    /// The same cut taken on a tilted plane has to produce the same 2D shape, since the frame
    /// follows the plane. This is what makes the projection independent of how the part is posed.
    #[test]
    fn a_tilted_section_has_the_same_shape_as_an_upright_one() {
        use crate::{Iso3, Mesh3};

        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let upright = mesh
            .section_with_plane(&Plane3::xy(), Some(1e-9), None)
            .unwrap()
            .to_2d_in_plane(&Plane3::xy())
            .unwrap();

        // Rotate the part, and cut with the correspondingly rotated plane.
        let iso = Iso3::new(Vector3::new(10.0, -5.0, 2.0), Vector3::new(0.3, -0.5, 0.2));
        let moved = mesh.transform_copy(&iso);
        let plane = Plane3::xy().transformed_by(&iso);

        let tilted = moved
            .section_with_plane(&plane, Some(1e-9), None)
            .unwrap()
            .to_2d_in_plane(&plane)
            .unwrap();

        assert_eq!(tilted.len(), upright.len());
        assert_relative_eq!(tilted.length(), upright.length(), epsilon = 1e-6);
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
