use crate::common::points::{dist, max_point_in_direction};
use crate::common::{Intersection, PCoords};
use crate::geom2::polyline2::polyline_intersections;
use crate::geom2::{Aabb2, LineOps2, intersection_param};
use crate::{Curve2, CurveStation2, Point2, SurfacePoint2, Vector2};
use parry2d_f64::partitioning::TraversalAction;
use parry2d_f64::query::{PointQueryWithLocation, Ray};
use crate::geom2::line2::slab_method2;

impl Curve2 {
    /// Returns the `CurveStation2` whose point is closest to `test_point`.
    ///
    /// # Arguments
    ///
    /// * `test_point`: The reference point to project onto the curve.
    ///
    /// returns: CurveStation2
    pub fn at_closest_to_point(&self, test_point: &impl PCoords<2>) -> CurveStation2<'_> {
        let p = Point2::from(test_point.coords());
        let (prj, loc) = self.shape.project_local_point_and_get_location(&p, false);
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

    /// Returns the minimum distance from any point on the curve to `test_point`.
    ///
    /// # Arguments
    ///
    /// * `test_point`: The point to measure the distance to.
    ///
    /// returns: f64
    pub fn dist_to_point(&self, test_point: &Point2) -> f64 {
        let (prj, _) = self
            .shape
            .project_local_point_and_get_location(test_point, false);
        dist(&prj.point, test_point)
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
        max_point_in_direction(self.points(), vector)
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
    #[deprecated]
    pub fn ray_intersections(&self, ray: &Ray) -> Vec<(f64, usize)> {
        polyline_intersections(&self.shape, ray)
    }

    /// Performs an intersection between the curve and any entity implementing `LineOps2`. Returns
    /// a sorted, deduplicated list of intersections, each as a pair of values representing the
    /// distance along the line and the index of the edge where the intersection occurs.
    ///
    /// # Arguments
    ///
    /// * `line`: The line to intersect with the curve.
    ///
    /// returns: Vec<(f64, usize), Global>
    pub fn intersections_with_line(&self, line: &impl LineOps2) -> Vec<(f64, usize)> {
        let mut candidates = Vec::new();
        let r_inv = Vector2::new(1.0 / line.dir().x, 1.0 / line.dir().y);

        self.shape.bvh().traverse(|node| {
            if let Some(data) = node.leaf_data() {
                candidates.push(data)
            };

            if !slab_method2(&node.aabb(), &line.origin(), &r_inv) {
                TraversalAction::Prune
            } else {
                TraversalAction::Continue
            }
        });

        let mut results = Vec::new();
        for i in candidates.iter() {
            if let Some(t) = intersect_with_edge(self, line, *i as usize) {
                results.push((t, *i as usize));
            }
        }

        results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        results.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-8 && a.1 == b.1);

        results
    }
}

impl<T: LineOps2> Intersection<T, Vec<(f64, usize)>> for Curve2 {
    fn intersection(&self, other: T) -> Vec<(f64, usize)> {
        self.intersections_with_line(&other)
    }
}

fn intersect_with_edge(curve: &Curve2, line: &impl LineOps2, edge_index: usize) -> Option<f64> {
    let v0 = curve.points()[edge_index];
    let v1 = curve.points()[edge_index + 1];
    let dir = v1 - v0;
    if let Some((t0, t1)) = intersection_param(&line.origin(), &line.dir(), &v0, &dir) {
        if (0.0..=1.0).contains(&t1) {
            Some(t0)
        } else {
            None
        }
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::*;
    use super::*;
    use crate::Line2;
    use approx::assert_relative_eq;

    #[test]
    fn intersection_simple() {
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, true).unwrap();

        let line = Line2::new([0.5, 0.0].into(), [0.0, 1.0].into());
        let intersections = curve
            .intersections_with_line(&line)
            .iter()
            .map(|i| i.0)
            .collect::<Vec<_>>();
        assert_eq!(intersections.len(), 2);

        assert_relative_eq!(intersections[0], 0.0, epsilon = 1e-6);
        assert_relative_eq!(intersections[1], 1.0, epsilon = 1e-6);
    }

    #[test]
    fn intersection_edge() {
        // Verifies that the edge of the curve still intersects the line. If this fails, it may
        // indicate an issue with the slab method of the AABB when taking an intersection with
        // a line coincident with the bounding volume's edge.
        let curve = Curve2::from_points(&sample_points(&sample1()), 1e-6, false).unwrap();

        let line = Line2::new([0.0, 0.0].into(), [0.0, 1.0].into());
        let intersections = curve
            .intersections_with_line(&line)
            .iter()
            .map(|i| i.0)
            .collect::<Vec<_>>();
        assert_eq!(intersections.len(), 2);

        assert_relative_eq!(intersections[0], 0.0, epsilon = 1e-6);
        assert_relative_eq!(intersections[1], 1.0, epsilon = 1e-6);
    }
}
