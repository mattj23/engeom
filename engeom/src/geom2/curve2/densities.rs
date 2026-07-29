//! This module has some special density-like functions for 2D curves.

use crate::common::DiscreteDomain;
use crate::geom2::LineOps2;
use crate::{Curve2, Line2, Series1, Vector2};
use parry2d_f64::partitioning::TraversalAction;

impl Curve2 {
    /// Calculates a `Series1` of the density of directional traversal along a given direction.
    ///
    /// Provide a parameterized line that defines the direction of interest. Imagine that the curve
    /// is reoriented so that the line origin is now at (0, 0) and the line direction is +Y. This
    /// method will then sweep through the X direction based on the provided spacing, looking at
    /// the bin from `x - bin_width / 2.0` to `x + bin_width / 2.0`, and calculating the amount of
    /// Y distance covered by the curve within that rectangular slice of space.
    ///
    /// The result will be a `Series1` struct containing the result. The domain (x values) of
    /// the series will correspond with the X values in the description above. The range (y values)
    /// will be the amount of Y distance covered by the curve in the associated bin.
    ///
    /// To transform an X position in the `Series1` domain back to a line in the original space,
    /// create a line at (x, 0) with a normal (0, 1) and transform it by the result of the
    /// `.to_iso_from_y()` method. Or, if you are using a `Line2` for the original line, use the
    /// `.offset_by(x)` method directly.
    ///
    /// # Arguments
    ///
    /// * `line`: A parameterized line that defines the direction of interest.
    /// * `spacing`: The spacing between bins in the X direction.
    /// * `bin_width`: The width of the bins in the X direction. The bin will range from `x -
    ///   bin_width / 2.0` to `x + bin_width / 2.0`.
    ///
    /// returns: Series1
    pub fn directional_density(
        &self,
        line: &impl LineOps2,
        spacing: f64,
        bin_width: f64,
    ) -> Series1 {
        let line = Line2::from(line);
        let iso = line.to_iso_from_y().inverse();
        let moved = self.new_transformed_by(&iso);

        let x0 = (moved.aabb().mins.x / spacing).floor() * spacing;
        let x1 = (moved.aabb().maxs.x / spacing).ceil() * spacing;
        let n = ((x1 - x0) / spacing).ceil() as usize + 1;
        let domain = DiscreteDomain::linear(x0, x1, n);
        let mut range = Vec::with_capacity(n);
        for x in domain.iter() {
            range.push(y_dist(&moved, x - bin_width / 2.0, x + bin_width / 2.0));
        }

        Series1::new(domain, range)
    }
}

fn y_dist(curve: &Curve2, x0: f64, x1: f64) -> f64 {
    let mut candidates = Vec::new();
    curve.shape.bvh().traverse(|node| {
        if let Some(data) = node.leaf_data() {
            candidates.push(data)
        };

        if node.aabb().maxs.x <= x0 || node.aabb().mins.x >= x1 {
            TraversalAction::Prune
        } else {
            TraversalAction::Continue
        }
    });

    let l0 = Line2::new([x0, 0.0].into(), Vector2::y());
    let l1 = Line2::new([x1, 0.0].into(), Vector2::y());

    let mut sum = 0.0;
    for i in candidates.iter() {
        let edge = curve.shape.indices()[*i as usize];
        let p0 = curve.shape.vertices()[edge[0] as usize];
        let p1 = curve.shape.vertices()[edge[1] as usize];

        // Check for complete miss
        if p0.x.max(p1.x) <= x0 || p0.x.min(p1.x) >= x1 {
            continue;
        }

        let l = Line2::from_points(&p0, &p1);

        let c0 = if x0 <= p0.x && p0.x <= x1 {
            p0
        } else {
            l0.intersect(&l).unwrap()
        };

        let c1 = if x0 <= p1.x && p1.x <= x1 {
            p1
        } else {
            l1.intersect(&l).unwrap()
        };

        sum += (c1.y - c0.y).abs();
    }

    sum
}

#[cfg(test)]
mod tests {
    use super::super::tests::sample_points;
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn y_dist_vertical() {
        let points = sample_points(&[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (2.0, 2.0)]);
        let curve = Curve2::from_points(&points, 0.01, false).unwrap();
        assert_abs_diff_eq!(y_dist(&curve, 0.5, 1.5), 2.0);
    }

    #[test]
    fn y_dist_diagonal() {
        let points = sample_points(&[(0.0, 0.0), (2.0, 2.0)]);
        let curve = Curve2::from_points(&points, 0.01, false).unwrap();
        assert_abs_diff_eq!(y_dist(&curve, 0.5, 1.5), 1.0);
    }

    #[test]
    fn directional_density() {
        let points = sample_points(&[(0.0, 0.0), (0.0, 1.0), (2.0, 1.0), (2.0, 2.0)]);
        let curve = Curve2::from_points(&points, 0.01, false).unwrap();
        let line = Line2::new([0.0, 0.0].into(), Vector2::x());
        let density = curve.directional_density(&line, 0.1, 0.1);

        let (x, y) = density.global_maxima_xy();
        assert_abs_diff_eq!(y, 2.0);
        assert_abs_diff_eq!(x, -1.0);
    }
}
