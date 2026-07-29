//! Construction methods for `BoundaryData2`

use crate::common::PCoords;
use crate::geom2::BoundaryData2;
use crate::geom2::boundary2::data::{BData, BoundaryAddData};
use crate::{AngleDir, Circle2, Line2, Point2, Result};

// ===============================================================================================
//  Common boundary construction tools
// ===============================================================================================

pub trait BoundaryEditor: BoundaryAddData {
    fn add_seg_xy(&mut self, x: f64, y: f64) -> u32 {
        let seg = BData::Seg((x, y));
        self.add_data(seg)
    }

    fn add_seg(&mut self, p: &impl PCoords<2>) -> u32 {
        self.add_seg_xy(p.coords().x, p.coords().y)
    }

    fn add_arc_xy(
        &mut self,
        center_x: f64,
        center_y: f64,
        end_x: f64,
        end_y: f64,
        clockwise: bool,
    ) -> u32 {
        let arc = BData::Arc((center_x, center_y, end_x, end_y, clockwise));
        self.add_data(arc)
    }

    fn add_arc(&mut self, center: &impl PCoords<2>, end: &impl PCoords<2>, clockwise: bool) -> u32 {
        self.add_arc_xy(
            center.coords().x,
            center.coords().y,
            end.coords().x,
            end.coords().y,
            clockwise,
        )
    }

    /// This method will create a sequence of line segments which meet at corner fillets to the
    /// boundary at the current cursor position.
    ///
    /// The geometry will be defined through a sequence of points and a common fillet radius. For
    /// `n` points, there will be `n` line segments and `n-1` fillets created. You must provide at
    /// least two points in the `corners` argument.
    ///
    /// # Arguments
    ///
    /// * `corners`: the endpoints of the line segments
    /// * `radius`: the radius of the fillets joining the line segments
    ///
    /// returns: Result<Vec<u32, Global>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn add_corner_fillets(&mut self, corners: &[impl PCoords<2>], radius: f64) -> Result<Vec<u32>> {
        if corners.len() < 2 {
            return Err("Must have at least 2 corners to fillet".into());
        }
        let mut ids = Vec::new();

        for i in 0..corners.len() - 1 {
            let c0 = &corners[i];
            let c1 = &corners[i + 1];

            let last_point = self
                .last_point()
                .ok_or("Cannot add corner fillet to empty boundary")?;

            // First we'll find the circle at the tangent of the corner
            let v0 = last_point.coords() - c0.coords();
            let v1 = c1.coords() - c0.coords();
            let c = Circle2::from_tangent_to_corner(c0, &v0, &v1, radius)?;

            // Now we need to find the two tangent endpoints
            let l0 = Line2::from_points(&last_point, c0).normalized();
            let (e0, e1) = c.tangent_points_to(c0).unwrap();
            let (e0, e1) = if l0.distance_to(&e0) < l0.distance_to(&e1) {
                (e0, e1)
            } else {
                (e1, e0)
            };

            // Check the length and figure out the direction
            if l0.scalar_project(&e0) < 0.0 {
                return Err("Invalid fillet radius, too large for corner".into());
            }
            let clockwise = l0.signed_distance_to(&e1).is_sign_positive();

            ids.push(self.add_seg(&e0));
            ids.push(self.add_arc(&c.center, &e1, clockwise));
        }

        let cf1 = &corners[corners.len() - 1];
        let cf0 = &corners[corners.len() - 2];
        let lf = Line2::from_points(cf0, cf1).normalized();
        let end = self
            .last_point()
            .ok_or("Cannot add corner to empty boundary")?;
        if lf.scalar_project(&end) < 0.0 {
            return Err("Invalid fillet radius, too large for corner".into());
        }
        ids.push(self.add_seg(cf1));

        Ok(ids)
    }

    /// Create a full round joining two line segments.  The first line segment starts at the last
    /// working point and goes to the tangency of a circle located at `center` of size `radius`.
    /// Then it adds an arc that following the perimeter of the circle until its tangency points at
    /// `end`. Then the final arc segment goes from the tangent point at the end of the arc to
    /// the `end` point.
    ///
    /// This method will pick the arc direction (clockwise/counter-clockwise) that results in
    /// the shorter arc length. This will prevent the side which intersects from being used, but
    /// if you need more control than this method offers, you will likely need to generate the
    /// geometry yourself.
    ///
    /// # Arguments
    ///
    /// * `center`: the center of the arc
    /// * `radius`: the radius of the arc
    /// * `end`: the final end point
    ///
    /// returns: Result<(u32, u32, u32), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn add_full_round(
        &mut self,
        center: &impl PCoords<2>,
        radius: f64,
        end: &impl PCoords<2>,
        direction: AngleDir,
    ) -> Result<(u32, u32, u32)> {
        let c = Circle2::new(center.coords().x, center.coords().y, radius);
        let start = self
            .last_point()
            .ok_or("Cannot add full round to empty boundary")?;

        let Some((start_cw, start_ccw)) = c.tangent_points_to(&start) else {
            return Err("Invalid full round, start point is inside the circle radius".into());
        };

        let Some((end_ccw, end_cw)) = c.tangent_points_to(end) else {
            return Err("Invalid full round, end point is inside the circle radius".into());
        };

        // Get the clockwise arc
        let (i0, i1) = match direction {
            AngleDir::Ccw => (
                self.add_seg(&start_ccw),
                self.add_arc(&c.center, &end_ccw, false),
            ),
            AngleDir::Cw => (
                self.add_seg(&start_cw),
                self.add_arc(&c.center, &end_cw, true),
            ),
        };

        let i2 = self.add_seg(end);
        Ok((i0, i1, i2))
    }
}

// ===============================================================================================
//  Boundary cursor
// ===============================================================================================

pub struct BCursor<'a> {
    data: &'a mut BoundaryData2,
    node_id: u32,
}

impl<'a> BCursor<'a> {
    pub fn new(data: &'a mut BoundaryData2, node_id: u32) -> Self {
        Self { data, node_id }
    }

    pub fn move_to(&mut self, id: u32) -> Result<()> {
        if self.data.get_node(id).is_none() {
            Err("Invalid cursor position".into())
        } else {
            self.node_id = id;
            Ok(())
        }
    }
}

impl BoundaryAddData for BCursor<'_> {
    fn add_data(&mut self, data: BData) -> u32 {
        let next_id = if self.data.is_empty() {
            self.data.insert_first(data).unwrap()
        } else {
            self.data.insert_after(self.node_id, data).unwrap()
        };

        self.node_id = next_id;
        next_id
    }

    fn last_point(&self) -> Option<Point2> {
        if self.data.is_empty() {
            if self.data.is_closed() {
                None
            } else {
                self.data.start
            }
        } else {
            self.data
                .get_node(self.node_id)
                .map(|node| node.data.end_point())
        }
    }
}

impl BoundaryEditor for BCursor<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AngleDir::{Ccw, Cw};
    use crate::common::points::to_points;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn corner_fillet_single() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let mut cursor = data.get_cursor(None);
        // cursor.add_corner_fillets(&[Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)], 0.25)?;
        cursor.add_corner_fillets(&to_points(&[[1.0, 0.0], [1.0, 1.0]]), 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 * 2.0 + 0.25 * PI / 2.0;

        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }

    #[test]
    fn corner_fillet_single_inverted() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let mut cursor = data.get_cursor(None);
        cursor.add_corner_fillets(&to_points(&[[1.0, 0.0], [1.0, -1.0]]), 0.25)?;

        // data.add_corner_fillets(&[Point2::new(1.0, 0.0), Point2::new(1.0, -1.0)], 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 * 2.0 + 0.25 * PI / 2.0;

        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }

    #[test]
    fn corner_fillet_double() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let mut cursor = data.get_cursor(None);
        let points = to_points(&[[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        cursor.add_corner_fillets(&points, 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 + 0.5 + 0.75 + 2.0 * 0.25 * PI / 2.0;
        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }

    #[test]
    fn create_simple_rectangle() -> Result<()> {
        let mut data = BoundaryData2::new_closed();
        let mut cursor = data.get_cursor(None);
        cursor.add_seg_xy(1.0, 0.0);
        cursor.add_seg_xy(1.0, 1.0);
        cursor.add_seg_xy(0.0, 1.0);
        cursor.add_seg_xy(0.0, 0.0);

        let boundary = data.try_to_boundary()?;
        assert_relative_eq!(boundary.length(), 4.0, epsilon = 1e-12);

        let aabb = boundary.aabb();
        assert_relative_eq!(aabb.mins, [0.0, 0.0].into(), epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs, [1.0, 1.0].into(), epsilon = 1e-12);

        assert_relative_eq!(
            boundary.at_start().point,
            [0.0, 0.0].into(),
            epsilon = 1e-12
        );
        assert_relative_eq!(boundary.at_end().point, [0.0, 0.0].into(), epsilon = 1e-12);

        Ok(())
    }

    #[test]
    fn full_round_ccw() -> Result<()> {
        let mut data = BoundaryData2::new_open_xy(0.0, 0.0);
        data.add_full_round(&Point2::new(1.0, 1.0), 1.0, &Point2::new(0.0, 2.0), Ccw)?;
        let boundary = data.try_to_boundary()?;
        assert_relative_eq!(boundary.length(), 2.0 + PI);
        Ok(())
    }

    #[test]
    fn full_round_cw() -> Result<()> {
        let mut data = BoundaryData2::new_open_xy(0.0, 2.0);
        data.add_full_round(&Point2::new(1.0, 1.0), 1.0, &Point2::new(0.0, 0.0), Cw)?;
        let boundary = data.try_to_boundary()?;
        assert_relative_eq!(boundary.length(), 2.0 + PI);
        Ok(())
    }
}
