//! Construction methods for `BoundaryData2`

use crate::common::PCoords;
use crate::geom2::BoundaryData2;
use crate::geom2::boundary::data::BData;
use crate::{Circle2, Line2, Point2, Result};

pub struct BCursor<'a> {
    data: &'a mut BoundaryData2,
    node_id: u32,
}


impl BoundaryData2 {
    pub fn add_seg_xy(&mut self, x: f64, y: f64) {
        self.push_data(BData::Seg((x, y))).unwrap();
    }

    pub fn add_seg(&mut self, p: &impl PCoords<2>) {
        self.add_seg_xy(p.coords().x, p.coords().y);
    }

    pub fn add_arc_xy(
        &mut self,
        center_x: f64,
        center_y: f64,
        end_x: f64,
        end_y: f64,
        clockwise: bool,
    ) {
        self.push_data(BData::Arc((center_x, center_y, end_x, end_y, clockwise))).unwrap();
    }

    pub fn add_arc(&mut self, center: &impl PCoords<2>, end: &impl PCoords<2>, clockwise: bool) {
        self.add_arc_xy(
            center.coords().x,
            center.coords().y,
            end.coords().x,
            end.coords().y,
            clockwise,
        );
    }

    pub fn add_corner_fillets(&mut self, corners: &[impl PCoords<2>], radius: f64) -> Result<()> {
        if corners.len() < 2 {
            return Err("Must have at least 2 corners to fillet".into());
        }

        for i in 0..corners.len() - 1 {
            let c0 = &corners[i];
            let c1 = &corners[i + 1];

            // First we'll find the circle at the tangent of the corner
            let v0 = self.last_point().coords() - c0.coords();
            let v1 = c1.coords() - c0.coords();
            let c = Circle2::tangent_to_corner(c0, &v0, &v1, radius)?;

            // Now we need to find the two tangent endpoints
            let l0 = Line2::from_points(&self.last_point(), c0).normalized();
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

            self.add_seg(&e0);
            self.add_arc(&c.center, &e1, clockwise);
        }

        let cf1 = &corners[corners.len() - 1];
        let cf0 = &corners[corners.len() - 2];
        let lf = Line2::from_points(cf0, cf1).normalized();
        if lf.scalar_project(&self.last_point()) < 0.0 {
            return Err("Invalid fillet radius, too large for corner".into());
        }
        self.add_seg(cf1);

        Ok(())
    }

    pub fn last_point(&self) -> Point2 {
        match self.last_data() {
            None => self.start.expect("open boundary must have a start point"),
            Some(BData::Seg((x, y))) => Point2::new(*x, *y),
            Some(BData::Arc((_cx, _cy, ex, ey, _))) => Point2::new(*ex, *ey),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn corner_fillet_single() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        data.add_corner_fillets(&[Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)], 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 * 2.0 + 0.25 * PI / 2.0;

        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }

    #[test]
    fn corner_fillet_single_inverted() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        data.add_corner_fillets(&[Point2::new(1.0, 0.0), Point2::new(1.0, -1.0)], 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 * 2.0 + 0.25 * PI / 2.0;

        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }

    #[test]
    fn corner_fillet_double() -> Result<()> {
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let points = vec![[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
            .iter()
            .map(|p| Point2::from(*p))
            .collect::<Vec<_>>();
        data.add_corner_fillets(&points, 0.25)?;

        let geom = data.try_to_boundary()?;
        let expected_len = 0.75 + 0.5 + 0.75 + 2.0 * 0.25 * PI / 2.0;
        assert_relative_eq!(geom.length(), expected_len);
        Ok(())
    }
}
