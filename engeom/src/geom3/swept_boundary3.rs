use crate::common::PCoords;
use crate::geom2::Boundary2;
use crate::geom3::align3::{AlignSurfMatch3, SurfaceTarget3};
use crate::{Iso3, Point3, To2D, To3D, UnitVec3};

pub struct ExtrudedBoundary3 {
    /// A boundary in the X-Y plane which is the basis of the sweep. It does not matter if the
    /// boundary is closed or open.
    shape: Boundary2,

    /// A transformation that brings the X-Y boundary to an arbitrary position in 3D space.
    start: Iso3,

    /// The inverse of the start isometry
    start_inv: Iso3,

    /// The length that the boundary is extruded along the `start` isometry's Z axis in the 3D
    /// world coordinates.
    length: f64,
}

impl ExtrudedBoundary3 {
    pub fn new(shape: Boundary2, start: Iso3, length: f64) -> Self {
        let start_inv = start.inverse();
        Self {
            shape,
            start,
            start_inv,
            length,
        }
    }
}

impl SurfaceTarget3 for ExtrudedBoundary3 {
    fn align_surf_closest_to(&self, p: &impl PCoords<3>) -> AlignSurfMatch3 {
        let p = Point3::from(p.coords());

        // First we want to bring the test point into the local coordinates of the start isometry
        // so we can work with the 2D boundary.
        let lp3 = self.start_inv * p;
        let lp2 = lp3.to_2d();

        // Now we get the closest point on the boundary and create the local 3d point and normal
        let (_, m) = self.shape.at_closest_to_point(&lp2);
        let lc3 = Point3::new(m.point.x, m.point.y, lp3.z);
        let ln3 = UnitVec3::new_normalize(m.normal.into_inner().to_3d());

        // Now we figure out if we're off the boundary by either missing the z interval or by
        // projecting onto a boundary end.
        let is_off = lp3.z < 0.0
            || lp3.z > self.length
            || (!self.shape.is_closed()
                && (m.l <= f64::EPSILON || m.l >= self.shape.length() - f64::EPSILON));

        AlignSurfMatch3::new(self.start * lc3, self.start * ln3, !is_off, 1.0)
    }
}

pub struct PolarBoundary3 {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::geom2::{BoundaryData2, BoundaryEditor};
    use crate::geom3::tests::RandomGeometry;
    use approx::assert_relative_eq;

    fn open_extruded() -> ExtrudedBoundary3 {
        let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
        data.add_seg_xy(0.0, 0.0);
        data.add_seg_xy(0.0, 1.0);
        let b = data.try_to_boundary().unwrap();
        ExtrudedBoundary3::new(b, Iso3::identity(), 1.0)
    }

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn extruded_points_off() {
        let target = open_extruded();

        assert_eq!(false, target.align_surf_closest_to(&p(2.0, 0.0, 0.5)).is_on);
        assert_eq!(false, target.align_surf_closest_to(&p(0.0, 2.0, 0.5)).is_on);
        assert_eq!(
            false,
            target.align_surf_closest_to(&p(0.0, 0.0, -1.0)).is_on
        );
        assert_eq!(false, target.align_surf_closest_to(&p(0.0, 0.0, 2.0)).is_on);
    }

    #[test]
    fn extruded_points_on() {
        let target = open_extruded();
        assert_eq!(true, target.align_surf_closest_to(&p(0.5, 0.1, 0.5)).is_on);
        assert_eq!(true, target.align_surf_closest_to(&p(0.1, 0.5, 0.5)).is_on);
    }

    #[test]
    fn extruded_known_points() {
        let pairs = vec![
            ([0.0, 2.0, 0.5], [0.0, 1.0, 0.5], None),
            ([2.0, 0.0, 0.5], [1.0, 0.0, 0.5], None),
            ([0.5, 0.1, 0.3], [0.5, 0.0, 0.3], Some([0.0, 1.0, 0.0])),
            ([0.1, 0.5, 0.3], [0.0, 0.5, 0.3], Some([1.0, 0.0, 0.0])),
            ([0.0, 2.0, -1.0], [0.0, 1.0, -1.0], None),
            ([2.0, 0.0, -1.0], [1.0, 0.0, -1.0], None),
            ([0.0, 2.0, 2.0], [0.0, 1.0, 2.0], None),
            ([2.0, 0.0, 2.0], [1.0, 0.0, 2.0], None),
            ([-1.0, -1.0, 0.5], [0.0, 0.0, 0.5], None),
            ([-1.0, -1.0, -1.0], [0.0, 0.0, -1.0], None),
            ([-1.0, -1.0, 2.0], [0.0, 0.0, 2.0], None),
        ];

        let mut rg = RandomGeometry::new();

        for _ in 0..100 {
            let iso = rg.iso3(10.0);

            let mut data = BoundaryData2::new_open_xy(1.0, 0.0);
            data.add_seg_xy(0.0, 0.0);
            data.add_seg_xy(0.0, 1.0);
            let b = data.try_to_boundary().unwrap();
            let target = ExtrudedBoundary3::new(b, iso.clone(), 1.0);

            for (t, e, n) in pairs.iter() {
                let t = iso * Point3::from(*t);
                let e = iso * Point3::from(*e);

                let c = target.align_surf_closest_to(&t);
                assert_relative_eq!(c.point, e, epsilon = 1.0e-6);

                if let Some(n) = n {
                    let n = iso * UnitVec3::new_normalize(Vector3::from_column_slice(n));
                    assert_relative_eq!(c.normal, n, epsilon = 1.0e-6);
                }
            }
        }
    }
}
