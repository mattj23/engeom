use crate::common::PCoords;
use crate::geom2::signed_angle;
use crate::geom3::IsoExtensions3;
use crate::{Iso3, Point3, Result, To2D, Vector2, Vector3};

/// A transformation between a general 3D cartesian space and a 3D cartesian space representing the
/// surface of a cylinder. The process of going from general 3D space to the cylinder surface is
/// analogous to unwrapping the surface of the cylinder onto the X-Y plane, but preserving the
/// relief in the radial direction as Z. The inverse process wraps the area around the X-Y plane
/// back onto the cylinder surface.
pub struct UnrollTransform {
    /// Transforms from cylinder aligned space to general 3D space
    iso: Iso3,

    /// Transforms from general 3D space to cylinder aligned space
    inv: Iso3,

    /// The radius of the cylinder surface
    pub radius: f64,
}

impl UnrollTransform {
    pub fn try_new(
        center_axis: Vector3,
        polar_axis: Vector3,
        radius: f64,
        origin: Option<Point3>,
    ) -> Result<UnrollTransform> {
        let iso = Iso3::try_from_basis_zx(&center_axis, &polar_axis, origin)?;
        let inv = iso.inverse();
        Ok(UnrollTransform { iso, inv, radius })
    }

    pub fn to_world(&self, point: &impl PCoords<3>) -> Point3 {
        // Find the angle around the cylinder based on the x coordinate divided by the radius
        let theta = point.coords().x / self.radius;
        let r = self.radius + point.coords().z;
        let p = Point3::new(r * theta.cos(), r * theta.sin(), point.coords().y);
        self.iso * p
    }

    pub fn to_surface(&self, point: &impl PCoords<3>) -> Point3 {
        let p = self.inv * Point3::from(point.coords());
        // The y coordinate is the transformed z, the z is the distance from the z axis minus the
        // radius, and the x value is the arc length around the cylinder at the nominal radius
        let xy = p.to_2d().coords;
        let theta = signed_angle(&Vector2::x(), &xy);

        Point3::new(theta * self.radius, p.z, xy.norm() - self.radius)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::XyzWpr;
    use crate::{Iso3, Vector3};
    use approx::assert_relative_eq;
    use rand::Rng;
    use std::f64::consts::PI;

    #[test]
    fn origin_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(1.0, 0.0, 0.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::origin(), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_reverse() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(1.0, 0.0, 0.0);
        let surface = unroll.to_world(&Point3::origin());
        assert_relative_eq!(world, surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_z_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(1.1, 0.0, 0.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::new(0.0, 0.0, 0.1), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_z_reverse() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let surface = Point3::new(0.0, 0.0, 0.1);
        let world = unroll.to_world(&surface);
        assert_relative_eq!(Point3::new(1.1, 0.0, 0.0), world, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_minus_z_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(0.9, 0.0, 0.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::new(0.0, 0.0, -0.1), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_x_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(0.0, 1.0, 0.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::new(PI / 2.0, 0.0, 0.0), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_x_reverse() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let surface = Point3::new(PI / 2.0, 0.0, 0.0);
        let world = unroll.to_world(&surface);
        assert_relative_eq!(Point3::new(0.0, 1.0, 0.0), world, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_minus_x_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(0.0, -1.0, 0.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::new(-PI / 2.0, 0.0, 0.0), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_y_forward() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let world = Point3::new(1.0, 0.0, 1.0);
        let surface = unroll.to_surface(&world);
        assert_relative_eq!(Point3::new(0.0, 1.0, 0.0), surface, epsilon = 1.0e-5);
    }

    #[test]
    fn origin_plus_y_reverse() {
        let unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let surface = Point3::new(0.0, 1.0, 0.0);
        let world = unroll.to_world(&surface);
        assert_relative_eq!(Point3::new(1.0, 0.0, 1.0), world, epsilon = 1.0e-5);
    }

    #[test]
    fn transformed_round_trip() {
        let std_unroll = UnrollTransform::try_new(Vector3::z(), Vector3::x(), 1.0, None).unwrap();
        let std_surface = Point3::origin();
        let std_world = std_unroll.to_world(&std_surface);

        assert_relative_eq!(Point3::new(1.0, 0.0, 0.0), std_world, epsilon = 1.0e-5);

        let transform = Iso3::from_rz(PI / 4.0);
        let unroll = UnrollTransform::try_new(
            transform * Vector3::z(),
            transform * Vector3::x(),
            1.0,
            Some(transform * Point3::origin()),
        )
        .unwrap();

        let exp_world = transform * std_world;
        let world = unroll.to_world(&std_surface);

        assert_relative_eq!(exp_world, world, epsilon = 1.0e-5);
    }

    #[test]
    fn stress_test_round_trip() {
        let mut rng = rand::rng();

        for _ in 0..1000 {
            let radius = rng.random_range(0.5..3.0);
            let max_x = radius * PI * 0.99;
            let max_z = radius * 0.9;
            let std_unroll =
                UnrollTransform::try_new(Vector3::z(), Vector3::x(), radius, None).unwrap();

            // Generate a random surface point
            let surface = Point3::new(
                rng.random_range(-max_x..max_x),
                rng.random_range(-5.0..5.0),
                rng.random_range(-max_z..max_z),
            );

            // Figure out what it is in the standard world
            let std_world = std_unroll.to_world(&surface);

            // Now bring it back and verify we get the same surface point
            let std_surface = std_unroll.to_surface(&std_world);
            assert_relative_eq!(surface, std_surface, epsilon = 1.0e-5);

            // Now we'll create a randomly positioned unroll transform and verify the same round
            // trip. The surface points should remain the same, but the world points will change.
            let p2 = PI / 2.0;
            let values = XyzWpr::new(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-p2..p2),
                rng.random_range(-p2..p2),
                rng.random_range(-p2..p2),
            );
            let iso = Iso3::from(&values);

            let rnd_unroll = UnrollTransform::try_new(
                iso * Vector3::z(),
                iso * Vector3::x(),
                radius,
                Some(iso * Point3::origin()),
            )
            .unwrap();

            let rnd_world_ex = iso * std_world;
            let rnd_world = rnd_unroll.to_world(&surface);

            assert_relative_eq!(rnd_world_ex, rnd_world, epsilon = 1.0e-5);

            let rnd_surface = rnd_unroll.to_surface(&rnd_world);
            assert_relative_eq!(surface, rnd_surface, epsilon = 1.0e-5);
        }
    }
}
