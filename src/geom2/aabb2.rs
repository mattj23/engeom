//! Tools for working with axis-aligned bounding boxes in 2D space.

use crate::{AngleInterval, Circle2, Point2, Vector2};
use num_traits::FloatConst;

pub type Aabb2 = parry2d_f64::bounding_volume::Aabb;

pub fn circle_aabb2(center: &Point2, radius: f64) -> Aabb2 {
    let half_extents = Vector2::new(radius, radius);
    Aabb2::new(center - half_extents, center + half_extents)
}

pub fn arc_aabb2(circle: &Circle2, angle0: f64, angle: f64) -> Aabb2 {
    let mut a = vec![angle0, angle0 + angle];
    let check = AngleInterval::new(angle0, angle);

    for i in 0..4 {
        // The angle in global space
        let t_g = i as f64 * f64::FRAC_PI_2();
        if check.contains(t_g) {
            a.push(t_g);
        }
    }

    Aabb2::from_points(
        &a.iter()
            .map(|&t| circle.point_at_angle(t))
            .collect::<Vec<_>>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::linear_space;
    use crate::common::points::dist;
    use approx::assert_relative_eq;
    use rand::Rng;

    #[test]
    fn stress_test_arc_aabb2() {
        let pi = std::f64::consts::PI;
        let mut rng = rand::rng();

        for _ in 0..1000 {
            // let center = Point2::new(rng.random_range(-10.0..10.0), rng.random_range(-10.0..10.0));
            let center = Point2::new(0.0, 0.0);
            let radius = rng.random_range(0.0..10.0);
            let circle = Circle2::from_point(center, radius);
            let angle0 = rng.random_range(-pi..pi);
            let angle = rng.random_range(-2.0 * pi..2.0 * pi);

            let aabb = arc_aabb2(&circle, angle0, angle);

            let thetas = linear_space(angle0, angle0 + angle, 1000);
            let points = thetas
                .iter()
                .map(|&t| circle.point_at_angle(t))
                .collect::<Vec<_>>();
            let expected = Aabb2::from_points(&points);
            assert!(
                dist(&aabb.mins, &expected.mins) < 1.0e-4,
                "radius: {}, angle0: {}, angle: {}",
                radius,
                angle0,
                angle
            );
            assert!(
                dist(&aabb.maxs, &expected.maxs) < 1.0e-4,
                "radius: {}, angle0: {}, angle: {}",
                radius,
                angle0,
                angle
            );

            // assert_relative_eq!(aabb.mins, expected.mins, epsilon = 1.0e-4);
            // assert_relative_eq!(aabb.maxs, expected.maxs, epsilon = 1.0e-4);
        }
    }

    #[test]
    fn arc_aabb_failure_case1() {
        let r = 4.865468058323412;
        let a0 = -2.3109122718910333;
        let a = 1.806427857721145;

        // Construct the expected AABB
        let circle = Circle2::from_point(Point2::new(0.0, 0.0), r);
        let thetas = linear_space(a0, a0 + a, 1000);
        let points = thetas
            .iter()
            .map(|&t| circle.point_at_angle(t))
            .collect::<Vec<_>>();
        let expected = Aabb2::from_points(&points);

        // Construct the AABB using the function
        let aabb = arc_aabb2(&circle, a0, a);

        assert_relative_eq!(aabb.mins, expected.mins, epsilon = 1.0e-4);
        assert_relative_eq!(aabb.maxs, expected.maxs, epsilon = 1.0e-4);
    }

    #[test]
    fn arc_aabb_failure_case2() {
        let r = 2.3918324202302066;
        let a0 = -3.1342936470156273;
        let a = -4.057966682848482;

        // Construct the expected AABB
        let circle = Circle2::from_point(Point2::new(0.0, 0.0), r);
        let thetas = linear_space(a0, a0 + a, 1000);
        let points = thetas
            .iter()
            .map(|&t| circle.point_at_angle(t))
            .collect::<Vec<_>>();
        let expected = Aabb2::from_points(&points);

        // Construct the AABB using the function
        let aabb = arc_aabb2(&circle, a0, a);

        assert_relative_eq!(aabb.mins, expected.mins, epsilon = 1.0e-4);
        assert_relative_eq!(aabb.maxs, expected.maxs, epsilon = 1.0e-4);
    }
}
