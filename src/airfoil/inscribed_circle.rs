//! This module contains a structure for representing an inscribed circle within an airfoil
//! section, consisting of a center point and a radius, and having the condition of having
//! two points on the circle which are tangent to the airfoil section while being otherwise
//! contained within the section.  Inscribed circles are used to calculate the camber line.

use crate::common::angle_in_direction;
use crate::common::points::linear_interpolation_error;
use crate::geom2::polyline2::SpanningRay;
use crate::geom2::{rot90, UnitVec2};
use crate::{AngleDir, Arc2, Circle2, Point2, SurfacePoint2};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
pub struct InscribedCircle {
    /// The spanning ray which crosses the airfoil section, on which the circle center
    /// is located.
    pub spanning_ray: SpanningRay,

    /// The contact point of the circle with the upper surface of the airfoil section.
    pub upper: Point2,

    /// The contact point of the circle with the lower surface of the airfoil section.
    pub lower: Point2,

    /// The circle that is inscribed within the airfoil section
    pub circle: Circle2,
    // pub thk: f64,
}

impl InscribedCircle {
    pub fn new(
        spanning_ray: SpanningRay,
        upper: Point2,
        lower: Point2,
        circle: Circle2,
    ) -> InscribedCircle {
        // let thk = dist(&upper, &lower);
        InscribedCircle {
            spanning_ray,
            upper,
            lower,
            circle,
            // thk,
        }
    }

    pub fn reversed(&self) -> InscribedCircle {
        InscribedCircle::new(
            self.spanning_ray.reversed(),
            self.lower,
            self.upper,
            self.circle,
        )
    }

    /// Return the clockwise arc that connects the upper and lower contact points of the inscribed
    /// circle, with the arc half returned being the one which has its center angle in the
    /// direction of the given vector.
    ///
    /// # Arguments
    ///
    /// * `direction`:
    ///
    /// returns: Arc2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn contact_arc(&self, direction: &UnitVec2) -> Arc2 {
        let angle_u = self.circle.angle_of_point(&self.upper);
        let angle_l = self.circle.angle_of_point(&self.lower);

        // Arc from upper to lower
        let arc0 = self
            .circle
            .to_partial_arc(angle_u, angle_in_direction(angle_u, angle_l, AngleDir::Ccw));

        // Arc from lower to upper
        let arc1 = self
            .circle
            .to_partial_arc(angle_l, angle_in_direction(angle_l, angle_u, AngleDir::Ccw));

        // Which of the two arcs will we use?
        let v0 = arc0.point_at_fraction(0.5) - self.circle.center;
        let v1 = arc1.point_at_fraction(0.5) - self.circle.center;

        if direction.dot(&v0) > direction.dot(&v1) {
            arc0
        } else {
            arc1
        }
    }

    pub fn radius(&self) -> f64 {
        self.circle.r()
    }

    pub fn center(&self) -> Point2 {
        self.circle.center
    }

    /// Calculates a point at the inscribed circle's center facing in the direction of the camber
    /// line. The direction is found by noting that the vector from the upper to lower contact
    /// points is perpendicular to the direction of the camber line at the inscribed circle's center.
    pub fn camber_point(&self) -> SurfacePoint2 {
        let dir = rot90(AngleDir::Cw) * (self.upper - self.lower);
        SurfacePoint2::new_normalize(self.center(), dir)
    }

    /// Calculates the interpolation error of the inscribed circle with respect to its two
    /// neighbors. The interpolation error is the distance that the inscribed circle deviates
    /// from the line segment connecting its two neighbors, and can be thought of as the error that
    /// would be present if the test circle did not exist.
    ///
    /// The error being checked is the maximum of the interpolation errors of the circle's center,
    /// upper contact point, and lower contact point.
    ///
    /// # Arguments
    ///
    /// * `s0`: the previous inscribed circle
    /// * `s1`: the next inscribed circle
    ///
    /// returns: f64
    pub fn interpolation_error(&self, s0: &Self, s1: &Self) -> f64 {
        let upper = linear_interpolation_error(&s0.upper, &s1.upper, &self.upper);
        let lower = linear_interpolation_error(&s0.lower, &s1.lower, &self.lower);
        let center = linear_interpolation_error(&s0.center(), &s1.center(), &self.center());

        upper.max(lower).max(center)
    }

    /// Reverses the direction of the spanning ray and swaps the upper and lower contact points.
    pub fn reverse_in_place(&mut self) {
        self.spanning_ray = self.spanning_ray.reversed();
        std::mem::swap(&mut self.upper, &mut self.lower);
    }
}
