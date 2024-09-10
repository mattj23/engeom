//! This module contains a structure for representing an inscribed circle within an airfoil
//! section, consisting of a center point and a radius, and having the condition of having
//! two points on the circle which are tangent to the airfoil section while being otherwise
//! contained within the section.  Inscribed circles are used to calculate the camber line.

use serde::Serialize;
use crate::{Circle2, Point2};
use crate::geom2::polyline2::SpanningRay;

#[derive(Serialize)]
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
    pub fn new(spanning_ray: SpanningRay, upper: Point2, lower: Point2, circle: Circle2) -> InscribedCircle {
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

    pub fn radius(&self) -> f64 {
        self.circle.r()
    }

    pub fn center(&self) -> Point2 {
        self.circle.center
    }

}