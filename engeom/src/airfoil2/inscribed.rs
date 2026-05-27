//! Common tools for finding and working with inscribed circles on airfoil sections.

use std::ops::Index;
use crate::airfoil2::SectionInput;
use crate::common::dist;
use crate::geom2::{rot90, LineOps2};
use crate::AngleDir::Ccw;
use crate::{AngleDir, Circle2, Line2, Point2, SurfacePoint2};
use serde::{Deserialize, Serialize};

/// Represents an inscribed circle inside an airfoil section. Contains the circle itself (center
/// and radius) and the two contact point with the perimeter of the section. The circle center is
/// definitionally on the mean camber line. The two points are not guaranteed to be in any order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inscribed {
    /// The inscribed circle, centered on the camber line with a radius that causes it to contact
    /// the section at two contact points on opposite sides of the camber line.
    pub c: Circle2,

    /// The first contact point, no guarantee on whether this is on the upper or lower surface, but
    /// it will be on the opposite side of the camber line from `p1`.
    pub p0: Point2,

    /// The second contact point, no guarantee on whether this is on the upper or lower surface, but
    /// it will be on the opposite side of the camber line from `p0`.
    pub p1: Point2,
}

impl Inscribed {
    pub fn new(c: Circle2, p0: Point2, p1: Point2) -> Self {
        Self { c, p0, p1 }
    }

    /// Create a SurfacePoint2 located at the inscribed circle center and pointing along the
    /// estimated direction of the camber line. The camber line direction is estimated by taking
    /// the unit vector from `p0` to `p1` and rotating it 90 degrees clockwise.
    pub fn camber_point(&self) -> SurfacePoint2 {
        let d = rot90(AngleDir::Cw) * (self.p1 - self.p0);
        SurfacePoint2::new_normalize(self.c.center, d)
    }

    /// Swap the points in this circle, so that `p0` becomes the old `p1` and vice versa. This will
    /// change the direction of the `.camber_point()` function.
    pub fn reverse_points(&mut self) {
        (self.p0, self.p1) = (self.p1, self.p0)
    }
}

/// This is a convenience container for a sequence of [`Inscribed`] circle entities, providing tools
/// for performing some collection-specific operations.
pub struct InscribedVec {
    items: Vec<Inscribed>,
}

impl InscribedVec {
    /// Create an empty collection
    pub fn empty() -> Self {
        Self { items: vec![] }
    }

    /// Create a new collection by taking ownership of a Vec of inscribed circles. The order is
    /// preserved exactly as it's received.
    pub fn new(items: Vec<Inscribed>) -> Self {
        Self { items }
    }

    /// Goes through all the inscribed circles in the collection and swaps their `p0` and `p1`
    /// points. This can be used to change the orientation of the upper/lower surface points once
    /// the direction has been established.
    pub fn reverse_points(&mut self) {
        for item in &mut self.items {
            item.reverse_points()
        }
    }

    /// Reverses the order of the inscribed circles in the collection, and also reverses the points
    /// in each circle so that the camber point direction is preserved. This is useful for switching
    /// between leading-to-trailing and trailing-to-leading orderings of the circles.
    pub fn reverse_order(&mut self) {
        self.items.reverse();
        self.reverse_points();
    }

    /// Remove the last item from the collection and return it, or return `None` if the collection
    /// is empty.
    pub fn pop(&mut self) -> Option<Inscribed> {
        self.items.pop()
    }

    /// Take ownership of the internal vector of inscribed circles, consuming the collection. The
    /// order is preserved exactly as it is in the collection.
    pub fn take_vec(self) -> Vec<Inscribed> {
        self.items
    }

    /// Return a reference to the last item in the collection, or return `None` if the collection
    /// is empty.
    pub fn last(&self) -> Option<&Inscribed> {
        self.items.last()
    }

    /// Return a reference to the first item in the collection, or return `None` if the collection
    /// is empty.
    pub fn first(&self) -> Option<&Inscribed> {
        self.items.first()
    }

    /// Extend the collection by taking ownership of another collection and appending all of its
    /// items onto the end of this collection. The order of both collections are preserved.
    pub fn extend(&mut self, other: Self) {
        self.items.extend(other.items)
    }

    /// Add a new inscribed circle to the collection, performing any necessary refinement between
    /// the current end of the collection and the new item.
    pub fn refine_and_push(&mut self, item: Inscribed, input: &SectionInput) {
        if self.items.is_empty() {
            self.items.push(item);
        } else {
            let mut stack = vec![item];
            while let Some(top) = stack.pop() {
                let last = self.items.last().unwrap().clone();
                if let Some(to_push) = input.inscribed_refined(&last, &top) {
                    stack.push(top);
                    stack.push(to_push);
                } else {
                    self.items.push(top);
                }
            }
        }
    }

    /// Find the index of the largest inscribed circle in the collection
    pub fn index_of_tmax(&self) -> usize {
        self.items
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.c.r().partial_cmp(&b.c.r()).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}

impl Index<usize> for InscribedVec {
    type Output = Inscribed;
    fn index(&self, index: usize) -> &Self::Output {
        &self.items[index]
    }
}

impl<'a> SectionInput<'a> {

    /// Try to find an inscribed circle using the position and orientation of a search line. This is
    /// the equivalent of first finding a crossing line using `.crossing_line(...)` and then
    /// using `.inscribed_from_crossing(...)` to get the inscribed circle. If the first step fails
    /// this function returns a `None`, doing both in a single step.
    pub fn try_inscribed(&self, line: &impl LineOps2) -> Option<Inscribed> {
        let donor = Line2::from(line);
        let crossing = self.crossing_line(&donor)?;
        Some(self.inscribed_from_crossing(&crossing))
    }

    /// Checks if a refinement step is needed between inscribed circles `a` and `b` to bring the
    /// linear interpolation of radius and circle center to within the supplied tolerance. If
    /// refinement is necessary, the inscribed circle between the two stations is returned. If it
    /// isn't necessary, a `None` is returned instead.
    pub fn inscribed_refined(&self, a: &Inscribed, b: &Inscribed) -> Option<Inscribed> {
        let ca = a.camber_point();
        let cb = b.camber_point();

        // Get the interpolation fraction (may need to be varied)
        let f = 0.5;

        // Get the interpolated point and radius
        let cp = ca.slerp(&cb, f);
        let cr = (b.c.r() - a.c.r()) * f + a.c.r();

        // Now find the measured circle
        let line = Line2::from(&cp.rot_normal_90(Ccw));
        let circle = self.try_inscribed(&line)?;
        if (circle.c.r() - cr).abs() < self.general_tol
            && dist(&cp, &circle.c.center) < self.general_tol
        {
            None
        } else {
            Some(circle)
        }
    }

    /// From a donor line, try to find a valid line which crosses the airfoil section. That means
    /// that it has exactly two intersections with the section, and is returned such that t=0 and
    /// t=1 are both exactly on the section perimeter.
    ///
    /// If the donor line does not result in exactly two intersections this function returns a
    /// `None`
    ///
    /// # Arguments
    ///
    /// * `donor`: Any parametric line coincident with the desired crossing line, does not need to
    ///   start or end on the section perimeter itself.
    ///
    /// returns: Option<Line2>
    pub fn crossing_line(&self, donor: &Line2) -> Option<Line2> {
        let ts = self
            .section
            .intersections_with_line(&donor)
            .iter()
            .map(|(t, _)| *t)
            .collect::<Vec<_>>();
        if ts.len() != 2 {
            None
        } else {
            Some(Line2::from_points(&donor.at(ts[0]), &donor.at(ts[1])))
        }
    }

    /// Given a valid crossing line, find the inscribed circle center along the line with a
    /// binary search, resolving the position down to the `resolve_tol` of the input entity.
    ///
    /// # Arguments
    ///
    /// * `line`: A valid crossing line for the airfoil (t=0 and t=1 must be on the section
    ///   perimeter)
    ///
    /// returns: Inscribed
    pub fn inscribed_from_crossing(&self, line: &Line2) -> Inscribed {
        struct InscribedEnd {
            f: f64,
            d: f64,
            p: Point2,
        }

        impl InscribedEnd {
            pub fn new(f: f64, d: f64, p: Point2) -> Self {
                InscribedEnd { f, d, p }
            }
        }

        let mut pos = InscribedEnd::new(1.0, 0.0, line.at(1.0));
        let mut neg = InscribedEnd::new(0.0, 0.0, line.at(0.0));
        let mut working;

        // While the distance between the positive and negative search bounds is greater than the
        // tolerance, continue to search for the inscribed circle center.
        while (pos.f - neg.f) * line.direction.norm() > self.resolve_tol {
            // TODO: Once we are close, there may be a geometric way to skip a bunch of iterations
            // We will update the working point to be right in the middle of the positive and negative
            // direction limits.
            let fraction = (pos.f + neg.f) * 0.5;
            working = line.at(fraction);

            // Now we find the closest position on the curve to the working point, and calculate the
            // distance and direction to that point. The direction will be used to determine which side
            // of the limits we will adjust.
            let closest = self.section.at_closest_to_point(&working);
            let to_closest = closest.point() - working; // The direction vector to the closest point
            let distance = dist(&working, &closest.point());
            let update = InscribedEnd::new(fraction, distance, closest.point());

            // If the direction vector to the closest point is in the positive direction of the ray,
            // then we will adjust the positive limit.  Otherwise, we will adjust the negative limit.
            if to_closest.dot(&line.direction) > 0.0 {
                pos = update;
            } else {
                neg = update;
            }
        }

        // Finally, we will put the center of the inscribed circle at the midpoint of the positive and
        // negative limits, splitting the difference one last time, and we will set the radius to be
        // the average of the positive and negative distances. By this point the difference will be
        // below the tolerance value.
        let c = Circle2::from_point(line.at((pos.f + neg.f) * 0.5), (pos.d + neg.d) * 0.5);

        Inscribed::new(c, neg.p, pos.p)
    }
}
