use super::helpers::{curve_from_inscribed_circles, find_tmax_circle, reverse_inscribed_circles};
use crate::airfoil::InscribedCircle;
use crate::{Curve2, Result, Vector2};

/// This trait defines an interface to perform orientation of the camber line of an airfoil
/// section, specifically referring to its order and its relationship to the leading edge. A
/// camber line begins at the leading edge and ends at the trailing edge.
pub trait CamberOrientation {
    /// Orient the camber line of an airfoil section based on the method specified by the
    /// implementation. The camber line at this stage is a series of inscribed circles whose
    /// adjacency in space is coupled to their ordering in the container. However, the orientation
    /// (whether the first circle is at the leading edge or the trailing edge) is not yet known.
    ///
    /// This method will return a new container of inscribed circles with the camber line oriented
    /// so that the first circle is closest to the leading edge and the last circle is closest to
    /// the trailing edge. The order of the circles will be preserved.
    ///
    /// This method will take ownership of the input container and return a new container with the
    /// circles in the correct order.
    ///
    /// # Arguments
    ///
    /// * `section`: the airfoil section curve used to generate the inscribed circles
    /// * `stations`: the inscribed circles in the airfoil section
    ///
    /// returns: Vec<InscribedCircle, Global>
    fn orient_camber_line(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
    ) -> Result<Vec<InscribedCircle>>;
}

/// This struct implements the `CamberOrientation` trait and orients the camber line of an airfoil
/// based on the principle that for most subsonic airfoils the point of max thickness is closer to
/// the leading edge than to the trailing edge along the camber line.  It will scan the inscribed
/// circles for the largest diameter, and evaluate the distance of the circle center along the
/// known camber line.  If that distance is closer to the end than the beginning, the order of the
/// circles will be reversed.
pub struct TMaxFwd {}

impl TMaxFwd {
    pub fn new() -> Self {
        TMaxFwd {}
    }

    /// Create a new boxed instance of the `TMaxFwd` struct.
    pub fn make() -> Box<dyn CamberOrientation> {
        Box::new(TMaxFwd::new())
    }
}

impl CamberOrientation for TMaxFwd {
    fn orient_camber_line(
        &self,
        _section: &Curve2,
        stations: Vec<InscribedCircle>,
    ) -> Result<Vec<InscribedCircle>> {
        // We're going to edit this in place, so we'll move it to a mutable variable.
        let mut moved = stations;

        // Construct the camber line curve
        let camber = curve_from_inscribed_circles(&moved, 1e-4)?;

        // Find the largest diameter circle
        let tmax = find_tmax_circle(&moved)
            .ok_or("Empty inscribed circles container during tmax forward orientation.")?;

        // How far along the camber line is the center of the circle?
        let l = camber.at_closest_to_point(&tmax.center()).length_along();

        // What fraction of the camber line is this distance?
        let fraction = l / camber.length();

        if fraction > 0.5 {
            // The circle is closer to the trailing edge than the leading edge, so we need to
            // reverse the order of the circles.
            reverse_inscribed_circles(&mut moved);
        }

        Ok(moved)
    }
}

/// This struct implements the `CamberOrientation` trait and orients the camber line of an airfoil
/// such that the leading edge is the camber line's end that is more in the direction of a
/// specified vector. This is useful if you know the direction of the leading edge already.
pub struct DirectionFwd {
    pub direction: Vector2,
}

impl DirectionFwd {
    pub fn new(direction: Vector2) -> Self {
        DirectionFwd { direction }
    }

    /// Create a new boxed instance of the `DirectionFwd` struct.
    pub fn make(direction: Vector2) -> Box<dyn CamberOrientation> {
        Box::new(DirectionFwd::new(direction))
    }
}

impl CamberOrientation for DirectionFwd {
    fn orient_camber_line(
        &self,
        _section: &Curve2,
        stations: Vec<InscribedCircle>,
    ) -> Result<Vec<InscribedCircle>> {
        // We're going to edit this in place, so we'll move it to a mutable variable.
        let mut moved = stations;

        let c0 = moved
            .first()
            .ok_or("Empty inscribed circles container during camber orientation.")?
            .circle
            .center;

        let c1 = moved
            .last()
            .ok_or("Empty inscribed circles container during camber orientation.")?
            .circle
            .center;

        if self.direction.dot(&c0.coords) < self.direction.dot(&c1.coords) {
            // The circle is closer to the trailing edge than the leading edge, so we need to
            // reverse the order of the circles.
            reverse_inscribed_circles(&mut moved);
        }

        Ok(moved)
    }
}
