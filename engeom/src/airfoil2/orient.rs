use crate::Resample::ByCount;
use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::{Line2, Result, Vector2};

/// An enum that specifies a method for determining the orientation of the forward/aft direction of
/// an airfoil after the inscribed circles are detected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrientFwdAft {
    /// Detect the forward/aft orientation of the airfoil by detecting which end of the camber line
    /// the largest inscribed circle (the maximum thickness) is closest to. That end will become
    /// the leading edge.
    TmaxFwd,

    /// Detect the forward/aft orientation of the airfoil by specifying the airflow direction. The
    /// end of the camber line that is more negative in this direction will become the leading
    /// edge.
    Airflow(Vector2),

    /// Detect the forward/aft orientation of the airfoil by specifying a vector pointing in the
    /// direction of the leading edge. For example, if an airfoil section is oriented so that the
    /// +X direction goes from leading edge to trailing edge (trailing edge points are at a higher
    /// X value than the leading edge points), specifying `Fwd(-Vector2::x())` would correctly
    /// identify the leading edge.
    Fwd(Vector2),
}

impl OrientFwdAft {
    /// Apply the orientation method to a vec of inscribed circles, returning a vec where the first
    /// element is towards the leading edge of the airfoil and the last element is towards the
    /// trailing edge.
    ///
    /// # Arguments
    ///
    /// * `circles`: the vec of inscribed circles to consume and re-order
    ///
    /// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
    pub fn apply(&self, circles: Vec<Inscribed>) -> Result<Vec<Inscribed>> {
        let mut working = InscribedVec::new(circles);
        working.throw_if_less_than(2)?;

        if match self {
            OrientFwdAft::TmaxFwd => check_tmax_fwd(&working)?,
            OrientFwdAft::Airflow(dir) => check_dir_fwd(&working, -dir)?,
            OrientFwdAft::Fwd(dir) => check_dir_fwd(&working, *dir)?,
        } {
            working.reverse_order_only();
        }

        Ok(working.take_vec())
    }
}

/// An enum that specifies a method for determining the orientation of the upper (suction) and lower
/// (pressure) sides of an airfoil after the inscribed circles are detected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrientUpperLower {
    /// Attempt to automatically detect the upper and lower sides of the airfoil by looking at the
    /// curvature of the camber line. The more concave direction of the camber line will become
    /// the lower side of the airfoil.
    Curvature,

    /// Use the direction provided to identify the upper side of the airfoil. The side whose points
    /// are more positive in this direction will become the upper side.
    Upper(Vector2),

    /// Use the direction provided to identify the lower side of the airfoil. The side whose points
    /// are more positive in this direction will become the lower side.
    Lower(Vector2),
}

impl OrientUpperLower {
    /// Apply the orientation method to a vec of inscribed circles, returning a vec where the `p0`
    /// points of each circle are on the lower (pressure) side and the `p1` point are on the
    /// upper (suction) side of the airfoil.
    ///
    /// # Arguments
    ///
    /// * `circles`: the vec of inscribed circles to consume and re-order
    ///
    /// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
    pub fn apply(&self, circles: Vec<Inscribed>) -> Result<Vec<Inscribed>> {
        let mut working = InscribedVec::new(circles);
        working.throw_if_less_than(2)?;

        if match self {
            OrientUpperLower::Curvature => check_curvature_upper(&working)?,
            OrientUpperLower::Upper(dir) => check_dir_upper(&working, *dir)?,
            OrientUpperLower::Lower(dir) => check_dir_upper(&working, -dir)?,
        } {
            working.reverse_points();
        }

        Ok(working.take_vec())
    }
}

fn check_dir_upper(working: &InscribedVec, dir: Vector2) -> Result<bool> {
    let mut vote_flip = 0;
    for c in working.iter() {
        if !dir.dot(&c.contact_dir()).is_sign_positive() {
            vote_flip += 1;
        }
    }

    Ok(vote_flip > working.len() - vote_flip)
}

/// Try to detect the upper and lower directions of the airfoil by seeing how many points fall to
/// which side of a chord line.
fn check_curvature_upper(working: &InscribedVec) -> Result<bool> {
    let tol = working.average_spacing()? * 5e-2;
    let camber = working.camber_curve(tol)?;
    let points = camber.resample(ByCount(50))?.points().to_vec();
    let line = Line2::from_points(&points[0], &points[points.len() - 1]);

    let mut to_right = 0;
    for p in points.iter() {
        if line.signed_distance_to(p).is_sign_positive() {
            to_right += 1;
        }
    }

    // The side which more points fall onto is the direction of the upper side
    let upper_to_right = to_right > points.len() - to_right;
    let mut vote_flip = 0;
    for c in working.iter() {
        let sp = camber.at_closest_to_point(&c.center()).surface_point();

        // If we are expecting the upper side to be on the right, the `p1` point should have a
        // positive signed distance to the surface point of the camber line. Otherwise, it should
        // be negative.
        if upper_to_right != sp.scalar_projection(&c.p1).is_sign_positive() {
            vote_flip += 1;
        }
    }

    Ok(vote_flip > working.len() - vote_flip)
}

fn check_dir_fwd(working: &InscribedVec, dir: Vector2) -> Result<bool> {
    let (front, back) = working.front_and_back()?;
    Ok(dir.dot(&front.center().coords) < dir.dot(&back.center().coords))
}

fn check_tmax_fwd(working: &InscribedVec) -> Result<bool> {
    let tol = working.average_spacing()? * 5e-2;
    let camber = working.camber_curve(tol)?;
    let i_max = working.index_of_tmax();

    let c = camber.at_closest_to_point(&working[i_max].c.center);
    if c.length_along() < camber.length() * 0.5 {
        Ok(false)
    } else {
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Circle2, Point2};

    fn ic(x: f64, y: f64, r: f64) -> Inscribed {
        Inscribed::new(
            Circle2::new(x, y, r),
            Point2::new(x, y - r),
            Point2::new(x, y + r),
        )
    }

    fn quick() -> InscribedVec {
        InscribedVec::new(vec![
            ic(-1.0, 0.0, 1.0),
            ic(0.0, 0.0, 2.0),
            ic(2.0, 0.0, 1.0),
        ])
    }

    #[test]
    fn tmax_fwd() -> Result<()> {
        let mut data = quick();
        assert!(!check_tmax_fwd(&data)?);

        data.reverse_order();
        assert!(check_tmax_fwd(&data)?);

        Ok(())
    }

    #[test]
    fn dir_fwd() -> Result<()> {
        let mut data = quick();
        let v = Vector2::x();
        assert!(check_dir_fwd(&data, v)?);

        data.reverse_order();
        assert!(!check_dir_fwd(&data, v)?);
        Ok(())
    }

    #[test]
    fn dir_upper() -> Result<()> {
        let mut data = quick();
        let v = Vector2::y();
        assert!(!check_dir_upper(&data, v)?);

        data.reverse_order();
        assert!(check_dir_upper(&data, v)?);
        Ok(())
    }
}
