use crate::common::Interval;
use crate::geom2::UnitVec2;
use crate::metrology::{SurfaceDeviation2, SurfaceDeviationSet2};
/// This module contains tools for perform GD&T line profiles on 2D curves. It relies on the
/// `Curve2` construct from the `geom2` module as a way of representing the nominal geometry, and
/// can represent actual geometry using either `Curve2` or a slice of `Point2` points.
/// Currently, this module is based on the ASME Y14.5 standard for line profiles, but eventually
/// if I acquire a copy of the ISO specification, I will implement that as well.
use crate::{Curve2, CurveStation2, Point2, SurfacePoint2};

pub fn point_curve2_deviation(station: &CurveStation2, point: &Point2) -> SurfaceDeviation2 {
    // TODO: is there a better way to handle corners?
    let sp = station.surface_point();
    let vector = point - station.point();
    let normal = if vector.norm() < 1e-6 {
        sp.normal
    } else {
        if vector.dot(&sp.normal) < 0.0 {
            UnitVec2::new_normalize(-vector)
        } else {
            UnitVec2::new_normalize(vector)
        }
    };

    SurfaceDeviation2::new(
        SurfacePoint2::new(station.point(), normal),
        vector.dot(&normal),
    )
}

pub fn line_surface_deviations(
    nominal: &Curve2,
    actual: &[Point2],
    interval: Option<Interval>,
) -> SurfaceDeviationSet2 {
    let mut result = SurfaceDeviationSet2::default();
    for p in actual {
        let closest = nominal.at_closest_to_point(p);
        if let Some(i) = interval {
            if !i.contains(closest.length_along()) {
                continue;
            }
        }

        result.push(point_curve2_deviation(&closest, p));
    }

    result
}
