use crate::common::SurfacePoint;
use parry3d_f64::na::Point;
use serde::{Deserialize, Serialize};
use std::ops::{Deref, Index};

/// A `SurfaceDeviation` is a struct which is used to represent a point on a reference surface
/// (n-1 dimensional) in n-dimensional space and a deviation from that surface in the direction
/// of the surface normal. The deviation is a scalar value which is the distance from the surface
/// point to the actual point.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SurfaceDeviation<const D: usize> {
    pub surface: SurfacePoint<D>,
    pub deviation: f64,
}

impl<const D: usize> SurfaceDeviation<D> {
    pub fn new(surface: SurfacePoint<D>, deviation: f64) -> Self {
        Self { surface, deviation }
    }

    pub fn actual_point(&self) -> Point<f64, D> {
        self.surface.at_distance(self.deviation)
    }
}

/// A `SurfaceDeviationSet` is a struct which is used to represent a collection of
/// `SurfaceDeviation` items, such as a set of deviations from a reference surface. This struct
/// implements methods to find the maximum and minimum deviations in the set.
#[derive(Clone, Serialize, Deserialize)]
pub struct SurfaceDeviationSet<const D: usize> {
    values: Vec<SurfaceDeviation<D>>,
    max_index: Option<usize>,
    min_index: Option<usize>,
}

impl<const D: usize> Default for SurfaceDeviationSet<D> {
    fn default() -> Self {
        Self {
            values: Vec::new(),
            max_index: None,
            min_index: None,
        }
    }
}

impl<const D: usize> Index<usize> for SurfaceDeviationSet<D> {
    type Output = SurfaceDeviation<D>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.values[index]
    }
}

impl<const D: usize> Deref for SurfaceDeviationSet<D> {
    type Target = [SurfaceDeviation<D>];

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl<const D: usize> SurfaceDeviationSet<D> {
    pub fn new(values: Vec<SurfaceDeviation<D>>) -> Self {
        let max_index = values
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.deviation.partial_cmp(&b.deviation).unwrap())
            .map(|(i, _)| i);
        let min_index = values
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.deviation.partial_cmp(&b.deviation).unwrap())
            .map(|(i, _)| i);
        Self {
            values,
            max_index,
            min_index,
        }
    }

    pub fn push_new(&mut self, surface: SurfacePoint<D>, deviation: f64) {
        self.push(SurfaceDeviation::new(surface, deviation));
    }

    pub fn push(&mut self, deviation: SurfaceDeviation<D>) {
        if self.max_index.is_none()
            || deviation.deviation > self.values[self.max_index.unwrap()].deviation
        {
            self.max_index = Some(self.values.len());
        }

        if self.min_index.is_none()
            || deviation.deviation < self.values[self.min_index.unwrap()].deviation
        {
            self.min_index = Some(self.values.len());
        }

        self.values.push(deviation);
    }

    pub fn iter(&self) -> impl Iterator<Item = &SurfaceDeviation<D>> {
        self.values.iter()
    }

    pub fn max(&self) -> Option<&SurfaceDeviation<D>> {
        self.max_index.map(|i| &self.values[i])
    }

    pub fn min(&self) -> Option<&SurfaceDeviation<D>> {
        self.min_index.map(|i| &self.values[i])
    }

    pub fn to_with_stats(&self) -> SurfaceDeviationSetWithStats<D> {
        todo!("Implement this method")
    }

    /// Calculates the symmetrical GD&T zone size of the set of deviations. This is the minimum
    /// sized symmetrical zone centered on the nominal value which contains all the deviations.
    ///
    /// **Note**: this is a GD&T specific calculation and is not necessarily a useful or intuitive
    /// measure outside of that context. It is not the same as the range of the deviations.
    /// Instead, it takes the larger of the absolute values of the maximum and minimum deviations
    /// and returns the double of that value.
    pub fn symmetrical_zone_size(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let v0 = self.max().unwrap().deviation.abs();
        let v1 = self.min().unwrap().deviation.abs();

        v0.max(v1) * 2.0
    }
}

/// This struct is used to represent a set of `SurfaceDeviation` items along with calculated
/// statistics about the set. It is created by transforming a `SurfaceDeviationSet` into a
/// `SurfaceDeviationSetWithStats` using the `to_with_stats` method, at which point the statistics
/// are calculated and stored in the struct along with the original data.  There is no way to
/// create an instance of this struct directly, and once created it is immutable.  This struct
/// can be converted back to a `SurfaceDeviationSet` using the `to_without_stats` method, but
/// the statistics will be lost.
pub struct SurfaceDeviationSetWithStats<const D: usize> {
    // TODO: Implement this struct
}

impl<const D: usize> SurfaceDeviationSetWithStats<D> {
    pub fn to_without_stats(&self) -> SurfaceDeviationSet<D> {
        todo!("Implement this method")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::SurfacePoint;
    use crate::{Point2, Vector2};

    fn make_dev(x: f64, y: f64, d: f64) -> SurfaceDeviation<2> {
        SurfaceDeviation::new(
            SurfacePoint::new_normalize(Point2::new(x, y), Vector2::new(0.0, 1.0)),
            d,
        )
    }

    #[test]
    fn surface_deviation_set_new_empty() {
        let sds: SurfaceDeviationSet<2> = SurfaceDeviationSet::default();
        assert!(sds.max().is_none());
        assert!(sds.min().is_none());
        assert!(sds.is_empty())
    }

    #[test]
    fn surface_deviation_set_push_single_element() {
        let mut sds = SurfaceDeviationSet::default();
        let sd = make_dev(1.0, 2.0, 0.5);
        sds.push(sd);
        assert_eq!(sds.max().unwrap().deviation, 0.5);
        assert_eq!(sds.min().unwrap().deviation, 0.5);
    }

    #[test]
    fn surface_deviation_set_push_multiple_elements() {
        let mut sds = SurfaceDeviationSet::default();
        let sd1 = make_dev(1.0, 2.0, 0.5);
        let sd2 = make_dev(3.0, 4.0, 1.5);
        sds.push(sd1);
        sds.push(sd2);
        assert_eq!(sds.max().unwrap().deviation, 1.5);
        assert_eq!(sds.min().unwrap().deviation, 0.5);
    }

    #[test]
    fn surface_deviation_set_push_multiple_updates() {
        let mut sds = SurfaceDeviationSet::default();
        sds.push(make_dev(1.0, 2.0, 0.5));
        sds.push(make_dev(3.0, 4.0, 1.5));
        sds.push(make_dev(5.0, 6.0, 1.0));
        sds.push(make_dev(7.0, 8.0, 2.0));
        sds.push(make_dev(9.0, 10.0, -0.5));
        sds.push(make_dev(11.0, 12.0, 10.0));
        assert_eq!(sds.max().unwrap().deviation, 10.0);
        assert_eq!(sds.min().unwrap().deviation, -0.5);
    }
}
