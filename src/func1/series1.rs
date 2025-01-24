use crate::common::{DiscreteDomain, Interval};
use crate::func1::{Func1, Line1};
use crate::{Point2, Result};
use serde::{Deserialize, Serialize};
use std::iter::Zip;
use std::ops;
use std::slice::Iter;

/// Represents a contiguous series of data points in a 2d plane where the x values go from
/// smallest to largest and the y values are associated with the x value at the same index.
/// TODO: Assert that x is always sorted
/// TODO: Assert that x and y are the same length
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Series1 {
    pub x: DiscreteDomain,
    pub y: Vec<f64>,
}

impl Series1 {
    /// Creates a new series from a pair of vectors, one for x values and one for y values. The two
    /// vectors must be the same length, and the x values must be sorted from smallest to largest.
    pub fn new(x: DiscreteDomain, y: Vec<f64>) -> Self {
        Self { x, y }
    }

    /// Try to create a new series from a pair of vectors, one for x values and one for y values.
    /// If  the two vectors are not the same length, an error will be returned. If the x values are
    /// not sorted from smallest to largest, an error will be returned.
    ///
    /// # Arguments
    ///
    /// * `x`: a vector of x values, matched in order to the y values
    /// * `y`: a vector of y values, matched in order to the x values
    ///
    /// returns: Result<Series1, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn try_new(x: Vec<f64>, y: Vec<f64>) -> Result<Self> {
        if x.len() != y.len() {
            Err("x and y must be the same length".into())
        } else {
            Ok(Self {
                x: DiscreteDomain::try_from(x)?,
                y,
            })
        }
    }

    /// Create a new series from a `Func1` by sampling it at a series of domain values and
    /// combining the result into a series.
    ///
    /// # Arguments
    ///
    /// * `f`:
    /// * `xs`:
    ///
    /// returns: Series1
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_sampled(f: &dyn Func1, xs: DiscreteDomain) -> Self {
        let ys = f.fs(&xs);
        Self::new(xs, ys)
    }

    /// Returns the smallest x value in the series, which also happens to be the first x value
    pub fn x_min(&self) -> f64 {
        self.x[0]
    }

    /// Finds and returns the largest x value in the series, which also happens to be last x value
    pub fn x_max(&self) -> f64 {
        self.x[self.x.len() - 1]
    }

    /// Finds and returns the smallest y value in the series
    pub fn y_min(&self) -> f64 {
        *self
            .y
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
    }

    /// Finds and returns the largest y value in the series
    pub fn y_max(&self) -> f64 {
        *self
            .y
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
    }

    /// Returns a new series where the x values are the same as the original series, but the y
    /// values are the absolute value of the original y values
    pub fn abs(&self) -> Self {
        let ys = self.y.iter().map(|v| v.abs()).collect::<Vec<f64>>();
        Self::new(self.x.clone(), ys)
    }

    /// Return a new series with the x values scaled by `scale_x` and the y values scaled by
    /// `scale_y`
    pub fn scaled_by(&self, scale_x: f64, scale_y: f64) -> Self {
        let xs = self
            .x
            .values()
            .iter()
            .map(|v| v * scale_x)
            .collect::<Vec<f64>>();
        let ys = self.y.iter().map(|v| v * scale_y).collect::<Vec<f64>>();

        if scale_x < 0.0 {
            let new_xs =
                DiscreteDomain::try_from(xs.into_iter().rev().collect::<Vec<_>>()).unwrap();
            Self::new(new_xs, ys.into_iter().rev().collect())
        } else {
            let new_xs = DiscreteDomain::try_from(xs).unwrap();
            Self::new(new_xs, ys)
        }
    }

    /// Returns a new series with the y values scaled by the function `fy` evaluated at the
    /// corresponding x value
    pub fn scaled_y(&self, fy: &dyn Func1) -> Self {
        let scale = fy.fs(&self.x);
        let ys = self
            .y
            .iter()
            .zip(scale.iter())
            .map(|(y, s)| y * s)
            .collect::<Vec<f64>>();
        Self::new(self.x.clone(), ys)
    }

    /// Return a new series with the x values shifted by `shift_x` and the y values shifted by
    /// `shift_y`
    pub fn shift_by(&self, shift_x: f64, shift_y: f64) -> Self {
        let xs = self.x.iter().map(|v| v + shift_x).collect::<Vec<f64>>();
        let ys = self.y.iter().map(|v| v + shift_y).collect::<Vec<f64>>();
        let new_xs = DiscreteDomain::try_from(xs).unwrap();
        Self::new(new_xs, ys)
    }

    pub fn remove_nan(&self) -> Self {
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for (x, y) in self.x.iter().zip(self.y.iter()) {
            if x.is_nan() || y.is_nan() {
                continue;
            }
            xs.push(*x);
            ys.push(*y);
        }
        let new_xs = DiscreteDomain::try_from(xs).unwrap();
        Self::new(new_xs, ys)
    }

    pub fn has_nan_between(&self, x0: f64, x1: f64) -> bool {
        for (x, y) in self.x.iter().zip(self.y.iter()) {
            if y.is_nan() && *x >= x0 && *x <= x1 {
                return true;
            }
        }
        false
    }

    pub fn has_nan(&self) -> bool {
        for y in self.y.iter() {
            if y.is_nan() {
                return true;
            }
        }
        false
    }

    /// Calculate and return the first derivative of the series
    pub fn dydx(&self) -> Self {
        let mut ys = Vec::new();
        for j in 0..self.y.len() {
            if j == 0 {
                ys.push((self.y[1] - self.y[0]) / (self.x[1] - self.x[0]));
            } else if j == self.y.len() - 1 {
                ys.push((self.y[j] - self.y[j - 1]) / (self.x[j] - self.x[j - 1]));
            } else {
                ys.push((self.y[j + 1] - self.y[j - 1]) / (self.x[j + 1] - self.x[j - 1]));
            }
        }

        Self::new(self.x.clone(), ys)
    }

    /// Perform a Savitzky-Golay smoothing on the series with a 5-point window
    pub fn savitzky_golay(&self) -> Self {
        savitzky_golay_5(self)
    }

    /// Return the x values of the places where the series crosses the y value `y_equals`. The
    /// values will be sorted from smallest to largest and will be unique.
    ///
    /// # Arguments
    ///
    /// * `y_equals`: the y value at which to find the corresponding x values
    ///
    /// returns: Vec<f64, Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn y_crossings(&self, y_equals: f64) -> Vec<f64> {
        let mut crossings = Vec::new();
        for j in 0..self.y.len() - 1 {
            let v0 = self.y[j];
            let v1 = self.y[j + 1];
            if v0 <= y_equals && v1 >= y_equals || v0 >= y_equals && v1 <= y_equals {
                let x0 = self.x[j];
                let x1 = self.x[j + 1];
                let m = (v1 - v0) / (x1 - x0);
                if !m.is_finite() {
                    continue;
                }
                let x = x0 + (y_equals - v0) / m;
                crossings.push(x);
            }
        }

        sort_and_dedup(&mut crossings);
        crossings
    }

    /// Return the index of the x value in the series which is either equal to or immediately
    /// after the given x value.  If the given x value is before the first x value in the series,
    /// 0 will be returned.  If the given x value is after the last x value in the series, the
    /// length of the series will be returned.
    pub fn index_of_x_after(&self, x: f64) -> usize {
        let search_result = self.x.binary_search_by(|v| v.partial_cmp(&x).unwrap());
        match search_result {
            Ok(actual_index) => actual_index,
            Err(index_after) => index_after,
        }
    }

    /// Return the interpolated y value at x. If the x value is outside the range of the series,
    /// NAN will be returned.
    pub fn interpolate(&self, x: f64) -> f64 {
        if x < self.x[0] || x > self.x[self.x.len() - 1] {
            f64::NAN
        } else {
            let search_result = self.x.binary_search_by(|v| v.partial_cmp(&x).unwrap());
            match search_result {
                Ok(actual_index) => self.y[actual_index],
                Err(index_after) => {
                    if index_after == 0 {
                        panic!("This should never happen")
                    } else {
                        let x0 = self.x[index_after - 1];
                        let x1 = self.x[index_after];
                        let y0 = self.y[index_after - 1];
                        let y1 = self.y[index_after];
                        let m = (y1 - y0) / (x1 - x0);
                        y0 + m * (x - x0)
                    }
                }
            }
        }
    }

    /// Return the y value at index i. If i is outside the range of the series, None will be
    /// returned instead.
    fn y_at_index(&self, i: i32) -> Option<f64> {
        if i < 0 || i >= self.y.len() as i32 {
            None
        } else {
            Some(self.y[i as usize])
        }
    }

    pub fn is_ordered(&self) -> bool {
        for i in 1..self.x.len() {
            if self.x[i] <= self.x[i - 1] {
                return false;
            }
        }
        true
    }

    /// Returns a new series that contains the region of the series between x0 and x1, including
    /// the end points. This breaks the regularity of the spacing
    pub fn between(&self, x0: f64, x1: f64) -> Self {
        let mut xs = Vec::new();
        let mut ys = Vec::new();

        // Find the first index that is >= x0
        let search_result = self.x.binary_search_by(|a| a.partial_cmp(&x0).unwrap());
        let mut i = match search_result {
            Ok(actual_index) => actual_index,
            Err(next_index) => {
                if next_index == 0 {
                    0
                } else {
                    xs.push(x0);
                    ys.push(self.interpolate(x0));
                    next_index
                }
            }
        };

        while i < self.x.len() && self.x[i] <= x1 {
            xs.push(self.x[i]);
            ys.push(self.y[i]);
            i += 1;
        }

        if xs[xs.len() - 1] < x1 {
            xs.push(x1);
            ys.push(self.interpolate(x1));
        }

        let new_xs = DiscreteDomain::try_from(xs).unwrap();
        Self::new(new_xs, ys)
    }

    pub fn interval(&self) -> Interval {
        Interval::new(self.x_min(), self.x_max())
    }

    pub fn in_interval(&self, interval: Interval) -> Self {
        self.between(interval.min, interval.max)
    }

    /// Calculates a set of intervals which represent the regions of the series that are either
    /// above the y=0 line or below the y=0 line.
    ///
    /// returns: a vector of `Interval` structs, which will be ordered from the lowest x value to
    /// the highest x value. The intervals will be non-overlapping and will cover the entire
    /// domain of the series.
    pub fn bounds_at_y0(&self) -> Vec<Interval> {
        let mut x_bounds = [self.y_crossings(0.0), vec![self.x_min(), self.x_max()]].concat();
        sort_and_dedup(&mut x_bounds);

        let mut bounds = Vec::new();
        for w in x_bounds.windows(2) {
            bounds.push(Interval::new(w[0], w[1]));
        }

        bounds
    }

    pub fn middle_reiemann_areas(&self) -> Vec<(f64, f64)> {
        let mut result = Vec::new();
        for i in 0..self.x.len() - 1 {
            let x0 = self.x[i];
            let x1 = self.x[i + 1];
            let y0 = self.y[i];
            let y1 = self.y[i + 1];
            result.push(((x1 + x0) * 0.5, (x1 - x0) * (y0 + y1) * 0.5));
        }
        result
    }

    /// Computes the area under the curve
    pub fn area_under(&self) -> f64 {
        self.middle_reiemann_areas()
            .into_iter()
            .map(|(_, a)| a)
            .sum()
    }

    /// Computes the least-squares line of best fit returning the slope and y-intercept (m and b)
    /// in that order as a tuple.
    pub fn best_fit_line(&self) -> Line1 {
        // TODO: Replace this with the polynomial least squares passthrough
        let n = self.x.len() as f64;
        let sum_x = self.x.iter().sum::<f64>();
        let sum_y = self.y.iter().sum::<f64>();
        let sum_xx = self.x.iter().map(|x| x * x).sum::<f64>();
        let sum_xy = self
            .x
            .iter()
            .zip(self.y.iter())
            .map(|(x, y)| x * y)
            .sum::<f64>();
        let m = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
        let b = (sum_y - m * sum_x) / n;
        Line1::new_mxb(m, b)
    }

    /// Calculates and returns a vector of the x (domain) location of all the local maxima in
    /// the series. A local maximum is defined as a point whose neighbors are both less than it.
    pub fn local_maxima_xs(&self) -> Vec<f64> {
        let mut maxima = Vec::new();
        for i in 0..self.x.len() {
            if i == 0 {
                if self.y[0] > self.y[1] {
                    maxima.push(self.x[0]);
                }
            } else if i == self.y.len() - 1 {
                if self.y[i] > self.y[i - 1] {
                    maxima.push(self.x[i]);
                }
            } else {
                let y0 = self.y[i - 1];
                let y1 = self.y[i];
                let y2 = self.y[i + 1];
                if y1 > y0 && y1 > y2 {
                    maxima.push(self.x[i]);
                }
            }
        }
        maxima
    }

    /// Calculates and returns the x and y location of the global maxima in the series. This is
    /// the x position of the y with the highest value. If there are multiple y values with the
    /// same maximum value, the one with the lowest x value is returned.
    pub fn global_maxima_xy(&self) -> (f64, f64) {
        let mut max_x = self.x[0];
        let mut max_y = self.y[0];
        for i in 1..self.x.len() {
            if self.y[i] > max_y {
                max_x = self.x[i];
                max_y = self.y[i];
            }
        }
        (max_x, max_y)
    }

    /// Calculates and returns the x and y location of the global minima in the series. This is
    /// the x position of the y with the lowest value. If there are multiple y values with the
    /// same minimum value, the one with the lowest x value is returned.
    pub fn global_minima_xy(&self) -> (f64, f64) {
        let mut min_x = self.x[0];
        let mut min_y = self.y[0];
        for i in 1..self.x.len() {
            if self.y[i] < min_y {
                min_x = self.x[i];
                min_y = self.y[i];
            }
        }
        (min_x, min_y)
    }

    /// Calculates and returns the x and y location of the global minima in the series. This is the
    /// x position of the y with the lowest value. If there are multiple y values with the same
    /// minimum value, the one with the lowest x value is returned.
    pub fn global_minima_x(&self) -> (f64, f64) {
        let mut min_x = self.x[0];
        let mut min_y = self.y[0];
        for i in 1..self.x.len() {
            if self.y[i] < min_y {
                min_x = self.x[i];
                min_y = self.y[i];
            }
        }
        (min_x, min_y)
    }

    /// Attempts to split the series into two series at the given x value, with the first value
    /// in the return tuple consisting of the series below the x value and the second value in the
    /// tuple consisting of the series above the x value. If the x value is not in the domain of
    /// the series, then the corresponding value in the tuple will be None.
    pub fn split_at_x(&self, x: f64) -> (Option<Self>, Option<Self>) {
        if x > self.x_max() {
            (Some(self.clone()), None)
        } else if x < self.x_min() {
            (None, Some(self.clone()))
        } else {
            (
                Some(self.between(self.x_min(), x)),
                Some(self.between(x, self.x_max())),
            )
        }
    }

    /// Calculates the lower and upper x bounds of a plateau based on an x value assumed to be
    /// of a local maxima.  This works by finding the value of the maxima, subtracting the
    /// tolerance, and then finding the first x value to the left and right of the maxima that has
    /// the y value of the maxima - tolerance.
    pub fn plateau_at_maxima(&self, x: f64, tol: f64) -> Option<Interval> {
        let v = self.interpolate(x);
        if v.is_nan() {
            return None;
        }
        let y = v - tol;
        let mut crossings = self.y_crossings(y);

        // If the first y value is above the threshold or the last y value is below the threshold,
        // then those points are crossings
        if self.y[0] > y {
            crossings.insert(0, self.x[0]);
        }

        if self.y[self.y.len() - 1] > y {
            crossings.push(self.x[self.x.len() - 1]);
        }

        // TODO: This is a linear search, but it could be a binary search
        for xs in crossings.windows(2) {
            if xs[0] <= x && x <= xs[1] {
                return Some(Interval::new(xs[0], xs[1]));
            }
        }
        None
    }

    /// Returns a zipped iterator over the x and y values of the series
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Series1;
    /// let series = Series1::try_new(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0])
    ///     .expect("Failed to create series");
    /// for (x, y) in series.xys() {
    ///    println!("x: {}, y: {}", x, y);
    /// }
    /// ```
    pub fn xys(&self) -> Zip<Iter<'_, f64>, Iter<'_, f64>> {
        self.x.values().iter().zip(self.y.iter())
    }

    /// Returns the series as a vector of `Point2` structs by zipping the x and y values together
    /// and converting them.  The order of the points will be from the smallest x value to the
    /// largest x value.
    pub fn as_points(&self) -> Vec<Point2> {
        self.xys()
            .map(|(x, y)| Point2::new(*x, *y))
            .collect::<Vec<_>>()
    }

    /// Resamples the series to have `n` evenly spaced points in the domain. The first and last
    /// points will be the same as the first and last points of the original series.
    ///
    /// # Arguments
    ///
    /// * `n`:
    ///
    /// returns: Series1
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Series1;
    /// let series = Series1::try_new(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0])
    ///     .expect("Failed to create series");
    /// let resampled = series.resampled_n(5);
    /// assert_eq!(resampled.x.values(), vec![0.0, 0.5, 1.0, 1.5, 2.0]);
    /// assert_eq!(resampled.y, vec![0.0, 0.5, 1.0, 1.5, 2.0]);
    /// ```
    pub fn resampled_n(&self, n: usize) -> Self {
        let step_size = (self.x_max() - self.x_min()) / (n as f64 - 1.0);
        // TODO: Make this linear rather than using the binary search
        let xs = (0..n)
            .map(|i| (self.x_min() + (i as f64) * step_size).min(self.x_max()))
            .collect::<Vec<_>>();
        let new_xs = DiscreteDomain::try_from(xs).unwrap();
        let ys = self.fs(&new_xs);

        Self::new(new_xs, ys)
    }

    /// Resamples the series to have roughly `x_spacing` between each point in the domain. The
    /// first and last points will be the same as the first and last points of the original series,
    /// and a number of evenly spaced points will be added in between.  The total number of points
    /// will be the length of the curve divided by `x_spacing` plus 1, rounded up to the nearest
    /// whole integer.
    ///
    /// # Arguments
    ///
    /// * `x_spacing`: the approximate spacing between points in the domain
    ///
    /// returns: Series1
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Series1;
    /// let series = Series1::try_new(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0])
    ///     .expect("Failed to create series");
    /// let resampled = series.resampled_x(0.5);
    /// assert_eq!(resampled.x.values(), vec![0.0, 0.5, 1.0, 1.5, 2.0]);
    /// assert_eq!(resampled.y, vec![0.0, 0.5, 1.0, 1.5, 2.0]);
    /// ```
    pub fn resampled_x(&self, x_spacing: f64) -> Self {
        let n = 1.0 + (self.x_max() - self.x_min()) / x_spacing;
        self.resampled_n(n.ceil() as usize)
    }

    /// Calculates and returns clusters of contiguous indices where the y value is within `d_tol`
    /// of the reference function and the x span of the cluster is at least `x_span`.  This is a
    /// specialized method used for identifying regions of a series that match to another series,
    /// but directly using the indices of the series rather than the x values.
    ///
    /// # Arguments
    ///
    /// * `reference`: a reference `Series1` to compare against
    /// * `d_tol`: the maximum allowable y difference between the reference and the series at any
    ///  point
    /// * `x_span`: the minimum allowable x span of any cluster
    ///
    /// returns: Vec<Vec<usize, Global>, Global>
    pub fn index_clusters_in_tol(
        &self,
        reference: &Series1,
        d_tol: f64,
        x_span: f64,
    ) -> Vec<Vec<usize>> {
        // First pass just collects everything within the tolerance
        let mut indices = Vec::new();
        for (i, (x, y)) in self.xys().enumerate() {
            if (reference.f(*x) - *y).abs() < d_tol {
                indices.push(i);
            }
        }

        let mut clusters = Vec::new();
        let mut current = Vec::new();
        for i in indices {
            if current.is_empty() {
                current.push(i);
            } else {
                let last = current.last().unwrap();
                if i - last == 1 {
                    current.push(i);
                } else {
                    clusters.push(current);
                    current = vec![i];
                }
            }
        }
        if !current.is_empty() {
            clusters.push(current);
        }

        // Now go through each cluster and if it is large enough, add all of the indices to the
        // result vector
        let mut results = Vec::new();
        for cluster in clusters {
            if self.x[*cluster.last().unwrap()] - self.x[cluster[0]] >= x_span {
                results.push(cluster);
            }
        }

        results
    }
}

impl Func1 for Series1 {
    fn f(&self, x: f64) -> f64 {
        self.interpolate(x)
    }
}

impl ops::Add<&dyn Func1> for &Series1 {
    type Output = Series1;

    fn add(self, rhs: &dyn Func1) -> Self::Output {
        let mut ys = Vec::new();
        for (x, y) in self.x.iter().zip(self.y.iter()) {
            ys.push(y + rhs.f(*x));
        }
        Self::Output::new(self.x.clone(), ys)
    }
}

impl ops::Sub<&dyn Func1> for &Series1 {
    type Output = Series1;

    fn sub(self, rhs: &dyn Func1) -> Self::Output {
        let mut ys = Vec::new();
        for (x, y) in self.x.iter().zip(self.y.iter()) {
            ys.push(y - rhs.f(*x));
        }
        Self::Output::new(self.x.clone(), ys)
    }
}

fn sort_and_dedup(xs: &mut Vec<f64>) {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
}

fn savitzky_golay_5(series: &Series1) -> Series1 {
    let mut ys = Vec::new();
    for j in 0..series.y.len() {
        let mut sum = 0.0;
        let mut total = 0.0;

        if let Some(v) = series.y_at_index(j as i32 - 2) {
            sum += v * -3.0;
            total += -3.0;
        }
        if let Some(v) = series.y_at_index(j as i32 - 1) {
            sum += v * 12.0;
            total += 12.0;
        }
        sum += series.y[j] * 17.0;
        total += 17.0;
        if let Some(v) = series.y_at_index(j as i32 + 1) {
            sum += v * 12.0;
            total += 12.0;
        }
        if let Some(v) = series.y_at_index(j as i32 + 2) {
            sum += v * -3.0;
            total += -3.0;
        }
        ys.push(sum / total);
    }

    Series1::new(series.x.clone(), ys)
}
