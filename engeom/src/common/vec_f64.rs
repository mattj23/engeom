//! This module contains a set of common functions for working with vectors and slices of f64 values

/// Compute the arithmetic mean of a slice of f64 values.
///
/// Returns `f64::NAN` if the slice is empty.
///
/// # Arguments
///
/// * `values`: the slice of f64 values to average
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::mean_value;
/// let values = vec![1.0, 2.0, 3.0];
/// assert_eq!(mean_value(&values), 2.0);
///
/// assert!(mean_value(&[]).is_nan());
/// ```
pub fn mean_value(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN; // Return NaN if the slice is empty
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Compute the mean and population standard deviation of a slice of f64 values.
///
/// Returns `None` if the slice is empty.
///
/// # Arguments
///
/// * `values`: the slice of f64 values to summarize
///
/// returns: `Option<(f64, f64)>` — `Some((mean, stdev))` or `None`
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::mean_and_stdev;
/// use approx::assert_relative_eq;
///
/// let (mean, stdev) = mean_and_stdev(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap();
/// assert_relative_eq!(mean, 5.0, epsilon = 1e-12);
/// assert_relative_eq!(stdev, 2.0, epsilon = 1e-12);
///
/// assert!(mean_and_stdev(&[]).is_none());
/// ```
pub fn mean_and_stdev(values: &[f64]) -> Option<(f64, f64)> {
    if values.is_empty() {
        return None;
    }
    let mean = mean_value(values);
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
    Some((mean, variance.sqrt()))
}

/// Checks if a slice contains any NaN values
///
/// # Arguments
///
/// * `values`: the slice of f64 values to check for NaN
///
/// returns: bool
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::has_nan;
/// let values = vec![1.0, f64::NAN, 3.0];
/// assert!(has_nan(&values));
///
/// let values = vec![1.0, 2.0, 3.0];
/// assert!(!has_nan(&values));
/// ```
pub fn has_nan(values: &[f64]) -> bool {
    values.iter().any(|v| v.is_nan())
}

/// Checks if all values in the slice are finite
///
/// # Arguments
///
/// * `values`: a slice of f64 values to test for finiteness
///
/// returns: bool
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::are_all_finite;
/// assert!(are_all_finite(&[1.0, 2.0, 3.0]));
/// assert!(!are_all_finite(&[1.0, f64::INFINITY, 3.0]));
/// assert!(!are_all_finite(&[1.0, f64::NAN, 3.0]));
/// ```
pub fn are_all_finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}

/// Checks if all values in the slice are sorted in ascending order
///
/// # Arguments
///
/// * `values`: a slice of f64 values to test for ascending order
///
/// returns: bool
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::are_in_ascending_order;
/// let values = vec![1.0, 2.0, 3.0];
/// assert!(are_in_ascending_order(&values));
///
/// let values = vec![1.0, 3.0, 2.0];
/// assert!(!are_in_ascending_order(&values));
/// ```
pub fn are_in_ascending_order(values: &[f64]) -> bool {
    values.windows(2).all(|w| w[0] <= w[1])
}

/// Checks if all values in the slice are sorted in descending order
///
/// # Arguments
///
/// * `values`: a slice of f64 values to test for descending order
///
/// returns: bool
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::are_in_descending_order;
/// let values = vec![3.0, 2.0, 1.0];
/// assert!(are_in_descending_order(&values));
///
/// let values = vec![1.0, 3.0, 2.0];
/// assert!(!are_in_descending_order(&values));
/// ```
pub fn are_in_descending_order(values: &[f64]) -> bool {
    values.windows(2).all(|w| w[0] >= w[1])
}

/// Sorts the slice in ascending order, NaN values are allowed and are sorted to the end of the
/// slice.
///
/// # Arguments
///
/// * `values`: a mutable slice of f64 values to sort
///
/// returns: ()
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::sort_with_nan;
/// let mut values = vec![3.0, f64::NAN, 1.0];
/// sort_with_nan(&mut values);
/// assert_eq!(values[0], 1.0);
/// assert_eq!(values[1], 3.0);
/// assert!(values[2].is_nan());
/// ```
pub fn sort_with_nan(values: &mut [f64]) {
    values.sort_by(|a, b| {
        if a.is_nan() {
            std::cmp::Ordering::Greater
        } else if b.is_nan() {
            std::cmp::Ordering::Less
        } else {
            a.partial_cmp(b).unwrap()
        }
    });
}

/// Sorts the slice in ascending order, NaN values are not allowed and will cause a panic. If you
/// aren't absolutely sure that the slice does not contain any NaN values, use the
/// `has_nan` function to check for NaN values before calling this function, or use the
/// `sort_with_nan` function in this module instead.
///
/// # Arguments
///
/// * `values`: a mutable slice of f64 values to sort, must not contain any NaN values
///
/// returns: ()
///
/// # Examples
///
/// ```
/// use engeom::common::vec_f64::sort_nan_panics;
/// let mut values = vec![3.0, 5.0, 1.0];
/// sort_nan_panics(&mut values);
/// assert_eq!(values, vec![1.0, 3.0, 5.0]);
/// ```
pub fn sort_nan_panics(values: &mut [f64]) {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn mean_value_empty_returns_nan() {
        assert!(mean_value(&[]).is_nan());
    }

    #[test]
    fn mean_value_single_element() {
        assert_relative_eq!(mean_value(&[7.0]), 7.0, epsilon = 1e-12);
    }

    #[test]
    fn mean_value_uniform_values() {
        assert_relative_eq!(mean_value(&[3.0, 3.0, 3.0]), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn mean_value_typical() {
        assert_relative_eq!(mean_value(&[1.0, 2.0, 3.0, 4.0]), 2.5, epsilon = 1e-12);
    }

    #[test]
    fn mean_and_stdev_empty_returns_none() {
        assert!(mean_and_stdev(&[]).is_none());
    }

    #[test]
    fn mean_and_stdev_single_element_has_zero_stdev() {
        let (mean, stdev) = mean_and_stdev(&[5.0]).unwrap();
        assert_relative_eq!(mean, 5.0, epsilon = 1e-12);
        assert_relative_eq!(stdev, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn mean_and_stdev_uniform_values_have_zero_stdev() {
        let (mean, stdev) = mean_and_stdev(&[4.0, 4.0, 4.0]).unwrap();
        assert_relative_eq!(mean, 4.0, epsilon = 1e-12);
        assert_relative_eq!(stdev, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn mean_and_stdev_known_values() {
        // Population stdev of [2, 4, 4, 4, 5, 5, 7, 9] = 2.0, mean = 5.0
        let (mean, stdev) = mean_and_stdev(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap();
        assert_relative_eq!(mean, 5.0, epsilon = 1e-12);
        assert_relative_eq!(stdev, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn has_nan_test() {
        let values = vec![1.0, f64::NAN, 3.0];
        assert!(has_nan(&values));

        let values = vec![1.0, 2.0, 3.0];
        assert!(!has_nan(&values));
    }

    #[test]
    fn are_all_finite_test() {
        let values = vec![1.0, 2.0, 3.0];
        assert!(are_all_finite(&values));

        let values = vec![1.0, f64::INFINITY, 3.0];
        assert!(!are_all_finite(&values));
    }

    #[test]
    fn are_in_ascending_order_test() {
        let values = vec![1.0, 2.0, 3.0];
        assert!(are_in_ascending_order(&values));

        let values = vec![1.0, 3.0, 2.0];
        assert!(!are_in_ascending_order(&values));
    }

    #[test]
    fn are_in_descending_order_test() {
        let values = vec![3.0, 2.0, 1.0];
        assert!(are_in_descending_order(&values));

        let values = vec![1.0, 3.0, 2.0];
        assert!(!are_in_descending_order(&values));
    }

    #[test]
    fn sort_with_nan_test() {
        let mut values = vec![3.0, f64::NAN, 1.0];
        sort_with_nan(&mut values);
        assert_eq!(values[0], 1.0);
        assert_eq!(values[1], 3.0);
        assert!(values[2].is_nan());
    }

    #[test]
    fn sort_nan_panics_test() {
        let mut values = vec![3.0, 5.0, 1.0];
        sort_nan_panics(&mut values);
        assert_eq!(values, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    #[should_panic]
    fn sort_nan_panics_with_nan_test() {
        let mut values = vec![3.0, f64::NAN, 1.0];
        sort_nan_panics(&mut values);
    }
}
