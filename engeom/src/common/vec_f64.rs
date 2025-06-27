//! This module contains a set of common functions for working with vectors and slices of f64 values

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
///
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
