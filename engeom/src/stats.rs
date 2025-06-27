//! This module contains simple statistical tools

use crate::Result;

/// This function computes the mean of a slice of f64 values.
///
/// # Arguments
///
/// * `values`:
///
/// returns: Result<f64, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn compute_mean(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        Err("Cannot compute mean of empty slice".into())
    } else {
        Ok(values.iter().sum::<f64>() / values.len() as f64)
    }
}

/// This function computes the variance of a slice of f64 values.
///
/// # Arguments
///
/// * `values`:
///
/// returns: Result<f64, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn compute_variance(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        Err("Cannot compute variance of empty slice".into())
    } else {
        let mean = compute_mean(values)?;
        let mut sum = 0.0;
        for v in values {
            sum += (v - mean).powi(2);
        }
        Ok(sum / values.len() as f64)
    }
}

/// This function computes the standard deviation of a slice of f64 values.
///
/// # Arguments
///
/// * `values`:
///
/// returns: Result<f64, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn compute_st_dev(values: &[f64]) -> Result<f64> {
    compute_variance(values).map(|v| v.sqrt())
}

/// This function computes the median of a slice of f64 values.
/// The median is the value separating the higher half from the lower half of a data sample.
/// For a data set, it may be thought of as the "middle" value.
///
/// # Arguments
///
/// * `values`:
///
/// returns: Result<f64, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn compute_median(values: &[f64]) -> Result<f64> {
    if values.is_empty() {
        Err("Cannot compute median of empty slice".into())
    } else {
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            Ok((sorted[mid] + sorted[mid - 1]) / 2.0)
        } else {
            Ok(sorted[mid])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_test() {
        let values = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(compute_mean(&values).unwrap(), 2.5);
    }

    #[test]
    fn variance_test() {
        let values = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(compute_variance(&values).unwrap(), 1.25);
    }

    #[test]
    fn std_dev_test() {
        let values = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(compute_st_dev(&values).unwrap(), 1.118033988749895);
    }
}
