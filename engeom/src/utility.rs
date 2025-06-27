//! Utility functions to help with serialization, transformation, and other tasks
//!
//!

use crate::errors::FailedConversion;
use crate::Result;
use parry3d_f64::na::{Point, SVector};

/// Converts a slice of arrays of floating point values into a vector of points with the specified
/// dimensionality.
///
/// # Arguments
///
/// * `slice`: the slice of [f64; D] arrays of floating point values to be converted
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
/// use engeom::Point2;
/// use engeom::utility::slice_to_points;
/// let slice = vec![[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]];
/// let points = slice_to_points(&slice);
/// assert_eq!(points.len(), 3);
/// assert_eq!(points[0], Point2::new(0.0, 0.5));
/// assert_eq!(points[1], Point2::new(1.0, 1.5));
/// assert_eq!(points[2], Point2::new(2.0, 2.5));
/// ```
pub fn slice_to_points<const D: usize>(slice: &[[f64; D]]) -> Vec<Point<f64, D>> {
    slice.iter().map(|p| Point::from(*p)).collect()
}

/// Converts a slice of arrays of floating point values into a vector of vectors with the specified
/// dimensionality.
///
/// # Arguments
///
/// * `slice`: the slice of [f64; D] arrays of floating point values to be converted
///
/// returns: Vec<Matrix<f64, Const<{ D }>, Const<1>, ArrayStorage<f64, { D }, 1>>, Global>
///
/// # Examples
///
/// ```
/// use engeom::Vector2;
/// use engeom::utility::slice_to_vectors;
/// let slice = vec![[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]];
/// let vectors = slice_to_vectors(&slice);
/// assert_eq!(vectors.len(), 3);
/// assert_eq!(vectors[0], Vector2::new(0.0, 0.5));
/// assert_eq!(vectors[1], Vector2::new(1.0, 1.5));
/// assert_eq!(vectors[2], Vector2::new(2.0, 2.5));
/// ```
pub fn slice_to_vectors<const D: usize>(slice: &[[f64; D]]) -> Vec<SVector<f64, D>> {
    slice.iter().map(|p| SVector::from(*p)).collect()
}

/// Takes a slice of floating point values and returns a vector of points with the specified
/// dimensionality.  The number of values must be a multiple of the dimensionality, otherwise an
/// `InvalidCount` error will be returned.
///
/// # Arguments
///
/// * `values`: the slice of floating point values to be unflattened, must be a multiple of `D`
///
/// returns: Result<Vec<OPoint<f64, Const<{ D }>>, Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
/// use engeom::Point2;
/// use engeom::utility::unflatten_points;
/// let values = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5];
/// let points = unflatten_points::<2>(&values).unwrap();
///
/// assert_eq!(points.len(), 3);
/// assert_eq!(points[0], Point2::new(0.0, 0.5));
/// assert_eq!(points[1], Point2::new(1.0, 1.5));
/// assert_eq!(points[2], Point2::new(2.0, 2.5));
/// ```
pub fn unflatten_points<const D: usize>(values: &[f64]) -> Result<Vec<Point<f64, D>>> {
    if values.len() % D != 0 {
        return Err(Box::new(FailedConversion::InvalidCount));
    }

    let mut points = Vec::new();
    for i in 0..values.len() / D {
        let mut coords = [0.0; D];
        for j in 0..D {
            coords[j] = values[i * D + j];
        }
        points.push(Point::from(coords));
    }

    Ok(points)
}

/// Takes a slice of floating point values and returns a `Vec` of vectors with the specified
/// dimensionality.  The number of values must be a multiple of the dimensionality, otherwise an
/// `InvalidCount` error will be returned.
///
/// # Arguments
///
/// * `values`: the slice of floating point values to be unflattened, must be a multiple of `D`
///
/// returns: Result<Vec<Matrix<f64, Const<{ D }>, Const<1>, ArrayStorage<f64, { D }, 1>>, Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
/// use engeom::Vector2;
/// use engeom::utility::unflatten_vectors;
/// let values = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5];
/// let points = unflatten_vectors::<2>(&values).unwrap();
///
/// assert_eq!(points.len(), 3);
/// assert_eq!(points[0], Vector2::new(0.0, 0.5));
/// assert_eq!(points[1], Vector2::new(1.0, 1.5));
/// assert_eq!(points[2], Vector2::new(2.0, 2.5));
/// ```
pub fn unflatten_vectors<const D: usize>(values: &[f64]) -> Result<Vec<SVector<f64, D>>> {
    if values.len() % D != 0 {
        return Err(Box::new(FailedConversion::InvalidCount));
    }

    let mut vectors = Vec::new();
    for i in 0..values.len() / D {
        let mut coords = [0.0; D];
        for j in 0..D {
            coords[j] = values[i * D + j];
        }
        vectors.push(SVector::from(coords));
    }

    Ok(vectors)
}

/// Flattens a slice of points into a vector of floating point values.
///
/// # Arguments
///
/// * `points`: the slice of points to be flattened
///
/// returns: Vec<f64, Global>
///
/// # Examples
///
/// ```
/// use engeom::Point2;
/// use engeom::utility::flatten_points;
/// let points = vec![Point2::new(0.0, 0.5), Point2::new(1.0, 1.5), Point2::new(2.0, 2.5)];
/// let values = flatten_points(&points);
///
/// assert_eq!(values, vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5]);
/// ```
pub fn flatten_points<const D: usize>(points: &[Point<f64, D>]) -> Vec<f64> {
    let mut values = Vec::new();
    for p in points.iter() {
        for i in 0..D {
            values.push(p[i]);
        }
    }
    values
}

/// Flattens a slice of vectors into a vec of floating point values.
///
/// # Arguments
///
/// * `points`: the slice of vectors to be flattened
///
/// returns: Vec<f64, Global>
///
/// # Examples
///
/// ```
/// use engeom::Vector2;
/// use engeom::utility::flatten_vectors;
/// let vectors = vec![Vector2::new(0.0, 0.5), Vector2::new(1.0, 1.5), Vector2::new(2.0, 2.5)];
/// let values = flatten_vectors(&vectors);
///
/// assert_eq!(values, vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5]);
/// ```
pub fn flatten_vectors<const D: usize>(vectors: &[SVector<f64, D>]) -> Vec<f64> {
    let mut values = Vec::new();
    for p in vectors.iter() {
        for i in 0..D {
            values.push(p[i]);
        }
    }
    values
}
