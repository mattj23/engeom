
use parry3d_f64::na::{DMatrix, Point, SVector, Unit};
use super::points::{mean_point, mean_point_weighted};

/// This structure contains the results of using singular value decomposition to determine the
/// basis vectors of a set of points and their singular values (scales). This can be used to roughly
/// estimate if a set of points in D-dimensional space falls along a point, line, or plane.
#[derive(Debug)]
pub struct SvdBasis<const D: usize> {
    /// The resultant basis vectors, sorted by their corresponding singular values so that the
    /// first vector is the most significant. These are given as unit vectors.
    pub basis: [SVector<f64, D>; D],

    /// The raw singular values associated with each basis vector. The singular values are the
    /// square root of the eigenvalues of the covariance matrix of the point set.  By squaring them
    /// and dividing by the number of points used to compute the basis, the variance accounted for
    /// by each basis vector can be determined.
    pub sv: [f64; D],

    /// The center of the original point set used to compute the basis. The SVD was computed by
    /// calculating this center (mean point) and then subtracting it from each point in the set.
    /// The basis vectors represent vectors relative to this center as their origin.
    pub center: Point<f64, D>,

    /// The number of points used to compute the basis
    pub n: usize,
}

impl<const D: usize> SvdBasis<D> {
    pub fn largest(&self) -> Unit<SVector<f64, D>> {
        Unit::new_unchecked(self.basis[0])
    }

    pub fn smallest(&self) -> Unit<SVector<f64, D>> {
        Unit::new_unchecked(self.basis[D - 1])
    }

    /// Calculates and returns the variance accounted for by each basis vector. The variance is
    /// calculated by squaring the singular value of each basis vector and dividing by the number
    /// of points used in the original decomposed matrix.
    ///
    /// The variance is the expected value of the squared deviation from the mean. It will have
    /// units of the distance squared.  For a measure of dispersion with units of distance, use
    /// the standard deviation.
    pub fn basis_variances(&self) -> [f64; D] {
        let mut result = [0.0; D];
        for (r, s) in result.iter_mut().zip(self.sv.iter()) {
            *r = s.powi(2) / (self.n as f64);
        }
        result
    }

    /// Calculates and returns the standard deviation of the point dispersion along each basis
    /// vector. These values have the same units as the original points (i.e. if the points were in
    /// millimeters, the standard deviations will be in millimeters).
    ///
    /// These are useful for interpreting the geometric results of the basis vectors. The actual
    /// dispersion of the points along each basis vector will vary depending on how the points are
    /// distributed in space, but it is reasonable to assume that almost all the points will be
    /// within a few standard deviations of the mean.  Thus, standard deviations can be used to
    /// estimate roughly how important a basis vector is in actual world units.
    pub fn basis_stdevs(&self) -> [f64; D] {
        let mut result = self.basis_variances();
        for r in result.iter_mut() {
            *r = r.sqrt();
        }
        result
    }

    /// Compute the basis vectors of a set of points using singular value decomposition.  This uses
    /// `nalgebra`'s SVD implementation. The basis vectors are sorted by their corresponding
    /// singular values so that the first vector is the most significant. The basis vectors are
    /// returned as unit vectors.
    ///
    /// The result struct can be used to compute the variance and/or standard deviation of the
    /// point set along each basis vector.  This can be used to quickly estimate if a set of points
    /// falls into a shape with a lower dimensionality than the space they are in, such as lying
    /// along a line or plane.
    ///
    /// This is a relatively fast operation on a set of points. Computation on a set of 500 points
    /// was measured to take about 25us on a Xeon E-2286M (passmark 15,120) and 70us on a Core
    /// i7-10510U (passmark 6,700).  The time appears to be linearly proportional to the number of
    /// points.
    ///
    /// # Arguments
    ///
    /// * `points`:
    /// * `weights`:
    ///
    /// returns: SvdBasis<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_points(points: &[Point<f64, D>], weights: Option<&[f64]>) -> Self {
        if let Some(w) = weights {
            let center = mean_point_weighted(points, w);
            let vectors = points
                .iter()
                .zip(w)
                .map(|(p, w)| p - center * *w)
                .collect::<Vec<_>>();
            svd_from_vectors(&vectors, Some(center))
        } else {
            let center = mean_point(points);
            let vectors = points.iter().map(|p| p - center).collect::<Vec<_>>();
            svd_from_vectors(&vectors, Some(center))
        }
    }

    /// Retrieve the rank of the decomposition by counting the number of singular values that are
    /// greater than the provided tolerance.  A rank of 0 indicates that all singular values are
    /// less than the tolerance, and thus the point set is essentially a single point. A rank of 1
    /// indicates that the point set is essentially a line. A rank of 2 indicates that the point
    /// set exists roughly in a plane.  The maximum rank is D, which indicates that the point set
    /// cannot be reduced to a lower dimension.
    ///
    /// The singular values do not directly have a clear physical meaning. They are square roots of
    /// the variance multiplied by the number of points used to compute the basis.  Thus they can
    /// be interpreted in relation to each other, and when they are very small.
    ///
    /// This method should be used either when you know roughly what a cutoff tolerance for the
    /// problem you're working on should be, or when you know the cutoff value should be very
    /// small.  Otherwise, consider examining the standard deviations of the basis vectors
    /// instead, as they will be easier to interpret (`basis_stdevs()`).
    ///
    /// # Arguments
    ///
    /// * `tol`: the largest value that a singular value can have and still be considered zero.
    ///
    /// returns: usize
    pub fn rank(&self, tol: f64) -> usize {
        let mut rank = 0;
        for s in self.sv.iter() {
            if *s > tol {
                rank += 1;
            }
        }
        rank
    }

    /// Given a point in the global coordinate system, return the coordinates of the point in the
    /// basis coordinate system.  This is done by subtracting the center of the basis from the
    /// point and then projecting the result onto each basis vector in sequence
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: OPoint<f64, Const<{ D }>>
    pub fn point_to_basis(&self, point: &Point<f64, D>) -> Point<f64, D> {
        let mut result = Point::<f64, D>::origin();
        for i in 0..D {
            let as_vec = point - self.center;
            result[i] = self.basis[i].dot(&as_vec);
        }
        result
    }

    pub fn vec_to_basis(&self, vector: &SVector<f64, D>) -> SVector<f64, D> {
        let mut result = SVector::<f64, D>::zeros();
        for i in 0..D {
            result[i] = self.basis[i].dot(vector);
        }
        result
    }

    /// Given a point in the basis coordinate system, return the coordinates of the point in the
    /// global coordinate system.  This is done by summing the product of the basis vectors and the
    /// corresponding coordinate of the point, and then adding the center of the basis.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: OPoint<f64, Const<{ D }>>
    pub fn point_from_basis(&self, point: &Point<f64, D>) -> Point<f64, D> {
        let mut result = Point::<f64, D>::origin();
        for i in 0..D {
            result += self.basis[i] * point[i];
        }
        result + self.center.coords
    }
}

fn svd_from_vectors<const D: usize>(
    vecs: &[SVector<f64, D>],
    center: Option<Point<f64, D>>,
) -> SvdBasis<D> {
    let n = vecs.len();
    let mut matrix = DMatrix::zeros(n, D);
    for (i, p) in vecs.iter().enumerate() {
        for j in 0..D {
            matrix[(i, j)] = p[j];
        }
    }

    let result = matrix.svd(false, true);
    let v_t = result.v_t.unwrap();

    let mut basis = [SVector::<f64, D>::zeros(); D];
    let mut scales = [0.0; D];
    for i in 0..D {
        for j in 0..D {
            basis[i][j] = v_t[(i, j)];
        }
        scales[i] = result.singular_values[i];
    }

    SvdBasis {
        basis,
        sv: scales,
        center: center.unwrap_or(Point::<f64, D>::origin()),
        n,
    }
}

impl From<&SvdBasis3> for Iso3 {
    fn from(value: &SvdBasis3) -> Self {
        iso3_from_basis(&value.basis, &value.center)
    }
}

impl From<&SvdBasis2> for Iso2 {
    fn from(value: &SvdBasis2) -> Self {
        iso2_from_basis(&value.basis, &value.center)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::geom3::{Point3, Vector3};

    #[test]
    fn from_points_perfect() {
        let points = vec![
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
        ];

        let result = SvdBasis3::from_points(&points, None);
        assert_relative_eq!(result.center, Point3::origin());
        assert_relative_eq!(result.basis[0], Vector3::x_axis());
        assert_relative_eq!(result.basis[1], Vector3::y_axis());
        assert_relative_eq!(result.basis[2], Vector3::z_axis());
        assert_eq!(result.n, 4);
    }


}
