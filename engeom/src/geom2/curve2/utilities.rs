use crate::{Curve2, Result};
use crate::common::dist;

impl Curve2 {

    /// This utility function takes ownership of a vector of `Curve2` instances and recursively
    /// merges them, end to front, until either there is only a single curve left, or there is no
    /// merge left that doesn't exceed the optional maximum distance threshold provided.
    ///
    /// Merges happen one at a time. The algorithm checks the end-to-start distance of every pair
    /// of curves in the collection.  If the distance of the shortest merge is below the optional
    /// `max_dist` threshold, then the two curves are merged and the process repeats. The lower
    /// tolerance of the two curves is used for the new curve.
    ///
    /// Note:
    ///  - **Closed curves are ignored in this process!**
    ///  - Merges only occur with the end of one curve joined to the front of the next, meaning that
    ///    the order of the points in both curves will stay the same. There are no front-to-front or
    ///    end-to-end merges.
    ///
    /// # Arguments
    ///
    /// * `curves`: a vector of `Curve2` instances to perform the chain merge on
    /// * `max_dist`: an optional maximum merge distance between the end of one curve and the front
    ///   of the candidate merge.
    ///
    /// returns: Result<Vec<Curve2, Global>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn chain_merge(mut curves: Vec<Curve2>, max_dist: Option<f64>) -> Result<Vec<Curve2>> {
        let mut no_op = false;
        let max_dist = max_dist.unwrap_or(f64::MAX);

        while !no_op {
            no_op = true;
            let mut distances = Vec::new();
            for i in 0..curves.len() {
                // Check that the front candidate isn't closed
                let front_candidate = &curves[i];
                if front_candidate.is_closed {
                    continue;
                }

                for k in 0..curves.len() {
                    if i == k {
                        continue;
                    }
                    let back_candidate = &curves[k];
                    if back_candidate.is_closed {
                        continue;
                    }

                    let d = dist(&front_candidate.at_back(), &back_candidate.at_front());
                    distances.push((d, i, k));
                }
            }

            if distances.is_empty() {
                continue;
            }

            let (d, i, k) = distances
                .iter()
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
                .unwrap();

            if *d <= max_dist {
                no_op = false;
                let points = [curves[*i].points(), curves[*k].points()].concat();
                let tol = curves[*i].tol.min(curves[*k].tol);
                let merged = Curve2::from_points(&points, tol, false)?;

                curves.remove(*i.max(k));
                curves.remove(*i.min(k));
                curves.push(merged);
            }
        }

        Ok(curves)
    }
}


#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use rand::prelude::SliceRandom;
    use rand::rng;
    use crate::curve2;
    use super::*;

    fn curve0() -> Curve2 {
        // Open curve 9 long
        curve2!((2, 2), (7, 2), (7, 6)).unwrap()
    }

    fn curve1() -> Curve2 {
        // Open curve 7 long
        curve2!(tol: 1e-4, closed: false; (5, 6), (1, 6), (1, 9)).unwrap()
    }

    fn curve2() -> Curve2 {
        // Open curve 6 long
        curve2!((7, 9), (13, 9)).unwrap()
    }

    fn curve4() -> Curve2 {
        // Closed curve 10 long
        curve2!(tol: 1e-4, closed: true; (8, 6), (8, 2), (9, 2), (9, 6)).unwrap()
    }

    fn curve5() -> Curve2 {
        // Closed curve 12 long
        curve2!(tol: 1e-4, closed: true; (-1, 9), (-2, 9), (-2, 4), (-1, 4)).unwrap()
    }

    #[test]
    fn check_sample_curves() {
        // This is just a pre-test to verify that I constructed the sample geometry correctly
        assert_relative_eq!(9.0, curve0().length(), epsilon = 1e-10);
        assert_relative_eq!(7.0, curve1().length(), epsilon = 1e-10);
        assert_relative_eq!(6.0, curve2().length(), epsilon = 1e-10);
        assert_relative_eq!(10.0, curve4().length(), epsilon = 1e-10);
        assert_relative_eq!(12.0, curve5().length(), epsilon = 1e-10);
        assert!(!curve0().is_closed);
        assert!(!curve1().is_closed);
        assert!(!curve2().is_closed);
        assert!(curve4().is_closed);
        assert!(curve5().is_closed);
    }

    fn find_curve(items: &[Curve2], len: f64, closed: bool) -> Option<&Curve2> {
        items.iter().find(|c| (c.length() - len).abs() < 1e-6 && c.is_closed == closed)
    }

    #[test]
    fn chain_merge_shuffle_no_max() {
        let mut rng = rng();
        let samples = vec![curve0(), curve1(), curve2(), curve4(), curve5()];
        for _ in 0..100 {
            let mut shuffled = samples.clone();
            shuffled.shuffle(&mut rng);

            let result = Curve2::chain_merge(shuffled, None).unwrap();
            assert_eq!(result.len(), 3);

            assert!(find_curve(&result, 30.0, false).is_some(), "Didn't find single merged curve");
            assert!(find_curve(&result, 10.0, true).is_some(), "Didn't find first closed curve");
            assert!(find_curve(&result, 12.0, true).is_some(), "Didn't find second closed curve");

            let merged = find_curve(&result, 30.0, false).unwrap();
            assert_relative_eq!(merged.tol, 1e-6, epsilon = 1e-10);
        }
    }

    #[test]
    fn chain_merge_shuffle_with_max() {
        let mut rng = rng();
        let samples = vec![curve0(), curve1(), curve2(), curve4(), curve5()];
        for _ in 0..100 {
            let mut shuffled = samples.clone();
            shuffled.shuffle(&mut rng);

            let result = Curve2::chain_merge(shuffled, Some(2.0)).unwrap();
            assert_eq!(result.len(), 4);

            assert!(find_curve(&result, 18.0, false).is_some(), "Didn't find single merged curve");
            assert!(find_curve(&result, 6.0, false).is_some(), "Didn't find single un-merged curve");
            assert!(find_curve(&result, 10.0, true).is_some(), "Didn't find first closed curve");
            assert!(find_curve(&result, 12.0, true).is_some(), "Didn't find second closed curve");

            let merged = find_curve(&result, 18.0, false).unwrap();
            assert_relative_eq!(merged.tol, 1e-6, epsilon = 1e-10);
        }
    }


}