use kiddo::immutable::float::kdtree;
use kiddo::SquaredEuclidean;
use parry3d_f64::na::Point;
use std::num::NonZero;

pub trait KdTreeSearch<const D: usize> {
    fn nearest_one(&self, point: &Point<f64, D>) -> (usize, f64);
    fn nearest(&self, point: &Point<f64, D>, count: NonZero<usize>) -> Vec<(usize, f64)>;
    fn within(&self, point: &Point<f64, D>, radius: f64) -> Vec<(usize, f64)>;

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// An immutable k-dimensional tree for fast searches on points in D dimensions
pub struct KdTree<const D: usize> {
    tree: kdtree::ImmutableKdTree<f64, usize, D, 32>,
}

impl<const D: usize> KdTree<D> {
    /// Create a new immutable kd-tree from a list of points.
    ///
    /// # Arguments
    ///
    /// * `points`: A slice of points.
    ///
    /// returns: KdTree<{ D }>
    pub fn new(points: &[Point<f64, D>]) -> Self {
        let mut entries: Vec<[f64; D]> = Vec::with_capacity(points.len());
        for p in points {
            entries.push(p.coords.into());
        }
        Self {
            tree: kdtree::ImmutableKdTree::new_from_slice(&entries),
        }
    }
}

impl<const D: usize> KdTreeSearch<D> for KdTree<D> {
    /// Find the nearest point in the kd-tree to a given test point, returning the index of the
    /// nearest point and the distance to it.
    ///
    /// # Arguments
    ///
    /// * `point`: A test point to find the nearest point to.
    ///
    /// returns: (usize, f64)
    ///
    /// # Examples
    ///
    /// ```
    ///
    ///
    /// ```
    fn nearest_one(&self, point: &Point<f64, D>) -> (usize, f64) {
        let result = self
            .tree
            .nearest_one::<SquaredEuclidean>(&point.coords.into());
        (result.item, result.distance.sqrt())
    }

    /// Find the nearest `count` points in the kd-tree to a given point.
    ///
    /// # Arguments
    ///
    /// * `point`:
    /// * `count`:
    ///
    /// returns: Vec<(usize, f64), Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn nearest(&self, point: &Point<f64, D>, count: NonZero<usize>) -> Vec<(usize, f64)> {
        let result = self
            .tree
            .nearest_n::<SquaredEuclidean>(&point.coords.into(), count);
        result
            .iter()
            .map(|r| (r.item, r.distance.sqrt()))
            .collect::<Vec<_>>()
    }

    /// Find all points within a given radius of a point.
    ///
    /// # Arguments
    ///
    /// * `point`:
    /// * `radius`:
    ///
    /// returns: Vec<(usize, f64), Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn within(&self, point: &Point<f64, D>, radius: f64) -> Vec<(usize, f64)> {
        let result = self
            .tree
            .within::<SquaredEuclidean>(&point.coords.into(), radius * radius);
        result
            .iter()
            .map(|r| (r.item, r.distance.sqrt()))
            .collect::<Vec<_>>()
    }

    /// Get the number of points in the kd-tree.
    fn len(&self) -> usize {
        self.tree.size()
    }
}

/// A wrapper around a KdTree and a list of indices into an original list of points. This allows
/// for searches to be performed on a subset of the original points with indices returned which
/// correspond to the original list.
pub struct PartialKdTree<const D: usize> {
    tree: KdTree<D>,
    index_map: Vec<usize>,
}

impl<const D: usize> PartialKdTree<D> {
    /// Create a new partial kd-tree from a list of points and a list of indices into the original.
    ///
    /// The `indices` array should be a list of indices into the `all_points` array. The tree will
    /// be built using only the points at those indices, *however*, the indices returned by the
    /// search methods will be the indices into the original `all_points` array.
    ///
    /// # Arguments
    ///
    /// * `all_points`: A slice of points.
    /// * `indices`: A slice of indices into the `all_points` array to use for the tree.
    ///
    /// returns: PartialKdTree<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new(all_points: &[Point<f64, D>], indices: &[usize]) -> Self {
        let points = indices.iter().map(|i| all_points[*i]).collect::<Vec<_>>();
        let tree = KdTree::new(&points);
        let index_map = indices.to_vec();
        Self { tree, index_map }
    }
}

impl<const D: usize> KdTreeSearch<D> for PartialKdTree<D> {
    fn nearest_one(&self, point: &Point<f64, D>) -> (usize, f64) {
        let (i, d) = self.tree.nearest_one(point);
        (self.index_map[i], d)
    }

    fn nearest(&self, point: &Point<f64, D>, count: NonZero<usize>) -> Vec<(usize, f64)> {
        let result = self.tree.nearest(point, count);
        result
            .iter()
            .map(|(i, d)| (self.index_map[*i], *d))
            .collect::<Vec<_>>()
    }

    fn within(&self, point: &Point<f64, D>, radius: f64) -> Vec<(usize, f64)> {
        let result = self.tree.within(point, radius);
        result
            .iter()
            .map(|(i, d)| (self.index_map[*i], *d))
            .collect::<Vec<_>>()
    }

    /// Get the number of points in the kd-tree.
    fn len(&self) -> usize {
        self.tree.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::indices::index_vec;
    use crate::Point2;
    use approx::assert_relative_eq;
    use rand::prelude::SliceRandom;

    #[test]
    fn kd_tree_build() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::new(&points);
        assert_eq!(tree.len(), 3);
    }

    #[test]
    fn kd_tree_nearest() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::new(&points);
        let (i, d) = tree.nearest_one(&Point2::new(1.25, 1.25));
        assert_eq!(i, 1);
        assert_relative_eq!(d, (0.25 * 0.25 * 2.0_f64).sqrt(), epsilon = 1e-6);
    }

    #[test]
    fn kd_tree_check_distances_within() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
        ];
        let tree = KdTree::new(&points);
        let within = tree.within(&Point2::new(3.5, 0.0), 2.0);
        assert_eq!(within.len(), 1);
        assert_eq!(within[0].0, 2);
        assert_relative_eq!(within[0].1, 1.5, epsilon = 1e-6);
    }

    #[test]
    fn partial_kd_tree_maps() {
        let points = (0..20)
            .flat_map(|i| (0..20).map(move |j| Point2::new(i as f64, j as f64)))
            .collect::<Vec<_>>();

        let fixed_tree = KdTree::new(&points);

        for _ in 0..1000 {
            let mut indices = index_vec(None, points.len());
            indices.shuffle(&mut rand::rng());
            let indices = &indices[0..(points.len() / 2)];

            let partial_tree = PartialKdTree::new(&points, &indices);

            for &i in indices.iter() {
                let p = &points[i];
                let (j, d) = fixed_tree.nearest_one(p);
                let (k, e) = partial_tree.nearest_one(p);
                assert_eq!(i, j);
                assert_eq!(i, k);
                assert_relative_eq!(0.0, e, epsilon = 1e-6);
                assert_relative_eq!(0.0, e, epsilon = 1e-6);
            }
        }
    }
}
