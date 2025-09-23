use crate::Result;
use crate::common::{IndexMask, PCoords};
use kdtree::KdTree as KdTreeInner;
use kdtree::distance::squared_euclidean;
use uuid::Uuid;

// fn check_tree() {
//     let a: ([f64; 2], usize) = ([0f64, 0f64], 0);
//     let b: ([f64; 2], usize) = ([1f64, 1f64], 1);
//     let c: ([f64; 2], usize) = ([2f64, 2f64], 2);
//     let d: ([f64; 2], usize) = ([3f64, 3f64], 3);
//
//     let dimensions = 2;
//     let mut kdtree = KdTreeInner::new(dimensions);
//
//     kdtree.add(&a.0, a.1).unwrap();
//     kdtree.add(&b.0, b.1).unwrap();
//     kdtree.add(&c.0, c.1).unwrap();
//     kdtree.add(&d.0, d.1).unwrap();
// }

/// A KD tree associated with a unique UUID, such that it can be checked to be matched against a
/// specific entity.  The idea is that the UUID of the associated object will change if the object
/// points are modified, and thus the matched tree can be validated before use.
pub struct MatchedTree<const D: usize> {
    tree_uuid: Uuid,
    tree: KdTree<D>,
}

impl<const D: usize> MatchedTree<D> {
    pub fn new(tree_uuid: Uuid, tree: KdTree<D>) -> Self {
        Self { tree_uuid, tree }
    }

    pub fn tree_uuid(&self) -> Uuid {
        self.tree_uuid
    }

    pub fn tree(&self) -> &KdTree<D> {
        &self.tree
    }
}

pub trait KdTreeSearch<const D: usize> {
    fn nearest_one(&self, point: &impl PCoords<D>) -> (usize, f64);
    fn nearest(&self, point: &impl PCoords<D>, count: usize) -> Vec<(usize, f64)>;
    fn within(&self, point: &impl PCoords<D>, radius: f64) -> Vec<(usize, f64)>;

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// An immutable k-dimensional tree for fast searches on points in D dimensions
pub struct KdTree<const D: usize> {
    tree: KdTreeInner<f64, usize, [f64; D]>,
}

impl<const D: usize> KdTree<D> {
    /// Create a new immutable kd-tree from a list of points.
    ///
    /// # Arguments
    ///
    /// * `points`: A slice of points.
    ///
    /// returns: KdTree<{ D }>
    pub fn new(points: &[impl PCoords<D>]) -> Result<Self> {
        let mut entries: Vec<[f64; D]> = Vec::with_capacity(points.len());
        for p in points {
            entries.push(p.coords().into());
        }
        let mut tree = KdTreeInner::new(D);
        for (i, e) in entries.iter().enumerate() {
            tree.add(*e, i)?;
        }
        Ok(Self { tree })
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
    fn nearest_one(&self, point: &impl PCoords<D>) -> (usize, f64) {
        if let Ok(item) = self
            .tree
            .nearest(point.coords().as_slice(), 1, &squared_euclidean)
        {
            let (d, u) = &item[0];
            (**u, d.sqrt())
        } else {
            (usize::MAX, f64::INFINITY)
        }
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
    fn nearest(&self, point: &impl PCoords<D>, count: usize) -> Vec<(usize, f64)> {
        let result = self
            .tree
            .nearest(point.coords().as_slice(), count, &squared_euclidean);

        if let Ok(neighbors) = result {
            neighbors.iter().map(|(d, u)| (**u, d.sqrt())).collect()
        } else {
            Vec::new()
        }
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
    fn within(&self, point: &impl PCoords<D>, radius: f64) -> Vec<(usize, f64)> {
        if let Ok(result) = self.tree.within(
            point.coords().as_slice(),
            radius * radius,
            &squared_euclidean,
        ) {
            result.iter().map(|(d, u)| (**u, d.sqrt())).collect()
        } else {
            Vec::new()
        }
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
    pub fn new(all_points: &[impl PCoords<D>], mask: &IndexMask) -> Result<Self> {
        if mask.len() != all_points.len() {
            panic!("Mask length must match the length of all_points");
        }

        let mut points = Vec::new();
        let mut index_map = Vec::new();
        for i in mask.iter_true() {
            points.push(all_points[i].coords());
            index_map.push(i);
        }
        let tree = KdTree::new(&points)?;
        Ok(Self { tree, index_map })
    }
}

impl<const D: usize> KdTreeSearch<D> for PartialKdTree<D> {
    fn nearest_one(&self, point: &impl PCoords<D>) -> (usize, f64) {
        let (i, d) = self.tree.nearest_one(point);
        (self.index_map[i], d)
    }

    fn nearest(&self, point: &impl PCoords<D>, count: usize) -> Vec<(usize, f64)> {
        let result = self.tree.nearest(point, count);
        result
            .iter()
            .map(|(i, d)| (self.index_map[*i], *d))
            .collect::<Vec<_>>()
    }

    fn within(&self, point: &impl PCoords<D>, radius: f64) -> Vec<(usize, f64)> {
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
    use crate::Point2;
    use crate::common::indices::index_vec;
    use approx::assert_relative_eq;
    use rand::prelude::SliceRandom;

    #[test]
    fn kd_tree_build() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::new(&points).expect("KD tree creation failed");
        assert_eq!(tree.len(), 3);
    }

    #[test]
    fn kd_tree_nearest() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::new(&points).expect("KD tree creation failed");
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
        let tree = KdTree::new(&points).expect("KD tree creation failed");
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

        let fixed_tree = KdTree::new(&points).expect("KD tree creation failed");

        for _ in 0..1000 {
            let mut test_select = index_vec(None, points.len());
            test_select.shuffle(&mut rand::rng());

            // Use only half of the points for the partial tree
            let mut mask = IndexMask::new(points.len(), false);
            for &i in test_select.iter().take(points.len() / 2) {
                mask.set(i, true);
            }
            let indices = mask.to_indices();

            let partial_tree = PartialKdTree::new(&points, &mask).expect("KD tree creation failed");

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
