//! A thin wrapper over the k-d tree backend, exposing the three queries `engeom` currently uses.
//! The backend used to be `kd-tree`, but I switched back to `kiddo` after doing some benchmarks and
//! seeing how much faster it is.  There are some additional tests in the library now to catch the
//! correctness bugs I was finding last year, but it looks like they got fixed.
//!
//! The backend is `kiddo`'s [`ImmutableKdTree`], which bulk-loads a point set that is known 100% up
//! front and is then only queried, which is every use in this crate. In addition to speed, bulk
//! loading picks split planes from the whole set at the beginning, so the tree ends up balanced
//! regardless of point order.  This avoids some of the dengenerate orderings that caused issues on
//! the `kd-tree` backend.
//!
//! <div class="warning">
//!
//! Important to note, distances cross this boundary as normal, non-squared Euclidean distances. For
//! optimization purposes, `kiddo` ingests and reports squared values, so distances are squared on
//! the way in and square-rooted (is that a word?) on the way out.
//!
//! </div>

use crate::Result;
use crate::common::{IndexMask, PCoords};
use kiddo::{ImmutableKdTree, SquaredEuclidean};
use std::num::NonZeroUsize;

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
    tree: ImmutableKdTree<f64, D>,

    /// The point count, kept alongside the tree because an empty tree is legal here and has to be
    /// detected before a query, not during one.
    len: usize,
}

impl<const D: usize> KdTree<D> {
    /// Build a k-d tree over a set of points.
    ///
    /// The returned indices of every query are indices into `points`. An empty set is allowed and
    /// produces a tree whose queries all come back empty.
    pub fn try_new(points: &[impl PCoords<D>]) -> Result<Self> {
        let entries: Vec<[f64; D]> = points.iter().map(|p| p.coords().into()).collect();

        // An empty slice has no split plane to choose, so the backend is not asked for one.
        if entries.is_empty() {
            return Ok(Self {
                tree: ImmutableKdTree::new_from_slice(&[])
                    .map_err(|e| format!("Failed to build an empty k-d tree: {e}"))?,
                len: 0,
            });
        }

        let tree = ImmutableKdTree::new_from_slice(&entries).map_err(|e| {
            format!(
                "Failed to build a k-d tree over {} points: {e}",
                entries.len()
            )
        })?;

        Ok(Self {
            tree,
            len: entries.len(),
        })
    }

    /// The query point as the fixed-size array the backend wants.
    fn key(point: &impl PCoords<D>) -> [f64; D] {
        point.coords().into()
    }
}

impl<const D: usize> KdTreeSearch<D> for KdTree<D> {
    /// Find the single nearest point to `point`, as an index into the set the tree was built over
    /// and the distance to it.
    ///
    /// An empty tree reports `(usize::MAX, f64::INFINITY)`, which is a sentinel rather than an
    /// index and will panic if used to subscript the point set. Check [`KdTreeSearch::is_empty`]
    /// first where an empty tree is possible.
    fn nearest_one(&self, point: &impl PCoords<D>) -> (usize, f64) {
        if self.len == 0 {
            return (usize::MAX, f64::INFINITY);
        }

        let found = self
            .tree
            .query(&Self::key(point))
            .nearest_one::<SquaredEuclidean<f64>>()
            .execute();

        (found.item as usize, found.distance.sqrt())
    }

    /// Find the `count` nearest points to `point`, ordered nearest first.
    ///
    /// Returns fewer than `count` entries when the tree holds fewer points, and none when `count`
    /// is zero.
    fn nearest(&self, point: &impl PCoords<D>, count: usize) -> Vec<(usize, f64)> {
        let Some(count) = NonZeroUsize::new(count.min(self.len)) else {
            return Vec::new();
        };

        self.tree
            .query(&Self::key(point))
            .nearest_n::<SquaredEuclidean<f64>>(count)
            .execute()
            .into_iter()
            .map(|f| (f.item as usize, f.distance.sqrt()))
            .collect()
    }

    /// Find every point within `radius` of `point`, ordered nearest first. The bound is inclusive.
    fn within(&self, point: &impl PCoords<D>, radius: f64) -> Vec<(usize, f64)> {
        if self.len == 0 {
            return Vec::new();
        }

        self.tree
            .query(&Self::key(point))
            .within::<SquaredEuclidean<f64>>(radius * radius)
            .execute()
            .into_iter()
            .map(|f| (f.item as usize, f.distance.sqrt()))
            .collect()
    }

    /// Get the number of points in the kd-tree.
    fn len(&self) -> usize {
        self.len
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
    /// The tree is built over only the points `mask` selects, *however* the indices returned by
    /// every query are indices into the full `all_points` array, not into the selected subset.
    ///
    /// # Panics
    ///
    /// If `mask` is not the same length as `all_points`.
    pub fn try_new(all_points: &[impl PCoords<D>], mask: &IndexMask) -> Result<Self> {
        if mask.len() != all_points.len() {
            panic!("Mask length must match the length of all_points");
        }

        let mut points = Vec::new();
        let mut index_map = Vec::new();
        for i in mask.iter_true() {
            points.push(all_points[i].coords());
            index_map.push(i);
        }
        let tree = KdTree::try_new(&points)?;
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
    use crate::common::indices::index_vec;
    use crate::common::points::dist;
    use crate::tests::stanford_bun_2;
    use crate::{Point2, Point3};
    use approx::assert_relative_eq;
    use rand::prelude::SliceRandom;

    #[test]
    fn kd_tree_build() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::try_new(&points).expect("KD tree creation failed");
        assert_eq!(tree.len(), 3);
    }

    #[test]
    fn kd_tree_nearest() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        ];
        let tree = KdTree::try_new(&points).expect("KD tree creation failed");
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
        let tree = KdTree::try_new(&points).expect("KD tree creation failed");
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

        let fixed_tree = KdTree::try_new(&points).expect("KD tree creation failed");

        for _ in 0..1000 {
            let mut test_select = index_vec(None, points.len());
            test_select.shuffle(&mut rand::rng());

            // Use only half of the points for the partial tree
            let mut mask = IndexMask::new(points.len(), false);
            for &i in test_select.iter().take(points.len() / 2) {
                mask.set(i, true);
            }
            let indices = mask.to_indices();

            let partial_tree =
                PartialKdTree::try_new(&points, &mask).expect("KD tree creation failed");

            for &i in indices.iter() {
                let p = &points[i];
                let (j, _d) = fixed_tree.nearest_one(p);
                let (k, e) = partial_tree.nearest_one(p);
                assert_eq!(i, j);
                assert_eq!(i, k);
                assert_relative_eq!(0.0, e, epsilon = 1e-6);
                assert_relative_eq!(0.0, e, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn partial_kd_tree_nearest() -> Result<()> {
        let mesh = stanford_bun_2();
        let even = (0..mesh.points().len()).step_by(2).collect::<Vec<_>>();
        let mask = IndexMask::try_from_indices(&even, mesh.points().len())?;
        let reduced_points = mask.clone_indices_of(mesh.points())?;
        let reduced_tree = KdTree::try_new(&reduced_points)?;

        let partial_tree = PartialKdTree::try_new(mesh.points(), &mask)?;
        let n = 3;

        for p in mesh.points().iter() {
            let mut result = partial_tree
                .nearest(p, n)
                .iter()
                .map(|x| x.0)
                .collect::<Vec<_>>();
            let mut reduced_result = reduced_tree
                .nearest(p, n)
                .iter()
                .map(|x| x.0)
                .collect::<Vec<_>>();
            result.sort();
            reduced_result.sort();

            assert_eq!(result.len(), reduced_result.len());
            for (a, b) in result.iter().zip(reduced_result.iter()) {
                assert_relative_eq!(mesh.points()[*a], reduced_points[*b], epsilon = 1e-6);
            }
        }

        Ok(())
    }

    #[test]
    fn kd_tree_nearest_one_matches_brute_force_on_stanford_bun_vertices() {
        let mesh = stanford_bun_2();
        let vertices = mesh.points().to_vec();
        let tree = KdTree::try_new(&vertices).expect("KD tree creation failed");

        for (expected_index, query) in vertices.iter().enumerate() {
            let (actual_index, actual_distance) = tree.nearest_one(query);
            assert_eq!(actual_index, expected_index);
            assert_relative_eq!(actual_distance, 0.0, epsilon = 1e-6);

            let expected = brute_force_nearest_n(&vertices, query, 1);
            assert_eq!(expected.len(), 1);
            assert_eq!(expected[0].0, expected_index);
            assert_relative_eq!(expected[0].1, 0.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn kd_tree_nearest_matches_brute_force_on_stanford_bun_vertices() {
        let mesh = stanford_bun_2();
        let vertices = mesh.points().to_vec();
        let tree = KdTree::try_new(&vertices).expect("KD tree creation failed");

        for query in &vertices {
            let expected = brute_force_nearest_n(&vertices, query, 3);
            let actual = tree.nearest(query, 3);

            assert_eq!(actual.len(), expected.len());
            for ((actual_index, actual_distance), (expected_index, expected_distance)) in
                actual.iter().zip(expected.iter())
            {
                assert_eq!(actual_index, expected_index);
                assert_relative_eq!(actual_distance, expected_distance, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn kd_tree_within_matches_brute_force_on_stanford_bun_vertices() {
        let mesh = stanford_bun_2();
        let vertices = mesh.points().to_vec();
        let tree = KdTree::try_new(&vertices).expect("KD tree creation failed");

        let search_dist = mesh.aabb().extents().x * 0.1;

        for query in vertices.iter() {
            let expected = vertices
                .iter()
                .filter(|p| dist(query, *p) < search_dist)
                .count();
            let tested = tree.within(query, search_dist).len();
            assert_eq!(tested, expected);
        }
    }

    /// A radius query must return every point inside the radius at a point count large enough to
    /// matter, checked against brute force.
    ///
    /// This exists because of the specific way the previous backend failed. `engeom` moved off
    /// `kiddo` in Sept 2025 over a radius query that did not return all of its neighbors *once the
    /// point count grew large*, and every other correctness test in this crate runs at a few
    /// thousand points: the Stanford bunny is about 8k, and `check_kiddo_bug` queries a Poisson
    /// sample of a couple of thousand. None of them reach the regime the bug lived in, so none of
    /// them would have caught it coming back.
    ///
    /// Brute force over every query would be O(N^2), so a small number of query points is checked
    /// exhaustively against the whole set instead. The radius is chosen to return a substantial
    /// neighborhood, since the failure mode was omission rather than a wrong answer.
    #[test]
    fn kd_tree_within_matches_brute_force_at_scale() {
        const N: usize = 500_000;
        const QUERIES: usize = 100;

        let points = fibonacci_sphere(N, 100.0);
        let tree = KdTree::try_new(&points).expect("KD tree creation failed");

        // Roughly a few hundred neighbors at this density.
        let radius = 30.0 * 200.0 / (N as f64).sqrt();

        for q in (0..QUERIES).map(|i| points[(i * (N / QUERIES)) % N]) {
            let mut expected: Vec<usize> = points
                .iter()
                .enumerate()
                .filter(|(_, p)| dist(&q, *p) <= radius)
                .map(|(i, _)| i)
                .collect();

            let mut actual: Vec<usize> = tree.within(&q, radius).iter().map(|(i, _)| *i).collect();

            assert!(
                expected.len() > 50,
                "radius {radius} returned only {} neighbors, too few to be a meaningful check",
                expected.len()
            );

            expected.sort_unstable();
            actual.sort_unstable();
            assert_eq!(
                actual, expected,
                "radius query missed or invented neighbors"
            );
        }
    }

    /// Deterministic surface-distributed points, so a failure is reproducible.
    fn fibonacci_sphere(n: usize, radius: f64) -> Vec<Point3> {
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        (0..n)
            .map(|i| {
                let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0;
                let r = (1.0 - y * y).max(0.0).sqrt();
                let theta = golden * i as f64;
                Point3::new(
                    radius * theta.cos() * r,
                    radius * y,
                    radius * theta.sin() * r,
                )
            })
            .collect()
    }

    /// Returns the `count` nearest points to `query` as `(index, distance)` pairs,
    /// sorted from nearest to farthest.
    ///
    /// If `count` is larger than the number of points, all points are returned.
    ///
    /// This is brute force: it checks every point in the collection.
    pub fn brute_force_nearest_n(
        points: &[Point3],
        query: &Point3,
        count: usize,
    ) -> Vec<(usize, f64)> {
        let mut neighbors: Vec<(usize, f64)> = points
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let dx = p.x - query.x;
                let dy = p.y - query.y;
                let dz = p.z - query.z;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                (i, dist)
            })
            .collect();

        neighbors.sort_by(|a, b| a.1.total_cmp(&b.1));
        neighbors.truncate(count);
        neighbors
    }
}
