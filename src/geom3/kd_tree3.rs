use crate::geom3::Point3;
use kiddo::SquaredEuclidean;

const BUCKET_SIZE: usize = 256;
pub type KdTree3 = kiddo::float::kdtree::KdTree<f64, usize, 3, BUCKET_SIZE, u32>;

pub fn to_kd_tree3(points: &[Point3]) -> KdTree3 {
    let entries = points.iter().map(|p| [p.x, p.y, p.z]).collect::<Vec<_>>();
    (&entries).into()
}

pub fn kd_tree_nearest_3d(tree: &KdTree3, point: &Point3) -> usize {
    let result = tree.nearest_one::<SquaredEuclidean>(&[point.x, point.y, point.z]);
    result.item
}

pub fn kd_tree_within_3d(tree: &KdTree3, point: &Point3, radius_squared: f64) -> Vec<usize> {
    let result = tree.within::<SquaredEuclidean>(&[point.x, point.y, point.z], radius_squared);
    result.iter().map(|r| r.item).collect::<Vec<_>>()
}

/// A wrapper around a KdTree3 and a list of indices into an original list of points. This allows
/// for searches to be performed on a subset of the original points with indices returned which
/// correspond to the original list.
pub struct PartialKdTree3 {
    tree: KdTree3,
    index_map: Vec<usize>,
}

impl PartialKdTree3 {
    pub fn new(all_points: &[Point3], indices: &[usize]) -> Self {
        let points = indices.iter().map(|i| all_points[*i]).collect::<Vec<_>>();
        let tree = to_kd_tree3(&points);
        let index_map = indices.to_vec();
        Self { tree, index_map }
    }

    pub fn nearest(&self, point: &Point3) -> usize {
        let result = kd_tree_nearest_3d(&self.tree, point);
        self.index_map[result]
    }

    pub fn nearest_as_new(&self, point: &Point3) -> usize {
        kd_tree_nearest_3d(&self.tree, point)
    }
}
