use crate::geom2::Point2;
use kiddo::SquaredEuclidean;

pub type KdTree2 = kiddo::float::kdtree::KdTree<f64, usize, 2, 128, u32>;

pub fn to_kd_tree2(points: &[Point2]) -> KdTree2 {
    let entries = points.iter().map(|p| [p.x, p.y]).collect::<Vec<_>>();
    (&entries).into()
}

pub fn kd_tree_nearest_2d(tree: &KdTree2, point: &Point2) -> usize {
    let result = tree.nearest_one::<SquaredEuclidean>(&[point.x, point.y]);
    result.item
}

pub fn kd_tree_within_2d(tree: &KdTree2, point: &Point2, radius_squared: f64) -> Vec<usize> {
    let result = tree.within::<SquaredEuclidean>(&[point.x, point.y], radius_squared);
    result.iter().map(|r| r.item).collect::<Vec<_>>()
}

pub fn sample_poisson_disk_2d(points: &[Point2], radius: f64, indices: &[usize]) -> Vec<usize> {
    let mut results = Vec::new();
    let mut tree = KdTree2::new();
    let threshold = radius * radius;
    for i in indices {
        let p = points[*i];
        let result = tree.nearest_one::<SquaredEuclidean>(&[p.x, p.y]);
        if result.distance > threshold {
            results.push(*i);
            tree.add(&[p.x, p.y], *i);
        }
    }

    results
}
