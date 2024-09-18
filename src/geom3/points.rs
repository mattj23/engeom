use kiddo::SquaredEuclidean;
use crate::geom3::kd_tree3::KdTree3;
use crate::Point3;

pub fn points_sample_poisson_disk(points: &[Point3], radius: f64, indices: &[usize]) -> Vec<usize> {
    let mut results = Vec::new();
    let mut tree = KdTree3::new();
    let threshold = radius * radius;
    for i in indices {
        let p = points[*i];
        let r = tree.nearest_one::<SquaredEuclidean>(&[p.x, p.y, p.z]);
        if r.distance > threshold {
            results.push(*i);
            tree.add(&[p.x, p.y, p.z], *i);
        }
    }

    results
}

