//! Build and query benchmarks for the in-tree k-d tree backend.
//!
//! Points are Fibonacci-sphere distributed rather than uniform in a volume. Every point cloud this
//! library sees is a scan of a surface, and a k-d tree over a 2-manifold embedded in 3-space
//! behaves differently from one over a solid.
//!
//! This exists because the library originally used `kiddo`, then I ran into some issues with data
//! packed into a single axis, then again with correctness, and I switched to `kd-tree`.  Then, when
//! revisiting the point cloud part of the library, I discovered that `kd-tree` is comparatively
//! very slow.  This benchmark is to leave some tooling up any time I need to re-verify the kd tree
//! performance.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use engeom::common::kd_tree::KdTreeSearch;
use engeom::{KdTree3, Point3};
use std::hint::black_box;

const SIZES: [usize; 3] = [10_000, 100_000, 1_000_000];

/// Queries per timed batch, large enough to swamp timer resolution.
const QUERIES: usize = 1_000;

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

/// A radius holding a handful of neighbors at the given density, matching how `within` is actually
/// used (a small neighborhood for normal estimation) rather than one returning thousands of points,
/// which would measure allocation instead of search.
fn neighborhood_radius(n: usize) -> f64 {
    3.0 * 200.0 / (n as f64).sqrt()
}

fn build(c: &mut Criterion) {
    let mut group = c.benchmark_group("kd_tree build");
    group.sample_size(10);

    for &n in SIZES.iter() {
        let points = fibonacci_sphere(n, 100.0);
        group.bench_with_input(BenchmarkId::from_parameter(n), &points, |b, points| {
            b.iter(|| KdTree3::try_new(black_box(points)).expect("build failed"))
        });
    }

    group.finish();
}

fn query(c: &mut Criterion) {
    let mut group = c.benchmark_group("kd_tree query");

    for &n in SIZES.iter() {
        let points = fibonacci_sphere(n, 100.0);
        let tree = KdTree3::try_new(&points).expect("build failed");
        let radius = neighborhood_radius(n);
        let queries: Vec<Point3> = (0..QUERIES)
            .map(|i| points[(i * (n / QUERIES).max(1)) % n])
            .collect();

        group.bench_with_input(
            BenchmarkId::new("nearest_one", n),
            &queries,
            |b, queries| b.iter(|| queries.iter().map(|q| tree.nearest_one(q).0).sum::<usize>()),
        );

        group.bench_with_input(BenchmarkId::new("nearest_7", n), &queries, |b, queries| {
            b.iter(|| {
                queries
                    .iter()
                    .map(|q| tree.nearest(q, 7).len())
                    .sum::<usize>()
            })
        });

        group.bench_with_input(BenchmarkId::new("within", n), &queries, |b, queries| {
            b.iter(|| {
                queries
                    .iter()
                    .map(|q| tree.within(q, radius).len())
                    .sum::<usize>()
            })
        });
    }

    group.finish();
}

criterion_group!(benches, build, query);
criterion_main!(benches);
