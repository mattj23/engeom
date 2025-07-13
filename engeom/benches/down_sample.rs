use criterion::{Criterion, criterion_group, criterion_main};
use engeom::Mesh;
use std::hint::black_box;

const N: usize = 10_000_000;

fn poisson_downsample(c: &mut Criterion) {
    c.bench_function("downsample mesh_poisson", |b| {
        let mesh = Mesh::create_sphere(100.0, 500, 500);
        b.iter(|| {
            let _results = black_box(&mesh).sample_poisson(5.0);
        })
    });
}

criterion_group!(benches, poisson_downsample,);
criterion_main!(benches);
