use criterion::{Criterion, criterion_group, criterion_main};
use engeom::Mesh3;
use std::hint::black_box;

#[allow(dead_code)] // retained bench parameter
const N: usize = 10_000_000;

fn poisson_downsample(c: &mut Criterion) {
    c.bench_function("downsample mesh_poisson", |b| {
        let mesh = Mesh3::create_sphere(100.0, 4.0e-3).unwrap();
        b.iter(|| {
            let _results = black_box(&mesh).sample_surface_poisson(5.0, None);
        })
    });
}

criterion_group!(benches, poisson_downsample,);
criterion_main!(benches);
