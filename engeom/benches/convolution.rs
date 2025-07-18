use criterion::{Criterion, criterion_group, criterion_main};
use engeom::na::DMatrix;
use engeom::raster2::{RasterKernel, ScalarRaster};
use std::hint::black_box;

fn kernel_convolution(c: &mut Criterion) {
    let target = striped_target();

    let max = target.max() * 2.0;
    let raster = ScalarRaster::from_matrix(&target, 1.0, -max, max);

    let kernel = RasterKernel::gaussian(1.0);

    c.bench_function("kernel convolution", |b| {
        b.iter(|| {
            let _result = black_box(&kernel).convolve(&raster, false, true);
        })
    });
}

fn striped_target() -> DMatrix<f64> {
    // Create a striped target matrix with some values
    let mut target = DMatrix::zeros(400, 600);
    for i in 0..target.nrows() {
        for j in 0..target.ncols() {
            target[(i, j)] = ((i + j) % 50) as f64 / 100.0;
        }
    }
    for i in 0..target.nrows() {
        for j in 0..target.ncols() {
            if (i + j) % 10 == 0 {
                target[(i, j)] = f64::NAN; // Mask some pixels
            }
        }
    }
    target
}

criterion_group!(benches, kernel_convolution,);
criterion_main!(benches);
