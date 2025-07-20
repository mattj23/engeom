use colorgrad::preset::turbo;
use criterion::{Criterion, criterion_group, criterion_main};
use engeom::na::DMatrix;
use engeom::raster2::{
    Point2I, Point2IIndexAccess, RasterKernel, RasterMask, ScalarRaster, SizeForIndex,
    ball_rolling_background,
};
use std::hint::black_box;
use std::path::Path;

fn raster_ball_rolling(c: &mut Criterion) {
    let target = build_target();
    target
        .render_with_cmap(
            &Path::new("D:/temp/k/rbtest.png"),
            &turbo(),
            Some((-2.0, 2.0)),
        )
        .unwrap();

    c.bench_function("raster2 ball_rolling", |b| {
        b.iter(|| {
            let _result = ball_rolling_background(black_box(&target), 2.0);
        })
    });
}

fn build_target() -> ScalarRaster {
    let mut target = DMatrix::zeros(300, 300);
    for p in target.iter_indices() {
        let v0 = (p.x as f64 * 50.0).sin();
        let v1 = (p.y as f64 * 50.0).cos();
        target.set_at(p, v0 + v1).unwrap();
    }
    let mut raster = ScalarRaster::from_matrix(&target, 0.1, -2.0, 2.0);
    let mut mask = RasterMask::empty_like(&raster.mask.buffer);
    mask.draw_circle_mut(Point2I::new(150, 150), 130, true, true);

    for p in mask.iter_all() {
        if !mask.get_point(p) {
            raster.set_f_at(p, f64::NAN).unwrap();
        }
    }

    raster
}

criterion_group!(benches, raster_ball_rolling);
criterion_main!(benches);
