use criterion::{Criterion, criterion_group, criterion_main};
use engeom::common::IndexMask;
use std::hint::black_box;

const N: usize = 10_000_000;

fn index_mask_set(c: &mut Criterion) {
    let to_set = (0..N).filter(|i| i % 2 == 0).collect::<Vec<_>>();
    let mut mask = IndexMask::new(N, false);

    c.bench_function("index_mask set", |b| {
        b.iter(|| {
            for x in to_set.iter() {
                black_box(&mut mask).set(*x, true);
            }
        })
    });
}

fn index_mask_get(c: &mut Criterion) {
    let mask = prep_mask(2);

    c.bench_function("index_mask get", |b| {
        b.iter(|| {
            for i in 0..N {
                let _x = black_box(&mask).get(i);
            }
        })
    });
}

fn index_mask_flip(c: &mut Criterion) {
    let mut mask = prep_mask(2);
    c.bench_function("index_mask flip", |b| {
        b.iter(|| {
            black_box(&mut mask).not_mut();
        })
    });
}

fn index_mask_to_vec(c: &mut Criterion) {
    let mask = prep_mask(2);

    c.bench_function("index_mask to vec", |b| {
        b.iter(|| {
            let _v = black_box(&mask).to_indices();
        })
    });
}

fn prep_mask(m: usize) -> IndexMask {
    let to_set = (0..N).filter(|i| i % m == 0).collect::<Vec<_>>();
    let mut mask = IndexMask::new(N, false);
    for x in to_set.iter() {
        mask.set(*x, true);
    }

    mask
}

criterion_group!(
    benches,
    index_mask_set,
    index_mask_get,
    index_mask_flip,
    index_mask_to_vec
);
criterion_main!(benches);
