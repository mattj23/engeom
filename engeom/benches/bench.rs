use std::collections::HashSet;
use engeom::common::IndexMask;
use divan::{Bencher, black_box};

const N: usize = 10_000_000;

fn main() {
    divan::main();
}

/// Test the speed at which an index mask can set half of its entities
#[divan::bench]
fn index_mask_set(bencher: Bencher) {
    let to_set = (0..N).filter(|i| i % 2 == 0).collect::<Vec<_>>();
    let mut mask = IndexMask::new(N, false);

    bencher.bench_local(move || {
        for x in to_set.iter() {
            black_box(&mut mask).set(*x, true);
        }
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

#[divan::bench]
fn index_mask_get(bencher: Bencher) {
    let mask = prep_mask(2);
    bencher.bench_local(move || {
        for i in 0..N {
            let x = black_box(&mask).get(i);
        }
    });
}

#[divan::bench]
fn index_mask_flip(bencher: Bencher) {
    let mut mask = prep_mask(2);
    bencher.bench_local(move || {
        black_box(&mut mask).flip();
    });
}

#[divan::bench]
fn index_mask_to_vec(bencher: Bencher) {
    let mask = prep_mask(2);

    bencher.bench_local(move || {
        let v = black_box(&mask).to_indices();
    });
}
