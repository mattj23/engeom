use std::collections::HashSet;
use engeom::common::IndexMask;
use divan::{Bencher, black_box};

const N: usize = 1_000_000;

fn main() {
    // Run registered benchmarks.
    divan::main();
}

/// Test the speed at which an index mask can set half of its entities
#[divan::bench]
fn index_mask_set(bencher: Bencher) {
    let to_set = (0..N).filter(|i| i % 2 == 0).collect::<Vec<_>>();
    let mut mask = HashSet::<usize>::new();

    bencher.bench_local(move || {
        for x in to_set.iter() {
            black_box(&mut mask).insert(*x);
        }
    });
}

fn prep_mask(m: usize) -> HashSet<usize> {
    let to_set = (0..N).filter(|i| i % m == 0).collect::<Vec<_>>();
    let mut mask = HashSet::<usize>::new();
    for x in to_set.iter() {
        mask.insert(*x);
    }

    mask
}

#[divan::bench]
fn index_mask_get(bencher: Bencher) {
    let mask = prep_mask(2);
    bencher.bench_local(move || {
        for i in 0..N {
            let x = black_box(&mask).contains(&i);
        }
    });
}

#[divan::bench]
fn index_mask_to_vec(bencher: Bencher) {
    let mask = prep_mask(2);

    bencher.bench_local(move || {
        let mut v = Vec::new();
        for i in black_box(&mask).iter() {
            v.push(black_box(*i));
        }
    });
}
