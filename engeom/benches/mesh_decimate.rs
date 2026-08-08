//! Benchmarks for tolerance-bounded mesh decimation.
//!
//! Decimation has two things worth measuring and they trade against each other: how long a run
//! takes, and how much of the mesh it manages to remove. A change which halves the runtime by
//! refusing more collapses is not an improvement, so every case here prints the surviving face
//! count and the work counters alongside the timing. Read them together.
//!
//! Run `cargo bench -p engeom --bench mesh_decimate` for the timings. For the work counters and
//! the deviation of each run, use the profile harness in the library instead:
//!
//! ```text
//! cargo test -r -p engeom --lib decimate_profile -- --ignored --nocapture
//! ```

use criterion::{Criterion, criterion_group, criterion_main};
use engeom::Mesh3;
use engeom::geom3::HalfEdgeMesh3;
use engeom::geom3::half_edge3::{DecimateOpts, Placement};
use engeom::io::read_tc_mesh_from;
use std::hint::black_box;

/// The engine blade fixture: 21795 vertices, 43586 faces, millimetres, watertight, stored at a one
/// micron tolerance. Embedded rather than loaded from a path so the bench does not depend on the
/// working directory.
fn engine_blade() -> Mesh3 {
    let bytes = include_bytes!("../tests/data/engine-blade.tcmesh");
    let data = read_tc_mesh_from(&mut { bytes.as_slice() }).unwrap();
    Mesh3::from_data(data, false).unwrap()
}

fn decimate_tolerance(c: &mut Criterion) {
    let mesh = engine_blade();
    let mut group = c.benchmark_group("decimate_to_tolerance");

    // The cost per run scales with the number of collapses, so these bracket the useful range: at
    // 0.001 the blade barely moves, at 0.05 it is down to under a third of its faces. Read these
    // alongside the faces remaining, since a change which buys time by decimating less is not an
    // improvement; the profile harness in the library prints both together.
    for tol in [0.001, 0.01, 0.05] {
        group.bench_function(format!("tol_{tol}"), |b| {
            b.iter(|| {
                let mut he = HalfEdgeMesh3::try_from(black_box(&mesh)).unwrap();
                he.decimate_guaranteed(&DecimateOpts::to_tolerance(tol))
                    .unwrap()
            })
        });
    }

    group.finish();
}

fn decimate_placement(c: &mut Criterion) {
    let mesh = engine_blade();
    let mut group = c.benchmark_group("decimate_placement");

    // Worth watching under the error volume: placement is no longer only a question of surface fit,
    // since where the merged vertex lands changes how far the two stars are from each other and so
    // how much of the tolerance a collapse spends.
    for (name, placement) in [
        ("optimal", Placement::Optimal),
        ("midpoint", Placement::Midpoint),
        ("endpoint", Placement::Endpoint),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                let mut he = HalfEdgeMesh3::try_from(black_box(&mesh)).unwrap();
                he.decimate_guaranteed(&DecimateOpts::to_tolerance(0.05).with_placement(placement))
                    .unwrap()
            })
        });
    }

    group.finish();
}

/// The ingest and extraction either side of a decimation run, which a caller pays once per
/// session. Worth keeping visible so it is clear how much of an end-to-end time is conversion.
fn half_edge_round_trip(c: &mut Criterion) {
    let mesh = engine_blade();

    c.bench_function("half_edge_round_trip", |b| {
        b.iter(|| {
            let he = HalfEdgeMesh3::try_from(black_box(&mesh)).unwrap();
            Mesh3::try_from(&he).unwrap()
        })
    });
}

criterion_group!(
    benches,
    decimate_tolerance,
    decimate_placement,
    half_edge_round_trip,
);
criterion_main!(benches);
