//! Encode and decode throughput.
//!
//! Size is the headline number for this crate and the size report covers it. This measures the
//! other axis, because a format that halves a file but takes ten times as long to read is a
//! different product. Reported as elements per second so cases of different sizes compare directly.
//!
//! ```text
//! cargo bench -p tol-compress --features corpus
//! cargo bench -p tol-compress --features corpus -- --save-baseline before
//! ```
//!
//! The `bits` group measures the primitive everything else sits on. Increment 5 adds per-block
//! adaptive widths on that same foundation, so a regression there would show up as a slowdown
//! everywhere at once and is worth being able to see in isolation.

#[path = "../tests/support/ply.rs"]
mod ply;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::io::Cursor;
use tol_compress::bits::{BitReader, BitWriter};
use tol_compress::corpus;
use tol_compress::testgen::Rng;
use tol_compress::{read_indices, read_points, write_indices, write_points};

struct Workload {
    name: String,
    points: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    tol: f64,
}

/// The cases benchmarked, chosen to span small real data, mid-size generated data, and a set large
/// enough that per-call overhead disappears.
fn workloads() -> Vec<Workload> {
    let smooth = corpus::smooth_surface();
    let clusters = corpus::distant_clusters();
    let scan = ply::load("scan-chunk.ply");

    let mut rng = Rng::new(4242);

    vec![
        Workload {
            name: "scan-chunk".to_string(),
            points: scan.points,
            faces: scan.faces,
            tol: 1e-3,
        },
        Workload {
            name: "smooth_surface".to_string(),
            points: smooth.points,
            faces: smooth.faces,
            tol: smooth.tol,
        },
        Workload {
            name: "distant_clusters".to_string(),
            points: clusters.points,
            faces: clusters.faces,
            tol: clusters.tol,
        },
        Workload {
            name: "large_random".to_string(),
            points: rng.points(500_000, -100.0, 100.0),
            faces: Vec::new(),
            tol: 1e-3,
        },
    ]
}

fn points_encode(c: &mut Criterion) {
    let mut group = c.benchmark_group("points/encode");

    for w in workloads() {
        let (points, tol) = (w.points, w.tol);
        group.throughput(Throughput::Elements(points.len() as u64));
        group.bench_function(&w.name, |b| {
            let mut buf = Vec::with_capacity(points.len() * 8);
            b.iter(|| {
                buf.clear();
                write_points(&mut buf, black_box(&points), tol).unwrap();
                black_box(buf.len())
            })
        });
    }

    group.finish();
}

fn points_decode(c: &mut Criterion) {
    let mut group = c.benchmark_group("points/decode");

    for w in workloads() {
        let (points, tol) = (w.points, w.tol);
        let mut buf = Vec::new();
        write_points(&mut buf, &points, tol).unwrap();

        group.throughput(Throughput::Elements(points.len() as u64));
        group.bench_function(&w.name, |b| {
            b.iter(|| {
                let mut cursor = Cursor::new(black_box(&buf));
                let out: Vec<[f64; 3]> = read_points(&mut cursor).unwrap();
                black_box(out.len())
            })
        });
    }

    group.finish();
}

fn indices_encode(c: &mut Criterion) {
    let mut group = c.benchmark_group("indices/encode");

    for w in workloads() {
        let faces = w.faces;
        if faces.is_empty() {
            continue;
        }
        let limit = w.points.len() as u32;

        group.throughput(Throughput::Elements(faces.len() as u64));
        group.bench_function(&w.name, |b| {
            let mut buf = Vec::with_capacity(faces.len() * 8);
            b.iter(|| {
                buf.clear();
                write_indices(&mut buf, black_box(&faces), limit).unwrap();
                black_box(buf.len())
            })
        });
    }

    group.finish();
}

fn indices_decode(c: &mut Criterion) {
    let mut group = c.benchmark_group("indices/decode");

    for w in workloads() {
        let faces = w.faces;
        if faces.is_empty() {
            continue;
        }
        let limit = w.points.len() as u32;
        let mut buf = Vec::new();
        write_indices(&mut buf, &faces, limit).unwrap();

        group.throughput(Throughput::Elements(faces.len() as u64));
        group.bench_function(&w.name, |b| {
            b.iter(|| {
                let mut cursor = Cursor::new(black_box(&buf));
                let out: Vec<[u32; 3]> = read_indices(&mut cursor, limit).unwrap();
                black_box(out.len())
            })
        });
    }

    group.finish();
}

/// The bit packer in isolation, at a byte-aligned width and at two that are not.
///
/// Widths that straddle byte boundaries are the ones this crate actually uses, so the comparison
/// against 16 shows what the packing costs relative to plain byte writes.
fn bits(c: &mut Criterion) {
    const COUNT: usize = 1_000_000;
    let values: Vec<u64> = (0..COUNT as u64)
        .map(|i| i.wrapping_mul(2_654_435_761))
        .collect();

    let mut group = c.benchmark_group("bits");
    group.throughput(Throughput::Elements(COUNT as u64));

    for width in [16u8, 17, 20] {
        let mask = (1u64 << width) - 1;
        let masked: Vec<u64> = values.iter().map(|v| v & mask).collect();

        group.bench_function(format!("write/{width}"), |b| {
            let mut buf = Vec::with_capacity(COUNT * 8);
            b.iter(|| {
                buf.clear();
                let mut w = BitWriter::new(&mut buf);
                for &v in black_box(&masked) {
                    w.write_bits(v, width).unwrap();
                }
                w.finish().unwrap();
                drop(w);
                black_box(buf.len())
            })
        });

        let mut buf = Vec::new();
        {
            let mut w = BitWriter::new(&mut buf);
            for &v in &masked {
                w.write_bits(v, width).unwrap();
            }
            w.finish().unwrap();
        }

        group.bench_function(format!("read/{width}"), |b| {
            b.iter(|| {
                let mut r = BitReader::new(Cursor::new(black_box(&buf)));
                let mut acc = 0u64;
                for _ in 0..COUNT {
                    acc = acc.wrapping_add(r.read_bits(width).unwrap());
                }
                black_box(acc)
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    points_encode,
    points_decode,
    indices_encode,
    indices_decode,
    bits,
);
criterion_main!(benches);
