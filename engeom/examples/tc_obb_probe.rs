//! Measure what an oriented box would be worth to the tcmesh point block, without changing anything.
//!
//! The format's partitions are axis-aligned. A scanned surface is a hollow shell, so a patch of it
//! sitting obliquely to the world axes fills its axis-aligned box almost entirely, while the same
//! patch in its own tangent frame is thin on one axis. This prices that difference on real data
//! before any of it is built.
//!
//! ```text
//! cargo run -r -p engeom --features ply --example tc_obb_probe -- <path> <tolerance> [label]
//! ```
//!
//! `.ply`, `.g3d` and `.tcmesh` are read by extension. What the number leans on, in both
//! directions, so it can be read for what it is:
//!
//! - Optimistic: the oriented widths come from an exact eigendecomposition, where a real encoder
//!   would store a quantized rotation and get a slightly fatter box. A twelve bit component is
//!   about 3e-4 radians, so over a patch a few millimetres across the inflation is far below a bit.
//! - Pessimistic: a covariance fit is not the minimum-volume box, only a cheap and decent proxy, so
//!   a better fit can only improve on what is measured here.
//! - Pessimistic: the partitions priced here are the ones the current planner chose while pricing
//!   axis-aligned boxes. A planner that knew about orientation would cut elsewhere. The coarser and
//!   finer rows bracket how much that is worth, though neither moves *where* a cut sits.
//! - Assumed: an oriented partition is charged the same corner bits as an axis-aligned one, on the
//!   grounds that a centre and three extents cost what two corners cost. Its rotation is charged on
//!   top of that, and the last columns show what happens when each partition may decline it.

use engeom::io::{load_g3d_mesh_data, read_ply_mesh_data, read_tc_mesh_file};
use engeom::na::Matrix3;
use engeom::{MeshData3, Result};
use std::io::BufReader;
use std::path::Path;
use tol_compress::{Effort, bits_for_tol, points::partition_bytes, reorder, segment};

/// The bits a quantized rotation would occupy: two to say which quaternion component is largest,
/// then the other three at twelve bits each.
///
/// It can be this coarse because rotation precision does not enter the tolerance budget. The
/// decoder applies the same quantized rotation the encoder used, and a rotation is an isometry, so
/// a coarse one costs a slightly looser box rather than any accuracy.
const ROTATION_BITS: u64 = 2 + 3 * 12;

fn load(path: &Path) -> Result<MeshData3> {
    match path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .as_deref()
    {
        Some("g3d") => load_g3d_mesh_data(path),
        Some("ply") => read_ply_mesh_data(BufReader::new(std::fs::File::open(path)?)),
        Some("tcmesh") => read_tc_mesh_file(path),
        other => Err(format!("unsupported extension {other:?}").into()),
    }
}

/// The per-axis extents of the tightest box around `points` in the frame that fits them best.
///
/// The frame comes from the eigenvectors of the covariance, which is the standard cheap stand-in
/// for a minimum-volume box and is a good one for a surface patch: the smallest eigenvector is the
/// patch normal, so the third extent collapses to how far the patch bends away from its own plane.
fn oriented_extents(points: &[[f64; 3]]) -> [f64; 3] {
    let n = points.len() as f64;
    let mut mean = [0.0; 3];
    for p in points {
        for i in 0..3 {
            mean[i] += p[i] / n;
        }
    }

    let mut cov = Matrix3::<f64>::zeros();
    for p in points {
        let d = [p[0] - mean[0], p[1] - mean[1], p[2] - mean[2]];
        for i in 0..3 {
            for j in 0..3 {
                cov[(i, j)] += d[i] * d[j] / n;
            }
        }
    }

    let axes = cov.symmetric_eigen().eigenvectors;

    let mut mins = [f64::MAX; 3];
    let mut maxs = [f64::MIN; 3];
    for p in points {
        let d = [p[0] - mean[0], p[1] - mean[1], p[2] - mean[2]];
        for i in 0..3 {
            let q = axes[(0, i)] * d[0] + axes[(1, i)] * d[1] + axes[(2, i)] * d[2];
            mins[i] = mins[i].min(q);
            maxs[i] = maxs[i].max(q);
        }
    }

    std::array::from_fn(|i| maxs[i] - mins[i])
}

/// The even split's widths for a set of extents, which is what the writer uses today.
fn widths_for(extents: &[f64; 3], tol: f64) -> [u8; 3] {
    let axis_tol = tol / 3f64.sqrt();
    std::array::from_fn(|i| bits_for_tol(extents[i], axis_tol).unwrap())
}

fn percent(actual: u64, baseline: u64) -> String {
    if baseline == 0 {
        "-".to_string()
    } else {
        format!("{:+.2}%", 100.0 * (actual as f64 / baseline as f64 - 1.0))
    }
}

/// Price one set of partitions three ways: as written today, oriented, and oriented only where that
/// pays for its own rotation.
struct Priced {
    parts: usize,
    flat: u64,
    oriented: u64,
    chosen: u64,
    turned: usize,
}

fn price(
    points: &[[f64; 3]],
    spans: &[(usize, usize)],
    anchor_widths: &[u8; 3],
    tol: f64,
) -> Priced {
    let mut out = Priced {
        parts: spans.len(),
        flat: 0,
        oriented: 0,
        chosen: 0,
        turned: 0,
    };

    for &(start, end) in spans {
        let run = &points[start..end];
        let count = end - start;

        let mut mins = run[0];
        let mut maxs = run[0];
        for p in run {
            for i in 0..3 {
                mins[i] = mins[i].min(p[i]);
                maxs[i] = maxs[i].max(p[i]);
            }
        }
        let flat_extents: [f64; 3] = std::array::from_fn(|i| maxs[i] - mins[i]);

        let flat = partition_bytes(count, &widths_for(&flat_extents, tol), anchor_widths, None);
        let obb = widths_for(&oriented_extents(run), tol);

        // The same accounting `partition_bytes` does, with the rotation added to the bitstream.
        let corners: u64 = anchor_widths.iter().map(|&w| 2 * u64::from(w)).sum();
        let per_point: u64 = obb.iter().map(|&w| u64::from(w)).sum();
        let oriented = 4 + 1 + 3 + (corners + ROTATION_BITS + count as u64 * per_point).div_ceil(8);

        out.flat += flat;
        out.oriented += oriented;
        out.chosen += flat.min(oriented);
        if oriented < flat {
            out.turned += 1;
        }
    }

    out
}

/// Cut each span into `pieces` contiguous chunks, which is what a planner that knew orientation was
/// available would be free to do and the current one is not.
fn subdivide(spans: &[(usize, usize)], pieces: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for &(start, end) in spans {
        let n = end - start;
        let per = n.div_ceil(pieces).max(1);
        let mut at = start;
        while at < end {
            out.push((at, (at + per).min(end)));
            at += per;
        }
    }
    out
}

/// Glue `runs` adjacent spans into one, which is the direction an orientation-aware planner is more
/// likely to move: cheaper points make a fixed header relatively more expensive, so the runs it
/// wants should be longer than the ones an axis-aligned planner settled on.
fn merge(spans: &[(usize, usize)], runs: usize) -> Vec<(usize, usize)> {
    spans
        .chunks(runs)
        .map(|c| (c[0].0, c[c.len() - 1].1))
        .collect()
}

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("path to a mesh file");
    let tol: f64 = args.next().expect("tolerance").parse().expect("a number");
    let path = Path::new(&path);
    let label = args
        .next()
        .unwrap_or_else(|| path.file_name().unwrap().to_string_lossy().into_owned());

    let data = load(path)?;
    let points: Vec<[f64; 3]> = data.points().iter().map(|p| [p.x, p.y, p.z]).collect();

    // The writer's own vertex order, so the runs measured here are the runs it would produce.
    let order = reorder::optimize(data.faces(), points.len())?;
    let points = reorder::permute(&points, &order.vertex_order);

    let plan = segment::plan(&points, tol, Effort::Balanced)?;
    let spans: Vec<(usize, usize)> = plan.runs.iter().map(|r| (r.start, r.end)).collect();

    println!("## {label}: {} points at tol {tol:.0e}\n", points.len());
    println!(
        "| cuts | parts | pts/part | flat | oriented | vs flat | best of each | vs flat | turned |"
    );
    println!("| --- | --: | --: | --: | --: | --: | --: | --: | --: |");

    // Negative is coarser, positive finer, zero the planner's own cuts.
    for step in [-8i32, -4, -2, 0, 2, 4, 8, 16] {
        let (name, spans) = match step {
            0 => ("as planned".to_string(), spans.clone()),
            s if s < 0 => (format!("x{}", -s), merge(&spans, (-s) as usize)),
            s => (format!("/{s}"), subdivide(&spans, s as usize)),
        };
        let p = price(&points, &spans, &plan.anchor_widths, tol);

        println!(
            "| {} | {} | {:.0} | {} | {} | {} | {} | {} | {:.0}% |",
            name,
            p.parts,
            points.len() as f64 / p.parts as f64,
            p.flat,
            p.oriented,
            percent(p.oriented, p.flat),
            p.chosen,
            percent(p.chosen, p.flat),
            100.0 * p.turned as f64 / p.parts as f64,
        );
    }

    Ok(())
}
