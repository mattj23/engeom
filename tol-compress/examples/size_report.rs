//! Prints what the format currently costs, per case, as a markdown table.
//!
//! Every increment after the bit-packing foundation is an optimization, and this is the artifact
//! that says whether one was worth doing. Rerun it after each and paste the delta into the commit
//! message.
//!
//! ```text
//! cargo run -r -p tol-compress --features corpus --example size_report
//! ```
//!
//! Two baselines are reported:
//!
//! - `vs old` is the previous whole-byte scheme, reproduced here from its width rules and its
//!   header layout rather than by depending on `engeom`. Those fully determine its output, so the
//!   comparison is exact. Headers are counted on both sides, which matters on small inputs: the
//!   old header carried a 56 byte quaternion that the current one replaces with a single flag byte.
//! - `vs f32` is uncompressed single-precision storage with `u32` indices, which is what the source
//!   PLY files actually cost. It is the absolute anchor.
//!
//! The `file` column is the real container output, so `over` is what the framing costs on top of
//! the raw blocks. Cases with faces go through the mesh container, point-only cases through the
//! cloud container.
//!
//! Meshes are measured **as the container writes them**, which since increment 4 means with their
//! vertices and faces renumbered by `reorder`. The `mode` column says which index coding won: `hw`
//! for high-water-mark deltas, `abs` for plain absolute indices, which is what a block too small to
//! pay for the block coder's overheads falls back to.

#[path = "../tests/support/ply.rs"]
mod ply;

use tol_compress::corpus::{self, Case};
use tol_compress::{Cloud3, Mesh3, cloud, mesh, reorder, write_indices, write_points};

fn main() {
    println!("# tol-compress size report\n");

    let mut rows: Vec<(String, Case)> = corpus::all()
        .into_iter()
        .map(|c| (c.name.to_string(), c))
        .collect();
    rows.extend(fixture_rows());

    println!("## Corpus at nominal tolerance\n");
    print_table(&rows);

    println!("\n## Tolerance sweep\n");
    println!("The same surface across four orders of magnitude of tolerance.\n");
    let sweep: Vec<(String, Case)> = [1e-1, 1e-2, 1e-3, 1e-4]
        .iter()
        .map(|&tol| {
            (
                format!("smooth_surface @ {tol:.0e}"),
                Case {
                    tol,
                    ..corpus::smooth_surface()
                },
            )
        })
        .collect();
    print_table(&sweep);

    println!("\n## Noise sweep\n");
    println!(
        "Coordinate prediction was measured and dropped, so noise costs nothing in the points \
         block and these rows stay identical. Noise does not disturb connectivity either, so the \
         index blocks match too.\n"
    );
    let noise: Vec<(String, Case)> = [0.0, 0.005, 0.05, 0.5]
        .iter()
        .map(|&s| {
            let case = if s == 0.0 {
                corpus::smooth_surface()
            } else {
                corpus::noisy_surface(s)
            };
            (format!("noise sigma {s}"), case)
        })
        .collect();
    print_table(&noise);
}

fn print_table(rows: &[(String, Case)]) {
    println!(
        "| case | verts | faces | tol | coords | indices | mode | file | over | B/vert | vs old | vs f32 |"
    );
    println!("| --- | --: | --: | --: | --: | --: | :-: | --: | --: | --: | --: | --: |");

    for (label, case) in rows {
        let Some(m) = measure(case) else {
            println!(
                "| {label} | {} | {} | {:.0e} | | | | | | | not representable |",
                case.points.len(),
                case.faces.len(),
                case.tol
            );
            continue;
        };

        let per_vertex = if case.points.is_empty() {
            0.0
        } else {
            m.file as f64 / case.points.len() as f64
        };

        println!(
            "| {label} | {} | {} | {:.0e} | {} | {} | {} | {} | {} | {:.2} | {} | {} |",
            case.points.len(),
            case.faces.len(),
            case.tol,
            m.coords,
            m.indices,
            m.mode,
            m.file,
            m.file - m.total(),
            per_vertex,
            percent(m.file, m.old_total()),
            percent(m.file, m.f32_total()),
        );
    }
}

/// Relative change against a baseline, or a dash where the baseline is meaningless.
fn percent(actual: usize, baseline: usize) -> String {
    if baseline == 0 {
        "-".to_string()
    } else {
        format!("{:+.1}%", 100.0 * (actual as f64 / baseline as f64 - 1.0))
    }
}

struct Measurement {
    coords: usize,
    indices: usize,
    /// Which index coding won, or a dash where there are no indices to code.
    mode: &'static str,
    /// The real container output, framing included.
    file: usize,
    old_coords: usize,
    old_indices: usize,
    f32_total: usize,
}

impl Measurement {
    fn total(&self) -> usize {
        self.coords + self.indices
    }

    fn old_total(&self) -> usize {
        self.old_coords + self.old_indices
    }

    fn f32_total(&self) -> usize {
        self.f32_total
    }
}

fn measure(case: &Case) -> Option<Measurement> {
    // Measured the way the container writes it, so the block columns and the file column describe
    // the same bytes.
    let (points, faces) = if case.faces.is_empty() {
        (case.points.clone(), case.faces.clone())
    } else {
        let plan = reorder::optimize(&case.faces, case.points.len()).ok()?;
        (
            reorder::permute(&case.points, &plan.vertex_order),
            plan.faces,
        )
    };

    let mut coord_buf = Vec::new();
    write_points(&mut coord_buf, &points, case.tol).ok()?;

    let mut index_buf = Vec::new();
    let limit = case.points.len() as u32;
    write_indices(&mut index_buf, &faces, limit).ok()?;

    // The mode byte sits after the arity and the simplex count.
    let mode = match index_buf.get(5) {
        Some(&tol_compress::indices::MODE_HIGH_WATER) => "hw",
        Some(_) => "abs",
        None => "-",
    };

    // Point-only cases have no business paying for an index block, so they go through the cloud
    // container and the framing comparison stays honest.
    let mut file_buf = Vec::new();
    if case.faces.is_empty() {
        cloud::write_one_to(&mut file_buf, &Cloud3::new(case.points.clone()), case.tol).ok()?;
    } else {
        mesh::write_one_to(
            &mut file_buf,
            &Mesh3::new(case.points.clone(), case.faces.clone()),
            case.tol,
        )
        .ok()?;
    }

    let extents = case.extents();
    let axis_tol = case.tol / 3f64.sqrt();

    // The old point block: partition count, point count, a 7 by f64 isometry, a 6 by f64 box, and
    // one width byte per axis. It had no empty-input case, so the header was always written.
    let old_coords = OLD_POINT_HEADER
        + case.points.len()
            * extents
                .iter()
                .map(|&e| old_bytes_for_tol(e, axis_tol) as usize)
                .sum::<usize>();

    // The old index block: just a count, then whole-byte indices.
    let old_indices = OLD_INDEX_HEADER + case.faces.len() * 3 * old_bytes_for_count(limit) as usize;

    Some(Measurement {
        coords: coord_buf.len(),
        indices: index_buf.len(),
        mode,
        file: file_buf.len(),
        old_coords,
        old_indices,
        f32_total: case.points.len() * 3 * 4 + case.faces.len() * 3 * 4,
    })
}

const OLD_POINT_HEADER: usize = 4 + 4 + 7 * 8 + 6 * 8 + 3;
const OLD_INDEX_HEADER: usize = 4;

/// The two real meshes, wrapped as corpus cases so they share the same reporting path.
fn fixture_rows() -> Vec<(String, Case)> {
    let bunny = ply::load("bunny.ply");
    let scan = ply::load("scan-chunk.ply");

    vec![
        (
            "bunny.ply".to_string(),
            Case {
                name: "bunny.ply",
                stresses: "real scanned geometry, smooth and low noise; units m",
                points: bunny.points,
                faces: bunny.faces,
                tol: 1e-5,
            },
        ),
        (
            "scan-chunk.ply".to_string(),
            Case {
                name: "scan-chunk.ply",
                stresses: "real scanner output with real measurement noise; units mm",
                points: scan.points,
                faces: scan.faces,
                tol: 1e-3,
            },
        ),
    ]
}

fn old_bytes_for_tol(range: f64, tol: f64) -> u32 {
    let mut bytes = 2u32;
    while range / 2f64.powi((bytes as i32) * 8) > tol {
        bytes += 1;
    }
    bytes
}

fn old_bytes_for_count(total_items: u32) -> u32 {
    let mut bytes = 2u32;
    while (1u64 << (bytes * 8)) <= u64::from(total_items) {
        bytes += 1;
    }
    bytes
}
