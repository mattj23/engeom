//! Measure what the tcmesh format costs on a mesh that is too large to keep in the repository.
//!
//! The corpus and the two committed PLY fixtures are small by design, so they say little about how
//! the format behaves on a real unthinned scan. This runs the same accounting as
//! `tol-compress`'s `size_report` example over a file given on the command line, so the numbers can
//! be taken without the file ever entering the tree.
//!
//! ```text
//! cargo run -r -p engeom --features ply --example tc_size_row -- <path> <tolerance> [label]
//! ```
//!
//! `.ply` and `.g3d` are read by extension, so the `ply` feature has to be on. The tolerance is in
//! the same units as the file's coordinates.

use engeom::io::{load_g3d_mesh_data, read_ply_mesh_data};
use engeom::{MeshData3, Result};
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;
use tol_compress::{Effort, Mesh3, mesh, reorder};

fn load(path: &Path) -> Result<MeshData3> {
    match path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .as_deref()
    {
        Some("g3d") => load_g3d_mesh_data(path),
        Some("ply") => read_ply_mesh_data(BufReader::new(std::fs::File::open(path)?)),
        other => Err(format!("unsupported extension {other:?}").into()),
    }
}

fn percent(actual: usize, baseline: usize) -> String {
    if baseline == 0 {
        "-".to_string()
    } else {
        format!("{:+.1}%", 100.0 * (actual as f64 / baseline as f64 - 1.0))
    }
}

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("path to a mesh file");
    let tol: f64 = args.next().expect("tolerance").parse().expect("a number");
    let path = Path::new(&path);
    let label = args
        .next()
        .unwrap_or_else(|| path.file_name().unwrap().to_string_lossy().into_owned());

    let started = Instant::now();
    let data = load(path)?;
    let loaded = started.elapsed();

    let points: Vec<[f64; 3]> = data.points().iter().map(|p| [p.x, p.y, p.z]).collect();
    let faces = data.faces().to_vec();

    let mut mins = points[0];
    let mut maxs = points[0];
    for p in &points {
        for i in 0..3 {
            mins[i] = mins[i].min(p[i]);
            maxs[i] = maxs[i].max(p[i]);
        }
    }
    let extents: Vec<f64> = (0..3).map(|i| maxs[i] - mins[i]).collect();

    // Measured the way the container writes it, so the block figures and the file figure describe
    // the same bytes.
    let plan = reorder::optimize(&faces, points.len())?;
    let ordered = reorder::permute(&points, &plan.vertex_order);

    let mut quick = Vec::new();
    tol_compress::write_points_with(&mut quick, &ordered, tol, Effort::Quick)?;

    let started = Instant::now();
    let mut coords = Vec::new();
    tol_compress::write_points_with(&mut coords, &ordered, tol, Effort::Balanced)?;
    let planned = started.elapsed();
    let parts = u32::from_le_bytes(coords[1..5].try_into().unwrap());

    let mut indices = Vec::new();
    tol_compress::write_indices(&mut indices, &plan.faces, points.len() as u32)?;

    let item = Mesh3::new(points.clone(), faces.clone());
    let started = Instant::now();
    let mut file = Vec::new();
    mesh::write_one_to(&mut file, &item, tol)?;
    let wrote = started.elapsed();

    let f32_total = points.len() * 3 * 4 + faces.len() * 3 * 4;

    println!(
        "| case | verts | faces | tol | extents | coords | parts | vs quick | indices | file | B/vert | vs f32 |"
    );
    println!("| --- | --: | --: | --: | --- | --: | --: | --: | --: | --: | --: | --: |");
    println!(
        "| {label} | {} | {} | {:.0e} | {:.1} x {:.1} x {:.1} | {} | {} | {} | {} | {} | {:.2} | {} |",
        points.len(),
        faces.len(),
        tol,
        extents[0],
        extents[1],
        extents[2],
        coords.len(),
        parts,
        percent(coords.len(), quick.len()),
        indices.len(),
        file.len(),
        file.len() as f64 / points.len() as f64,
        percent(file.len(), f32_total),
    );
    println!(
        "\nload {:.2}s, plan+encode points {:.2}s, whole file {:.2}s",
        loaded.as_secs_f64(),
        planned.as_secs_f64(),
        wrote.as_secs_f64()
    );

    Ok(())
}
