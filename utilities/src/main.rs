use clap::{Parser, Subcommand};
use engeom::io::{read_tc_mesh_file, write_tc_mesh_file};
use engeom::{Mesh3, MeshData3, Result};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "engeom-utils")]
#[command(about = "Engeom geometry utilities")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert an STL or PLY mesh file to the tcmesh format
    ToTcmesh {
        /// Input STL or PLY file
        input: PathBuf,
        /// Output .tcmesh file (defaults to input path with .tcmesh extension)
        output: Option<PathBuf>,
        /// Round-trip position tolerance in model units
        #[arg(long, default_value_t = 1e-6)]
        tol: f64,
    },
}

/// Load a mesh from whichever of the supported input formats the extension names.
///
/// The STL path goes through `Mesh3` only to reach the point-merging and degenerate-triangle
/// cleanup, which lives on the accelerated constructor because it renumbers points.
fn load_input_mesh(input: &Path, ext: Option<&str>) -> Result<MeshData3> {
    match ext {
        Some("stl") => Ok(Mesh3::load_stl(input, false, true, true)?.into_data()),
        Some("ply") => MeshData3::load_ply(input),
        _ => Err("Input file must have a .stl or .ply extension".into()),
    }
}

fn cmd_to_tcmesh(input: &Path, output: Option<&Path>, tol: f64) -> Result<()> {
    let output = output
        .map(PathBuf::from)
        .unwrap_or_else(|| input.with_extension("tcmesh"));
    let ext = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    let mesh = load_input_mesh(input, ext.as_deref())?;

    write_tc_mesh_file(&output, &mesh, tol)?;

    let recovered = read_tc_mesh_file(&output)?;
    if mesh.points().len() != recovered.points().len() {
        return Err("Vertex count mismatch after round-trip".into());
    }

    // Writing renumbers vertices, so the two point arrays cannot be compared position by position.
    // The ordering is derived from connectivity, so recomputing the same plan says which original
    // vertex each recovered one came from. Zipping the arrays directly would compare unrelated
    // points and report a "deviation" the size of the part.
    let plan = engeom::tol_compress::reorder::optimize(mesh.faces(), mesh.points().len())?;
    let deviations = plan
        .vertex_order
        .iter()
        .enumerate()
        .map(|(new, &old)| (mesh.points()[old as usize] - recovered.points()[new]).norm())
        .collect::<Vec<_>>();

    let max_dev = deviations
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let avg_dev = deviations.iter().sum::<f64>() / deviations.len() as f64;

    println!("Saved tcmesh to {}", output.display());
    println!(
        " > {} vertices, {} faces",
        mesh.points().len(),
        mesh.faces().len()
    );
    println!(" > Tolerance: {tol}");
    println!(" > Max deviation: {max_dev}");
    println!(" > Average deviation: {avg_dev}");

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::ToTcmesh { input, output, tol } => cmd_to_tcmesh(input, output.as_deref(), *tol),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
