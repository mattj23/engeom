use clap::{Parser, Subcommand};
use engeom::Result;
use engeom::io::{load_ply_mesh, read_mesh_stl, u_bytes_to_mesh, u_mesh_to_bytes};
use std::fs;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "engeom-utils")]
#[command(about = "Engeom geometry utilities")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert an STL or PLY mesh file to the micro mesh binary format
    ToUmesh {
        /// Input STL or PLY file
        input: PathBuf,
        /// Output binary file
        output: PathBuf,
    },
}

fn cmd_to_umesh(input: &PathBuf, output: &PathBuf) -> Result<()> {
    let ext = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    let (vertices, triangles) = match ext.as_deref() {
        Some("stl") => {
            let mesh = read_mesh_stl(input, true, true)?;
            let verts = mesh.vertices().to_vec();
            let tris = mesh.faces().to_vec();
            (verts, tris)
        }
        Some("ply") => {
            let mesh = load_ply_mesh(input)?;
            let verts = mesh.vertices().to_vec();
            let tris = mesh.faces().to_vec();
            (verts, tris)
        }
        _ => return Err("Input file must have a .stl or .ply extension".into()),
    };

    let bytes = u_mesh_to_bytes(&vertices, &triangles)?;

    // Load it back again and check the deviation
    let (u_vert, _) = u_bytes_to_mesh(&bytes)?;

    // Verify that the number of vertices is the same
    if (vertices.len() as u32) != u_vert.len() as u32 {
        return Err(
            "Number of vertices in the original mesh and the micro mesh are not the same".into(),
        );
    }

    let deviations = vertices
        .iter()
        .zip(u_vert.iter())
        .map(|(a, b)| (a - b).norm())
        .collect::<Vec<_>>();

    // Compute min, max, and average deviation
    let max_dev = deviations
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let avg_dev = deviations.iter().sum::<f64>() / deviations.len() as f64;

    println!("Saved micro mesh to {}", output.to_str().unwrap());
    println!(" > {} vertices, {} faces", vertices.len(), triangles.len());
    println!(" > Max deviation: {}", max_dev);
    println!(" > Average deviation: {}", avg_dev);

    fs::write(output, bytes)?;
    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::ToUmesh { input, output } => cmd_to_umesh(input, output),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
