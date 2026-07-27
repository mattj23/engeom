use clap::{Parser, Subcommand};
use engeom::io::{
    read_tc_mesh_file, u_bytes_to_mesh_data, u_mesh_data_to_bytes, write_tc_mesh_file,
};
use engeom::{Mesh3, MeshData3, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs;
use std::io::Write;
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
    /// Convert an STL or PLY mesh file to the micro mesh binary format
    ToUmesh {
        /// Input STL or PLY file
        input: PathBuf,
        /// Output binary file
        output: PathBuf,
        /// Compress the output with gzip
        #[arg(long)]
        compress: bool,
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

    write_tc_mesh_file(&output, &mesh, tol, false)?;

    let recovered = read_tc_mesh_file(&output)?;
    if mesh.points().len() != recovered.points().len() {
        return Err("Vertex count mismatch after round-trip".into());
    }

    let deviations = mesh
        .points()
        .iter()
        .zip(recovered.points().iter())
        .map(|(a, b)| (a - b).norm())
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

fn cmd_to_umesh(input: &Path, output: &PathBuf, compress: bool) -> Result<()> {
    let ext = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    let mesh = load_input_mesh(input, ext.as_deref())?;
    let vertices = mesh.points().to_vec();

    let bytes = u_mesh_data_to_bytes(&mesh, false)?;

    // Load it back again and check the deviation
    let u_vert = u_bytes_to_mesh_data(&bytes)?.points().to_vec();

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
    println!(
        " > {} vertices, {} faces",
        vertices.len(),
        mesh.faces().len()
    );
    println!(" > Max deviation: {}", max_dev);
    println!(" > Average deviation: {}", avg_dev);

    if compress {
        let file = fs::File::create(output)?;
        let mut encoder = GzEncoder::new(file, Compression::best());
        encoder.write_all(&bytes)?;
        encoder.finish()?;
    } else {
        fs::write(output, bytes)?;
    }
    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::ToTcmesh { input, output, tol } => cmd_to_tcmesh(input, output.as_deref(), *tol),
        Commands::ToUmesh {
            input,
            output,
            compress,
        } => cmd_to_umesh(input, output, *compress),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
