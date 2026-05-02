use clap::{Parser, Subcommand};
use engeom::Result;
use engeom::io::{
    load_ply_mesh, read_mesh_stl, read_tc_mesh_from, u_bytes_to_mesh_data, u_mesh_data_to_bytes,
    write_tc_mesh_to,
};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs;
use std::io::{Cursor, Write};
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
        /// Output .tcmesh file
        output: PathBuf,
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

fn cmd_to_tcmesh(input: &Path, output: &Path, tol: f64) -> Result<()> {
    let ext = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    let mesh = match ext.as_deref() {
        Some("stl") => read_mesh_stl(input, true, true)?,
        Some("ply") => load_ply_mesh(input)?,
        _ => return Err("Input file must have a .stl or .ply extension".into()),
    };

    let mut buf = Vec::new();
    write_tc_mesh_to(&mut buf, &mesh, tol)?;

    let recovered = read_tc_mesh_from(&mut Cursor::new(&buf))?;
    if mesh.vertices().len() != recovered.vertices().len() {
        return Err("Vertex count mismatch after round-trip".into());
    }

    let deviations = mesh
        .vertices()
        .iter()
        .zip(recovered.vertices().iter())
        .map(|(a, b)| (a - b).norm())
        .collect::<Vec<_>>();

    let max_dev = deviations
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let avg_dev = deviations.iter().sum::<f64>() / deviations.len() as f64;

    fs::write(output, &buf)?;

    println!("Saved tcmesh to {}", output.to_str().unwrap());
    println!(" > {} vertices, {} faces", mesh.vertices().len(), mesh.faces().len());
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

    let bytes = u_mesh_data_to_bytes(&vertices, &triangles)?;

    // Load it back again and check the deviation
    let (u_vert, _) = u_bytes_to_mesh_data(&bytes)?;

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
        Commands::ToTcmesh { input, output, tol } => cmd_to_tcmesh(input, output, *tol),
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
