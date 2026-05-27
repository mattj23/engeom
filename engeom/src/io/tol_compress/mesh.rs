//! Mesh serialization for the practical tol-compression mesh format. This uses the algorithms
//! in the `engeom::io::tol_compress::core` module and its submodules to dynamically adjust the 
//! byte width and store mesh information in an efficient binary representation.
//!
//! The recommended format extension is `.tcmesh`

use crate::{Mesh, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use crate::io::tol_compress::core::{read_indices, read_tc_points3, write_indices, write_tc_points3};

const MAGIC: &[u8; 6] = b"TCMESH";

pub fn write_tc_mesh_to<W: Write>(writer: &mut W, mesh: &Mesh, tol: f64) -> Result<()> {
    writer.write_all(MAGIC)?;
    write_tc_points3(writer, mesh.vertices(), tol)?;
    write_indices(writer, mesh.faces(), mesh.vertices().len() as u32)?;
    Ok(())
}

pub fn read_tc_mesh_from<R: Read>(reader: &mut R) -> Result<Mesh> {
    let mut magic = [0u8; 6];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err("Not a tcmesh file: invalid magic bytes".into());
    }

    let vertices = read_tc_points3(reader)?;
    let faces = read_indices::<_, 3>(reader, vertices.len() as u32)?;
    Ok(Mesh::new(vertices, faces, false))
}

pub fn write_tc_mesh_file(path: &Path, mesh: &Mesh, tol: f64) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_tc_mesh_to(&mut writer, mesh, tol)
}

pub fn read_tc_mesh_file(path: &Path) -> Result<Mesh> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    read_tc_mesh_from(&mut reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::{stanford_bun_2, stanford_bun_4};
    use approx::assert_relative_eq;
    use std::io::Cursor;

    fn check_round_trip(original: &Mesh, recovered: &Mesh, tol: f64) {
        assert_eq!(original.vertices().len(), recovered.vertices().len());
        assert_eq!(original.faces().len(), recovered.faces().len());
        for (a, b) in original.vertices().iter().zip(recovered.vertices().iter()) {
            assert_relative_eq!(a, b, epsilon = tol);
        }
        for (a, b) in original.faces().iter().zip(recovered.faces().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn round_trip_bytes() {
        let mesh = stanford_bun_2();
        let tol = 1e-4;
        let mut buf = Vec::new();
        write_tc_mesh_to(&mut buf, &mesh, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_mesh_from(&mut cursor).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    #[test]
    fn round_trip_file() {
        let mesh = stanford_bun_2();
        let tol = 1e-4;
        let path = std::env::temp_dir().join("stanford_bun_2_round_trip.tcmesh");
        write_tc_mesh_file(&path, &mesh, tol).unwrap();
        let recovered = read_tc_mesh_file(&path).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    #[test]
    fn round_trip_high_res() {
        // stanford_bun_4 has more vertices, exercising a larger index byte width
        let mesh = stanford_bun_4();
        let tol = 1e-4;
        let mut buf = Vec::new();
        write_tc_mesh_to(&mut buf, &mesh, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_mesh_from(&mut cursor).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    #[test]
    fn invalid_magic_rejected() {
        let bad = b"NOPE\x00\x01\x00\x00\x00";
        let mut cursor = Cursor::new(bad);
        assert!(read_tc_mesh_from(&mut cursor).is_err());
    }
}
