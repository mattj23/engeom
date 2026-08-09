//! Mesh serialization for the practical tol-compression format.
//!
//! This is a thin adapter over the [`tol_compress`] crate: it converts between engeom's
//! [`MeshData3`] and that crate's mesh container, and adds nothing of its own. The `tol` parameter
//! throughout is the maximum acceptable round-trip position error for any vertex, in the same units
//! as the coordinates.
//!
//! The recommended format extension is `.tcmesh`.
//!
//! # Vertices are renumbered
//!
//! Writing reorders vertices and faces so that the connectivity compresses, which is where most of
//! the format's advantage over a naive quantized mesh comes from. A mesh read back therefore
//! describes the **same surface** but not necessarily with the same indices, and code that holds
//! per-vertex data outside the file cannot assume its indices still line up. The ordering is derived
//! entirely from connectivity, so nothing about it is stored and a caller who needs to follow it can
//! compute the same plan up front; see the `tol_compress::mesh` module documentation for that
//! recipe.
//!
//! # Attributes
//!
//! The format stores geometry only. Writing a mesh that carries attributes is an **error** naming
//! them, rather than a silent drop or an opt-in discard, because the format will grow attribute
//! support later and a file written today should not have quietly lost data by then.

use crate::{MeshData3, Point3, Result};
use std::io::{Read, Write};
use std::path::Path;
use tol_compress::{Mesh3 as TcMesh, mesh as tc};

/// Refuse a mesh carrying attributes, naming what would be lost.
///
/// Unlike the other geometry-only formats in this module there is no `allow_attribute_loss` escape
/// hatch, so the message does not offer one.
fn refuse_attributes(mesh: &MeshData3) -> Result<()> {
    if mesh.attrs().is_empty() {
        return Ok(());
    }

    let mut lost = mesh.attrs().point_attr_labels();
    lost.extend(mesh.attrs().face_attr_labels());

    Err(format!(
        "The tcmesh format cannot yet represent the attributes on this mesh ({}), so writing it \
         would discard them. Remove them from the mesh explicitly if that is what you want.",
        lost.join(", ")
    )
    .into())
}

/// A mesh as the storable container, refusing anything the format cannot represent.
fn to_container(mesh: &MeshData3) -> Result<TcMesh> {
    refuse_attributes(mesh)?;
    let points = mesh.points().iter().map(|p| [p.x, p.y, p.z]).collect();
    Ok(TcMesh::new(points, mesh.faces().to_vec()))
}

/// The inverse of [`to_container`].
fn from_container(item: TcMesh) -> Result<MeshData3> {
    let points = item
        .points
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect();
    MeshData3::new(points, item.faces)
}

/// Serialize a mesh into the tcmesh format, writing to any [`Write`] sink.
///
/// Vertex positions are quantized to the narrowest bit width per axis that keeps every one of them
/// within `tol` of where it started, and the connectivity is stored exactly. Smaller `tol` values
/// produce more accurate output at the cost of more bytes per vertex.
///
/// Vertices are renumbered; see the module documentation.
///
/// # Errors
///
/// Returns an error naming the attributes if the mesh carries any, since the format stores geometry
/// only.
///
/// Use [`write_tc_mesh_file`] for the common case of writing directly to a file path.
pub fn write_tc_mesh_to<W: Write>(writer: &mut W, mesh: &MeshData3, tol: f64) -> Result<()> {
    tc::write_one_to(writer, &to_container(mesh)?, tol)?;
    Ok(())
}

/// Deserialize a mesh from a tcmesh-format byte stream.
///
/// The recovered vertex positions are guaranteed to be within the tolerance that was supplied at
/// write time, and the connectivity is exact.
///
/// Use [`read_tc_mesh_file`] for the common case of reading from a file path.
pub fn read_tc_mesh_from<R: Read>(reader: &mut R) -> Result<MeshData3> {
    from_container(tc::read_one_from(reader)?)
}

/// Write a mesh to a tcmesh file at the given path. See [`write_tc_mesh_to`] for format details.
pub fn write_tc_mesh_file(path: &Path, mesh: &MeshData3, tol: f64) -> Result<()> {
    tc::write_one_file(path, &to_container(mesh)?, tol)?;
    Ok(())
}

/// Read a mesh from a tcmesh file at the given path. See [`read_tc_mesh_from`] for format details.
pub fn read_tc_mesh_file(path: &Path) -> Result<MeshData3> {
    from_container(tc::read_one_file(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::{stanford_bun_2, stanford_bun_4};
    use std::io::Cursor;
    use tol_compress::reorder;

    /// Assert the recovered mesh describes the same surface as the original.
    ///
    /// Array equality is the wrong assertion now that writing renumbers vertices, and quietly
    /// weakening it to a count comparison would be the way this stops catching anything. Instead the
    /// same connectivity-derived plan the encoder used is recomputed here, and every face corner is
    /// checked against the position it came from.
    fn check_round_trip(original: &MeshData3, recovered: &MeshData3, tol: f64) {
        assert_eq!(original.points().len(), recovered.points().len());
        assert_eq!(original.faces().len(), recovered.faces().len());

        let plan = reorder::optimize(original.faces(), original.points().len()).unwrap();

        for (new, &old) in plan.face_order.iter().enumerate() {
            for corner in 0..3 {
                let got = recovered.points()[recovered.faces()[new][corner] as usize];
                let want = original.points()[original.faces()[old as usize][corner] as usize];
                let d = (got - want).norm();
                assert!(d <= tol, "face {new} corner {corner} moved {d}, tol {tol}");
            }
        }

        for (i, &old) in plan.vertex_order.iter().enumerate() {
            let d = (recovered.points()[i] - original.points()[old as usize]).norm();
            assert!(d <= tol, "vertex {i} moved {d}, tol {tol}");
        }
    }

    #[test]
    fn round_trip_bytes() {
        let mesh = stanford_bun_2().to_data();
        let tol = 1e-4;
        let mut buf = Vec::new();
        write_tc_mesh_to(&mut buf, &mesh, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_mesh_from(&mut cursor).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    #[test]
    fn round_trip_file() {
        let mesh = stanford_bun_2().to_data();
        let tol = 1e-4;
        let path = std::env::temp_dir().join("stanford_bun_2_round_trip.tcmesh");
        write_tc_mesh_file(&path, &mesh, tol).unwrap();
        let recovered = read_tc_mesh_file(&path).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    #[test]
    fn round_trip_high_res() {
        // stanford_bun_4 has more vertices, exercising a larger index width
        let mesh = stanford_bun_4().to_data();
        let tol = 1e-4;
        let mut buf = Vec::new();
        write_tc_mesh_to(&mut buf, &mesh, tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered = read_tc_mesh_from(&mut cursor).unwrap();
        check_round_trip(&mesh, &recovered, tol);
    }

    /// The format stores geometry only, and there is no longer a flag that accepts the loss, so an
    /// attributed mesh has to be refused outright.
    #[test]
    fn writing_refuses_a_mesh_carrying_attributes() {
        let mut mesh = stanford_bun_2().to_data();
        mesh.set_face_labels(Some(vec![1; mesh.face_count()]))
            .unwrap();

        let mut buf = Vec::new();
        let err = match write_tc_mesh_to(&mut buf, &mesh, 1e-4) {
            Err(e) => e.to_string(),
            Ok(()) => panic!("an attributed mesh must not be silently written"),
        };
        assert!(
            err.contains("face_labels"),
            "the error must name what would be lost, got: {err}"
        );
    }

    #[test]
    fn invalid_magic_rejected() {
        let bad = b"NOPE\x00\x01\x00\x00\x00";
        let mut cursor = Cursor::new(bad);
        assert!(read_tc_mesh_from(&mut cursor).is_err());
    }

    /// Renumbering is the reason this format beats a naive quantized mesh, so it is worth pinning
    /// that it actually happens rather than trusting the default never changes.
    #[test]
    fn writing_renumbers_a_badly_ordered_mesh() {
        // The embedded fixtures are themselves tcmesh files, so they arrive already in a good order
        // and would make a poor subject. Reversing the faces gives something that certainly is not.
        let mesh = stanford_bun_2().to_data();
        let mut faces = mesh.faces().to_vec();
        faces.reverse();
        let scrambled = MeshData3::new(mesh.points().to_vec(), faces).unwrap();

        let mut buf = Vec::new();
        write_tc_mesh_to(&mut buf, &scrambled, 1e-4).unwrap();
        let recovered = read_tc_mesh_from(&mut Cursor::new(&buf)).unwrap();

        assert_ne!(
            recovered.faces(),
            scrambled.faces(),
            "a reversed face order cannot also be a good one"
        );
    }

    /// Writing a mesh that is already in the order the encoder would choose has to reproduce that
    /// order, and rewriting a file has to come back to somewhere it has already been rather than
    /// walking away from the geometry a tolerance at a time.
    ///
    /// It is not a fixed point and cannot be. The encoder decides where to partition a point block
    /// from where the points are, the quantizer moves them, and a mesh read back out therefore
    /// plans slightly differently than the one that went in. Some inputs settle on a single file
    /// and some alternate between two.
    ///
    /// What matters is that the orbit is short and closed, because a repeat means the decoded
    /// points repeated exactly. Error cannot accumulate through a pipeline that stores and restores.
    #[test]
    fn rewriting_a_mesh_reaches_a_cycle() {
        let mesh = stanford_bun_2().to_data();
        let tol = 3e-6;

        let mut current = mesh.clone();
        let mut seen: Vec<Vec<u8>> = Vec::new();
        let mut closed = None;

        for round in 0..8 {
            let mut buf = Vec::new();
            write_tc_mesh_to(&mut buf, &current, tol).unwrap();
            if seen.contains(&buf) {
                closed = Some(round);
                break;
            }
            seen.push(buf.clone());
            current = read_tc_mesh_from(&mut Cursor::new(&buf)).unwrap();

            assert_eq!(
                current.faces(),
                mesh.faces(),
                "connectivity must not move, round {round}"
            );
        }

        assert!(
            matches!(closed, Some(r) if r <= 4),
            "rewriting must close its orbit, got {closed:?}"
        );
    }
}
