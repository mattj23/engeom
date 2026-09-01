//! The tolerance guarantee, checked against real measured geometry rather than generated shapes.
//!
//! The synthetic corpus covers the awkward cases deliberately. These two fixtures cover the
//! ordinary one: data that came off actual hardware, with whatever noise, tessellation and
//! numerical history it happens to carry.
//!
//! Both fixtures store coordinates as `f32`, so their own resolution is about seven significant
//! digits. Tolerances below that are measuring the source file's rounding rather than ours, which
//! is why the sweeps here stop where they do.

#[path = "support/ply.rs"]
mod ply;

use tol_compress::{read_indices, read_points, write_indices, write_points};

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}

fn check_points(name: &str, points: &[[f64; 3]], tol: f64) -> usize {
    let mut buf = Vec::new();
    write_points(&mut buf, points, tol).unwrap_or_else(|e| panic!("{name} at {tol}: {e}"));

    let mut cursor = std::io::Cursor::new(&buf);
    let recovered: Vec<[f64; 3]> = read_points(&mut cursor).unwrap();

    assert_eq!(
        cursor.position() as usize,
        buf.len(),
        "{name}: decoder left bytes unread"
    );
    assert_eq!(recovered.len(), points.len(), "{name}");

    for (i, (o, r)) in points.iter().zip(recovered.iter()).enumerate() {
        let d = distance(o, r);
        assert!(
            d <= tol,
            "{name} at tol {tol}: vertex {i} recovered {d} away"
        );
    }

    buf.len()
}

#[test]
fn bunny_round_trips_within_tolerance() {
    let mesh = ply::load("bunny.ply");
    assert_eq!(mesh.points.len(), 453);
    assert_eq!(mesh.faces.len(), 948);

    // Units are metres and the model is about 150 mm across.
    for tol in [1e-3, 1e-4, 1e-5, 1e-6] {
        check_points("bunny", &mesh.points, tol);
    }
}

#[test]
fn scan_chunk_round_trips_within_tolerance() {
    let mesh = ply::load("scan-chunk.ply");
    assert_eq!(mesh.points.len(), 8091);

    // Units are millimetres over roughly a 20 mm region, so a micron is a realistic floor.
    for tol in [1e-1, 1e-2, 1e-3, 1e-4] {
        check_points("scan-chunk", &mesh.points, tol);
    }
}

#[test]
fn real_connectivity_round_trips_exactly() {
    for name in ["bunny.ply", "scan-chunk.ply"] {
        let mesh = ply::load(name);
        let limit = mesh.points.len() as u32;

        let mut buf = Vec::new();
        write_indices(&mut buf, &mesh.faces, limit).unwrap();

        let mut cursor = std::io::Cursor::new(&buf);
        let recovered: Vec<[u32; 3]> = read_indices(&mut cursor, limit).unwrap();

        assert_eq!(cursor.position() as usize, buf.len(), "{name}");
        assert_eq!(recovered, mesh.faces, "{name}");
    }
}

/// Real data should not need the widths a whole-byte scheme charged. Asserted rather than merely
/// reported, so a regression in width selection fails the build instead of quietly showing up as a
/// worse number in the size report.
#[test]
fn real_meshes_beat_the_whole_byte_scheme() {
    for (name, tol) in [("bunny.ply", 1e-5), ("scan-chunk.ply", 1e-3)] {
        let mesh = ply::load(name);
        let extents = mesh.extents();
        let axis_tol = tol / 3f64.sqrt();

        let old_coords: usize = mesh.points.len()
            * extents
                .iter()
                .map(|&e| old_bytes_for_tol(e, axis_tol) as usize)
                .sum::<usize>();
        let old_indices =
            mesh.faces.len() * 3 * old_bytes_for_count(mesh.points.len() as u32) as usize;

        let coords = check_points(name, &mesh.points, tol);
        let mut index_buf = Vec::new();
        write_indices(&mut index_buf, &mesh.faces, mesh.points.len() as u32).unwrap();

        let new_total = coords + index_buf.len();
        let old_total = old_coords + old_indices;

        assert!(
            new_total < old_total,
            "{name}: {new_total} bytes is no better than the old {old_total}"
        );
    }
}

/// The container path, exercised through an actual file rather than an in-memory buffer, because
/// `BufReader`/`BufWriter` behaviour and path handling are part of what the container adds.
#[test]
fn real_meshes_round_trip_through_a_tcmesh_file() {
    for (name, tol) in [("bunny.ply", 1e-5), ("scan-chunk.ply", 1e-3)] {
        let mesh = ply::load(name);
        let item = tol_compress::Mesh3::new(mesh.points.clone(), mesh.faces.clone())
            .named(name.trim_end_matches(".ply"));

        let path = std::env::temp_dir().join(format!(
            "tol-compress-fixture-{}-{}.tcmesh",
            std::process::id(),
            name.trim_end_matches(".ply")
        ));

        tol_compress::mesh::write_one_file(&path, &item, tol).unwrap();
        let back = tol_compress::mesh::read_one_file(&path).unwrap();
        let size = std::fs::metadata(&path).unwrap().len();
        let _ = std::fs::remove_file(&path);

        assert_eq!(back.name.as_deref(), Some(name.trim_end_matches(".ply")));
        assert_eq!(back.points.len(), mesh.points.len(), "{name}");
        assert_eq!(back.faces.len(), mesh.faces.len(), "{name}");

        // `write_one_file` renumbers vertices and faces, so the arrays are not the ones that went
        // in. What must hold is that it is the same mesh: every triangle joining the same three
        // positions it joined before, which is what a points array permuted out of step with its
        // faces would break.
        let plan = tol_compress::reorder::optimize(&mesh.faces, mesh.points.len()).unwrap();

        for (new, &old) in plan.face_order.iter().enumerate() {
            for corner in 0..3 {
                let got = back.points[back.faces[new][corner] as usize];
                let want = mesh.points[mesh.faces[old as usize][corner] as usize];

                let d = distance(&got, &want);
                assert!(
                    d <= tol,
                    "{name}: face {new} corner {corner} landed {d} away"
                );
            }
        }

        // And the order-preserving path must hand the arrays back untouched.
        tol_compress::mesh::write_one_file_preserving_order(&path, &item, tol).unwrap();
        let exact = tol_compress::mesh::read_one_file(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(exact.faces, mesh.faces, "{name}: preserved connectivity");
        for (i, (o, r)) in mesh.points.iter().zip(exact.points.iter()).enumerate() {
            let d = distance(o, r);
            assert!(d <= tol, "{name}: vertex {i} recovered {d} away");
        }

        // Container framing is 11 bytes plus a per-item preamble, so it must be negligible against
        // the geometry rather than a meaningful share of the file.
        let raw = mesh.points.len() * 3 * 4 + mesh.faces.len() * 3 * 4;
        assert!(
            size < raw as u64,
            "{name}: {size} bytes is no better than {raw} raw f32 bytes"
        );
    }
}

/// The previous whole-byte width rules, reproduced so the comparison is exact and needs no
/// dependency on `engeom`.
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
