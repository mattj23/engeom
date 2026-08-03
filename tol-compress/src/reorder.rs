//! Reordering a mesh's faces and vertices so its indices compress.
//!
//! An index block costs `log2(V)` bits per corner no matter what, unless the indices themselves can
//! be made small. Two reorderings together make them small, and both are derived entirely from
//! connectivity, so **no permutation is ever stored**: an encoder applies them, and a decoder
//! recovers the same arrangement for free by reading the streams in the order they were written.
//! Decoding costs nothing extra at all.
//!
//! This is the cheap majority of what Edgebreaker buys without its manifold requirement or its
//! traversal machinery, which matters because metrology scan meshes are frequently non-manifold or
//! boundary-heavy.
//!
//! # The two orderings
//!
//! **Faces** are visited breadth-first over adjacency, so a face is followed by faces near it. On
//! `stanford_bun_2`, a zippered range-image merge whose face order jumps across the surface,
//! reordering vertices gives -5.5% and adding the face traversal gives -39.6%.
//!
//! Adjacency here is by shared vertex, not shared edge. It needs no half-edge structure and
//! makes no manifold assumption, and it is what the measurements were taken with.
//!
//! **Vertices** are then numbered by first use across that face order, so every index is at most
//! the running count of vertices introduced so far. That is what lets [`crate::indices`] code a
//! corner as a small distance back from a high-water mark instead of an absolute index. Vertices no
//! face references keep their relative order and are appended after the referenced ones.
//!
//! # Applying it
//!
//! [`Reordering`] hands back both permutations as well as the remapped faces, because a caller has
//! to move _everything_ that is indexed by vertex or by face, not just the positions.

use crate::error::{Error, Result};
use std::collections::VecDeque;

/// A mesh's faces and vertices renumbered for compact index coding.
///
/// Both orders map **new position to old position**, so `old_points[vertex_order[i]]` is the new
/// vertex `i`. [`permute`] applies that.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reordering<const N: usize> {
    /// New vertex index to old vertex index. Always covers every vertex, referenced or not.
    pub vertex_order: Vec<u32>,
    /// New face index to old face index.
    pub face_order: Vec<u32>,
    /// The faces in the new order, with corners renumbered into the new vertex numbering.
    pub faces: Vec<[u32; N]>,
}

/// Rearrange `items` into a new order, where `order` maps new position to old.
///
/// This is what keeps per-vertex and per-face data lined up with the geometry it belongs to.
///
/// # Panics
///
/// Panics if any entry of `order` is out of range for `items`.
pub fn permute<T: Clone>(items: &[T], order: &[u32]) -> Vec<T> {
    order.iter().map(|&i| items[i as usize].clone()).collect()
}

/// Order faces breadth-first over shared-vertex adjacency, returning new face index to old.
///
/// Faces belonging to separate connected components follow one another in the order their first
/// face appeared, so a disconnected mesh is handled without any special case.
///
/// # Errors
///
/// [`Error::Malformed`] if a face refers to a vertex at or beyond `vertex_count`, or if there are
/// more faces than a `u32` can count.
pub fn greedy_face_order<const N: usize>(
    faces: &[[u32; N]],
    vertex_count: usize,
) -> Result<Vec<u32>> {
    if u32::try_from(faces.len()).is_err() {
        return Err(Error::Malformed(
            "mesh holds more faces than a u32 can count",
        ));
    }

    let (starts, incident) = vertex_to_faces(faces, vertex_count)?;

    let mut visited = vec![false; faces.len()];
    let mut order = Vec::with_capacity(faces.len());
    let mut queue: VecDeque<u32> = VecDeque::new();

    for seed in 0..faces.len() {
        if visited[seed] {
            continue;
        }

        // Marking on push rather than on pop is what keeps a face out of the queue twice.
        visited[seed] = true;
        queue.push_back(seed as u32);

        while let Some(face) = queue.pop_front() {
            order.push(face);

            for &v in &faces[face as usize] {
                let v = v as usize;
                for &neighbour in &incident[starts[v]..starts[v + 1]] {
                    if !visited[neighbour as usize] {
                        visited[neighbour as usize] = true;
                        queue.push_back(neighbour);
                    }
                }
            }
        }
    }

    Ok(order)
}

/// Reorder a mesh's faces and then number its vertices by first use.
///
/// # Errors
///
/// [`Error::Malformed`] if a face refers to a vertex at or beyond `vertex_count`, or if either
/// count is too large for a `u32`.
pub fn optimize<const N: usize>(faces: &[[u32; N]], vertex_count: usize) -> Result<Reordering<N>> {
    if u32::try_from(vertex_count).is_err() {
        return Err(Error::Malformed(
            "mesh holds more vertices than a u32 can index",
        ));
    }

    let face_order = greedy_face_order(faces, vertex_count)?;

    const UNSEEN: u32 = u32::MAX;
    let mut new_id = vec![UNSEEN; vertex_count];
    let mut vertex_order: Vec<u32> = Vec::with_capacity(vertex_count);
    let mut out_faces: Vec<[u32; N]> = Vec::with_capacity(faces.len());

    for &old_face in &face_order {
        let mut face = [0u32; N];

        for (corner, &v) in faces[old_face as usize].iter().enumerate() {
            let v = v as usize;
            if new_id[v] == UNSEEN {
                new_id[v] = vertex_order.len() as u32;
                vertex_order.push(v as u32);
            }
            face[corner] = new_id[v];
        }

        out_faces.push(face);
    }

    // Vertices no face touches still have to survive the round trip, so they follow the referenced
    // ones with their relative order intact.
    for (v, id) in new_id.iter_mut().enumerate() {
        if *id == UNSEEN {
            *id = vertex_order.len() as u32;
            vertex_order.push(v as u32);
        }
    }

    Ok(Reordering {
        vertex_order,
        face_order,
        faces: out_faces,
    })
}

/// Which faces touch each vertex, as a flat array plus per-vertex start offsets.
///
/// A compressed adjacency layout rather than a vector of vectors: one allocation instead of one per
/// vertex, and the neighbours of a vertex are contiguous.
fn vertex_to_faces<const N: usize>(
    faces: &[[u32; N]],
    vertex_count: usize,
) -> Result<(Vec<usize>, Vec<u32>)> {
    let mut starts = vec![0usize; vertex_count + 1];

    for face in faces {
        for &v in face {
            let v = v as usize;
            if v >= vertex_count {
                return Err(Error::Malformed(
                    "simplex refers to an index at or beyond the item count",
                ));
            }
            starts[v + 1] += 1;
        }
    }

    for i in 0..vertex_count {
        starts[i + 1] += starts[i];
    }

    let mut incident = vec![0u32; starts[vertex_count]];
    let mut cursor = starts.clone();

    for (index, face) in faces.iter().enumerate() {
        for &v in face {
            let v = v as usize;
            incident[cursor[v]] = index as u32;
            cursor[v] += 1;
        }
    }

    Ok((starts, incident))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corpus;
    use crate::testgen::Rng;

    /// The property the index coder depends on: every corner is at most the number of vertices
    /// introduced before it.
    fn assert_first_use_order<const N: usize>(faces: &[[u32; N]]) {
        let mut high_water = 0u32;
        for (f, face) in faces.iter().enumerate() {
            for (corner, &v) in face.iter().enumerate() {
                assert!(
                    v <= high_water,
                    "face {f} corner {corner} is {v}, past the high water mark {high_water}"
                );
                if v == high_water {
                    high_water += 1;
                }
            }
        }
    }

    /// A reordering has to describe the same mesh, which means the same triangles over the same
    /// vertices once both permutations are undone.
    fn assert_same_mesh<const N: usize>(
        original: &[[u32; N]],
        vertex_count: usize,
        r: &Reordering<N>,
    ) {
        assert_eq!(r.faces.len(), original.len(), "face count");
        assert_eq!(r.face_order.len(), original.len(), "face order length");
        assert_eq!(r.vertex_order.len(), vertex_count, "vertex order length");

        let mut seen = vec![false; vertex_count];
        for &v in &r.vertex_order {
            assert!(!seen[v as usize], "vertex {v} appears twice in the order");
            seen[v as usize] = true;
        }

        let mut seen_face = vec![false; original.len()];
        for &f in &r.face_order {
            assert!(
                !seen_face[f as usize],
                "face {f} appears twice in the order"
            );
            seen_face[f as usize] = true;
        }

        for (new, &old) in r.face_order.iter().enumerate() {
            for (corner, &want) in original[old as usize].iter().enumerate() {
                let recovered = r.vertex_order[r.faces[new][corner] as usize];
                assert_eq!(
                    recovered, want,
                    "face {new} corner {corner} does not map back"
                );
            }
        }
    }

    #[test]
    fn every_corpus_mesh_reorders_into_first_use_order() {
        for case in corpus::all() {
            if case.faces.is_empty() {
                continue;
            }

            let r = optimize(&case.faces, case.points.len()).unwrap();
            assert_first_use_order(&r.faces);
            assert_same_mesh(&case.faces, case.points.len(), &r);
        }
    }

    /// The whole scheme has to depend on connectivity and not on the order the caller happened to
    /// supply, which is exactly what `shuffled` exists to check.
    #[test]
    fn a_shuffled_mesh_reorders_just_as_well() {
        let ordered = corpus::smooth_surface();
        let shuffled = corpus::shuffled();

        let a = optimize(&ordered.faces, ordered.points.len()).unwrap();
        let b = optimize(&shuffled.faces, shuffled.points.len()).unwrap();

        assert_first_use_order(&a.faces);
        assert_first_use_order(&b.faces);

        // The measure that matters is how far back the corners reach, since that is what the index
        // coder pays for. Shuffling the input must not make it materially worse.
        let span_a = mean_high_water_delta(&a.faces);
        let span_b = mean_high_water_delta(&b.faces);
        assert!(
            span_b < span_a * 1.5,
            "shuffled input reached back {span_b} on average against {span_a} for ordered input"
        );
    }

    fn mean_high_water_delta<const N: usize>(faces: &[[u32; N]]) -> f64 {
        let mut high_water = 0u32;
        let mut total = 0u64;
        let mut count = 0u64;

        for face in faces {
            for &v in face {
                total += u64::from(high_water - v);
                count += 1;
                if v == high_water {
                    high_water += 1;
                }
            }
        }

        total as f64 / count.max(1) as f64
    }

    #[test]
    fn a_disconnected_mesh_keeps_every_component() {
        // Two separate triangles sharing nothing.
        let faces = [[0u32, 1, 2], [3, 4, 5]];
        let r = optimize(&faces, 6).unwrap();

        assert_first_use_order(&r.faces);
        assert_same_mesh(&faces, 6, &r);
    }

    #[test]
    fn unreferenced_vertices_survive() {
        // Vertices 1 and 4 are touched by nothing.
        let faces = [[0u32, 2, 3]];
        let r = optimize(&faces, 5).unwrap();

        assert_eq!(r.vertex_order.len(), 5);
        assert_same_mesh(&faces, 5, &r);

        // The referenced ones come first, and the rest keep their relative order behind them.
        assert_eq!(&r.vertex_order[..3], &[0, 2, 3]);
        assert_eq!(&r.vertex_order[3..], &[1, 4]);
    }

    #[test]
    fn a_mesh_with_no_faces_is_the_identity() {
        let faces: Vec<[u32; 3]> = Vec::new();
        let r = optimize(&faces, 4).unwrap();

        assert!(r.faces.is_empty());
        assert!(r.face_order.is_empty());
        assert_eq!(r.vertex_order, vec![0, 1, 2, 3]);
    }

    #[test]
    fn an_empty_mesh_is_handled() {
        let faces: Vec<[u32; 3]> = Vec::new();
        let r = optimize(&faces, 0).unwrap();

        assert!(r.vertex_order.is_empty());
        assert!(r.faces.is_empty());
    }

    /// Non-manifold and boundary-heavy geometry must go through unremarkably, since shared-vertex
    /// adjacency was chosen precisely so that scan data needs no special handling.
    #[test]
    fn boundary_heavy_geometry_reorders() {
        let case = corpus::boundary_heavy();
        let r = optimize(&case.faces, case.points.len()).unwrap();

        assert_first_use_order(&r.faces);
        assert_same_mesh(&case.faces, case.points.len(), &r);
    }

    /// A vertex used three times by one face is degenerate but not malformed, and must not confuse
    /// the adjacency build.
    #[test]
    fn a_degenerate_face_is_tolerated() {
        let faces = [[0u32, 0, 0], [0, 1, 2]];
        let r = optimize(&faces, 3).unwrap();

        assert_first_use_order(&r.faces);
        assert_same_mesh(&faces, 3, &r);
    }

    #[test]
    fn an_out_of_range_index_is_refused() {
        let faces = [[0u32, 1, 9]];

        assert!(matches!(optimize(&faces, 3), Err(Error::Malformed(_))));
        assert!(matches!(
            greedy_face_order(&faces, 3),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn permute_moves_data_with_its_vertices() {
        let values = ['a', 'b', 'c', 'd'];
        let order = vec![2u32, 0, 3, 1];

        assert_eq!(permute(&values, &order), vec!['c', 'a', 'd', 'b']);
    }

    /// The reordering must line data up with the geometry: after permuting, a vertex's payload has
    /// to still belong to the position at the same index.
    #[test]
    fn permuting_data_keeps_it_with_its_vertex() {
        let mut rng = Rng::new(4_100);
        let case = corpus::smooth_surface();
        let count = case.points.len();

        // A distinct tag per vertex, so a misplacement cannot go unnoticed.
        let tags: Vec<u64> = (0..count).map(|_| rng.next_u64()).collect();

        let r = optimize(&case.faces, count).unwrap();
        let moved_points = permute(&case.points, &r.vertex_order);
        let moved_tags = permute(&tags, &r.vertex_order);

        for i in 0..count {
            let old = r.vertex_order[i] as usize;
            assert_eq!(moved_points[i], case.points[old]);
            assert_eq!(moved_tags[i], tags[old]);
        }
    }
}
