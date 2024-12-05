use crate::geom3::Mesh;
use crate::Point3;
use std::collections::{HashMap, HashSet};

fn edge_key(i: usize, f: &[u32; 3]) -> (u32, u32) {
    (f[i], f[(i + 1) % 3])
}

fn make_sym(k: &(u32, u32)) -> (u32, u32) {
    if k.0 < k.1 {
        (k.0, k.1)
    } else {
        (k.1, k.0)
    }
}

pub fn compute_boundary_edges(mesh: &Mesh, patch: &[usize]) -> Vec<(u32, u32)> {
    // Boundary edges are ones which belong only to one face in the patch.  We will start by
    // identifying all the edges in the patch, and then we will count the number of times each
    // edge appears in the patch.  If an edge appears only once, it is a boundary edge.
    let mut edge_counts = HashMap::new();
    let mut keys = HashSet::new();

    for face_index in patch {
        let face = &mesh.triangles()[*face_index];
        for i in 0..3 {
            let rk = edge_key(i, face);
            let sk = make_sym(&rk);

            keys.insert(rk);
            *edge_counts.entry(sk).or_insert(0) += 1;
        }
    }

    keys.into_iter()
        .filter(|k| edge_counts[&make_sym(k)] == 1)
        .collect()
}

fn take_one_boundary(order: &mut HashMap<u32, u32>) -> Option<Vec<u32>> {
    let start = *order.keys().next()?;
    let mut sequence = vec![start];
    let mut next = order.remove(&start)?;

    while next != start {
        sequence.push(next);
        next = order.remove(&next)?;
    }

    Some(sequence)
}

pub fn compute_boundary_points(mesh: &Mesh, patch: &[usize]) -> Vec<Vec<Point3>> {
    let mut order = HashMap::new();
    let mut remaining = HashSet::new();
    let edges = compute_boundary_edges(&mesh, patch);
    for (v0, v1) in edges {
        if order.contains_key(&v0) {
            panic!("Duplicate boundary point found!");
        } else {
            order.insert(v0, v1);
        }
        remaining.insert(v0);
    }

    let mut sequences = Vec::new();
    while let Some(sequence) = take_one_boundary(&mut order) {
        sequences.push(sequence);
    }

    sequences
        .into_iter()
        .map(|seq| {
            seq.into_iter()
                .map(|v| mesh.vertices()[v as usize])
                .collect()
        })
        .collect()
}

/// Uses the adjacency information of the mesh to compute a list of patches of connected faces. This
/// function will not work with non-manifold meshes!.
///
/// # Arguments
///
/// * `mesh`:
///
/// returns: Vec<Vec<usize, Global>, Global>
///
/// # Examples
pub fn compute_patch_indices(mesh: &Mesh) -> Vec<Vec<usize>> {
    // We will start by computing a table of edges and the faces that contain them. Then we will
    // create a bin of remaining faces, pick one at random, and add it to a new patch, putting all
    // of its edges in a working queue.  We will then retrieve each remaining face associated with
    // each edge in the working queue, adding it to the patch and adding its edges to the working
    // queue.  When the working queue is empty, we will have completed one patch.  We will then
    // repeat the process until all faces are accounted for.

    // Compute the edge table
    let mut edge_table = HashMap::new();
    let mut remaining_faces = HashSet::new();

    for (i, face) in mesh.triangles().iter().enumerate() {
        edge_table.insert((face[0], face[1]), i);
        edge_table.insert((face[1], face[2]), i);
        edge_table.insert((face[2], face[0]), i);
        remaining_faces.insert(i);
    }

    // Create the patches
    let mut patches = Vec::new();

    while !remaining_faces.is_empty() {
        let mut working_queue = Vec::new();

        // Pick any face from the remaining faces and remove it from the set
        let face_index = *remaining_faces.iter().next().unwrap();
        remaining_faces.remove(&face_index);

        // Add the edges of the face to the working queue
        working_queue.push((
            mesh.triangles()[face_index][0],
            mesh.triangles()[face_index][1],
        ));
        working_queue.push((
            mesh.triangles()[face_index][1],
            mesh.triangles()[face_index][2],
        ));
        working_queue.push((
            mesh.triangles()[face_index][2],
            mesh.triangles()[face_index][0],
        ));

        let mut patch = vec![face_index];

        while !working_queue.is_empty() {
            let (v0, v1) = working_queue.pop().unwrap();
            let e0 = (v0, v1);
            let e1 = (v1, v0);

            if let Some(f0) = edge_table.get(&e0) {
                if remaining_faces.contains(f0) {
                    patch.push(*f0);
                    remaining_faces.remove(f0);
                    working_queue.push((mesh.triangles()[*f0][0], mesh.triangles()[*f0][1]));
                    working_queue.push((mesh.triangles()[*f0][1], mesh.triangles()[*f0][2]));
                    working_queue.push((mesh.triangles()[*f0][2], mesh.triangles()[*f0][0]));
                }
            }

            if let Some(f1) = edge_table.get(&e1) {
                if remaining_faces.contains(f1) {
                    patch.push(*f1);
                    remaining_faces.remove(f1);
                    working_queue.push((mesh.triangles()[*f1][0], mesh.triangles()[*f1][1]));
                    working_queue.push((mesh.triangles()[*f1][1], mesh.triangles()[*f1][2]));
                    working_queue.push((mesh.triangles()[*f1][2], mesh.triangles()[*f1][0]));
                }
            }
        }

        patches.push(patch);
    }

    patches
}
