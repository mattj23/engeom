//! Conformal mapping

use super::MeshEdges;
use crate::common::points::dist;
use crate::{Point2, Point3, Result};
use faer::Mat;
use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::Lu;
use faer::sparse::{SparseColMat, Triplet};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

type SparseMat = SparseColMat<u32, f64>;

impl MeshEdges<'_> {
    pub fn boundary_first_flatten(&self) -> Result<Vec<Point2>> {
        let n_vert = self.vertices().len();

        // Get a single boundary loop or return an error
        if self.boundary_loops.len() != 1 {
            return Err("Mesh must have a single boundary loop".into());
        }
        let i_bound = self.boundary_loops[0].as_slice();

        // Get the inner vertices
        let i_inner = inner_vertices(self, i_bound)?;

        // Calculate the face angles, which will be used both to calculate the angle defects for the
        // dirichlet boundary condition and to calculate the cotangent weights for the laplacian.
        let face_angles = calc_face_angles(self)?;

        // Calculate the cluster of values for the cotangent laplacian matrix and its sub-matrices.
        // We'll create the A matrix (the full cotangent laplacian) and the A_ii, A_ib, and A_bb
        // sub-matrices, which are the rows of the inner vertices against the inner vertices, the rows
        // of the inner vertices against the boundary vertices, and the rows of the boundary vertices
        // against the boundary vertices, respectively.
        let triplets =
            cotan_laplacian_triplets(&face_angles, n_vert, &self.edges, &self.face_edges)?;
        let (a, aii, aib, abb) = laplacian_set(n_vert, &i_inner, i_bound, &triplets)?;

        // We'll pre-factor the A matrix and the A_ii matrix, the latter which will be used while
        // setting the dirichlet boundary condition and determining the positions of the boundary
        // vertices in the final layout, and the former which will be used to extend the outer layout
        // boundary vertices into the interior of the layout.
        let a_lu = a.sp_lu()?;
        let aii_lu = aii.sp_lu()?;

        // Calculate the angle defects to be used in the dirichlet boundary condition.
        let angle_defects = calc_angle_defects(n_vert, i_bound, &face_angles, self.faces())?;

        // The ub vector is set to zeroes when using the minimum distortion boundary condition.
        let ub = Mat::<f64>::zeros(i_bound.len(), 1);

        // Set the target im_k vector for the boundary vertices.
        let im_k = dirichlet_boundary(&ub, &aii_lu, &aib, &abb, &i_inner, i_bound, &angle_defects)?;

        // Calculate the uv positions of the boundary vertices in the final layout
        let boundary_edge_len = boundary_edge_lengths(self.vertices(), i_bound);
        let uvb = best_fit_curve(&ub, &im_k, &boundary_edge_len)?;

        // Finally, extend the boundary vertices into the interior of the layout.
        let uv = extend_curve(&a_lu, &aii_lu, &aib, &uvb, n_vert, i_bound, &i_inner)?;

        Ok(uv.iter().map(|row| Point2::new(row[0], row[1])).collect())
    }
}

fn inner_vertices(mesh: &MeshEdges, boundary_vertices: &[u32]) -> Result<Vec<u32>> {
    let boundary_set: HashSet<u32> = boundary_vertices.iter().cloned().collect();
    let inner: Vec<u32> = (0..mesh.vertices().len() as u32)
        .filter(|&i| !boundary_set.contains(&i))
        .collect();
    Ok(inner)
}

// Quick helper for making a single column matrix from slice
fn single_col_matrix(data: &[f64]) -> Mat<f64> {
    Mat::from_fn(data.len(), 1, |i, _| data[i])
}

fn invert_2x2(m: &Mat<f64>) -> Result<Mat<f64>> {
    if m.nrows() != 2 || m.ncols() != 2 {
        return Err("Matrix is not 2x2".into());
    }

    let det = m[(0, 0)] * m[(1, 1)] - m[(0, 1)] * m[(1, 0)];
    if det == 0.0 {
        return Err("Matrix is singular".into());
    }

    let inv_det = 1.0 / det;
    let mut result = Mat::zeros(2, 2);
    result[(0, 0)] = m[(1, 1)] * inv_det;
    result[(0, 1)] = -m[(0, 1)] * inv_det;
    result[(1, 0)] = -m[(1, 0)] * inv_det;
    result[(1, 1)] = m[(0, 0)] * inv_det;

    Ok(result)
}

fn cumulative_sum(a: &[f64], scale: f64) -> Vec<f64> {
    let mut sum = 0.0;
    a.iter()
        .map(|&x| {
            sum += x * scale;
            sum
        })
        .collect()
}

fn zip_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect()
}

///
///
/// # Arguments
///
/// * `ub`: an array with a single row and `b_lengths.len()` columns
/// * `b_lengths`:
///
/// returns: Mat<f64, usize, usize>
fn calc_im_elen(ub: &Mat<f64>, b_lengths: &[f64]) -> Mat<f64> {
    let ub_slice = ub.col_as_slice(0);
    let values = b_lengths
        .iter()
        .enumerate()
        .map(|(i, &l)| {
            let a = (ub_slice[i] + ub_slice[(i + 1) % ub_slice.len()]) / 2.0;
            l * a.exp()
        })
        .collect::<Vec<_>>();

    single_col_matrix(&values)
}

fn calc_best_fit_tangents(im_k: &[f64]) -> Mat<f64> {
    let phi = cumulative_sum(im_k, -1.0);
    Mat::from_fn(phi.len(), 2, |i, j| match j {
        0 => phi[i].cos(),
        1 => phi[i].sin(),
        _ => unreachable!(),
    })
}

fn calc_boundary_vertex_masses(b_lengths: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; b_lengths.len()];
    for (i, l) in b_lengths.iter().enumerate() {
        let ni = (i + 1) % b_lengths.len();
        result[ni] = (l + b_lengths[ni]) / 2.0;
    }
    result
}

fn calc_best_fit_n1(b_lengths: &[f64]) -> Result<SparseMat> {
    let bvm = calc_boundary_vertex_masses(b_lengths);
    SparseMat::try_new_from_triplets(
        bvm.len(),
        bvm.len(),
        &bvm.iter()
            .enumerate()
            .map(|(i, &v)| Triplet::new(i as u32, i as u32, v))
            .collect::<Vec<_>>(),
    )
    .map_err(Into::into)
}

fn modified_im_elen(im_elen: &Mat<f64>, n1: &SparseMat, tangents: &Mat<f64>) -> Mat<f64> {
    let core = invert_2x2(&(&tangents.transpose() * n1 * tangents)).unwrap();
    let im_elen_sub = n1 * tangents * core * tangents.transpose() * im_elen;
    im_elen - im_elen_sub
}

fn best_fit_curve(ub: &Mat<f64>, im_k: &[f64], b_lengths: &[f64]) -> Result<Mat<f64>> {
    let tangents = calc_best_fit_tangents(im_k);
    let im_elen = calc_im_elen(ub, b_lengths);
    let n1 = calc_best_fit_n1(b_lengths)?;

    // Modify the im_elen matrix
    let im_elen = modified_im_elen(&im_elen, &n1, &tangents);

    // Any negative values in im_elen indicate that there is a boundary edge with a
    // negative length.
    if im_elen.col_as_slice(0).iter().any(|&x| x < 0.0) {
        return Err("Negative values in im_elen".into());
    }

    // We'll do a row-wise multiplication of im_elen against the tangents matrix, so that the
    // first column of the result is im_elen[(i, 0)] * tangents[(i, 0)] and the second is
    // im_elen[(i, 0)] * tangents[(i, 1)]
    let col0 = zip_product(im_elen.col_as_slice(0), tangents.col_as_slice(0));
    let col1 = zip_product(im_elen.col_as_slice(0), tangents.col_as_slice(1));

    // We'll take the cumulative sum of each
    let col0 = cumulative_sum(&col0, 1.0);
    let col1 = cumulative_sum(&col1, 1.0);

    // Finally, the result will be combined into a 2-column matrix and rolled forward by one
    Ok(Mat::from_fn(col0.len(), 2, |i, j| {
        let insert_i = ((i + col0.len()) - 1) % col0.len();
        match j {
            0 => col0[insert_i],
            1 => col1[insert_i],
            _ => unreachable!(),
        }
    }))
}

fn calc_extend_h(n_vert: usize, i_bound: &[u32], uvb: &Mat<f64>) -> Mat<f64> {
    let mut h = Mat::zeros(n_vert, 1);
    let bn = uvb.shape().0;
    for (i_b, &i_all) in i_bound.iter().enumerate() {
        // The value `i_b` is the index of the boundary vertex in the `i_bound` array, and `i_all`
        // is the index of the same boundary vertex in the `vertices` array.
        let i_b_prev = (i_b + bn - 1) % bn;
        let i_b_next = (i_b + 1) % bn;
        h[(i_all as usize, 0)] = 0.5 * (uvb[(i_b_prev, 0)] - uvb[(i_b_next, 0)]);
    }

    h
}

fn calc_extend_uv_xs(
    n_vert: usize,
    i_bound: &[u32],
    i_inner: &[u32],
    uvb: &Mat<f64>,
    aib: &SparseMat,
    aii_lu: &Lu<u32, f64>,
) -> Mat<f64> {
    let mut uv = Mat::zeros(n_vert, 2);

    // Copy the boundary vertex x values into the first column of the uv matrix from the ubv matrix
    for (&i_b, &v) in i_bound.iter().zip(uvb.col_as_slice(0)) {
        uv[(i_b as usize, 0)] = v;
    }

    let core = aib * &uvb.subcols(0, 1);
    let uv_inner = aii_lu.solve(-core);

    for (&i, &v) in i_inner.iter().zip(uv_inner.col_as_slice(0)) {
        uv[(i as usize, 0)] = v;
    }

    uv
}

fn extend_curve(
    a_lu: &Lu<u32, f64>,
    aii_lu: &Lu<u32, f64>,
    aib: &SparseMat,
    uvb: &Mat<f64>,
    n_vert: usize,
    i_bound: &[u32],
    i_inner: &[u32],
) -> Result<Vec<[f64; 2]>> {
    // Calculate the x values for all vertices
    let uv = calc_extend_uv_xs(n_vert, i_bound, i_inner, uvb, aib, aii_lu);

    // Solve for the y values
    let h = calc_extend_h(n_vert, i_bound, uvb);
    let y = a_lu.solve(-&h);

    Ok(uv
        .col_as_slice(0)
        .iter()
        .zip(y.col_as_slice(0).iter())
        .map(|(&x, &y)| [x, y])
        .collect())
}

fn boundary_edge_lengths(vertices: &[Point3], i_bound: &[u32]) -> Vec<f64> {
    i_bound
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let next = i_bound[(i + 1) % i_bound.len()];
            dist(&vertices[v as usize], &vertices[next as usize])
        })
        .collect()
}

/// Calculate the angles of the faces in the mesh. Each element `i` of the returned list corresponds
/// with the face at `mesh.faces[i]`, while the value at `j` for that face corresponds with the
/// angle formed at the vertex opposite to the `j`th edge `mesh.face_edges[i][j]`.
fn calc_face_angles(mesh: &MeshEdges) -> Result<Vec<[f64; 3]>> {
    use std::f64::consts::PI;
    let mut angles: Vec<[f64; 3]> = Vec::with_capacity(mesh.faces().len());

    for face_indices in mesh.face_edges.iter() {
        let a = mesh.edge_lengths[face_indices[0] as usize];
        let b = mesh.edge_lengths[face_indices[1] as usize];
        let c = mesh.edge_lengths[face_indices[2] as usize];

        // Check for degenerate faces
        let face_angles = if a > b + c {
            [PI, 0.0, 0.0]
        } else if b > a + c {
            [0.0, PI, 0.0]
        } else if c > a + b {
            [0.0, 0.0, PI]
        } else {
            let cos_a = (b.powi(2) + c.powi(2) - a.powi(2)) / (2.0 * b * c);
            let cos_b = (a.powi(2) + c.powi(2) - b.powi(2)) / (2.0 * a * c);
            let cos_c = (a.powi(2) + b.powi(2) - c.powi(2)) / (2.0 * a * b);
            [cos_a.acos(), cos_b.acos(), cos_c.acos()]
        };

        angles.push(face_angles);
    }

    Ok(angles)
}

fn calc_angle_defects(
    n: usize,
    i_bound: &[u32],
    face_angles: &[[f64; 3]],
    faces: &[[u32; 3]],
) -> Result<Vec<f64>> {
    let mut thetas = vec![2.0 * PI; n];
    for &i in i_bound {
        thetas[i as usize] = PI;
    }

    for (face, angles) in faces.iter().zip(face_angles.iter()) {
        thetas[face[0] as usize] -= angles[0];
        thetas[face[1] as usize] -= angles[1];
        thetas[face[2] as usize] -= angles[2];
    }

    Ok(thetas)
}

fn cotan_laplacian_triplets(
    face_angles: &[[f64; 3]],
    n_vert: usize,
    edges: &[[u32; 2]],
    face_edges: &[[u32; 3]],
) -> Result<Vec<Triplet<u32, u32, f64>>> {
    let cotans = face_angles
        .iter()
        .map(|angles| {
            [
                1.0 / angles[0].tan(),
                1.0 / angles[1].tan(),
                1.0 / angles[2].tan(),
            ]
        })
        .collect::<Vec<[f64; 3]>>();

    let mut values = vec![0.0; edges.len()];
    for (face, cotan) in face_edges.iter().zip(cotans.iter()) {
        for (i, &edge) in face.iter().enumerate() {
            values[edge as usize] += cotan[i];
        }
    }

    // Multiply by 0.5 to account for the fact that each edge is shared by two faces
    for value in values.iter_mut() {
        *value *= 0.5;
    }

    // Prepare the diagonal values
    let mut diagonals = vec![0.0; n_vert];
    for (edge, &value) in edges.iter().zip(values.iter()) {
        diagonals[edge[0] as usize] += value;
        diagonals[edge[1] as usize] += value;
    }

    // Build the sparse matrix
    let mut triplets = Vec::new();
    for (i, &value) in diagonals.iter().enumerate() {
        // The 1e-8 is added for stability (ensures the matrix is positive definite?)
        triplets.push(Triplet::new(i as u32, i as u32, value + 1e-8));
    }

    for (edge, &value) in edges.iter().zip(values.iter()) {
        triplets.push(Triplet::new(edge[0], edge[1], -value));
        triplets.push(Triplet::new(edge[1], edge[0], -value));
    }

    Ok(triplets)
}

fn slice_triplets_to_sparse(
    rows: &[u32],
    cols: &[u32],
    triplets: &[Triplet<u32, u32, f64>],
) -> Result<SparseMat> {
    let row_check: HashMap<u32, u32> = rows
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i as u32))
        .collect();
    let col_check: HashMap<u32, u32> = cols
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i as u32))
        .collect();

    let updated = triplets
        .iter()
        .filter_map(|t| {
            if let (Some(&row), Some(&col)) = (row_check.get(&t.row), col_check.get(&t.col)) {
                Some(Triplet::new(row, col, t.val))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    SparseColMat::try_new_from_triplets(rows.len(), cols.len(), &updated).map_err(|e| e.into())
}

/// Calculate the set of sparse laplacian matrices A, AII, AIB, and ABB for the given data.
///
/// # Arguments
///
/// * `n`: the number of vertices in the mesh
/// * `i_inner`: the indices of the inner vertices, in their original order
/// * `i_bound`: the indices of the boundary vertices, in their original order
/// * `triplets`: the total set of triplets for the cotangent laplacian matrix
///
/// returns: Result<(SparseColMat<u32, f64, usize, usize>, SparseColMat<u32, f64, usize, usize>, SparseColMat<u32, f64, usize, usize>, SparseColMat<u32, f64, usize, usize>), Box<dyn Error, Global>>
fn laplacian_set(
    n: usize,
    i_inner: &[u32],
    i_bound: &[u32],
    triplets: &[Triplet<u32, u32, f64>],
) -> Result<(SparseMat, SparseMat, SparseMat, SparseMat)> {
    let a = SparseColMat::try_new_from_triplets(n, n, triplets)?;
    let aii = slice_triplets_to_sparse(i_inner, i_inner, triplets)?;
    let aib = slice_triplets_to_sparse(i_inner, i_bound, triplets)?;
    let abb = slice_triplets_to_sparse(i_bound, i_bound, triplets)?;

    Ok((a, aii, aib, abb))
}

fn dirichlet_boundary(
    ub: &Mat<f64>,
    aii_lu: &Lu<u32, f64>,
    aib: &SparseMat,
    abb: &SparseMat,
    i_inner: &[u32],
    i_bounds: &[u32],
    angle_defects: &[f64],
) -> Result<Vec<f64>> {
    let inner_defects = i_inner
        .iter()
        .map(|i| angle_defects[*i as usize])
        .collect::<Vec<_>>();

    let defects = Mat::from_fn(inner_defects.len(), 1, |i, _| inner_defects[i]);

    let value = &defects + aib * ub;
    let ui = -aii_lu.solve(&value);

    let h = -aib.transpose() * &ui - abb * ub;
    let h = h.row_iter().map(|r| r[0]).collect::<Vec<f64>>();

    let im_k = i_bounds
        .iter()
        .zip(h.iter())
        .map(|(&vi, &hv)| angle_defects[vi as usize] - hv)
        .collect::<Vec<f64>>();

    Ok(im_k)
}
