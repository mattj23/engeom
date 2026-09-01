//! A benchmark corpus of deliberately awkward geometry.
//!
//! Everything the format does after basic bit packing is an optimization, and an optimization
//! measured only on data that suits it is not measured at all. These generators exist so that each
//! later technique can be shown the case designed to defeat it as well as the case designed to
//! flatter it.
//!
//! Nothing here is committed as a data file. Every case is generated from a fixed seed, so the
//! corpus is reproducible without putting meshes in the repository.
//!
//! Enable with the `corpus` feature. It is not part of the default build.

use crate::testgen::Rng;

/// One corpus entry: a point set, optionally with connectivity, and the tolerance it is intended
/// to be measured at.
#[derive(Debug, Clone)]
pub struct Case {
    /// Short identifier, used as a row label in reports.
    pub name: &'static str,
    /// What this case is designed to stress, and which technique it is evidence for or against.
    pub stresses: &'static str,
    /// The vertices.
    pub points: Vec<[f64; 3]>,
    /// Triangles, empty for point-only cases.
    pub faces: Vec<[u32; 3]>,
    /// The tolerance this case is meant to be encoded at.
    pub tol: f64,
}

impl Case {
    /// Whether this case carries connectivity as well as positions.
    pub fn is_mesh(&self) -> bool {
        !self.faces.is_empty()
    }

    /// The span of the point set along each axis.
    pub fn extents(&self) -> [f64; 3] {
        if self.points.is_empty() {
            return [0.0; 3];
        }
        let mut mins = self.points[0];
        let mut maxs = self.points[0];
        for p in &self.points {
            for i in 0..3 {
                mins[i] = mins[i].min(p[i]);
                maxs[i] = maxs[i].max(p[i]);
            }
        }
        std::array::from_fn(|i| maxs[i] - mins[i])
    }
}

/// Every case in the corpus.
pub fn all() -> Vec<Case> {
    let mut cases = vec![
        smooth_surface(),
        noisy_surface(0.05),
        shuffled(),
        irregular_tessellation(),
        boundary_heavy(),
        planar(),
        distant_clusters(),
        oblique_plate(),
        wide_dynamic_range(),
    ];
    cases.extend(degenerate());
    cases
}

// ---------------------------------------------------------------------------------------------
// Building blocks
// ---------------------------------------------------------------------------------------------

/// A triangulated height field over a regular grid.
///
/// Vertices run row-major, so neighbouring triangles reference nearby indices. That locality is
/// what index delta coding will later exploit, and what [`shuffled`] destroys on purpose.
fn grid(nx: usize, ny: usize, extent: f64, amplitude: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut points = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = extent * (i as f64 / (nx - 1) as f64 - 0.5);
            let y = extent * (j as f64 / (ny - 1) as f64 - 0.5);
            let z = amplitude
                * (3.0 * std::f64::consts::PI * x / extent).sin()
                * (2.0 * std::f64::consts::PI * y / extent).cos();
            points.push([x, y, z]);
        }
    }

    (points, grid_faces(nx, ny))
}

/// The triangle pairs for an `nx` by `ny` vertex grid.
fn grid_faces(nx: usize, ny: usize) -> Vec<[u32; 3]> {
    let mut faces = Vec::with_capacity((nx - 1) * (ny - 1) * 2);
    for j in 0..ny - 1 {
        for i in 0..nx - 1 {
            let a = (j * nx + i) as u32;
            let b = a + 1;
            let c = a + nx as u32;
            let d = c + 1;
            faces.push([a, b, d]);
            faces.push([a, d, c]);
        }
    }
    faces
}

/// Rodrigues rotation of a point about a unit axis.
fn rotate(p: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let k = [axis[0] / n, axis[1] / n, axis[2] / n];
    let (s, c) = angle.sin_cos();

    let dot = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
    let cross = [
        k[1] * p[2] - k[2] * p[1],
        k[2] * p[0] - k[0] * p[2],
        k[0] * p[1] - k[1] * p[0],
    ];

    std::array::from_fn(|i| p[i] * c + cross[i] * s + k[i] * dot * (1.0 - c))
}

// ---------------------------------------------------------------------------------------------
// Cases
// ---------------------------------------------------------------------------------------------

/// The favourable case: a smooth, regularly tessellated surface with well-ordered vertices.
///
/// Every later technique should do well here. It is the baseline the awkward cases are compared
/// against, not evidence of anything on its own.
pub fn smooth_surface() -> Case {
    let (points, faces) = grid(120, 120, 100.0, 8.0);
    Case {
        name: "smooth_surface",
        stresses: "nothing; the flattering baseline every other case is measured against",
        points,
        faces,
        tol: 1e-3,
    }
}

/// A smooth surface buried in measurement noise of the given standard deviation.
///
/// Parallelogram prediction assumes the surface is locally close to flat. Noise breaks that
/// assumption in proportion to its amplitude, so sweeping `sigma` finds the point where prediction
/// stops paying for itself. Real scan data sits somewhere on that sweep rather than at either end.
pub fn noisy_surface(sigma: f64) -> Case {
    let (mut points, faces) = grid(120, 120, 100.0, 8.0);
    let mut rng = Rng::new(9_001);
    for p in points.iter_mut() {
        for c in p.iter_mut() {
            *c += rng.gaussian(sigma);
        }
    }

    Case {
        name: "noisy_surface",
        stresses: "parallelogram prediction, which degrades as noise swamps local flatness",
        points,
        faces,
        tol: 1e-3,
    }
}

/// The smooth surface with its vertices randomly permuted and faces remapped to match.
///
/// Geometrically identical to [`smooth_surface`], so the coordinate block should be the same size.
/// Any difference between the two is entirely attributable to the loss of index locality, which
/// makes this the direct measurement of what a spatial reorder is worth.
pub fn shuffled() -> Case {
    let base = smooth_surface();
    let mut rng = Rng::new(9_002);

    let n = base.points.len();
    let mut order: Vec<usize> = (0..n).collect();
    rng.shuffle(&mut order);

    // `order[new] = old`, so the inverse is what remaps a face's old index to its new slot.
    let mut inverse = vec![0u32; n];
    for (new, &old) in order.iter().enumerate() {
        inverse[old] = new as u32;
    }

    let points = order.iter().map(|&old| base.points[old]).collect();
    let faces = base
        .faces
        .iter()
        .map(|f| {
            [
                inverse[f[0] as usize],
                inverse[f[1] as usize],
                inverse[f[2] as usize],
            ]
        })
        .collect();

    Case {
        name: "shuffled",
        stresses: "index delta coding and spatial reordering, by destroying vertex locality",
        points,
        faces,
        tol: 1e-3,
    }
}

/// A surface whose triangles vary in size by two orders of magnitude across the sheet.
///
/// Adaptive tessellation is normal in real CAD output. Prediction schemes that assume neighbouring
/// triangles are similar in size have their residuals blow up where the size changes fastest.
pub fn irregular_tessellation() -> Case {
    let (nx, ny) = (120usize, 120usize);
    let extent = 100.0;

    // Geometric spacing in x compresses one edge and stretches the other.
    let mut points = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let t = i as f64 / (nx - 1) as f64;
            let x = extent * ((10.0f64.powf(2.0 * t) - 1.0) / 99.0 - 0.5);
            let y = extent * (j as f64 / (ny - 1) as f64 - 0.5);
            let z = 4.0 * (2.0 * std::f64::consts::PI * x / extent).sin();
            points.push([x, y, z]);
        }
    }

    Case {
        name: "irregular_tessellation",
        stresses: "prediction's assumption that neighbouring triangles are similar in size",
        points,
        faces: grid_faces(nx, ny),
        tol: 1e-3,
    }
}

/// A scan-like sheet with holes punched through it and non-manifold fins attached.
///
/// Metrology meshes are routinely like this, which is why full Edgebreaker connectivity coding is
/// out of scope. Nothing in the format may assume two-manifold topology, and this is the case that
/// catches such an assumption if one creeps in.
pub fn boundary_heavy() -> Case {
    let (points, mut faces) = grid(100, 100, 100.0, 6.0);
    let mut rng = Rng::new(9_003);

    // Punch roughly a third of the faces out, producing many interior boundaries.
    faces.retain(|_| rng.next_f64() > 0.33);

    // Attach fins: triangles sharing an edge that already has two incident faces, which makes
    // those edges non-manifold.
    let n = points.len() as u32;
    let fin_count = faces.len() / 50;
    for _ in 0..fin_count {
        let f = faces[rng.below(faces.len())];
        let far = rng.below(n as usize) as u32;
        if far != f[0] && far != f[1] {
            faces.push([f[0], f[1], far]);
        }
    }

    Case {
        name: "boundary_heavy",
        stresses: "any assumption of manifold topology; holes and non-manifold edges throughout",
        points,
        faces,
        tol: 1e-3,
    }
}

/// A surface lying exactly in a plane, so one axis has zero extent.
///
/// The previous implementation charged two bytes per point for that axis. It should now cost
/// nothing at all, and every point should come back on the plane exactly rather than near it.
pub fn planar() -> Case {
    let (mut points, faces) = grid(120, 120, 100.0, 0.0);
    for p in points.iter_mut() {
        p[2] = 12.5;
    }

    Case {
        name: "planar",
        stresses: "the zero-bit axis path; a degenerate dimension must be free",
        points,
        faces,
        tol: 1e-3,
    }
}

/// Small dense clusters scattered across a large empty volume.
///
/// The global bounding box is almost entirely empty, so every point pays for a range it does not
/// occupy. This is the case that spatial partitioning exists to fix, and the one that should show
/// the least benefit from everything before it.
pub fn distant_clusters() -> Case {
    let mut rng = Rng::new(9_004);
    let mut points = Vec::with_capacity(8 * 500);

    for _ in 0..8 {
        let centre = rng.point::<3>(-500.0, 500.0);
        for _ in 0..500 {
            points.push(std::array::from_fn(|i| centre[i] + rng.gaussian(0.5)));
        }
    }

    Case {
        name: "distant_clusters",
        stresses: "spatial partitioning; the global box is mostly empty space",
        points,
        faces: Vec::new(),
        tol: 1e-3,
    }
}

/// A thin plate rotated so that it is oblique to all three axes.
///
/// Axis-aligned bounds cannot describe it tightly: the plate is 0.02 thick but every axis-aligned
/// extent is large. Reorientation should collapse one axis to nearly nothing, and until it lands
/// this case pays for a dimension that carries almost no information.
pub fn oblique_plate() -> Case {
    let (points, faces) = grid(120, 120, 100.0, 0.0);
    let axis = [1.0, 1.0, 1.0];
    let angle = 0.6;

    let points = points
        .into_iter()
        .map(|mut p| {
            p[2] = 0.01 * (p[0] * 0.3).sin();
            rotate(p, axis, angle)
        })
        .collect();

    Case {
        name: "oblique_plate",
        stresses: "reorientation; thin in one direction but wide along every axis",
        points,
        faces,
        tol: 1e-3,
    }
}

/// Two metres of part measured to a micron.
///
/// A range-to-tolerance ratio around two million, which is about as demanding as real metrology
/// gets. Exercises the wide end of the width table without reaching the representable limit.
pub fn wide_dynamic_range() -> Case {
    let (points, faces) = grid(100, 100, 2000.0, 40.0);

    Case {
        name: "wide_dynamic_range",
        stresses: "the wide end of the width table; a two million to one range/tolerance ratio",
        points,
        faces,
        tol: 1e-3,
    }
}

/// The trivially small inputs, where off-by-one handling tends to live.
pub fn degenerate() -> Vec<Case> {
    vec![
        Case {
            name: "empty",
            stresses: "empty input; no bounding box exists",
            points: Vec::new(),
            faces: Vec::new(),
            tol: 1e-3,
        },
        Case {
            name: "single_point",
            stresses: "a point set with zero extent on every axis",
            points: vec![[1.5, -2.5, 3.5]],
            faces: Vec::new(),
            tol: 1e-3,
        },
        Case {
            name: "two_points",
            stresses: "the narrowest non-degenerate range",
            points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            faces: Vec::new(),
            tol: 1e-3,
        },
        Case {
            name: "single_triangle",
            stresses: "the smallest possible mesh; index widths at their floor",
            points: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            faces: vec![[0, 1, 2]],
            tol: 1e-3,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Error, read_indices, read_points, write_indices, write_points};
    use std::io::Cursor;

    fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    /// The guarantee, checked on every case at every tolerance, for every point. This is the
    /// product; the compression ratio is a side effect.
    #[test]
    fn every_case_round_trips_within_tolerance() {
        for case in all() {
            for scale in [0.1, 1.0, 10.0] {
                let tol = case.tol * scale;

                let mut buf = Vec::new();
                write_points(&mut buf, &case.points, tol)
                    .unwrap_or_else(|e| panic!("{} at tol {tol}: {e}", case.name));

                let mut cursor = Cursor::new(&buf);
                let recovered: Vec<[f64; 3]> = read_points(&mut cursor).unwrap();

                assert_eq!(
                    cursor.position() as usize,
                    buf.len(),
                    "{}: decoder left bytes unread",
                    case.name
                );
                assert_eq!(recovered.len(), case.points.len(), "{}", case.name);

                for (i, (o, r)) in case.points.iter().zip(recovered.iter()).enumerate() {
                    let d = distance(o, r);
                    assert!(
                        d <= tol,
                        "{} at tol {tol}: point {i} recovered {d} away",
                        case.name
                    );
                }
            }
        }
    }

    /// Connectivity is not lossy at all, so anything short of exact equality is a bug.
    #[test]
    fn every_mesh_case_round_trips_connectivity_exactly() {
        for case in all() {
            if !case.is_mesh() {
                continue;
            }
            let limit = case.points.len() as u32;

            let mut buf = Vec::new();
            write_indices(&mut buf, &case.faces, limit)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));

            let mut cursor = Cursor::new(&buf);
            let recovered: Vec<[u32; 3]> = read_indices(&mut cursor, limit).unwrap();

            assert_eq!(cursor.position() as usize, buf.len(), "{}", case.name);
            assert_eq!(recovered, case.faces, "{}", case.name);
        }
    }

    /// A generator that quietly stops being adversarial makes every later measurement taken
    /// against it meaningless, and nothing else would notice. These assert the defining property
    /// of each awkward case.
    #[test]
    fn the_awkward_cases_are_actually_awkward() {
        let clusters = distant_clusters();
        let ext = clusters.extents();
        assert!(
            ext.iter().all(|&e| e > 300.0),
            "clusters should span a large volume, got {ext:?}"
        );

        let plate = oblique_plate();
        let ext = plate.extents();
        assert!(
            ext.iter().all(|&e| e > 40.0),
            "an oblique plate should be wide on every axis, got {ext:?}"
        );

        let flat = planar();
        assert_eq!(
            flat.extents()[2],
            0.0,
            "planar case must have a zero extent"
        );

        let wide = wide_dynamic_range();
        assert!(
            wide.extents()[0] / wide.tol > 1e6,
            "wide_dynamic_range should have a demanding range/tolerance ratio"
        );

        // Shuffled must be geometrically identical to the surface it came from, so that any size
        // difference between them is attributable to index locality alone.
        let smooth = smooth_surface();
        let shuf = shuffled();
        assert_eq!(smooth.points.len(), shuf.points.len());
        assert_eq!(smooth.faces.len(), shuf.faces.len());
        for i in 0..3 {
            assert!((smooth.extents()[i] - shuf.extents()[i]).abs() < 1e-9);
        }

        // And its index locality must actually be worse, which is the whole point.
        assert!(
            mean_index_spread(&shuf.faces) > 10.0 * mean_index_spread(&smooth.faces),
            "shuffling should destroy index locality"
        );

        let boundary = boundary_heavy();
        assert!(
            boundary.faces.len() < grid_faces(100, 100).len(),
            "boundary_heavy should have holes punched in it"
        );
    }

    /// Mean spread between the smallest and largest index within a face, as a proxy for how much
    /// delta coding will be able to do.
    fn mean_index_spread(faces: &[[u32; 3]]) -> f64 {
        if faces.is_empty() {
            return 0.0;
        }
        let total: u64 = faces
            .iter()
            .map(|f| {
                let max = f.iter().max().unwrap();
                let min = f.iter().min().unwrap();
                u64::from(max - min)
            })
            .sum();
        total as f64 / faces.len() as f64
    }

    /// The flat axis must cost nothing and come back exactly, not merely within tolerance.
    #[test]
    fn the_planar_case_is_recovered_exactly_on_its_flat_axis() {
        let case = planar();
        let mut buf = Vec::new();
        write_points(&mut buf, &case.points, case.tol).unwrap();

        let mut cursor = Cursor::new(&buf);
        let recovered: Vec<[f64; 3]> = read_points(&mut cursor).unwrap();

        for r in &recovered {
            assert_eq!(r[2], 12.5);
        }

        let widths_at = 1 + 4 + 4 + 1 + 8 * 2 * 3;
        assert_eq!(buf[widths_at + 2], 0, "flat axis should take zero bits");
    }

    /// Past the representable limit the encoder must refuse rather than quietly return values
    /// outside the tolerance it promised.
    #[test]
    fn an_impossible_ratio_is_refused() {
        let case = wide_dynamic_range();
        let mut buf = Vec::new();

        assert!(matches!(
            write_points(&mut buf, &case.points, 1e-15),
            Err(Error::ToleranceNotRepresentable { .. })
        ));
    }

    /// Noise should cost real bits: a noisy surface cannot compress as well as a clean one at the
    /// same tolerance, and if it ever does, the noise is not reaching the encoder.
    #[test]
    fn noise_costs_bits() {
        let clean = smooth_surface();
        let noisy = noisy_surface(0.05);

        let mut a = Vec::new();
        let mut b = Vec::new();
        write_points(&mut a, &clean.points, clean.tol).unwrap();
        write_points(&mut b, &noisy.points, noisy.tol).unwrap();

        // The two have the same tolerance and near-identical extents, so a single bounding box
        // charges them identically and this pair used to come out byte for byte the same. What
        // separates them now is partitioning: noise pushes points out of the tight local boxes a
        // clean surface falls into, so fewer cuts pay and the ones that do save less.
        assert!(
            a.len() < b.len(),
            "noise should cost bits: clean {} against noisy {}",
            a.len(),
            b.len()
        );
    }
}
