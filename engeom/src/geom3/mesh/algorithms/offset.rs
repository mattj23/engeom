//! This module offsets a triangle mesh, moving each point as if the faces it is attached to were
//! offset by their face normal direction by a given distance.  It works using the point and face
//! buffers alone, and for points that belong to three or fewer faces it provides an exact solution.
//! For points with more than three faces it provides a reasonable estimate, unless the faces are
//! symmetric about the point in which case it's the exact solution.
//!
//! # How it works
//!
//! Write the new position of a point as `p + u`. Every incident face plane, moved out by `d`, is
//! the constraint `n_i . u = d`. For a point where `k` faces meet, that is `k` linear equations in
//! three unknowns:
//!
//! - `k = 3` with independent normals is square and has an exact unique solution.
//! - `k > 3` is overdetermined, and is only consistent when the incident normals lie on a common
//!   cone. Symmetric features satisfy this exactly: a regular pyramid apex with four faces solves
//!   with zero residual. Generic geometry does not, and the least squares answer is the best a
//!   single displacement can do.
//! - `k < 3` is underdetermined: a point in a flat region constrains only one direction, a point on
//!   a straight edge only two. The minimum norm solution is the right answer in both cases, moving
//!   the point along the normal and along the edge bisector respectively without wandering
//!   sideways.
//!
//! Minimizing `sum_i w_i (n_i . u - d)^2` gives the normal equations
//! `(sum_i w_i n_i n_i^T) u = d (sum_i w_i n_i)`. The matrix on the left is the Garland-Heckbert
//! quadric, the same one behind quadric error simplification and dual contouring, and it
//! accumulates in a single pass over the faces with no adjacency structure needed. Its rank is
//! exactly the local feature dimension: 3 at a corner, 2 on an edge, 1 on a flat region, which is
//! why a pseudoinverse handles all three cases without any of them being special-cased.
//!
//! Weights are interior angles, for the same reason `compute_point_normals` uses them: a flat
//! quadrilateral split into two triangles must not vote twice at the corners the diagonal touches.
//!
//! Things to be aware of when usinig this:
//!
//! - Self intersections aren't detected; if you offset inward on a convex region or outward on a
//!   concave region by more than the local radius of curvature, the surface will end up folding
//!   through itself.  That's not detected.
//! - The topology isn't changed.  A mathematically complete offset would grow new faces at convex
//!   edges and remove them at concave ones, but this is just moving existing points.

use crate::geom3::mesh::algorithms::normals::compute_face_normal;
use crate::na::Matrix3;
use crate::{Point3, Result, UnitVec3, Vector3};

/// Options controlling how the offset handles the two ill-conditioned cases.
///
/// Marked `#[non_exhaustive]` with a `Default`, so construct it as
/// `OffsetOpts { max_ratio: 4.0, ..Default::default() }` and new fields will not break your code.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OffsetOpts {
    /// The relative singular value below which a direction is treated as unconstrained, as a
    /// fraction of the largest singular value at that point.
    ///
    /// This is what decides whether a point sits on a genuine corner or on a smooth surface. The
    /// small singular values of the quadric scale as the square of the angular spread of the
    /// incident normals, so a cutoff of `t` treats a spread below roughly `sqrt(t)` radians as
    /// flat. The default of `1e-6` puts that boundary at about 0.06 degrees.
    ///
    /// Raising it makes gently curved regions offset along their average normal, which is smoother
    /// but under-offsets shallow features. Lowering it lets nearly-flat vertices behave as corners,
    /// where the solution is real but poorly conditioned, and the displacement grows quickly.
    pub rank_tolerance: f64,

    /// The largest displacement any point may take, as a multiple of the offset distance.
    ///
    /// A displacement of `d/cos(theta)` diverges as the half-angle of the cone of incident normals
    /// approaches 90 degrees, so a near-degenerate spike in the mesh would otherwise be flung an
    /// arbitrary distance. A point whose solution exceeds this is scaled back to it, keeping its
    /// direction. The default of `10.0` corresponds to a cone half-angle of about 84 degrees, which
    /// no well-formed surface reaches.
    ///
    /// Note that `rank_tolerance` already puts a ceiling of roughly `1/sqrt(rank_tolerance)` on the
    /// amplification, since a feature sharp enough to exceed that has an unresolvable axis and
    /// falls back to the average normal. This cap is the one that does the useful work, catching
    /// spikes well before they reach that fallback, where the result would be arbitrary rather than
    /// merely large.
    pub max_ratio: f64,
}

impl Default for OffsetOpts {
    fn default() -> Self {
        Self {
            rank_tolerance: 1.0e-6,
            max_ratio: 10.0,
        }
    }
}

/// Compute the offset positions of every point of a mesh, moving the surface `distance` along its
/// own normal direction.
///
/// Each point is placed at the least squares intersection of its incident face planes, each moved
/// out by `distance`. This is the geometrically correct offset: a convex corner travels farther
/// than `distance` so that every face touching it lands exactly `distance` from where it was. See
/// the module documentation for why displacing along point normals instead is wrong, and for the
/// two things this deliberately does not handle.
///
/// A positive distance moves the surface along its normals, which for a closed mesh wound outward
/// means outward. A negative distance moves it the other way.
///
/// Triangles with no well-defined normal contribute nothing and are skipped.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
/// * `distance`: how far to move the surface, in the mesh's length units
/// * `opts`: how to handle unconstrained directions and near-degenerate spikes
///
/// returns: `Result<Vec<Point3>>`, one position per point, in the same order. Fails if a point
/// belongs to no face or only to degenerate ones, since there is nothing to offset it against
pub fn compute_face_offset_points(
    points: &[Point3],
    faces: &[[u32; 3]],
    distance: f64,
    opts: &OffsetOpts,
) -> Result<Vec<Point3>> {
    if !distance.is_finite() {
        return Err(format!("The offset distance must be finite, but it was {distance}").into());
    }
    if !(opts.rank_tolerance > 0.0 && opts.rank_tolerance < 1.0) {
        return Err("The rank tolerance must lie strictly between 0 and 1".into());
    }
    if !(opts.max_ratio >= 1.0 && opts.max_ratio.is_finite()) {
        return Err("The maximum displacement ratio must be finite and at least 1".into());
    }

    let quadrics = accumulate_quadrics(points, faces)?;

    let limit = opts.max_ratio * distance.abs();
    let mut out = Vec::with_capacity(points.len());

    for (i, q) in quadrics.iter().enumerate() {
        if !q.has_faces {
            return Err(format!(
                "Point {i} cannot be offset. It belongs to no face or only to degenerate faces, so \
                 there are no surfaces to offset it away from."
            )
            .into());
        }

        let mut u = solve_min_norm(&q.a, &(q.b * distance), opts.rank_tolerance);

        // A cone of incident normals approaching a half-angle of 90 degrees sends the exact
        // solution to infinity, so cap the travel while keeping its direction.
        let norm = u.norm();
        if norm > limit && norm > 0.0 {
            u *= limit / norm;
        }

        out.push(points[i] + u);
    }

    Ok(out)
}

/// Displace every point along its own point normal, which is **not** a surface offset.
///
/// This under-offsets every convex feature and over-offsets every concave one, by the cosine of the
/// angle between the point normal and each incident face normal: a cube corner ends up at
/// `distance/sqrt(3)` from where the three offset planes actually meet. It is here because it is
/// cheap, because it never produces the large displacements a true offset does at a sharp feature,
/// and because it is what some callers actually mean when they want to nudge a surface. It is named
/// so that nobody reaches for it expecting a real offset.
///
/// Use `compute_offset_points` unless you specifically want this.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `normals`: one unit normal per point, in the same order, as from `compute_point_normals`
/// * `distance`: how far to move each point along its normal
///
/// returns: `Result<Vec<Point3>>`
pub fn compute_normal_displaced_points(
    points: &[Point3],
    normals: &[UnitVec3],
    distance: f64,
) -> Result<Vec<Point3>> {
    if points.len() != normals.len() {
        return Err(format!(
            "There are {} points but {} normals, and displacement needs one normal per point",
            points.len(),
            normals.len()
        )
        .into());
    }
    if !distance.is_finite() {
        return Err(format!("The offset distance must be finite, but it was {distance}").into());
    }

    Ok(points
        .iter()
        .zip(normals.iter())
        .map(|(p, n)| *p + n.into_inner() * distance)
        .collect())
}

/// The accumulated least squares system at a single point.
struct Quadric {
    /// `sum_i w_i n_i n_i^T`, the Garland-Heckbert quadric.
    a: Matrix3<f64>,

    /// `sum_i w_i n_i`, which becomes the right hand side once scaled by the distance.
    b: Vector3,

    /// Whether any non-degenerate face touched this point. A point with none has an all-zero
    /// system, which is indistinguishable from a solved one without this flag.
    has_faces: bool,
}

/// Accumulate the per-point least squares systems in a single pass over the faces.
fn accumulate_quadrics(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<Quadric>> {
    let mut out = (0..points.len())
        .map(|_| Quadric {
            a: Matrix3::zeros(),
            b: Vector3::zeros(),
            has_faces: false,
        })
        .collect::<Vec<_>>();

    for face in faces.iter() {
        let p = [
            *points
                .get(face[0] as usize)
                .ok_or("Face refers to a point which does not exist")?,
            *points
                .get(face[1] as usize)
                .ok_or("Face refers to a point which does not exist")?,
            *points
                .get(face[2] as usize)
                .ok_or("Face refers to a point which does not exist")?,
        ];

        let Some(normal) = compute_face_normal(&p) else {
            continue;
        };
        let n = normal.into_inner();
        let outer = n * n.transpose();

        for k in 0..3 {
            let e1 = p[(k + 1) % 3] - p[k];
            let e2 = p[(k + 2) % 3] - p[k];

            // The interior angle at this corner, which is the weight. Angle weighting makes the
            // result independent of how a flat polygon was split into triangles.
            let w = e1.cross(&e2).norm().atan2(e1.dot(&e2));

            let q = &mut out[face[k] as usize];
            q.a += outer * w;
            q.b += n * w;
            q.has_faces = true;
        }
    }

    Ok(out)
}

/// Solve `a u = c` for the minimum norm `u`, treating singular values below `tol` times the largest
/// as zero.
///
/// The matrix is symmetric positive semi-definite by construction, and its rank is the dimension of
/// the local feature: 3 at a corner, 2 on an edge, 1 on a flat region. Discarding the null space
/// rather than trying to invert through it is what makes the underdetermined cases come out right,
/// since the minimum norm solution is the one that does not move the point in a direction nothing
/// asked it to move.
fn solve_min_norm(a: &Matrix3<f64>, c: &Vector3, tol: f64) -> Vector3 {
    let svd = a.svd(true, true);
    let (Some(u), Some(v_t)) = (svd.u, svd.v_t) else {
        return Vector3::zeros();
    };

    let largest = svd
        .singular_values
        .iter()
        .fold(0.0_f64, |acc, s| acc.max(*s));
    if largest <= 0.0 {
        return Vector3::zeros();
    }

    let cutoff = largest * tol;
    let mut result = Vector3::zeros();
    for i in 0..3 {
        let s = svd.singular_values[i];
        if s <= cutoff {
            continue;
        }
        // Column i of U and row i of V^T, which is column i of V.
        let projection = u.column(i).dot(c) / s;
        result += v_t.row(i).transpose() * projection;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::mesh::algorithms::normals::compute_point_normals;
    use approx::assert_relative_eq;

    /// A cube of side 2 centred on the origin, with each square face split into two triangles.
    fn cube() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut points = Vec::new();
        for x in [-1.0, 1.0] {
            for y in [-1.0, 1.0] {
                for z in [-1.0, 1.0] {
                    points.push(Point3::new(x, y, z));
                }
            }
        }

        let quads = [
            [1, 5, 7, 3], // +z
            [0, 2, 6, 4], // -z
            [0, 4, 5, 1], // -y
            [2, 3, 7, 6], // +y
            [4, 6, 7, 5], // +x
            [0, 1, 3, 2], // -x
        ];

        let mut faces = Vec::new();
        for q in quads {
            faces.push([q[0], q[1], q[2]]);
            faces.push([q[0], q[2], q[3]]);
        }

        (points, faces)
    }

    /// A grid of points in the xy plane, split into triangles.
    fn flat_grid(n: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut points = Vec::new();
        for j in 0..n {
            for i in 0..n {
                points.push(Point3::new(i as f64, j as f64, 0.0));
            }
        }

        let mut faces = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let a = j * n + i;
                faces.push([a, a + 1, a + n + 1]);
                faces.push([a, a + n + 1, a + n]);
            }
        }

        (points, faces)
    }

    /// The defining property of a correct offset: every face of the offset mesh lies exactly
    /// `distance` from the plane of the face it came from. Returns the worst violation.
    fn worst_plane_error(
        points: &[Point3],
        faces: &[[u32; 3]],
        moved: &[Point3],
        distance: f64,
    ) -> f64 {
        let mut worst: f64 = 0.0;
        for face in faces {
            let p = face.map(|i| points[i as usize]);
            let Some(n) = compute_face_normal(&p) else {
                continue;
            };

            for i in face {
                let u = moved[*i as usize] - points[*i as usize];
                worst = worst.max((n.dot(&u) - distance).abs());
            }
        }
        worst
    }

    /// A cube corner has three independent normals, so the system is square and exact. The corner
    /// must travel `d*sqrt(3)`, not `d`.
    #[test]
    fn a_cube_corner_reaches_the_true_plane_intersection() -> Result<()> {
        let (points, faces) = cube();
        let d = 0.25;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        // The offset of a cube of half-side 1 is a cube of half-side 1 + d.
        for (before, after) in points.iter().zip(moved.iter()) {
            let expected = before + before.coords.map(|c| c.signum()) * d;
            assert_relative_eq!(after.coords, expected.coords, epsilon = 1.0e-12);
        }

        // Which means each corner moved d*sqrt(3), well past the d that displacing along the point
        // normal would have given.
        for (before, after) in points.iter().zip(moved.iter()) {
            assert_relative_eq!(
                (after - before).norm(),
                d * 3.0_f64.sqrt(),
                epsilon = 1.0e-12
            );
        }

        assert_relative_eq!(
            worst_plane_error(&points, &faces, &moved, d),
            0.0,
            epsilon = 1.0e-12
        );

        Ok(())
    }

    /// The comparison that justifies the whole module. Displacing along point normals moves every
    /// cube face by only `d/sqrt(3)`.
    #[test]
    fn displacing_along_point_normals_under_offsets_a_corner() -> Result<()> {
        let (points, faces) = cube();
        let d = 0.25;

        let normals = compute_point_normals(&points, &faces)?;
        let naive = compute_normal_displaced_points(&points, &normals, d)?;

        // Every face ends up short by the cosine between the point normal and the face normal,
        // which on a cube is 1/sqrt(3).
        let shortfall = d - d / 3.0_f64.sqrt();
        assert_relative_eq!(
            worst_plane_error(&points, &faces, &naive, d),
            shortfall,
            epsilon = 1.0e-12
        );

        // That is a 42% error, not a rounding difference.
        assert!(shortfall / d > 0.42);

        Ok(())
    }

    /// A flat region constrains only one direction, so the system has rank 1 and the minimum norm
    /// solution must move the point along the normal and nowhere else.
    #[test]
    fn a_flat_sheet_offsets_along_its_normal() -> Result<()> {
        let (points, faces) = flat_grid(4);
        let d = 2.5;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        for (before, after) in points.iter().zip(moved.iter()) {
            assert_relative_eq!(
                after.coords,
                before.coords + Vector3::z() * d,
                epsilon = 1.0e-12
            );
        }

        Ok(())
    }

    /// A straight edge constrains two directions, so the system has rank 2 and the point must
    /// travel along the bisector by `d/cos(half angle)`, staying in the plane of the two normals.
    #[test]
    fn a_valley_floor_offsets_along_the_bisector() -> Result<()> {
        // Two flat strips meeting along the y axis at 90 degrees, wound so that both normals point
        // upward: a valley seen from above, with its floor along the y axis.
        let points = vec![
            Point3::new(-1.0, 0.0, 1.0),
            Point3::new(-1.0, 1.0, 1.0),
            Point3::new(0.0, 0.0, 0.0), // on the floor
            Point3::new(0.0, 1.0, 0.0), // on the floor
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
        ];
        let faces = vec![[0, 2, 3], [0, 3, 1], [2, 4, 5], [2, 5, 3]];

        let d = 0.5;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        // The two normals are (1,0,1)/sqrt(2) and (-1,0,1)/sqrt(2), 90 degrees apart, so their
        // bisector is +z and the floor travels d/cos(45) = d*sqrt(2) along it. Nothing along x or
        // y, which is what the minimum norm solution buys: the rank 2 system says nothing about
        // the direction along the valley, and the point must not wander down it.
        for i in [2, 3] {
            let u = moved[i] - points[i];
            assert_relative_eq!(u.x, 0.0, epsilon = 1.0e-12);
            assert_relative_eq!(u.y, 0.0, epsilon = 1.0e-12);
            assert_relative_eq!(u.z, d * 2.0_f64.sqrt(), epsilon = 1.0e-12);
        }

        assert_relative_eq!(
            worst_plane_error(&points, &faces, &moved, d),
            0.0,
            epsilon = 1.0e-12
        );

        Ok(())
    }

    /// Four faces at a point is overdetermined, but a regular pyramid's normals lie on a common
    /// cone so the system stays consistent and solves exactly.
    #[test]
    fn a_regular_pyramid_apex_solves_exactly() -> Result<()> {
        let points = vec![
            Point3::new(0.0, 0.0, 1.0), // apex
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];

        let d = 0.3;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        // Only the apex is checked for exactness: the base points touch two side faces and no
        // base face, which is a genuine edge rather than a symmetric cone.
        let mut worst: f64 = 0.0;
        for face in faces.iter() {
            let p = face.map(|i| points[i as usize]);
            let n = compute_face_normal(&p).unwrap();
            let u = moved[0] - points[0];
            worst = worst.max((n.dot(&u) - d).abs());
        }
        assert_relative_eq!(worst, 0.0, epsilon = 1.0e-12);

        // By symmetry the apex travels straight up.
        let u = moved[0] - points[0];
        assert_relative_eq!(u.x, 0.0, epsilon = 1.0e-12);
        assert_relative_eq!(u.y, 0.0, epsilon = 1.0e-12);
        assert!(u.z > d, "the apex should travel farther than the distance");

        Ok(())
    }

    /// A sphere is the case where the answer is known in closed form: offsetting by `d` should
    /// grow the radius by `d` everywhere.
    #[test]
    fn a_sphere_grows_by_the_offset_distance() -> Result<()> {
        let sphere = crate::Mesh::create_sphere(10.0, 64, 64);
        let d = 1.5;
        let moved =
            compute_face_offset_points(sphere.vertices(), sphere.faces(), d, &OffsetOpts::default())?;

        let mut worst: f64 = 0.0;
        for p in moved.iter() {
            worst = worst.max((p.coords.norm() - 11.5).abs());
        }

        // The discretization means the faceted surface sits slightly inside the true sphere, so
        // this cannot be exact, but it should be far tighter than the facet size.
        assert!(worst < 0.02, "worst radius error was {worst}");

        Ok(())
    }

    #[test]
    fn a_negative_distance_offsets_the_other_way() -> Result<()> {
        let (points, faces) = cube();
        let d = -0.25;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        // The cube shrinks to half-side 0.75.
        for after in moved.iter() {
            assert_relative_eq!(after.coords.abs().max(), 0.75, epsilon = 1.0e-12);
        }

        assert_relative_eq!(
            worst_plane_error(&points, &faces, &moved, d),
            0.0,
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn a_zero_distance_leaves_every_point_alone() -> Result<()> {
        let (points, faces) = cube();
        let moved = compute_face_offset_points(&points, &faces, 0.0, &OffsetOpts::default())?;

        for (before, after) in points.iter().zip(moved.iter()) {
            assert_relative_eq!(after.coords, before.coords, epsilon = 1.0e-15);
        }

        Ok(())
    }

    /// The whole computation is linear in the mesh scale, and the degeneracy and rank tests are
    /// both relative, so a mesh in microns must offset the same way as one in meters.
    #[test]
    fn scale_does_not_change_the_result() -> Result<()> {
        let (big, faces) = cube();
        let small: Vec<Point3> = big.iter().map(|p| *p * 1.0e-6).collect();

        let a = compute_face_offset_points(&big, &faces, 0.25, &OffsetOpts::default())?;
        let b = compute_face_offset_points(&small, &faces, 0.25e-6, &OffsetOpts::default())?;

        for (x, y) in a.iter().zip(b.iter()) {
            assert_relative_eq!(x.coords * 1.0e-6, y.coords, epsilon = 1.0e-15);
        }

        Ok(())
    }

    /// A point with nothing to offset against used to be left silently in place, which produced a
    /// mesh that was offset almost everywhere.
    #[test]
    fn a_point_with_no_faces_is_an_error() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(5.0, 5.0, 5.0), // orphan
        ];
        let faces = vec![[0, 1, 2]];

        let err = compute_face_offset_points(&points, &faces, 1.0, &OffsetOpts::default())
            .unwrap_err()
            .to_string();
        assert!(err.contains("Point 3"), "{err}");
    }

    #[test]
    fn an_out_of_range_face_index_is_an_error() {
        let (points, mut faces) = flat_grid(3);
        faces.push([0, 1, 99]);

        assert!(compute_face_offset_points(&points, &faces, 1.0, &OffsetOpts::default()).is_err());
    }

    /// A tall thin needle: three faces whose normals are nearly perpendicular to its axis, so the
    /// exact solution runs away and the cap has to catch it.
    fn needle(base_radius: f64, height: f64) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let r = base_radius;
        let s = 3.0_f64.sqrt() / 2.0;
        let points = vec![
            Point3::new(0.0, 0.0, height), // the tip
            Point3::new(r, 0.0, 0.0),
            Point3::new(-r * 0.5, r * s, 0.0),
            Point3::new(-r * 0.5, -r * s, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3], [0, 3, 1]];

        (points, faces)
    }

    #[test]
    fn a_needle_tip_is_clamped() -> Result<()> {
        // The inradius of the equilateral base is r/2, so the faces sit at an angle whose sine is
        // (r/2)/h, and the tip wants to travel d divided by that. At r = 0.04 and h = 1 that is
        // fifty times the offset distance.
        let (points, faces) = needle(0.04, 1.0);
        let d = 1.0;

        let loose = OffsetOpts {
            max_ratio: 1.0e6,
            ..Default::default()
        };
        let unclamped = compute_face_offset_points(&points, &faces, d, &loose)?;
        let far = (unclamped[0] - points[0]).norm();
        assert!(
            (40.0..60.0).contains(&far),
            "expected roughly 50x the distance, got {far}"
        );

        // With the default cap it goes exactly as far as it is allowed to, and no farther.
        let clamped = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;
        assert_relative_eq!(
            (clamped[0] - points[0]).norm(),
            OffsetOpts::default().max_ratio * d,
            epsilon = 1.0e-9
        );

        // The direction is preserved, only the magnitude is cut.
        let a = (unclamped[0] - points[0]).normalize();
        let b = (clamped[0] - points[0]).normalize();
        assert_relative_eq!(a, b, epsilon = 1.0e-9);

        Ok(())
    }

    /// The rank cutoff puts its own ceiling on the amplification, well before a needle gets thin
    /// enough to overflow anything: once the incident normals deviate from coplanar by less than
    /// roughly `sqrt(rank_tolerance)`, the axis reads as unconstrained and the tip simply offsets
    /// along the average normal instead.
    #[test]
    fn an_extremely_thin_needle_falls_back_to_the_average_normal() -> Result<()> {
        // A million to one aspect ratio, far past where the quadric can resolve the axis.
        let (points, faces) = needle(1.0e-6, 1.0);
        let d = 1.0;

        let loose = OffsetOpts {
            max_ratio: 1.0e12,
            ..Default::default()
        };
        let moved = compute_face_offset_points(&points, &faces, d, &loose)?;

        // Bounded by the rank cutoff at roughly 1/sqrt(tolerance), not by max_ratio.
        let travel = (moved[0] - points[0]).norm();
        assert!(
            travel < 1.0 / OffsetOpts::default().rank_tolerance.sqrt(),
            "travel of {travel} was not bounded by the rank cutoff"
        );

        Ok(())
    }

    /// The rank tolerance is what separates a corner from a smooth surface. A vertex whose normals
    /// deviate by far less than the cutoff must behave as if it were flat.
    #[test]
    fn a_nearly_flat_vertex_is_treated_as_flat() -> Result<()> {
        // A grid with its centre lifted by a nanometre, which is a corner in principle and a flat
        // sheet in every way that matters.
        let (mut points, faces) = flat_grid(3);
        points[4].z = 1.0e-9;

        let d = 1.0;
        let moved = compute_face_offset_points(&points, &faces, d, &OffsetOpts::default())?;

        // Without the cutoff the ill-conditioned solve at the centre would throw it a long way.
        let u = moved[4] - points[4];
        assert_relative_eq!(u.x, 0.0, epsilon = 1.0e-6);
        assert_relative_eq!(u.y, 0.0, epsilon = 1.0e-6);
        assert_relative_eq!(u.z, d, epsilon = 1.0e-4);

        Ok(())
    }

    #[test]
    fn degenerate_faces_are_skipped() -> Result<()> {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3], [0, 1, 1]];

        let moved = compute_face_offset_points(&points, &faces, 1.0, &OffsetOpts::default())?;
        for (before, after) in points.iter().zip(moved.iter()) {
            assert_relative_eq!(
                after.coords,
                before.coords + Vector3::z(),
                epsilon = 1.0e-12
            );
        }

        Ok(())
    }

    #[test]
    fn invalid_arguments_are_rejected() {
        let (points, faces) = cube();

        assert!(compute_face_offset_points(&points, &faces, f64::NAN, &OffsetOpts::default()).is_err());
        assert!(
            compute_face_offset_points(
                &points,
                &faces,
                1.0,
                &OffsetOpts {
                    rank_tolerance: 0.0,
                    ..Default::default()
                }
            )
            .is_err()
        );
        assert!(
            compute_face_offset_points(
                &points,
                &faces,
                1.0,
                &OffsetOpts {
                    max_ratio: 0.5,
                    ..Default::default()
                }
            )
            .is_err()
        );

        let normals = compute_point_normals(&points, &faces).unwrap();
        assert!(compute_normal_displaced_points(&points, &normals[..4], 1.0).is_err());
    }
}
