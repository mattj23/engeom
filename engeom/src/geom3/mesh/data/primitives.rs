//! This module contains the constructors which build a `MeshData3` for a simple geometric shape,
//! along with the sample meshes embedded in the binary.
//!
//! These constructors live on the data container rather than the accelerated one because a
//! tessellation produces only points and faces. `Mesh3` provides a pass-through for each
//! constructor and adds the `is_solid` flag, which controls query behavior rather than the
//! geometry itself.
//!
//! # Conventions
//!
//! Every shape is built at the origin. A body of revolution puts its axis on `z`, so a cylinder
//! and a cone run along `z`, a sphere's poles are at `+/-z`, and a flat circle's normal is `+z`.
//! The two shapes that span a pair of points first build their canonical form and then transform it
//! onto the segment.
//!
//! Curvature is tessellated to a chordal tolerance rather than a fixed division count. The points
//! lie on the true surface, and the facets sag inward by at most `tol`. The
//! `arc_segments_for_tol` function calculates the division count and applies the minimum segment
//! count. A shape curved in two directions splits the tolerance equally between them because the
//! two deviations are additive. The ruled surface of a cone contributes no deviation, so only its
//! base circle consumes the tolerance.
//!
//! Closed shapes are watertight and manifold: a seam ring is one set of shared vertices, not two
//! coincident ones, and every face is wound counterclockwise as viewed from outside, so the signed
//! volume is positive. None of the constructors attach attributes; a tessellation produces only
//! points and faces.

use super::MeshData3;
use crate::common::arc_segments_for_tol;
use crate::geom3::IsoExtensions3;
use crate::io::read_tc_mesh_from;
use crate::{Iso3, Point3, Result, Vector3};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

// ===============================================================================================
// Tessellation helpers
// ===============================================================================================

/// Append `n` evenly spaced points around a circle of `radius` in the plane at `z`, starting at
/// +X and proceeding counterclockwise as viewed from +Z. Returns the index of the first point
/// appended.
fn push_ring(points: &mut Vec<Point3>, radius: f64, z: f64, n: usize) -> u32 {
    let start = points.len() as u32;
    for i in 0..n {
        let angle = TAU * (i as f64) / (n as f64);
        points.push(Point3::new(radius * angle.cos(), radius * angle.sin(), z));
    }
    start
}

/// Append the 2n triangles that join two rings of `n` points each, both wound in the same
/// direction and starting at the same angle. The winding is counterclockwise as viewed from
/// outside a surface that runs from the `lower` ring up to the `upper` one.
fn band_faces(faces: &mut Vec<[u32; 3]>, lower: u32, upper: u32, n: usize) {
    let n = n as u32;
    for i in 0..n {
        let j = (i + 1) % n;
        faces.push([lower + i, lower + j, upper + j]);
        faces.push([lower + i, upper + j, upper + i]);
    }
}

/// Append the `n` triangles that fan a ring of `n` points to a single center point. When `up` is
/// `true`, the fan is wound counterclockwise as viewed from +Z, which faces outward for the top cap
/// of a body. When `up` is `false`, the winding is reversed for a bottom cap.
fn fan_faces(faces: &mut Vec<[u32; 3]>, center: u32, ring: u32, n: usize, up: bool) {
    let n = n as u32;
    for i in 0..n {
        let j = (i + 1) % n;
        if up {
            faces.push([center, ring + i, ring + j]);
        } else {
            faces.push([center, ring + j, ring + i]);
        }
    }
}

// ===============================================================================================
// Shapes centered at the origin
// ===============================================================================================

impl MeshData3 {
    /// Create a box mesh with the given dimensions, centered at the origin.
    ///
    /// # Arguments
    ///
    /// * `length`: the dimension of the box in the x direction
    /// * `width`: the dimension of the box in the y direction
    /// * `height`: the dimension of the box in the z direction
    ///
    /// returns: `MeshData3`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let mesh = MeshData3::create_box(2.0, 4.0, 6.0);
    ///
    /// for p in mesh.points() {
    ///     assert_relative_eq!(p.x.abs(), 1.0);
    ///     assert_relative_eq!(p.y.abs(), 2.0);
    ///     assert_relative_eq!(p.z.abs(), 3.0);
    /// }
    /// ```
    pub fn create_box(length: f64, width: f64, height: f64) -> Self {
        let hx = length / 2.0;
        let hy = width / 2.0;
        let hz = height / 2.0;

        // Two rings of four corners, both wound counterclockwise as viewed from +Z, with the
        // bottom ring first. All eight vertices are shared between their three faces, so the box
        // is watertight.
        let points = vec![
            Point3::new(-hx, -hy, -hz),
            Point3::new(hx, -hy, -hz),
            Point3::new(hx, hy, -hz),
            Point3::new(-hx, hy, -hz),
            Point3::new(-hx, -hy, hz),
            Point3::new(hx, -hy, hz),
            Point3::new(hx, hy, hz),
            Point3::new(-hx, hy, hz),
        ];

        // Each face pair is wound counterclockwise as viewed from outside the box.
        let faces = vec![
            [0, 2, 1], // -z
            [0, 3, 2],
            [4, 5, 6], // +z
            [4, 6, 7],
            [0, 1, 5], // -y
            [0, 5, 4],
            [1, 2, 6], // +x
            [1, 6, 5],
            [2, 3, 7], // +y
            [2, 7, 6],
            [3, 0, 4], // -x
            [3, 4, 7],
        ];

        Self::new(points, faces).expect("A box tessellation is always a valid mesh")
    }

    /// Create a closed spherical mesh centered at the origin, with its poles on the local `z`
    /// axis. This is a UV sphere: rings of evenly spaced points at evenly spaced polar angles,
    /// with a triangle fan closing each pole.
    ///
    /// Every point lies on the true sphere, so the facets sag inward. A facet's deviation comes
    /// from curvature in both the azimuthal and polar directions. Because these contributions are
    /// additive, each is held to half of `tol`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the sphere, which must be positive and finite
    /// * `tol`: the maximum allowed deviation of a facet from the true sphere, which must be
    ///   positive and finite. No matter how loose this is, a full turn around the equator never
    ///   gets fewer than 8 segments and a pole-to-pole sweep never gets fewer than 4.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let sphere = MeshData3::create_sphere(1.0, 1.0e-3).unwrap();
    ///
    /// // Verify that the points are on the surface of the unit sphere.
    /// for p in sphere.points() {
    ///     assert_relative_eq!(p.coords.norm(), 1.0)
    /// }
    /// ```
    pub fn create_sphere(radius: f64, tol: f64) -> Result<Self> {
        let n_theta = arc_segments_for_tol(radius, TAU, tol / 2.0)?;
        let n_phi = arc_segments_for_tol(radius, PI, tol / 2.0)?;

        let mut points = Vec::with_capacity(n_theta * (n_phi - 1) + 2);
        let mut faces = Vec::with_capacity(2 * n_theta * (n_phi - 1));

        let north = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, radius));

        // The rings proceed down from the north pole, so the second ring in each consecutive pair
        // is the lower one.
        let mut rings = Vec::with_capacity(n_phi - 1);
        for k in 1..n_phi {
            let phi = PI * (k as f64) / (n_phi as f64);
            let (sin, cos) = phi.sin_cos();
            rings.push(push_ring(&mut points, radius * sin, radius * cos, n_theta));
        }

        let south = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, -radius));

        fan_faces(&mut faces, north, rings[0], n_theta, true);
        for w in rings.windows(2) {
            band_faces(&mut faces, w[1], w[0], n_theta);
        }
        fan_faces(&mut faces, south, rings[rings.len() - 1], n_theta, false);

        Ok(Self::new(points, faces).expect("A sphere tessellation is always a valid mesh"))
    }

    /// Create a closed cylindrical mesh centered at the origin and aligned with the local `z`
    /// axis. The wall points lie on the true cylinder, so the chords sag inward by at most `tol`,
    /// and each flat cap is a triangle fan from its own center point.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the cylinder, which must be positive and finite
    /// * `height`: full height of the cylinder, along the z-axis
    /// * `tol`: the maximum allowed chordal deviation of the wall, which must be positive and
    ///   finite. The circumference never gets fewer than 8 segments, no matter how loose this is.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let cyl = MeshData3::create_cylinder(1.0, 4.0, 1.0e-3).unwrap();
    ///
    /// for p in cyl.points() {
    ///     // Every point lies on one of the two cap planes, either out on the wall or at the
    ///     // center of the cap.
    ///     assert_relative_eq!(p.z.abs(), 2.0);
    ///
    ///     let r = p.x.hypot(p.y);
    ///     assert!(r == 0.0 || (r - 1.0).abs() < 1.0e-12);
    /// }
    /// ```
    pub fn create_cylinder(radius: f64, height: f64, tol: f64) -> Result<Self> {
        let n = arc_segments_for_tol(radius, TAU, tol)?;
        let hz = height / 2.0;

        let mut points = Vec::with_capacity(2 * n + 2);
        let mut faces = Vec::with_capacity(4 * n);

        // The wall's two rings are shared with the caps, so the result is watertight.
        let bottom = push_ring(&mut points, radius, -hz, n);
        let top = push_ring(&mut points, radius, hz, n);
        band_faces(&mut faces, bottom, top, n);

        let bottom_center = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, -hz));
        fan_faces(&mut faces, bottom_center, bottom, n, false);

        let top_center = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, hz));
        fan_faces(&mut faces, top_center, top, n, true);

        Ok(Self::new(points, faces).expect("A cylinder tessellation is always a valid mesh"))
    }

    /// Create a closed conical mesh centered at the origin and aligned with the local `z` axis,
    /// with its apex at `+height/2` and its base at `-height/2`.
    ///
    /// The lateral surface is ruled, so only the base circle carries curvature and `tol` is the
    /// chordal deviation of that circle.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the base of the cone, which must be positive and finite
    /// * `height`: full height of the cone, along the z-axis
    /// * `tol`: the maximum allowed chordal deviation of the base circle, which must be positive
    ///   and finite. The base never gets fewer than 8 segments, no matter how loose this is.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let cone = MeshData3::create_cone(2.0, 10.0, 1.0e-3).unwrap();
    ///
    /// // The apex sits at +height/2 and the base ring at -height/2, at the given radius.
    /// let apex = cone.points().iter().map(|p| p.z).fold(f64::MIN, f64::max);
    /// assert_relative_eq!(apex, 5.0);
    ///
    /// let widest = cone.points().iter().map(|p| p.x.hypot(p.y)).fold(f64::MIN, f64::max);
    /// assert_relative_eq!(widest, 2.0);
    /// ```
    pub fn create_cone(radius: f64, height: f64, tol: f64) -> Result<Self> {
        let n = arc_segments_for_tol(radius, TAU, tol)?;
        let hz = height / 2.0;

        let mut points = Vec::with_capacity(n + 2);
        let mut faces = Vec::with_capacity(2 * n);

        // The base ring is shared between the lateral surface and the base fan.
        let base = push_ring(&mut points, radius, -hz, n);

        let apex = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, hz));
        for i in 0..n as u32 {
            let j = (i + 1) % n as u32;
            faces.push([base + i, base + j, apex]);
        }

        let base_center = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, -hz));
        fan_faces(&mut faces, base_center, base, n, false);

        Ok(Self::new(points, faces).expect("A cone tessellation is always a valid mesh"))
    }

    /// Create a flat, filled circle mesh lying in the XY plane, centered at the origin, with the
    /// normal pointing along +Z. The mesh is a triangle fan from the center to evenly spaced
    /// perimeter points, which lie on the true circle so the chords sag inward by at most `tol`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the circle, which must be positive and finite
    /// * `tol`: the maximum allowed chordal deviation of the perimeter, which must be positive
    ///   and finite. A full circle never gets fewer than 8 segments, no matter how loose this is.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    ///
    /// let circle = MeshData3::create_circle(1.0, 1.0e-3).unwrap();
    /// assert_eq!(circle.point_count(), circle.face_count() + 1); // center + perimeter
    /// assert!(circle.face_count() >= 8);
    /// ```
    pub fn create_circle(radius: f64, tol: f64) -> Result<Self> {
        let n = arc_segments_for_tol(radius, TAU, tol)?;
        let mut points = Vec::with_capacity(n + 1);
        let mut faces = Vec::with_capacity(n);

        points.push(Point3::origin());
        let ring = push_ring(&mut points, radius, 0.0, n);
        fan_faces(&mut faces, 0, ring, n, true);

        Ok(Self::new(points, faces).expect("A triangle fan is always a valid mesh"))
    }
}

// ===============================================================================================
// Shapes spanning two points
// ===============================================================================================

impl MeshData3 {
    /// Create a capsule mesh, which is a cylinder with a hemispherical cap on each end, spanning
    /// the segment between two points. The two points are the centers of the caps, so the capsule
    /// extends by `radius` beyond each endpoint of the segment.
    ///
    /// Every point lies on the true surface, so the facets sag inward. The caps are curved in two
    /// directions, and their contributions are additive, so each is held to half of `tol`. The tube
    /// inherits the circumferential density of the cap equators because they share the same rings.
    ///
    /// # Arguments
    ///
    /// * `p0`: one end of the segment
    /// * `p1`: the other end of the segment
    /// * `radius`: radius of the cylinder and of the caps, which must be positive and finite
    /// * `tol`: the maximum allowed deviation of a facet from the true surface, which must be
    ///   positive and finite. No matter how loose this is, a full turn around the tube never gets
    ///   fewer than 8 segments and each cap never gets fewer than 2 rows.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid, or if the two points
    /// are coincident and so do not define a direction
    pub fn create_capsule(p0: &Point3, p1: &Point3, radius: f64, tol: f64) -> Result<Self> {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;

        // As with the cylinder, the capsule is built along its own z-axis and then moved. The second
        // basis vector is arbitrary for a body of revolution.
        let transform = Iso3::from_basis_zx(&v, &Vector3::x(), Some(pc))
            .or_else(|_| Iso3::from_basis_zy(&v, &Vector3::y(), Some(pc)))
            .map_err(|_| "The two ends of a capsule must not be coincident")?;

        let n_theta = arc_segments_for_tol(radius, TAU, tol / 2.0)?;
        let n_cap = arc_segments_for_tol(radius, FRAC_PI_2, tol / 2.0)?;
        let hz = v.norm() / 2.0;

        let mut points = Vec::with_capacity(2 * n_cap * n_theta + 2);
        let mut faces = Vec::with_capacity(4 * n_cap * n_theta);

        let north = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, hz + radius));

        // The rings proceed from the north pole down over the top cap, across the tube, and then
        // continue down over the bottom cap. The last ring of the top cap and the first ring of the
        // bottom cap form the two ends of the tube, so banding every consecutive pair builds the
        // tube as well.
        let mut rings = Vec::with_capacity(2 * n_cap);
        for k in 1..=n_cap {
            let phi = FRAC_PI_2 * (k as f64) / (n_cap as f64);
            let (sin, cos) = phi.sin_cos();
            rings.push(push_ring(
                &mut points,
                radius * sin,
                hz + radius * cos,
                n_theta,
            ));
        }
        for k in (1..=n_cap).rev() {
            let phi = FRAC_PI_2 * (k as f64) / (n_cap as f64);
            let (sin, cos) = phi.sin_cos();
            rings.push(push_ring(
                &mut points,
                radius * sin,
                -hz - radius * cos,
                n_theta,
            ));
        }

        let south = points.len() as u32;
        points.push(Point3::new(0.0, 0.0, -hz - radius));

        fan_faces(&mut faces, north, rings[0], n_theta, true);
        for w in rings.windows(2) {
            band_faces(&mut faces, w[1], w[0], n_theta);
        }
        fan_faces(&mut faces, south, rings[rings.len() - 1], n_theta, false);

        let mut mesh =
            Self::new(points, faces).expect("A capsule tessellation is always a valid mesh");
        mesh.transform_in_place(&transform);
        Ok(mesh)
    }

    /// Create a cylindrical mesh spanning the segment between two points, with the given radius.
    ///
    /// # Arguments
    ///
    /// * `p0`: one end of the segment
    /// * `p1`: the other end of the segment
    /// * `radius`: radius of the cylinder, which must be positive and finite
    /// * `tol`: the maximum allowed chordal deviation of the wall, which must be positive and
    ///   finite. The circumference never gets fewer than 8 segments, no matter how loose this is.
    ///
    /// returns: `Result<MeshData3>`, failing if `radius` or `tol` is invalid, or if the two points
    /// are coincident and so do not define a direction
    pub fn create_cylinder_between(
        p0: &Point3,
        p1: &Point3,
        radius: f64,
        tol: f64,
    ) -> Result<Self> {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;

        // The cylinder is built along its own z-axis, so the segment direction becomes z. The
        // world axis that supplies the second basis vector does not matter for a body of
        // revolution, so x is tried first and y is the fallback when the segment is already
        // parallel to x. Both can fail only when the segment has no direction at all.
        let transform = Iso3::from_basis_zx(&v, &Vector3::x(), Some(pc))
            .or_else(|_| Iso3::from_basis_zy(&v, &Vector3::y(), Some(pc)))
            .map_err(|_| "The two ends of a cylinder must not be coincident")?;

        let mut mesh = Self::create_cylinder(radius, v.norm(), tol)?;
        mesh.transform_in_place(&transform);
        Ok(mesh)
    }

    /// Create a rectangular beam spanning the segment between two points, with its cross section
    /// oriented by an up direction.
    ///
    /// # Arguments
    ///
    /// * `p0`: one end of the segment
    /// * `p1`: the other end of the segment
    /// * `width`: the cross section dimension perpendicular to both the segment and `up`
    /// * `height`: the cross section dimension along `up`
    /// * `up`: the direction the beam's height is measured along, which must not be parallel to the
    ///   segment
    ///
    /// returns: `Result<MeshData3>`, failing if `up` is parallel to the segment, since the cross
    /// section would have no defined orientation
    pub fn create_rect_beam_between(
        p0: &Point3,
        p1: &Point3,
        width: f64,
        height: f64,
        up: &Vector3,
    ) -> Result<Self> {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;

        // The box is built with its long dimension along its own z-axis, so the segment direction
        // becomes z and `up` becomes y.
        let transform = Iso3::from_basis_zy(&v, up, Some(pc))?;

        let mut mesh = Self::create_box(width, height, v.norm());
        mesh.transform_in_place(&transform);
        Ok(mesh)
    }
}

// ===============================================================================================
// Sample meshes embedded in the binary
// ===============================================================================================

impl MeshData3 {
    /// Load a Stanford bunny mesh embedded in the binary with 453 points and 948 faces. This mesh
    /// has been quantized to 16 bits per axis. The mesh structure is the same as the corresponding
    /// `bun_zipper_res4.ply` mesh, but some precision has been lost in the conversion. The maximum
    /// point deviation from the original is 0.00000189 meters.
    ///
    /// Note that the vertex order is not the same as the original ply file's.
    pub fn stanford_bunny_res4() -> Self {
        embedded_mesh(include_bytes!(
            "../../../../tests/data/stanford_bun_4.tcmesh"
        ))
    }

    /// Load a Stanford bunny mesh embedded in the binary with 1889 points and 3851 faces. This mesh
    /// has been quantized to 16 bits per axis. The mesh structure is the same as the corresponding
    /// `bun_zipper_res3.ply` mesh, but some precision has been lost in the conversion. The maximum
    /// point deviation from the original is 0.00000189 meters.
    ///
    /// Note that the vertex order is not the same as the original ply file's.
    pub fn stanford_bunny_res3() -> Self {
        embedded_mesh(include_bytes!(
            "../../../../tests/data/stanford_bun_3.tcmesh"
        ))
    }

    /// Load a Stanford bunny mesh embedded in the binary with 8171 points and 16301 faces. This
    /// mesh has been quantized to 16 bits per axis. The mesh structure is the same as the
    /// corresponding `bun_zipper_res2.ply` mesh, but some precision has been lost in the conversion.
    /// The maximum point deviation from the original is 0.00000189 meters.
    ///
    /// Note that the vertex order is not the same as the original ply file's.
    pub fn stanford_bunny_res2() -> Self {
        embedded_mesh(include_bytes!(
            "../../../../tests/data/stanford_bun_2.tcmesh"
        ))
    }
}

/// Decode one of the embedded fixtures.
///
/// These were stored in a 16-bit micro-mesh format until that was retired in favour of tcmesh. The
/// re-encoding was deliberately done at a tolerance that lands on the same 16 bits per axis over the
/// same bounds, so the coordinates are bit-identical to what the old format held and the deviation
/// figures above did not have to be recomputed. Only the vertex numbering changed, because writing
/// a tcmesh reorders vertices so the connectivity compresses.
fn embedded_mesh(bytes: &[u8]) -> MeshData3 {
    read_tc_mesh_from(&mut { bytes }).expect("embedded fixture must decode")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Mesh3;
    use crate::na;
    use approx::assert_relative_eq;

    /// Signed volume by the divergence theorem: positive when the winding is counterclockwise as
    /// viewed from outside everywhere, and equal to the enclosed volume when the mesh is closed.
    fn signed_volume(data: &MeshData3) -> f64 {
        data.faces()
            .iter()
            .map(|f| {
                let p0 = data.points()[f[0] as usize].coords;
                let p1 = data.points()[f[1] as usize].coords;
                let p2 = data.points()[f[2] as usize].coords;
                p0.dot(&p1.cross(&p2)) / 6.0
            })
            .sum()
    }

    /// A closed primitive must produce a manifold mesh with no boundary loops.
    fn assert_watertight(data: &MeshData3) {
        let mesh = Mesh3::from_data(data.clone(), true).expect("primitive must build");
        let edges = mesh.compute_edges().expect("primitive must be manifold");
        assert!(edges.boundary_loops.is_empty());
    }

    /// Every perimeter chord's midpoint must sag inward from the true circle by no more than
    /// the requested tolerance.
    #[test]
    fn a_circle_respects_the_chordal_tolerance() {
        for &tol in &[0.02, 1.0e-3] {
            let radius = 2.5;
            let circle = MeshData3::create_circle(radius, tol).unwrap();
            for f in circle.faces() {
                let a = circle.points()[f[1] as usize];
                let b = circle.points()[f[2] as usize];
                let mid = na::center(&a, &b);
                let sag = radius - mid.coords.norm();
                assert!(sag >= 0.0, "chord midpoint outside the circle");
                assert!(sag <= tol + 1.0e-12, "sag {sag} exceeds tol {tol}");
            }
        }
    }

    #[test]
    fn a_circle_never_gets_fewer_than_eight_segments() {
        let circle = MeshData3::create_circle(1.0, 10.0).unwrap();
        assert_eq!(circle.face_count(), 8);
        assert_eq!(circle.point_count(), 9);
    }

    #[test]
    fn a_circle_is_flat_open_and_faces_up() {
        let circle = MeshData3::create_circle(1.0, 1.0e-3).unwrap();

        // One boundary loop: the perimeter
        let mesh = Mesh3::from_data(circle.clone(), false).unwrap();
        let edges = mesh.compute_edges().unwrap();
        assert_eq!(edges.boundary_loops.len(), 1);

        // Every face normal points along +Z
        for n in circle.compute_face_normals().unwrap() {
            assert_relative_eq!(n.dot(&Vector3::z()), 1.0, epsilon = 1.0e-12);
        }
    }

    #[test]
    fn a_circle_rejects_invalid_arguments() {
        assert!(MeshData3::create_circle(1.0, 0.0).is_err());
        assert!(MeshData3::create_circle(1.0, f64::NAN).is_err());
        assert!(MeshData3::create_circle(0.0, 0.1).is_err());
    }

    #[test]
    fn a_box_has_positive_and_exact_volume() {
        let bx = MeshData3::create_box(2.0, 3.0, 4.0);
        assert_relative_eq!(signed_volume(&bx), 24.0, epsilon = 1.0e-12);
    }

    #[test]
    fn a_box_is_watertight() {
        assert_watertight(&MeshData3::create_box(2.0, 3.0, 4.0));
    }

    #[test]
    fn a_beam_between_is_watertight_and_spans_the_points() {
        let p0 = Point3::new(1.0, 2.0, 3.0);
        let p1 = Point3::new(4.0, 6.0, 3.0);
        let beam = MeshData3::create_rect_beam_between(&p0, &p1, 1.0, 2.0, &Vector3::z()).unwrap();

        assert_watertight(&beam);
        assert_relative_eq!(signed_volume(&beam), 1.0 * 2.0 * 5.0, epsilon = 1.0e-12);

        // The vertex centroid sits at the segment midpoint
        let center = beam
            .points()
            .iter()
            .fold(Vector3::zeros(), |acc, p| acc + p.coords)
            / beam.point_count() as f64;
        assert_relative_eq!(center.x, 2.5, epsilon = 1.0e-12);
        assert_relative_eq!(center.y, 4.0, epsilon = 1.0e-12);
        assert_relative_eq!(center.z, 3.0, epsilon = 1.0e-12);
    }

    /// A cone's apex must sit at +height/2 on the z-axis and its base ring at -height/2, at the
    /// requested radius. The two dimension arguments used to be swapped on the way through to
    /// `parry3d`, so this pins which is which.
    #[test]
    fn a_cone_has_the_requested_radius_and_full_height() {
        let cone = MeshData3::create_cone(2.0, 10.0, 1.0e-3).unwrap();

        let max_z = cone.points().iter().map(|p| p.z).fold(f64::MIN, f64::max);
        let min_z = cone.points().iter().map(|p| p.z).fold(f64::MAX, f64::min);
        assert_relative_eq!(max_z, 5.0, epsilon = 1.0e-12);
        assert_relative_eq!(min_z, -5.0, epsilon = 1.0e-12);

        let max_r = cone
            .points()
            .iter()
            .map(|p| p.x.hypot(p.y))
            .fold(f64::MIN, f64::max);
        assert_relative_eq!(max_r, 2.0, epsilon = 1.0e-12);

        // The apex is a single point on the axis, not a ring.
        let on_axis = cone
            .points()
            .iter()
            .filter(|p| p.x.hypot(p.y) < 1.0e-12 && p.z > 0.0)
            .count();
        assert_eq!(on_axis, 1);
    }

    /// Every base chord's midpoint must sag inward from the true base circle by no more than the
    /// requested tolerance.
    #[test]
    fn a_cone_respects_the_chordal_tolerance() {
        for &tol in &[0.02, 1.0e-3] {
            let radius = 2.5;
            let cone = MeshData3::create_cone(radius, 4.0, tol).unwrap();

            // The base ring comes first, followed by the apex and the base center.
            let n = cone.point_count() - 2;
            for i in 0..n {
                let a = cone.points()[i];
                let b = cone.points()[(i + 1) % n];
                let mid = na::center(&a, &b);
                let sag = radius - mid.x.hypot(mid.y);
                assert!(sag >= 0.0, "chord midpoint outside the base circle");
                assert!(sag <= tol + 1.0e-12, "sag {sag} exceeds tol {tol}");
            }
        }
    }

    #[test]
    fn a_cone_never_gets_fewer_than_eight_segments() {
        let cone = MeshData3::create_cone(1.0, 2.0, 10.0).unwrap();
        assert_eq!(cone.point_count(), 8 + 2);
        assert_eq!(cone.face_count(), 8 * 2);
    }

    #[test]
    fn a_cone_is_watertight_and_wound_outward() {
        let cone = MeshData3::create_cone(2.0, 10.0, 1.0e-4).unwrap();
        assert_watertight(&cone);

        // Inscribed, so the tessellation slightly under-fills the true volume.
        let volume = signed_volume(&cone);
        let exact = std::f64::consts::PI * 4.0 * 10.0 / 3.0;
        assert!(volume > 0.0, "winding is inside out");
        assert!(volume <= exact, "volume {volume} exceeds the true {exact}");
        assert_relative_eq!(volume, exact, max_relative = 1.0e-3);
    }

    #[test]
    fn a_cone_rejects_invalid_arguments() {
        assert!(MeshData3::create_cone(1.0, 2.0, 0.0).is_err());
        assert!(MeshData3::create_cone(0.0, 2.0, 0.1).is_err());
    }

    /// The cylinder takes a full height too, so the two round bodies agree on what the argument
    /// means, and its axis runs along z.
    #[test]
    fn a_cylinder_has_the_requested_radius_and_full_height() {
        let cyl = MeshData3::create_cylinder(2.0, 10.0, 1.0e-3).unwrap();

        for p in cyl.points() {
            assert_relative_eq!(p.z.abs(), 5.0, epsilon = 1.0e-12);

            // Wall points are out at the radius, the two cap centers are on the axis.
            let r = p.x.hypot(p.y);
            assert!(
                r == 0.0 || (r - 2.0).abs() < 1.0e-12,
                "unexpected radius {r}"
            );
        }
    }

    /// Every wall chord's midpoint must sag inward from the true cylinder by no more than the
    /// requested tolerance.
    #[test]
    fn a_cylinder_respects_the_chordal_tolerance() {
        for &tol in &[0.02, 1.0e-3] {
            let radius = 2.5;
            let cyl = MeshData3::create_cylinder(radius, 4.0, tol).unwrap();

            // Adjacent points on a wall ring define the chords, so walk the bottom ring.
            let n = (cyl.point_count() - 2) / 2;
            for i in 0..n {
                let a = cyl.points()[i];
                let b = cyl.points()[(i + 1) % n];
                let mid = na::center(&a, &b);
                let sag = radius - mid.x.hypot(mid.y);
                assert!(sag >= 0.0, "chord midpoint outside the cylinder");
                assert!(sag <= tol + 1.0e-12, "sag {sag} exceeds tol {tol}");
            }
        }
    }

    #[test]
    fn a_cylinder_never_gets_fewer_than_eight_segments() {
        let cyl = MeshData3::create_cylinder(1.0, 2.0, 10.0).unwrap();
        assert_eq!(cyl.point_count(), 8 * 2 + 2);
        assert_eq!(cyl.face_count(), 8 * 4);
    }

    #[test]
    fn a_cylinder_is_watertight_and_wound_outward() {
        let cyl = MeshData3::create_cylinder(2.0, 10.0, 1.0e-4).unwrap();
        assert_watertight(&cyl);

        // Inscribed, so the tessellation slightly under-fills the true volume.
        let volume = signed_volume(&cyl);
        let exact = std::f64::consts::PI * 4.0 * 10.0;
        assert!(volume > 0.0, "winding is inside out");
        assert!(volume <= exact, "volume {volume} exceeds the true {exact}");
        assert_relative_eq!(volume, exact, max_relative = 1.0e-3);
    }

    #[test]
    fn a_cylinder_rejects_invalid_arguments() {
        assert!(MeshData3::create_cylinder(1.0, 2.0, 0.0).is_err());
        assert!(MeshData3::create_cylinder(0.0, 2.0, 0.1).is_err());
    }

    #[test]
    fn a_cylinder_between_spans_the_two_points() {
        let p0 = Point3::new(1.0, 2.0, 3.0);
        let p1 = Point3::new(1.0, 2.0, 9.0);
        let cyl = MeshData3::create_cylinder_between(&p0, &p1, 0.5, 1.0e-3).unwrap();

        let min_z = cyl.points().iter().map(|p| p.z).fold(f64::MAX, f64::min);
        let max_z = cyl.points().iter().map(|p| p.z).fold(f64::MIN, f64::max);
        assert_relative_eq!(min_z, 3.0, epsilon = 1.0e-12);
        assert_relative_eq!(max_z, 9.0, epsilon = 1.0e-12);

        // Every point is either on the surface of the tube around the segment or at the center of
        // one of the two end caps.
        for p in cyl.points() {
            let r = (p.x - 1.0).hypot(p.y - 2.0);
            assert!(
                r < 1.0e-12 || (r - 0.5).abs() < 1.0e-12,
                "off the tube: {r}"
            );
        }

        assert_watertight(&cyl);
    }

    #[test]
    fn a_cylinder_between_needs_two_distinct_points() {
        let p0 = Point3::new(1.0, 2.0, 3.0);
        assert!(MeshData3::create_cylinder_between(&p0, &p0, 0.5, 1.0e-3).is_err());
    }

    #[test]
    fn a_beam_between_needs_an_up_which_is_not_parallel_to_the_segment() {
        let p0 = Point3::origin();
        let p1 = Point3::new(0.0, 0.0, 10.0);

        assert!(MeshData3::create_rect_beam_between(&p0, &p1, 1.0, 2.0, &Vector3::y()).is_ok());
        assert!(MeshData3::create_rect_beam_between(&p0, &p1, 1.0, 2.0, &Vector3::z()).is_err());
    }

    /// The primitives carry geometry and nothing else, so none of them may invent an attribute.
    #[test]
    fn no_primitive_arrives_with_attributes() {
        assert!(MeshData3::create_box(1.0, 1.0, 1.0).attrs().is_empty());
        assert!(
            MeshData3::create_sphere(1.0, 1.0e-3)
                .unwrap()
                .attrs()
                .is_empty()
        );
        assert!(
            MeshData3::create_cylinder(1.0, 2.0, 1.0e-3)
                .unwrap()
                .attrs()
                .is_empty()
        );
        assert!(
            MeshData3::create_cone(1.0, 2.0, 1.0e-3)
                .unwrap()
                .attrs()
                .is_empty()
        );
        assert!(
            MeshData3::create_circle(1.0, 1.0e-3)
                .unwrap()
                .attrs()
                .is_empty()
        );
        assert!(MeshData3::stanford_bunny_res4().attrs().is_empty());
    }

    /// No point of any facet may sit outside the true sphere, and none may sag further inside it
    /// than the requested tolerance.
    #[test]
    fn a_sphere_respects_the_tolerance() {
        for &tol in &[0.05, 1.0e-3] {
            let radius = 2.5;
            let sphere = MeshData3::create_sphere(radius, tol).unwrap();

            let mut worst: f64 = 0.0;
            for f in sphere.faces() {
                let a = sphere.points()[f[0] as usize];
                let b = sphere.points()[f[1] as usize];
                let c = sphere.points()[f[2] as usize];

                let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.0);
                for p in [
                    centroid,
                    na::center(&a, &b),
                    na::center(&b, &c),
                    na::center(&c, &a),
                ] {
                    let sag = radius - p.coords.norm();
                    assert!(sag >= -1.0e-12, "a facet point sits outside the sphere");
                    worst = worst.max(sag);
                }
            }

            assert!(
                worst <= tol + 1.0e-12,
                "worst sag {worst} exceeds tol {tol}"
            );
        }
    }

    #[test]
    fn a_sphere_never_gets_fewer_than_the_floor_counts() {
        let sphere = MeshData3::create_sphere(1.0, 10.0).unwrap();

        // 8 segments around the equator and 4 from pole to pole produce 3 rings plus the two poles.
        assert_eq!(sphere.point_count(), 8 * 3 + 2);
        assert_eq!(sphere.face_count(), 8 * 2 + 2 * 8 * 2);
    }

    #[test]
    fn a_sphere_has_its_poles_on_z() {
        let radius = 3.0;
        let sphere = MeshData3::create_sphere(radius, 1.0e-3).unwrap();

        for p in sphere.points() {
            assert_relative_eq!(p.coords.norm(), radius, epsilon = 1.0e-12);
        }

        // One point at each pole, and nothing else on the axis.
        let on_axis = sphere
            .points()
            .iter()
            .filter(|p| p.x.hypot(p.y) < 1.0e-12)
            .count();
        assert_eq!(on_axis, 2);

        let max_z = sphere.points().iter().map(|p| p.z).fold(f64::MIN, f64::max);
        let min_z = sphere.points().iter().map(|p| p.z).fold(f64::MAX, f64::min);
        assert_relative_eq!(max_z, radius, epsilon = 1.0e-12);
        assert_relative_eq!(min_z, -radius, epsilon = 1.0e-12);
    }

    #[test]
    fn a_sphere_is_watertight_and_wound_outward() {
        let radius = 2.0;
        let sphere = MeshData3::create_sphere(radius, 1.0e-4).unwrap();
        assert_watertight(&sphere);

        // Inscribed, so the tessellation slightly under-fills the true volume.
        let volume = signed_volume(&sphere);
        let exact = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
        assert!(volume > 0.0, "winding is inside out");
        assert!(volume <= exact, "volume {volume} exceeds the true {exact}");
        assert_relative_eq!(volume, exact, max_relative = 1.0e-3);
    }

    #[test]
    fn a_sphere_rejects_invalid_arguments() {
        assert!(MeshData3::create_sphere(1.0, 0.0).is_err());
        assert!(MeshData3::create_sphere(1.0, f64::NAN).is_err());
        assert!(MeshData3::create_sphere(0.0, 0.1).is_err());
    }

    /// No point of any facet may sit outside the true capsule surface, and none may sag further
    /// inside it than the requested tolerance.
    #[test]
    fn a_capsule_respects_the_tolerance() {
        for &tol in &[0.05, 1.0e-3] {
            let radius = 2.5;
            let hz = 4.0;
            let p0 = Point3::new(0.0, 0.0, -hz);
            let p1 = Point3::new(0.0, 0.0, hz);
            let capsule = MeshData3::create_capsule(&p0, &p1, radius, tol).unwrap();

            // Distance from a point to the capsule's axis segment, which is `radius` everywhere on
            // the true surface.
            let surface_distance = |p: &Point3| {
                let z = p.z.clamp(-hz, hz);
                (p - Point3::new(0.0, 0.0, z)).norm()
            };

            let mut worst: f64 = 0.0;
            for f in capsule.faces() {
                let a = capsule.points()[f[0] as usize];
                let b = capsule.points()[f[1] as usize];
                let c = capsule.points()[f[2] as usize];

                let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.0);
                for p in [
                    centroid,
                    na::center(&a, &b),
                    na::center(&b, &c),
                    na::center(&c, &a),
                ] {
                    let sag = radius - surface_distance(&p);
                    assert!(sag >= -1.0e-12, "a facet point sits outside the capsule");
                    worst = worst.max(sag);
                }
            }

            assert!(
                worst <= tol + 1.0e-12,
                "worst sag {worst} exceeds tol {tol}"
            );
        }
    }

    #[test]
    fn a_capsule_never_gets_fewer_than_the_floor_counts() {
        let p0 = Point3::origin();
        let p1 = Point3::new(0.0, 0.0, 4.0);
        let capsule = MeshData3::create_capsule(&p0, &p1, 1.0, 10.0).unwrap();

        // 8 segments around, 2 rows per cap, so 4 rings plus the two poles.
        assert_eq!(capsule.point_count(), 8 * 4 + 2);
        assert_eq!(capsule.face_count(), 8 * 2 + 3 * 8 * 2);
    }

    /// The two points are the cap centers, so the capsule runs half a radius past each of them,
    /// and its poles land on the segment direction.
    #[test]
    fn a_capsule_spans_the_two_points_plus_a_radius_at_each_end() {
        let p0 = Point3::new(1.0, 2.0, 3.0);
        let p1 = Point3::new(1.0, 2.0, 9.0);
        let radius = 0.5;
        let capsule = MeshData3::create_capsule(&p0, &p1, radius, 1.0e-3).unwrap();

        let min_z = capsule
            .points()
            .iter()
            .map(|p| p.z)
            .fold(f64::MAX, f64::min);
        let max_z = capsule
            .points()
            .iter()
            .map(|p| p.z)
            .fold(f64::MIN, f64::max);
        assert_relative_eq!(min_z, 3.0 - radius, epsilon = 1.0e-12);
        assert_relative_eq!(max_z, 9.0 + radius, epsilon = 1.0e-12);

        assert_watertight(&capsule);
    }

    /// A capsule is a sphere with a cylinder spliced into its equator, so its volume is the sum of
    /// the two.
    #[test]
    fn a_capsule_is_watertight_and_wound_outward() {
        let radius = 2.0;
        let length = 10.0;
        let p0 = Point3::new(0.0, 0.0, -length / 2.0);
        let p1 = Point3::new(0.0, 0.0, length / 2.0);
        let capsule = MeshData3::create_capsule(&p0, &p1, radius, 1.0e-4).unwrap();
        assert_watertight(&capsule);

        let volume = signed_volume(&capsule);
        let exact = std::f64::consts::PI * radius.powi(2) * (length + 4.0 / 3.0 * radius);
        assert!(volume > 0.0, "winding is inside out");
        assert!(volume <= exact, "volume {volume} exceeds the true {exact}");
        assert_relative_eq!(volume, exact, max_relative = 1.0e-3);
    }

    #[test]
    fn a_capsule_rejects_invalid_arguments() {
        let p0 = Point3::origin();
        let p1 = Point3::new(0.0, 0.0, 4.0);

        assert!(MeshData3::create_capsule(&p0, &p1, 1.0, 0.0).is_err());
        assert!(MeshData3::create_capsule(&p0, &p1, 0.0, 0.1).is_err());
        assert!(MeshData3::create_capsule(&p0, &p0, 1.0, 0.1).is_err());
    }

    #[test]
    fn the_embedded_bunnies_have_their_documented_sizes() {
        assert_eq!(MeshData3::stanford_bunny_res4().point_count(), 453);
        assert_eq!(MeshData3::stanford_bunny_res4().face_count(), 948);
        assert_eq!(MeshData3::stanford_bunny_res3().point_count(), 1889);
        assert_eq!(MeshData3::stanford_bunny_res2().point_count(), 8171);
    }
}
