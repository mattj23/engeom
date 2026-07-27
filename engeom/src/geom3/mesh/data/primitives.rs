//! This module contains the constructors which build a `MeshData3` for a simple geometric shape,
//! along with the sample meshes embedded in the binary.
//!
//! These live on the data container rather than the accelerated one because a tessellation produces
//! points and faces and nothing else. `Mesh3` carries a pass-through for each, adding the `is_solid`
//! flag, which is a query behavior and not part of the geometry.
//!
//! # These are not our tessellators yet
//!
//! Most of the bodies below still call `parry3d`'s `to_trimesh`, which is why the tessellation
//! density arguments read the way they do and why the axis conventions vary between shapes. The
//! public signatures here are the ones we intend to keep, so each body can be replaced with our own
//! implementation one at a time without another breaking change. What we gain by writing them
//! ourselves is control over winding, seam handling, and whether point normals come out of the
//! construction rather than having to be computed afterwards.

use super::MeshData3;
use crate::geom3::IsoExtensions3;
use crate::io::{deflate_bytes, u_bytes_to_mesh_data};
use crate::{Iso3, Point3, Result, Vector3};
use parry3d_f64::shape;

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
        let bx = shape::Cuboid::new(Vector3::new(length / 2.0, width / 2.0, height / 2.0));
        let (points, faces) = bx.to_trimesh();
        Self::new(points, faces).expect("A cuboid tessellation is always a valid mesh")
    }

    /// Create a spherical mesh centered at the origin. The `n_theta` and `n_phi` parameters control
    /// the tessellation density.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the sphere
    /// * `n_theta`: number of subdivisions around the polar direction
    /// * `n_phi`: number of subdivisions around the azimuthal direction
    ///
    /// returns: `MeshData3`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let n_t = 14;
    /// let n_p = 15;
    /// let sphere = MeshData3::create_sphere(1.0, n_t, n_p);
    ///
    /// assert_eq!(sphere.point_count(), n_t * (n_p - 1) + 2);
    ///
    /// // Verify that the points are on the surface of the unit sphere.
    /// for p in sphere.points() {
    ///     assert_relative_eq!(p.coords.norm(), 1.0)
    /// }
    /// ```
    pub fn create_sphere(radius: f64, n_theta: usize, n_phi: usize) -> Self {
        let sphere = shape::Ball::new(radius);
        let (points, faces) = sphere.to_trimesh(n_theta as u32, n_phi as u32);
        Self::new(points, faces).expect("A ball tessellation is always a valid mesh")
    }

    /// Create a cylindrical mesh centered at the origin and aligned with the local `y` axis.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the cylinder
    /// * `height`: full height of the cylinder, along the y axis
    /// * `steps`: number of subdivisions around the circumference
    ///
    /// returns: `MeshData3`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{MeshData3, Point3};
    /// use approx::assert_relative_eq;
    ///
    /// let cyl = MeshData3::create_cylinder(1.0, 4.0, 16);
    ///
    /// for p in cyl.points() {
    ///     assert_relative_eq!(p.y.abs(), 2.0);
    ///
    ///     let proj = Point3::new(p.x, 0.0, p.z);
    ///     assert_relative_eq!(proj.coords.norm(), 1.0);
    /// }
    /// ```
    pub fn create_cylinder(radius: f64, height: f64, steps: usize) -> Self {
        let cyl = shape::Cylinder::new(height / 2.0, radius);
        let (points, faces) = cyl.to_trimesh(steps as u32);
        Self::new(points, faces).expect("A cylinder tessellation is always a valid mesh")
    }

    /// Create a conical mesh centered at the origin and aligned with the local `y` axis, with its
    /// apex at `+height/2` and its base at `-height/2`.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the base of the cone
    /// * `height`: full height of the cone, along the y axis
    /// * `steps`: number of subdivisions around the circumference
    ///
    /// returns: `MeshData3`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    /// use approx::assert_relative_eq;
    ///
    /// let cone = MeshData3::create_cone(2.0, 10.0, 32);
    ///
    /// // The apex sits at +height/2 and the base ring at -height/2, at the given radius.
    /// let apex = cone.points().iter().map(|p| p.y).fold(f64::MIN, f64::max);
    /// assert_relative_eq!(apex, 5.0);
    ///
    /// let widest = cone.points().iter().map(|p| p.x.hypot(p.z)).fold(f64::MIN, f64::max);
    /// assert_relative_eq!(widest, 2.0);
    /// ```
    pub fn create_cone(radius: f64, height: f64, steps: usize) -> Self {
        let cone = shape::Cone::new(height / 2.0, radius);
        let (points, faces) = cone.to_trimesh(steps as u32);
        Self::new(points, faces).expect("A cone tessellation is always a valid mesh")
    }

    /// Create a flat, filled circle mesh lying in the XY plane, centered at the origin, with the
    /// normal pointing along +Z. The mesh is a triangle fan from the center to `segments` evenly
    /// spaced perimeter points.
    ///
    /// # Arguments
    ///
    /// * `radius`: radius of the circle
    /// * `segments`: number of perimeter points, and of triangles. Must be at least 3.
    ///
    /// returns: `MeshData3`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::MeshData3;
    ///
    /// let circle = MeshData3::create_circle(1.0, 32);
    /// assert_eq!(circle.point_count(), 33); // center + 32 perimeter
    /// assert_eq!(circle.face_count(), 32);
    /// ```
    pub fn create_circle(radius: f64, segments: usize) -> Self {
        use std::f64::consts::TAU;
        let mut points = Vec::with_capacity(segments + 1);
        let mut faces = Vec::with_capacity(segments);

        points.push(Point3::origin());
        for i in 0..segments {
            let angle = TAU * (i as f64) / (segments as f64);
            points.push(Point3::new(radius * angle.cos(), radius * angle.sin(), 0.0));
        }

        for i in 0..segments {
            let a = (i + 1) as u32;
            let b = ((i + 1) % segments + 1) as u32;
            faces.push([0u32, a, b]);
        }

        Self::new(points, faces).expect("A triangle fan is always a valid mesh")
    }
}

// ===============================================================================================
// Shapes spanning two points
// ===============================================================================================

impl MeshData3 {
    /// Create a capsule mesh, which is a cylinder with a hemispherical cap on each end, spanning
    /// the segment between two points.
    ///
    /// # Arguments
    ///
    /// * `p0`: one end of the segment
    /// * `p1`: the other end of the segment
    /// * `radius`: radius of the cylinder and of the caps
    /// * `n_theta`: number of subdivisions around the circumference
    /// * `n_phi`: number of subdivisions over each cap
    ///
    /// returns: `MeshData3`
    pub fn create_capsule(
        p0: &Point3,
        p1: &Point3,
        radius: f64,
        n_theta: usize,
        n_phi: usize,
    ) -> Self {
        let capsule = shape::Capsule::new(*p0, *p1, radius);
        let (points, faces) = capsule.to_trimesh(n_theta as u32, n_phi as u32);
        Self::new(points, faces).expect("A capsule tessellation is always a valid mesh")
    }

    /// Create a cylindrical mesh spanning the segment between two points, with the given radius.
    ///
    /// # Arguments
    ///
    /// * `p0`: one end of the segment
    /// * `p1`: the other end of the segment
    /// * `radius`: radius of the cylinder
    /// * `steps`: number of subdivisions around the circumference
    ///
    /// returns: `MeshData3`
    pub fn create_cylinder_between(p0: &Point3, p1: &Point3, radius: f64, steps: usize) -> Self {
        let v = *p1 - *p0;
        let pc = *p0 + v / 2.0;
        let cyl = shape::Cylinder::new(v.norm() / 2.0, radius);

        // The cylinder is built along its own y axis, so the segment direction becomes y. Which
        // world axis supplies the second basis vector does not matter for a body of revolution, so
        // z is tried first and x is the fallback for a segment which is already parallel to z.
        let transform = Iso3::from_basis_yz(&v, &Vector3::z(), Some(pc))
            .unwrap_or(Iso3::from_basis_yx(&v, &Vector3::x(), Some(pc)).unwrap());

        let (points, faces) = cyl.to_trimesh(steps as u32);
        let mut mesh =
            Self::new(points, faces).expect("A cylinder tessellation is always a valid mesh");
        mesh.transform_in_place(&transform);
        mesh
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
        let box_geom = shape::Cuboid::new(Vector3::new(width / 2.0, height / 2.0, v.norm() / 2.0));

        // The box is built with its long dimension along its own z axis, so the segment direction
        // becomes z and `up` becomes y.
        let transform = Iso3::from_basis_zy(&v, up, Some(pc))?;

        let (points, faces) = box_geom.to_trimesh();
        let mut mesh = Self::new(points, faces).expect("A cuboid tessellation is always valid");
        mesh.transform_in_place(&transform);
        Ok(mesh)
    }
}

// ===============================================================================================
// Sample meshes embedded in the binary
// ===============================================================================================

impl MeshData3 {
    /// Load a Stanford bunny mesh embedded in the binary with 453 points and 948 faces. This mesh
    /// has been compressed into the 16-bit micro mesh format. The mesh structure is the same as the
    /// corresponding `bun_zipper_res4.ply` mesh, but some precision has been lost in the conversion.
    /// The maximum point deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res4() -> Self {
        let bytes = include_bytes!("../../../../tests/data/stanford_bun_4.umesh.gz");
        u_bytes_to_mesh_data(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 1889 points and 3851 faces. This mesh
    /// has been compressed into the 16-bit micro mesh format. The mesh structure is the same as the
    /// corresponding `bun_zipper_res3.ply` mesh, but some precision has been lost in the conversion.
    /// The maximum point deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res3() -> Self {
        let bytes = include_bytes!("../../../../tests/data/stanford_bun_3.umesh.gz");
        u_bytes_to_mesh_data(&deflate_bytes(bytes).unwrap()).unwrap()
    }

    /// Load a Stanford bunny mesh embedded in the binary with 8171 points and 16301 faces. This
    /// mesh has been compressed into the 16-bit micro mesh format. The mesh structure is the same as
    /// the corresponding `bun_zipper_res2.ply` mesh, but some precision has been lost in the
    /// conversion. The maximum point deviation from the original is 0.00000189 meters.
    pub fn stanford_bunny_res2() -> Self {
        let bytes = include_bytes!("../../../../tests/data/stanford_bun_2.umesh.gz");
        u_bytes_to_mesh_data(&deflate_bytes(bytes).unwrap()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// A cone's apex must sit at +height/2 and its base ring at -height/2, with the base at the
    /// requested radius. The two arguments used to be swapped on the way through to `parry3d`, so
    /// this pins which is which.
    #[test]
    fn a_cone_has_the_requested_radius_and_full_height() {
        let cone = MeshData3::create_cone(2.0, 10.0, 32);

        let max_y = cone.points().iter().map(|p| p.y).fold(f64::MIN, f64::max);
        let min_y = cone.points().iter().map(|p| p.y).fold(f64::MAX, f64::min);
        assert_relative_eq!(max_y, 5.0, epsilon = 1.0e-12);
        assert_relative_eq!(min_y, -5.0, epsilon = 1.0e-12);

        let max_r = cone
            .points()
            .iter()
            .map(|p| p.x.hypot(p.z))
            .fold(f64::MIN, f64::max);
        assert_relative_eq!(max_r, 2.0, epsilon = 1.0e-12);
    }

    /// The cylinder takes a full height too, so the two round bodies agree on what the argument
    /// means.
    #[test]
    fn a_cylinder_has_the_requested_radius_and_full_height() {
        let cyl = MeshData3::create_cylinder(2.0, 10.0, 32);

        let max_y = cyl.points().iter().map(|p| p.y).fold(f64::MIN, f64::max);
        assert_relative_eq!(max_y, 5.0, epsilon = 1.0e-12);

        for p in cyl.points() {
            assert_relative_eq!(p.x.hypot(p.z), 2.0, epsilon = 1.0e-12);
        }
    }

    #[test]
    fn a_cylinder_between_spans_the_two_points() {
        let p0 = Point3::new(1.0, 2.0, 3.0);
        let p1 = Point3::new(1.0, 2.0, 9.0);
        let cyl = MeshData3::create_cylinder_between(&p0, &p1, 0.5, 16);

        let min_z = cyl.points().iter().map(|p| p.z).fold(f64::MAX, f64::min);
        let max_z = cyl.points().iter().map(|p| p.z).fold(f64::MIN, f64::max);
        assert_relative_eq!(min_z, 3.0, epsilon = 1.0e-12);
        assert_relative_eq!(max_z, 9.0, epsilon = 1.0e-12);

        // Every point is on the surface of the tube around the segment.
        for p in cyl.points() {
            assert_relative_eq!((p.x - 1.0).hypot(p.y - 2.0), 0.5, epsilon = 1.0e-12);
        }
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
        assert!(MeshData3::create_sphere(1.0, 8, 8).attrs().is_empty());
        assert!(MeshData3::create_cylinder(1.0, 2.0, 8).attrs().is_empty());
        assert!(MeshData3::create_cone(1.0, 2.0, 8).attrs().is_empty());
        assert!(MeshData3::create_circle(1.0, 8).attrs().is_empty());
        assert!(MeshData3::stanford_bunny_res4().attrs().is_empty());
    }

    #[test]
    fn the_embedded_bunnies_have_their_documented_sizes() {
        assert_eq!(MeshData3::stanford_bunny_res4().point_count(), 453);
        assert_eq!(MeshData3::stanford_bunny_res4().face_count(), 948);
        assert_eq!(MeshData3::stanford_bunny_res3().point_count(), 1889);
        assert_eq!(MeshData3::stanford_bunny_res2().point_count(), 8171);
    }
}
