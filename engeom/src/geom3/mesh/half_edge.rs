//! This module provides an implementation of a half edge mesh structure using the `alum` library.
//! The main type is `HalfEdgeMesh`, which is a polyhedral mesh that can be converted to a `Mesh`
//! type using the `TryFrom` trait. The `NaAdaptor` struct provides the necessary adaptors for
//! `alum` to work with `nalgebra` types.
//!
//! Submodules of this module wrap and compose `alum` functionality to provide tools for editing
//! meshes and their topology.

mod smoothing;

use crate::{Mesh, Point3, Vector3};
use alum;
use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use smoothing::*;
use std::error::Error;

pub trait HalfEdgeCloneOps {
    fn clone_vertices(&self) -> crate::Result<Vec<Point3>>;

    fn clone_faces(&self) -> crate::Result<Vec<[u32; 3]>>;
}

impl HalfEdgeCloneOps for HalfEdgeMesh {
    fn clone_vertices(&self) -> crate::Result<Vec<Point3>> {
        let borrow_vert = self.points();
        let borrow_vert = borrow_vert
            .try_borrow()
            .map_err(|_| "Failed to borrow points")?;
        Ok(borrow_vert.iter().map(|v| Point3::from(*v)).collect())
    }

    fn clone_faces(&self) -> crate::Result<Vec<[u32; 3]>> {
        let f_status = self.face_status_prop();
        let f_status = f_status
            .try_borrow()
            .map_err(|_| "Failed to borrow face status")?;
        Ok(self
            .triangulated_vertices(&f_status)
            .map(|f| [f[0].index(), f[1].index(), f[2].index()])
            .collect())
    }
}

impl TryFrom<&HalfEdgeMesh> for Mesh {
    type Error = Box<dyn Error>;

    fn try_from(value: &HalfEdgeMesh) -> Result<Self, Self::Error> {
        let vertices = value.clone_vertices()?;
        let faces = value.clone_faces()?;
        Ok(Mesh::new(vertices, faces, false))
    }
}

impl TryFrom<&Mesh> for HalfEdgeMesh {
    type Error = Box<dyn Error>;

    fn try_from(value: &Mesh) -> Result<Self, Self::Error> {
        let mut result = HalfEdgeMesh::new();
        let mut indices = Vec::new();
        for v in value.vertices() {
            let handle = result
                .add_vertex(v.coords)
                .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
            indices.push(handle);
        }

        for f in value.faces() {
            result
                .add_tri_face(
                    indices[f[0] as usize],
                    indices[f[1] as usize],
                    indices[f[2] as usize],
                )
                .map_err(|e| format!("Failed to add face: {:?}", e))?;
        }

        Ok(result)
    }
}

pub struct NaAdaptor {}

impl Adaptor<3> for NaAdaptor {
    type Vector = Vector3;
    type Scalar = f64;

    fn vector(coords: [Self::Scalar; 3]) -> Self::Vector {
        Vector3::new(coords[0], coords[1], coords[2])
    }

    fn zero_vector() -> Self::Vector {
        Vector3::zeros()
    }

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
    }
}

impl VectorLengthAdaptor<3> for NaAdaptor {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        v.norm()
    }
}

impl VectorNormalizeAdaptor<3> for NaAdaptor {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for NaAdaptor {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(&b)
    }
}

impl VectorAngleAdaptor for NaAdaptor {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle(&b)
    }
}

impl CrossProductAdaptor for NaAdaptor {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(&b)
    }
}

impl FloatScalarAdaptor<3> for NaAdaptor {
    fn scalarf32(val: f32) -> Self::Scalar {
        val as f64
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val
    }
}

pub type HalfEdgeMesh = alum::PolyMeshT<3, NaAdaptor>;
