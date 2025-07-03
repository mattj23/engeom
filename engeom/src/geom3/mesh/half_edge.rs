use crate::{Mesh, Point3, Vector3};
use alum;
use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
use std::error::Error;

impl TryFrom<HalfEdgeMesh> for Mesh {
    type Error = Box<dyn Error>;

    fn try_from(value: HalfEdgeMesh) -> Result<Self, Self::Error> {
        let borrow_vert = value.points();
        let borrow_vert = borrow_vert
            .try_borrow()
            .map_err(|_| "Failed to borrow points")?;
        let vertices = borrow_vert
            .iter()
            .map(|v| Point3::from(*v))
            .collect::<Vec<_>>();

        let f_status = value.face_status_prop();
        let f_status = f_status
            .try_borrow()
            .map_err(|_| "Failed to borrow face status")?;
        let faces = value
            .triangulated_vertices(&f_status)
            .map(|f| [f[0].index(), f[1].index(), f[2].index()])
            .collect::<Vec<_>>();

        Ok(Mesh::new(vertices, faces, false))
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
