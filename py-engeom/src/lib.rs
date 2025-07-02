mod airfoil;
pub mod alignments;
mod bounding;
mod common;
mod conversions;
mod geom2;
mod geom3;
mod mesh;
mod metrology;
mod point_cloud;
mod raster;
mod ray_casting;
mod sensors;
mod svd_basis;

use pyo3::prelude::*;

/// Geometry in 2D space.
fn register_geom2(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_geom2")?;
    // Primitive geometry types
    child.add_class::<geom2::Iso2>()?;
    child.add_class::<geom2::Vector2>()?;
    child.add_class::<geom2::Point2>()?;
    child.add_class::<geom2::SurfacePoint2>()?;
    child.add_class::<geom2::Circle2>()?;

    // Curves and other complex geometries
    child.add_class::<geom2::Curve2>()?;
    child.add_class::<geom2::CurveStation2>()?;

    // Bounding and tools
    child.add_class::<bounding::Aabb2>()?;
    child.add_class::<svd_basis::SvdBasis2>()?;

    parent_module.add_submodule(&child)
}

/// Geometry in 3D space.
fn register_geom3(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_geom3")?;

    // Primitive geometry types
    child.add_class::<geom3::Iso3>()?;
    child.add_class::<geom3::Vector3>()?;
    child.add_class::<geom3::Point3>()?;
    child.add_class::<geom3::Plane3>()?;
    child.add_class::<geom3::SurfacePoint3>()?;

    // Mesh, curves, other complex geometries
    child.add_class::<mesh::Mesh>()?;
    child.add_class::<mesh::MeshCollisionSet>()?;
    child.add_class::<mesh::FaceFilterHandle>()?;
    child.add_class::<geom3::Curve3>()?;
    child.add_class::<geom3::CurveStation3>()?;
    child.add_class::<point_cloud::PointCloud>()?;

    // Bounding and tools
    child.add_class::<bounding::Aabb3>()?;
    child.add_class::<svd_basis::SvdBasis3>()?;

    // Intersection and ray casting
    child.add_class::<ray_casting::RayBundle3>()?;

    parent_module.add_submodule(&child)
}

fn register_align_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_align")?;
    child.add_function(wrap_pyfunction!(alignments::points_to_mesh, &child)?)?;
    child.add_function(wrap_pyfunction!(alignments::points_to_cloud, &child)?)?;
    parent_module.add_submodule(&child)
}

fn register_airfoil_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_airfoil")?;

    child.add_class::<airfoil::MclOrient>()?;
    child.add_class::<airfoil::FaceOrient>()?;
    child.add_class::<airfoil::EdgeFind>()?;
    child.add_class::<airfoil::AfGage>()?;

    child.add_class::<airfoil::InscribedCircle>()?;
    child.add_class::<airfoil::AirfoilGeometry>()?;

    child.add_function(wrap_pyfunction!(
        airfoil::compute_inscribed_circles,
        &child
    )?)?;

    parent_module.add_submodule(&child)
}

fn register_metrology_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_metrology")?;
    child.add_class::<metrology::Distance2>()?;
    child.add_class::<metrology::Distance3>()?;

    parent_module.add_submodule(&child)
}

fn register_raster3_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_raster3")?;

    child.add_function(wrap_pyfunction!(raster::clusters_from_sparse, &child)?)?;

    parent_module.add_submodule(&child)
}

fn register_sensor_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_sensors")?;

    child.add_class::<sensors::LaserProfile>()?;
    child.add_class::<sensors::PanningLaserProfile>()?;

    parent_module.add_submodule(&child)
}

/// Engeom is a library for geometric operations in 2D and 3D space.
#[pymodule(name = "engeom")]
fn py_engeom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // 2D geometry submodule
    register_geom2(m)?;

    // 3D geometry submodule
    register_geom3(m)?;

    // 3D raster module
    register_raster3_module(m)?;

    // Alignment submodule
    register_align_module(m)?;

    // Airfoil submodule
    register_airfoil_module(m)?;

    // Metrology submodule
    register_metrology_module(m)?;

    // Sensor submodule
    register_sensor_module(m)?;

    // Common features and primitives
    m.add_class::<common::DeviationMode>()?;
    m.add_class::<common::Resample>()?;
    m.add_class::<common::SelectOp>()?;

    Ok(())
}
