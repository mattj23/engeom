mod airfoil2;
mod align3;
pub mod alignments;
mod boundary2;
mod bounding;
mod common;
mod conversions;
mod geom2;
mod geom3;
mod mesh;
mod metrology;
mod point_cloud;
mod raster;
mod raster2;
mod ray_casting;
mod sensors;
mod svd_basis;

use pyo3::prelude::*;

/// Raster in 2D space.
fn register_raster2(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_raster2")?;

    // Primitive geometry types
    child.add_class::<raster2::ScalarRaster>()?;

    parent_module.add_submodule(&child)
}

/// Geometry in 2D space.
fn register_geom2(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_geom2")?;
    // Primitive geometry types
    child.add_class::<geom2::Iso2>()?;
    child.add_class::<geom2::Vector2>()?;
    child.add_class::<geom2::Point2>()?;
    child.add_class::<geom2::SurfacePoint2>()?;
    child.add_class::<geom2::Circle2>()?;
    child.add_class::<geom2::Arc2>()?;
    child.add_class::<geom2::Segment2>()?;
    child.add_class::<geom2::Line2>()?;
    child.add_class::<geom2::CubicSpline2>()?;
    child.add_class::<geom2::CubicSplineQueries2>()?;
    child.add_class::<geom2::SplineProjection>()?;
    child.add_function(wrap_pyfunction!(geom2::fit_spline_to_points, &child)?)?;

    // Angle functions
    child.add_function(wrap_pyfunction!(geom2::rot90, &child)?)?;
    child.add_function(wrap_pyfunction!(geom2::rot270, &child)?)?;
    child.add_function(wrap_pyfunction!(geom2::signed_angle, &child)?)?;
    child.add_function(wrap_pyfunction!(geom2::directed_angle, &child)?)?;
    child.add_function(wrap_pyfunction!(geom2::convex_hull_2d, &child)?)?;
    child.add_function(wrap_pyfunction!(geom2::convex_hull_idx, &child)?)?;

    // Curves and other complex geometries
    child.add_class::<geom2::Curve2>()?;
    child.add_class::<geom2::CurveStation2>()?;

    // Boundary geometry
    child.add_class::<boundary2::Manifold1Pos2>()?;
    child.add_class::<boundary2::BoundaryData2>()?;
    child.add_class::<boundary2::Boundary2>()?;
    child.add_function(wrap_pyfunction!(boundary2::fit_boundary_to_points, &child)?)?;
    child.add_function(wrap_pyfunction!(
        boundary2::fit_boundary_to_surface_points,
        &child
    )?)?;

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
    child.add_class::<geom3::Line3>()?;
    child.add_class::<geom3::Segment3>()?;
    child.add_class::<geom3::Sphere3>()?;
    child.add_class::<geom3::Manifold1Pos3>()?;
    child.add_class::<geom3::Circle3>()?;
    child.add_class::<geom3::Cylinder3>()?;
    child.add_class::<geom3::Cone3>()?;

    // Mesh, curves, other complex geometries
    child.add_class::<mesh::Mesh>()?;
    child.add_class::<mesh::MeshCollisionSet>()?;
    child.add_class::<mesh::FaceFilterHandle>()?;
    child.add_class::<geom3::Curve3>()?;
    child.add_class::<geom3::CurveStation3>()?;
    child.add_class::<geom3::CubicSpline3>()?;
    child.add_class::<geom3::CubicSplineQueries3>()?;
    child.add_class::<geom2::SplineProjection>()?;
    child.add_function(wrap_pyfunction!(geom3::fit_spline_to_points, &child)?)?;
    child.add_class::<point_cloud::PointCloud>()?;

    // Bounding and tools
    child.add_class::<bounding::Aabb3>()?;
    child.add_class::<svd_basis::SvdBasis3>()?;

    // Intersection and ray casting
    child.add_class::<ray_casting::RayBundle3>()?;

    parent_module.add_submodule(&child)
}

fn register_align3_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_align3")?;
    child.add_class::<align3::Dof6>()?;
    child.add_class::<align3::AlignParams3>()?;
    child.add_class::<align3::Alignment3>()?;
    child.add_function(wrap_pyfunction!(align3::points_to_mesh, &child)?)?;
    parent_module.add_submodule(&child)
}

fn register_airfoil2_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_airfoil2")?;

    child.add_class::<airfoil2::Inscribed>()?;
    child.add_class::<airfoil2::AfEdgeGeometry>()?;
    child.add_class::<airfoil2::AfEdge>()?;
    child.add_class::<airfoil2::AfEdgeFit>()?;
    child.add_class::<airfoil2::OrientFwdAft>()?;
    child.add_class::<airfoil2::OrientUpperLower>()?;
    child.add_class::<airfoil2::AfEdgeSearch>()?;
    child.add_class::<airfoil2::AfGeometry>()?;

    child.add_function(wrap_pyfunction!(
        airfoil2::extract_inscribed_circles,
        &child
    )?)?;
    child.add_function(wrap_pyfunction!(airfoil2::fit_square_edge, &child)?)?;
    child.add_function(wrap_pyfunction!(airfoil2::fit_rounded_square_edge, &child)?)?;
    child.add_function(wrap_pyfunction!(airfoil2::fit_sharp_edge, &child)?)?;
    child.add_function(wrap_pyfunction!(airfoil2::fit_full_round_edge, &child)?)?;
    child.add_function(wrap_pyfunction!(airfoil2::fit_blended_round_edge, &child)?)?;

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

fn register_common_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let child = PyModule::new(parent_module.py(), "_common")?;

    child.add_class::<common::AngleDir>()?;
    child.add_class::<common::AngleInterval>()?;
    child.add_function(wrap_pyfunction!(common::angle_in_direction, &child)?)?;
    child.add_function(wrap_pyfunction!(common::shortest_angle_between, &child)?)?;
    child.add_function(wrap_pyfunction!(common::angle_signed_pi, &child)?)?;
    child.add_function(wrap_pyfunction!(common::angle_to_2pi, &child)?)?;
    child.add_function(wrap_pyfunction!(common::signed_compliment_2pi, &child)?)?;

    parent_module.add_submodule(&child)
}

/// Engeom is a library for geometric operations in 2D and 3D space.
#[pymodule(name = "engeom")]
fn py_engeom(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // 2D geometry submodule
    register_geom2(m)?;

    // 2D raster submodule
    register_raster2(m)?;

    // 3D geometry submodule
    register_geom3(m)?;

    // 3D raster module
    register_raster3_module(m)?;

    // Alignment submodule
    register_align3_module(m)?;

    // Airfoil2 submodule
    register_airfoil2_module(m)?;

    // Metrology submodule
    register_metrology_module(m)?;

    // Sensor submodule
    register_sensor_module(m)?;

    // Common submodule
    register_common_module(m)?;

    // Common features and primitives
    m.add_class::<common::DeviationMode>()?;
    m.add_class::<common::Resample>()?;
    m.add_class::<common::SelectOp>()?;
    m.add_class::<common::VecDot>()?;

    Ok(())
}
