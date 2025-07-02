//! This module contains tools for simulating sensors and sensor data

mod laser_profile;

use crate::{Iso3, Mesh, PointCloud};
use parry3d_f64::query::RayCast;

pub use laser_profile::{LaserProfile, PanningLaserProfile};

pub trait SimulatedPointSensor {
    /// Simulate the measurement of a sensor which produces a point cloud. The exact way that the
    /// sensor works is defined by the specific implementation of this trait, and the parameters of
    /// the sensor...including the meaning of the sensor's coordinate system...will vary by the
    /// type of the sensor.
    ///
    /// This function requires a target mesh which is the body which the sensor is intended to
    /// measure. The sensor will produce a point cloud of points which are the result of the
    /// simulated interaction with the body. The sensor may also take an optional obstruction mesh,
    /// which is a mesh which will be included in the scene and is capable of blocking the sensor's
    /// view of the target mesh. The way that the obstruction mesh interferes with the sensor will
    /// vary by the type of sensor and how it captures data; for example, a laser triangulation
    /// sensor will be unable to capture points in which either the detector or the emitter is
    /// blocked by the obstruction mesh, while a laser time of flight sensor will only be blocked
    /// if the laser rays strike the obstruction mesh before they reach the target mesh.
    ///
    /// Finally, the sensor will take an isometry which defines the position and orientation of the
    /// sensor in the world. This isometry moves the sensor coordinate system so that it overlaps
    /// with the area of the world being measured.
    ///
    /// The return value of this function is a tuple containing (1) the simulated point cloud
    /// sampled by the sensor, _in the sensor's coordinate system_, and (2) an optional vector of
    /// standard deviations for each point in the point cloud.
    ///
    /// To get the position of the points in the world coordinate system, you can apply the
    /// isometry to the points in the point cloud.
    ///
    /// The standard deviations are the sigma value for a Gaussian distribution which can be used
    /// to simulate where the point would sample on repeated measurements. This is only provided by
    /// sensors which can model their uncertainty.
    ///
    /// # Arguments
    ///
    /// * `target`: a mesh which the sensor is intended to measure. The sensor will produce a point
    ///   cloud of points and normals which lie on the surface of this mesh.
    /// * `obstruction`: an optional mesh which will be included in the scene and can block the
    ///   sensor's view of the target mesh. The way that the obstruction mesh interferes with the
    ///   sensor will vary by the type of sensor and how it captures data.
    /// * `iso`: an isometry which defines the position and orientation of the sensor in the world.
    ///
    /// returns: (PointCloud, Option<Vec<f64, Global>>)
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>);
}
