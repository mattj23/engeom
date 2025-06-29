//! This module contains tools for simulating sensors and sensor data

mod laser_triangulation;

use crate::{Iso3, Mesh, PointCloud};
use parry3d_f64::query::RayCast;

pub use laser_triangulation::{LaserLine, PanningLaserLine};

pub trait SimulatedPointSensor {
    fn get_points(
        &self,
        target: &Mesh,
        obstruction: Option<&Mesh>,
        iso: &Iso3,
    ) -> (PointCloud, Option<Vec<f64>>);
}
