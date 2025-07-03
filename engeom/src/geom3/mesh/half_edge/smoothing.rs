//! This module implements simple smoothing algorithms for half-edge meshes.

use super::HalfEdgeMesh;
use crate::common::points::dist;
use crate::{Iso3, Point3, Result, SvdBasis3};
use alum::{Handle, HasIterators, HasTopology};

pub trait HalfEdgeSmoothing {
    fn neighborhood_smooth(&mut self) -> Result<()>;
}

impl HalfEdgeSmoothing for HalfEdgeMesh {
    fn neighborhood_smooth(&mut self) -> Result<()> {
        let vertices = self.vertices().collect::<Vec<_>>();
        let mut adjusted = Vec::new();

        for vh in vertices {
            let this_point: Point3 = self
                .point(vh)
                .map_err(|e| format!("Failed to get point for vertex {}: {:?}", vh.index(), e))?
                .into();

            let mut neighbors = Vec::new();
            for he in self.voh_ccw_iter(vh) {
                let neighbor_point: Point3 = self
                    .point(he.head(self))
                    .map_err(|e| {
                        format!("Failed to get point for half-edge {}: {:?}", he.index(), e)
                    })?
                    .into();

                neighbors.push((neighbor_point, dist(&this_point, &neighbor_point)));
            }

            if neighbors.len() < 3 {
                continue; // No neighbors to smooth with
            }

            // Naive smoothing to start with
            let mut n_points: Vec<Point3> = neighbors.iter().map(|(p, _)| *p).collect();
            n_points.push(this_point); // Include the current point

            let basis = SvdBasis3::from_points(&n_points, None);
            let t = Iso3::from(&basis);

            let transformed = n_points.iter().map(|p| t * p).collect::<Vec<_>>();

            // Average the z values
            let avg_z = transformed.iter().map(|p| p.z).sum::<f64>() / transformed.len() as f64;
            let t_point = t * this_point;
            let new_point = t.inverse() * Point3::new(t_point.x, t_point.y, avg_z);
            adjusted.push((vh, new_point));
        }

        // Apply the adjustments
        for (vh, new_point) in adjusted {
            self.set_point(vh, new_point.coords)
                .map_err(|e| format!("Failed to set point for vertex {}: {:?}", vh.index(), e))?;
        }

        Ok(())
    }
}
