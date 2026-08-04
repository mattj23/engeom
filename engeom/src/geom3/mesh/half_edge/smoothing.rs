//! This module implements simple smoothing algorithms for half-edge meshes.

use super::HalfEdgeMesh;
use crate::common::points::dist;
use crate::{Iso3, Point3, Result, SvdBasis3};
use alum::{Handle, HasIterators, HasTopology};

impl HalfEdgeMesh {
    /// Smooth the mesh by fitting a plane to each vertex's one-ring and moving the vertex onto the
    /// average height of that ring within the fitted frame.
    ///
    /// A vertex with fewer than three neighbors is left where it is, since no plane can be fitted
    /// to it. The pass is two-phase, collecting every adjustment before applying any of them, so
    /// the result does not depend on the order the vertices are visited.
    ///
    /// This moves measured points, so it changes the geometry rather than only the topology.
    pub fn neighborhood_smooth(&mut self) -> Result<()> {
        let mesh = &mut self.inner;
        let vertices = mesh.vertices().collect::<Vec<_>>();
        let mut adjusted = Vec::new();

        for vh in vertices {
            let this_point: Point3 = mesh
                .point(vh)
                .map_err(|e| format!("Failed to get point for vertex {}: {:?}", vh.index(), e))?
                .into();

            let mut neighbors = Vec::new();
            for he in mesh.voh_ccw_iter(vh) {
                let neighbor_point: Point3 = mesh
                    .point(he.head(mesh))
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

            let Some(basis) = SvdBasis3::from_points(&n_points, None) else {
                continue;
            };

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
            mesh.set_point(vh, new_point.coords)
                .map_err(|e| format!("Failed to set point for vertex {}: {:?}", vh.index(), e))?;
        }

        Ok(())
    }
}
