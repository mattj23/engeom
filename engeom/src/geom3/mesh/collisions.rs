//! Mesh collision detection and distance checks

use crate::{Iso3, Mesh, Result};
use parry3d_f64::query::intersection_test;
use std::cmp::PartialEq;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq)]
enum MeshType {
    Stationary,
    Moving,
}

struct MeshItem {
    mesh: Mesh,
    mesh_type: MeshType,
}

pub struct MeshCollisionSet {
    meshes: HashMap<usize, MeshItem>,
    exceptions: HashSet<(usize, usize)>,
}

impl MeshCollisionSet {
    pub fn new() -> Self {
        Self {
            meshes: HashMap::new(),
            exceptions: HashSet::new(),
        }
    }

    pub fn add_exception(&mut self, id1: usize, id2: usize) {
        let lower = id1.min(id2);
        let upper = id1.max(id2);
        self.exceptions.insert((lower, upper));
    }

    fn skip_collision(&self, id1: usize, id2: usize) -> bool {
        let lower = id1.min(id2);
        let upper = id1.max(id2);
        self.exceptions.contains(&(lower, upper))
    }

    fn add_mesh(&mut self, mesh: Mesh, mesh_type: MeshType) -> usize {
        let id = self.meshes.len();
        self.meshes.insert(id, MeshItem { mesh, mesh_type });

        id
    }

    pub fn add_stationary(&mut self, mesh: Mesh) -> usize {
        self.add_mesh(mesh, MeshType::Stationary)
    }

    pub fn add_moving(&mut self, mesh: Mesh) -> usize {
        self.add_mesh(mesh, MeshType::Moving)
    }

    /// This function will check for all collisions between the meshes in the set, according to the
    /// following rules:
    ///
    /// - Moving meshes will be checked against all meshes that don't contain an exception,
    ///   including both stationary and other moving meshes
    /// - Stationary meshes will not be checked against any other meshes, and so a collision will
    ///   only be reported if it is with a stationary mesh
    ///
    /// # Arguments
    ///
    /// * `transforms`: transforms for the moving meshes
    /// * `stop_at_first`: If true, the function will stop at the first collision found for each
    ///   moving mesh. If false, it will check all collisions.
    ///
    /// returns: Vec<(usize, usize), Global>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn check_all(
        &self,
        transforms: &[(usize, Iso3)],
        stop_at_first: bool,
    ) -> Result<Vec<(usize, usize)>> {
        // Create the fast isometry lookup:
        let lookups = self.quick_lookups(transforms)?;

        // We'll iterate through all the moving meshes (this can be parallelized later).
        // For each moving mesh, we'll iterate through all stationary meshes and then all moving
        // meshes. Whether a collision check occurs depends on the following:
        // - Is the current mesh id lower than the other mesh id?
        // - Is there an exception for the current pair of meshes?
        let mut pairs = Vec::new();

        for (&id1, mesh1) in self.meshes.iter() {
            if mesh1.mesh_type == MeshType::Stationary {
                continue;
            }

            let iso1 = &lookups[id1];

            for (&id2, mesh2) in self.meshes.iter() {
                if mesh2.mesh_type == MeshType::Moving && id1 >= id2 {
                    continue;
                }

                if self.skip_collision(id1, id2) {
                    continue;
                }

                let iso2 = &lookups[id2];

                // Check for collision
                if let Ok(check) =
                    intersection_test(iso1, mesh1.mesh.tri_mesh(), iso2, mesh2.mesh.tri_mesh())
                {
                    if check {
                        pairs.push((id1, id2));
                        if stop_at_first {
                            break;
                        }
                    }
                }
            }
        }

        Ok(pairs)
    }

    fn quick_lookups(&self, transforms: &[(usize, Iso3)]) -> Result<Vec<Iso3>> {
        let mut lookups = vec![Iso3::identity(); self.meshes.len()];
        for &(id, iso) in transforms.iter() {
            if id >= self.meshes.len() {
                return Err(format!("Transform id {} out of bounds", id).into());
            }

            lookups[id] = iso.clone();
        }

        Ok(lookups)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Iso3;
    use crate::geom3::Mesh;

    #[test]
    fn collision_set() {
        let mut set = MeshCollisionSet::new();
        let mesh1 = Mesh::create_box(1.0, 1.0, 1.0, true);
        let mesh2 = Mesh::create_box(1.0, 1.0, 1.0, true);

        let id1 = set.add_stationary(mesh1);
        let id2 = set.add_moving(mesh2);

        let transforms = vec![(id2, Iso3::translation(0.5, 0.5, 0.5))];

        let pairs = set.check_all(&transforms, false).unwrap();
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn collision_set_exception_skips() {
        let mut set = MeshCollisionSet::new();
        let mesh1 = Mesh::create_box(1.0, 1.0, 1.0, true);
        let mesh2 = Mesh::create_box(1.0, 1.0, 1.0, true);

        let id1 = set.add_stationary(mesh1);
        let id2 = set.add_moving(mesh2);
        set.add_exception(id2, id1);

        let transforms = vec![(id2, Iso3::translation(0.5, 0.5, 0.5))];

        let pairs = set.check_all(&transforms, false).unwrap();
        assert_eq!(pairs.len(), 0);
    }

    #[test]
    fn collision_set_misses() {
        let mut set = MeshCollisionSet::new();
        let mesh1 = Mesh::create_box(1.0, 1.0, 1.0, true);
        let mesh2 = Mesh::create_box(1.0, 1.0, 1.0, true);

        let id1 = set.add_stationary(mesh1);
        let id2 = set.add_moving(mesh2);

        let transforms = vec![(id2, Iso3::translation(2.5, 0.5, 0.5))];

        let pairs = set.check_all(&transforms, false).unwrap();
        assert_eq!(pairs.len(), 0);
    }
}
