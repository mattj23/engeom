use crate::common::kd_tree::{KdTreeSearch, MatchedTree};
use crate::common::points::dist;
use crate::common::poisson_disk::sample_poisson_disk;
use crate::{Iso3, KdTree3, Mesh, Point3, Result, SurfacePoint3, UnitVec3};
use bounding_volume::Aabb;
use parry3d_f64::bounding_volume;
use uuid::Uuid;

pub trait PointCloudOverlap<TOther> {
    fn overlap_by_reciprocity(&self, other: &TOther, max_distance: f64) -> Vec<usize>;
}

pub trait PointCloudFeatures {
    fn points(&self) -> &[Point3];
    fn normals(&self) -> Option<&[UnitVec3]>;
    fn colors(&self) -> Option<&[[u8; 3]]>;

    fn is_empty(&self) -> bool {
        self.points().is_empty()
    }

    fn len(&self) -> usize {
        self.points().len()
    }

    fn aabb(&self) -> Aabb {
        Aabb::from_points(self.points())
    }

    fn create_from_indices(&self, indices: &[usize]) -> Result<PointCloud> {
        // Verify that all indices are valid
        if indices.iter().any(|&i| i >= self.len()) {
            return Err("Index out of bounds".into());
        }

        let points = self.points();
        let normals = self.normals();
        let colors = self.colors();

        let points = indices.iter().map(|i| points[*i]).collect();
        let normals = normals.map(|n| indices.iter().map(|i| n[*i]).collect());
        let colors = colors.map(|c| indices.iter().map(|i| c[*i]).collect());

        PointCloud::try_new(points, normals, colors)
    }
}

/// A mutable point cloud with optional normals and colors.
#[derive(Clone)]
pub struct PointCloud {
    tree_uuid: Uuid,
    points: Vec<Point3>,
    normals: Option<Vec<UnitVec3>>,
    colors: Option<Vec<[u8; 3]>>,
}

impl PointCloud {
    /// Create a new point cloud from points and, optionally, normals and colors.
    ///
    /// # Arguments
    ///
    /// * `points`: The points in the point cloud.
    /// * `normals`: Optional normals to be associated with the points. If provided, the number of
    ///   normals must match the number of points.
    /// * `colors`: Optional colors to be associated with the points. If provided, the number of
    ///   colors must match the number of points.
    ///
    /// returns: PointCloud
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn try_new(
        points: Vec<Point3>,
        normals: Option<Vec<UnitVec3>>,
        colors: Option<Vec<[u8; 3]>>,
    ) -> Result<Self> {
        if let Some(normals) = &normals {
            if normals.len() != points.len() {
                return Err("normals must have the same length as points".into());
            }
        }

        if let Some(colors) = &colors {
            if colors.len() != points.len() {
                return Err("colors must have the same length as points".into());
            }
        }
        Ok(Self {
            tree_uuid: Uuid::new_v4(),
            points,
            normals,
            colors,
        })
    }

    pub fn from_surface_points(points: &[SurfacePoint3]) -> Self {
        let normals = points.iter().map(|p| p.normal).collect::<Vec<_>>();
        let points = points.iter().map(|p| p.point).collect();
        Self::try_new(points, Some(normals), None).unwrap()
    }

    /// Merges another point cloud into this one, modifying this point cloud in place and
    /// consuming the other. The two point clouds must either both have normals or both not have
    /// normals, and either both have colors or both not have colors.
    ///
    /// If the point clouds' normal or color data is inconsistent, an error will be returned before
    /// any data is merged, however the other point cloud will still have been moved. Thus, it is
    /// recommended to check the normal and color data of both point clouds before calling this
    /// method.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn merge(&mut self, other: PointCloud) -> Result<()> {
        // Pre-merge checks to ensure that the colors and normals are both either present or absent
        // in both point clouds.
        if self.normals.is_some() != other.normals.is_some() {
            return Err("Cannot merge point clouds with inconsistent normal data".into());
        }
        if self.colors.is_some() != other.colors.is_some() {
            return Err("Cannot merge point clouds with inconsistent color data".into());
        }

        // Merge the points
        self.points.extend(other.points);

        // Merge the normals if they are present
        if let Some(normals) = other.normals {
            self.normals.as_mut().unwrap().extend(normals);
        }

        // Merge the colors if they are present
        if let Some(colors) = other.colors {
            self.colors.as_mut().unwrap().extend(colors);
        }

        self.tree_uuid = Uuid::new_v4();

        Ok(())
    }

    /// Add a single point to the point cloud, along with optional normal and color data. If the
    /// point cloud already has normals the new point must have a normal, and the same goes for
    /// colors. If this consistency check fails, an error will be returned.
    ///
    /// # Arguments
    ///
    /// * `point`: The point to add to the cloud
    /// * `normal`: An optional normal to add to the point, this must be provided if the point
    ///   cloud already has normals and excluded if it does not.
    /// * `color`: An optional color to add to the point, this must be provided if the point cloud
    ///   already has colors and excluded if it does not.
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn append(
        &mut self,
        point: Point3,
        normal: Option<UnitVec3>,
        color: Option<[u8; 3]>,
    ) -> Result<()> {
        // Check that the normal and color data is consistent with the existing point cloud
        if self.normals.is_some() != normal.is_some() {
            return Err("Cannot append point with inconsistent normal data".into());
        }

        if self.colors.is_some() != color.is_some() {
            return Err("Cannot append point with inconsistent color data".into());
        }

        self.points.push(point);
        if let Some(normal) = normal {
            self.normals.as_mut().unwrap().push(normal);
        }
        if let Some(color) = color {
            self.colors.as_mut().unwrap().push(color);
        }

        self.tree_uuid = Uuid::new_v4();
        Ok(())
    }

    /// Create an empty point cloud with the specified normal and color data. The point cloud will
    /// initialize with an empty vector for the points. If `has_normals` is true, an empty vector
    /// will be created for the normals, and the same goes for colors.  Any data appended or merged
    /// into this point cloud must be consistent with the presence/absence of normal and color data.
    ///
    /// # Arguments
    ///
    /// * `has_normals`: if true, the point cloud will have normals
    /// * `has_colors`: if true, the point cloud will have colors
    ///
    /// returns: PointCloud
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn empty(has_normals: bool, has_colors: bool) -> Self {
        Self {
            tree_uuid: Uuid::new_v4(),
            points: Vec::new(),
            normals: if has_normals { Some(Vec::new()) } else { None },
            colors: if has_colors { Some(Vec::new()) } else { None },
        }
    }

    /// Transform the point cloud by applying a transformation to all points and normals. This
    /// modifies the point cloud in place.
    ///
    /// # Arguments
    ///
    /// * `transform`: The transformation to apply to the point cloud.
    ///
    /// returns: ()
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn transform_by(&mut self, transform: &Iso3) {
        for p in &mut self.points {
            *p = transform * *p;
        }

        if let Some(normals) = &mut self.normals {
            for n in normals {
                *n = transform * *n;
            }
        }

        self.tree_uuid = Uuid::new_v4();
    }

    pub fn create_matched_tree(&self) -> Result<MatchedTree<3>> {
        if self.points.is_empty() {
            return Err("Cannot create a KD tree from an empty point cloud".into());
        }
        let tree = KdTree3::new(&self.points);
        Ok(MatchedTree::new(self.tree_uuid, tree))
    }
}

impl TryFrom<(&[Point3], &[UnitVec3])> for PointCloud {
    type Error = Box<dyn std::error::Error>;

    fn try_from(value: (&[Point3], &[UnitVec3])) -> Result<Self> {
        let (points, normals) = value;
        if points.len() != normals.len() {
            return Err("points and normals must have the same length".into());
        }

        Self::try_new(points.to_vec(), Some(normals.to_vec()), None)
    }
}

impl From<&[Point3]> for PointCloud {
    fn from(points: &[Point3]) -> Self {
        Self::try_new(points.to_vec(), None, None)
            .expect("Failed to create point cloud from points, this should not happen")
    }
}

impl From<&[SurfacePoint3]> for PointCloud {
    fn from(points: &[SurfacePoint3]) -> Self {
        let normals = points.iter().map(|p| p.normal).collect::<Vec<_>>();
        let points = points.iter().map(|p| p.point).collect();
        Self::try_new(points, Some(normals), None)
            .expect("Points and normals must have the same length, this should not have happened")
    }
}

impl PointCloudFeatures for PointCloud {
    fn points(&self) -> &[Point3] {
        &self.points
    }

    fn normals(&self) -> Option<&[UnitVec3]> {
        self.normals.as_deref()
    }

    fn colors(&self) -> Option<&[[u8; 3]]> {
        self.colors.as_deref()
    }
}

pub struct PointCloudKdTree<'a> {
    cloud: &'a PointCloud,
    tree: &'a MatchedTree<3>,
}

impl<'a> PointCloudKdTree<'a> {
    pub fn try_new(cloud: &'a PointCloud, tree: &'a MatchedTree<3>) -> Result<Self> {
        if cloud.tree_uuid != tree.tree_uuid() {
            return Err("The point cloud and the KD tree do not match".into());
        }
        Ok(Self { cloud, tree })
    }

    pub fn tree(&self) -> &KdTree3 {
        self.tree.tree()
    }

    /// Performs a Poisson disk sampling on the point cloud, returning a vector of indices of points
    /// which are at least `radius` distance apart from each other.
    ///
    /// # Arguments
    ///
    /// * `radius`: The minimum distance between sampled points.
    ///
    /// returns: Vec<usize, Global>
    pub fn sample_poisson_disk(&self, radius: f64) -> Vec<usize> {
        // This extra work is because it seems that kiddo does not necessarily return all points
        // within the radius.
        let mut mask = vec![true; self.cloud.points.len()];
        let mut previous_pass = Vec::new();
        for i in 0..self.cloud.points.len() {
            if !mask[i] {
                continue;
            }
            previous_pass.push(i);
            let neighbors = self.tree().within(&self.cloud.points[i], radius);
            for (n, _) in neighbors {
                mask[n] = false;
            }
        }

        let mut next_pass = sample_poisson_disk(self.cloud.points(), &previous_pass, radius);
        while next_pass.len() != previous_pass.len() {
            previous_pass = next_pass;
            next_pass = sample_poisson_disk(self.cloud.points(), &previous_pass, radius);
        }

        next_pass
    }

    /// Create a new point cloud from a Poisson disk sampling of the original point cloud. The new
    /// point cloud will contain only points that are at least `radius` distance apart from each
    /// other.
    ///
    /// This is the equivalent of performing `create_from_indices` using the indices returned by
    /// `sample_poisson_disk`.
    ///
    /// # Arguments
    ///
    /// * `radius`: The minimum distance between sampled points.
    ///
    /// returns: Result<PointCloud, Box<dyn Error, Global>>
    pub fn create_from_poisson_sample(&self, radius: f64) -> Result<PointCloud> {
        let indices = self.sample_poisson_disk(radius);
        self.cloud.create_from_indices(&indices)
    }
}

impl PointCloudOverlap<Mesh> for PointCloudKdTree<'_> {
    /// Find the indices of points in this point cloud that "overlap" with a mesh by looking for
    /// reciprocity in the closest point in each direction.
    ///
    /// For each point in this point cloud "p_this", we will find the closest point in the mesh
    /// "p_other".  Then we take "p_other" and find the closest point to it in this point cloud,
    /// "p_recip".
    ///
    /// In an ideally overlapping point cloud, "p_recip" should be the same as "p_this".  We will
    /// use a maximum distance tolerance to determine if "p_recip" is close enough to "p_this" that
    /// "p_this" is considered to be overlapping with the mesh.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh to check for overlap with
    /// * `max_distance`: the maximum distance between "p_this" and "p_recip" for them to be
    ///   considered overlapping
    ///
    /// returns: Vec<usize, Global>
    fn overlap_by_reciprocity(&self, mesh: &Mesh, max_distance: f64) -> Vec<usize> {
        let mut result = Vec::new();
        for (i, p_this) in self.cloud.points.iter().enumerate() {
            // Find the closest point in the mesh
            let p_other = mesh.point_closest_to(p_this);

            // Find the reciprocal point in this point cloud
            let (k, _) = self.tree().nearest_one(&p_other);
            let p_recip = self.cloud.points()[k];

            if dist(&p_this, &p_recip) < max_distance {
                result.push(i);
            }
        }

        result
    }
}

impl PointCloudOverlap<PointCloudKdTree<'_>> for PointCloudKdTree<'_> {
    /// Find the indices of points in this point cloud that "overlap" with points in another point
    /// cloud by looking for reciprocity in the closest point in each direction.
    ///
    /// For each point in this point cloud "p_this", we will find the closest point in the other
    /// point cloud "p_other".  Then we take "p_other" and find the closest point to it in this
    /// point cloud, "p_recip".
    ///
    /// In an ideally overlapping point cloud, "p_recip" should be the same as "p_this".  We will
    /// use a maximum distance tolerance to determine if "p_recip" is close enough to "p_this" that
    /// "p_this" is considered to be overlapping with the other point cloud.
    ///
    /// # Arguments
    ///
    /// * `other`: the other point cloud to check for overlap with
    /// * `max_distance`: the maximum distance between "p_this" and "p_recip" for them to be
    ///   considered overlapping
    ///
    /// returns: Vec<usize, Global>
    fn overlap_by_reciprocity(&self, other: &PointCloudKdTree, max_distance: f64) -> Vec<usize> {
        let mut result = Vec::new();
        for (i, p_this) in self.cloud.points.iter().enumerate() {
            // Find the closest point in the other point cloud
            let (j, _) = other.tree().nearest_one(p_this);
            let p_other = other.cloud.points()[j];

            // Find the reciprocal point in this point cloud
            let (k, _) = self.tree().nearest_one(&p_other);
            let p_recip = self.cloud.points()[k];

            if dist(&p_this, &p_recip) < max_distance {
                result.push(i);
            }
        }

        result
    }
}

impl PointCloudFeatures for PointCloudKdTree<'_> {
    fn points(&self) -> &[Point3] {
        &self.cloud.points
    }

    fn normals(&self) -> Option<&[UnitVec3]> {
        self.cloud.normals.as_deref()
    }

    fn colors(&self) -> Option<&[[u8; 3]]> {
        self.cloud.colors.as_deref()
    }
}
