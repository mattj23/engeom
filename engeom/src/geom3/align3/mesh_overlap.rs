// //! This module contains a set of helpful tools to determine the overlap between sampled alignment
// //! points on meshes.
//
// use crate::Result;
// use crate::na::DMatrix;
// use std::collections::HashMap;
//
// pub struct CloudOverlap {
//     pub matrix: DMatrix<f64>,
//     pub indices: HashMap<(usize, usize), Vec<usize>>,
// }
//
// impl CloudOverlap {
//     fn new(matrix: DMatrix<f64>, indices: HashMap<(usize, usize), Vec<usize>>) -> Self {
//         Self { matrix, indices }
//     }
// }
//
// pub fn compute_sampled_overlap_matrix(
//     clouds: &[PointCloudWithTree],
//     max_distance: f64,
//     sample_radius: f64,
//     transforms: Option<&[Iso3]>,
// ) -> Result<CloudOverlap> {
//     // Compute the poisson disk sampling of each cloud, resulting in a collection of vecs containing
//     // sampled indices for each cloud.
//     let sampled = clouds
//         .iter()
//         .map(|cloud| cloud.sample_poisson_disk(sample_radius, None))
//         .collect::<Vec<_>>();
//
//     compute_overlap_matrix(clouds, max_distance, Some(&sampled), transforms)
// }
//
// pub fn compute_overlap_matrix(
//     clouds: &[PointCloudWithTree],
//     max_distance: f64,
//     indices: Option<&[Vec<usize>]>,
//     transforms: Option<&[Iso3]>,
// ) -> Result<CloudOverlap> {
//     // If indices to use were provided then we will retain them, otherwise we will generate
//     // temporary indices for each cloud.
//     let indices = match indices {
//         Some(indices) => indices.to_vec(),
//         None => clouds
//             .iter()
//             .map(|cloud| (0..cloud.len()).collect())
//             .collect::<Vec<_>>(),
//     };
//
//     // If transforms were provided then we will retain them, otherwise we will use the identity
//     // transform for each cloud.
//     let transforms = match transforms {
//         Some(transforms) => transforms.to_vec(),
//         None => vec![Iso3::identity(); clouds.len()],
//     };
//
//     let mut result = CloudOverlap::new(DMatrix::zeros(clouds.len(), clouds.len()), HashMap::new());
//
//     let threshold = max_distance * max_distance;
//
//     // We are going to iterate through every combination of clouds and find the indices of the
//     // points in the test cloud that are within the threshold distance of the reference cloud. The
//     // outer loop iterating through the reference clouds, while the inner loop is iterating through
//     // the collections of test indices in the test clouds.
//     for (i, reference) in clouds.iter().enumerate() {
//         for (j, test_indices) in indices.iter().enumerate() {
//             if i == j {
//                 continue;
//             }
//
//             // At this point `i` is the index of the reference cloud and `j` is the index of the
//             // test cloud. The `reference` variable is the `PointCloudWithTree` entity itself and
//             // the `test_indices` is a vector of `usize` indices into the test cloud representing
//             // the points that were specified to look for a match at.
//             let test_cloud = &clouds[j];
//
//             // Get the relative transform between the clouds, allowing a point in the test cloud
//             // to be transformed to the reference cloud with the same distances between points as
//             // if both clouds were transformed by their respective transforms.
//             let t = transforms[i].inv_mul(&transforms[j]);
//
//             let mut indices_with_matches = Vec::new();
//             for &index in test_indices {
//                 let test_point = t * test_cloud.points()[index];
//                 let (d, _) = reference.nearest(&test_point);
//                 if d < threshold {
//                     indices_with_matches.push(index);
//                 }
//             }
//
//             result.matrix[(i, j)] = indices_with_matches.len() as f64;
//             result.indices.insert((i, j), indices_with_matches);
//         }
//     }
//
//     // println!("overlap matrix: {}", result.matrix);
//
//     Ok(result)
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::geom3::point_cloud::test::stanford_bunny;
//     use approx::assert_relative_eq;
//     use parry3d_f64::na::{Translation3, UnitQuaternion};
//
//     /// This is the expected overlap matrix for the bunny clouds at 0.05m max distance and 0.001m
//     /// sample radius.  Because the sampling is random, the exact values will vary slightly
//     fn expected_overlap() -> DMatrix<f64> {
//         DMatrix::from_row_slice(
//             3,
//             3,
//             &[
//                 0.0, 11515.0, 8730.0, 11896.0, 0.0, 9096.0, 11743.0, 11515.0, 0.0,
//             ],
//         )
//     }
//
//     #[test]
//     fn test_bunny_overlaps() {
//         let bunny = stanford_bunny();
//         let meshes = vec![
//             bunny["bun000"].clone().into_with_tree(),
//             bunny["bun045"].clone().into_with_tree(),
//             bunny["bun090"].clone().into_with_tree(),
//         ];
//
//         let overlap = compute_sampled_overlap_matrix(&meshes, 0.05, 0.001, None).unwrap();
//         println!("{}", overlap.matrix);
//
//         let expected = expected_overlap();
//         assert_relative_eq!(overlap.matrix, expected, epsilon = 100.0);
//     }
//
//     #[test]
//     fn test_bunny_with_transforms() {
//         let bunny = stanford_bunny();
//         let mut meshes = vec![
//             bunny["bun000"].clone(),
//             bunny["bun045"].clone(),
//             bunny["bun090"].clone(),
//         ];
//
//         let transforms = vec![
//             Iso3::from_parts(
//                 Translation3::new(0.12, -0.34, 0.21),
//                 UnitQuaternion::from_euler_angles(-0.03, 0.11, 0.5),
//             ),
//             Iso3::from_parts(
//                 Translation3::new(-0.23, 0.12, 0.34),
//                 UnitQuaternion::from_euler_angles(0.12, -0.5, 0.11),
//             ),
//             Iso3::from_parts(
//                 Translation3::new(0.34, 0.21, -0.12),
//                 UnitQuaternion::from_euler_angles(0.5, 0.11, -0.03),
//             ),
//         ];
//
//         for (mesh, transform) in meshes.iter_mut().zip(transforms.iter()) {
//             mesh.transform(&transform.inverse());
//         }
//
//         let meshes = meshes
//             .into_iter()
//             .map(|mesh| mesh.into_with_tree())
//             .collect::<Vec<_>>();
//         let overlap =
//             compute_sampled_overlap_matrix(&meshes, 0.05, 0.001, Some(&transforms)).unwrap();
//
//         let expected = expected_overlap();
//         assert_relative_eq!(overlap.matrix, expected, epsilon = 100.0);
//     }
// }
