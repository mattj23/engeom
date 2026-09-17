//! Choosing which way each normal of a point cloud faces.
//!
//! A plane fit through a neighborhood gives an axis without a direction, so newly estimated
//! normals can point arbitrarily inward or outward. Any downstream operation that uses a normal to
//! identify the outside of a surface requires consistent directions. This is especially important
//! for a signed distance field, where a reversed normal can create a region with the wrong sign and
//! an unwanted surface in the mesh.
//!
//! Direction can come from a known viewpoint or from propagation between neighboring points. A
//! known viewpoint provides the stronger constraint and should be preferred when available.
//!
//! # From a viewpoint, when one is known
//!
//! A measured surface point faces the sensor because the sensor observes the visible side of the
//! surface. If the sensor position is known, [`orient_by_viewpoint`] can determine the direction
//! with one comparison per point. [`orient_by_viewpoints`] applies the same rule to a moving sensor
//! by accepting one sensor position per point.
//!
//! # By propagation, when one is not
//!
//! Without a viewpoint, orientation relies on local consistency. Neighboring points on a smooth
//! surface have nearly parallel normals, so the direction selected for one point constrains its
//! neighbors. [`orient_by_propagation`] implements the method from Hoppe's 1992 paper.
//!
//! Propagation contains two heuristic choices. First, a minimum spanning tree prioritizes the most
//! nearly parallel pairs. This routes propagation around sharp edges, where neighboring normals
//! genuinely disagree and orientation errors are more likely. Second, each connected component
//! needs a seed with no previously oriented neighbor. The algorithm uses the component's topmost
//! point and directs its normal upward. This typically produces outward normals for a closed shape
//! in front of a scanner, but an arbitrary point set can receive either global orientation. An
//! incorrectly oriented seed reverses the entire component.
//!
//! Check the component count in [`OrientationReport3`]. Each component has an independently chosen
//! seed, so a large component count increases the chance that some components have the wrong
//! global orientation.
//!
//! Propagation also requires sampling that is fine relative to the surface curvature. If adjacent
//! normals differ by a large angle, there is no consistent direction to propagate. The result can
//! then be nearly arbitrary even when the algorithm reports one connected component and no error.
//! In the tests, a sphere sampled at one-tenth of its radius is oriented correctly, while the same
//! sphere sampled at one-half of its radius is not.

use crate::common::kd_tree::KdTreeSearch;
use crate::{Point3, Result, UnitVec3, Vector3};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Where the direction of a set of estimated normals should come from.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum NormalOrientation3 {
    /// Face every normal toward one fixed sensor position.
    Viewpoint(Point3),

    /// Face each normal toward its corresponding sensor position for a scan taken while the sensor
    /// was moving.
    /// The array must have one entry per point.
    Viewpoints(Vec<Point3>),

    /// Settle directions by propagating agreement between neighbors when no viewpoint is known.
    /// `k` is the number of neighbors linked to each point. Twelve is a reasonable default; see
    /// [`orient_by_propagation`].
    Propagate { k: usize },
}

/// Statistics from a normal-orientation pass.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OrientationReport3 {
    /// The number of normals that were reversed.
    pub flipped: usize,

    /// The number of connected components that propagation seeded independently, or `None` when the
    /// directions came from a viewpoint and no graph was built. See the module documentation for
    /// why a large number here is a warning.
    pub components: Option<usize>,
}

/// The vector from each point toward a viewpoint, in the form
/// [`super::normal_estimation::estimate_by_neighborhood`] wants.
pub fn must_match_toward(points: &[Point3], viewpoint: &Point3) -> Vec<Vector3> {
    points.iter().map(|p| viewpoint - p).collect()
}

/// Face every normal toward a viewpoint and return the number reversed.
///
/// A normal perpendicular to the line of sight is unchanged because both directions have the same
/// dot product with the view vector.
///
/// # Arguments
///
/// * `points`: the positions the normals belong to
/// * `normals`: the normals to orient, modified in place
/// * `viewpoint`: the position the surface was measured from
///
/// returns: the number of normals reversed
///
/// # Panics
///
/// If `normals` is not the same length as `points`.
pub fn orient_by_viewpoint(
    points: &[Point3],
    normals: &mut [UnitVec3],
    viewpoint: &Point3,
) -> usize {
    assert_eq!(
        points.len(),
        normals.len(),
        "there must be one normal per point"
    );

    let mut flipped = 0;
    for (p, n) in points.iter().zip(normals.iter_mut()) {
        if n.dot(&(viewpoint - p)) < 0.0 {
            *n = -*n;
            flipped += 1;
        }
    }
    flipped
}

/// Face each normal toward its corresponding viewpoint for a scan taken while the sensor moved.
///
/// # Arguments
///
/// * `points`: the positions the normals belong to
/// * `normals`: the normals to orient, modified in place
/// * `viewpoints`: the position each point was measured from, one per point
///
/// returns: `Result<usize>`, the number of normals reversed, failing if the lengths disagree
pub fn orient_by_viewpoints(
    points: &[Point3],
    normals: &mut [UnitVec3],
    viewpoints: &[Point3],
) -> Result<usize> {
    if points.len() != normals.len() {
        return Err(format!(
            "there are {} points but {} normals",
            points.len(),
            normals.len()
        )
        .into());
    }

    if viewpoints.len() != points.len() {
        return Err(format!(
            "there are {} points but {} viewpoints",
            points.len(),
            viewpoints.len()
        )
        .into());
    }

    let mut flipped = 0;
    for ((p, n), v) in points.iter().zip(normals.iter_mut()).zip(viewpoints.iter()) {
        if n.dot(&(v - p)) < 0.0 {
            *n = -*n;
            flipped += 1;
        }
    }
    Ok(flipped)
}

/// Settle the direction of every normal by propagating agreement outward from a seed.
///
/// Points are linked to their `k` nearest neighbors, each link weighted by `1 - |na . nb|` so that
/// pairs that already nearly agree weigh least. A minimum spanning tree over that graph is walked
/// from a seed, and each normal is reversed where it disagrees with the one it was reached from.
/// Taking the cheapest available link at each step delays crossing a sharp edge until no smoother
/// route remains, which reduces the risk of propagating an incorrect direction across the edge.
///
/// The seed of each connected piece is its topmost point, with its normal turned to face up. See
/// the module documentation for the limitations of this seed rule.
///
/// # Arguments
///
/// * `points`: the positions the normals belong to
/// * `normals`: the normals to orient, modified in place
/// * `tree`: a k-d tree built over `points`
/// * `k`: how many neighbors to link each point to. Larger values connect sparse data at the cost
///   of memory and of links long enough to bridge between separate surfaces. Internally, the value
///   is clamped to `max(1, min(k, point_count - 1))`, with saturating subtraction. Values of zero
///   and values larger than the cloud are therefore accepted.
///
/// returns: `Result<OrientationReport3>`, failing if the normal count differs from the point count
/// or the cloud contains more points than a `u32` index can represent
///
/// # Memory
///
/// The neighbor graph holds roughly `3 * k` four-byte indices per point while it is being built.
/// At the default `k` that is around a gigabyte for ten million points, so a cloud of that size is
/// better reduced first with `PointCloud3::reduce_by_voxel`, which also improves the normals.
pub fn orient_by_propagation(
    points: &[Point3],
    normals: &mut [UnitVec3],
    tree: &(impl KdTreeSearch<3> + Sync),
    k: usize,
) -> Result<OrientationReport3> {
    if points.len() != normals.len() {
        return Err(format!(
            "there are {} points but {} normals",
            points.len(),
            normals.len()
        )
        .into());
    }

    if points.is_empty() {
        return Ok(OrientationReport3 {
            flipped: 0,
            components: Some(0),
        });
    }

    if points.len() > u32::MAX as usize {
        return Err(format!(
            "normal propagation indexes points with u32 and cannot take {} of them",
            points.len()
        )
        .into());
    }

    let graph = build_neighbor_graph(points, tree, k);

    // Points from the top down, which is the order seeds are taken in.
    let mut by_height: Vec<u32> = (0..points.len() as u32).collect();
    by_height.sort_unstable_by(|&a, &b| {
        points[b as usize]
            .z
            .total_cmp(&points[a as usize].z)
            .then_with(|| a.cmp(&b))
    });

    let mut visited = vec![false; points.len()];
    let mut heap: BinaryHeap<Link> = BinaryHeap::new();
    let mut report = OrientationReport3 {
        flipped: 0,
        components: Some(0),
    };
    let mut cursor = 0;

    loop {
        // The topmost point not yet reached seeds the next piece.
        while cursor < by_height.len() && visited[by_height[cursor] as usize] {
            cursor += 1;
        }
        if cursor >= by_height.len() {
            break;
        }

        let seed = by_height[cursor] as usize;
        report.components = report.components.map(|c| c + 1);

        if normals[seed].z < 0.0 {
            normals[seed] = -normals[seed];
            report.flipped += 1;
        }

        visited[seed] = true;
        heap.clear();
        push_links(&graph, normals, seed, &mut heap);

        while let Some(link) = heap.pop() {
            let node = link.node as usize;
            if visited[node] {
                continue;
            }
            visited[node] = true;

            if normals[link.parent as usize].dot(&normals[node]) < 0.0 {
                normals[node] = -normals[node];
                report.flipped += 1;
            }

            push_links(&graph, normals, node, &mut heap);
        }
    }

    Ok(report)
}

/// A candidate step of the spanning tree, ordered so that the cheapest comes off the heap first.
struct Link {
    cost: f64,
    parent: u32,
    node: u32,
}

impl PartialEq for Link {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Link {}

impl PartialOrd for Link {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Link {
    /// Reversed, because `BinaryHeap` is a max-heap and this wants the smallest cost. The node
    /// index breaks ties so that equal costs come off in a fixed order, making the result
    /// independent of the heap's internal arrangement.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.node.cmp(&self.node))
            .then_with(|| other.parent.cmp(&self.parent))
    }
}

fn push_links(
    graph: &NeighborGraph,
    normals: &[UnitVec3],
    from: usize,
    heap: &mut BinaryHeap<Link>,
) {
    for &to in graph.neighbors(from) {
        heap.push(Link {
            // Parallel axes have low cost, and perpendicular axes have high cost. The absolute
            // value makes the cost independent of the directions being selected.
            cost: 1.0 - normals[from].dot(&normals[to as usize]).abs(),
            parent: from as u32,
            node: to,
        });
    }
}

/// The k-nearest-neighbor graph, made symmetric and stored as one flat index array with per-point
/// offsets.
struct NeighborGraph {
    offsets: Vec<u32>,
    indices: Vec<u32>,
}

impl NeighborGraph {
    fn neighbors(&self, i: usize) -> &[u32] {
        &self.indices[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }
}

/// Build the symmetric k-nearest-neighbor graph.
///
/// Nearest-neighbor is not a symmetric relation: a point in a dense cluster can be the neighbor of
/// a distant outlier without the outlier being one of its own. Using the one-directional graph
/// would leave such a point unreachable and create a separate component with an independently
/// chosen seed. Adding the reverse of every link prevents this at the cost of storing both
/// directions.
fn build_neighbor_graph(
    points: &[Point3],
    tree: &(impl KdTreeSearch<3> + Sync),
    k: usize,
) -> NeighborGraph {
    let n = points.len();
    let k = k.min(n.saturating_sub(1)).max(1);

    // The nearest neighbors of each point, padded to a fixed stride so that the whole thing is one
    // allocation. A query can return fewer than asked for, and it returns the point itself, which
    // is dropped here.
    let mut raw = vec![u32::MAX; n * k];
    raw.par_chunks_mut(k).enumerate().for_each(|(i, slot)| {
        let found = tree.nearest(&points[i], k + 1);
        let mut written = 0;
        for (j, _) in found {
            if j == i || written >= k {
                continue;
            }
            slot[written] = j as u32;
            written += 1;
        }
    });

    let mut counts = vec![0u32; n + 1];
    for (i, chunk) in raw.chunks(k).enumerate() {
        for &j in chunk {
            if j == u32::MAX {
                continue;
            }
            counts[i] += 1;
            counts[j as usize] += 1;
        }
    }

    let mut offsets = vec![0u32; n + 1];
    let mut total = 0u32;
    for i in 0..n {
        offsets[i] = total;
        total += counts[i];
    }
    offsets[n] = total;

    let mut indices = vec![0u32; total as usize];
    let mut fill = offsets.clone();
    for (i, chunk) in raw.chunks(k).enumerate() {
        for &j in chunk {
            if j == u32::MAX {
                continue;
            }
            indices[fill[i] as usize] = j;
            fill[i] += 1;
            indices[fill[j as usize] as usize] = i as u32;
            fill[j as usize] += 1;
        }
    }
    drop(raw);

    // A link present in both directions was added twice. Duplicates only cost redundant heap
    // pushes, but there are a lot of them, so they are worth dropping while the graph is laid out.
    let mut compacted: Vec<u32> = Vec::with_capacity(indices.len());
    let mut new_offsets = vec![0u32; n + 1];
    for i in 0..n {
        new_offsets[i] = compacted.len() as u32;

        let start = offsets[i] as usize;
        let end = offsets[i + 1] as usize;
        let slice = &mut indices[start..end];
        slice.sort_unstable();

        let mut last = u32::MAX;
        for &j in slice.iter() {
            if j != last {
                compacted.push(j);
                last = j;
            }
        }
    }
    new_offsets[n] = compacted.len() as u32;

    NeighborGraph {
        offsets: new_offsets,
        indices: compacted,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{KdTree3, Mesh3, PointCloud3};

    /// A deterministic generator, so that a failure can be reproduced.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        fn coin(&mut self) -> bool {
            self.next_u64() & 1 == 1
        }
    }

    /// A sphere sampled with its true outward normals.
    fn sphere_cloud(radius: f64, spacing: f64) -> PointCloud3 {
        Mesh3::create_sphere(radius, radius * 0.002)
            .expect("sphere creation failed")
            .sample_poisson(spacing, None)
            .expect("sampling failed")
    }

    /// Reverse a deterministic scattering of normals, returning how many were turned.
    fn scramble(normals: &mut [UnitVec3], seed: u64) -> usize {
        let mut rng = Rng(seed);
        let mut turned = 0;
        for n in normals.iter_mut() {
            if rng.coin() {
                *n = -*n;
                turned += 1;
            }
        }
        turned
    }

    /// How many normals agree in direction with a reference set.
    fn agreement(a: &[UnitVec3], b: &[UnitVec3]) -> f64 {
        let matching = a
            .iter()
            .zip(b.iter())
            .filter(|(x, y)| x.dot(y) > 0.0)
            .count();
        matching as f64 / a.len() as f64
    }

    // ===========================================================================================
    // Viewpoint orientation
    // ===========================================================================================

    #[test]
    fn a_viewpoint_turns_every_normal_towards_it() {
        let cloud = sphere_cloud(5.0, 0.4);
        let truth = cloud
            .point_normals()
            .expect("the sample carries normals")
            .to_vec();

        let mut normals = truth.clone();
        let turned = scramble(&mut normals, 12345);
        assert!(turned > 0, "the scramble did nothing");

        // The viewpoint is far outside the sphere. Normals face the viewpoint rather than pointing
        // outward across the whole cloud, so only the near side should match the true normals.
        let viewpoint = Point3::new(0.0, 0.0, 100.0);
        orient_by_viewpoint(cloud.points(), &mut normals, &viewpoint);

        for (p, n) in cloud.points().iter().zip(normals.iter()) {
            assert!(
                n.dot(&(viewpoint - p)) >= 0.0,
                "a normal at {p:?} still faces away from the viewpoint"
            );
        }
    }

    /// A viewpoint at the center of a sphere faces every normal inward, so the reversed count must
    /// equal the number that initially pointed outward.
    #[test]
    fn the_reversed_count_is_the_number_that_disagreed() {
        let cloud = sphere_cloud(5.0, 0.5);
        let truth = cloud.point_normals().expect("normals").to_vec();

        let mut normals = truth.clone();
        let center = Point3::origin();
        let expected = normals
            .iter()
            .zip(cloud.points())
            .filter(|(n, p)| n.dot(&(center - *p)) < 0.0)
            .count();

        let flipped = orient_by_viewpoint(cloud.points(), &mut normals, &center);
        assert_eq!(flipped, expected);

        // Every normal now points toward the center, opposite to the sampled normals.
        assert_eq!(agreement(&normals, &truth), 0.0);
    }

    #[test]
    fn per_point_viewpoints_must_match_the_cloud() {
        let cloud = sphere_cloud(3.0, 1.0);
        let mut normals = cloud.point_normals().expect("normals").to_vec();

        let short = vec![Point3::origin(); cloud.points().len() - 1];
        assert!(orient_by_viewpoints(cloud.points(), &mut normals, &short).is_err());

        let long = vec![Point3::origin(); cloud.points().len() + 1];
        assert!(orient_by_viewpoints(cloud.points(), &mut normals, &long).is_err());

        let right = vec![Point3::origin(); cloud.points().len()];
        assert!(orient_by_viewpoints(cloud.points(), &mut normals, &right).is_ok());
    }

    /// Per-point viewpoints can orient an entire closed shape outward when every viewpoint lies
    /// outside its corresponding point.
    #[test]
    fn per_point_viewpoints_can_orient_a_closed_shape() {
        let cloud = sphere_cloud(5.0, 0.4);
        let truth = cloud.point_normals().expect("normals").to_vec();

        let mut normals = truth.clone();
        scramble(&mut normals, 777);

        // Each viewpoint is twice its point's distance from the center, outside the sphere and on
        // the same side as the point.
        let viewpoints: Vec<Point3> = cloud.points().iter().map(|p| *p * 2.0).collect();
        orient_by_viewpoints(cloud.points(), &mut normals, &viewpoints)
            .expect("orientation failed");

        assert_eq!(agreement(&normals, &truth), 1.0);
    }

    #[test]
    fn must_match_toward_points_at_the_viewpoint() {
        let points = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(0.0, 0.0, -3.0),
        ];
        let viewpoint = Point3::new(0.0, 0.0, 10.0);

        let vectors = must_match_toward(&points, &viewpoint);
        assert_eq!(vectors.len(), points.len());
        for (p, v) in points.iter().zip(vectors.iter()) {
            assert_eq!(*v, viewpoint - p);
        }
    }

    // ===========================================================================================
    // Propagation
    // ===========================================================================================

    #[test]
    fn propagation_over_nothing_is_nothing() {
        let points: Vec<Point3> = Vec::new();
        let mut normals: Vec<UnitVec3> = Vec::new();
        let tree = KdTree3::try_new(&points).expect("tree failed");

        let report = orient_by_propagation(&points, &mut normals, &tree, 12).expect("failed");
        assert_eq!(report.flipped, 0);
        assert_eq!(report.components, Some(0));
    }

    #[test]
    fn propagation_rejects_a_normal_count_that_does_not_match() {
        let points = vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)];
        let mut normals = vec![UnitVec3::new_normalize(Vector3::z())];
        let tree = KdTree3::try_new(&points).expect("tree failed");

        assert!(orient_by_propagation(&points, &mut normals, &tree, 4).is_err());
    }

    /// A neighbor count larger than the cloud must be clamped because a caller using a default
    /// cannot know the size of every cloud.
    ///
    /// This test does not assert orientation accuracy. The cloud is small enough that adjacent
    /// normals turn sharply, so it does not satisfy propagation's smoothness precondition. Measured
    /// agreement on this fixture is 0.54, close to random orientation. This is expected behavior;
    /// [`propagation_recovers_a_scrambled_sphere`] covers data that satisfies the precondition.
    #[test]
    fn propagation_clamps_a_neighbor_count_larger_than_the_cloud() {
        let cloud = sphere_cloud(2.0, 1.2);
        let truth = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        let mut normals = truth.clone();
        scramble(&mut normals, 99);

        let report = orient_by_propagation(cloud.points(), &mut normals, &tree, 100_000)
            .expect("propagation failed");

        assert_eq!(report.components, Some(1));
        assert_eq!(normals.len(), truth.len());

        // Whatever it decided, every normal is still a unit vector along the original axis.
        for (n, t) in normals.iter().zip(truth.iter()) {
            assert!((n.dot(t).abs() - 1.0).abs() < 1e-12);
        }
    }

    /// The main claim: propagation recovers a consistent orientation over a smooth closed surface
    /// whose normals were scrambled.
    #[test]
    fn propagation_recovers_a_scrambled_sphere() {
        let cloud = sphere_cloud(5.0, 0.3);
        let truth = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        for seed in [1u64, 2, 3, 4, 5] {
            let mut normals = truth.clone();
            scramble(&mut normals, seed);

            let report = orient_by_propagation(cloud.points(), &mut normals, &tree, 12)
                .expect("propagation failed");

            assert_eq!(
                report.components,
                Some(1),
                "seed {seed}: a sampled sphere should be one connected piece"
            );

            // The seed is the topmost point with its normal turned up, which on a sphere is the
            // outward direction, so the recovered set should match the sampled normals rather than
            // their reverse.
            // Measured at 1.0000 for every one of twenty seeds on this fixture, where the
            // sampling is fine enough relative to the curvature for propagation.
            let matched = agreement(&normals, &truth);
            assert!(
                matched > 0.995,
                "seed {seed}: only {:.4} of the normals came back right",
                matched
            );
        }
    }

    /// Two separated shapes cannot share neighbor links, so propagation must seed and report each
    /// component independently.
    #[test]
    fn separate_shapes_are_reported_as_separate_pieces() {
        let near = sphere_cloud(2.0, 0.3);
        let mut points: Vec<Point3> = near.points().to_vec();
        let mut truth: Vec<UnitVec3> = near.point_normals().expect("normals").to_vec();

        // The same sphere again, far enough away that no neighbor link can reach across.
        let offset = Vector3::new(1000.0, 0.0, 0.0);
        points.extend(near.points().iter().map(|p| p + offset));
        truth.extend(near.point_normals().expect("normals").iter().copied());

        let tree = KdTree3::try_new(&points).expect("tree failed");
        let mut normals = truth.clone();
        scramble(&mut normals, 4242);

        let report =
            orient_by_propagation(&points, &mut normals, &tree, 12).expect("propagation failed");

        assert_eq!(report.components, Some(2));

        // Each component must be internally consistent. Their global directions can differ because
        // their seeds are selected independently, so this test does not compare the components.
        let half = near.points().len();
        for (piece, range) in [(0, 0..half), (1, half..points.len())] {
            let matched = agreement(&normals[range.clone()], &truth[range]);
            assert!(
                matched > 0.995 || matched < 0.005,
                "piece {piece} came back mixed, at {matched}"
            );
        }
    }

    /// A scanned shape with sharp creases and thin features that are absent from a sphere.
    #[test]
    fn propagation_recovers_the_stanford_bunny() {
        let cloud = crate::tests::stanford_bun_3()
            .sample_poisson(0.002, None)
            .expect("sampling failed");
        let truth = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        let mut normals = truth.clone();
        scramble(&mut normals, 20260916);

        let report = orient_by_propagation(cloud.points(), &mut normals, &tree, 12)
            .expect("propagation failed");

        let matched = agreement(&normals, &truth);

        // Measured at 0.9825 on this fixture, over a single connected piece. The points that come
        // back wrong occur where a neighbor across a thin gap is as close as a neighbor along the
        // surface, particularly at the ears and the underside of the base.
        //
        // Agreement decreases as sampling becomes coarser, as predicted by the module's sampling
        // precondition. The same bunny at twice the spacing reaches only 0.9454. Use a viewpoint
        // when every normal must have the correct global direction.
        assert!(
            matched > 0.97,
            "only {matched:.4} of the bunny normals came back right, over {:?} pieces",
            report.components
        );
        assert_eq!(report.components, Some(1));
    }

    /// Propagation must be independent of the heap's ordering of equal-cost links.
    #[test]
    fn propagation_is_reproducible() {
        let cloud = sphere_cloud(4.0, 0.35);
        let truth = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        let run = || {
            let mut normals = truth.clone();
            scramble(&mut normals, 31337);
            let report = orient_by_propagation(cloud.points(), &mut normals, &tree, 12)
                .expect("propagation failed");
            (normals, report)
        };

        let (first, first_report) = run();
        let (second, second_report) = run();

        assert_eq!(first_report, second_report);
        assert_eq!(first, second);
    }

    // ===========================================================================================
    // The combined entry point
    // ===========================================================================================

    #[test]
    fn estimating_with_a_viewpoint_gives_outward_normals_on_the_near_side() {
        let cloud = sphere_cloud(5.0, 0.25);
        let index = cloud.compute_index().expect("index failed");

        let viewpoint = Point3::new(0.0, 0.0, 50.0);
        let (estimates, report) = index
            .estimate_normals_oriented(0.6, &NormalOrientation3::Viewpoint(viewpoint))
            .expect("estimation failed");

        assert_eq!(estimates.normals.len(), cloud.points().len());
        assert_eq!(report.components, None);

        for (p, n) in cloud.points().iter().zip(estimates.normals.iter()) {
            assert!(n.dot(&(viewpoint - p)) >= 0.0);
        }
    }

    #[test]
    fn estimating_with_propagation_agrees_with_the_sampled_normals() {
        let cloud = sphere_cloud(5.0, 0.25);
        let truth = cloud.point_normals().expect("normals").to_vec();
        let index = cloud.compute_index().expect("index failed");

        let (estimates, report) = index
            .estimate_normals_oriented(0.6, &NormalOrientation3::Propagate { k: 12 })
            .expect("estimation failed");

        assert_eq!(report.components, Some(1));

        // These are fitted normals rather than the sampled ones, so they differ a little in angle;
        // this test verifies that they all have the same orientation.
        let matched = agreement(&estimates.normals, &truth);
        assert!(matched > 0.995, "only {matched:.4} agreed");
    }

    #[test]
    fn estimating_with_bad_viewpoints_fails() {
        let cloud = sphere_cloud(3.0, 1.0);
        let index = cloud.compute_index().expect("index failed");

        let wrong = vec![Point3::origin(); cloud.points().len() + 3];
        assert!(
            index
                .estimate_normals_oriented(1.5, &NormalOrientation3::Viewpoints(wrong))
                .is_err()
        );
    }

    /// The older entry point takes directions from a caller-supplied array. Separating the axis fit
    /// must preserve that behavior.
    #[test]
    fn the_must_match_entry_point_still_follows_the_vectors_it_is_given() {
        let cloud = sphere_cloud(5.0, 0.4);
        let index = cloud.compute_index().expect("index failed");

        let outward: Vec<Vector3> = cloud.points().iter().map(|p| p.coords).collect();
        let estimates = index
            .estimate_normals(&outward, 0.8)
            .expect("estimation failed");

        assert_eq!(estimates.normals.len(), cloud.points().len());
        assert_eq!(estimates.confidence.len(), cloud.points().len());

        for (n, m) in estimates.normals.iter().zip(outward.iter()) {
            assert!(
                n.dot(m) >= 0.0,
                "a normal did not follow its must_match vector"
            );
        }

        // Reversing the requested directions reverses every normal while preserving confidence.
        let inward: Vec<Vector3> = outward.iter().map(|v| -v).collect();
        let reversed = index
            .estimate_normals(&inward, 0.8)
            .expect("estimation failed");

        for (a, b) in estimates.normals.iter().zip(reversed.normals.iter()) {
            assert!(
                (a.dot(b) + 1.0).abs() < 1e-12,
                "the two runs are not opposites"
            );
        }
        assert_eq!(estimates.confidence, reversed.confidence);
    }

    #[test]
    fn estimate_normals_rejects_a_must_match_of_the_wrong_length() {
        let cloud = sphere_cloud(3.0, 1.0);
        let index = cloud.compute_index().expect("index failed");

        let wrong = vec![Vector3::z(); cloud.points().len() + 1];
        assert!(index.estimate_normals(&wrong, 1.5).is_err());
    }
}
