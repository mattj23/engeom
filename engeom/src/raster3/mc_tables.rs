//! The marching cubes case table, and the cube numbering it is built on.
//!
//! # Why this is generated rather than written down
//!
//! Marching cubes implementations usually contain a transcribed 256-row table of edge indices.
//! Such a table has approximately four thousand manually placed numbers. One incorrect number can
//! create a hole in a single configuration while the other configurations remain visibly correct.
//! This module derives the table from a rule whose consistency can be demonstrated.
//!
//! # The rule
//!
//! Take a cube whose eight corners are each inside or outside the surface. A cube edge is crossed
//! when its two corners disagree, and the surface meets each of the six faces in some set of
//! segments joining the crossings on that face.
//!
//! Walk the boundary of one face counter-clockwise as seen from outside the cube. Each crossing is
//! either an *entry*, where the walk passes from an outside corner to an inside one, or an *exit*,
//! where it does the reverse; around a closed face the two alternate. Join every entry to the next
//! exit the walk reaches, and direct the segment that way.
//!
//! Two properties establish the consistency of this rule:
//!
//! 1. **Neighboring cells agree on their shared face.** The neighbor walks that same face in the
//!    opposite direction, because consistently oriented faces of a polyhedron traverse a shared
//!    edge in opposite directions. Reversing the walk swaps every entry with every exit, so
//!    "entry joined to the next exit" picks out the same pairs and directs them backwards. The two
//!    cells therefore produce identical segments with opposite winding, as required for a
//!    crack-free surface with consistent orientation. The cells need only agree on their four
//!    shared corner values, which they do by construction.
//!
//! 2. **The segments close into loops.** Every crossed edge belongs to two faces and is traversed
//!    once by each, in opposite directions, so it is an exit on one face and an entry on the other.
//!    That makes it the head of one directed segment and the tail of another, so the segments form
//!    disjoint directed cycles with no loose ends.
//!
//! This module returns the directed cycles. The extractor triangulates them because a fan from the
//! first vertex can fail when neighboring cells share the resulting edges; see
//! [`super::extract_isosurface`]. Any triangulation of a cycle inherits its direction. The winding
//! makes the surface face the *positive* side of the field, matching the sign convention on
//! [`super::SdfVoxel`]: positive is the side the surface normals point toward.
//!
//! # Topology guarantees
//!
//! The rule separates the two diagonal inside corners of an ambiguous face, which is the same
//! choice the original Lorensen and Cline table makes, so the output agrees with the classic
//! tables on the cases where the classic tables are unambiguous.
//!
//! The rule does not guarantee the topologically correct result inside the cube. A face with four
//! crossings does not determine whether the two inside corners join within the cell. Resolving that
//! ambiguity requires the field's behavior in the interior, which marching cubes does not inspect.
//! The rule guarantees a closed, consistently wound surface that separates inside corners from
//! outside corners. A future topology-preserving variant can be implemented by changing only this
//! file.

use std::sync::LazyLock;

/// The position of each cube corner, as a 0 or 1 offset on each axis from the cell's lowest corner.
///
/// This is the numbering used by every published marching cubes table, where the lower and upper
/// faces are each walked as a ring rather than counted in binary.
pub(crate) const CORNER_OFFSETS: [[i32; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];

/// The two corners joined by each of the twelve cube edges.
pub(crate) const EDGE_CORNERS: [(u8, u8); 12] = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

/// Each cube edge as the offset of its lower endpoint from the cell's lowest corner, paired with
/// the axis it runs along.
///
/// This gives an edge an identity independent of the cell that references it. Two cells sharing an
/// edge compute the same offset and axis, so both find the same vertex.
pub(crate) const EDGE_LOWER: [([i32; 3], usize); 12] = [
    ([0, 0, 0], 0),
    ([1, 0, 0], 1),
    ([0, 1, 0], 0),
    ([0, 0, 0], 1),
    ([0, 0, 1], 0),
    ([1, 0, 1], 1),
    ([0, 1, 1], 0),
    ([0, 0, 1], 1),
    ([0, 0, 0], 2),
    ([1, 0, 0], 2),
    ([1, 1, 0], 2),
    ([0, 1, 0], 2),
];

/// The six cube faces, each as its four corners in counter-clockwise order seen from outside the
/// cube. The pairing rule in the module documentation depends on this orientation, so tests verify
/// the orderings.
const FACES: [[u8; 4]; 6] = [
    [0, 3, 2, 1], // z = 0, outward normal -Z
    [4, 5, 6, 7], // z = 1, outward normal +Z
    [0, 1, 5, 4], // y = 0, outward normal -Y
    [2, 3, 7, 6], // y = 1, outward normal +Y
    [0, 4, 7, 3], // x = 0, outward normal -X
    [1, 2, 6, 5], // x = 1, outward normal +X
];

/// The closed loops of crossing edges produced by each of the 256 corner sign configurations.
///
/// The table contains loops without triangulating them because the extractor has the geometric
/// information needed to choose how to fill each loop. See
/// [`super::extract_isosurface`] for why a loop longer than three is not simply fanned.
///
/// Stored flat, with one level of offsets for the loops of each case and another for the edges of
/// each loop, so a lookup during extraction requires two index reads instead of a
/// pointer chase through nested vectors.
pub(crate) struct CaseTable {
    edges: Vec<u8>,
    loop_starts: Vec<u32>,
    case_starts: [u32; 257],
}

/// The largest number of separate loops any corner configuration produces, which bounds the
/// per-cell bookkeeping the extractor allocates. Checked by test.
pub(crate) const MAX_LOOPS: usize = 4;

impl CaseTable {
    /// The number of loops for a corner configuration, where bit `i` of `case` is set when corner
    /// `i` is on the negative side of the surface.
    pub(crate) fn loop_count(&self, case: u8) -> usize {
        (self.case_starts[case as usize + 1] - self.case_starts[case as usize]) as usize
    }

    /// One loop of a corner configuration, as the cube edge indices it passes through in order.
    ///
    /// The loop is directed so that a surface filling it faces the positive side of the field.
    pub(crate) fn loop_edges(&self, case: u8, index: usize) -> &[u8] {
        let slot = self.case_starts[case as usize] as usize + index;
        let start = self.loop_starts[slot] as usize;
        let end = self.loop_starts[slot + 1] as usize;
        &self.edges[start..end]
    }
}

/// The table, built once on first use.
pub(crate) fn case_table() -> &'static CaseTable {
    static TABLE: LazyLock<CaseTable> = LazyLock::new(build_table);
    &TABLE
}

fn build_table() -> CaseTable {
    let mut edges = Vec::new();
    let mut loop_starts = Vec::new();
    let mut case_starts = [0u32; 257];

    for (case, start) in case_starts.iter_mut().enumerate().take(256) {
        *start = loop_starts.len() as u32;

        for cycle in build_case(case as u8) {
            loop_starts.push(edges.len() as u32);
            edges.extend_from_slice(&cycle);
        }
    }
    case_starts[256] = loop_starts.len() as u32;
    loop_starts.push(edges.len() as u32);

    CaseTable {
        edges,
        loop_starts,
        case_starts,
    }
}

/// The index of the edge joining two corners, which must be adjacent on the cube.
fn edge_between(a: u8, b: u8) -> u8 {
    for (i, &(x, y)) in EDGE_CORNERS.iter().enumerate() {
        if (x == a && y == b) || (x == b && y == a) {
            return i as u8;
        }
    }
    panic!("corners {a} and {b} are not joined by a cube edge");
}

/// The loops for one corner configuration, following the rule in the module documentation.
fn build_case(case: u8) -> Vec<Vec<u8>> {
    let inside = |i: u8| (case >> i) & 1 == 1;

    // For each crossed edge, the edge its directed segment runs to. A crossed edge is an entry on
    // one of its two faces and an exit on the other, so each one is written exactly once here.
    let mut next = [u8::MAX; 12];

    for face in FACES.iter() {
        // The crossings met while walking this face counter-clockwise from outside, each tagged
        // with whether the walk was entering the inside region at that point.
        let mut walk: Vec<(u8, bool)> = Vec::new();
        for k in 0..4 {
            let a = face[k];
            let b = face[(k + 1) % 4];
            if inside(a) != inside(b) {
                walk.push((edge_between(a, b), !inside(a) && inside(b)));
            }
        }

        let n = walk.len();
        debug_assert!(n == 0 || n == 2 || n == 4, "a face had {n} crossings");

        for i in 0..n {
            if !walk[i].1 {
                continue;
            }
            for j in 1..=n {
                let (edge, is_entry) = walk[(i + j) % n];
                if !is_entry {
                    next[walk[i].0 as usize] = edge;
                    break;
                }
            }
        }
    }

    let mut loops = Vec::new();
    let mut visited = [false; 12];

    for start in 0..12u8 {
        if visited[start as usize] || next[start as usize] == u8::MAX {
            continue;
        }

        let mut cycle = Vec::new();
        let mut current = start;
        while !visited[current as usize] {
            visited[current as usize] = true;
            cycle.push(current);
            debug_assert!(
                next[current as usize] != u8::MAX,
                "case {case} left a segment with no continuation"
            );
            current = next[current as usize];
        }

        debug_assert!(
            current == start,
            "case {case} produced a path which did not close"
        );
        debug_assert!(cycle.len() >= 3, "case {case} produced a degenerate loop");

        loops.push(cycle);
    }

    loops
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};

    /// The corner offsets of a face, as f64 vectors.
    fn corner(i: u8) -> [f64; 3] {
        let o = CORNER_OFFSETS[i as usize];
        [o[0] as f64, o[1] as f64, o[2] as f64]
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    /// The whole pairing rule rests on the face corner lists running counter-clockwise as seen
    /// from outside the cube. If one of them were wound the other way, that face would hand back
    /// segments pointing backwards and the surface would come out with a reversed patch.
    #[test]
    fn faces_are_wound_counter_clockwise_seen_from_outside() {
        let center = [0.5, 0.5, 0.5];

        for face in FACES.iter() {
            let a = corner(face[0]);
            let b = corner(face[1]);
            let c = corner(face[2]);
            let normal = cross(sub(b, a), sub(c, b));

            // The face plane is the axis on which all four corners agree, and "outside" is away
            // from the middle of the cube.
            let outward = sub(a, center);
            assert!(
                dot(normal, outward) > 0.0,
                "face {face:?} is wound towards the inside of the cube"
            );

            // All four corners must lie in one plane of the cube, which is what makes it a face.
            let fixed = (0..3)
                .filter(|&d| face.iter().all(|&i| corner(i)[d] == a[d]))
                .count();
            assert_eq!(fixed, 1, "face {face:?} is not a single plane of the cube");
        }
    }

    #[test]
    fn faces_cover_every_cube_edge_twice() {
        let mut seen: HashMap<u8, usize> = HashMap::new();

        for face in FACES.iter() {
            for k in 0..4 {
                let edge = edge_between(face[k], face[(k + 1) % 4]);
                *seen.entry(edge).or_insert(0) += 1;
            }
        }

        assert_eq!(seen.len(), 12, "the faces did not touch all twelve edges");
        for edge in 0..12u8 {
            assert_eq!(
                seen.get(&edge).copied().unwrap_or(0),
                2,
                "edge {edge} is not shared by two faces"
            );
        }
    }

    /// Two faces sharing an edge must walk it in opposite directions. That is the property which
    /// makes a crossing an entry on one face and an exit on the other, so the segments link up
    /// head to tail instead of leaving loose ends.
    #[test]
    fn shared_edges_are_walked_in_opposite_directions() {
        let mut seen: HashMap<u8, Vec<(u8, u8)>> = HashMap::new();

        for face in FACES.iter() {
            for k in 0..4 {
                let a = face[k];
                let b = face[(k + 1) % 4];
                seen.entry(edge_between(a, b)).or_default().push((a, b));
            }
        }

        for (edge, walks) in seen.iter() {
            assert_eq!(walks.len(), 2, "edge {edge} was not walked twice");
            assert_eq!(
                walks[0],
                (walks[1].1, walks[1].0),
                "edge {edge} was walked the same way by both faces"
            );
        }
    }

    #[test]
    fn edge_lower_agrees_with_the_corner_numbering() {
        for (edge, &(a, b)) in EDGE_CORNERS.iter().enumerate() {
            let oa = CORNER_OFFSETS[a as usize];
            let ob = CORNER_OFFSETS[b as usize];

            let differing: Vec<usize> = (0..3).filter(|&d| oa[d] != ob[d]).collect();
            assert_eq!(
                differing.len(),
                1,
                "edge {edge} joins corners which are not adjacent"
            );

            let axis = differing[0];
            let lower = if oa[axis] < ob[axis] { oa } else { ob };

            assert_eq!(EDGE_LOWER[edge], (lower, axis), "edge {edge} is mislabeled");
        }
    }

    // ===============================================================================================
    // The generated cases
    // ===============================================================================================

    /// Which cube edges are crossed by a given corner configuration.
    fn crossing_edges(case: u8) -> HashSet<u8> {
        let inside = |i: u8| (case >> i) & 1 == 1;
        (0..12u8)
            .filter(|&e| {
                let (a, b) = EDGE_CORNERS[e as usize];
                inside(a) != inside(b)
            })
            .collect()
    }

    #[test]
    fn a_uniform_cube_produces_nothing() {
        assert_eq!(case_table().loop_count(0), 0);
        assert_eq!(case_table().loop_count(255), 0);
    }

    /// Every loop runs over crossing edges, and consecutive edges of a loop share a cube face,
    /// because each step of a loop is a segment lying in one face.
    #[test]
    fn loops_step_between_edges_which_share_a_face() {
        let face_edges: Vec<Vec<u8>> = FACES
            .iter()
            .map(|f| (0..4).map(|k| edge_between(f[k], f[(k + 1) % 4])).collect())
            .collect();

        for case in 0..=255u8 {
            let crossing = crossing_edges(case);

            for l in 0..case_table().loop_count(case) {
                let cycle = case_table().loop_edges(case, l);
                assert!(
                    cycle.len() >= 3,
                    "case {case} has a loop of {} edges",
                    cycle.len()
                );

                for k in 0..cycle.len() {
                    let a = cycle[k];
                    let b = cycle[(k + 1) % cycle.len()];

                    assert!(crossing.contains(&a), "case {case} loop used edge {a}");
                    assert!(
                        face_edges.iter().any(|f| f.contains(&a) && f.contains(&b)),
                        "case {case}: the step from edge {a} to {b} is not within a face"
                    );
                }
            }
        }
    }

    #[test]
    fn no_case_makes_more_loops_than_the_bookkeeping_allows() {
        for case in 0..=255u8 {
            assert!(
                case_table().loop_count(case) <= MAX_LOOPS,
                "case {case} has {} loops, over the limit of {MAX_LOOPS}",
                case_table().loop_count(case)
            );
        }
    }

    /// The loops of each case must partition the crossing edges, using every one exactly once.
    ///
    /// This is the structural claim the whole construction makes. If it holds for all 256 cases
    /// then every crossing edge is met by one piece of surface, which is what lets neighboring
    /// cells weld into a closed sheet with no edge left hanging and none covered twice.
    #[test]
    fn the_loops_partition_the_crossing_edges() {
        for case in 0..=255u8 {
            let mut seen: Vec<u8> = Vec::new();
            for l in 0..case_table().loop_count(case) {
                seen.extend_from_slice(case_table().loop_edges(case, l));
            }
            seen.sort_unstable();

            let mut expected: Vec<u8> = crossing_edges(case).into_iter().collect();
            expected.sort_unstable();

            assert_eq!(
                seen, expected,
                "case {case}: the loops did not cover the crossing edges once each"
            );
        }
    }

    /// A single inside corner must give one loop of three whose winding faces away from that
    /// corner, which means toward the positive side of the field.
    #[test]
    fn one_inside_corner_gives_one_outward_loop() {
        for c in 0..8u8 {
            let case = 1u8 << c;
            assert_eq!(
                case_table().loop_count(case),
                1,
                "case {case} is not one loop"
            );

            let cycle = case_table().loop_edges(case, 0);
            assert_eq!(cycle.len(), 3, "case {case} should cut three edges");

            // The midpoint of each crossed edge stands in for its vertex.
            let midpoint = |edge: u8| {
                let (a, b) = EDGE_CORNERS[edge as usize];
                let (pa, pb) = (corner(a), corner(b));
                [
                    (pa[0] + pb[0]) * 0.5,
                    (pa[1] + pb[1]) * 0.5,
                    (pa[2] + pb[2]) * 0.5,
                ]
            };

            let (p0, p1, p2) = (midpoint(cycle[0]), midpoint(cycle[1]), midpoint(cycle[2]));
            let normal = cross(sub(p1, p0), sub(p2, p0));

            let away = sub(p0, corner(c));
            assert!(
                dot(normal, away) > 0.0,
                "case {case} produced a loop wound towards the inside corner"
            );
        }
    }
}
