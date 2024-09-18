//! This module contains tools for computing hulls (convex and otherwise) around sets of 2d points
//!

use crate::common::points::dist;
use crate::geom2::{directed_angle, signed_angle, Arc2, Circle2, Iso2, Point2, Vector2};

use crate::common::AngleDir;
use crate::geom2::kd_tree2::{kd_tree_nearest_2d, kd_tree_within_2d, to_kd_tree2};
use crate::Result;
use parry2d_f64::shape::ConvexPolygon;
use parry2d_f64::transformation::convex_hull_idx;
use serde::Serialize;
use std::collections::HashSet;
use std::f64::consts::{FRAC_PI_2, PI};

/// Computes the convex hull of a set of 2d points, returning a vector of `usize` elements that
/// specify the indices of the points in the original set which make up the hull. The indices are
/// ordered in a counter-clockwise direction.
///
/// This is a direct wrapper around the `convex_hull_idx` function from the parry2d crate.
///
/// # Arguments
///
/// * `points`: The 2D points on which to compute the convex hull
///
/// returns: Vec<usize, Global>
pub fn convex_hull_2d(points: &[Point2]) -> Vec<usize> {
    convex_hull_idx(points)
}

/// Finds the indices of the two points in a convex hull which are farthest apart. This is done by
/// calculating the distance between every pair of points in the hull and returning the indices of
/// the pair with the greatest distance. Needs to be replaced with the rotating caliper algorithm.
///
/// # Arguments
///
/// * `hull`: the convex hull for which to find the farthest pair of points
///
/// returns: (usize, usize)
///
/// # Examples
///
/// ```
///
/// ```
pub fn farthest_pair_indices(hull: &ConvexPolygon) -> (usize, usize) {
    // TODO: Replace this with the rotating calipers algorithm
    let mut max_dist = 0.0;
    let mut max_pair = (0, 0);
    for i in 0..hull.points().len() {
        for j in i + 1..hull.points().len() {
            let d = dist(&hull.points()[i], &hull.points()[j]);
            if d > max_dist {
                max_dist = d;
                max_pair = (i, j);
            }
        }
    }

    max_pair
}

/// Estimate the direction of ordering of a set of points.  This is done by calculating and
/// comparing to the convex hull of the points.
///
/// This function computes the convex hull of the points, and then checks if more points on the
/// hull have an index (from the original list) greater than or less than the index of their
/// immediate neighbor. Because the convex hull is always oriented counter-clockwise, ascending
/// indices indicate that the points are also ordered counter-clockwise, and descending indices
/// indicate that the points are ordered clockwise.
///
/// For the result of this function to mean anything, the points ordering within the slice should
/// be meaningful.  Clockwise or counter-clockwise will be in reference to this order and this
/// order alone.
///
/// # Arguments
///
/// * `points`: a slice of points to check for clockwise or counterclockwise order.
///
/// returns: RotationDirection
///
/// # Examples
///
/// ```
///
/// ```
pub fn point_order_direction(points: &[Point2]) -> AngleDir {
    let mut d_sum = 0;
    let hull = convex_hull_2d(points);
    for i in 0..hull.len() {
        let j = (i + 1) % hull.len();
        let d = hull[j] as i32 - hull[i] as i32;
        d_sum += d.signum();
    }

    if d_sum > 0 {
        AngleDir::Ccw
    } else {
        AngleDir::Cw
    }
}

#[derive(Copy, Clone, Debug)]
pub enum BallPivotStart {
    StartOnIndex(usize),
    StartOnIndexDir(usize, Vector2),
    StartOnConvex,
}

#[derive(Copy, Clone, Debug)]
pub enum BallPivotEnd {
    EndOnIndex(usize),
    EndOnRepeat,
}

/// Runs the ball pivoting algorithm on a set of 2d points, generating a list of indices into the
/// points which constitute the outer hull and a list of the center points of the balls which
/// contact each pair of indices. The list of centers is one element shorter than the list of
/// indices.  The center at index `i` is the center of the ball which contacts the points at
/// indices `i` and `i+1`.
///
/// # Arguments
///
/// * `points`:
/// * `start`:
/// * `end`:
/// * `pivot_direction`:
/// * `radius`:
///
/// returns: Result<(Vec<usize, Global>, Vec<OPoint<f64, Const<2>>, Global>)>
///
/// # Examples
///
/// ```
///
/// ```
pub fn ball_pivot_with_centers_2d(
    points: &[Point2],
    start: BallPivotStart,
    end: BallPivotEnd,
    pivot_direction: AngleDir,
    radius: f64,
) -> Result<(Vec<usize>, Vec<Point2>)> {
    // To prepare, we create a KD tree of the points to speed up searching, and a vector of circles
    // with the ball radius. These circles represent the arc of the ball's movement around
    // each point, not the actual ball itself. We use these to find the intersections between the
    // ball trajectories.
    let tree = to_kd_tree2(points);
    let circles = points
        .iter()
        .map(|p| Circle2::new(p.x, p.y, radius))
        .collect::<Vec<_>>();

    // To start, we must have a working index and a vector with the direction from the point at
    // that index to the center of the ball in the starting location.  The vector does not need to
    // be normalized, but it must be non-zero.
    let (start_index, start_direction) = match start {
        BallPivotStart::StartOnIndex(i) => find_start_on_index(points, i, radius)?,
        BallPivotStart::StartOnIndexDir(i, v) => (i, v),
        BallPivotStart::StartOnConvex => {
            let convex = convex_hull_2d(points);
            let v = points[convex[1]] - points[convex[0]];
            (convex[0], Iso2::rotation(-FRAC_PI_2) * v)
        }
    };

    // The `working_index` is the index of the point which currently has the ball, and will be
    // updated as we traverse around the hull. The `results` vector is the final list of indices
    // being returned. The `direction` vector is the direction from the working point to the center
    // of the ball where it currently is, and will also be updated as we traverse around the hull.
    let mut working_index = start_index;
    let mut results = vec![working_index];
    let mut centers = Vec::new();
    let mut direction = start_direction.normalize();

    // We also track a set of completed indices, which we will use if the stop condition is based
    // on the ball returning to a previously visited point.
    let mut completed = HashSet::new();
    completed.insert(working_index);

    // The distance we must search for points is double the radius of the ball, then squared for
    // convenience when working with the kd-tree's squared Euclidean distance function.
    let search2 = (2.0 * radius).powi(2);

    // let mut count = 0;
    loop {
        // Get the neighborhood of points within 2x the radius
        let neighbors = kd_tree_within_2d(&tree, &points[working_index], search2);

        // let mut debug_output = DebugOutput {
        //     points: points.to_vec(),
        //     working_index,
        //     working_direction: direction,
        //     radius,
        //     results: results.clone(),
        //     neighbors: neighbors.clone(),
        //     pivots: Vec::new(),
        // };

        if results.len() > completed.len() * 3 {
            if let BallPivotEnd::EndOnIndex(end_index) = end {
                println!("Should have ended on {}", end_index);
                println!("Total: {}", points.len());
            }
            // json_elements_save("ball-pivot-debug.json".as_ref(), &debug_output).unwrap();
            println!("Results: {:?}", results);
            return Err("Loop detected".into());
        }

        // Now we go through every possible circle in the neighborhood and check all intersections
        // between circles. The one with the intersection that has the smallest positive angle to
        // the last ball contact point is the one we choose to pivot on
        let mut best: Option<PivotPoint> = None;
        for ni in neighbors.iter() {
            // We want to skip the neighbor two elements back, because that's the one we just came
            // from and it will otherwise have a perfect intersection at 0 degrees.
            if results.len() >= 2 && *ni == results[results.len() - 2] {
                continue;
            }

            for pi in circles[working_index].intersections_with(&circles[*ni]) {
                let di = pi - points[working_index];
                let angle = directed_angle(&direction, &di, pivot_direction);
                if angle < 1e-6 {
                    continue;
                }

                let pivot = PivotPoint::new(*ni, pi, angle);
                // debug_output.pivots.push(pivot);
                best = pivot.better_of(best);
            }
        }

        // println!("Best: {:?}", best);

        if let Some(best_item) = best {
            // if best_item.point_index < working_index {
            //     // Temp debugging on known list
            //     println!("{:?}", best_item);
            //     json_elements_save("ball-pivot-debug.json".as_ref(), &debug_output).unwrap();
            //     println!("Results: {:?}", results);
            //     panic!("Loop detected");
            // }
            working_index = best_item.point_index;
            // println!("  * Adding point {}", working_index);
            direction = (best_item.point - points[working_index]).normalize();
            results.push(working_index);
            centers.push(best_item.point);
        } else {
            println!("No intersections");
            println!("Neighbors: {:?}", neighbors);
            break;
        }

        // json_elements_save(
        //     format!("debug_output_{}.json", count).as_ref(),
        //     &debug_output,
        // )
        // .unwrap();
        // count += 1;

        // Finally we check if the end condition is met
        match end {
            BallPivotEnd::EndOnIndex(end_index) => {
                if working_index == end_index {
                    break;
                }
            }
            BallPivotEnd::EndOnRepeat => {
                if completed.contains(&working_index) {
                    break;
                }
            }
        }

        completed.insert(working_index);
    }

    Ok((results, centers))
}

/// Runs the ball pivoting algorithm on a set of 2d points, generating a list of indices into the
/// points which constitute the outer hull. There are several different ways to start and stop the
/// algorithm, and the caller can specify the direction of rotation for the ball to pivot.
///
/// The algorithm is a simplified version of the common ball pivoting algorithm used for creating
/// triangle meshes from point clouds.  The algorithm works by starting with a ball contacting a
/// specified point and at a specified direction from that point.  It then pivots in the specified
/// direction (clockwise or counter-clockwise) until it contacts another point.  The algorithm
/// adds these contact points to the result vector in the order it contacts them, and repeats until
/// the stopping criteria is met.
///
/// # Arguments
///
/// * `points`: the set of points on which the algorithm is run. The result vector will contain
/// indices into this vector.
/// * `start`: The starting condition for the algorithm. This can either be a specific index and
/// direction, or it can be a request to start on one of the convex hull points with a default
/// direction pointing out from the hull.
/// * `stop`: The stopping condition for the algorithm. This can either be a specific index to
/// stop on, or it can be a request to stop when the ball returns to a previously visited point.
/// * `pivot_direction`: The direction in which the ball should pivot. This can either be
/// counter-clockwise or clockwise.
/// * `radius`: The radius of the ball to use for the algorithm. The maximum gap which the
/// algorithm will be able to traverse is twice this radius.
///
/// returns: Vec<usize, Global>
pub fn ball_pivot_2d(
    points: &[Point2],
    start: BallPivotStart,
    end: BallPivotEnd,
    pivot_direction: AngleDir,
    radius: f64,
) -> Result<Vec<usize>> {
    let (result, _) = ball_pivot_with_centers_2d(points, start, end, pivot_direction, radius)?;
    Ok(result)
}

pub fn ball_pivot_fill_gaps_2d(
    points: &[Point2],
    start: BallPivotStart,
    end: BallPivotEnd,
    pivot_direction: AngleDir,
    radius: f64,
    max_spacing: f64,
) -> Result<Vec<Point2>> {
    let (indices, centers) =
        ball_pivot_with_centers_2d(points, start, end, pivot_direction, radius)?;

    // The arc between any two points is a circle, centered on the center point, with the radius
    // equal to the ball radius, starting at the first point and ending at the second point in
    // the opposite direction of the pivot direction.  We will fill in gaps greater than the
    // max gap value by adding points along this arc.

    let mut result = Vec::new();
    for i in 0..centers.len() {
        let p0 = points[indices[i]];
        let p1 = points[indices[i + 1]];
        result.push(p0);
        if dist(&p0, &p1) > max_spacing {
            // Fill in the points
            let v0 = p0 - centers[i];
            let v1 = p1 - centers[i];
            // let angle = directed_angle(&v0, &v1, pivot_direction.opposite());
            let angle = signed_angle(&v0, &v1);
            let arc = Arc2::circle_point_angle(centers[i], radius, p0, angle);

            // How many times do we need to break it up?
            let n = (arc.length() / max_spacing).ceil() as usize;
            let f = 1.0 / n as f64;

            // We skip the first point and the last point
            for j in 1..n {
                let p = arc.point_at_fraction(f * j as f64);
                result.push(p);
            }
        }
    }

    result.push(points[indices[indices.len() - 1]]);

    Ok(result)
}

fn find_start_on_index(points: &[Point2], index: usize, radius: f64) -> Result<(usize, Vector2)> {
    let tree_points = points
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != index)
        .map(|(_, p)| *p)
        .collect::<Vec<_>>();
    let tree = to_kd_tree2(&tree_points);

    let circle = Circle2::from_point(points[index], radius);
    let mut best_distance = 0.0;
    let mut best_angle = 0.0;

    for i in 0..1000 {
        let angle = (i as f64 / 1000.0) * 2.0 * PI;
        let p = circle.point_at_angle(angle);
        let ni = kd_tree_nearest_2d(&tree, &p);
        let d = dist(&p, &tree_points[ni]);
        if d > best_distance {
            best_distance = d;
            best_angle = angle;
        }
    }

    if best_distance < radius {
        Err("Couldn't find a start point".into())
    } else {
        Ok((index, Iso2::rotation(best_angle) * Vector2::new(1.0, 0.0)))
    }
}

// #[derive(Serialize)]
// struct DebugOutput {
//     points: Vec<Point2>,
//     working_index: usize,
//     working_direction: Vector2,
//     radius: f64,
//     results: Vec<usize>,
//     neighbors: Vec<usize>,
//     pivots: Vec<PivotPoint>,
// }

#[derive(Copy, Clone, Serialize, Debug)]
struct PivotPoint {
    point_index: usize,
    point: Point2,
    angle: f64,
}

impl PivotPoint {
    fn new(point_index: usize, point: Point2, angle: f64) -> Self {
        Self {
            point_index,
            point,
            angle,
        }
    }

    fn better_of(&self, other: Option<Self>) -> Option<Self> {
        if let Some(oth) = other {
            if self.angle < oth.angle {
                Some(*self)
            } else {
                Some(oth)
            }
        } else {
            Some(*self)
        }
    }
}
