//! This module contains tools for computing hulls (convex and otherwise) around sets of 2d points
//!

use crate::common::points::dist;
use crate::geom2::{Arc2, Circle2, Iso2, Point2, Vector2, directed_angle, signed_angle};

use crate::common::AngleDir;
use crate::common::kd_tree::KdTreeSearch;
use crate::{KdTree2, Result};
use parry2d_f64::shape::ConvexPolygon;
use parry2d_f64::transformation::convex_hull_idx;
use serde::{Deserialize, Serialize};
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

pub fn farthest_pair_on_hull(all_points: &[Point2], hull_indices: &[usize]) -> (usize, usize) {
    // Find the farthest pair of points on the convex hull
    let mut max_dist = 0.0;
    let mut max_pair = (0, 0);
    for i in 0..hull_indices.len() {
        for j in i + 1..hull_indices.len() {
            let d = dist(&all_points[hull_indices[i]], &all_points[hull_indices[j]]);
            if d > max_dist {
                max_dist = d;
                max_pair = (hull_indices[i], hull_indices[j]);
            }
        }
    }

    max_pair
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
    let tree = KdTree2::new(points);
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

    // The distance we must search for points is double the radius of the ball. The ball center's
    // unimpeded trajectory lies on the circle with radius r centered around the working point, and
    // the farthest point that the ball can sweep is that radius again, so the total distance we
    // need to search for neighbors is twice the radius.
    let search2 = radius * 2.0 + f64::EPSILON;

    // let mut count = 0;
    loop {
        // Get the neighborhood of points within 2x the radius
        let neighbors = tree.within(&points[working_index], search2);

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
        for (ni, _) in neighbors.iter() {
            if ni == &working_index {
                // We don't want to pivot on the same point we are currently at
                continue;
            }

            // We want to skip the neighbor two elements back, because that's the one we just came
            // from, and it will otherwise have a perfect intersection at 0 degrees.
            if results.len() >= 2 && *ni == results[results.len() - 2] {
                continue;
            }

            for pi in circles[working_index].intersections_with(&circles[*ni]) {
                let di = pi - points[working_index];
                let angle = directed_angle(&direction, &di, pivot_direction);
                if angle < f64::EPSILON {
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
///   indices into this vector.
/// * `start`: The starting condition for the algorithm. This can either be a specific index and
///   direction, or it can be a request to start on one of the convex hull points with a default
///   direction pointing out from the hull.
/// * `stop`: The stopping condition for the algorithm. This can either be a specific index to
///   stop on, or it can be a request to stop when the ball returns to a previously visited point.
/// * `pivot_direction`: The direction in which the ball should pivot. This can either be
///   counter-clockwise or clockwise.
/// * `radius`: The radius of the ball to use for the algorithm. The maximum gap which the
///   algorithm will be able to traverse is twice this radius.
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
    let tree = KdTree2::new(&tree_points);

    let circle = Circle2::from_point(points[index], radius);
    let mut best_distance = 0.0;
    let mut best_angle = 0.0;

    for i in 0..1000 {
        let angle = (i as f64 / 1000.0) * 2.0 * PI;
        let p = circle.point_at_angle(angle);
        let (_ni, d) = tree.nearest_one(&p);

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

#[derive(Serialize, Deserialize)]
pub struct DebugBallData {
    pub points: Vec<Point2>,
    pub start_index: usize,
    pub end_index: usize,
    pub direction: AngleDir,
    pub radius: f64,
    pub max_spacing: f64,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_bad_case() {
        let points = known_bad_points();
        let indices = ball_pivot_2d(
            &points,
            BallPivotStart::StartOnIndexDir(0, Vector2::new(-1.0, 0.0)),
            BallPivotEnd::EndOnIndex(points.len() - 1),
            AngleDir::Ccw,
            6.0,
        )
        .unwrap();

        // First point should be zero, the last point should be the last index, and every index
        // should be greater than the previous index.
        assert_eq!(0, indices[0]);
        assert_eq!(points.len() - 1, indices[indices.len() - 1]);

        for i in 1..indices.len() {
            assert!(
                indices[i] > indices[i - 1],
                "Index {} is not greater than previous index {}",
                indices[i],
                indices[i - 1]
            );
        }
    }

    fn known_bad_points() -> Vec<Point2> {
        let copied = vec![
            [-3.0, 0.011533450828884677],
            [-2.9931149848644925, 0.011610177002391808],
            [-2.9432264768244427, 0.012108803082076472],
            [-2.893337968784407, 0.012865090394418202],
            [-2.843449460744357, 0.013755138728188954],
            [-2.793560952704318, 0.013544406581500773],
            [-2.743672444664279, 0.0105300017968558],
            [-2.693783936624225, 0.008028708380988125],
            [-2.6438955439843173, 0.007781413112970176],
            [-2.5940082076871853, 0.011075848108987735],
            [-2.544120871390071, 0.014452923948401364],
            [-2.4942335350929428, 0.010191328654021309],
            [-2.444346198795829, 0.0010423045622110494],
            [-2.3944588624986936, -0.006330606206679162],
            [-2.3445726282685593, -0.003050054725069279],
            [-2.294686568644922, 0.005217642116127255],
            [-2.244800509021263, 0.00918227580768019],
            [-2.194914449397626, 0.005420166699465601],
            [-2.1450283897739744, -0.004586111128981107],
            [-2.0951423301503405, -0.013897272675003208],
            [-2.045256270526689, -0.024676815496321336],
            [-1.9953702109030407, -0.01632336855917753],
            [-1.9454841512794, -0.0058278268102799485],
            [-1.895598091655752, 7.058320980651844e-5],
            [-1.845712032032111, 0.0025999804896405937],
            [-1.7958268410917246, 0.014139947171216083],
            [-1.7459418342170188, 0.020314149788346975],
            [-1.6960568273422991, 0.01917579530621044],
            [-1.6461718204675968, 0.013865636242524987],
            [-1.596286813592891, 0.007982364577035037],
            [-1.5464018067181784, 0.0009267556525584114],
            [-1.4965167998434796, -0.0010509568045465462],
            [-1.4466317929687666, 0.0009217117943644965],
            [-1.396747421605193, 0.007816736948706382],
            [-1.3468636803567602, 0.010406960548908314],
            [-1.2969799391083274, 0.007751817873913477],
            [-1.2470961978599013, 0.007919302395563765],
            [-1.197212456611465, 0.00726236985739924],
            [-1.1473287153630323, 0.005953219623361733],
            [-1.097444974114603, 0.0031078336380426647],
            [-1.047561232866181, 0.003436412103385106],
            [-0.9976780703089783, 0.0007020551574275224],
            [-0.9477952573796458, 0.0018948526720762443],
            [-0.8979124444503102, 0.005201389804361925],
            [-0.848029631520971, 0.008642351666349385],
            [-0.7981468185916487, 0.003923463368246255],
            [-0.7482640056623202, -0.003786303237351643],
            [-0.6983811927329882, -0.00025548576631315156],
            [-0.6484983798036597, 0.007609548951302933],
            [-0.5986155668743303, 0.014658791855104922],
            [-0.5487327539449982, 0.01694921974417797],
            [-0.49884994101565905, 0.025981518754550376],
            [-0.44896726887375227, 0.03779249779993181],
            [-0.39908569721380616, 0.055377070889870435],
            [-0.3492041255538494, 0.07112458674456873],
            [-0.29932255389390683, 0.08248156072899855],
            [-0.2494409822339536, 0.08718439060809569],
            [-0.19956025490192886, 0.0872003335694326],
            [-0.14968119003076197, 0.0857094943076752],
            [-0.09980212515960663, 0.08600812572555475],
            [-0.04992306028844329, 0.09382169294225874],
            [-4.399541728705714e-5, 0.10285869399819128],
            [0.049835069453878944, 0.10545253606650348],
            [0.09971413432503518, 0.09482820921553503],
            [0.14959319919620206, 0.08259851808351136],
            [0.19947226406736096, 0.06978675613971425],
            [0.24935132893852074, 0.05933905132192448],
            [0.299230393809677, 0.048702060955515644],
            [0.3491094586808394, 0.03918785546817684],
            [0.39898852355199566, 0.028509016089007334],
            [0.44886758842316965, 0.016733277758132853],
            [0.49874665329432144, 0.008375830413864499],
            [0.5486257181654919, 0.004384429891011274],
            [0.5985047830366446, 0.006198707624819319],
            [0.6483838479078106, 0.007870290986872877],
            [0.6982629127789739, 0.007058751585470244],
            [0.7481419776501266, -0.00023070977722045027],
            [0.7980210425212926, -0.00891748302942287],
            [0.8479001073924453, -0.013739126437549182],
            [0.8977791722636219, -0.012353303807385363],
            [0.9476582371347746, -0.006991210766181935],
            [0.9975373020059308, -0.0017734416035632],
            [1.047416907674978, 0.004111849062055847],
            [1.0972966681199257, 0.0099146764268815],
            [1.1471764285648742, 0.015186679597431946],
            [1.1970561890098148, 0.013606759162732462],
            [1.2469359494547767, 0.00999729297573669],
            [1.2968157098997288, 0.005926584067701177],
            [1.3466954703446765, 0.0059281093927150835],
            [1.3965752307896349, 0.0014417863367119504],
            [1.4464516549110567, -0.006738314204891807],
            [1.496327993014634, -0.009799971589715947],
            [1.5462043311182008, -0.006894562997243614],
            [1.5960806654413764, -0.0034472546335543724],
            [1.6459566612819563, -0.003429782573289089],
            [1.695832657122546, -0.00024743672278720874],
            [1.7457086529631223, -0.0014260752120152098],
            [1.7955846488037057, -0.004551456930650191],
            [1.8454606446442998, -0.00822393127012791],
            [1.8953366404848762, -0.007025576265045463],
            [1.9452126363254703, -0.008470813109135888],
            [1.9950886321660466, -0.007606837155512414],
            [2.0449646280066327, -0.004463903851894116],
            [2.094840623847216, 0.0003315382136287043],
            [2.144716619687796, -0.0013963791991287722],
            [2.1945926155283937, -0.005867469025177249],
            [2.24446861136897, -0.008199107524823654],
            [2.294344607209564, -0.0034618582317463287],
            [2.344220603050144, 0.003534732496884449],
            [2.394096598890723, 0.010275410632399007],
            [2.4439725947313136, 0.01388569764505404],
            [2.4938485905718935, 0.011780595973343078],
            [2.5437245864124804, 0.004127206034884326],
            [2.5936005822530603, 0.0005184877739736332],
            [2.643476578093651, -0.00048808226547177036],
            [2.6933525739342308, -0.0033951318837023153],
            [2.7432285697748178, -0.004062875510231559],
            [2.793104565615401, -0.0026590133915761344],
            [2.8429805614559838, -0.0017361585014103453],
            [2.8928565572965743, -0.004047216769118039],
            [2.9427325531371507, -0.0006016559385569192],
            [2.9926085489777483, 0.0016243854590760748],
            [3.0, 0.001803219787837871],
        ];

        copied
            .into_iter()
            .map(|p| Point2::new(p[0], p[1]))
            .collect::<Vec<_>>()
    }
}
