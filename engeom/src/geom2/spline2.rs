use crate::AngleDir::Cw;
use crate::common::PCoords;
use crate::geom2::{Aabb2, BoundaryElement, ManifoldPosition2, rot90};
use crate::{Point2, UnitVec2, Vector2};
use serde::{Deserialize, Serialize};

// 10-point Gauss-Legendre quadrature nodes and weights mapped to [0, 1].
// Derived from the standard [-1, 1] rule via x = (t + 1) / 2, w_new = w / 2.
const GL_NODES: [f64; 10] = [
    0.013046735741414139,
    0.067468316655507744,
    0.160295215850487796,
    0.283302302935376404,
    0.425562830509184394,
    0.574437169490815606,
    0.716697697064623596,
    0.839704784149512204,
    0.932531683344492256,
    0.986953264258585861,
];

const GL_WEIGHTS: [f64; 10] = [
    0.033335672154344069,
    0.074725674575290296,
    0.109543181257991022,
    0.134633359654998177,
    0.147762112357376435,
    0.147762112357376435,
    0.134633359654998177,
    0.109543181257991022,
    0.074725674575290296,
    0.033335672154344069,
];

/// One cubic segment of the spline: S(τ) = a + b·τ + c·τ² + d·τ³ for τ ∈ [0, 1].
#[derive(Clone, Debug, Serialize, Deserialize)]
struct SplineSegment {
    a: Vector2,
    b: Vector2,
    c: Vector2,
    d: Vector2,
}

impl SplineSegment {
    fn eval(&self, tau: f64) -> Point2 {
        let tau2 = tau * tau;
        let tau3 = tau2 * tau;
        Point2::from(self.a + self.b * tau + self.c * tau2 + self.d * tau3)
    }

    fn deriv(&self, tau: f64) -> Vector2 {
        self.b + self.c * (2.0 * tau) + self.d * (3.0 * tau * tau)
    }

    fn deriv2(&self, tau: f64) -> Vector2 {
        self.c * 2.0 + self.d * (6.0 * tau)
    }

    /// Arc length from τ=0 to τ=`end` via scaled Gauss-Legendre quadrature.
    fn arc_length_to(&self, end: f64) -> f64 {
        GL_NODES
            .iter()
            .zip(GL_WEIGHTS.iter())
            .map(|(&x, &w)| w * self.deriv(x * end).norm() * end)
            .sum()
    }

    /// Total arc length of this segment (τ ∈ [0, 1]).
    fn total_arc_length(&self) -> f64 {
        GL_NODES
            .iter()
            .zip(GL_WEIGHTS.iter())
            .map(|(&x, &w)| w * self.deriv(x).norm())
            .sum()
    }
}

/// A natural cubic spline interpolating through a set of 2D control points.
///
/// Given n+1 control points the spline is composed of n cubic segments with C2 continuity at each
/// interior knot and zero second-derivative boundary conditions at the two endpoints (natural
/// spline). The segments are parameterised uniformly in the global parameter `t ∈ [0, 1]`, where
/// `t = 0` is the first control point and `t = 1` is the last.
///
/// Arc-length queries are available through the [`BoundaryElement`] trait; the `at_t` method
/// provides direct access to the normalized parameter.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CubicSpline2 {
    segments: Vec<SplineSegment>,
    /// Cumulative arc length at the start of each segment, plus the total at the end.
    /// Length = n_segments + 1, with `cumulative_lengths[0] == 0`.
    cumulative_lengths: Vec<f64>,
    total_length: f64,
}

impl CubicSpline2 {
    /// Create a natural cubic spline interpolating through `points`.
    ///
    /// Returns `None` if fewer than 2 points are supplied.
    pub fn try_new(points: &[Point2]) -> Option<Self> {
        let n = points.len();
        if n < 2 {
            return None;
        }

        // Compute second derivatives σ at each knot (natural boundary: σ_0 = σ_{n-1} = 0).
        let sigma = compute_sigma(points);

        // Build segment coefficients.
        // For segment i: S_i(τ) = a + b·τ + c·τ² + d·τ³
        //   a = p_i
        //   b = p_{i+1} - p_i - (2σ_i + σ_{i+1}) / 6
        //   c = σ_i / 2
        //   d = (σ_{i+1} - σ_i) / 6
        let n_segs = n - 1;
        let segments: Vec<SplineSegment> = (0..n_segs)
            .map(|i| {
                let a = points[i].coords;
                let b = points[i + 1].coords - points[i].coords
                    - (sigma[i] * 2.0 + sigma[i + 1]) / 6.0;
                let c = sigma[i] / 2.0;
                let d = (sigma[i + 1] - sigma[i]) / 6.0;
                SplineSegment { a, b, c, d }
            })
            .collect();

        // Precompute cumulative arc lengths.
        let mut cumulative_lengths = vec![0.0_f64; n_segs + 1];
        for i in 0..n_segs {
            cumulative_lengths[i + 1] =
                cumulative_lengths[i] + segments[i].total_arc_length();
        }
        let total_length = *cumulative_lengths.last().unwrap();

        Some(Self {
            segments,
            cumulative_lengths,
            total_length,
        })
    }

    /// Number of cubic segments.
    fn n_segs(&self) -> usize {
        self.segments.len()
    }

    /// Map global `t ∈ [0, 1]` to `(segment_index, local_τ ∈ [0, 1])`.
    fn t_to_seg(&self, t: f64) -> (usize, f64) {
        let n = self.n_segs();
        let scaled = t * n as f64;
        let idx = (scaled.floor() as usize).min(n - 1);
        (idx, scaled - idx as f64)
    }

    /// Evaluate the spline at normalized parameter `t ∈ [0, 1]`.
    ///
    /// `t = 0` corresponds to the first control point; `t = 1` to the last. The `l` field of the
    /// returned [`ManifoldPosition2`] holds the arc length at `t`.
    pub fn at_t(&self, t: f64) -> ManifoldPosition2 {
        let t = t.clamp(0.0, 1.0);
        let (idx, tau) = self.t_to_seg(t);
        self.position_at(idx, tau)
    }

    /// Build a [`ManifoldPosition2`] for segment `idx` at local parameter `tau`.
    fn position_at(&self, idx: usize, tau: f64) -> ManifoldPosition2 {
        let seg = &self.segments[idx];
        let point = seg.eval(tau);
        let deriv = seg.deriv(tau);

        // Guard against degenerate (zero-speed) points.
        let direction = if deriv.norm_squared() > 1e-24 {
            UnitVec2::new_normalize(deriv)
        } else {
            // Fall back to the segment chord direction.
            let chord = Point2::from(self.segments[idx].a + self.segments[idx].b)
                - Point2::from(self.segments[idx].a);
            UnitVec2::new_normalize(chord)
        };

        let normal = rot90(Cw) * direction;
        let l = self.cumulative_lengths[idx] + seg.arc_length_to(tau);
        ManifoldPosition2::new(l, point, direction, normal)
    }

    /// Within segment `seg_idx`, find the local `τ` whose accumulated arc length from `τ = 0`
    /// equals `target` using Newton's method.
    fn tau_at_seg_length(&self, seg_idx: usize, target: f64) -> f64 {
        let seg_len =
            self.cumulative_lengths[seg_idx + 1] - self.cumulative_lengths[seg_idx];
        if target <= 0.0 {
            return 0.0;
        }
        if target >= seg_len {
            return 1.0;
        }

        let seg = &self.segments[seg_idx];
        let mut tau = target / seg_len; // linear initial guess

        for _ in 0..20 {
            let current_len = seg.arc_length_to(tau);
            let speed = seg.deriv(tau).norm();
            if speed < 1e-12 {
                break;
            }
            let delta = (target - current_len) / speed;
            tau = (tau + delta).clamp(0.0, 1.0);
            if delta.abs() < 1e-12 {
                break;
            }
        }

        tau.clamp(0.0, 1.0)
    }

    /// Refine a candidate `(seg_idx, tau)` using Newton's method on the squared-distance function.
    /// Returns the refined `(seg_idx, tau, dist_sq)`.
    fn refine_projection(
        &self,
        seg_idx: usize,
        mut tau: f64,
        p: &Vector2,
    ) -> (usize, f64, f64) {
        let seg = &self.segments[seg_idx];

        for _ in 0..20 {
            let s = seg.eval(tau);
            let sp = seg.deriv(tau);
            let spp = seg.deriv2(tau);
            let diff = s.coords - p;

            // f'(τ) = diff · S'   f''(τ) = ||S'||² + diff · S''
            let f_prime = diff.dot(&sp);
            let f_double = sp.norm_squared() + diff.dot(&spp);

            if f_double.abs() < 1e-12 {
                break;
            }
            let delta = -f_prime / f_double;
            tau = (tau + delta).clamp(0.0, 1.0);
            if delta.abs() < 1e-12 {
                break;
            }
        }

        let dist_sq = (self.segments[seg_idx].eval(tau).coords - p).norm_squared();
        (seg_idx, tau, dist_sq)
    }
}

/// Solve the natural cubic spline second-derivative system via the Thomas algorithm.
///
/// Returns a `Vec<Vector2>` of length `n` with `σ[0] == σ[n-1] == 0`.
fn compute_sigma(points: &[Point2]) -> Vec<Vector2> {
    let n = points.len();
    let m = n - 2; // number of interior knots

    if m == 0 {
        // Two points: linear — all second derivatives are zero.
        return vec![Vector2::zeros(); n];
    }

    // Build RHS: r_i = 6 * (p_{i-1} - 2·p_i + p_{i+1}) for i = 1..n-2.
    let mut rhs_x = vec![0.0_f64; m];
    let mut rhs_y = vec![0.0_f64; m];
    for i in 0..m {
        let r = 6.0
            * (points[i].coords - 2.0 * points[i + 1].coords + points[i + 2].coords);
        rhs_x[i] = r.x;
        rhs_y[i] = r.y;
    }

    // Thomas algorithm on the tridiagonal system: sub=1, main=4, super=1.
    let mut c_prime = vec![0.0_f64; m];
    let mut dx = vec![0.0_f64; m];
    let mut dy = vec![0.0_f64; m];

    c_prime[0] = 1.0 / 4.0;
    dx[0] = rhs_x[0] / 4.0;
    dy[0] = rhs_y[0] / 4.0;

    for i in 1..m {
        let denom = 4.0 - c_prime[i - 1];
        c_prime[i] = if i + 1 < m { 1.0 / denom } else { 0.0 };
        dx[i] = (rhs_x[i] - dx[i - 1]) / denom;
        dy[i] = (rhs_y[i] - dy[i - 1]) / denom;
    }

    let mut sigma_inner = vec![Vector2::zeros(); m];
    sigma_inner[m - 1] = Vector2::new(dx[m - 1], dy[m - 1]);
    for i in (0..m - 1).rev() {
        sigma_inner[i] = Vector2::new(
            dx[i] - c_prime[i] * sigma_inner[i + 1].x,
            dy[i] - c_prime[i] * sigma_inner[i + 1].y,
        );
    }

    // Assemble full vector with zero boundary conditions.
    let mut sigma = vec![Vector2::zeros(); n];
    for i in 0..m {
        sigma[i + 1] = sigma_inner[i];
    }
    sigma
}

impl BoundaryElement for CubicSpline2 {
    fn length(&self) -> f64 {
        self.total_length
    }

    fn at_length(&self, length: f64) -> ManifoldPosition2 {
        let length = length.clamp(0.0, self.total_length);
        let n_segs = self.n_segs();

        // Binary search: find segment whose interval contains `length`.
        let raw = self.cumulative_lengths.partition_point(|&l| l <= length);
        let idx = raw.saturating_sub(1).min(n_segs - 1);

        let seg_offset = length - self.cumulative_lengths[idx];
        let tau = self.tau_at_seg_length(idx, seg_offset);
        self.position_at(idx, tau)
    }

    fn closest_to_point(&self, point: &impl PCoords<2>) -> ManifoldPosition2 {
        let p = point.coords();

        // Coarse grid: 20 uniform samples per segment.
        const N_SAMPLES: usize = 20;

        let mut best_dist_sq = f64::INFINITY;
        let mut best_idx = 0;
        let mut best_tau = 0.0_f64;

        for (seg_idx, seg) in self.segments.iter().enumerate() {
            for k in 0..=N_SAMPLES {
                let tau = k as f64 / N_SAMPLES as f64;
                let dist_sq = (seg.eval(tau).coords - p).norm_squared();
                if dist_sq < best_dist_sq {
                    best_dist_sq = dist_sq;
                    best_idx = seg_idx;
                    best_tau = tau;
                }
            }
        }

        // Newton refinement on the best candidate.
        let (idx, tau, _) = self.refine_projection(best_idx, best_tau, &p);
        self.position_at(idx, tau)
    }

    fn aabb(&self) -> Aabb2 {
        const N: usize = 20;
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for seg in &self.segments {
            for k in 0..=N {
                let tau = k as f64 / N as f64;
                let pt = seg.eval(tau);
                min_x = min_x.min(pt.x);
                min_y = min_y.min(pt.y);
                max_x = max_x.max(pt.x);
                max_y = max_y.max(pt.y);
            }
        }

        Aabb2::new(Point2::new(min_x, min_y), Point2::new(max_x, max_y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn line_points() -> Vec<Point2> {
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 0.0),
        ]
    }

    #[test]
    fn interpolates_control_points() {
        let pts = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 1.0),
        ];
        let spline = CubicSpline2::try_new(&pts).unwrap();

        // The spline must pass through every control point.
        let n = pts.len();
        for (i, p) in pts.iter().enumerate() {
            let t = i as f64 / (n - 1) as f64;
            let pos = spline.at_t(t);
            assert_relative_eq!(pos.point, *p, epsilon = 1e-9);
        }
    }

    #[test]
    fn line_arc_length() {
        // Four collinear points spaced 1 apart → total length must be 3.
        let spline = CubicSpline2::try_new(&line_points()).unwrap();
        assert_relative_eq!(spline.length(), 3.0, epsilon = 1e-6);
    }

    #[test]
    fn at_length_midpoint_on_line() {
        let spline = CubicSpline2::try_new(&line_points()).unwrap();
        let mid = spline.at_length(1.5);
        assert_relative_eq!(mid.point, Point2::new(1.5, 0.0), epsilon = 1e-6);
    }

    #[test]
    fn project_point_onto_line() {
        let spline = CubicSpline2::try_new(&line_points()).unwrap();
        let query = Point2::new(1.5, 1.0);
        let pos = spline.closest_to_point(&query);
        assert_relative_eq!(pos.point, Point2::new(1.5, 0.0), epsilon = 1e-5);
    }

    #[test]
    fn two_point_spline_is_linear() {
        let pts = vec![Point2::new(0.0, 0.0), Point2::new(4.0, 3.0)];
        let spline = CubicSpline2::try_new(&pts).unwrap();
        assert_relative_eq!(spline.length(), 5.0, epsilon = 1e-6);
        let mid = spline.at_t(0.5);
        assert_relative_eq!(mid.point, Point2::new(2.0, 1.5), epsilon = 1e-9);
    }
}
