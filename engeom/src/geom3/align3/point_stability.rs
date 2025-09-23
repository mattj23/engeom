//! This module is for testing the alignment stability contributions of points

use crate::common::points::mean_point;
use crate::geom3::align3::RcParams3;
use crate::geom3::align3::jacobian::{copy_jacobian, point_plane_jacobian};
use crate::geom3::mesh::MeshSurfPoint;
use crate::na::{Dyn, Matrix, Owned, U1, U6, Vector};
use crate::{Iso3, Result, SurfacePoint3};
use itertools::all;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::utils::hashset::HashSet;

#[derive(Debug, Clone)]
pub struct StabilityResult {
    pub tx: Vec<f64>,
    pub ty: Vec<f64>,
    pub tz: Vec<f64>,
    pub rx: Vec<f64>,
    pub ry: Vec<f64>,
    pub rz: Vec<f64>,
}

impl StabilityResult {
    pub fn summary(&self) -> [f64; 6] {
        [
            sq_resid(&self.tx),
            sq_resid(&self.ty),
            sq_resid(&self.tz),
            sq_resid(&self.rx),
            sq_resid(&self.ry),
            sq_resid(&self.rz),
        ]
    }
}

fn sq_resid(res: &[f64]) -> f64 {
    res.iter().map(|r| r * r).sum()
}

pub fn point_stability_reduce(
    fraction: f64,
    points: &[MeshSurfPoint],
    dt: f64,
    dr: Option<f64>,
) -> Result<(StabilityResult, Vec<MeshSurfPoint>)> {
    let dr = if let Some(v) = dr {
        v
    } else {
        dr_from_dt(points, dt)
    };

    let baseline = point_stability(points, dt, Some(dr))?;
    let reference = baseline
        .summary()
        .into_iter()
        .map(|v| v * fraction)
        .collect::<Vec<_>>();

    let ordered = {
        let mut tx_i = baseline.tx.iter().enumerate().collect::<Vec<_>>();
        let mut ty_i = baseline.ty.iter().enumerate().collect::<Vec<_>>();
        let mut tz_i = baseline.tz.iter().enumerate().collect::<Vec<_>>();
        let mut rx_i = baseline.rx.iter().enumerate().collect::<Vec<_>>();
        let mut ry_i = baseline.ry.iter().enumerate().collect::<Vec<_>>();
        let mut rz_i = baseline.rz.iter().enumerate().collect::<Vec<_>>();
        tx_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        ty_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        tz_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        rx_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        ry_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        rz_i.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        vec![tx_i, ty_i, tz_i, rx_i, ry_i, rz_i]
    };

    let mut pop_order = Vec::new();
    let mut passed = vec![HashSet::new(); 6];

    for i in 0..points.len() {
        let mut added = HashSet::new();
        for (j, p) in passed.iter_mut().take(6).enumerate() {
            p.insert(ordered[j][i].0);
            added.insert(ordered[j][i].0);
        }
        let mut to_remove = Vec::new();
        for idx in added.iter() {
            if all(&passed, |s| s.contains(idx)) {
                pop_order.push(*idx);
                to_remove.push(*idx);
            }
        }
        for idx in to_remove.iter() {
            for p in passed.iter_mut() {
                p.remove(idx);
            }
            // for j in 0..6 {
            //     passed[j].remove(idx);
            // }
        }
    }

    let mut lower = 0;
    let mut upper = points.len();
    let mut last_good = baseline.clone();

    let mut to_index = (upper - lower) / 2 + lower;
    loop {
        let skip = pop_order.iter().take(to_index).collect::<HashSet<_>>();
        let reduced = points
            .iter()
            .enumerate()
            .filter(|(i, _)| !skip.contains(i))
            .map(|(_, p)| *p)
            .collect::<Vec<_>>();
        let reduced_stability = point_stability(&reduced, dt, Some(dr))?;
        let reduced_summary = reduced_stability.summary();
        let is_ok = all(0..6, |i| reduced_summary[i] >= reference[i]);
        if is_ok {
            lower = to_index;
            last_good = reduced_stability;
        } else {
            upper = to_index;
        }
        let next_index = (upper - lower) / 2 + lower;
        if next_index == to_index {
            return Ok((last_good, reduced));
        }
        to_index = next_index;
    }
}

fn dr_from_dt(points: &[MeshSurfPoint], dt: f64) -> f64 {
    let mean = mean_point(points);

    let mut mx = 0.0f64;
    let mut my = 0.0f64;
    let mut mz = 0.0f64;
    for p in points {
        mx = mx.max((p.point().x - mean.x).abs());
        my = my.max((p.point().y - mean.y).abs());
        mz = mz.max((p.point().z - mean.z).abs());
    }
    let m = mx.max(my).max(mz);
    dt / m
}

pub fn point_stability(
    points: &[MeshSurfPoint],
    dt: f64,
    dr: Option<f64>,
) -> Result<StabilityResult> {
    let dr = if let Some(v) = dr {
        v
    } else {
        dr_from_dt(points, dt)
    };

    Ok(StabilityResult {
        tx: sub_problem(points, 0, dt)?,
        ty: sub_problem(points, 1, dt)?,
        tz: sub_problem(points, 2, dt)?,
        rx: sub_problem(points, 3, dr)?,
        ry: sub_problem(points, 4, dr)?,
        rz: sub_problem(points, 5, dr)?,
    })
}

fn sub_problem(points: &[MeshSurfPoint], fixed_dim: usize, fixed_value: f64) -> Result<Vec<f64>> {
    let points = points.iter().map(|p| p.sp).collect::<Vec<_>>();
    let problem = PointStability::new(&points, fixed_dim, fixed_value);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    if report.termination.was_successful() {
        Ok(result.residuals().unwrap().as_slice().to_vec())
    } else {
        Err(format!(
            "Failed to solve sub-problem for fixed_dim {} fixed_value {}, \
        problem may be poorly conditioned",
            fixed_dim, fixed_value
        )
        .into())
    }
}

struct PointStability {
    fixed_dim: usize,
    fixed_value: f64,
    points: Vec<SurfacePoint3>,
    moved: Vec<SurfacePoint3>,
    params: RcParams3,
}

impl PointStability {
    fn new(points: &[SurfacePoint3], fixed_dim: usize, fixed_value: f64) -> Self {
        let mean = mean_point(points);
        let params = RcParams3::from_initial(&Iso3::identity(), &mean);

        let mut item = Self {
            fixed_dim,
            fixed_value,
            points: points.to_vec(),
            moved: points.to_vec(),
            params,
        };

        item.move_points();
        item
    }

    fn move_points(&mut self) {
        let mut copied = self.params.clone();
        copied.set_index(self.fixed_dim, self.fixed_value);

        for (i, p) in self.points.iter().enumerate() {
            self.moved[i] = copied.transform() * p;
        }
    }
}

impl LeastSquaresProblem<f64, Dyn, U6> for PointStability {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params.set(x);
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params.x
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.points.iter()).enumerate() {
            res[i] = c.scalar_projection(&p.point);
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let _center = self.params.transform() * self.params.rc;
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.points.iter()).enumerate() {
            let mut values = point_plane_jacobian(&p.point, c, &self.params);
            values[self.fixed_dim] = 0.0;
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}
