//! This module has an abstraction to help manage optimization parameters for multiple relative
//! entities

use crate::geom3::align3::{RcParams3, T3Storage};
use crate::geom3::{Iso3, Point3};
use parry3d_f64::na::{DVector, Dyn, Matrix, Owned};

type Jacobian = Matrix<f64, Dyn, Dyn, Owned<f64, Dyn, Dyn>>;

fn enforce_initial(count: usize, initial: Option<&[Iso3]>) -> Vec<Iso3> {
    if let Some(initial) = initial {
        assert_eq!(initial.len(), count);
        initial.to_vec()
    } else {
        vec![Iso3::identity(); count]
    }
}

pub struct ParamHandler {
    pub static_i: usize,
    pub params: Vec<RcParams3>,
    raw_params: DVector<f64>,
    count: usize,
}

impl ParamHandler {
    pub fn new(static_i: usize, mean_points: Vec<Point3>, initial: Option<&[Iso3]>) -> Self {
        let count = mean_points.len();
        let initial = enforce_initial(count, initial);
        let raw_params = DVector::zeros((count - 1) * 6);
        assert_eq!(count, initial.len());

        let params = initial
            .iter()
            .zip(mean_points.iter())
            .map(|(t, p)| RcParams3::from_initial(t, p))
            .collect();

        let mut item = Self {
            static_i,
            params,
            raw_params,
            count,
        };
        item.compute();
        item
    }

    pub fn params(&self) -> &DVector<f64> {
        &self.raw_params
    }

    pub fn p_index(&self, cloud_i: usize) -> usize {
        if cloud_i > self.static_i {
            cloud_i - 1
        } else {
            cloud_i
        }
    }

    pub fn get_transform(&self, cloud_i: usize) -> Iso3 {
        *self.params[cloud_i].transform()
    }

    pub fn set_param(&mut self, x: &DVector<f64>) {
        self.raw_params.copy_from(x);
        self.compute();
    }

    fn compute(&mut self) {
        for i in 0..self.count {
            if i != self.static_i {
                let p_index = self.p_index(i);
                let param = self.raw_params.fixed_rows::<6>(p_index * 6);
                self.params[i].set(&param.into_owned());
            }
        }
    }

    pub fn relative_transform(&self, test_i: usize, ref_i: usize) -> Iso3 {
        self.params[ref_i].inverse() * self.params[test_i].transform()
    }

    pub fn set_jacobian(
        &self,
        matrix: &mut Jacobian,
        row: usize,
        cloud_index: usize,
        values: &T3Storage,
    ) {
        if cloud_index != self.static_i {
            let start_col = self.p_index(cloud_index) * 6;
            matrix[(row, start_col)] = values[0];
            matrix[(row, start_col + 1)] = values[1];
            matrix[(row, start_col + 2)] = values[2];
            matrix[(row, start_col + 3)] = values[3];
            matrix[(row, start_col + 4)] = values[4];
            matrix[(row, start_col + 5)] = values[5];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geom3::Vector3;

    fn gen_mean_pts() -> Vec<Point3> {
        vec![
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(4.0, 5.0, 6.0),
            Point3::new(7.0, 8.0, 9.0),
        ]
    }

    fn gen_transform(x: f64, y: f64, z: f64, rx: f64, ry: f64, rz: f64) -> Iso3 {
        Iso3::new(Vector3::new(x, y, z), Vector3::new(rx, ry, rz))
    }

    // #[test]
    // fn test_param_round_trip() {
    //     let mut handler = ParamHandler::new(1, gen_mean_pts());
    //
    //     let t0 = gen_transform(1.0, 2.0, 3.0, 0.5, 0.6, 0.7);
    //     let t1 = gen_transform(4.0, 5.0, 6.0, 0.8, 0.9, 1.0);
    //
    //     handler.initialize_transform(0, &t0);
    //     handler.initialize_transform(2, &t1);
    //
    //     let t0_ = handler.get_transform(0);
    //     let t1_ = handler.get_transform(2);
    //     assert_relative_eq!(t0, t0_, epsilon = 1e-10);
    //     assert_relative_eq!(t1, t1_, epsilon = 1e-10);
    // }
    //
    // #[test]
    // fn test_relative_transform() {
    //     let mut handler = ParamHandler::new(1, gen_mean_pts());
    //     handler.initialize_transform(0, &gen_transform(1.0, 2.0, 3.0, 0.5, 0.6, 0.7));
    //     handler.initialize_transform(2, &gen_transform(4.0, 5.0, 6.0, 0.8, 0.9, 1.0));
    //
    //     let p0 = Point3::new(10.0, 11.0, 12.0);
    //     let p1 = Point3::new(-5.0, -4.0, -3.0);
    //
    //     let expected = dist(
    //         &(handler.get_transform(0) * p0),
    //         &(handler.get_transform(2) * p1),
    //     );
    //     let test = dist(&(handler.relative_transform(0, 2) * p0), &p1);
    //
    //     assert_relative_eq!(expected, test, epsilon = 1e-10);
    // }
}
