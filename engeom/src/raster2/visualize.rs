//! Quick plotting of raster and matrix data

use crate::Result;
use crate::image::{Rgba, RgbaImage};
use crate::na::DMatrix;
use crate::raster2::{Point2IIndexAccess, SizeForIndex, d_matrix_min_max};
use colorgrad::Gradient;
use std::path::Path;

pub fn render_d_matrix(
    path: &Path,
    matrix: &DMatrix<f64>,
    gradient: &dyn Gradient,
    limits: Option<(f64, f64)>,
) -> Result<()> {
    let (min_z, max_z) = if let Some((min, max)) = limits {
        (min, max)
    } else {
        d_matrix_min_max(matrix)
    };

    let mut image = RgbaImage::new(matrix.ncols() as u32, matrix.nrows() as u32);
    for p in matrix.iter_indices() {
        let v = matrix.get_at(p).ok_or("Failed to get value at point")?;

        if v.is_finite() {
            let f = (v - min_z) / (max_z - min_z);
            let color = gradient.at(f as f32).to_rgba8();
            image.put_pixel(p.x as u32, p.y as u32, Rgba(color));
        }
    }

    image.save(path).map_err(|e| e.into())
}
