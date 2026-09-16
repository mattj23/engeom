mod binary_mesh;
mod g3d;
mod lptf3;
mod pcd;
mod ply;
pub mod row_scan;
mod stl;
pub mod tol_compress;
pub use tol_compress::curve::*;
pub use tol_compress::mesh::*;
pub use tol_compress::row_points::*;

use crate::{Point3, Result, Vector3};
pub use binary_mesh::*;
pub use g3d::*;
pub use lptf3::*;
pub use pcd::*;
pub use row_scan::{RowAxis, RowScan3, compute_row_scan_mesh, smooth_row_scan, thin_row_scan};
use serde::Serialize;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

#[cfg(feature = "ply")]
pub use ply::*;

#[cfg(feature = "stl")]
pub use stl::*;

pub fn write_xyz(path: &Path, points: &[Point3]) -> Result<()> {
    if path.exists() {
        std::fs::remove_file(path)?;
    }

    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;

    for point in points {
        writeln!(&mut file, "{} {} {}", point.x, point.y, point.z)?;
    }

    Ok(())
}

pub fn write_xyzn(path: &Path, points: &[Point3], normals: &[Vector3]) -> Result<()> {
    if path.exists() {
        std::fs::remove_file(path)?;
    }

    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut buffered = BufWriter::new(file);

    for (point, normal) in points.iter().zip(normals.iter()) {
        writeln!(
            &mut buffered,
            "{} {} {} {} {} {}",
            point.x, point.y, point.z, normal.x, normal.y, normal.z
        )?;
    }

    Ok(())
}

pub fn json_elements_save<T>(path: &Path, item: &T) -> Result<()>
where
    T: Serialize + ?Sized,
{
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let bytes = serde_json::to_vec_pretty(item)?;
    writer.write_all(&bytes)?;
    Ok(())
}
