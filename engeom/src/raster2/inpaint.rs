use crate::image::{GrayImage, Luma};
use crate::raster2::raster_mask::RasterMask;
use crate::raster2::{Point2I, ScalarImage};
use parry3d_f64::na::DMatrix;

const KNOWN: u8 = 0;
const BAND: u8 = 1;
const UNKNOWN: u8 = 2;

/// This function performs an inpainting operation on a scalar valued image using Alexander Telea's
/// algorithm based on the Fast Marching Method, using a fill mask to determine which pixels should
/// be inpainted and an image mask to determine which pixels are part of the image.  Pixels which
/// are not on the image mask are not used to determine the inpainted value of neighboring pixels.
///
/// # Arguments
///
/// * `values`: The scalar image to inpaint.
/// * `fill_mask`: A mask image where pixels with a value of 255 are identified as pixels to be
///   inpainted, and pixels with a value of 0 are not to be inpainted.
/// * `image_mask`: A mask image where pixels with a value of 255 are part of the image and will
///   contribute to determining the inpainted value of neighboring fill regions, and pixels with a
///   value of 0 are not part of the image and are ignored in the inpainting process.
/// * `radius`: The radius of the neighborhood to consider when inpainting pixels. This determines
///   how far from the fill mask pixels the algorithm will look to find neighboring pixels to
///   compute the inpainted value.
///
/// returns: ImageBuffer<Luma<u16>, Vec<u16, Global>>
pub fn inpaint(
    depth: &ScalarImage<u16>,
    fill_mask: &GrayImage,
    image_mask: &RasterMask,
    radius: usize,
) -> ScalarImage<u16> {
    let mut fill = Fill::new(depth, fill_mask, image_mask, radius);

    while let Some(px) = fill.band.pop() {
        fill.flags[(px.y, px.x)] = KNOWN;

        let neighbors = [
            (px.x as i32 - 1, px.y as i32),
            (px.x as i32 + 1, px.y as i32),
            (px.x as i32, px.y as i32 - 1),
            (px.x as i32, px.y as i32 + 1),
        ];
        for (nx, ny) in neighbors.iter() {
            if *nx < 0 || *ny < 0 || *nx >= depth.width() as i32 || *ny >= depth.height() as i32 {
                continue;
            }

            if fill.flags[(*ny as usize, *nx as usize)] != UNKNOWN {
                continue;
            }

            let ns = [
                solve_eikonal(
                    ny - 1,
                    *nx,
                    *ny,
                    nx - 1,
                    depth.height() as usize,
                    depth.width() as usize,
                    &fill.distances,
                    &fill.flags,
                ),
                solve_eikonal(
                    ny + 1,
                    *nx,
                    *ny,
                    nx + 1,
                    depth.height() as usize,
                    depth.width() as usize,
                    &fill.distances,
                    &fill.flags,
                ),
                solve_eikonal(
                    ny - 1,
                    *nx,
                    *ny,
                    nx + 1,
                    depth.height() as usize,
                    depth.width() as usize,
                    &fill.distances,
                    &fill.flags,
                ),
                solve_eikonal(
                    ny + 1,
                    *nx,
                    *ny,
                    nx - 1,
                    depth.height() as usize,
                    depth.width() as usize,
                    &fill.distances,
                    &fill.flags,
                ),
            ];
            let d = ns.iter().fold(f32::INFINITY, |a, b| a.min(*b));
            fill.distances[(*ny as usize, *nx as usize)] = d;

            // Fill the pixel
            let val = inpaint_pixel(*ny, *nx, &fill);
            fill.values.put_pixel(*nx as u32, *ny as u32, Luma([val]));
            // fill.image_mask.set(*nx as u32, *ny as u32, true);
            fill.image_mask
                .set_point_unchecked(Point2I::new(*nx, *ny), true);

            fill.flags[(*ny as usize, *nx as usize)] = BAND;

            // Sort in reverse order so that the last element is the one with the smallest distance
            fill.band.push(Pixel::new(*nx as usize, *ny as usize, d));
            fill.band
                .sort_by(|a, b| b.distance.partial_cmp(&a.distance).unwrap());
        }
    }

    fill.values
}

#[derive(Copy, Clone, Debug)]
struct Pixel {
    x: usize,
    y: usize,
    distance: f32,
}

impl Pixel {
    fn new(x: usize, y: usize, distance: f32) -> Self {
        Self { x, y, distance }
    }
}

struct Fill {
    values: ScalarImage<u16>,
    distances: DMatrix<f32>,
    flags: DMatrix<u8>,
    band: Vec<Pixel>,
    radius: usize,
    image_mask: RasterMask,
}

impl Fill {
    fn new(
        values: &ScalarImage<u16>,
        fill_mask: &GrayImage,
        image_mask: &RasterMask,
        radius: usize,
    ) -> Self {
        let mut distances = DMatrix::zeros(values.height() as usize, values.width() as usize);

        // Fill the initial distances with the infinity distance value
        distances.fill(f32::INFINITY);

        // Turn the mask pixels into UNKNOWN flags
        let mut flags = DMatrix::zeros(values.height() as usize, values.width() as usize);
        for (x, y, pixel) in fill_mask.enumerate_pixels() {
            if pixel[0] > 0 {
                flags[(y as usize, x as usize)] = UNKNOWN;
            }
        }

        // Create the narrow band of pixels on the border of the mask
        let mut band = Vec::new();
        for (x, y, pixel) in fill_mask.enumerate_pixels() {
            if pixel[0] == 0 {
                continue;
            }

            let neighbors = [
                (x as i32 - 1, y as i32),
                (x as i32 + 1, y as i32),
                (x as i32, y as i32 - 1),
                (x as i32, y as i32 + 1),
            ];
            for (nx, ny) in neighbors.iter() {
                if *nx < 0
                    || *ny < 0
                    || *nx >= values.width() as i32
                    || *ny >= values.height() as i32
                {
                    continue;
                }

                if fill_mask.get_pixel(*nx as u32, *ny as u32)[0] == 0 {
                    distances[(*ny as usize, *nx as usize)] = 0.0;
                    flags[(*ny as usize, *nx as usize)] = BAND;
                    band.push(Pixel::new(*nx as usize, *ny as usize, 0.0));
                }
            }
        }

        // Compute the distances between the initial mask contour and the pixels outside the mask
        // using the Fast Marching Method
        compute_outside_distances(
            values.height() as usize,
            values.width() as usize,
            &mut distances,
            &flags,
            &band,
            radius,
        );

        Self {
            values: values.clone(),
            distances,
            flags,
            band,
            radius,
            image_mask: image_mask.clone(),
        }
    }
}

fn compute_outside_distances(
    h: usize,
    w: usize,
    distances: &mut DMatrix<f32>,
    flags: &DMatrix<u8>,
    band: &[Pixel],
    radius: usize,
) {
    let mut work_band = band.to_vec();
    let mut work_flags = flags.clone();

    // Swap the values in work_flags such that KNOWN becomes UNKNOWN and vice versa
    for i in 0..work_flags.nrows() {
        for j in 0..work_flags.ncols() {
            if flags[(i, j)] == KNOWN {
                work_flags[(i, j)] = UNKNOWN;
            } else if flags[(i, j)] == UNKNOWN {
                work_flags[(i, j)] = KNOWN;
            }
        }
    }

    let mut last_distance = 0.0;
    while !work_band.is_empty() {
        if last_distance > (radius as f32) * 2.0 {
            break;
        }

        let px = work_band.pop().unwrap();
        work_flags[(px.y, px.x)] = KNOWN;
        let neighbors = [
            (px.x as i32 - 1, px.y as i32),
            (px.x as i32 + 1, px.y as i32),
            (px.x as i32, px.y as i32 - 1),
            (px.x as i32, px.y as i32 + 1),
        ];
        for (nx, ny) in neighbors.iter() {
            if *nx < 0 || *ny < 0 || *nx >= w as i32 || *ny >= h as i32 {
                continue;
            }

            if work_flags[(*ny as usize, *nx as usize)] != UNKNOWN {
                continue;
            }

            let ns = [
                solve_eikonal(ny - 1, *nx, *ny, nx - 1, h, w, distances, flags),
                solve_eikonal(ny + 1, *nx, *ny, nx + 1, h, w, distances, flags),
                solve_eikonal(ny - 1, *nx, *ny, nx + 1, h, w, distances, flags),
                solve_eikonal(ny + 1, *nx, *ny, nx - 1, h, w, distances, flags),
            ];
            last_distance = ns.iter().fold(f32::INFINITY, |a, b| a.min(*b));
            distances[(*ny as usize, *nx as usize)] = last_distance;
            work_flags[(*ny as usize, *nx as usize)] = BAND;

            // Sort in reverse order so that the last element is the one with the smallest distance
            work_band.push(Pixel::new(*nx as usize, *ny as usize, last_distance));
            work_band.sort_by(|a, b| b.distance.partial_cmp(&a.distance).unwrap());
        }
    }

    // Distances are opposite of the actual propagation distance, so invert them
    for i in 0..distances.nrows() {
        for j in 0..distances.ncols() {
            distances[(i, j)] *= -1.0;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_eikonal(
    y1: i32,
    x1: i32,
    y2: i32,
    x2: i32,
    h: usize,
    w: usize,
    distances: &DMatrix<f32>,
    flags: &DMatrix<u8>,
) -> f32 {
    if y1 < 0 || y1 >= h as i32 || x1 < 0 || x1 >= w as i32 {
        return f32::INFINITY;
    }

    if y2 < 0 || y2 >= h as i32 || x2 < 0 || x2 >= w as i32 {
        return f32::INFINITY;
    }

    let f1 = flags[(y1 as usize, x1 as usize)];
    let f2 = flags[(y2 as usize, x2 as usize)];

    if f1 == KNOWN && f2 == KNOWN {
        let d1 = distances[(y1 as usize, x1 as usize)];
        let d2 = distances[(y2 as usize, x2 as usize)];
        let d = 2.0 - (d1 - d2) * (d1 - d2);
        if d > 0.0 {
            let r = d.sqrt();
            let s = (d1 + d2 - r) / 2.0;
            if s >= d1 && s >= d2 {
                return s;
            }
            let s = s + r;
            if s >= d1 && s >= d2 {
                return s;
            }
            // Unsolvable
            return f32::INFINITY;
        }
    }

    if f1 == KNOWN {
        let d1 = distances[(y1 as usize, x1 as usize)];
        return 1.0 + d1;
    }

    if f2 == KNOWN {
        let d2 = distances[(y2 as usize, x2 as usize)];
        return 1.0 + d2;
    }

    f32::INFINITY
}

fn inpaint_pixel(ny: i32, nx: i32, fill: &Fill) -> u16 {
    let d = fill.distances[(ny as usize, nx as usize)];
    let (dgy, dgx) = pixel_gradient(ny as u32, nx as u32, fill);
    let mut sum: f32 = 0.0;
    let mut weight_sum: f32 = 0.0;

    for nby in (ny - fill.radius as i32)..(ny + fill.radius as i32 + 1) {
        if nby < 0 || nby >= fill.values.height() as i32 {
            continue;
        }

        for nbx in (nx - fill.radius as i32)..(nx + fill.radius as i32 + 1) {
            if nbx < 0
                || nbx >= fill.values.width() as i32
                // || !fill.image_mask.get(nbx as u32, nby as u32)
                || !fill.image_mask.get_point(Point2I::new(nbx, nby))
            {
                continue;
            }

            if fill.flags[(nby as usize, nbx as usize)] == UNKNOWN {
                continue;
            }

            let dir_y = (ny as f32) - (nby as f32);
            let dir_x = (nx as f32) - (nbx as f32);
            let dir_len_sq = dir_y * dir_y + dir_x * dir_x;
            let dir_len = dir_len_sq.sqrt();
            if dir_len > fill.radius as f32 {
                continue;
            }

            // Compute weight
            let mut dir_factor = (dir_y * dgy + dir_x * dgx).abs();
            if dir_factor == 0.0 {
                dir_factor = 1e-6;
            }

            let nb_dist = fill.distances[(nby as usize, nbx as usize)];
            let level_factor = 1.0 / (1.0 + (nb_dist - d).abs());
            let dist_factor = 1.0 / (dir_len * dir_len_sq);

            let weight = (dir_factor * level_factor * dist_factor).abs();
            sum += weight * (fill.values.get_pixel(nbx as u32, nby as u32)[0] as f32);
            weight_sum += weight;
        }
    }

    (sum / weight_sum) as u16
}

fn pixel_gradient(y: u32, x: u32, fill: &Fill) -> (f32, f32) {
    let v = fill.values.get_pixel(x, y)[0] as f32;

    let prev_y = (y as i32) - 1;
    let next_y = (y as i32) + 1;

    let grad_y = if prev_y < 0 || next_y >= fill.values.height() as i32 {
        f32::INFINITY
    } else {
        let flag_prev_y = fill.flags[(prev_y as usize, x as usize)];
        let flag_next_y = fill.flags[(next_y as usize, x as usize)];
        if flag_prev_y != UNKNOWN && flag_next_y != UNKNOWN {
            (fill.values.get_pixel(x, next_y as u32)[0] as f32
                - fill.values.get_pixel(x, prev_y as u32)[0] as f32)
                / 2.0
        } else if flag_prev_y != UNKNOWN {
            v - fill.values.get_pixel(x, prev_y as u32)[0] as f32
        } else if flag_next_y != UNKNOWN {
            fill.values.get_pixel(x, next_y as u32)[0] as f32 - v
        } else {
            0.0
        }
    };

    let prev_x = (x as i32) - 1;
    let next_x = (x as i32) + 1;
    let grad_x = if prev_x < 0 || next_x >= fill.values.width() as i32 {
        f32::INFINITY
    } else {
        let flag_prev_x = fill.flags[(y as usize, prev_x as usize)];
        let flag_next_x = fill.flags[(y as usize, next_x as usize)];
        if flag_prev_x != UNKNOWN && flag_next_x != UNKNOWN {
            (fill.values.get_pixel(next_x as u32, y)[0] as f32
                - fill.values.get_pixel(prev_x as u32, y)[0] as f32)
                / 2.0
        } else if flag_prev_x != UNKNOWN {
            v - fill.values.get_pixel(prev_x as u32, y)[0] as f32
        } else if flag_next_x != UNKNOWN {
            fill.values.get_pixel(next_x as u32, y)[0] as f32 - v
        } else {
            0.0
        }
    };

    (grad_y, grad_x)
}
