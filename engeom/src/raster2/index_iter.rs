use crate::image::Pixel;
use crate::na::DMatrix;
use crate::raster2::{Point2I, ScalarRaster};
use imageproc::definitions::Image;

pub trait SizeForIndex {
    fn iter_indices(&self) -> IndexIter;
}

pub struct IndexIter {
    pub x: usize,
    pub y: usize,
    pub width: usize,
    pub height: usize,
}

impl IndexIter {
    pub fn new(width: usize, height: usize) -> Self {
        IndexIter {
            x: 0,
            y: 0,
            width,
            height,
        }
    }
}

impl Iterator for IndexIter {
    type Item = Point2I;

    fn next(&mut self) -> Option<Self::Item> {
        if self.y >= self.height {
            return None;
        }

        let point = Point2I::new(self.x as i32, self.y as i32);
        self.x += 1;

        if self.x >= self.width {
            self.x = 0;
            self.y += 1;
        }

        Some(point)
    }
}

impl<T: Pixel> SizeForIndex for Image<T> {
    fn iter_indices(&self) -> IndexIter {
        IndexIter::new(self.width() as usize, self.height() as usize)
    }
}

impl<T> SizeForIndex for DMatrix<T> {
    fn iter_indices(&self) -> IndexIter {
        IndexIter::new(self.ncols(), self.nrows())
    }
}

impl SizeForIndex for ScalarRaster {
    fn iter_indices(&self) -> IndexIter {
        IndexIter::new(self.width() as usize, self.height() as usize)
    }
}
