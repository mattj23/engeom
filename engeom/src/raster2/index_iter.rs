use crate::image::Pixel;
use crate::raster2::Point2I;
use imageproc::definitions::Image;

pub trait SizeForIndex {
    fn height(&self) -> usize;
    fn width(&self) -> usize;

    fn iter_indices(&self) -> IndexIter<'_, Self>
    where
        Self: Sized,
    {
        IndexIter {
            size: self,
            x: 0,
            y: 0,
        }
    }
}

pub struct IndexIter<'a, T: SizeForIndex + Sized> {
    pub size: &'a T,
    pub x: usize,
    pub y: usize,
}

impl <'a, T: SizeForIndex> IndexIter<'a, T> {
    pub fn new(size: &'a T) -> Self {
        IndexIter {
            size,
            x: 0,
            y: 0,
        }
    }
}

impl<'a, T: SizeForIndex> Iterator for IndexIter<'a, T> {
    type Item = Point2I;

    fn next(&mut self) -> Option<Self::Item> {
        if self.y >= self.size.height() {
            return None;
        }

        let point = Point2I::new(self.x as i32, self.y as i32);
        self.x += 1;

        if self.x >= self.size.width() {
            self.x = 0;
            self.y += 1;
        }

        Some(point)
    }
}

impl<T: Pixel> SizeForIndex for Image<T> {
    fn height(&self) -> usize {
        self.height() as usize
    }

    fn width(&self) -> usize {
        self.width() as usize
    }
}
