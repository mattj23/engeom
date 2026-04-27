use std::collections::HashMap;
use crate::geom2::{Boundary2, BoundaryElement, Segment2};
use crate::{Arc2, Point2};


#[derive(Debug, Clone)]
enum BData {
    /// A line segment, containing the end point
    Seg((f64, f64)),

    /// An arc, containing the center point, end point, and whether the arc is clockwise
    Arc((f64, f64, f64, f64, bool)),
}

#[derive(Debug, Clone)]
struct BNode {
    id: u32,
    next_id: Option<u32>,
    prev_id: Option<u32>,
    data: BData
}

impl BNode {
    pub fn new(id: u32, next_id: Option<u32>, prev_id: Option<u32>, data: BData) -> BNode {
        Self {
            id, next_id, prev_id, data
        }
    }

}

pub struct BCursor<'a> {
    data: &'a mut BoundaryData2,
    node_id: u32,
}

#[derive(Debug, Clone)]
pub struct BoundaryData2 {
    pub start: Option<Point2>,
    elements: HashMap<u32, BNode>,
}

impl BoundaryData2 {
    pub fn new_open(start: Point2) -> Self {
        Self {
            start: Some(start),
            elements: HashMap::new(),
        }
    }

    pub fn new_closed() -> Self {
        Self {
            start: None,
            elements: HashMap::new(),
        }
    }

    pub fn is_closed(&self) -> bool {
        self.start.is_none()
    }

    pub fn try_to_boundary(&self) -> crate::Result<Boundary2> {
        // let mut elements: Vec<Box<dyn BoundaryElement>> = Vec::new();
        // let mut last_point = self.start;
        //
        // for e in self.elements.iter() {
        //     match e {
        //         BData::Seg((x, y)) => {
        //             let end = Point2::new(*x, *y);
        //             let seg = Segment2::try_new(&last_point, &end)?;
        //             last_point = end;
        //             elements.push(Box::new(seg));
        //         }
        //         BData::Arc((cx, cy, ex, ey, cw)) => {
        //             let end = Point2::new(*ex, *ey);
        //             let center = Point2::new(*cx, *cy);
        //             let arc = Arc2::try_new_ends(&last_point, &end, &center, *cw)?;
        //             last_point = end;
        //             elements.push(Box::new(arc));
        //         }
        //     }
        // }
        //
        // Ok(Boundary2::new(elements))
        todo!()
    }
}

