use std::collections::HashMap;
use crate::geom2::{Boundary2, BoundaryElement, Segment2};
use crate::{Arc2, Point2, Result};


#[derive(Debug, Clone)]
pub(super) enum BData {
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

#[derive(Debug, Clone)]
pub struct BoundaryData2 {
    pub start: Option<Point2>,
    nodes: HashMap<u32, BNode>,
    next_unique_id: u32,
    head_id: u32,
}

impl BoundaryData2 {
    pub fn new_open(start: Point2) -> Self {
        Self {
            start: Some(start),
            nodes: HashMap::new(),
            next_unique_id: 0,
            head_id: u32::MAX,
        }
    }

    pub fn new_closed() -> Self {
        Self {
            start: None,
            nodes: HashMap::new(),
            next_unique_id: 0,
            head_id: u32::MAX,
        }
    }

    pub fn is_closed(&self) -> bool {
        self.start.is_none()
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub(super) fn insert_first(&mut self, data: BData) -> Result<u32> {
        if self.len() > 0 {
            return Err("Boundary already has a first element".into());
        };

        let id = self.next_unique_id;
        self.next_unique_id = self.next_unique_id.wrapping_add(1);

        self.head_id = id;
        if self.is_closed() {
            let node = BNode::new(id, Some(id), Some(id), data);
            self.nodes.insert(id, node);
        } else {
            let node = BNode::new(id, None, None, data);
            self.nodes.insert(id, node);
        }

        Ok(id)
    }

    pub(super) fn insert_after(&mut self, after_this_id: u32, data: BData) -> Result<u32> {
        if !self.nodes.contains_key(&after_this_id) {
            return Err("Node not found".into());
        }

        let old_next_id = self.nodes[&after_this_id].next_id;

        let new_id = self.next_unique_id;
        self.next_unique_id = self.next_unique_id.wrapping_add(1);

        self.nodes.insert(new_id, BNode::new(new_id, old_next_id, Some(after_this_id), data));
        self.nodes.get_mut(&after_this_id).unwrap().next_id = Some(new_id);

        if let Some(next_id) = old_next_id {
            self.nodes.get_mut(&next_id).unwrap().prev_id = Some(new_id);
        }

        Ok(new_id)
    }

    pub(super) fn insert_before(&mut self, before_this_id: u32, data: BData) -> Result<u32> {
        if !self.nodes.contains_key(&before_this_id) {
            return Err("Node not found".into());
        }

        let old_prev_id = self.nodes[&before_this_id].prev_id;

        let new_id = self.next_unique_id;
        self.next_unique_id = self.next_unique_id.wrapping_add(1);

        self.nodes.insert(new_id, BNode::new(new_id, Some(before_this_id), old_prev_id, data));
        self.nodes.get_mut(&before_this_id).unwrap().prev_id = Some(new_id);

        if let Some(prev_id) = old_prev_id {
            self.nodes.get_mut(&prev_id).unwrap().next_id = Some(new_id);
        }

        // In an open boundary, inserting before the head makes the new node the new head
        if !self.is_closed() && before_this_id == self.head_id {
            self.head_id = new_id;
        }

        Ok(new_id)
    }

    pub(super) fn try_remove(&mut self, id: u32) -> Result<()> {
        if !self.nodes.contains_key(&id) {
            return Err("Node not found".into());
        }
        let node = self.nodes.remove(&id).unwrap();

        // Filter out self-references so a single-node closed boundary collapses cleanly
        let eff_prev = node.prev_id.filter(|&p| p != id);
        let eff_next = node.next_id.filter(|&n| n != id);

        if let Some(prev) = eff_prev {
            self.nodes.get_mut(&prev).unwrap().next_id = eff_next;
        }
        if let Some(next) = eff_next {
            self.nodes.get_mut(&next).unwrap().prev_id = eff_prev;
        }

        if id == self.head_id {
            self.head_id = eff_next.unwrap_or(u32::MAX);
        }

        Ok(())
    }

    pub(super) fn tail_id(&self) -> Option<u32> {
        if self.len() == 0 {
            return None;
        }
        if self.is_closed() {
            self.nodes[&self.head_id].prev_id
        } else {
            let mut cur = self.head_id;
            loop {
                match self.nodes[&cur].next_id {
                    Some(next) => cur = next,
                    None => return Some(cur),
                }
            }
        }
    }

    pub(super) fn push_data(&mut self, data: BData) -> Result<()> {
        match self.tail_id() {
            None => { self.insert_first(data)?; }
            Some(tail) => { self.insert_after(tail, data)?; }
        }
        Ok(())
    }

    pub(super) fn last_data(&self) -> Option<&BData> {
        self.tail_id().map(|id| &self.nodes[&id].data)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(x: f64, y: f64) -> BData { BData::Seg((x, y)) }

    /// Walk an open boundary from head to tail, returning node ids in order.
    fn walk_open(bd: &BoundaryData2) -> Vec<u32> {
        if bd.len() == 0 { return vec![]; }
        let mut ids = Vec::new();
        let mut cur = bd.head_id;
        loop {
            ids.push(cur);
            match bd.nodes[&cur].next_id {
                Some(next) => cur = next,
                None => break,
            }
        }
        ids
    }

    /// Walk a closed boundary one full revolution from head, returning node ids in order.
    fn walk_closed(bd: &BoundaryData2) -> Vec<u32> {
        if bd.len() == 0 { return vec![]; }
        let mut ids = Vec::new();
        let start = bd.head_id;
        let mut cur = start;
        loop {
            ids.push(cur);
            let next = bd.nodes[&cur].next_id.expect("closed node must have next");
            if next == start { break; }
            cur = next;
        }
        ids
    }

    /// Verify every next/prev back-link is mutually consistent.
    fn check_integrity(bd: &BoundaryData2) {
        for (id, node) in &bd.nodes {
            if let Some(next) = node.next_id {
                assert_eq!(
                    bd.nodes[&next].prev_id,
                    Some(*id),
                    "node {id} -> next {next}, but {next}.prev != {id}"
                );
            }
            if let Some(prev) = node.prev_id {
                assert_eq!(
                    bd.nodes[&prev].next_id,
                    Some(*id),
                    "node {id} -> prev {prev}, but {prev}.next != {id}"
                );
            }
        }
    }

    #[test]
    fn insert_first_open() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let id = bd.insert_first(seg(1.0, 0.0)).unwrap();
        assert_eq!(bd.head_id, id);
        assert_eq!(bd.len(), 1);
        assert_eq!(bd.nodes[&id].next_id, None);
        assert_eq!(bd.nodes[&id].prev_id, None);
    }

    #[test]
    fn insert_first_closed() {
        let mut bd = BoundaryData2::new_closed();
        let id = bd.insert_first(seg(1.0, 0.0)).unwrap();
        assert_eq!(bd.head_id, id);
        assert_eq!(bd.len(), 1);
        // single-node closed boundary is self-referential
        assert_eq!(bd.nodes[&id].next_id, Some(id));
        assert_eq!(bd.nodes[&id].prev_id, Some(id));
    }

    #[test]
    fn insert_first_twice_is_error() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        bd.insert_first(seg(1.0, 0.0)).unwrap();
        assert!(bd.insert_first(seg(2.0, 0.0)).is_err());
    }

    #[test]
    fn insert_after_open_single() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_open(&bd), vec![a, b]);
        assert_eq!(bd.nodes[&b].next_id, None);
    }

    #[test]
    fn insert_after_open_middle() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let c = bd.insert_after(a, seg(3.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_open(&bd), vec![a, b, c]);
    }

    #[test]
    fn insert_after_closed_single() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_closed(&bd), vec![a, b]);
    }

    #[test]
    fn insert_after_closed_middle() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let c = bd.insert_after(a, seg(3.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_closed(&bd), vec![a, b, c]);
    }

    #[test]
    fn insert_after_invalid_id_is_error() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        assert!(bd.insert_after(99, seg(1.0, 0.0)).is_err());
    }

    #[test]
    fn insert_before_open_head_updates_head() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let b = bd.insert_first(seg(2.0, 0.0)).unwrap();
        let a = bd.insert_before(b, seg(1.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(bd.head_id, a);
        assert_eq!(walk_open(&bd), vec![a, b]);
        assert_eq!(bd.nodes[&a].prev_id, None);
    }

    #[test]
    fn insert_before_open_middle() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let c = bd.insert_after(a, seg(3.0, 0.0)).unwrap();
        let b = bd.insert_before(c, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_open(&bd), vec![a, b, c]);
    }

    #[test]
    fn insert_before_closed_head_does_not_change_head() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let original_head = bd.head_id;
        let z = bd.insert_before(a, seg(0.0, 0.0)).unwrap();
        check_integrity(&bd);
        // head_id must not change for a closed boundary
        assert_eq!(bd.head_id, original_head);
        // z should appear just before a in the walk
        let order = walk_closed(&bd);
        let pos_a = order.iter().position(|&x| x == a).unwrap();
        let pos_z = order.iter().position(|&x| x == z).unwrap();
        assert_eq!((pos_z + 1) % order.len(), pos_a);
    }

    #[test]
    fn insert_before_closed_middle() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let c = bd.insert_after(a, seg(3.0, 0.0)).unwrap();
        let b = bd.insert_before(c, seg(2.0, 0.0)).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_closed(&bd), vec![a, b, c]);
    }

    #[test]
    fn insert_before_invalid_id_is_error() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        assert!(bd.insert_before(99, seg(1.0, 0.0)).is_err());
    }

    #[test]
    fn remove_only_node_open() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        bd.try_remove(a).unwrap();
        assert_eq!(bd.len(), 0);
        assert_eq!(bd.head_id, u32::MAX);
    }

    #[test]
    fn remove_only_node_closed() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        bd.try_remove(a).unwrap();
        assert_eq!(bd.len(), 0);
        assert_eq!(bd.head_id, u32::MAX);
    }

    #[test]
    fn remove_head_open() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        let c = bd.insert_after(b, seg(3.0, 0.0)).unwrap();
        bd.try_remove(a).unwrap();
        check_integrity(&bd);
        assert_eq!(bd.head_id, b);
        assert_eq!(bd.nodes[&b].prev_id, None);
        assert_eq!(walk_open(&bd), vec![b, c]);
    }

    #[test]
    fn remove_tail_open() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        let c = bd.insert_after(b, seg(3.0, 0.0)).unwrap();
        bd.try_remove(c).unwrap();
        check_integrity(&bd);
        assert_eq!(bd.nodes[&b].next_id, None);
        assert_eq!(walk_open(&bd), vec![a, b]);
    }

    #[test]
    fn remove_middle_open() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        let c = bd.insert_after(b, seg(3.0, 0.0)).unwrap();
        bd.try_remove(b).unwrap();
        check_integrity(&bd);
        assert_eq!(walk_open(&bd), vec![a, c]);
    }

    #[test]
    fn remove_head_closed() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        let c = bd.insert_after(b, seg(3.0, 0.0)).unwrap();
        bd.try_remove(a).unwrap();
        check_integrity(&bd);
        assert_eq!(bd.head_id, b);
        assert_eq!(walk_closed(&bd), vec![b, c]);
    }

    #[test]
    fn remove_non_head_closed() {
        let mut bd = BoundaryData2::new_closed();
        let a = bd.insert_first(seg(1.0, 0.0)).unwrap();
        let b = bd.insert_after(a, seg(2.0, 0.0)).unwrap();
        let c = bd.insert_after(b, seg(3.0, 0.0)).unwrap();
        bd.try_remove(b).unwrap();
        check_integrity(&bd);
        assert_eq!(bd.head_id, a);
        assert_eq!(walk_closed(&bd), vec![a, c]);
    }

    #[test]
    fn remove_invalid_id_is_error() {
        let mut bd = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        assert!(bd.try_remove(99).is_err());
    }
}

