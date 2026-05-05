use super::IntervalOps;

pub struct MergeDomain<T: IntervalOps> {
    pub items: Vec<T>,
}

impl<T: IntervalOps> MergeDomain<T> {
    pub fn from_intervals(mut intervals: Vec<T>) -> Self {
        intervals.sort_by(|a, b| a.min().partial_cmp(&b.min()).unwrap());
        let mut items: Vec<T> = Vec::with_capacity(intervals.len());
        for x in intervals {
            match items.last_mut() {
                Some(last) if last.overlaps(x) => {
                    *last = x.new_containing(last);
                }
                _ => items.push(x),
            }
        }

        let mut constructed = MergeDomain { items };
        constructed.mut_resolve_wrap();
        constructed
    }

    pub fn new_empty() -> Self {
        MergeDomain { items: Vec::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn insert(&mut self, item: T) {
        let search = self
            .items
            .binary_search_by(|i| i.min().partial_cmp(&item.min()).unwrap());

        // The interval either has the same minimum value as the element at check_i, or is just
        // below the minimum value of the interval at check_i.
        let check_i = search.unwrap_or_else(|i| i);

        // Now we determine which intervals it overlaps with.
        let overlaps = self.find_overlaps_no_wrap(&item, check_i.saturating_sub(1));

        if overlaps.is_empty() {
            self.items.insert(check_i, item);
        } else {
            // If there are overlaps, we work in reverse order (to preserve earlier indices),
            // removing intervals and merging them into a working replacement interval.
            let mut working = item.clone();
            for i in overlaps.iter().rev() {
                let popped = self.items.remove(*i);
                working = working.new_containing(&popped);
            }

            // Finally, we insert the accumulated interval back into the collection.
            self.items.insert(overlaps[0], working);
        }

        self.mut_resolve_wrap();
    }

    fn mut_resolve_wrap(&mut self) {
        if T::wraps() && self.items.len() > 2 {
            // If the domain wraps, it is possible that the last and the first element may overlap.
            // There cannot be additional overlaps, because the item after the first element is
            // guaranteed to not overlap with the first element, and the item before the last
            // element is guaranteed to not overlap with the last element. Therefore, we only need
            // to check the first and last.
            //
            // If there is an overlap, the easiest thing to do is to pop the last element off of
            // the list and replace the first.
            if self.items[0].overlaps(self.items[self.items.len() - 1]) {
                let last = self.items.pop().unwrap();
                self.items[0] = last.new_containing(&self.items[0]);
            }
        }
    }

    /// Starting at `start_i`, will return a Vec of all indices which overlap with the given
    /// interval. These will inherently be arranged in ascending order.
    fn find_overlaps_no_wrap(&self, item: &T, start_i: usize) -> Vec<usize> {
        let mut overlaps = vec![];
        for i in start_i..self.items.len() {
            if self.items[i].min() > item.max() {
                break;
            }
            if item.overlaps(self.items[i]) {
                overlaps.push(i);
            }
        }
        overlaps
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::interval::{Interval, IntervalMerge};

    fn no_overlaps(x: &IntervalMerge) -> bool {
        for i in 0..x.items.len() {
            for k in i + 1..x.items.len() {
                if x.items[i].overlaps(x.items[k]) {
                    return false;
                }
            }
        }
        true
    }

    fn is_sorted(x: &IntervalMerge) -> bool {
        for i in 0..(x.items.len().saturating_sub(1)) {
            if x.items[i + 1].min < x.items[i].min {
                return false;
            }
        }
        true
    }

    fn interval(min: f64, max: f64) -> Interval {
        Interval::new(min, max)
    }

    // Insert a single interval into an empty domain.
    #[test]
    fn insert_into_empty() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(1.0, 3.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 1.0);
        assert_eq!(d.items[0].max, 3.0);
    }

    // Two disjoint intervals inserted in order : both should be kept, sorted.
    #[test]
    fn insert_two_disjoint_in_order() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        assert_eq!(d.items.len(), 2);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[1].min, 2.0);
    }

    // Two disjoint intervals inserted in reverse order : result must be sorted.
    #[test]
    fn insert_two_disjoint_reverse_order() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(2.0, 3.0));
        d.insert(interval(0.0, 1.0));
        assert_eq!(d.items.len(), 2);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[1].min, 2.0);
    }

    // Three disjoint intervals inserted in random order.
    #[test]
    fn insert_three_disjoint_random_order() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(4.0, 5.0));
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        assert_eq!(d.items.len(), 3);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
    }

    // Two overlapping intervals should merge into one.
    #[test]
    fn insert_overlapping_merges() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(0.0, 2.0));
        d.insert(interval(1.0, 3.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 3.0);
    }

    // Inserting a fully-contained interval leaves the existing one unchanged.
    #[test]
    fn insert_contained_interval() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(0.0, 10.0));
        d.insert(interval(3.0, 7.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 10.0);
    }

    // Inserting an interval that contains an existing one expands outward.
    #[test]
    fn insert_containing_interval() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(3.0, 7.0));
        d.insert(interval(0.0, 10.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 10.0);
    }

    // A new interval that bridges two existing ones should merge all three.
    #[test]
    fn insert_bridges_two_existing() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(3.0, 4.0));
        d.insert(interval(0.5, 3.5)); // overlaps both
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 4.0);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
    }

    // A new interval that bridges three existing ones should merge all four.
    #[test]
    fn insert_bridges_three_existing() {
        let mut d = IntervalMerge::new_empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        d.insert(interval(4.0, 5.0));
        d.insert(interval(0.5, 4.5)); // spans all three
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 5.0);
    }
}
