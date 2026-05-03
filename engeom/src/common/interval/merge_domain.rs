//! Experimental tools to work with merging intervals over a domain

use crate::common::Interval;

/// Represents a continuous domain storing non-overlapping intervals.
pub struct IntervalMergeDomain {
    /// The items in the domain. These must be kept sorted.
    items: Vec<Interval>,
}

impl<'a> IntoIterator for &'a IntervalMergeDomain {
    type Item = &'a Interval;
    type IntoIter = std::slice::Iter<'a, Interval>;

    fn into_iter(self) -> Self::IntoIter {
        self.items.iter()
    }
}

impl IntervalMergeDomain {
    pub fn empty() -> IntervalMergeDomain {
        IntervalMergeDomain { items: vec![] }
    }

    /// Build a domain from an unsorted, possibly-overlapping collection of intervals.
    /// Intervals are sorted by min, then merged in a single pass.
    pub fn from_intervals(mut intervals: Vec<Interval>) -> IntervalMergeDomain {
        intervals.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        let mut items: Vec<Interval> = Vec::with_capacity(intervals.len());
        for iv in intervals {
            match items.last_mut() {
                Some(last) if last.overlaps(&iv) => *last = Interval::new_contains(last, &iv),
                _ => items.push(iv),
            }
        }
        IntervalMergeDomain { items }
    }

    /// Insert an interval into the merge domain. If it overlaps with one or more existing interval,
    /// they will be replaced with the merged result.
    ///
    /// # Arguments
    ///
    /// * `item`: the new interval to insert
    ///
    /// returns: ()
    pub fn insert(&mut self, item: Interval) {
        let search = self
            .items
            .binary_search_by(|i| i.min.partial_cmp(&item.min).unwrap());

        // The interval either has the same minimum value as the element at check_i, or is just
        // below the minimum value of the interval at check_i.
        let check_i = search.unwrap_or_else(|i| i);

        // Now we determine which intervals it overlaps with.
        let overlaps = self.check_overlaps(&item, check_i.saturating_sub(1));

        if overlaps.is_empty() {
            self.items.insert(check_i, item);
        } else {
            // If there are overlaps, we work in reverse order (to preserve earlier indices),
            // removing intervals and merging them into a working replacement interval.
            let mut working = item.clone();
            for i in overlaps.iter().rev() {
                let popped = self.items.remove(*i);
                working = Interval::new_contains(&working, &popped);
            }

            // Finally, we insert the accumulated interval back into the collection.
            self.items.insert(overlaps[0], working);
        }
    }

    /// Starting at `start_i`, will return a Vec of all indices which overlap with the given
    /// interval. These will inherently be arranged in ascending order.
    fn check_overlaps(&self, item: &Interval, start_i: usize) -> Vec<usize> {
        let mut overlaps = vec![];
        for i in start_i..self.items.len() {
            if self.items[i].min > item.max {
                break;
            }
            if item.overlaps(&self.items[i]) {
                overlaps.push(i);
            }
        }
        overlaps
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Interval;

    fn no_overlaps(x: &IntervalMergeDomain) -> bool {
        for i in 0..x.items.len() {
            for k in i + 1..x.items.len() {
                if x.items[i].overlaps(&x.items[k]) {
                    return false;
                }
            }
        }
        true
    }

    fn is_sorted(x: &IntervalMergeDomain) -> bool {
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
    fn test_insert_into_empty() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(1.0, 3.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 1.0);
        assert_eq!(d.items[0].max, 3.0);
    }

    // Two disjoint intervals inserted in order : both should be kept, sorted.
    #[test]
    fn test_insert_two_disjoint_in_order() {
        let mut d = IntervalMergeDomain::empty();
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
    fn test_insert_two_disjoint_reverse_order() {
        let mut d = IntervalMergeDomain::empty();
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
    fn test_insert_three_disjoint_random_order() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(4.0, 5.0));
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        assert_eq!(d.items.len(), 3);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
    }

    // Two overlapping intervals should merge into one.
    #[test]
    fn test_insert_overlapping_merges() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 2.0));
        d.insert(interval(1.0, 3.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 3.0);
    }

    // Inserting a fully-contained interval leaves the existing one unchanged.
    #[test]
    fn test_insert_contained_interval() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 10.0));
        d.insert(interval(3.0, 7.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 10.0);
    }

    // Inserting an interval that contains an existing one expands outward.
    #[test]
    fn test_insert_containing_interval() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(3.0, 7.0));
        d.insert(interval(0.0, 10.0));
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 10.0);
    }

    // A new interval that bridges two existing ones should merge all three.
    #[test]
    fn test_insert_bridges_two_existing() {
        let mut d = IntervalMergeDomain::empty();
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
    fn test_insert_bridges_three_existing() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        d.insert(interval(4.0, 5.0));
        d.insert(interval(0.5, 4.5)); // spans all three
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 5.0);
    }

    // Iterator yields all intervals in sorted order as shared references.
    #[test]
    fn test_iter_yields_sorted_references() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(4.0, 5.0));
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));

        let collected: Vec<&Interval> = d.into_iter().collect();
        assert_eq!(collected.len(), 3);
        assert_eq!(collected[0].min, 0.0);
        assert_eq!(collected[1].min, 2.0);
        assert_eq!(collected[2].min, 4.0);

        // References point into the domain, domain is still usable after iteration.
        assert_eq!(d.items.len(), 3);
    }

    // from_intervals with an empty vec produces an empty domain.
    #[test]
    fn test_from_intervals_empty() {
        let d = IntervalMergeDomain::from_intervals(vec![]);
        assert_eq!(d.items.len(), 0);
    }

    // from_intervals sorts and keeps disjoint intervals.
    #[test]
    fn test_from_intervals_disjoint_unsorted() {
        let d = IntervalMergeDomain::from_intervals(vec![
            interval(4.0, 5.0),
            interval(0.0, 1.0),
            interval(2.0, 3.0),
        ]);
        assert_eq!(d.items.len(), 3);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[1].min, 2.0);
        assert_eq!(d.items[2].min, 4.0);
    }

    // from_intervals merges overlapping intervals.
    #[test]
    fn test_from_intervals_with_overlaps() {
        let d = IntervalMergeDomain::from_intervals(vec![
            interval(0.0, 2.0),
            interval(3.0, 5.0),
            interval(1.5, 4.0), // bridges the first two
        ]);
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 5.0);
    }

    // from_intervals handles a mix: some overlap, some don't.
    #[test]
    fn test_from_intervals_mixed() {
        let d = IntervalMergeDomain::from_intervals(vec![
            interval(0.0, 1.0),
            interval(0.5, 1.5), // overlaps first
            interval(3.0, 4.0),
            interval(3.5, 5.0), // overlaps previous
            interval(7.0, 8.0), // disjoint
        ]);
        assert_eq!(d.items.len(), 3);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 1.5);
        assert_eq!(d.items[1].min, 3.0);
        assert_eq!(d.items[1].max, 5.0);
        assert_eq!(d.items[2].min, 7.0);
        assert_eq!(d.items[2].max, 8.0);
    }

    // Stress test: many random-ish inserts; invariants must always hold.
    #[test]
    fn test_invariants_hold_after_many_inserts() {
        let mut d = IntervalMergeDomain::empty();
        let starts = [5.0, 0.0, 8.0, 3.0, 6.5, 1.0, 10.0, 2.5];
        let ends =   [6.0, 2.0, 9.0, 4.0, 7.5, 1.5, 11.0, 3.5];
        for (s, e) in starts.iter().zip(ends.iter()) {
            d.insert(interval(*s, *e));
            assert!(is_sorted(&d), "not sorted after inserting [{s}, {e}]");
            assert!(no_overlaps(&d), "overlaps after inserting [{s}, {e}]");
        }
    }
}
