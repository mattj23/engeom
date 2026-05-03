//! Experimental tools to work with merging intervals over a domain

use crate::common::Interval;

/// Represents a continuous domain storing non-overlapping intervals.
pub struct IntervalMergeDomain {
    /// The items in the domain. These must be kept sorted.
    items: Vec<Interval>,
}

impl IntervalMergeDomain {
    pub fn empty() -> IntervalMergeDomain {
        IntervalMergeDomain { items: vec![] }
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
            // If there are no overlaps, we simply insert into the items collection
            self.items.push(item);
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
        let mut i = start_i;
        let mut overlaps = vec![];
        while item.overlaps(&self.items[i]) {
            overlaps.push(i);
            i += 1;
        }
        overlaps
    }
}
