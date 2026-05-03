//! Experimental tools to work with merging intervals over a domain

use crate::common::Interval;

/// Represents a continuous domain storing non-overlapping intervals.
pub struct IntervalMergeDomain {
    /// The items in the domain. These must be kept sorted.
    items: Vec<Interval>,
}

/// A pending modification to an interval in the domain.
enum ModAction {
    Remove,
    Replace(Interval),
}

/// Allows queuing a modification to an interval in the domain.
/// Calling `set` with a zero-length interval is equivalent to `remove`.
pub trait IntervalEdit {
    fn remove(&mut self);
    fn set(&mut self, interval: Interval);
}

fn set_action(action: &mut Option<ModAction>, interval: Interval) {
    if interval.length() == 0.0 {
        *action = Some(ModAction::Remove);
    } else {
        *action = Some(ModAction::Replace(interval));
    }
}

fn apply_queued(domain: &mut IntervalMergeDomain, index: usize, action: ModAction) {
    domain.items.remove(index);
    if let ModAction::Replace(iv) = action {
        domain.insert(iv);
    }
}

/// A handle to a single interval in the domain. The queued modification is applied when
/// the handle is dropped.
pub struct IntervalHandle<'a> {
    domain: &'a mut IntervalMergeDomain,
    index: usize,
    action: Option<ModAction>,
}

impl<'a> IntervalHandle<'a> {
    pub fn interval(&self) -> &Interval {
        &self.domain.items[self.index]
    }
}

impl<'a> IntervalEdit for IntervalHandle<'a> {
    fn remove(&mut self) {
        self.action = Some(ModAction::Remove);
    }

    fn set(&mut self, interval: Interval) {
        set_action(&mut self.action, interval);
    }
}

impl<'a> Drop for IntervalHandle<'a> {
    fn drop(&mut self) {
        if let Some(action) = self.action.take() {
            apply_queued(self.domain, self.index, action);
        }
    }
}

/// An item yielded by [`IntervalModIter`] that allows queuing a modification.
pub struct IntervalModItem<'a> {
    pub interval: &'a Interval,
    action: &'a mut Option<ModAction>,
}

impl<'a> IntervalEdit for IntervalModItem<'a> {
    fn remove(&mut self) {
        *self.action = Some(ModAction::Remove);
    }

    fn set(&mut self, interval: Interval) {
        set_action(self.action, interval);
    }
}

/// An iterator-like type that yields each interval in the domain as an [`IntervalModItem`],
/// allowing modifications to be queued. All queued modifications are applied when the
/// iterator is dropped.
pub struct IntervalModIter<'a> {
    domain: &'a mut IntervalMergeDomain,
    index: usize,
    actions: Vec<Option<ModAction>>,
}

impl<'a> IntervalModIter<'a> {
    /// Advance to the next interval, returning a handle that allows queuing a modification.
    /// Works like `Iterator::next` but uses a lending lifetime so items may not outlive
    /// the call: use `while let Some(item) = iter.next_item()` rather than `for`.
    pub fn next_item(&mut self) -> Option<IntervalModItem<'_>> {
        if self.index < self.domain.items.len() {
            let i = self.index;
            self.index += 1;
            Some(IntervalModItem {
                interval: &self.domain.items[i],
                action: &mut self.actions[i],
            })
        } else {
            None
        }
    }
}

impl<'a> Drop for IntervalModIter<'a> {
    fn drop(&mut self) {
        // Collect replacement intervals before touching the domain.
        let mut replacements: Vec<Interval> = Vec::new();
        // Remove marked items in reverse order to keep earlier indices valid.
        for i in (0..self.actions.len()).rev() {
            if let Some(action) = self.actions[i].take() {
                self.domain.items.remove(i);
                if let ModAction::Replace(iv) = action {
                    replacements.push(iv);
                }
            }
        }
        // Re-insert replacements; insert() handles merging with remaining items.
        for iv in replacements {
            self.domain.insert(iv);
        }
    }
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

    /// Return a handle to the interval at `index`. The queued modification is applied on drop.
    pub fn modify_at(&mut self, index: usize) -> Option<IntervalHandle<'_>> {
        if index < self.items.len() {
            Some(IntervalHandle { domain: self, index, action: None })
        } else {
            None
        }
    }

    /// Return an iterator-like value over all intervals that allows queuing modifications.
    /// All queued modifications are applied when the returned value is dropped.
    pub fn modify_iter(&mut self) -> IntervalModIter<'_> {
        let actions = self.items.iter().map(|_| None).collect();
        IntervalModIter { domain: self, index: 0, actions }
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

    // modify_at: removing an interval shrinks the domain.
    #[test]
    fn test_modify_at_remove() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        d.insert(interval(4.0, 5.0));
        {
            let mut h = d.modify_at(1).unwrap();
            assert_eq!(h.interval().min, 2.0);
            h.remove();
        }
        assert_eq!(d.items.len(), 2);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[1].min, 4.0);
    }

    // modify_at: no queued action leaves the domain unchanged.
    #[test]
    fn test_modify_at_no_action_is_noop() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        { let _h = d.modify_at(0).unwrap(); } // dropped without queuing anything
        assert_eq!(d.items.len(), 1);
    }

    // modify_at: replacing an interval re-inserts and merges if needed.
    #[test]
    fn test_modify_at_replace_merges() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(5.0, 6.0));
        {
            let mut h = d.modify_at(0).unwrap();
            h.set(interval(0.0, 5.5)); // now overlaps the second interval
        }
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 6.0);
    }

    // modify_at: set with a zero-length interval acts as remove.
    #[test]
    fn test_modify_at_set_zero_length_removes() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        { d.modify_at(0).unwrap().set(interval(0.5, 0.5)); }
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 2.0);
    }

    // modify_at: out-of-bounds index returns None.
    #[test]
    fn test_modify_at_out_of_bounds() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        assert!(d.modify_at(1).is_none());
    }

    // modify_iter: removing items via the bulk iterator.
    #[test]
    fn test_modify_iter_remove_some() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        d.insert(interval(4.0, 5.0));
        {
            let mut it = d.modify_iter();
            while let Some(mut item) = it.next_item() {
                if item.interval.min >= 2.0 {
                    item.remove();
                }
            }
        }
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
    }

    // modify_iter: replacing items via the bulk iterator; invariants hold after.
    #[test]
    fn test_modify_iter_replace_and_invariants() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(3.0, 4.0));
        d.insert(interval(7.0, 8.0));
        {
            let mut it = d.modify_iter();
            while let Some(mut item) = it.next_item() {
                // Expand every interval by 0.5 on each side.
                let iv = item.interval;
                item.set(interval(iv.min - 0.5, iv.max + 0.5));
            }
        }
        // [−0.5, 1.5], [2.5, 4.5], [6.5, 8.5]: no overlaps.
        assert_eq!(d.items.len(), 3);
        assert!(is_sorted(&d));
        assert!(no_overlaps(&d));
    }

    // modify_iter: replacements that now overlap each other are merged.
    #[test]
    fn test_modify_iter_replacements_merge() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        {
            let mut it = d.modify_iter();
            while let Some(mut item) = it.next_item() {
                // Expand both so they overlap: [0, 2.5] and [1.5, 3].
                let iv = item.interval;
                item.set(interval(iv.min, iv.max + 1.5));
            }
        }
        assert_eq!(d.items.len(), 1);
        assert_eq!(d.items[0].min, 0.0);
        assert_eq!(d.items[0].max, 4.5);
    }

    // modify_iter: no queued actions leaves the domain unchanged.
    #[test]
    fn test_modify_iter_no_action_is_noop() {
        let mut d = IntervalMergeDomain::empty();
        d.insert(interval(0.0, 1.0));
        d.insert(interval(2.0, 3.0));
        { let _it = d.modify_iter(); } // iterate nothing, drop immediately
        assert_eq!(d.items.len(), 2);
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
