//! This module contains an efficient construct for iterating over the indices of a sized window
//! in a domain. For instance, when iterating over a fixed sized window on a manifold, where the
//! locations on the manifold are not evenly spaced.

///
///
/// # Arguments
///
/// * `domain`:
/// * `window_size`:
///
/// returns: DomainWindowIter
///
/// # Examples
///
/// ```
///
/// ```
pub fn iter_domain_window(domain: &[f64], window_size: f64) -> DomainWindowIter {
    DomainWindowIter::new(domain, window_size)
}

/// Represents the indices within a distance window of a central index in a domain. A "domain" in
/// this case is a scalar field with discrete values at each index, such as discrete distances
/// along a manifold or some other field. Each index in the domain has a value, and the values are
/// in order from smallest to largest.
#[derive(Debug, Copy, Clone)]
pub struct DomainWindow {
    /// The index of the central (focus) value in the window, this is the value which the window
    /// is computed around
    pub index: usize,

    /// The index of the first value in the domain which is within the distance window centered on
    /// the focus value.  The largest this value can be is identical to the `index` value.
    pub first: usize,

    /// The index of the last value in the domain which is within the distance window centered on
    /// the focus value.  The smallest this value can be is identical to the `index` value.
    pub last: usize,
}

impl DomainWindow {
    fn new(index: usize, first: usize, last: usize) -> Self {
        Self { index, first, last }
    }

    pub fn len(&self) -> usize {
        self.last - self.first + 1
    }

    pub fn is_empty(&self) -> bool {
        self.first > self.last
    }

    pub fn contains(&self, index: usize) -> bool {
        index >= self.first && index <= self.last
    }

    pub fn iter(&self) -> impl Iterator<Item = usize> {
        self.first..=self.last
    }
}

pub struct DomainWindowIter<'a> {
    domain: &'a [f64],
    half_size: f64,
    window: Option<DomainWindow>,
}

impl<'a> DomainWindowIter<'a> {
    pub fn new(domain: &'a [f64], window_size: f64) -> Self {
        Self {
            domain,
            half_size: (window_size / 2.0).abs(),
            window: None,
        }
    }
}

impl<'a> Iterator for DomainWindowIter<'a> {
    type Item = DomainWindow;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(window) = &self.window {
            // There is already window information in the iterator, so we need to advance it. If
            // the index is the last index in the domain we are done and return None, otherwise
            // we have to move index and then advance the first and last indices

            if window.index < self.domain.len() - 1 {
                // This is is the case where we are not at the end of the domain, so we advance
                // the index and the first and last indices

                let next_index = window.index + 1;
                let next_position = self.domain[next_index];

                let mut first = window.first;
                while self.domain[first] < next_position - self.half_size {
                    first += 1;
                }

                let mut last = window.last;
                while last < self.domain.len() - 1
                    && self.domain[last + 1] <= next_position + self.half_size
                {
                    last += 1;
                }

                // Now we update the window in the iterator and return it
                self.window = Some(DomainWindow::new(next_index, first, last));
                self.window
            } else {
                // This is the case where we are at the end of the domain, so we return None
                None
            }
        } else {
            // This is the case where the iterator has not yet been used and so the stored window
            // is still empty.  We initialize the window by setting the index to zero, and by
            // definition we know that the first index is also zero. The last index is found by
            // searching the domain to find the index of the last value that is less than or equal
            // to the half size.

            let result = self
                .domain
                .binary_search_by(|probe| probe.partial_cmp(&self.half_size).unwrap());
            let end = result.unwrap_or_else(|next_index| next_index - 1);
            self.window = Some(DomainWindow::new(0, 0, end));
            self.window
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_domain() -> Vec<f64> {
        (0..=100).map(|i| i as f64).collect()
    }

    #[test]
    fn test_window() {
        let domain = make_domain();
        let iter = DomainWindowIter::new(&domain, 6.5);
        let half_size = iter.half_size;

        for window in iter {
            let current_position = domain[window.index];
            let indices: Vec<usize> = domain
                .iter()
                .enumerate()
                .filter(|(_i, s)| {
                    **s >= current_position - half_size && **s <= current_position + half_size
                })
                .map(|(i, _s)| i)
                .collect();

            assert_eq!(window.first, indices[0]);
            assert_eq!(window.last, indices[indices.len() - 1]);
        }
    }
}
