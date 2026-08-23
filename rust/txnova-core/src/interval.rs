//! 1-based inclusive intervals.

#[derive(Clone, Debug)]
pub struct Interval<T> {
    pub start: u64,
    pub end: u64,
    pub data: T,
}

#[inline]
pub fn overlaps(a0: u64, a1: u64, b0: u64, b1: u64) -> bool {
    a0 <= b1 && b0 <= a1
}

#[inline]
pub fn gap(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    if overlaps(a0, a1, b0, b1) {
        0
    } else if a1 < b0 {
        b0 - a1 - 1
    } else {
        a0 - b1 - 1
    }
}

/// Sorted by start. `max_end[i]` is the max `end` among `items[0..=i]`,
/// so a query can skip the prefix whose every interval already ended
/// before `qstart`. Long intervals at small starts stay visible.
#[derive(Clone, Debug, Default)]
pub struct IvIndex<T> {
    items: Vec<Interval<T>>,
    max_end: Vec<u64>,
}

impl<T> IvIndex<T> {
    pub fn from_intervals(mut items: Vec<Interval<T>>) -> Self {
        items.sort_by_key(|iv| iv.start);
        let mut max_end = Vec::with_capacity(items.len());
        let mut m = 0u64;
        for iv in &items {
            m = m.max(iv.end);
            max_end.push(m);
        }
        Self { items, max_end }
    }

    pub fn overlapping(&self, start: u64, end: u64) -> impl Iterator<Item = &Interval<T>> {
        let right = self.items.partition_point(|iv| iv.start <= end);
        let left = self.max_end[..right].partition_point(|&m| m < start);
        self.items[left..right]
            .iter()
            .filter(move |iv| iv.end >= start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn idx(ivs: &[(u64, u64, &'static str)]) -> IvIndex<&'static str> {
        IvIndex::from_intervals(
            ivs.iter()
                .map(|&(s, e, d)| Interval {
                    start: s,
                    end: e,
                    data: d,
                })
                .collect(),
        )
    }

    fn hits<'a>(ix: &'a IvIndex<&'static str>, start: u64, end: u64) -> Vec<&'static str> {
        ix.overlapping(start, end).map(|iv| iv.data).collect()
    }

    #[test]
    fn long_interval_is_not_swallowed() {
        let ix = idx(&[
            (1, 1_000_000, "long"),
            (2, 3, "a"),
            (4, 5, "b"),
            (6, 7, "c"),
        ]);
        assert_eq!(hits(&ix, 6, 6), vec!["long", "c"]);
    }

    #[test]
    fn adjacent_closed_intervals_overlap() {
        let ix = idx(&[(10, 20, "left"), (30, 40, "right")]);
        assert!(hits(&ix, 21, 29).is_empty());
        assert_eq!(hits(&ix, 20, 21), vec!["left"]);
        assert_eq!(hits(&ix, 1, 10), vec!["left"]);
        assert_eq!(hits(&ix, 40, 50), vec!["right"]);
    }

    #[test]
    fn empty_index() {
        let ix: IvIndex<u8> = IvIndex::from_intervals(Vec::new());
        assert!(ix.overlapping(1, 10).next().is_none());
    }

    #[test]
    fn prefix_max_end_still_sees_long_span() {
        let mut ivs: Vec<(u64, u64, &'static str)> = vec![(1, 1_000_000, "long")];
        for i in 0..100u64 {
            let s = 10 + i * 10;
            ivs.push((s, s + 2, "s"));
        }
        let ix = idx(&ivs);
        let got = hits(&ix, 500_000, 500_010);
        assert!(got.contains(&"long"));
        assert!(!got.contains(&"s"));
    }
}
