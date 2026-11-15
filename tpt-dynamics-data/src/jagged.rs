//! Jagged (variable-length) arrays for per-event particle collections.
//!
//! A [`Jagged`] stores `n` events, each holding a variable number of values,
//! as one flat `values` buffer plus `n + 1` offsets. Event `i` is
//! `values[offsets[i]..offsets[i + 1]]`. This is the same layout ROOT and Arrow
//! use for list columns, and it keeps per-event slicing and reductions
//! allocation-free.
//!
//! ```
//! use tpt_dynamics_data::jagged::Jagged;
//!
//! // Three events with 2, 0 and 3 particle pT values (GeV).
//! let pt = Jagged::from_nested(&[vec![45.0, 12.5], vec![], vec![7.0, 3.0, 1.0]]);
//! assert_eq!(pt.len(), 3);
//! assert_eq!(pt.counts(), vec![2, 0, 3]);
//! assert_eq!(pt.event(2), &[7.0, 3.0, 1.0]);
//! assert_eq!(pt.sums(), vec![57.5, 0.0, 11.0]);
//! assert_eq!(pt.maxes(), vec![Some(45.0), None, Some(7.0)]);
//! ```

use core::fmt;

/// Errors from building a [`Jagged`] from raw offsets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum JaggedError {
    /// The offsets slice is empty; it must hold at least the leading `0`.
    NoOffsets,
    /// The first offset is not zero.
    BadFirstOffset(usize),
    /// Offsets decrease at the given index.
    NonMonotone(usize),
    /// The last offset does not equal the number of values.
    LengthMismatch {
        /// Last offset supplied.
        last_offset: usize,
        /// Number of values supplied.
        nvalues: usize,
    },
    /// An event index is out of range.
    OutOfRange {
        /// Requested event index.
        index: usize,
        /// Number of events.
        nevents: usize,
    },
}

impl fmt::Display for JaggedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoOffsets => write!(f, "offsets must contain at least one entry"),
            Self::BadFirstOffset(v) => write!(f, "first offset must be 0, found {v}"),
            Self::NonMonotone(i) => write!(f, "offsets decrease at index {i}"),
            Self::LengthMismatch { last_offset, nvalues } => {
                write!(f, "last offset {last_offset} does not match {nvalues} values")
            }
            Self::OutOfRange { index, nevents } => {
                write!(f, "event {index} out of range for {nevents} events")
            }
        }
    }
}

impl std::error::Error for JaggedError {}

/// Variable-length arrays: `n` events with flat storage and `n + 1` offsets.
#[derive(Clone, Debug, PartialEq)]
pub struct Jagged<T> {
    offsets: Vec<usize>,
    values: Vec<T>,
}

impl<T> Default for Jagged<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Jagged<T> {
    /// An empty jagged array with no events.
    pub fn new() -> Self {
        Self { offsets: vec![0], values: Vec::new() }
    }

    /// Build from raw parts, validating the offsets.
    ///
    /// The offsets must start at 0, be non-decreasing, and end at
    /// `values.len()`.
    pub fn from_parts(offsets: Vec<usize>, values: Vec<T>) -> Result<Self, JaggedError> {
        let Some(&first) = offsets.first() else {
            return Err(JaggedError::NoOffsets);
        };
        if first != 0 {
            return Err(JaggedError::BadFirstOffset(first));
        }
        for (i, pair) in offsets.windows(2).enumerate() {
            if pair[1] < pair[0] {
                return Err(JaggedError::NonMonotone(i + 1));
            }
        }
        let last = offsets[offsets.len() - 1];
        if last != values.len() {
            return Err(JaggedError::LengthMismatch { last_offset: last, nvalues: values.len() });
        }
        Ok(Self { offsets, values })
    }

    /// Number of events.
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no events.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Total number of values across all events.
    pub fn total_values(&self) -> usize {
        self.values.len()
    }

    /// The offsets (`len() + 1` entries, starting at 0).
    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// The flat value buffer.
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// Consume the array and return `(offsets, values)`.
    pub fn into_parts(self) -> (Vec<usize>, Vec<T>) {
        (self.offsets, self.values)
    }

    /// Values of event `i`, or `None` if out of range.
    pub fn get(&self, i: usize) -> Option<&[T]> {
        if i >= self.len() {
            return None;
        }
        Some(&self.values[self.offsets[i]..self.offsets[i + 1]])
    }

    /// Values of event `i`.
    ///
    /// # Panics
    /// Panics if `i >= self.len()`. Use [`Jagged::get`] for a checked lookup.
    pub fn event(&self, i: usize) -> &[T] {
        self.get(i).unwrap_or_else(|| panic!("event {i} out of range for {} events", self.len()))
    }

    /// Number of values in event `i`, or `None` if out of range.
    pub fn count_in(&self, i: usize) -> Option<usize> {
        if i >= self.len() {
            return None;
        }
        Some(self.offsets[i + 1] - self.offsets[i])
    }

    /// Iterate over events as slices.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &[T]> + '_ {
        self.offsets.windows(2).map(move |w| &self.values[w[0]..w[1]])
    }

    /// Append one event built from an iterator of values.
    pub fn push_event<I: IntoIterator<Item = T>>(&mut self, values: I) {
        self.values.extend(values);
        self.offsets.push(self.values.len());
    }

    /// Number of values per event.
    pub fn counts(&self) -> Vec<usize> {
        self.offsets.windows(2).map(|w| w[1] - w[0]).collect()
    }

    /// Reduce each event with `f`, producing one output per event.
    ///
    /// `f` receives the event's values and returns its reduction; this is the
    /// building block for `count`, `sum`, `max` and similar per-event reductions.
    pub fn reduce<U, F: FnMut(&[T]) -> U>(&self, f: F) -> Vec<U> {
        self.iter().map(f).collect()
    }

    /// Keep only the values for which `keep(event_index, value)` is true.
    ///
    /// Events are preserved even when all their values are removed, so the
    /// output has the same number of events as the input.
    #[must_use]
    pub fn filter<F: FnMut(usize, &T) -> bool>(&self, mut keep: F) -> Self
    where
        T: Clone,
    {
        let mut offsets = Vec::with_capacity(self.offsets.len());
        let mut values = Vec::new();
        offsets.push(0);
        for (i, ev) in self.iter().enumerate() {
            values.extend(ev.iter().filter(|v| keep(i, v)).cloned());
            offsets.push(values.len());
        }
        Self { offsets, values }
    }

    /// Keep only the events for which `keep(event_index, event_values)` is true.
    #[must_use]
    pub fn filter_events<F: FnMut(usize, &[T]) -> bool>(&self, mut keep: F) -> Self
    where
        T: Clone,
    {
        let mut out = Self::new();
        for (i, ev) in self.iter().enumerate() {
            if keep(i, ev) {
                out.push_event(ev.iter().cloned());
            }
        }
        out
    }

    /// Apply `f` to every value, keeping the event structure.
    pub fn map<U, F: FnMut(&T) -> U>(&self, f: F) -> Jagged<U> {
        Jagged { offsets: self.offsets.clone(), values: self.values.iter().map(f).collect() }
    }

    /// Flatten all events into one vector, dropping the event structure.
    pub fn flatten(&self) -> Vec<T>
    where
        T: Clone,
    {
        self.values.clone()
    }

    /// Build from nested vectors (one inner `Vec` per event).
    pub fn from_nested(events: &[Vec<T>]) -> Self
    where
        T: Clone,
    {
        let total: usize = events.iter().map(Vec::len).sum();
        let mut offsets = Vec::with_capacity(events.len() + 1);
        let mut values = Vec::with_capacity(total);
        offsets.push(0);
        for ev in events {
            values.extend_from_slice(ev);
            offsets.push(values.len());
        }
        Self { offsets, values }
    }

    /// Convert to nested vectors (one inner `Vec` per event).
    pub fn to_nested(&self) -> Vec<Vec<T>>
    where
        T: Clone,
    {
        self.iter().map(<[T]>::to_vec).collect()
    }
}

impl Jagged<f64> {
    /// Number of values per event as a reduction (alias of [`Jagged::counts`]).
    pub fn sums(&self) -> Vec<f64> {
        self.reduce(|ev| ev.iter().sum())
    }

    /// Largest value per event, or `None` for an empty event. NaN values are
    /// ignored.
    pub fn maxes(&self) -> Vec<Option<f64>> {
        self.reduce(|ev| {
            ev.iter()
                .copied()
                .filter(|v| !v.is_nan())
                .fold(None, |acc: Option<f64>, v| Some(acc.map_or(v, |a| a.max(v))))
        })
    }

    /// Smallest value per event, or `None` for an empty event. NaN values are
    /// ignored.
    pub fn mins(&self) -> Vec<Option<f64>> {
        self.reduce(|ev| {
            ev.iter()
                .copied()
                .filter(|v| !v.is_nan())
                .fold(None, |acc: Option<f64>, v| Some(acc.map_or(v, |a| a.min(v))))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Jagged<f64> {
        Jagged::from_nested(&[vec![45.0, 12.5], vec![], vec![7.0, 3.0, 1.0]])
    }

    #[test]
    fn nested_round_trip() {
        let nested = vec![vec![1.0], vec![], vec![2.0, 3.0]];
        let j = Jagged::from_nested(&nested);
        assert_eq!(j.to_nested(), nested);
        assert_eq!(j.offsets(), &[0, 1, 1, 3]);
        assert_eq!(j.total_values(), 3);
    }

    #[test]
    fn empty_and_default() {
        let j: Jagged<f64> = Jagged::new();
        assert_eq!(j.len(), 0);
        assert_eq!(j.offsets(), &[0]);
        assert_eq!(Jagged::<f64>::from_nested(&[]), j);
        assert_eq!(j.sums().len(), 0);
    }

    #[test]
    fn slicing_and_bounds() {
        let j = sample();
        assert_eq!(j.event(0), &[45.0, 12.5]);
        assert_eq!(j.event(1).len(), 0);
        assert_eq!(j.get(3), None);
        assert_eq!(j.count_in(2), Some(3));
        assert_eq!(j.count_in(9), None);
        let collected: Vec<&[f64]> = j.iter().collect();
        assert_eq!(collected.len(), 3);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn event_panics_out_of_range() {
        let _ = sample().event(3);
    }

    #[test]
    fn reductions() {
        let j = sample();
        assert_eq!(j.counts(), vec![2, 0, 3]);
        assert_eq!(j.sums(), vec![57.5, 0.0, 11.0]);
        assert_eq!(j.maxes(), vec![Some(45.0), None, Some(7.0)]);
        assert_eq!(j.mins(), vec![Some(12.5), None, Some(1.0)]);
        let sizes = j.reduce(<[f64]>::len);
        assert_eq!(sizes, vec![2, 0, 3]);
    }

    #[test]
    fn nan_ignored_by_max() {
        let j = Jagged::from_nested(&[vec![f64::NAN, 2.0]]);
        assert_eq!(j.maxes(), vec![Some(2.0)]);
    }

    #[test]
    fn filter_values_keeps_events() {
        let j = sample();
        let hard = j.filter(|_, &pt| pt > 5.0);
        assert_eq!(hard.len(), 3);
        assert_eq!(hard.to_nested(), vec![vec![45.0, 12.5], vec![], vec![7.0]]);
    }

    #[test]
    fn filter_events_drops_events() {
        let j = sample();
        let two_plus = j.filter_events(|_, ev| ev.len() >= 2);
        assert_eq!(two_plus.to_nested(), vec![vec![45.0, 12.5], vec![7.0, 3.0, 1.0]]);
    }

    #[test]
    fn filter_sees_event_index() {
        let j = sample();
        let only_second = j.filter(|i, _| i == 2);
        assert_eq!(only_second.to_nested(), vec![vec![], vec![], vec![7.0, 3.0, 1.0]]);
    }

    #[test]
    fn map_and_flatten() {
        let j = sample();
        let doubled = j.map(|v| v * 2.0);
        assert_eq!(doubled.event(0), &[90.0, 25.0]);
        assert_eq!(j.flatten(), vec![45.0, 12.5, 7.0, 3.0, 1.0]);
    }

    #[test]
    fn push_event_appends() {
        let mut j: Jagged<i64> = Jagged::new();
        j.push_event([1, 2]);
        j.push_event(core::iter::empty());
        j.push_event(vec![3]);
        assert_eq!(j.offsets(), &[0, 2, 2, 3]);
        assert_eq!(j.event(2), &[3]);
    }

    #[test]
    fn from_parts_validates() {
        assert_eq!(Jagged::<f64>::from_parts(vec![], vec![]), Err(JaggedError::NoOffsets));
        assert_eq!(
            Jagged::from_parts(vec![1, 2], vec![0.0, 0.0]),
            Err(JaggedError::BadFirstOffset(1))
        );
        assert_eq!(
            Jagged::from_parts(vec![0, 2, 1], vec![0.0, 0.0]),
            Err(JaggedError::NonMonotone(2))
        );
        assert_eq!(
            Jagged::from_parts(vec![0, 1], vec![0.0, 0.0]),
            Err(JaggedError::LengthMismatch { last_offset: 1, nvalues: 2 })
        );
        let ok = Jagged::from_parts(vec![0, 0, 2], vec![4.0, 5.0]).unwrap();
        assert_eq!(ok.to_nested(), vec![vec![], vec![4.0, 5.0]]);
    }

    #[test]
    fn into_parts_round_trip() {
        let j = sample();
        let (off, vals) = j.clone().into_parts();
        assert_eq!(Jagged::from_parts(off, vals).unwrap(), j);
    }

    #[test]
    fn error_display() {
        let msg = JaggedError::OutOfRange { index: 4, nevents: 2 }.to_string();
        assert!(msg.contains("event 4"));
    }
}
