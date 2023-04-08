pub trait SliceExt {
    /// Get median of an unsorted slice in O(n) time.
    /// Returns NaN if slice is empty.
    fn median(&mut self) -> f64;
}

impl SliceExt for [f64] {
    fn median(&mut self) -> f64 {
        if self.is_empty() {
            f64::NAN
        } else if self.len() % 2 == 1 {
            let right = *self
                .select_nth_unstable_by(self.len() / 2, f64::total_cmp)
                .1;
            let left = *self
                .select_nth_unstable_by(self.len() / 2 - 1, f64::total_cmp)
                .1;
            (right + left) * 0.5
        } else {
            *self
                .select_nth_unstable_by(self.len() / 2, f64::total_cmp)
                .1
        }
    }
}
