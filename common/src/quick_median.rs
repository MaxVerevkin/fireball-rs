pub trait SliceExt {
    /// Get median of an unsorted slice in O(n) time.
    /// Returns NaN if slice is empty.
    fn median(&mut self) -> f64;
}

pub trait SliceExt2<T> {
    fn median_by_key(&mut self, key_fn: impl Fn(T) -> f64) -> f64;
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

impl<T: Copy> SliceExt2<T> for [T] {
    fn median_by_key(&mut self, key_fn: impl Fn(T) -> f64) -> f64 {
        let cmp = |a: &T, b: &T| f64::total_cmp(&key_fn(*a), &key_fn(*b));
        if self.is_empty() {
            f64::NAN
        } else if self.len() % 2 == 1 {
            let right = *self.select_nth_unstable_by(self.len() / 2, cmp).1;
            let left = *self.select_nth_unstable_by(self.len() / 2 - 1, cmp).1;
            (key_fn(right) + key_fn(left)) * 0.5
        } else {
            key_fn(*self.select_nth_unstable_by(self.len() / 2, cmp).1)
        }
    }
}
