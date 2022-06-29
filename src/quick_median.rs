pub fn quick_median(list: &mut [f64]) -> f64 {
    if list.is_empty() {
        f64::NAN
    } else {
        *list
            .select_nth_unstable_by(list.len() / 2, |a, b| {
                a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Less)
            })
            .1
    }
}
