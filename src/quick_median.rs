use std::mem::swap;

pub fn quick_median(list: &mut [f64]) -> f64 {
    if list.is_empty() {
        f64::NAN
    } else {
        quick_sort(list, list.len() / 2);
        list[list.len() / 2]
    }
}

fn quick_sort(list: &mut [f64], target_index: usize) {
    match list {
        [] | [_] => (),
        [x, y] => {
            if *x > *y {
                swap(x, y);
            }
        }
        [x, y, z] => {
            if *x > *y {
                swap(x, y);
            }
            if *x > *z {
                swap(x, z);
            }
            if *y > *z {
                swap(y, z);
            }
        }
        [pivot, list @ ..] => {
            let mut min_left = None;
            let mut min_left_val = f64::INFINITY;
            let mut left = 0;
            let mut right = list.len();
            while left < right {
                if list[left] <= *pivot {
                    if list[left] < min_left_val {
                        min_left_val = list[left];
                        min_left = Some(left);
                    }
                    left += 1;
                } else {
                    right -= 1;
                    list.swap(left, right);
                }
            }

            if let Some(min_left) = min_left {
                swap(pivot, &mut list[min_left]);
            }

            if target_index > right {
                quick_sort(&mut list[right..], target_index - (left + 1));
            } else if target_index > 0 {
                quick_sort(&mut list[..left], target_index - 1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn quick_median_test(mut array in prop::collection::vec(-1e9f64..1e9, 1..=1000)) {
            let mut sorted = array.clone();
            sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            let median = quick_median(&mut array);
            prop_assert_eq!(median, sorted[sorted.len()/2]);
        }
    }
}
