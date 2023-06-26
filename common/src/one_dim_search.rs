pub fn one_dim_search<F>(mut eval: F, mut min: f64, mut max: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    // On each iteration: sample 200 points and half the search nage.
    // 15 iterations will shrink the search range by a factor of 2^15 ~ 32_000.
    const SAMPLES_PER_ITERATION: usize = 200;
    const ITERATIONS: usize = 15;

    let mut best_err = f64::INFINITY;
    let mut best_point = (max + min) / 2.0;

    for _ in 0..ITERATIONS {
        let search_range = max - min;
        let step = search_range / SAMPLES_PER_ITERATION as f64;

        let mut x = min;
        while x <= max {
            let err = eval(x);
            if err < best_err {
                best_err = err;
                best_point = x;
            }
            x += step;
        }

        let mut new_min = best_point - search_range / 2.0;
        let mut new_max = best_point + search_range / 2.0;

        if new_min < min {
            new_max += min - new_min;
            new_min = min;
        } else if new_max > max {
            new_min -= new_max - max;
            new_max = max;
        }

        min = new_min;
        max = new_max;
    }

    best_point
}

pub fn one_dim_gradint_descent<F>(eval: F, start: f64, d: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let diff = |err0: f64, x: f64| (eval(x + d) - err0) / d;

    let mut k = 1e10;
    let mut x = start;
    let mut err = eval(x);

    for _ in 0..10_000 {
        let diff = diff(err, x);
        let old_x = x;
        let old_err = err;

        x -= diff * k;
        err = eval(x);

        if err >= old_err {
            k /= 2.0;
            x = old_x;
            err = old_err;
        }

        if k < 1e-100 {
            break;
        }
    }

    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_quadratic() {
        let eval = |x: f64| x * (x - 1.0); // min in 0.5
        let min = one_dim_search(eval, 0.0, 1.0);
        assert!((min - 0.5).abs() <= 1e-6);
    }
}
