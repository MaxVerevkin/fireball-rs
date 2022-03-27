pub fn draw_hitogram(values: &[f64], lines: usize) {
    let (min, max) = min_max(values);
    let k = f64::sqrt(values.len() as f64).ceil() as usize;
    //let k = f64::log2(values.len() as f64).ceil() as usize + 1;
    let step = (max - min) / k as f64;
    let mut buckets = vec![0usize; k];

    for &x in values {
        let rel = x - min;
        let bucket = (rel / step) as usize;
        if bucket == k {
            buckets[k - 1] += 1;
        } else {
            buckets[bucket] += 1;
        }
    }

    let max_h = *buckets.iter().max().unwrap();

    let top_label = max_h.to_string();
    let bottom_label = format!("{}0", " ".repeat(top_label.len() - 1));
    let left_padding = " ".repeat(top_label.len());

    for line in 0..lines {
        let sub = max_h as f64 * (lines - line - 1) as f64 / lines as f64;

        let scaled: Vec<f64> = buckets
            .iter()
            .map(|&x| (x as f64 - sub) / (max_h / lines) as f64)
            .collect();

        let label = match line {
            0 => &top_label,
            x if x == lines - 1 => &bottom_label,
            _ => &left_padding,
        };
        println!("{}{}", label, draw_line(&scaled));
    }

    let left_label = format!("{:.1}", min);
    let mid_label = format!("{:.1}", (min + max) / 2.);
    let righ_label = format!("{:.1}", max);

    let padding = k
        .saturating_sub(left_label.len())
        .saturating_sub(mid_label.len())
        .saturating_sub(righ_label.len()) as f64
        / 2.;

    println!(
        "{}{}{}{}{}{}",
        left_padding,
        left_label,
        " ".repeat(padding.floor() as usize),
        mid_label,
        " ".repeat(padding.ceil() as usize),
        righ_label
    );
}

fn draw_line(values: &[f64]) -> String {
    values
        .iter()
        .map(|x| match (x * 9.).clamp(0., 9.) as u8 {
            0 => ' ',
            1 => '\u{2581}',
            2 => '\u{2582}',
            3 => '\u{2583}',
            4 => '\u{2584}',
            5 => '\u{2585}',
            6 => '\u{2586}',
            7 => '\u{2587}',
            _ => '\u{2588}',
        })
        .collect()
}

fn min_max(values: &[f64]) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = -f64::INFINITY;
    for &x in values {
        if x < min {
            min = x;
        }
        if x > max {
            max = x;
        }
    }
    (min, max)
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
