use std::f64::consts::*;
use std::fmt::Debug;
use std::path::Path;

use clap::Parser;

use common::obs_data::RawData;
use common::plot::plotters::style::BLACK;
use common::plot::{draw_plot_svg_with_named_axes, weight_to_rgb, RGBColor};
use common::quick_median::SliceExt;
use common::structs::UnitVec3;
use common::utils::DisplayWithSigmaPercents;
use common::{maths::*, one_dim_search::*};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    in_files: Vec<String>,
}

fn main() {
    let args = Args::parse();
    assert!(!args.in_files.is_empty());

    let mut samples: Vec<_> = args
        .in_files
        .iter()
        .filter_map(parse_file)
        .flatten()
        .collect();
    println!("#samples = {}", samples.len());

    let mut flipped_cnt = 0;
    for s in &mut samples {
        let flipped_reported = (s.reported + PI) % TAU;
        if (flipped_reported - s.actual).abs() < FRAC_PI_2 {
            s.reported = flipped_reported;
            flipped_cnt += 1;
        }
    }
    println!();
    println!(
        "flipped {flipped_cnt} ({:.0}%)",
        flipped_cnt as f64 / samples.len() as f64 * 100.0
    );
    println!();

    let (a, sigma) = lms(&samples);
    let weights = get_weights(&samples, a, sigma);
    let a = ls(&samples, &weights, true);
    let a = gradient(&samples, &weights, a, true);

    let sigma = calc_sigmas(&samples, &weights, a);
    println!("a = {}", DisplayWithSigmaPercents::new(a, sigma));

    plot("plot.svg", &samples, &weights, a, sigma);

    let mut flipped_cnt = 0;
    for s in &mut samples {
        let flipped_reported = (s.reported + PI) % TAU;
        if (flipped_reported - s.actual).abs() < FRAC_PI_2 {
            s.reported = flipped_reported;
            flipped_cnt += 1;
        }
    }
    println!();
    println!(
        "flipped {flipped_cnt} ({:.0}%)",
        flipped_cnt as f64 / samples.len() as f64 * 100.0
    );
    println!();

    let (a, sigma) = lms(&samples);
    let weights = get_weights(&samples, a, sigma);
    let a = ls(&samples, &weights, true);
    let a = gradient(&samples, &weights, a, true);

    let sigma = calc_sigmas(&samples, &weights, a);
    println!("a = {}", DisplayWithSigmaPercents::new(a, sigma));

    plot("plot-flipped.svg", &samples, &weights, a, sigma);
}

fn plot(name: &str, samples: &[Sample], _weights: &[f64], a: f64, sigma: f64) {
    draw_plot_svg_with_named_axes(
        name,
        &samples
            .iter()
            // .zip(weights)
            // .map(|(s, w)| {
            .map(|s| {
                (
                    s.reported.to_degrees(),
                    s.actual.to_degrees(),
                    // weight_to_rgb(*w),
                    BLACK,
                    1.0,
                )
            })
            .collect::<Vec<_>>(),
        &[],
        &[],
        &[
            // (&|x| f(x.to_radians(), a).to_degrees(), weight_to_rgb(1.0)),
            // (
            //     &|x| f(x.to_radians(), a + sigma * 3.0).to_degrees(),
            //     weight_to_rgb(0.0),
            // ),
            // (
            //     &|x| f(x.to_radians(), a - sigma * 3.0).to_degrees(),
            //     weight_to_rgb(0.0),
            // ),
            (&|x| x, RGBColor(0, 0, 255)),
        ],
        "Reported Descent Angle (degrees)",
        "Actual Descent Angle (degrees)",
    )
    .unwrap();
}

fn f(reported: f64, a: f64) -> f64 {
    if (FRAC_PI_2..(FRAC_PI_2 * 3.0)).contains(&reported) {
        reported - a * f64::sin(reported * 2.0)
    } else {
        reported
    }
}

fn f_inv(actual: f64, a: f64) -> Option<f64> {
    if a.abs() >= 0.5 {
        None
    } else if !(FRAC_PI_2..(FRAC_PI_2 * 3.0)).contains(&actual) {
        Some(actual)
    } else {
        let get_next_reported = |x, y| y + a * f64::sin(x * 2.0);
        let mut reported = actual;
        while (f(reported, a) - actual).abs() >= 1e-7 {
            reported = get_next_reported(reported, actual);
        }
        Some(reported)
    }
}

fn lms(points: &[Sample]) -> (f64, f64) {
    fn eval_a(points: &[Sample], a: f64) -> f64 {
        points
            .iter()
            .map(|s| angle_diff(f_inv(s.actual, a).unwrap(), s.reported))
            .map(|div| div * div)
            .collect::<Vec<_>>()
            .median()
    }

    let best_a = one_dim_search(|a| eval_a(points, a), -0.49, 0.49);
    let best_err = eval_a(points, best_a);
    let sigma = best_err.sqrt() * 1.483 * (1. + 5. / (points.len() - 6) as f64);

    println!("[lms] a = {best_a} err = {best_err} sigma = {sigma}");

    (best_a, sigma)
}

fn ls(points: &[Sample], weights: &[f64], print: bool) -> f64 {
    fn eval_a(points: &[Sample], weights: &[f64], a: f64) -> f64 {
        points
            .iter()
            .map(|s| angle_diff(f_inv(s.actual, a).unwrap(), s.reported))
            .zip(weights)
            .map(|(div, &w)| div * div * w)
            .sum::<f64>()
            / weights.iter().sum::<f64>()
    }

    let best_a = one_dim_search(|a| eval_a(points, weights, a), -0.49, 0.49);
    let best_err = eval_a(points, weights, best_a);

    if print {
        println!("[ls] a = {best_a} err = {best_err}");
    }

    best_a
}

fn gradient(points: &[Sample], weights: &[f64], a: f64, print: bool) -> f64 {
    fn eval_a(points: &[Sample], weights: &[f64], a: f64) -> f64 {
        if a.abs() >= 0.5 {
            f64::INFINITY
        } else {
            points
                .iter()
                .map(|s| angle_diff(f_inv(s.actual, a).unwrap(), s.reported))
                .zip(weights)
                .map(|(div, &w)| div * div * w)
                .sum::<f64>()
                / weights.iter().sum::<f64>()
        }
    }

    let a = one_dim_gradint_descent(|a| eval_a(points, weights, a), a, 1e-9);
    let best_err = eval_a(points, weights, a);

    if print {
        println!("[gd] a = {a} err = {best_err}");
    }

    a
}

fn get_weights(points: &[Sample], a: f64, sigma: f64) -> Vec<f64> {
    let get_weight = |err: f64| -> f64 {
        // Smooth transition from 1 to 0
        // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
        const O: f64 = 1.0;
        const F: f64 = 0.6;
        // const O: f64 = 1.5;
        // const F: f64 = 0.6;
        0.5 * (1.0 - ((err.abs() / sigma - O) / F).tanh())
    };
    points
        .iter()
        .map(|s| get_weight(angle_diff(s.actual, f(s.reported, a))))
        .collect()
}

fn calc_sigmas(points: &[Sample], weights: &[f64], orig_a: f64) -> f64 {
    assert!(!points.is_empty());
    let thread_cnt = std::thread::available_parallelism()
        .unwrap()
        .get()
        .min(points.len());

    const D_REPORTED: f64 = to_radians(2.0);

    std::thread::scope(|s| {
        (0..thread_cnt)
            .map(|thread_i| {
                let mut samples = points.to_vec();
                s.spawn(move || {
                    let mut sum = 0.0;
                    for i in (thread_i..samples.len()).step_by(thread_cnt) {
                        let sample = samples[i];

                        let eval_a = |a: f64, samples: &[Sample]| {
                            if a.abs() >= 0.5 {
                                f64::INFINITY
                            } else {
                                samples
                                    .iter()
                                    .map(|s| angle_diff(f_inv(s.actual, a).unwrap(), s.reported))
                                    .zip(weights)
                                    .map(|(div, &w)| div * div * w)
                                    .sum::<f64>()
                                    / weights.iter().sum::<f64>()
                            }
                        };

                        samples[i].reported += D_REPORTED;
                        let new_a = one_dim_gradint_descent(
                            |a| eval_a(a, &samples),
                            orig_a,
                            to_radians(0.01),
                        );
                        samples[i].reported = sample.reported;

                        let da_dr = (orig_a - new_a) / D_REPORTED;
                        let dr = sample.reported - f_inv(sample.actual, orig_a).unwrap();
                        sum += dr * dr * da_dr * da_dr;
                    }
                    sum
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .sum::<f64>()
            .sqrt()
    })
}

#[derive(Debug, Clone, Copy)]
struct Sample {
    reported: f64,
    actual: f64,
}

fn parse_file(in_file: impl AsRef<Path>) -> Option<impl Iterator<Item = Sample>> {
    let data = RawData::read_from_toml_file(in_file).finalize();
    let traj = data.answer?.traj?;

    Some(data.samples.into_iter().flat_map(move |s| {
        let reported = s.da?;
        let k = UnitVec3::new_normalize(traj.point - s.location);
        // let k = traj.point - (traj.direction * 35_000.0) - s.location;
        let actual = descent_angle(s.zenith_dir, k, traj.direction.into_inner()).rem_euclid(TAU);
        Some(Sample { reported, actual })
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use common::rand::random;

    #[test]
    fn f_inv_does_not_stall() {
        for _ in 0..100_000 {
            let a = loop {
                let a = random::<f64>() - 0.5;
                if a != -0.5 {
                    break a;
                };
            };
            let y = random::<f64>() * TAU;
            let x = f_inv(y, a);
            assert!((f(x.unwrap(), a) - y).abs() < 1e-4,);
        }
    }
}
