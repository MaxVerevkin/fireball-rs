use std::f64::consts::*;
use std::fmt::Debug;
use std::path::Path;

use clap::Parser;

use common::obs_data::RawData;
use common::plot::plotters::style::BLACK;
use common::plot::{draw_plot_svg_with_named_axes, weight_to_rgb, RGBColor};
use common::quick_median::SliceExt;
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

    let samples: Vec<_> = args
        .in_files
        .iter()
        .filter_map(parse_file)
        .flatten()
        .collect();

    draw_plot_svg_with_named_axes(
        "plot-clear.svg",
        &samples
            .iter()
            .map(|s| {
                (
                    s.reported.to_degrees(),
                    s.actual.to_degrees(),
                    RGBColor(0, 0, 0),
                    1.0,
                )
            })
            .collect::<Vec<_>>(),
        &[],
        &[],
        &[(&|x| x, RGBColor(0, 0, 255))],
        "Reported Azimuth (degrees)",
        "Actual Azimuth (degrees)",
    )
    .unwrap();

    let (a, sigma) = lms(&samples);
    let weights = get_weights(&samples, a, sigma);
    let a = ls(&samples, &weights, true);
    let a = gradient(&samples, &weights, a, true);

    let sigma = calc_sigmas(&samples, &weights, a);
    println!("a = {}", DisplayWithSigmaPercents::new(a, sigma));

    plot("plot.svg", &samples, &weights, a, sigma);
}

fn plot(name: &str, samples: &[Sample], weights: &[f64], a: f64, sigma: f64) {
    draw_plot_svg_with_named_axes(
        name,
        &samples
            .iter()
            // .zip(weights)
            // .map(|(s, w)| {
            .map(|s| (s.reported.to_degrees(), s.actual.to_degrees(), BLACK, 1.0))
            .collect::<Vec<_>>(),
        &[],
        &[],
        &[
            (&|x| f(x.to_radians(), a).to_degrees(), weight_to_rgb(1.0)),
            (
                &|x| f(x.to_radians(), a + sigma * 3.0).to_degrees(),
                weight_to_rgb(0.0),
            ),
            (
                &|x| f(x.to_radians(), a - sigma * 3.0).to_degrees(),
                weight_to_rgb(0.0),
            ),
            (&|x| x, RGBColor(0, 0, 255)),
        ],
        "Reported Azimuth (degrees)",
        "Actual Azimuth (degrees)",
    )
    .unwrap();
}

fn f(reported: f64, a: f64) -> f64 {
    assert!(reported >= 0.0);
    assert!(reported < TAU);
    reported + a * reported.sin()
}

fn f_inv(actual: f64, a: f64) -> Option<f64> {
    assert!(actual >= 0.0);
    assert!(actual < TAU);
    if a.abs() >= 1.0 {
        None
    } else {
        let get_next_reported =
            |actual: f64, reported: f64| (actual - a * reported.sin()).rem_euclid(TAU);
        let mut reported = actual;
        while (f(reported, a) - actual).abs() >= 1e-7 {
            reported = get_next_reported(actual, reported);
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

    let best_a = one_dim_search(|a| eval_a(points, a), -0.99, 0.99);
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

    let best_a = one_dim_search(|a| eval_a(points, weights, a), -0.99, 0.99);
    let best_err = eval_a(points, weights, best_a);

    if print {
        println!("[ls] a = {best_a} err = {best_err}");
    }

    best_a
}

fn gradient(points: &[Sample], weights: &[f64], a: f64, print: bool) -> f64 {
    fn eval_a(points: &[Sample], weights: &[f64], a: f64) -> f64 {
        if a.abs() >= 1.0 {
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

    const D_REPORTED: f64 = radians(0.1);

    std::thread::scope(|s| {
        (0..thread_cnt)
            .map(|thread_i| {
                let mut samples = points.to_vec();
                s.spawn(move || {
                    let mut sum = 0.0;
                    for i in (thread_i..samples.len()).step_by(thread_cnt) {
                        let sample = samples[i];

                        let eval_a = |a: f64, samples: &[Sample]| {
                            if a.abs() >= 1.0 {
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
                        let new_a =
                            one_dim_gradint_descent(|a| eval_a(a, &samples), orig_a, radians(0.01));
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
    let point = data.answer?.point?;

    Some(data.samples.into_iter().flat_map(move |s| {
        let reported = s.z0?.rem_euclid(TAU);
        let actual = s.calc_azimuth(point);
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
