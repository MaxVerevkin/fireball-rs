use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::path::Path;
use std::sync::Mutex;

use common::plot::{draw_plot_svg, weight_to_rgb};
use common::quick_median::SliceExt;
use common::{rand, rand_distr};

use rand::prelude::*;
use rand_distr::Normal;

use clap::Parser;

use common::maths::*;
use common::obs_data::Data;

const LS_ITERS: usize = 300_000;
const SIGMA_LS_ITERS: usize = 60_000;
const SIGMA_D_X_RAD: f64 = 0.05;
const SIGMA_SIGMA: f64 = 0.01; // FIXME: do not hard-code sigma

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    in_files: Vec<String>,
}

fn main() {
    let args = Args::parse();
    assert!(!args.in_files.is_empty());

    let mut samples: Vec<_> = args.in_files.iter().flat_map(parse_file).collect();

    println!("LS_ITERS: {LS_ITERS}");
    println!("SIGMA_LS_ITERS: {SIGMA_LS_ITERS}");
    println!("SIGMA_D_X_RAD: {SIGMA_D_X_RAD}");
    println!("SIGMA_SIGMA: {SIGMA_SIGMA}");

    let (a, sigma) = lms(&samples);
    let weights = get_weights(&samples, a, sigma);
    let a = ls(&samples, &weights, a, sigma, LS_ITERS, true);
    let comp_sigma = calc_sigmas(&mut samples, &weights, a);

    println!("comp_sigma = {comp_sigma}");

    draw_plot_svg(
        "plot.svg",
        &samples
            .iter()
            .zip(&weights)
            .map(|(s, w)| (s.reported, s.actual, weight_to_rgb(*w), 1.0))
            .collect::<Vec<_>>(),
        &[],
        &[],
        &[
            (&|x| f(x, a), weight_to_rgb(1.0)),
            (&|x| f(x, a + comp_sigma * 3.0), weight_to_rgb(0.0)),
            (&|x| f(x, a - comp_sigma * 3.0), weight_to_rgb(0.0)),
        ],
    )
    .iter();

    let mut flipped_cnt = 0;
    for s in &mut samples {
        let flipped_reported = (s.reported + PI) % TAU;
        if (flipped_reported - s.actual).abs() < FRAC_PI_2 {
            s.reported = flipped_reported;
            flipped_cnt += 1;
        }
    }
    println!("flipped_cnt = {flipped_cnt}");

    let (a, sigma) = lms(&samples);
    let weights = get_weights(&samples, a, sigma);
    let a = ls(&samples, &weights, a, sigma, LS_ITERS, true);
    let comp_sigma = calc_sigmas(&mut samples, &weights, a);

    println!("comp_sigma = {comp_sigma}");

    draw_plot_svg(
        "plot-flipped.svg",
        &samples
            .iter()
            .zip(&weights)
            .map(|(s, w)| (s.reported, s.actual, weight_to_rgb(*w), 1.0))
            .collect::<Vec<_>>(),
        &[],
        &[],
        &[
            (&|x| f(x, a), weight_to_rgb(1.0)),
            (&|x| f(x, a + comp_sigma * 3.0), weight_to_rgb(0.0)),
            (&|x| f(x, a - comp_sigma * 3.0), weight_to_rgb(0.0)),
        ],
    )
    .iter();

    println!();
}

fn f(x: f64, a: f64) -> f64 {
    if x < FRAC_PI_2 || x > (FRAC_PI_2 * 3.0) {
        x
    } else {
        x - a * f64::sin(x * 2.0)
    }
}

fn f_inv(y: f64, a: f64) -> f64 {
    assert!(a > 0.0);
    if y < FRAC_PI_2 || y > (FRAC_PI_2 * 3.0) {
        y
    } else {
        let get_next_x = |x, y| y + a * f64::sin(x * 2.0);
        let mut x = y;
        while (f(x, a) - y).abs() > 1e-5 {
            x = get_next_x(x, y);
        }
        x
    }
}

fn lms(points: &[Sample]) -> (f64, f64) {
    fn eval_a(points: &[Sample], a: f64) -> f64 {
        points
            .iter()
            .map(|s| s.actual - f(s.reported, a))
            .map(|div| div * div)
            .collect::<Vec<_>>()
            .median()
    }

    let mut best_err = f64::INFINITY;
    let mut best_a = -1.0;
    let mut cur_a = 0.0;

    while cur_a < 1.0 {
        let cur_err = eval_a(&points, cur_a);
        if cur_err < best_err {
            best_a = cur_a;
            best_err = cur_err;
        }
        cur_a += 0.001;
    }

    let sigma = best_err.sqrt() * 1.483 * (1. + 5. / (points.len() - 6) as f64);

    println!("[lms] a = {best_a} err = {best_err} sigma = {sigma}");

    (best_a, sigma)
}

fn ls(
    points: &[Sample],
    weights: &[f64],
    prev_a: f64,
    sigma: f64,
    iters: usize,
    print: bool,
) -> f64 {
    fn eval_a(points: &[Sample], weights: &[f64], a: f64) -> f64 {
        points
            .iter()
            .map(|s| s.actual - f(s.reported, a))
            .zip(weights)
            .map(|(div, &w)| div * div * w)
            .sum::<f64>()
            / weights.iter().sum::<f64>()
    }

    let mut best_err = eval_a(points, weights, prev_a);
    let mut best_a = prev_a;

    let mut rng = thread_rng();
    let distr = Normal::new(prev_a, sigma).unwrap();

    for _ in 0..iters {
        let cur_a = rng.sample(distr);
        let cur_err = eval_a(points, weights, cur_a);
        if cur_err < best_err {
            best_a = cur_a;
            best_err = cur_err;
        }
    }

    if print {
        println!("[ls] a = {best_a} err = {best_err}");
    }

    best_a
}

fn get_weights(points: &[Sample], a: f64, sigma: f64) -> Vec<f64> {
    let get_weight = |err: f64| -> f64 {
        // Smooth transition from 1 to 0
        // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
        const O: f64 = 1.5;
        const F: f64 = 0.6;
        0.5 * (1.0 - ((err.abs() / sigma - O) / F).tanh())
    };
    points
        .iter()
        .map(|s| get_weight(s.actual - f(s.reported, a)))
        .collect()
}

fn calc_sigmas(points: &mut [Sample], weights: &[f64], a: f64) -> f64 {
    let shared_sum = Mutex::new(0f64);

    std::thread::scope(|s| {
        for thread_i in 0..16 {
            let mut samples = points.to_vec();
            let shared_sum = &shared_sum;
            s.spawn(move || {
                let mut sum = 0.0;
                for i in (thread_i..).step_by(16) {
                    let Some(sample) = samples.get(i).copied()
                    else { break };

                    samples[i].reported += SIGMA_D_X_RAD;
                    let new_a = ls(&samples, weights, a, SIGMA_SIGMA, SIGMA_LS_ITERS, false);
                    samples[i].reported = sample.reported;

                    let da_dx = (a - new_a) / SIGMA_D_X_RAD;
                    let dx = sample.reported - f_inv(sample.actual, a);
                    sum += weights[i] * dx * dx * da_dx * da_dx;

                    if i % 100 == 0 {
                        eprintln!("{thread_i}: {i}/{}", samples.len());
                    }
                }
                *shared_sum.lock().unwrap() += sum;
            });
        }
    });

    let sum = *shared_sum.lock().unwrap();
    f64::sqrt(sum)
}

#[derive(Debug, Clone, Copy)]
struct Sample {
    reported: f64,
    actual: f64,
}

fn parse_file(in_file: impl AsRef<Path>) -> impl Iterator<Item = Sample> {
    let data = Data::read_from_toml(in_file);
    let answer = data
        .answer
        .expect("given toml file does not contain answer");
    data.samples.into_iter().flat_map(move |s| {
        let reported = s.da?;
        let k = answer.traj.point - s.location;
        let actual = descent_angle(s.location, k, answer.traj.direction.into_inner());
        Some(Sample { reported, actual })
    })
}
