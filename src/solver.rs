//! The implementation

use std::collections::HashMap;
use std::f64::consts::*;
use std::ops;
use std::time::Instant;

use common::histogram::draw_hitogram;
use common::maths::*;
use common::obs_data::{Data, DataSample, PartDiff};
use common::plot::{draw_plot_svg, plotters, weight_to_rgb};
use common::quick_median::SliceExt;
use common::rand::random;
use common::structs::*;
use common::{constants::*, Sigmas};
use common::{rand, rand_distr};

use image::{Rgb, RgbImage};

use plotters::style::colors;
use rand::prelude::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::StandardNormal;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

pub mod pair;
use pair::PairTrajectory;

const ERR_INCR_COLOR: yansi::Color = yansi::Color::Red;
const ERR_DECR_COLOR: yansi::Color = yansi::Color::Unset;
fn eval_err_colored(err: f64) -> yansi::Paint<f64> {
    yansi::Paint::new(err).fg(if err > f64::EPSILON * 16.0 {
        ERR_INCR_COLOR
    } else {
        ERR_DECR_COLOR
    })
}

const GD_FIRST_PARAMS: GradDescentParams = GradDescentParams {
    target_err_rel: 1e-5,
    debug: false,
    realtime_graph: true,
};

const GD_FINAL_PARAMS: GradDescentParams = GradDescentParams {
    target_err_rel: 1e-10,
    debug: false,
    realtime_graph: true,
};

const GD_SIGMA_PARAMS: GradDescentParams = GradDescentParams {
    target_err_rel: 1e-10,
    debug: false,
    realtime_graph: false,
};

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
    altitude_error_multiplier: f64,
}

/// Parameters to tweak the algorithm
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub no_da_flip: bool,
    pub no_altitudes: bool,
    pub no_azimuths: bool,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub flash: Geodetic,
    pub velocity: Vec3,
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver {
            data,
            params,
            altitude_error_multiplier: 1.0,
        }
    }

    /// The weighted sum of pairwise plane crossings
    pub fn pairwise(&self) -> (Line, Sigmas) {
        let mut point = Vec3::new(0.0, 0.0, 0.0);
        let mut dir = Vec3::new(0.0, 0.0, 0.0);
        let mut sum = 0.0;
        let mut results = HashMap::new();
        for (i, s1) in self.data.iter().enumerate() {
            for (j, s2) in self.data.iter().enumerate().skip(i + 1) {
                if let Some(pair) = PairTrajectory::calculate(s1, s2) {
                    point += pair.line.point * pair.weight;
                    dir += pair.line.direction * pair.weight;
                    sum += pair.weight;
                    results.insert((i, j), pair);
                }
            }
        }
        let ans = Line {
            point: point / sum,
            direction: UnitVec3::new_normalize(dir),
        };

        let mut sigmas = Sigmas::default();

        let mut push = |si: usize, err: f64, pd: PartDiff| {
            let (p_diff, d_diff) = self.pairwise_diff(si, &pd, &results, point, dir, sum);
            assert_approx_eq!(ans.direction.dot(d_diff), 0.0);
            sigmas.add_point(p_diff, err);
            sigmas.v_angle += d_diff.norm_squared() * err * err;
        };

        for (sample_i, sample) in self.data.iter().enumerate() {
            let plane_norm = (ans.point - sample.location)
                .cross(*ans.direction)
                .normalize();
            let semispace_k = ans.direction.cross(ans.point - sample.location);

            let project_k = |k: UnitVec3| -> Vec3 {
                let proj = if k.dot(semispace_k) >= 0.0 {
                    (k.into_inner() - plane_norm * k.dot(plane_norm)).normalize()
                } else {
                    ans.direction * k.dot(*ans.direction).signum()
                };
                assert_approx_eq!(proj.dot(plane_norm), 0.0);
                proj
            };

            if let Some(k_end) = sample.k_end {
                // Project `k_end` on the apparent arc as defined by `ans`.
                let k_end_projected = project_k(k_end);

                // Calculate z0 and h0 errors
                let projected_az =
                    Azimuthal::from(sample.geo_location.geocentric_to_local(k_end_projected));
                let h0_error = angle_diff(sample.h0.unwrap(), projected_az.h);

                // HACK: using this approximation because of zenith problems.
                // let z0_error = angle_diff(sample.z0.unwrap(), projected_az.z);
                let z0_error = (k_end.dot(k_end_projected).clamp(-1.0, 1.0).acos().powi(2)
                    - h0_error * h0_error)
                    .max(0.0)
                    .sqrt();

                push(
                    sample_i,
                    z0_error,
                    PartDiff {
                        dz0: 1.0,
                        dh0: 0.0,
                        dzb: 0.0,
                        dhb: 0.0,
                        dda: 0.0,
                    },
                );

                push(
                    sample_i,
                    h0_error,
                    PartDiff {
                        dz0: 0.0,
                        dh0: 1.0,
                        dzb: 0.0,
                        dhb: 0.0,
                        dda: 0.0,
                    },
                );
            }
            if let Some(k_start) = sample.k_start {
                // Project `k_start` on the apparent arc as defined by `ans`.
                let k_start_projected = project_k(k_start);

                // Calculate zb and hb errors
                let projected_az =
                    Azimuthal::from(sample.geo_location.geocentric_to_local(k_start_projected));
                let hb_error = angle_diff(sample.hb.unwrap(), projected_az.h);

                // HACK: using this approximation because of zenith problems.
                // let zb_error = angle_diff(sample.zb.unwrap(), projected_az.z);
                let zb_error = (k_start
                    .dot(k_start_projected)
                    .clamp(-1.0, 1.0)
                    .acos()
                    .powi(2)
                    - hb_error * hb_error)
                    .max(0.0)
                    .sqrt();

                push(
                    sample_i,
                    zb_error,
                    PartDiff {
                        dz0: 0.0,
                        dh0: 0.0,
                        dzb: 1.0,
                        dhb: 0.0,
                        dda: 0.0,
                    },
                );

                push(
                    sample_i,
                    hb_error,
                    PartDiff {
                        dz0: 0.0,
                        dh0: 0.0,
                        dzb: 0.0,
                        dhb: 1.0,
                        dda: 0.0,
                    },
                );
            }
        }

        (ans, sigmas.sqrt())
    }

    /// Returns (point_diff, direction_diff)
    pub fn pairwise_diff(
        &self,
        diff_i: usize,
        pd: &PartDiff,
        results: &HashMap<(usize, usize), PairTrajectory>,
        point: Vec3,
        dir: Vec3,
        sum: f64,
    ) -> (Vec3, Vec3) {
        let mut point_diff = Vec3::new(0.0, 0.0, 0.0);
        let mut dir_diff = Vec3::new(0.0, 0.0, 0.0);
        let mut sum_diff = 0.0;
        // let mut stat = Vec::new();
        for i in 0..self.data.len() {
            if i != diff_i {
                if let Some(pair) = results.get(&(i.min(diff_i), i.max(diff_i))) {
                    let (p_diff, d_diff, w_diff) =
                        PairTrajectory::diff(&self.data.samples[diff_i], &self.data.samples[i], pd)
                            .unwrap();
                    point_diff += pair.line.point * w_diff + p_diff * pair.weight;
                    dir_diff += pair.line.direction * w_diff + d_diff * pair.weight;
                    sum_diff += w_diff;
                    // if diff_i == 13 && pd.dzb == 1.0 {
                    //     stat.push((
                    //         (pair.line.point * w_diff + p_diff * pair.weight).norm(),
                    //         pair.weight,
                    //     ));
                    // }
                }
            }
        }
        // if diff_i == 13 && pd.dzb == 1.0 {
        //     stat.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
        //     dbg!(stat);
        // }
        (
            point_diff / sum - point / (sum * sum) * sum_diff,
            dir.normalize_diff(dir_diff),
        )
    }

    #[allow(dead_code)]
    fn draw_2d_image<F: Fn(Line) -> f64>(
        &self,
        name: &str,
        traj: Line,
        size_m: f64,
        contrast_k: f64,
        eval_fn: F,
    ) {
        fn error_to_color(err: f64, k: f64) -> Rgb<u8> {
            let err = FRAC_2_PI * (err * k).atan();
            Rgb([(err * 255.0) as u8, ((1.0 - err) * 255.0) as u8, 0])
        }

        let initial_spherical = Geodetic::from_geocentric_cartesian(traj.point, 10);
        let d_lat = size_m * 0.5 / (EARTH_R * initial_spherical.lat.cos());
        let d_lon = size_m * 0.5 / EARTH_R;
        let lat = initial_spherical.lat;
        let lon = initial_spherical.lon;

        let mut min = f64::INFINITY;
        let mut max = -f64::INFINITY;

        struct Flat2D<const W: usize, const H: usize> {
            mem: Vec<f64>,
        }
        impl<const W: usize, const H: usize> Flat2D<W, H> {
            fn new() -> Self {
                Self {
                    mem: vec![0.0; W * H],
                }
            }

            fn get(&mut self, x: usize, y: usize) -> &mut f64 {
                &mut self.mem[x + y * W]
            }
        }

        const DIM: usize = 4000;
        let mut err_map = Flat2D::<DIM, DIM>::new();

        for x in 0..DIM {
            let lon = lon + lerp(-d_lon, d_lon, x as f64 / DIM as f64);
            for y in 0..DIM {
                let lat = lat + lerp(d_lat, -d_lat, y as f64 / DIM as f64);

                let sph = Geodetic {
                    lat,
                    lon,
                    ..initial_spherical
                };
                let mut copy = traj;
                copy.point = sph.into_geocentric_cartesian();

                let error = eval_fn(copy);
                *err_map.get(x, y) = error;

                min = min.min(error);
                max = max.max(error);
            }
        }

        let mut img = RgbImage::new(DIM as u32, DIM as u32);

        const WHITE: Rgb<u8> = Rgb([u8::MAX, u8::MAX, u8::MAX]);
        const BLUE: Rgb<u8> = Rgb([0, 0, u8::MAX]);
        const BLACK: Rgb<u8> = Rgb([0, 0, 0]);

        for x in 0..DIM {
            for y in 0..DIM {
                let e = *err_map.get(x, y);
                let e_percent = (e - min) / min * 100.0;
                let color = if e - min < 1e-6 {
                    WHITE
                } else if e_percent % 1.0 < 0.075 {
                    BLUE
                } else {
                    error_to_color(e - min, contrast_k)
                };
                img.put_pixel(x as u32, y as u32, color);
            }
        }

        // Fill a pixel for each observer
        for s in &self.data.samples {
            let x = ((s.geo_location.lon - initial_spherical.lon) / d_lon * 500.0 + 500.0).round()
                as i32;
            let y = ((s.geo_location.lat - initial_spherical.lat) / d_lat * -500.0 + 500.0).round()
                as i32;
            let mut put_pixel = |x: i32, y: i32| {
                if (0..DIM as i32).contains(&x) && (0..DIM as i32).contains(&y) {
                    img.put_pixel(x as u32, y as u32, BLACK);
                }
            };
            put_pixel(x, y);
            put_pixel(x, y + 1);
            put_pixel(x + 1, y);
            put_pixel(x + 1, y + 1);
        }

        let path = format!(
            "plots/{name}-{}-2d-{}.png",
            self.data.name.as_deref().unwrap_or("<unknown>"),
            random::<u32>()
        );
        eprintln!("saving to {path}");
        img.save(path).unwrap();
    }

    fn generic_iterative_search<F: Fn(Line) -> f64 + Send + Copy>(
        name: &str,
        mut best_traj: Line,
        dir_iters: usize,
        dir_sigma: f64,
        pos_iters: usize,
        pos_sigma: f64,
        eval_fn: F,
    ) -> Line {
        let threads_cnt = std::thread::available_parallelism().unwrap().get();

        let mut best_error = eval_fn(best_traj);

        for iter_i in 1.. {
            eprintln!("{name} Iteration #{iter_i}");
            let prev_iter_err = best_error;
            let last_best = best_traj;

            // Adjust direction
            std::thread::scope(|s| {
                (0..threads_cnt)
                    .map(|_| {
                        s.spawn(move || {
                            let mut rng = SmallRng::from_entropy();
                            let mut local_best = best_traj;
                            let mut local_best_error = best_error;

                            for _ in 0..(dir_iters / threads_cnt) {
                                let mut copy = best_traj;
                                copy.direction = copy.direction.tilt_random(dir_sigma, &mut rng);
                                let error = eval_fn(copy);
                                if error < local_best_error {
                                    local_best_error = error;
                                    local_best = copy;
                                }
                            }

                            (local_best, local_best_error)
                        })
                    })
                    .collect::<Vec<_>>()
                    .into_iter()
                    .map(|handle| handle.join().unwrap())
                    .for_each(|(traj, error)| {
                        if error < best_error {
                            best_error = error;
                            best_traj = traj;
                        }
                    });
            });
            eprintln!(
                "  direction adjusted by {:.1}{DEGREE_SYM} ({:.1}% error)",
                last_best
                    .direction
                    .dot(*best_traj.direction)
                    .min(1.0)
                    .acos()
                    .to_degrees(),
                (best_error - prev_iter_err) / prev_iter_err * 100.0,
            );

            let initial_spherical = Geodetic::from_geocentric_cartesian(best_traj.point, 10);
            let last_best = best_traj;
            let d_lat = pos_sigma / (EARTH_R * initial_spherical.lat.cos());
            let d_lon = pos_sigma / EARTH_R;
            let err_before_pos_adj = best_error;

            // Adjust position
            std::thread::scope(|s| {
                (0..threads_cnt)
                    .map(|_| {
                        s.spawn(move || {
                            let mut rng = SmallRng::from_entropy();
                            let mut local_best = best_traj;
                            let mut local_best_error = best_error;

                            for _ in 0..(pos_iters / threads_cnt) {
                                let h_mid = lerp(40_000.0, initial_spherical.h, 0.4);
                                let dh = (h_mid - 40_000.0).abs() / 2.0 + 5_000.0;

                                let mut geodetic = initial_spherical;
                                geodetic.lat += d_lat * rng.sample::<f64, _>(StandardNormal);
                                geodetic.lon += d_lon * rng.sample::<f64, _>(StandardNormal);
                                geodetic.h = h_mid + dh * rng.sample::<f64, _>(StandardNormal);

                                let mut copy = best_traj;
                                copy.point = geodetic.into_geocentric_cartesian();

                                let error = eval_fn(copy);
                                if error < local_best_error {
                                    local_best_error = error;
                                    local_best = copy;
                                }
                            }

                            (local_best, local_best_error)
                        })
                    })
                    .collect::<Vec<_>>()
                    .into_iter()
                    .map(|handle| handle.join().unwrap())
                    .for_each(|(traj, error)| {
                        if error < best_error {
                            best_error = error;
                            best_traj = traj;
                        }
                    });
            });
            eprintln!(
                "  position adjusted by {:.1}km ({:.1}% error)",
                (last_best.point - best_traj.point).norm() * 1e-3,
                (best_error - err_before_pos_adj) / err_before_pos_adj * 100.0,
            );

            if (prev_iter_err - best_error) / best_error < 0.001 {
                break;
            }
        }

        best_traj
    }

    /// Perform the "least-median-of-squares"
    pub fn lms(&self, best: Line) -> Line {
        // self.draw_2d_image("lms", best, 600_000.0, 60.0, |traj| {
        //     self.eval_median_of_squares(traj)
        // });

        Self::generic_iterative_search(
            "LMS",
            best,
            3_000_000,
            35f64.to_radians(),
            6_000_000,
            150_000.0,
            |traj| self.eval_median_of_squares(traj),
        )
    }

    /// Perform the "weidghted least-squares"
    pub fn ls(&self, traj: Line, weights: Option<&[Weight]>) -> Line {
        // self.draw_2d_image("ls", traj, 500_000.0, 100.0, |traj| {
        //     self.eval_mean_of_squares(traj, weights)
        // });

        Self::generic_iterative_search(
            "LS",
            traj,
            5_000_000,
            15f64.to_radians(),
            6_000_000,
            50_000.0,
            |traj| self.eval_mean_of_squares(traj, weights),
        )
    }

    pub fn gradient_descent_point(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        iterations: u32,
    ) -> (Line, GradientDescentStats) {
        const K: f64 = -1e10;
        const MAX_JUMP: f64 = 1.0;

        let mut cycle = CycleDetector::default();
        let mut stats = GradientDescentStats::default();

        while stats.iters < iterations && !cycle.had(&traj.point) {
            stats.iters += 1;
            cycle.push(traj.point);

            let grad = self.eval_point_grad_mean_of_squares(traj, weights);

            let mut delta = grad * K;
            let delta_dist = delta.norm();
            stats.mean_unclipped_jump += delta_dist;
            if delta_dist > MAX_JUMP {
                stats.clipped_jumps_percent += 1.0;
                delta *= MAX_JUMP / delta_dist;
            }

            traj.point += delta;
        }

        stats.clipped_jumps_percent *= 100.0 / stats.iters as f64;
        stats.mean_unclipped_jump /= stats.iters as f64;

        (traj, stats)
    }

    pub fn gradient_descent_dir(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        iterations: u32,
    ) -> (Line, GradientDescentStats) {
        const K: f64 = -1e-2;
        const MAX_JUMP: f64 = 0.00002;

        let mut cycle = CycleDetector::default();
        let mut stats = GradientDescentStats::default();

        while stats.iters < iterations && !cycle.had(&*traj.direction) {
            cycle.push(*traj.direction);
            stats.iters += 1;

            let grad = self.eval_dir_grad_mean_of_squares(traj, weights);

            let mut delta = grad * K;
            let delta_dist = delta.norm();
            stats.mean_unclipped_jump += delta_dist;
            if delta_dist > MAX_JUMP {
                stats.clipped_jumps_percent += 1.0;
                delta *= MAX_JUMP / delta_dist;
            }
            traj.direction = UnitVec3::new_normalize(traj.direction + delta);
        }

        stats.clipped_jumps_percent *= 100.0 / stats.iters as f64;
        stats.mean_unclipped_jump /= stats.iters as f64;

        (traj, stats)
    }

    pub fn gradient_descent_complete(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        params: &GradDescentParams,
    ) -> (Line, usize) {
        const ITERS_POINT: u32 = 1_000;
        const ITERS_DIR: u32 = 1_000;

        let plotter = params
            .realtime_graph
            .then(|| plotter::Builder::new().build());
        let mut err = self.eval_mean_of_squares(traj, weights);

        let init_traj = traj;

        for runs in 1.. {
            let prev_iter_err = err;

            if let Some(p) = &plotter {
                p.set_meta(format!(
                    "{DELTA_SYM}p={:.3}km\n{DELTA_SYM}{ALPHA_SYM}={:.3}{DEGREE_SYM}\n{}",
                    (init_traj.point - traj.point).norm() * 1e-3,
                    init_traj
                        .direction
                        .dot(*traj.direction)
                        .clamp(-1.0, 1.0)
                        .acos()
                        .to_degrees(),
                    self.data.name.as_deref().unwrap_or("N/A"),
                ));
                p.push(err as f32);
            }

            {
                let (new_traj, stats) = self.gradient_descent_point(traj, weights, ITERS_POINT);
                let new_err = self.eval_mean_of_squares(new_traj, weights);
                if params.debug {
                    let g = Geodetic::from_geocentric_cartesian(new_traj.point, 20);
                    let p_grad = self.eval_point_grad_mean_of_squares(new_traj, weights);
                    let derr_percents = (new_err - err) / err * 100.0;
                    let mean_zenith_coef = self
                        .data
                        .samples
                        .iter()
                        .map(|s| {
                            let ek = UnitVec3::new_normalize(new_traj.point - s.location);
                            ZenithCoef::new(s, ek).coef
                        })
                        .sum::<f64>()
                        / self.data.samples.len() as f64;
                    let grad_magnitude_m = p_grad.norm();
                    eprintln!(
                        "PNT ({iters:4}) ({mean_jump:8.2e}m {clipped:5.1}% clipped) | Err: {derr_percents:15.12}% | Mean Z_K: {mean_zenith_coef:.2} | {DELTA_SYM}={:.6}m | |p_grad| = {grad_magnitude_m:.3e}m | {g}",
                        (new_traj.point - traj.point).norm(),
                        iters = stats.iters,
                        mean_jump = stats.mean_unclipped_jump,
                        clipped = stats.clipped_jumps_percent,
                        derr_percents = eval_err_colored(derr_percents),
                    );

                    // if derr_percents > 0.00000053 {
                    //     // let weight_out = std::fs::File::create("/tmp/w.json").unwrap();
                    //     // let w = weights.unwrap();
                    //     // serde_json::to_writer(weight_out, w).unwrap();
                    //     dbg!(traj, err, new_err);
                    //     panic!();
                    // }
                }

                traj = new_traj;
                err = new_err;
            }

            {
                let (new_traj, stats) = self.gradient_descent_dir(traj, weights, ITERS_DIR);
                let new_err = self.eval_mean_of_squares(new_traj, weights);
                if params.debug {
                    let d_grad = self.eval_dir_grad_mean_of_squares(new_traj, weights);
                    let derr_percents = (new_err - err) / err * 100.0;
                    let mean_zenith_coef = self
                        .data
                        .samples
                        .iter()
                        .map(|s| {
                            let ek = UnitVec3::new_normalize(new_traj.point - s.location);
                            ZenithCoef::new(s, ek).coef
                        })
                        .sum::<f64>()
                        / self.data.samples.len() as f64;
                    let grad_magnitude_m = d_grad.norm();
                    eprintln!(
                        "DIR ({iters:4}) ({mean_jump:8.2e}m {clipped:5.1}% clipped) | Err: {derr_percents:15.12}% | Mean Z_K: {mean_zenith_coef:.2} | {DELTA_SYM}={:.6}{DEGREE_SYM} | |d_grad| = {grad_magnitude_m:.3e}m",
                        new_traj
                            .direction
                            .dot(*traj.direction)
                            .min(1.0)
                            .acos()
                            .to_degrees(),
                        iters = stats.iters,
                        mean_jump = stats.mean_unclipped_jump,
                        clipped = stats.clipped_jumps_percent,
                        derr_percents = eval_err_colored(derr_percents),
                    );

                    // if derr_percents > 0.000003 {
                    //     let weight_out = std::fs::File::create("/tmp/w.json").unwrap();
                    //     let w = weights.unwrap();
                    //     serde_json::to_writer(weight_out, w).unwrap();
                    //     dbg!(traj, err, new_err);
                    //     panic!();
                    // }
                }

                traj = new_traj;
                err = new_err;
            }

            if err <= prev_iter_err
                && (prev_iter_err - err) / prev_iter_err <= params.target_err_rel
            {
                return (traj, runs);
            }
        }

        unreachable!()
    }

    pub fn compute_weights(&self, traj: Line) -> Vec<Weight> {
        let evals = self
            .data
            .samples
            .iter()
            .map(|s| self.evaluate_traj(s, traj))
            .collect::<Vec<Evaluation>>();

        let med_mul = 1.483 * (1. + 5. / (self.data.samples.len() - 6) as f64);

        let mut buf = Vec::with_capacity(self.data.samples.len());
        buf.extend(evals.iter().flat_map(|eval| eval.h_end.map(|x| x * x)));
        let h_end_med = buf.median().sqrt() * med_mul;
        buf.clear();
        buf.extend(evals.iter().flat_map(|eval| eval.z_end.map(|x| x * x)));
        let z_end_med = buf.median().sqrt() * med_mul;
        buf.clear();
        buf.extend(evals.iter().flat_map(|eval| eval.da.map(|x| x * x)));
        let da_med = buf.median().sqrt() * med_mul;

        let get_weight = |err: f64, best_err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const O: f64 = 1.50;
            const F: f64 = 0.40;
            0.5 * (1.0 - ((err.abs() / best_err - O) / F).tanh())
        };
        let get_weight_h_end = |err| get_weight(err, h_end_med);
        let get_weight_z_end = |err| get_weight(err, z_end_med);
        let get_weight_da = |err| get_weight(err, da_med);

        let mut weights: Vec<Weight> = Vec::with_capacity(self.data.samples.len());
        let mut hw = Vec::new();
        let mut zw = Vec::new();
        let mut daw = Vec::new();
        let mut he = Vec::new();
        let mut ze = Vec::new();
        let mut dae = Vec::new();
        for eval in &evals {
            let b1 = eval.h_end.map_or(0.0, get_weight_h_end);
            let b2 = eval.z_end.map_or(0.0, get_weight_z_end);
            let b3 = eval.da.map_or(0.0, get_weight_da);
            weights.push(Weight {
                h_end: b1,
                z_end: b2,
                da: b3,
            });
            eval.h_end.map(|x| he.push(x.abs() / h_end_med));
            eval.z_end.map(|x| ze.push(x.abs() / z_end_med));
            eval.da.map(|x| dae.push(x.abs() / da_med));
            eval.h_end.map(|_| hw.push(b1));
            eval.z_end.map(|_| zw.push(b2));
            eval.da.map(|_| daw.push(b3));
        }

        println!("Error distribution");
        println!();
        if !hw.is_empty() {
            println!("H-End ({}):", he.len());
            draw_hitogram(&he, 5);
        }
        if !zw.is_empty() {
            println!("Z-End ({}):", ze.len());
            draw_hitogram(&ze, 5);
        }
        if !daw.is_empty() {
            println!("DA ({}):", dae.len());
            draw_hitogram(&dae, 5);
        }
        println!();

        println!("Weight distribution");
        println!();
        if !hw.is_empty() {
            println!("H-End ({}):", hw.len());
            draw_hitogram(&hw, 5);
        }
        if !zw.is_empty() {
            println!("Z-End ({}):", zw.len());
            draw_hitogram(&zw, 5);
        }
        if !daw.is_empty() {
            println!("DA ({}):", daw.len());
            draw_hitogram(&daw, 5);
        }
        println!();

        weights
    }

    pub fn eval_median_of_squares(&self, traj: Line) -> f64 {
        let mut buf = Vec::with_capacity(self.data.samples.len() * 3);
        buf.extend(
            self.data
                .samples
                .iter()
                .flat_map(|s| self.evaluate_traj(s, traj).iter_squared()),
        );
        buf.median()
    }

    pub fn eval_mean_of_squares(&self, traj: Line, weights: Option<&[Weight]>) -> f64 {
        let mut sum = 0.0;
        let mut cnt = 0.0;

        match weights {
            Some(weights) => {
                for (s, w) in self.data.iter().zip(weights) {
                    self.evaluate_traj(s, traj)
                        .weighted_sum_sq(w, &mut sum, &mut cnt);
                }
            }
            None => {
                for s in self.data.iter() {
                    self.evaluate_traj(s, traj).sum_sq(&mut sum, &mut cnt);
                }
            }
        }

        sum / cnt
    }

    pub fn eval_point_grad_mean_of_squares(&self, traj: Line, weights: Option<&[Weight]>) -> Vec3 {
        let mut half_sum = Vec3::default();
        let mut cnt = 0.0;

        match weights {
            Some(weights) => {
                for (s, w) in self.data.iter().zip(weights) {
                    let (eval, grad) = self.evaluate_point_grad_traj(s, traj);
                    eval.weighted_sum_sq_grad(&grad, w, &mut half_sum, &mut cnt);
                }
            }
            None => {
                for s in self.data.iter() {
                    let (eval, grad) = self.evaluate_point_grad_traj(s, traj);
                    eval.sum_sq_grad(&grad, &mut half_sum, &mut cnt);
                }
            }
        }

        half_sum / cnt * 2.0
    }

    pub fn eval_dir_grad_mean_of_squares(&self, traj: Line, weights: Option<&[Weight]>) -> Vec3 {
        let mut half_sum = Vec3::default();
        let mut cnt = 0.0;

        match weights {
            Some(weights) => {
                for (s, w) in self.data.iter().zip(weights) {
                    let (eval, grad) = self.evaluate_dir_grad_traj(s, traj);
                    eval.weighted_sum_sq_grad(&grad, w, &mut half_sum, &mut cnt);
                }
            }
            None => {
                for s in self.data.iter() {
                    let (eval, grad) = self.evaluate_dir_grad_traj(s, traj);
                    eval.sum_sq_grad(&grad, &mut half_sum, &mut cnt);
                }
            }
        }

        half_sum / cnt * 2.0
    }

    fn _draw_plots(&self, traj: Line, weights: Option<&[(f64, f64, f64)]>, stage: &str) {
        let mut points = Vec::new();
        let mut buf = Vec::new();

        let get_w = |i: usize| weights.map(|w| w[i]).unwrap_or((0.0, 0.0, 0.0));

        points.clear();
        buf.clear();
        for (i, s) in self.data.samples.iter().enumerate() {
            if let Some(da) = self.evaluate_traj(s, traj).da {
                let x = (i + 1) as f64;
                let y = da.to_degrees();
                buf.push(y);
                points.push((x, y, weight_to_rgb(get_w(i).2), 2.5));
            }
        }
        let med = buf.median();
        eprintln!("DA_error med = {med}");
        draw_plot_svg(
            &format!("plots/da-errors-{:.0}-{stage}.svg", traj.point.x),
            &points,
            &[],
            &[(0.0, colors::BLACK), (med, colors::BLUE)],
            &[],
        )
        .unwrap();
    }

    fn calc_sigma_sample(
        &self,
        sample_i: usize,
        flash: Vec3,
        dir: Azimuthal,
        speed: f64,
        traj: Line,
        weights: Option<&[Weight]>,
    ) -> Sigmas {
        const DA_D: f64 = to_radians(0.1);
        const AZ_D: f64 = to_radians(0.1);
        const ALT_D: f64 = to_radians(0.1);

        let mut this = self.clone();

        let eval = this.evaluate_traj(&this.data.samples[sample_i], traj);
        let mut s = Sigmas::default();

        if let Some(da_err) = eval.da {
            let old_da = this.data.samples[sample_i].da;
            *this.data.samples[sample_i].da.as_mut().unwrap() += DA_D;

            let (new_traj, _gd_i) = this.gradient_descent_complete(traj, weights, &GD_SIGMA_PARAMS);
            let (new_flash, new_speed) = this.calc_flash_and_speed(new_traj, weights);

            this.data.samples[sample_i].da = old_da;

            let mul = da_err * da_err / DA_D / DA_D;
            s.x += mul * f64::powi(new_flash.x - flash.x, 2);
            s.y += mul * f64::powi(new_flash.y - flash.y, 2);
            s.z += mul * f64::powi(new_flash.z - flash.z, 2);
            s.v_angle +=
                mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
            s.speed += mul * f64::powi(new_speed - speed, 2);
        }

        if let Some(z0_err) = eval.z_end {
            let old_z0 = this.data.samples[sample_i].z0;
            let old_k_end = this.data.samples[sample_i].k_end;
            *this.data.samples[sample_i].z0.as_mut().unwrap() += AZ_D;
            this.data.samples[sample_i].update_k_end();

            let (new_traj, _gd_i) = this.gradient_descent_complete(traj, weights, &GD_SIGMA_PARAMS);
            let (new_flash, new_speed) = this.calc_flash_and_speed(new_traj, weights);

            this.data.samples[sample_i].z0 = old_z0;
            this.data.samples[sample_i].k_end = old_k_end;

            let mul = z0_err * z0_err / AZ_D / AZ_D;
            s.x += mul * f64::powi(new_flash.x - flash.x, 2);
            s.y += mul * f64::powi(new_flash.y - flash.y, 2);
            s.z += mul * f64::powi(new_flash.z - flash.z, 2);
            s.v_angle +=
                mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
            s.speed += mul * f64::powi(new_speed - speed, 2);
        }

        if let Some(h0_err) = eval.h_end {
            let old_h0 = this.data.samples[sample_i].h0;
            let old_k_end = this.data.samples[sample_i].k_end;
            *this.data.samples[sample_i].h0.as_mut().unwrap() += ALT_D;
            this.data.samples[sample_i].update_k_end();

            let (new_traj, _gd_i) = this.gradient_descent_complete(traj, weights, &GD_SIGMA_PARAMS);
            let (new_flash, new_speed) = this.calc_flash_and_speed(new_traj, weights);

            this.data.samples[sample_i].h0 = old_h0;
            this.data.samples[sample_i].k_end = old_k_end;

            let mul = h0_err * h0_err / AZ_D / AZ_D;
            s.x += mul * f64::powi(new_flash.x - flash.x, 2);
            s.y += mul * f64::powi(new_flash.y - flash.y, 2);
            s.z += mul * f64::powi(new_flash.z - flash.z, 2);
            s.v_angle +=
                mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
            s.speed += mul * f64::powi(new_speed - speed, 2);
        }

        s
    }

    fn calc_sigmas(&self, traj: Line, weights: Option<&[Weight]>) -> Sigmas {
        let (flash, speed) = self.calc_flash_and_speed(traj, weights);
        let dir: Azimuthal = traj.direction.into_inner().into();

        (0..self.data.samples.len())
            .into_par_iter()
            .map(|i| self.calc_sigma_sample(i, flash, dir, speed, traj, weights))
            .sum::<Sigmas>()
            .sqrt()
    }

    fn flip_da(&mut self, traj: Line) -> usize {
        let mut flipped = 0;

        for i in 0..self.data.samples.len() {
            let Some(dda) = self.evaluate_traj(&self.data.samples[i], traj).da else {
                continue;
            };

            if dda.abs() > f64::to_radians(120.0) {
                let da = self.data.samples[i].da.as_mut().unwrap();
                *da = (*da + PI) % TAU;
                flipped += 1;
            }
        }

        flipped
    }

    fn run_stage<F>(&mut self, f: F) -> Line
    where
        F: FnOnce(&mut Self) -> (String, Line, Option<Sigmas>),
    {
        let timer = Instant::now();
        let (title, result, sigmas) = f(self);
        let duration = timer.elapsed();

        self.data.compare(result, sigmas, &title);
        eprintln!("This step took {duration:?}\n");

        result
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        let traj = self.run_stage(|this| {
            let (result, sigmas) = this.pairwise();
            ("Initial guess (pairwise)".into(), result, Some(sigmas))
        });

        self.update_altidutes_error_k(traj);

        let traj = self.run_stage(|this| {
            let (traj, gd_i) = this.gradient_descent_complete(traj, None, &GD_FIRST_PARAMS);
            (format!("Gradient descent ({gd_i} runs)"), traj, None)
        });

        let traj = self.run_stage(|this| ("After LMS".into(), this.lms(traj), None));

        let mut weights = self.compute_weights(traj);

        let mut traj =
            self.run_stage(|this| ("After LS".into(), this.ls(traj, Some(&weights)), None));

        if !self.params.no_da_flip {
            eprintln!("flipped {} descent angles", self.flip_da(traj));
            eprintln!();
            weights = self.compute_weights(traj);
            traj = self.ls(traj, Some(&weights));
            self.data.compare(traj, None, "LS (after DA-flip)");
        }

        dbg!(traj);

        let traj = self.run_stage(|this| {
            let (traj, gd_i) =
                this.gradient_descent_complete(traj, Some(&weights), &GD_FINAL_PARAMS);
            (format!("Gradient descent ({gd_i} runs)"), traj, None)
        });

        dbg!(traj);
        dump_weights(&weights);

        let sigmas = self.calc_sigmas(traj, Some(&weights));
        self.data.compare(traj, Some(sigmas), "Final");

        let (flash, speed) = self.calc_flash_and_speed(traj, Some(&weights));

        Solution {
            flash: Geodetic::from_geocentric_cartesian(flash, 10),
            velocity: traj.direction.into_inner() * speed,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, traj: Line, _weights: Option<&[Weight]>) -> (Vec3, f64) {
        let mut lambdas = Vec::new();
        let mut speeds = Vec::new();

        for s in &self.data.samples {
            let Some(k_end) = s.k_end else { continue };

            // total_l_end += 1;

            // A vertor that describes the hemispace of "allowed" observations
            let perp = traj
                .direction
                .cross(traj.point - s.location)
                .cross(*traj.direction);
            if perp.dot(*k_end) <= 0.0 {
                // skipped_l_end_dir += 1;
                continue;
            }

            let l_end = lambda(s.location, k_end, traj.point, traj.direction);

            // if weights.is_some_and(|w| w[i].end > 0.5) {
            lambdas.push(l_end);
            // lambdas_all.push(l_end);
            // } else {
            // lambdas_all.push(l_end);
            // skipped_l_end_w += 1;
            // }

            if s.observation_matches(traj) {
                if let Some((duration, k_start)) = s.dur.zip(s.k_start) {
                    let l_start = lambda(s.location, k_start, traj.point, traj.direction);
                    let dist = l_end - l_start;
                    speeds.push(dist / duration);
                }
            }
        }

        let l_end = lambdas.median();

        // dbg!(lambdas);

        let speed = speeds.median();
        (traj.point + traj.direction * l_end, speed)
    }

    fn update_altidutes_error_k(&mut self, traj: Line) {
        self.altitude_error_multiplier = 1.0;
        let mut h_end = Vec::new();
        let mut z_end = Vec::new();
        let mut da = Vec::new();

        for s in &self.data.samples {
            let e = self.evaluate_traj(s, traj);
            if let Some(x) = e.h_end {
                h_end.push(x.abs());
            }
            if let Some(x) = e.z_end {
                z_end.push(x.abs());
            }
            if let Some(x) = e.da {
                da.push(x.abs());
            }
        }

        if z_end.is_empty() || h_end.is_empty() || da.is_empty() {
            println!("altitude_error_multiplier not updated");
            return;
        }

        let z_end = z_end.median();
        let h_end = h_end.median();
        let da = da.median();

        let avg = (z_end + da) / 2.0;

        self.altitude_error_multiplier = avg / h_end;
        println!(
            "altitude_error_multiplier = {:.1}",
            self.altitude_error_multiplier
        );
    }

    fn evaluate_traj_with(
        &self,
        sample: &DataSample,
        traj: Line,
        ek: UnitVec3,
        zenith_k: f64,
    ) -> Evaluation {
        Evaluation {
            h_end: if !self.params.no_altitudes {
                sample.h0.map(|h0| {
                    angle_diff(h0, sample.calc_altitude(ek))
                        * zenith_k
                        * self.altitude_error_multiplier
                })
            } else {
                None
            },
            z_end: if !self.params.no_azimuths {
                sample
                    .z0
                    .map(|z0| angle_diff(z0, sample.calc_azimuth(*ek)) * zenith_k)
            } else {
                None
            },
            da: sample.da.map(|da| {
                angle_diff(descent_angle(sample.zenith_dir, ek, *traj.direction), da) * zenith_k
            }),
        }
    }

    pub fn evaluate_traj(&self, sample: &DataSample, traj: Line) -> Evaluation {
        let ek = UnitVec3::new_normalize(traj.point - sample.location);
        let zenith_k = ZenithCoef::new(sample, ek).coef;
        self.evaluate_traj_with(sample, traj, ek, zenith_k)
    }

    pub fn evaluate_point_grad_traj(
        &self,
        s: &DataSample,
        traj: Line,
    ) -> (Evaluation, Evaluation<Vec3>) {
        let k = traj.point - s.location;
        let (ek, k_norm) = UnitVec3::new_and_get(k);
        let k_norm_recip = k_norm.recip();
        let zenith_coef = ZenithCoef::new(s, ek);

        let err0 = self.evaluate_traj_with(s, traj, ek, zenith_coef.coef);

        let end_diff_x = k.dot(*s.east_dir);
        let end_diff_y = k.dot(*s.north_dir);

        let make_diff_ek_and_diff_zenith_coef = |mask: Vec3| {
            // k_norm = (k_norm^2)^1/2
            // k_norm,x = 1/2 * (k_norm^2)^(-1/2) * 2kx = kx / k_norm
            let diff_k_norm = k.dot(mask) * k_norm_recip;
            // ek = k * k_norm^-1
            // ek,x = (k,x * k_norm^-1) + (k * -1 * k_norm^-2 * k_norm,x)
            // ek,x = (k,x - (k * k_norm,x * k_norm^-1)) * * k_norm^-1
            let diff_ek = (mask - (k * diff_k_norm * k_norm_recip)) * k_norm_recip;
            (diff_ek, zenith_coef.diff(s, diff_ek))
        };

        let make_h_end_diff = |h0: f64, mask: Vec3| {
            let (diff_ek, diff_zenith_coef) = make_diff_ek_and_diff_zenith_coef(mask);

            // E_0 = 4 * diff(h0, asin(ek*zenith))
            // E_0 = 4[asin(ek*zenith) - h0]
            // E_0,i = 4 * (diff_ek * zenith) / sqrt(1-(ek*zenith)^2)
            let e0 = angle_diff(h0, ek.dot(*s.zenith_dir).clamp(-1.0, 1.0).asin())
                * self.altitude_error_multiplier;
            let e0_diff = diff_ek.dot(*s.zenith_dir) * self.altitude_error_multiplier
                / f64::sqrt(1.0 - ek.dot(*s.zenith_dir) * ek.dot(*s.zenith_dir));

            // E = E_0 * zenith_coef
            // E,x = E_0,x * zenith_coef + E_0 * zenith_coef,x
            e0_diff * zenith_coef.coef + e0 * diff_zenith_coef
        };

        let make_z_end_diff = |z0: f64, mask: Vec3| {
            let (_, diff_zenith_coef) = make_diff_ek_and_diff_zenith_coef(mask);

            // E_0 = diff(z0, atan2(k*east, k*north))
            // E_0 = atan2(k*east, k*north) - z0
            // E_0,x = -(k*east) / (xx+yy) * north_x + (k*north) / (xx+yy) * east_x
            // E_0,x = (k*north * east_x - k*east * north_x) / (xx+yy)
            let e0 = angle_diff(z0, f64::atan2(end_diff_x, end_diff_y));
            let e0_diff = (end_diff_y * s.east_dir.dot(mask) - end_diff_x * s.north_dir.dot(mask))
                / (end_diff_x * end_diff_x + end_diff_y * end_diff_y);

            // E = E_0 * zenith_coef
            // E,x = E_0,x * zenith_coef + E_0 * zenith_coef,x
            e0_diff * zenith_coef.coef + e0 * diff_zenith_coef
        };

        let make_da_diff_xa = ek.cross(*s.zenith_dir);
        let make_da_diff_ya = make_da_diff_xa.cross(*ek);
        let make_da_diff_x = traj.direction.dot(make_da_diff_xa);
        let make_da_diff_y = traj.direction.dot(make_da_diff_ya);
        let make_da_diff_dac = f64::atan2(make_da_diff_x, make_da_diff_y);
        let make_da_diff = |da: f64, mask: Vec3| -> f64 {
            let (diff_ek, diff_zenith_coef) = make_diff_ek_and_diff_zenith_coef(mask);

            let diff_xa = diff_ek.cross(*s.zenith_dir);
            let diff_ya = diff_xa.cross(*ek) + make_da_diff_xa.cross(diff_ek);

            let diff_x = traj.direction.dot(diff_xa);
            let diff_y = traj.direction.dot(diff_ya);

            let diff_dac = (make_da_diff_y * diff_x - make_da_diff_x * diff_y)
                / (make_da_diff_x * make_da_diff_x + make_da_diff_y * make_da_diff_y);

            // Err = (da - dac) * zenith_coef = da*zenith_coef - dac*zenith_coef
            // Err,i = da*zenith_coef,i - dac*zenith_coef,i - dac,i*zenith_coef
            //       = zenith_coef,i*(da - dac) - dac,i*zenith_coef
            diff_zenith_coef * angle_diff(make_da_diff_dac, da) - diff_dac * zenith_coef.coef
        };

        let h_end = (s.h0.is_some() && !self.params.no_altitudes).then(|| {
            Vec3::new(
                make_h_end_diff(s.h0.unwrap(), Vec3::x()),
                make_h_end_diff(s.h0.unwrap(), Vec3::y()),
                make_h_end_diff(s.h0.unwrap(), Vec3::z()),
            )
        });

        let z_end = (s.z0.is_some() && !self.params.no_azimuths).then(|| {
            Vec3::new(
                make_z_end_diff(s.z0.unwrap(), Vec3::x()),
                make_z_end_diff(s.z0.unwrap(), Vec3::y()),
                make_z_end_diff(s.z0.unwrap(), Vec3::z()),
            )
        });

        let da = s.da.map(|da| {
            Vec3::new(
                make_da_diff(da, Vec3::x()),
                make_da_diff(da, Vec3::y()),
                make_da_diff(da, Vec3::z()),
            )
        });

        (err0, Evaluation { h_end, z_end, da })
    }

    pub fn evaluate_dir_grad_traj(
        &self,
        s: &DataSample,
        traj: Line,
    ) -> (Evaluation, Evaluation<Vec3>) {
        let k = traj.point - s.location;
        let ek = UnitVec3::new_normalize(k);
        let zenith_coef = ZenithCoef::new(s, ek).coef;

        let err0 = self.evaluate_traj_with(s, traj, ek, zenith_coef);

        let xa = ek.cross(*s.zenith_dir);
        let ya = xa.cross(*ek);

        let x = traj.direction.dot(xa);
        let y = traj.direction.dot(ya);

        let make_da_diff_dir = |mask: Vec3| -> f64 {
            // v = {x, y, z}
            // dx
            // v' = {x+dx, y, z} / (1 + dx*dx + 2x*dx)
            // dv = ...
            // dv/dx = {1-xx, -xy, -xz}
            let diff_vel = mask - traj.direction * traj.direction.dot(mask);

            let diff_x = diff_vel.dot(xa);
            let diff_y = diff_vel.dot(ya);

            (x * diff_y - y * diff_x) / (x * x + y * y) * zenith_coef
        };

        let da = s.da.map(|_| {
            Vec3::new(
                make_da_diff_dir(Vec3::x()),
                make_da_diff_dir(Vec3::y()),
                make_da_diff_dir(Vec3::z()),
            )
        });

        (
            err0,
            Evaluation {
                // h_end does not depend on the direction
                h_end: None,
                // z_end does not depend on the direction
                z_end: None,
                da,
            },
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Evaluation<T = f64> {
    pub h_end: Option<T>,
    pub z_end: Option<T>,
    pub da: Option<T>,
}

impl ops::Sub for Evaluation {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            z_end: self.z_end.zip(rhs.z_end).map(|(x, y)| x - y),
            h_end: self.h_end.zip(rhs.h_end).map(|(x, y)| x - y),
            da: self.da.zip(rhs.da).map(|(x, y)| x - y),
        }
    }
}

impl ops::Div<f64> for Evaluation {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            z_end: self.z_end.map(|x| x / rhs),
            h_end: self.h_end.map(|x| x / rhs),
            da: self.da.map(|x| x / rhs),
        }
    }
}

impl Evaluation {
    fn iter_squared(self) -> impl Iterator<Item = f64> {
        let sq = |x| x * x;
        self.h_end
            .map(sq)
            .into_iter()
            .chain(self.z_end.map(sq))
            .chain(self.da.map(sq))
    }

    fn sum_sq(&self, sum: &mut f64, cnt: &mut f64) {
        if let Some(x) = self.h_end {
            *sum += x * x;
            *cnt += 1.0;
        }
        if let Some(x) = self.z_end {
            *sum += x * x;
            *cnt += 1.0;
        }
        if let Some(x) = self.da {
            *sum += x * x;
            *cnt += 1.0;
        }
    }

    fn weighted_sum_sq(&self, w: &Weight, sum: &mut f64, cnt: &mut f64) {
        if let Some(x) = self.h_end {
            *sum += x * x * w.h_end;
            *cnt += w.h_end;
        }
        if let Some(x) = self.z_end {
            *sum += x * x * w.z_end;
            *cnt += w.z_end;
        }
        if let Some(x) = self.da {
            *sum += x * x * w.da;
            *cnt += w.da;
        }
    }

    fn sum_sq_grad(&self, grad: &Evaluation<Vec3>, half_sum: &mut Vec3, cnt: &mut f64) {
        if let (Some(x), Some(dx)) = (self.h_end, grad.h_end) {
            *half_sum += dx * x;
            *cnt += 1.0;
        }
        if let (Some(x), Some(dx)) = (self.z_end, grad.z_end) {
            *half_sum += dx * x;
            *cnt += 1.0;
        }
        if let (Some(x), Some(dx)) = (self.da, grad.da) {
            *half_sum += dx * x;
            *cnt += 1.0;
        }
    }

    fn weighted_sum_sq_grad(
        &self,
        grad: &Evaluation<Vec3>,
        w: &Weight,
        half_sum: &mut Vec3,
        cnt: &mut f64,
    ) {
        if let (Some(x), Some(dx)) = (self.h_end, grad.h_end) {
            *half_sum += dx * x * w.h_end;
            *cnt += w.h_end;
        }
        if let (Some(x), Some(dx)) = (self.z_end, grad.z_end) {
            *half_sum += dx * x * w.z_end;
            *cnt += w.z_end;
        }
        if let (Some(x), Some(dx)) = (self.da, grad.da) {
            *half_sum += dx * x * w.da;
            *cnt += w.da;
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Weight {
    pub h_end: f64,
    pub z_end: f64,
    pub da: f64,
}

fn lerp(a: f64, b: f64, w: f64) -> f64 {
    a + (b - a) * w
}

const fn sq(x: f64) -> f64 {
    x * x
}

struct ZenithCoef {
    pub to_zenith_cos: f64,
    pub to_zenith: f64,
    pub coef: f64,
}

impl ZenithCoef {
    const ZENITH_COEF_ANGLE: f64 = to_radians(40.0);
    const ZENITH_COEF_SIN_K: f64 = FRAC_PI_2 / Self::ZENITH_COEF_ANGLE;
    const ZENITH_COEF_K: f64 = 5.0;

    fn new(s: &DataSample, ek: UnitVec3) -> Self {
        let to_zenith_cos = s.zenith_dir.dot(*ek);
        let to_zenith = to_zenith_cos.min(1.0).acos();

        let coef = if to_zenith < Self::ZENITH_COEF_ANGLE {
            Self::ZENITH_COEF_K / f64::sin(to_zenith * Self::ZENITH_COEF_SIN_K) + 1.0
                - Self::ZENITH_COEF_K
        } else if to_zenith > PI - Self::ZENITH_COEF_ANGLE {
            Self::ZENITH_COEF_K / f64::sin((PI - to_zenith) * Self::ZENITH_COEF_SIN_K) + 1.0
                - Self::ZENITH_COEF_K
        } else {
            1.0
        };

        Self {
            to_zenith_cos,
            to_zenith,
            coef,
        }
    }

    fn diff(&self, s: &DataSample, ek_diff: Vec3) -> f64 {
        // k(x) = 1/sin(x * ZENITH_COEF_SIN_K)
        // k'(x) = a * ZENITH_COEF_SIN_K / (a^2 - 1), where a = cos(x * ZENITH_COEF_SIN_K)

        let to_zenith_cos_diff = s.zenith_dir.dot(ek_diff);

        let to_zenith_diff =
            -1.0 / f64::sqrt(1.0 - self.to_zenith_cos * self.to_zenith_cos) * to_zenith_cos_diff;

        if self.to_zenith < Self::ZENITH_COEF_ANGLE {
            let a = f64::cos(self.to_zenith * Self::ZENITH_COEF_SIN_K);
            a * Self::ZENITH_COEF_K * Self::ZENITH_COEF_SIN_K * to_zenith_diff / (a * a - 1.0)
        } else if self.to_zenith > PI - Self::ZENITH_COEF_ANGLE {
            let a = f64::cos((PI - self.to_zenith) * Self::ZENITH_COEF_SIN_K);
            a * Self::ZENITH_COEF_K * Self::ZENITH_COEF_SIN_K * to_zenith_diff / (1.0 - a * a)
        } else {
            0.0
        }
    }
}

#[derive(Debug, Default)]
pub struct GradientDescentStats {
    pub iters: u32,
    pub clipped_jumps_percent: f64,
    pub mean_unclipped_jump: f64,
}

/// Detect cycles of stats.
///
/// Current implementation only detects cycles of up to 10 states. It is efficient and in most cases
/// sufficient.
///
/// It is assumed that `T::default()` is not a "used" value.
pub struct CycleDetector<T> {
    history: [T; 10],
}

impl<T: Default + Copy> Default for CycleDetector<T> {
    fn default() -> Self {
        Self {
            history: [T::default(); 10],
        }
    }
}

impl<T: PartialEq + Copy> CycleDetector<T> {
    pub fn had(&self, v: &T) -> bool {
        self.history.contains(v)
    }

    pub fn push(&mut self, v: T) {
        self.history.rotate_right(1);
        self.history[0] = v;
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GradDescentParams {
    pub target_err_rel: f64,
    pub debug: bool,
    pub realtime_graph: bool,
}

fn dump_weights(weights: &[Weight]) {
    let weight_out = std::fs::File::create("./weights.json").unwrap();
    serde_json::to_writer(weight_out, weights).unwrap();
}

struct RangeIter {
    start: f64,
    end: f64,
    step: f64,
}

impl Iterator for RangeIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start <= self.end {
            let retval = self.start;
            self.start += self.step;
            Some(retval)
        } else {
            None
        }
    }
}

fn range(start: f64, end: f64, step: f64) -> RangeIter {
    RangeIter { start, end, step }
}
