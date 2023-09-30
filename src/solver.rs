//! The implementation

use std::f64::consts::*;
use std::ops;
use std::time::Instant;

use common::constants::*;
use common::histogram::draw_hitogram;
use common::maths::*;
use common::obs_data::{Data, DataSample};
use common::plot::{draw_plot_svg, plotters, weight_to_rgb};
use common::quick_median::{SliceExt, SliceExt2};
use common::rand::random;
use common::structs::*;
use common::{rand, rand_distr};

use image::{Rgb, RgbImage};

use plotters::style::colors;
use rand::prelude::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::StandardNormal;
use rayon::prelude::*;

pub mod pair;
use pair::PairTrajectory;

const GD_FIRST_ERR_DIFF: f64 = 1e-8;
const GD_TARGET_ERR_DIFF: f64 = 1e-11;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Parameters to tweak the algorithm
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub da_only: bool,
    pub no_da_flip: bool,
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
        Solver { data, params }
    }

    /// The weighted sum of pairwise plane crossings
    pub fn pairwise(&self) -> Line {
        let mut point = Vec3::new(0.0, 0.0, 0.0);
        let mut dir = Vec3::new(0.0, 0.0, 0.0);
        let mut sum = 0.0;
        for (i, s1) in self.data.iter().enumerate() {
            for s2 in self.data.iter().skip(i + 1) {
                if let Some(pair) = PairTrajectory::calculate(s1, s2) {
                    point += pair.line.point * pair.weight;
                    dir += pair.line.direction.into_inner() * pair.weight;
                    sum += pair.weight;
                }
            }
        }
        Line {
            point: point / sum,
            direction: UnitVec3::new_normalize(dir),
        }
    }

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
            let prev_iter_err = best_error;
            let last_best = best_traj;
            let d_lat = pos_sigma / (EARTH_R * initial_spherical.lat.cos());
            let d_lon = pos_sigma / EARTH_R;

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
                (best_error - prev_iter_err) / prev_iter_err * 100.0,
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
        iterations: usize,
    ) -> Line {
        const K: f64 = -1e9;
        const MAX_JUMP: f64 = 100.0;

        let mut maxed = 0;

        for _ in 0..iterations {
            let (_, grad) = self.eval_grad_mean_of_squares(traj, weights);

            let mut delta = grad.point * K;
            let delta_dist = delta.norm();
            if delta_dist > MAX_JUMP {
                maxed += 1;
                delta *= MAX_JUMP / delta_dist;
            }
            traj.point += delta;
        }

        eprint!("maxed {maxed} ");

        traj
    }

    pub fn gradient_descent_dir(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        iterations: usize,
    ) -> Line {
        const K: f64 = -0.001;
        // const K: f64 = -1e9;
        const MAX_JUMP: f64 = 0.01;

        let mut maxed = 0;

        for _ in 0..iterations {
            let (_, grad) = self.eval_grad_mean_of_squares(traj, weights);

            let mut delta = grad.direction * K;
            let delta_dist = delta.norm();
            if delta_dist > MAX_JUMP {
                maxed += 1;
                delta *= MAX_JUMP / delta_dist;
            }
            traj.direction = UnitVec3::new_normalize(traj.direction + delta);
        }

        eprint!("maxed {maxed} ");

        traj
    }

    pub fn gradient_descent_complete(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        target_err_diff: f64,
    ) -> (Line, usize) {
        const ITERS_POINT: usize = 100_000;
        const ITERS_DIR: usize = 50_000;

        let mut err = self.eval_mean_of_squares(traj, weights);
        dbg!(err);

        for runs in 1.. {
            let prev_iter_err = err;

            {
                let new_traj = self.gradient_descent_point(traj, weights, ITERS_POINT);
                let g = Geodetic::from_geocentric_cartesian(new_traj.point, 20);
                let (new_err, grad) = self.eval_grad_mean_of_squares(new_traj, weights);
                let derr_percents = (new_err - err) / err * 100.0;
                let mean_zenith_coef = self
                    .data
                    .samples
                    .iter()
                    .map(|s| {
                        let ek = UnitVec3::new_normalize(new_traj.point - s.location);
                        Self::calc_zenith_coef(s, ek)
                    })
                    .sum::<f64>()
                    / self.data.samples.len() as f64;
                let grad_magnitude_m = grad.point.norm();
                eprintln!(
                    "PNT {g} | Err: {derr_percents:.12}% | Mean Z_K: {mean_zenith_coef:.2} | {DELTA_SYM}={:.6}m | |p_grad| = {grad_magnitude_m:.3e}m",
                    (new_traj.point - traj.point).norm()
                );

                traj = new_traj;
                err = new_err;
            }

            {
                let new_traj = self.gradient_descent_dir(traj, weights, ITERS_DIR);
                let g = Geodetic::from_geocentric_cartesian(new_traj.point, 20);
                let (new_err, grad) = self.eval_grad_mean_of_squares(new_traj, weights);
                let derr_percents = (new_err - err) / err * 100.0;
                let mean_zenith_coef = self
                    .data
                    .samples
                    .iter()
                    .map(|s| {
                        let ek = UnitVec3::new_normalize(new_traj.point - s.location);
                        Self::calc_zenith_coef(s, ek)
                    })
                    .sum::<f64>()
                    / self.data.samples.len() as f64;
                let grad_magnitude_m = grad.direction.norm();
                eprintln!(
                    "DIR {g} | Err: {derr_percents:.12}% | Mean Z_K: {mean_zenith_coef:.2} | {DELTA_SYM}={:.6}{DEGREE_SYM} | |d_grad| = {grad_magnitude_m:.3e}m",
                    new_traj
                        .direction
                        .dot(*traj.direction)
                        .min(1.0)
                        .acos()
                        .to_degrees(),
                );

                traj = new_traj;
                err = new_err;
            }

            let new_eval = self.eval_mean_of_squares(traj, weights);
            dbg!(new_eval);
            if new_eval <= prev_iter_err
                && (prev_iter_err - new_eval) / prev_iter_err <= target_err_diff
            {
                return (traj, runs);
            }
        }

        unreachable!()
    }

    pub fn compute_weights(&self, traj: Line) -> Vec<Weight> {
        let best_err = self.eval_median_of_squares(traj).sqrt()
            * 1.483
            * (1. + 5. / (self.data.samples.len() - 6) as f64);
        dbg!(best_err.to_degrees());
        let get_weight = |err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const O: f64 = 0.90;
            const F: f64 = 0.25;
            0.5 * (1.0 - ((err.abs() / best_err - O) / F).tanh())
        };

        let mut weights: Vec<Weight> = Vec::with_capacity(self.data.samples.len());
        let mut sw = Vec::new();
        let mut ew = Vec::new();
        let mut daw = Vec::new();
        let mut se = Vec::new();
        let mut ee = Vec::new();
        let mut dae = Vec::new();
        for s in &self.data.samples {
            let eval = self.evaluate_traj(s, traj);
            let b1 = eval.start.map_or(0.0, get_weight);
            let b2 = eval.end.map_or(0.0, get_weight);
            let b3 = eval.da.map_or(0.0, get_weight);
            weights.push(Weight {
                start: b1,
                end: b2,
                da: b3,
            });
            eval.start.map(|x| se.push(x.abs() / best_err));
            eval.end.map(|x| ee.push(x.abs() / best_err));
            eval.da.map(|x| dae.push(x.abs() / best_err));
            eval.start.map(|_| sw.push(b1));
            eval.end.map(|_| ew.push(b2));
            eval.da.map(|_| daw.push(b3));
        }

        println!("Error distribution");
        println!();
        if !sw.is_empty() {
            println!("Start ({}):", se.len());
            draw_hitogram(&se, 5);
        }
        if !ew.is_empty() {
            println!("End ({}):", ee.len());
            draw_hitogram(&ee, 5);
        }
        if !daw.is_empty() {
            println!("DA ({}):", dae.len());
            draw_hitogram(&dae, 5);
        }
        println!();

        println!("Weight distribution");
        println!();
        if !sw.is_empty() {
            println!("Start ({}):", sw.len());
            draw_hitogram(&sw, 5);
        }
        if !ew.is_empty() {
            println!("End ({}):", ew.len());
            draw_hitogram(&ew, 5);
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
                    let eval = self.evaluate_traj(s, traj);
                    if let Some(x) = eval.start {
                        sum += x * x * w.start;
                        cnt += w.start;
                    }
                    if let Some(x) = eval.end {
                        sum += x * x * w.end;
                        cnt += w.end;
                    }
                    if let Some(x) = eval.da {
                        sum += x * x * w.da;
                        cnt += w.da;
                    }
                }
            }
            None => {
                for s in self.data.iter() {
                    let eval = self.evaluate_traj(s, traj);
                    if let Some(x) = eval.start {
                        sum += x * x;
                        cnt += 1.0;
                    }
                    if let Some(x) = eval.end {
                        sum += x * x;
                        cnt += 1.0;
                    }
                    if let Some(x) = eval.da {
                        sum += x * x;
                        cnt += 1.0;
                    }
                }
            }
        }

        sum / cnt
    }

    /// Compule the value and gradient of `eval_mean_of_squares`
    pub fn eval_grad_mean_of_squares(
        &self,
        traj: Line,
        weights: Option<&[Weight]>,
    ) -> (f64, LineGrad) {
        let mut grad_sum = LineGrad::default();
        let mut val_sum = 0.0;
        let mut cnt = 0.0;

        match weights {
            Some(weights) => {
                for (s, w) in self.data.iter().zip(weights) {
                    let (eval, grad) = self.evaluate_grad_traj(s, traj);
                    if let (Some(x), Some(dx)) = (eval.start, grad.start) {
                        grad_sum += dx * 2.0 * x * w.start;
                        val_sum += x * x * w.start;
                        cnt += w.start;
                    }
                    if let (Some(x), Some(dx)) = (eval.end, grad.end) {
                        grad_sum += dx * 2.0 * x * w.end;
                        val_sum += x * x * w.end;
                        cnt += w.end;
                    }
                    if let (Some(x), Some(dx)) = (eval.da, grad.da) {
                        grad_sum += dx * 2.0 * x * w.da;
                        val_sum += x * x * w.da;
                        cnt += w.da;
                    }
                }
            }
            None => {
                for s in self.data.iter() {
                    let (eval, grad) = self.evaluate_grad_traj(s, traj);
                    if let (Some(x), Some(dx)) = (eval.start, grad.start) {
                        grad_sum += dx * 2.0 * x;
                        val_sum += x * x;
                        cnt += 1.0;
                    }
                    if let (Some(x), Some(dx)) = (eval.end, grad.end) {
                        grad_sum += dx * 2.0 * x;
                        val_sum += x * x;
                        cnt += 1.0;
                    }
                    if let (Some(x), Some(dx)) = (eval.da, grad.da) {
                        grad_sum += dx * 2.0 * x;
                        val_sum += x * x;
                        cnt += 1.0;
                    }
                }
            }
        }

        (val_sum / cnt, grad_sum / cnt)
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
        const DA_D: f64 = to_radians(0.5);
        const AZ_D: f64 = to_radians(0.5);

        let mut this = self.clone();

        let eval = this.evaluate_traj(&this.data.samples[sample_i], traj);
        let mut s = Sigmas::default();

        if let Some(da_err) = eval.da {
            let old_da = this.data.samples[sample_i].da;
            *this.data.samples[sample_i].da.as_mut().unwrap() += DA_D;

            let (new_traj, gd_i) =
                this.gradient_descent_complete(traj, weights, GD_TARGET_ERR_DIFF);
            let (new_flash, new_speed) = this.calc_flash_and_speed(new_traj, weights);

            println!(
                "{sample_i}  [da] ({gd_i} gd runs): {:.0}m/{DEGREE_SYM} * {:.0}{DEGREE_SYM} = {:.3}km\n        {:.3}{DEGREE_SYM}/{DEGREE_SYM} * {:.0}{DEGREE_SYM} = {}{DEGREE_SYM}",
                (new_flash - flash).norm() / DA_D.to_degrees(),
                da_err.to_degrees().abs(),
                (new_flash - flash).norm() / DA_D * da_err.abs() * 1e-3,
                (new_traj.direction.dot(*traj.direction).min(1.0).acos() / DA_D).to_degrees(),
                da_err.to_degrees().abs(),
                (new_traj.direction.dot(*traj.direction).min(1.0).acos() * da_err.abs() / DA_D).to_degrees(),
            );

            let new_dir: Azimuthal = new_traj.direction.into_inner().into();
            this.data.samples[sample_i].da = old_da;

            let mul = da_err * da_err / DA_D / DA_D;
            s.x += mul * f64::powi(new_flash.x - flash.x, 2);
            s.y += mul * f64::powi(new_flash.y - flash.y, 2);
            s.z += mul * f64::powi(new_flash.z - flash.z, 2);
            s.v_z += mul * f64::powi(new_dir.z - dir.z, 2);
            s.v_h += mul * f64::powi(new_dir.h - dir.h, 2);
            s.v_angle +=
                mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
            s.speed += mul * f64::powi(new_speed - speed, 2);
        }

        if let Some(z0_err) = eval.end {
            let old_z0 = this.data.samples[sample_i].z0;
            *this.data.samples[sample_i].z0.as_mut().unwrap() += AZ_D;

            let (new_traj, gd_i) =
                this.gradient_descent_complete(traj, weights, GD_TARGET_ERR_DIFF);
            let (new_flash, new_speed) = this.calc_flash_and_speed(new_traj, weights);
            println!(
                "{sample_i}   [z0] ({gd_i} gd runs): {:.0}m/{DEGREE_SYM} * {:.0}{DEGREE_SYM} = {:.3}km\n        {:.3}{DEGREE_SYM}/{DEGREE_SYM} * {:.0}{DEGREE_SYM} = {}{DEGREE_SYM}",
                (new_flash - flash).norm() / DA_D.to_degrees(),
                z0_err.to_degrees().abs(),
                (new_flash - flash).norm() / DA_D * z0_err.abs() * 1e-3,
                (new_traj.direction.dot(*traj.direction).min(1.0).acos() / DA_D).to_degrees(),
                z0_err.to_degrees().abs(),
                (new_traj.direction.dot(*traj.direction).min(1.0).acos() * z0_err.abs() / DA_D).to_degrees(),
            );

            let new_dir: Azimuthal = new_traj.direction.into_inner().into();
            this.data.samples[sample_i].z0 = old_z0;

            let mul = z0_err * z0_err / AZ_D / AZ_D;
            s.x += mul * f64::powi(new_flash.x - flash.x, 2);
            s.y += mul * f64::powi(new_flash.y - flash.y, 2);
            s.z += mul * f64::powi(new_flash.z - flash.z, 2);
            s.v_z += mul * f64::powi(new_dir.z - dir.z, 2);
            s.v_h += mul * f64::powi(new_dir.h - dir.h, 2);
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
        F: FnOnce(&mut Self) -> (String, Line),
    {
        let timer = Instant::now();
        let (title, result) = f(self);
        let duration = timer.elapsed();

        self.data.compare(result, &title);
        eprintln!("This step took {duration:?}\n");

        result
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        let traj = self.run_stage(|this| ("Initial guess (pairwise)".into(), this.pairwise()));

        let traj = self.run_stage(|this| {
            let (traj, gd_i) = this.gradient_descent_complete(traj, None, GD_FIRST_ERR_DIFF);
            (format!("Gradient descent ({gd_i} runs)"), traj)
        });

        // {
        // let mut errs = Vec::new();
        // let mut errs_sq = Vec::new();
        // for s in &self.data.samples {
        //     let eval = self.evaluate_traj(s, traj);
        //     if let Some(da) = eval.da {
        //         let da = da.to_degrees();
        //         errs.push(da);
        //         errs_sq.push(da * da);
        //     }
        // }
        // draw_hitogram(&errs, 5);
        // println!();
        // draw_hitogram(&errs_sq, 5);
        // println!();

        // let mut list = Vec::new();
        // for (i, s) in self.data.samples.iter().enumerate() {
        //     let eval = self.evaluate_traj(s, traj);
        //     if let Some(da) = eval.da {
        //         list.push((da.to_degrees(), s.name.as_deref().unwrap_or("N/A"), i));
        //     }
        // }
        // list.sort_unstable_by(|a, b| a.0.abs().total_cmp(&b.0.abs()));
        // for entry in &list {
        //     dbg!(entry);
        // }
        // eprintln!("median: {:#?}", list[list.len() / 2]);
        // }

        // {
        //     let mut points = Vec::new();
        //     let mut p = Geodetic::from_geocentric_cartesian(traj.point, 10);
        //     p.h = 0.0;
        //     while p.h < 6_000_000.0 {
        //         let mut traj = traj;
        //         traj.point = p.into_geocentric_cartesian();
        //         points.push((
        //             p.h * 1e-3,
        //             self.eval_mean_of_squares(traj, None),
        //             BLACK,
        //             1.0,
        //         ));
        //         p.h += 10_000.0;
        //     }
        //     draw_plot_svg_with_named_axes(
        //         "plots/height-err-after.svg",
        //         &points,
        //         &[],
        //         &[],
        //         &[],
        //         "height",
        //         "error",
        //     )
        //     .unwrap();
        // }

        let traj = self.lms(traj);
        self.data.compare(traj, "After LMS");

        let mut weights = self.compute_weights(traj);

        let mut traj = self.ls(traj, Some(&weights));
        self.data.compare(traj, "After LS");

        if !self.params.no_da_flip {
            eprintln!("flipped {} descent angles", self.flip_da(traj));
            eprintln!();
            weights = self.compute_weights(traj);
            traj = self.ls(traj, Some(&weights));
            self.data.compare(traj, "LS (after DA-flip)");
        }

        // self.draw_2d_image("after-ls", traj, 100_000.0, 120.0, |traj| {
        //     self.eval_mean_of_squares(traj, Some(&weights))
        // });

        let (traj, gd_i) = self.gradient_descent_complete(traj, Some(&weights), GD_TARGET_ERR_DIFF);
        self.data
            .compare(traj, &format!("Gradient descent ({} runs)", gd_i));

        let sigmas = dbg!(self.calc_sigmas(traj, Some(&weights)));
        dbg!(Vec3::new(sigmas.x, sigmas.y, sigmas.z).norm() / 1000.0);
        dbg!(sigmas.v_angle.to_degrees());

        self.data
            .compare(traj, &format!("Gradient descent ({} runs)", gd_i));

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

    const ZENITH_COEF_ANGLE: f64 = to_radians(40.0);
    const ZENITH_COEF_SIN_K: f64 = FRAC_PI_2 / Self::ZENITH_COEF_ANGLE;
    const ZENITH_COEF_K: f64 = 5.0;

    fn calc_zenith_coef(sample: &DataSample, ek: UnitVec3) -> f64 {
        // k(x) = 1/sin(x * ZENITH_COEF_SIN_K)
        // k'(x) = -1 * a * ZENITH_COEF_SIN_K / (1 - a^2), where a = cos(x * ZENITH_COEF_SIN_K)
        let to_zenith = sample.zenith_dir.dot(*ek).min(1.0).acos();
        if to_zenith < Self::ZENITH_COEF_ANGLE {
            Self::ZENITH_COEF_K / f64::sin(to_zenith * Self::ZENITH_COEF_SIN_K)
        } else if to_zenith > PI - Self::ZENITH_COEF_ANGLE {
            Self::ZENITH_COEF_K / f64::sin((PI - to_zenith) * Self::ZENITH_COEF_SIN_K)
        } else {
            1.0
        }
        // 1.0
    }

    fn calc_zenith_coef_diff(sample: &DataSample, ek: UnitVec3, ek_diff: Vec3) -> f64 {
        // k(x) = 1/sin(x * ZENITH_COEF_SIN_K)
        // k'(x) = -1 * a * ZENITH_COEF_SIN_K / (1 - a^2), where a = cos(x * ZENITH_COEF_SIN_K)

        let to_zenith_cos = sample.zenith_dir.dot(*ek);
        let to_zenith_cos_diff = sample.zenith_dir.dot(ek_diff);

        let to_zenith = to_zenith_cos.min(1.0).acos();
        let to_zenith_diff =
            -1.0 / f64::sqrt(1.0 - to_zenith_cos * to_zenith_cos) * to_zenith_cos_diff;

        if to_zenith < Self::ZENITH_COEF_ANGLE {
            let a = f64::cos(to_zenith * Self::ZENITH_COEF_SIN_K);
            -1.0 * Self::ZENITH_COEF_K * a * Self::ZENITH_COEF_SIN_K * to_zenith_diff
                / (1.0 - a * a)
        } else if to_zenith > PI - Self::ZENITH_COEF_ANGLE {
            let a = f64::cos((PI - to_zenith) * Self::ZENITH_COEF_SIN_K);
            a * Self::ZENITH_COEF_K * Self::ZENITH_COEF_SIN_K * to_zenith_diff / (1.0 - a * a)
        } else {
            0.0
        }

        // 0.0
    }

    pub fn evaluate_traj(&self, sample: &DataSample, traj: Line) -> Evaluation {
        let ek = UnitVec3::new_normalize(traj.point - sample.location);
        let zenith_coef = Self::calc_zenith_coef(sample, ek);

        Evaluation {
            start: None,
            end: if !self.params.da_only {
                sample
                    .z0
                    .map(|z0| angle_diff(z0, sample.calc_azimuth(traj.point)))
            } else {
                None
            },
            da: sample.da.map(|da| {
                angle_diff(descent_angle(sample.zenith_dir, ek, *traj.direction), da) * zenith_coef
            }),
        }
    }

    /// Compute the value and gradient of `evaluate_traj`
    pub fn evaluate_grad_traj(
        &self,
        s: &DataSample,
        traj: Line,
    ) -> (Evaluation, Evaluation<LineGrad>) {
        let err0 = self.evaluate_traj(s, traj);

        let k = traj.point - s.location;
        let (ek, k_norm) = UnitVec3::new_and_get(k);
        let zenith_coef = Self::calc_zenith_coef(s, ek);

        fn make_end_diff(s: &DataSample, k: Vec3, mask: Vec3) -> f64 {
            // E = diff(z0, atan2(k*east, k*north))
            // E = atan2(k*east, k*north) - z0
            // E,x = -(k*east) / (xx+yy) * north_x + (k*north) / (xx+yy) * east_x
            // E,x = (k*north * east_x - k*east * north_x) / (xx+yy)
            let x = k.dot(*s.east_dir);
            let y = k.dot(*s.north_dir);
            (y * s.east_dir.dot(mask) - x * s.north_dir.dot(mask)) / (x * x + y * y)
        }

        let make_da_diff_point = |s: &DataSample, da: f64, vel: Vec3, mask: Vec3| -> f64 {
            // k_norm = (k_norm^2)^1/2
            // k_norm,x = 1/2 * (k_norm^2)^(-1/2) * 2kx = kx / k_norm
            let diff_k_norm = k.dot(mask) / k_norm;

            let diff_k = mask;

            // ek = k * k_norm^-1
            // ek,x = (k,x * k_norm^-1) + (k * -1 * k_norm^-2 * k_norm,x)
            let diff_ek = (diff_k / k_norm) - (k * diff_k_norm / k_norm / k_norm);

            let diff_zenith_coef = Self::calc_zenith_coef_diff(s, ek, diff_ek);

            let xa = ek.cross(*s.zenith_dir);
            let diff_xa = diff_ek.cross(*s.zenith_dir);

            let ya = xa.cross(*ek);
            let diff_ya = diff_xa.cross(*ek) + xa.cross(diff_ek);

            let x = vel.dot(xa);
            let diff_x = vel.dot(diff_xa);

            let y = vel.dot(ya);
            let diff_y = vel.dot(diff_ya);

            // let x = k.cross(*zenith); // -> +0
            // let y = x.cross(*k); // -> +0
            let dac = f64::atan2(x, y);
            let diff_dac = (y * diff_x - x * diff_y) / (x * x + y * y);

            // Err = (da - dac) * zenith_coef = da*zenith_coef - dac*zenith_coef
            // Err,i = da*zenith_coef,i - dac*zenith_coef,i - dac,i*zenith_coef
            //       = zenith_coef,i*(da - dac) - dac,i*zenith_coef
            diff_zenith_coef * angle_diff(dac, da) - diff_dac * zenith_coef
        };

        let make_da_diff_dir = |s: &DataSample, ek: UnitVec3, vel: Vec3, mask: Vec3| -> f64 {
            let xa = ek.cross(*s.zenith_dir);
            let ya = xa.cross(*ek);

            // v = {x, y, z}
            // dx
            // v' = {x+dx, y, z} / (1 + dx*dx + 2x*dx)
            // dv = ...
            // dv/dx = {1-xx, -xy, -xz}
            let diff_vel = mask - vel * mask.dot(vel);

            let x = vel.dot(xa);
            let diff_x = diff_vel.dot(xa);

            let y = vel.dot(ya);
            let diff_y = diff_vel.dot(ya);

            (x * diff_y - y * diff_x) / (x * x + y * y) * zenith_coef
        };

        let end = (s.z0.is_some() && !self.params.da_only).then(|| LineGrad {
            point: Vec3::new(
                make_end_diff(s, k, Vec3::x()),
                make_end_diff(s, k, Vec3::y()),
                make_end_diff(s, k, Vec3::z()),
            ),
            direction: Vec3::default(),
        });

        let da = s.da.map(|da| LineGrad {
            point: Vec3::new(
                make_da_diff_point(s, da, *traj.direction, Vec3::x()),
                make_da_diff_point(s, da, *traj.direction, Vec3::y()),
                make_da_diff_point(s, da, *traj.direction, Vec3::z()),
            ),
            direction: Vec3::new(
                make_da_diff_dir(s, ek, *traj.direction, Vec3::x()),
                make_da_diff_dir(s, ek, *traj.direction, Vec3::y()),
                make_da_diff_dir(s, ek, *traj.direction, Vec3::z()),
            ),
        });

        (
            err0,
            Evaluation {
                start: None,
                end,
                da,
            },
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Evaluation<T = f64> {
    pub start: Option<T>,
    pub end: Option<T>,
    pub da: Option<T>,
}

impl ops::Sub for Evaluation {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            start: self.start.zip(rhs.start).map(|(x, y)| x - y),
            end: self.end.zip(rhs.end).map(|(x, y)| x - y),
            da: self.da.zip(rhs.da).map(|(x, y)| x - y),
        }
    }
}

impl ops::Div<f64> for Evaluation {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            start: self.start.map(|x| x / rhs),
            end: self.end.map(|x| x / rhs),
            da: self.da.map(|x| x / rhs),
        }
    }
}

impl Evaluation {
    fn iter_squared(self) -> impl Iterator<Item = f64> {
        let sq = |x| x * x;
        self.start
            .map(sq)
            .into_iter()
            .chain(self.end.map(sq))
            .chain(self.da.map(sq))
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Weight {
    pub start: f64,
    pub end: f64,
    pub da: f64,
}

fn lerp(a: f64, b: f64, w: f64) -> f64 {
    a + (b - a) * w
}

const fn sq(x: f64) -> f64 {
    x * x
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Sigmas {
    /// Meters
    pub x: f64,
    /// Meters
    pub y: f64,
    /// Meters
    pub z: f64,
    /// Radians
    pub v_z: f64,
    /// Radians
    pub v_h: f64,
    /// Radians
    pub v_angle: f64,
    /// Meters/Seconds
    pub speed: f64,
}

impl Sigmas {
    fn sqrt(self) -> Self {
        Self {
            x: self.x.sqrt(),
            y: self.y.sqrt(),
            z: self.z.sqrt(),
            v_z: self.v_z.sqrt(),
            v_h: self.v_h.sqrt(),
            v_angle: self.v_angle.sqrt(),
            speed: self.speed.sqrt(),
        }
    }
}

impl ops::AddAssign for Sigmas {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
        self.v_z += rhs.v_z;
        self.v_h += rhs.v_h;
        self.v_angle += rhs.v_angle;
        self.speed += rhs.speed;
    }
}

impl ops::Add for Sigmas {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl std::iter::Sum for Sigmas {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |s, x| s + x)
    }
}
