//! The implementation

use std::f64::consts::*;
use std::intrinsics::unlikely;
use std::ops;

use common::constants::*;
use common::histogram::draw_hitogram;
use common::maths::*;
use common::obs_data::{Data, DataSample};
use common::plot::{draw_plot_svg, plotters, weight_to_rgb};
use common::quick_median::SliceExt;
use common::rand::random;
use common::structs::*;
use common::{rand, rand_distr};

use image::{Rgb, RgbImage};

use plotters::style::colors;
use rand::prelude::SmallRng;
use rand::{Rng, SeedableRng};
use rand_distr::StandardNormal;

pub mod pair;
use pair::PairTrajectory;

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
                "  direction adjusted by {:.1}{DEGREE} ({:.1}% error)",
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

    pub fn gradient_descent(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
        iterations: usize,
    ) -> Line {
        let diff_dir = |err0: f64, dir_a: Azimuthal, point: Vec3| {
            const D_V: f64 = radians(1e-6);

            let eval =
                |direction: UnitVec3| self.eval_mean_of_squares(Line { point, direction }, weights);

            let dir_z = UnitVec3::from(Azimuthal {
                z: dir_a.z + D_V,
                h: dir_a.h,
            });

            let dir_h = UnitVec3::from(Azimuthal {
                z: dir_a.z,
                h: dir_a.h + D_V,
            });

            let diff_vz = (eval(dir_z) - err0) / D_V;
            let diff_vh = (eval(dir_h) - err0) / D_V;

            (diff_vz, diff_vh)
        };
        let diff_p = |err0: f64, traj: Line| {
            const D_P: f64 = 1e-6;

            let eval = |direction: UnitVec3, point: Vec3| {
                self.eval_mean_of_squares(Line { point, direction }, weights)
            };
            let eval_with_offset = |offset: Vec3| eval(traj.direction, traj.point + offset);

            let diff_x = (eval_with_offset(Vec3::x() * D_P) - err0) / D_P;
            let diff_y = (eval_with_offset(Vec3::y() * D_P) - err0) / D_P;
            let diff_z = (eval_with_offset(Vec3::z() * D_P) - err0) / D_P;

            Vec3::new(diff_x, diff_y, diff_z)
        };

        let mut point_descent_coeff: f64 = 1e10;
        let mut direction_descent_coeff: f64 = 1e10;

        let mut cur_eval = self.eval_mean_of_squares(traj, weights);
        let mut cur_dir_a = Azimuthal::from(traj.direction);

        for _ in 0..iterations {
            if point_descent_coeff > 1e-100 {
                // for dp_exp in (-10..=4).rev() {
                //     let dp = 10f64.powi(dp_exp);
                //     eprintln!(
                //         "10^{dp_exp:+03} -> 10^{}",
                //         diff_p(cur_eval, traj, dp).norm().log10()
                //     );
                // }
                // eprintln!(
                //     "2.329*10^-10 -> 10^{}",
                //     diff_p(cur_eval, traj, 2.329e-10).norm().log10()
                // );
                // eprintln!();
                // let diff = diff_p(cur_eval, traj, 0.0001);
                let diff = diff_p(cur_eval, traj);
                let old_eval = cur_eval;
                let old = traj;

                traj.point -= diff * point_descent_coeff;
                cur_eval = self.eval_mean_of_squares(traj, weights);

                if cur_eval >= old_eval {
                    traj = old;
                    cur_eval = old_eval;
                    point_descent_coeff *= 0.5;
                }
            }

            //----------------------------------------------------

            if direction_descent_coeff > 1e-100 {
                let (vz, vh) = diff_dir(cur_eval, cur_dir_a, traj.point);
                let old_dir_a = cur_dir_a;
                let old_eval = cur_eval;
                let old_traj = traj;

                cur_dir_a.z -= vz * direction_descent_coeff;
                cur_dir_a.h -= vh * direction_descent_coeff;
                traj.direction = UnitVec3::from(cur_dir_a);
                cur_eval = self.eval_mean_of_squares(traj, weights);

                if cur_eval >= old_eval {
                    cur_dir_a = old_dir_a;
                    traj = old_traj;
                    cur_eval = old_eval;
                    direction_descent_coeff *= 0.5;
                }
            }

            //----------------------------------------------------

            if point_descent_coeff < 1e-100 && direction_descent_coeff < 1e-100 {
                break;
            }
        }

        traj
    }

    pub fn gradient_descent_complete(
        &self,
        mut traj: Line,
        weights: Option<&[Weight]>,
    ) -> (Line, usize) {
        const ITERS: usize = 10_000;
        const MIN_D_POINT_SQ: f64 = sq(0.0001);
        const MIN_D_VEL_SQ: f64 = sq(radians(0.001));

        for runs in 1.. {
            let new_traj = self.gradient_descent(traj, weights, ITERS);
            if (new_traj.point - traj.point).norm_squared() <= MIN_D_POINT_SQ
                && (*new_traj.direction - *traj.direction).norm_squared() <= MIN_D_VEL_SQ
            {
                // dbg!(self.eval_mean_of_squares(traj, weights));
                return (new_traj, runs);
            }
            traj = new_traj;
        }

        unreachable!()
    }

    pub fn compute_weights(&self, traj: Line) -> Vec<Weight> {
        let best_err = self.eval_median_of_squares(traj).sqrt()
            * 1.483
            * (1. + 5. / (self.data.samples.len() - 6) as f64);
        let get_weight = |err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const O: f64 = 0.95;
            const F: f64 = 0.36;
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

        match weights {
            Some(weights) => {
                let mut cnt = 0.0;
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
                sum / cnt
            }
            None => {
                let mut cnt = 0usize;
                for s in self.data.iter() {
                    let eval = self.evaluate_traj(s, traj);
                    if let Some(x) = eval.start {
                        sum += x * x;
                        cnt += 1;
                    }
                    if let Some(x) = eval.end {
                        sum += x * x;
                        cnt += 1;
                    }
                    if let Some(x) = eval.da {
                        sum += x * x;
                        cnt += 1;
                    }
                }
                sum / cnt as f64
            }
        }
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

    fn calc_sigmas(&mut self, traj: Line, weights: Option<&[Weight]>) -> Sigmas {
        let (flash, speed) = self.calc_flash_and_speed(traj, weights);
        let dir: Azimuthal = traj.direction.into_inner().into();

        // dbg!(flash);
        // dbg!(dir);

        let mut sigmas = Sigmas::default();

        const DA_D: f64 = radians(0.02);
        const AZ_D: f64 = radians(0.02);

        for i in 0..(self.data.samples.len()) {
            if i != 0 {
                print!("{ANSI_GOTO_PREV_LINE}{ANSI_CLEAR_LINE}");
            }

            println!(
                "Calculating sigmas: {} ({}/{})",
                &"##########----------"[(10 - (i + 1) * 10 / self.data.len())..][..10],
                i + 1,
                self.data.samples.len()
            );
            let eval = self.evaluate_traj(&self.data.samples[i], traj);

            if let Some(da_err) = eval.da {
                let old_da = self.data.samples[i].da;
                *self.data.samples[i].da.as_mut().unwrap() += DA_D;

                let (new_traj, _gd_i) = self.gradient_descent_complete(traj, weights);
                let (new_flash, new_speed) = self.calc_flash_and_speed(new_traj, weights);
                // println!(
                //     "  [da] ({}K): {:.0}m/{DEGREE} * {:.0}{DEGREE} = {:.3}km\n        {:.3}{DEGREE}/{DEGREE} * {:.0}{DEGREE} = {}",
                //     gd_i / 1000,
                //     (new_flash - flash).norm() / DA_D.to_degrees(),
                //     da_err.to_degrees().abs(),
                //     (new_flash - flash).norm() / DA_D * da_err.abs() * 1e-3,
                //     new_traj.direction.dot(*traj.direction).min(1.0).acos().to_degrees() / DA_D,
                //     da_err.to_degrees().abs(),
                //     new_traj.direction.dot(*traj.direction).min(1.0).acos() * da_err.abs() / DA_D,
                // );

                let new_dir: Azimuthal = new_traj.direction.into_inner().into();
                self.data.samples[i].da = old_da;

                let mul = da_err * da_err / DA_D / DA_D;
                sigmas.x += mul * f64::powi(new_flash.x - flash.x, 2);
                sigmas.y += mul * f64::powi(new_flash.y - flash.y, 2);
                sigmas.z += mul * f64::powi(new_flash.z - flash.z, 2);
                sigmas.v_z += mul * f64::powi(new_dir.z - dir.z, 2);
                sigmas.v_h += mul * f64::powi(new_dir.h - dir.h, 2);
                sigmas.v_angle +=
                    mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
                sigmas.speed += mul * f64::powi(new_speed - speed, 2);
            }

            // if let Some(z0_err) = eval.end {
            //     let old_z0 = self.data.samples[i].z0;
            //     *self.data.samples[i].z0.as_mut().unwrap() += AZ_D;

            //     let (new_traj, _gd_i) = self.gradient_descent_complete(traj, weights);
            //     let (new_flash, new_speed) = self.calc_flash_and_speed(new_traj, weights);
            //     // println!(
            //     //     "    [z0] ({}K): {:.0}m/{DEGREE} * {:.0}{DEGREE} = {:.3}km",
            //     //     gd_i / 1000,
            //     //     (new_flash - flash).norm() / DA_D.to_degrees(),
            //     //     z0_err.to_degrees().abs(),
            //     //     (new_flash - flash).norm() / DA_D * z0_err.abs() * 1e-3,
            //     // );

            //     let new_dir: Azimuthal = new_traj.direction.into_inner().into();
            //     self.data.samples[i].z0 = old_z0;

            //     let mul = z0_err * z0_err / AZ_D / AZ_D;
            //     sigmas.x += mul * f64::powi(new_flash.x - flash.x, 2);
            //     sigmas.y += mul * f64::powi(new_flash.y - flash.y, 2);
            //     sigmas.z += mul * f64::powi(new_flash.z - flash.z, 2);
            //     sigmas.v_z += mul * f64::powi(new_dir.z - dir.z, 2);
            //     sigmas.v_h += mul * f64::powi(new_dir.h - dir.h, 2);
            //     sigmas.v_angle +=
            //         mul * f64::powi(new_traj.direction.dot(*traj.direction).min(1.0).acos(), 2);
            //     sigmas.speed += mul * f64::powi(new_speed - speed, 2);
            // }
        }

        sigmas.x = sigmas.x.sqrt();
        sigmas.y = sigmas.y.sqrt();
        sigmas.z = sigmas.z.sqrt();
        sigmas.v_z = sigmas.v_z.sqrt();
        sigmas.v_h = sigmas.v_h.sqrt();
        sigmas.v_angle = sigmas.v_angle.sqrt();
        sigmas.speed = sigmas.speed.sqrt();

        sigmas
    }

    fn flip_da(&mut self, traj: Line) -> usize {
        let mut flipped = 0;

        for i in 0..self.data.samples.len() {
            let Some(dda) = self.evaluate_traj(&self.data.samples[i], traj).da
            else { continue };

            if dda.abs() > f64::to_radians(120.0) {
                let da = self.data.samples[i].da.as_mut().unwrap();
                *da = (*da + PI) % TAU;
                flipped += 1;
            }
        }

        flipped
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        let traj = self.pairwise();
        self.data.compare(traj, "Initial guess (pairwise)");

        let (traj, gd_i) = self.gradient_descent_complete(traj, None);
        self.data
            .compare(traj, &format!("Gradient descent ({} runs)", gd_i));

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

        let (traj, gd_i) = self.gradient_descent_complete(traj, Some(&weights));
        self.data
            .compare(traj, &format!("Gradient descent ({} runs)", gd_i));

        // let sigmas = dbg!(self.calc_sigmas(traj, Some(&weights)));
        // dbg!(Vec3::new(sigmas.x, sigmas.y, sigmas.z).norm() / 1000.0);
        // dbg!(sigmas.v_angle.to_degrees());

        self.draw_2d_image("ls", traj, 100_000.0, 120.0, |traj| {
            self.eval_mean_of_squares(traj, Some(&weights))
        });

        let (flash, speed) = self.calc_flash_and_speed(traj, Some(&weights));

        Solution {
            flash: Geodetic::from_geocentric_cartesian(flash, 10),
            velocity: traj.direction.into_inner() * speed,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, traj: Line, weights: Option<&[Weight]>) -> (Vec3, f64) {
        let mut lambdas = Vec::new();
        // let mut lambdas_all = Vec::new();
        let mut speeds = Vec::new();

        // let mut skipped_l_end_w = 0;
        // let mut skipped_l_end_dir = 0;
        // let mut total_l_end = 0;

        for (i, s) in self.data.samples.iter().enumerate() {
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

            if weights.is_some_and(|w| w[i].end > 0.5) {
                lambdas.push(l_end);
                // lambdas_all.push(l_end);
            } else {
                // lambdas_all.push(l_end);
                // skipped_l_end_w += 1;
            }

            if s.observation_matches(traj) {
                if let Some((duration, k_start)) = s.dur.zip(s.k_start) {
                    let l_start = lambda(s.location, k_start, traj.point, traj.direction);
                    let dist = l_end - l_start;
                    speeds.push(dist / duration);
                }
            }
        }

        // eprintln!("skipped {skipped_l_end_w}(weight) + {skipped_l_end_dir}(direction) l_end out of {total_l_end}");

        let l_end = lambdas.median();
        // let l_end_all = lambdas_all.median();

        // dbg!(l_end - l_end_all);

        let speed = speeds.median();
        (traj.point + traj.direction.into_inner() * l_end, speed)
    }

    pub fn evaluate_traj(&self, sample: &DataSample, traj: Line) -> Evaluation {
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
                angle_diff(
                    descent_angle(
                        sample.location,
                        traj.point - sample.location,
                        traj.direction.into_inner(),
                    ),
                    da,
                )
            }),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Evaluation {
    pub start: Option<f64>,
    pub end: Option<f64>,
    pub da: Option<f64>,
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

pub struct Weight {
    pub start: f64,
    pub end: f64,
    pub da: f64,
}

fn lerp(a: f64, b: f64, w: f64) -> f64 {
    a + (b - a) * w
}

#[derive(Debug, Default)]
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

const fn sq(x: f64) -> f64 {
    x * x
}
