//! The implementation

use std::f64::consts::*;
use std::intrinsics::unlikely;

use common::constants::EARTH_R;
use common::histogram::draw_hitogram;
use common::maths::*;
use common::obs_data::{Data, DataSample};
use common::plot::{draw_plot_svg, plotters, weight_to_rgb};
use common::quick_median::SliceExt;
use common::rand::random;
use common::structs::*;
use common::{nalgebra, rand, rand_distr};

use image::{Rgb, RgbImage};

use nalgebra::Unit;
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
    /// Whether to tilt the observation to better match given descent angle
    pub da_correction: f64,
    /// Max correction that DAC is allowed to apply
    pub dac_max: f64,
    /// `evaluate_traj` uses only descent angles
    pub da_only: bool,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub flash: Vec3,
    pub velocity: Vec3,
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver { data, params }
    }

    /// We assume that k_start & k_end are less accurate than descent_angle, so we tweak them to
    /// match with descent_angle
    pub fn enhance_data(&mut self) {
        if self.params.da_correction > 0.0 {
            for s in &mut self.data.samples {
                if let Some(((da, axis), comp_da)) = s.da.zip(s.axis).zip(s.calculated_da()) {
                    let angle = angle_diff(comp_da, da) * self.params.da_correction;
                    if angle.abs() <= self.params.dac_max {
                        let q = UnitQuaternion::from_axis_angle(&axis, angle);
                        s.plane.as_mut().map(|k| *k = q * *k);
                        s.k_start.as_mut().map(|k| *k = q * *k);
                        s.k_end.as_mut().map(|k| *k = q * *k);
                    }
                }
            }
        }
    }

    /// The weighted sum of pairwise plane crossings
    pub fn pairwise(&self) -> Line {
        let mut point = Vec3::default();
        let mut dir = Vec3::default();
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
            direction: Unit::new_normalize(dir),
        }
    }

    pub fn draw_2d_image<F: Fn(Line) -> f64>(
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

        let initial_spherical: Spherical = traj.point.into();
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

            fn get(&mut self, x: u32, y: u32) -> &mut f64 {
                &mut self.mem[x as usize + y as usize * W]
            }
        }

        let mut err_map = Flat2D::<1000, 1000>::new();

        for x in 0..1000 {
            let lon = lon + lerp(-d_lon, d_lon, x as f64 / 1000.0);
            for y in 0..1000 {
                let lat = lat + lerp(d_lat, -d_lat, y as f64 / 1000.0);

                let sph = Spherical {
                    lat,
                    lon,
                    ..initial_spherical
                };
                let mut copy = traj;
                copy.point = sph.into();

                let error = eval_fn(copy);
                *err_map.get(x, y) = error;

                min = min.min(error);
                max = max.max(error);
            }
        }

        let mut img = RgbImage::new(1000, 1000);

        const WHITE: Rgb<u8> = Rgb([u8::MAX, u8::MAX, u8::MAX]);
        const BLUE: Rgb<u8> = Rgb([0, 0, u8::MAX]);
        const BLACK: Rgb<u8> = Rgb([0, 0, 0]);

        for x in 0..1000 {
            for y in 0..1000 {
                let e = *err_map.get(x, y);
                let color = if unlikely(e == min) {
                    WHITE
                } else if unlikely(e - min < min * 0.01) {
                    BLUE
                } else {
                    error_to_color(e - min, contrast_k)
                };
                img.put_pixel(x, y, color);
            }
        }

        // Fill a pixel for each observer
        for s in &self.data.samples {
            let x = ((s.geo_location.lon - initial_spherical.lon) / d_lon * 500.0 + 500.0).round()
                as i32;
            let y = ((s.geo_location.lat - initial_spherical.lat) / d_lat * -500.0 + 500.0).round()
                as i32;
            let mut put_pixel = |x: i32, y: i32| {
                if (0..1000).contains(&x) && (0..1000).contains(&y) {
                    img.put_pixel(x as u32, y as u32, BLACK);
                }
            };
            put_pixel(x, y);
            put_pixel(x, y + 1);
            put_pixel(x + 1, y);
            put_pixel(x + 1, y + 1);
        }

        let path = format!("plots/{name}-2d-{}.png", random::<u32>());
        eprintln!("saving to {path}");
        img.save(path).unwrap();
    }

    fn generic_iterative_search<F: Fn(Line) -> f64 + Send + Copy>(
        mut best: Line,
        dir_iters: usize,
        dir_sigma: f64,
        pos_iters: usize,
        pos_sigma: f64,
        eval_fn: F,
    ) -> Line {
        let threads_cnt = std::thread::available_parallelism().unwrap().get();

        let mut best_error = eval_fn(best);
        let mut last_point = best.point;
        eprintln!("Initial point: {}", Spherical::from(last_point));
        eprintln!("initial best_error: {best_error}");

        for iter_i in 1.. {
            eprintln!("Iteration #{iter_i}");
            let prev_iter_err = best_error;

            // Adjust direction
            std::thread::scope(|s| {
                let handles: Vec<_> = (0..threads_cnt)
                    .map(|_| {
                        s.spawn(move || {
                            let mut rng = SmallRng::from_entropy();
                            let mut local_best = best;
                            let mut local_best_error = best_error;

                            for _ in 0..(dir_iters / threads_cnt) {
                                let mut copy = best;
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
                    .collect();

                for handle in handles {
                    let (t, te) = handle.join().unwrap();
                    if te < best_error {
                        best_error = te;
                        best = t;
                    }
                }
            });
            eprintln!("best_error after direction adjustments: {best_error}");

            let initial_spherical: Spherical = best.point.into();
            let d_lat = pos_sigma / (EARTH_R * initial_spherical.lat.cos());
            let d_lon = pos_sigma / EARTH_R;

            // Adjust position
            std::thread::scope(|s| {
                let handles: Vec<_> = (0..threads_cnt)
                    .map(|_| {
                        s.spawn(move || {
                            let mut rng = SmallRng::from_entropy();
                            let mut local_best = best;
                            let mut local_best_error = best_error;

                            for _ in 0..(pos_iters / threads_cnt) {
                                const R_MIN: f64 = EARTH_R + 40_000.0;
                                let r_mid = lerp(R_MIN, initial_spherical.r, 0.4);
                                let dr = (r_mid - R_MIN).abs() / 2.0 + 5_000.0;

                                let mut spherical = initial_spherical;
                                spherical.lat += d_lat * rng.sample::<f64, _>(StandardNormal);
                                spherical.lon += d_lon * rng.sample::<f64, _>(StandardNormal);
                                spherical.r = r_mid + dr * rng.sample::<f64, _>(StandardNormal);

                                let mut copy = best;
                                copy.point = spherical.into();

                                let error = eval_fn(copy);
                                if error < local_best_error {
                                    local_best_error = error;
                                    local_best = copy;
                                }
                            }

                            (local_best, local_best_error)
                        })
                    })
                    .collect();

                for handle in handles {
                    let (t, te) = handle.join().unwrap();
                    if te < best_error {
                        best_error = te;
                        best = t;
                    }
                }
            });

            if last_point == best.point {
                eprintln!("Point: N/C");
            } else {
                last_point = best.point;
                eprintln!("Point: {}", Spherical::from(last_point));
            }
            eprintln!("best_error after position adjustments: {best_error}");

            if (prev_iter_err - best_error) / best_error < 0.001 {
                break;
            }
        }

        best
    }

    /// Perform the "least-median-of-squares"
    pub fn lms(&self, best: Line) -> Line {
        self.draw_2d_image("lms", best, 600_000.0, 60.0, |traj| {
            self.eval_median_of_squares(traj)
        });

        Self::generic_iterative_search(
            best,
            1_000_000,
            40f64.to_radians(),
            15_000_000,
            300_000.0,
            |traj| self.eval_median_of_squares(traj),
        )
    }

    /// Perform the "weidghted least-squares"
    pub fn weighted_ls(&self, traj: Line, weights: &[(f64, f64, f64)]) -> Line {
        self.draw_2d_image("ls", traj, 500_000.0, 100.0, |traj| {
            self.eval_mean_of_squares(traj, weights)
        });

        Self::generic_iterative_search(
            traj,
            1_000_000,
            20f64.to_radians(),
            15_000_000,
            100_000.0,
            |traj| self.eval_mean_of_squares(traj, weights),
        )
    }

    pub fn gradient_descent(
        &self,
        mut traj: Line,
        weights: &[(f64, f64, f64)],
        iterations: usize,
    ) -> Line {
        let diff = |point: Vec3, dir: Azimuthal| {
            let d_p: f64 = 50.0;
            let d_v: f64 = 0.05f64.to_radians();

            let evall_ls = |dir: Azimuthal, point| {
                self.eval_mean_of_squares(
                    Line {
                        point,
                        direction: dir.into(),
                    },
                    weights,
                )
            };

            let err0 = evall_ls(dir, point).sqrt();

            let diff_x = (evall_ls(dir, point + Vec3::x() * d_p).sqrt() - err0) / d_p;
            let diff_y = (evall_ls(dir, point + Vec3::y() * d_p).sqrt() - err0) / d_p;
            let diff_z = (evall_ls(dir, point + Vec3::z() * d_p).sqrt() - err0) / d_p;

            let diff_vz = (evall_ls(
                Azimuthal {
                    z: dir.z + d_v,
                    ..dir
                },
                point,
            )
            .sqrt()
                - err0)
                / d_v;
            let diff_vh = (evall_ls(
                Azimuthal {
                    h: dir.h + d_v,
                    ..dir
                },
                point,
            )
            .sqrt()
                - err0)
                / d_v;

            (diff_x, diff_y, diff_z, diff_vz, diff_vh)
        };

        let mut point_descent_coeff: f64 = 1e6;
        let mut direction_descent_coeff: f64 = 1.0;

        for _ in 0..iterations {
            let mut dir_azimuthal: Azimuthal = traj.direction.into_inner().into();
            let (x, y, z, vz, vh) = diff(traj.point, dir_azimuthal);

            let old = traj;

            traj.point.x -= x * point_descent_coeff;
            traj.point.y -= y * point_descent_coeff;
            traj.point.z -= z * point_descent_coeff;

            if self.eval_mean_of_squares(traj, weights) > self.eval_mean_of_squares(old, weights) {
                traj = old;
                point_descent_coeff *= 0.5;
            }

            //----------------------------------------------------

            let mut dir_azimuthal: Azimuthal = traj.direction.into_inner().into();
            let (x, y, z, vz, vh) = diff(traj.point, dir_azimuthal);

            let old = traj;

            dir_azimuthal.z -= vz * direction_descent_coeff;
            dir_azimuthal.h -= vh * direction_descent_coeff;
            traj.direction = dir_azimuthal.into();

            if self.eval_mean_of_squares(traj, weights) > self.eval_mean_of_squares(old, weights) {
                traj = old;
                direction_descent_coeff *= 0.5;
            }
        }

        traj
    }

    pub fn compute_weights(&self, traj: Line) -> Vec<(f64, f64, f64)> {
        let best_err = self.eval_median_of_squares(traj).sqrt()
            * 1.483
            * (1. + 5. / (self.data.samples.len() - 6) as f64);
        let get_weight = |err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const O: f64 = 1.5;
            const F: f64 = 0.6;
            0.5 * (1.0 - ((err.abs() / best_err - O) / F).tanh())
        };

        let mut weights: Vec<(f64, f64, f64)> = Vec::with_capacity(self.data.samples.len());
        let mut sw = Vec::new();
        let mut ew = Vec::new();
        let mut daw = Vec::new();
        for s in &self.data.samples {
            let eval = self.evaluate_traj(s, traj);
            let b1 = eval.start.map_or(0.0, get_weight);
            let b2 = eval.end.map_or(0.0, get_weight);
            let b3 = eval.da.map_or(0.0, get_weight);
            weights.push((b1, b2, b3));
            eval.start.map(|_| sw.push(b1));
            eval.end.map(|_| ew.push(b2));
            eval.da.map(|_| daw.push(b3));
        }

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

    pub fn eval_mean_of_squares(&self, traj: Line, weights: &[(f64, f64, f64)]) -> f64 {
        let mut sum = 0.0;
        let mut cnt = 0.0;
        let mut push = |x, w| {
            sum += x * w;
            cnt += w;
        };
        for (s, &(w1, w2, w3)) in self.data.iter().zip(weights) {
            let eval = self.evaluate_traj(s, traj);
            eval.start.map(|x| push(x * x, w1));
            eval.end.map(|x| push(x * x, w2));
            eval.da.map(|x| push(x * x, w3));
        }
        sum / cnt
    }

    fn draw_plots(&self, traj: Line, weights: &[(f64, f64, f64)], stage: &str) {
        let mut points = Vec::new();
        let mut buf = Vec::new();
        // let mut buf2 = Vec::new();

        points.clear();
        buf.clear();
        for (i, s) in self.data.samples.iter().enumerate() {
            if let Some(da) = self.evaluate_traj(s, traj).da {
                let x = (i + 1) as f64;
                let y = da.to_degrees();
                buf.push(y);
                points.push((x, y, weight_to_rgb(weights[i].2), 2.5));
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

        // points.clear();
        // buf.clear();
        // for (i, s) in self.data.samples.iter().enumerate() {
        //     let Some(da_observed) = s.da else { continue };
        //     let k = self.data.answer.unwrap().traj.point - s.location;
        //     let da_comp = descent_angle(
        //         s.location,
        //         k,
        //         self.data.answer.unwrap().traj.direction.into_inner(),
        //     );
        //     // let x = i as f64;
        //     let x = da_observed.to_degrees();
        //     let y = da_comp.to_degrees();
        //     buf.push(y);
        //     buf2.push(x);
        //     points.push((x, y, weight_to_rgb(0.0), 2.5));
        //     // points.push((x, y2, weight_to_rgb(0.0), 2.5));
        // }
        // let med_comp = buf.median();
        // let med_observed = buf2.median();
        // eprintln!("da_comp med = {med_comp}");
        // eprintln!("da_observed med = {med_observed}");
        // draw_plot_svg(
        //     &format!("plots/da-{:.0}-{stage}.svg", traj.point.x),
        //     &points,
        //     &[],
        //     &[
        //         (0.0, colors::BLACK),
        //         // (med_comp, colors::GREEN),
        //         // (med_observed, colors::RED),
        //     ],
        // )
        // .unwrap();

        // use std::io::Write;
        // let color = match points.len() {
        //     203 => "Cyprus",
        //     18 => "Denmark",
        //     300 => "Florida",
        //     116 => "Italy",
        //     581 => "Netherlands",
        //     _ => unreachable!(),
        // };
        // let mut file = BufWriter::new(
        //     File::options()
        //         .create(true)
        //         .append(true)
        //         .open("da-points.txt")
        //         .unwrap(),
        // );
        // // writeln!(file, "x,y,color").unwrap();
        // for (x, y, _, _) in &points {
        //     writeln!(file, "{x},{y},{color}").unwrap();
        // }
        // eprintln!("Written {} points", points.len());

        // points.clear();
        // let mut buf = Vec::new();
        // for (i, s) in self.data.samples.iter().enumerate() {
        //     if let Some(error) = Self::evaluate_traj(s, traj).0 {
        //         if let Some(k_start) = s.k_start {
        //             let l_start = lambda(s.location, k_start, traj.point, traj.direction);
        //             let x = l_start * 1e-3;
        //             buf.push(x);
        //             let y = error.to_degrees();
        //             let weight = weights[i].0;
        //             if x.abs() > 6_000.0 {
        //                 eprintln!(
        //                     "Warning: point skipped. Resaon: x = {x:.00}km. Weight = {weight:.03}."
        //                 );
        //             } else {
        //                 points.push((x, y, weight_to_rgb(weight), 2.5));
        //             }
        //         }
        //     }
        // }
        // let med = buf.median();
        // draw_plot_svg(
        //     &format!("plots/ls-{:.0}.svg", traj.point.x),
        //     &points,
        //     &[(med, colors::BLUE)],
        //     &[(0.0, colors::BLACK)],
        // )
        // .unwrap();
        //
        // points.clear();
        // buf.clear();
        // for (i, s) in self.data.samples.iter().enumerate() {
        //     if let Some(error) = Self::evaluate_traj(s, traj).1 {
        //         if let Some(k_end) = s.k_end {
        //             let l_end = lambda(s.location, k_end, traj.point, traj.direction);
        //             let x = l_end * 1e-3;
        //             buf.push(x);
        //             let y = error.to_degrees();
        //             let weight = weights[i].1;
        //             if x.abs() > 6_000.0 {
        //                 eprintln!(
        //                     "Warning: point skipped. Resaon: x = {x:.00}km. Weight = {weight:.03}."
        //                 );
        //             } else {
        //                 points.push((x, y, weight_to_rgb(weight), 2.5));
        //             }
        //         }
        //     }
        // }
        // let med = buf.median();
        // draw_plot_svg(
        //     &format!("plots/le-{:.0}.svg", traj.point.x),
        //     &points,
        //     &[(med, colors::BLUE)],
        //     &[(0.0, colors::BLACK)],
        // )
        // .unwrap();
    }

    // pub fn altitude_correction_given_answer(&self, answer: Line) {
    //     let mut points = Vec::new();
    //
    //     let new_h = |s: &DataSample, k: UnitVec3| {
    //         let lm = lambda(s.location, k, answer.point, answer.direction);
    //         let k_p = answer.point + answer.direction.into_inner() * lm - s.location;
    //         let p: Azimuthal = Unit::new_normalize(k_p)
    //             .to_local(s.geo_location)
    //             .into_inner()
    //             .into();
    //         p.h
    //
    //         // let n = k.cross(&s.location);
    //         // let l = (s.location - answer.point).dot(&n) / answer.direction.dot(&n);
    //         // let p = answer.point + answer.direction.into_inner() * l;
    //         // let k_new = (p - s.location).normalize();
    //         // s.location.normalize().dot(&k_new).asin()
    //     };
    //
    //     points.clear();
    //     let mut buf_x = Vec::new();
    //     let mut buf_y = Vec::new();
    //     for s in self.data.samples.iter() {
    //         if let Some(k) = s.k_start {
    //             let az: Azimuthal = k.to_local(s.geo_location).into_inner().into();
    //             let x = az.h.to_degrees();
    //             let y = new_h(s, k).to_degrees();
    //             if x.abs() > 6_000.0 {
    //                 eprintln!("Warning: point skipped. Resaon: x = {x:.00}km.");
    //             } else {
    //                 buf_x.push(x);
    //                 buf_y.push(y);
    //                 points.push((x, y, weight_to_rgb(s.exp.unwrap_or(1.0) / 5.0), 2.5));
    //             }
    //         }
    //     }
    //     let med_x = buf_x.median();
    //     let med_y = buf_y.median();
    //
    //     draw_plot_svg(
    //         &format!("plots/ac-ls-{:.0}.svg", answer.point.x),
    //         &points,
    //         &[(med_x, colors::BLUE)],
    //         &[(0.0, colors::BLACK), (med_y, colors::BLUE)],
    //     )
    //     .unwrap();
    //
    //     points.clear();
    //     let mut buf_x = Vec::new();
    //     let mut buf_y = Vec::new();
    //     for s in self.data.samples.iter() {
    //         if let Some(k) = s.k_end {
    //             let az: Azimuthal = k.to_local(s.geo_location).into_inner().into();
    //             let x = az.h.to_degrees();
    //             let y = new_h(s, k).to_degrees();
    //             if x.abs() > 6_000.0 {
    //                 eprintln!("Warning: point skipped. Resaon: x = {x:.00}km.");
    //             } else {
    //                 buf_x.push(x);
    //                 buf_y.push(y);
    //                 points.push((x, y, weight_to_rgb(s.exp.unwrap_or(1.0) / 5.0), 2.5));
    //             }
    //         }
    //     }
    //     let med_x = buf_x.median();
    //     let med_y = buf_y.median();
    //     draw_plot_svg(
    //         &format!("plots/ac-le-{:.0}.svg", answer.point.x),
    //         &points,
    //         &[(med_x, colors::BLUE)],
    //         &[(0.0, colors::BLACK), (med_y, colors::BLUE)],
    //     )
    //     .unwrap();
    //     eprintln!("{med_x:.0} -> {med_y:.0}");
    // }

    fn calc_sigmas(&mut self, traj: Line, weights: &[(f64, f64, f64)]) -> [f64; 6] {
        let (flash, speed) = self.calc_flash_and_speed(traj, weights);
        let dir: Azimuthal = traj.direction.into_inner().into();

        dbg!(flash);
        dbg!(dir);

        let mut sum = [0.0; 6];

        const DA_D: f64 = 0.01;

        for i in 0..(self.data.iter().len()) {
            if self.data.samples[i].da.is_none() {
                continue;
            }

            if let Some(da_err) = self.evaluate_traj(&self.data.samples[i], traj).da {
                let old_da = self.data.samples[i].da;
                *self.data.samples[i].da.as_mut().unwrap() += DA_D;

                println!("{i}: flash.x = {} ", flash.x);
                let new_traj = self.gradient_descent(traj, weights, 1_000);
                let (new_flash, _new_speed) = self.calc_flash_and_speed(new_traj, weights);
                println!("{i}: flash.x(1_000) = {} ", new_flash.x);
                let new_traj = self.gradient_descent(new_traj, weights, 9_000);
                let (new_flash, _new_speed) = self.calc_flash_and_speed(new_traj, weights);
                println!("{i}: flash.x(10_000) = {} ", new_flash.x);
                let new_traj = self.gradient_descent(new_traj, weights, 40_000);
                let (new_flash, _new_speed) = self.calc_flash_and_speed(new_traj, weights);
                println!("{i}: flash.x(50_000) = {} ", new_flash.x);
                let new_traj = self.gradient_descent(new_traj, weights, 50_000);
                let (new_flash, new_speed) = self.calc_flash_and_speed(new_traj, weights);
                println!("{i}: flash.x(100_000) = {} ", new_flash.x);
                println!();

                let new_dir: Azimuthal = new_traj.direction.into_inner().into();
                self.data.samples[i].da = old_da;

                let mul = weights[i].2 * da_err * da_err / DA_D / DA_D;
                sum[0] += mul * f64::powi(new_flash.x - flash.x, 2);
                sum[1] += mul * f64::powi(new_flash.y - flash.y, 2);
                sum[2] += mul * f64::powi(new_flash.z - flash.z, 2);
                sum[3] += mul * f64::powi(new_dir.z - dir.z, 2);
                sum[4] += mul * f64::powi(new_dir.h - dir.h, 2);
                sum[5] += mul * f64::powi(new_speed - speed, 2);
            }
        }

        sum.map(f64::sqrt)
    }

    fn flip_da(&mut self, traj: Line) -> usize {
        let mut flipped = 0;

        for i in 0..self.data.samples.len() {
            let Some(dda) = self.evaluate_traj(&self.data.samples[i], traj).da
            else { continue };

            if dda.abs() > FRAC_PI_2 {
                let da = self.data.samples[i].da.as_mut().unwrap();
                *da = (*da + PI) % TAU;
                flipped += 1;
            }
        }

        flipped
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        // self.enhance_data();

        let traj = self.pairwise();
        self.data.compare(traj, "Initial guess (pairwise)");
        self.draw_plots(
            traj,
            &vec![(1.0, 1.0, 1.0); self.data.samples.len()],
            "initial",
        );

        let traj = self.lms(traj);
        self.data.compare(traj, "After LMS");
        self.draw_plots(
            traj,
            &vec![(1.0, 1.0, 1.0); self.data.samples.len()],
            "after-lms",
        );

        let weights = self.compute_weights(traj);
        self.draw_plots(traj, &weights, "with-weights");

        let traj = self.weighted_ls(traj, &weights);
        self.data.compare(traj, "After LS");

        self.draw_plots(traj, &weights, "after-ls");

        // let traj = self.gradient_descent(traj, &weights, 100_000);
        // self.data.compare(traj, "After Gradient Descent");
        // self.draw_plots(traj, &weights, "final");

        dbg!(self.flip_da(traj));

        let weights = self.compute_weights(traj);
        self.draw_plots(traj, &weights, "with-weights");

        let traj = self.weighted_ls(traj, &weights);
        self.data.compare(traj, "After LS");

        self.draw_plots(traj, &weights, "after-ls");

        // let traj = self.gradient_descent(traj, &weights, 100_000);
        // self.data.compare(traj, "After Gradient Descent");
        // self.draw_plots(traj, &weights, "final");

        // dbg!(self.calc_sigmas(traj, &weights));

        let (flash, speed) = self.calc_flash_and_speed(traj, &weights);

        Solution {
            flash,
            velocity: traj.direction.into_inner() * speed,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, traj: Line, weights: &[(f64, f64, f64)]) -> (Vec3, f64) {
        let mut lambdas = Vec::new();
        let mut lambdas_all = Vec::new();
        let mut speeds = Vec::new();

        let mut skipped_l_end_w = 0;
        let mut skipped_l_end_dir = 0;
        let mut total_l_end = 0;

        for (i, s) in self.data.samples.iter().enumerate() {
            if s.k_end.is_some() {
                total_l_end += 1;
            }

            if !s.observation_matches(traj) {
                skipped_l_end_dir += 1;
                continue;
            }

            if let Some(k_end) = s.k_end {
                let l_end = lambda(s.location, k_end, traj.point, traj.direction);

                if weights[i].1 > 0.5 {
                    lambdas.push(l_end);
                    lambdas_all.push(l_end);
                } else {
                    lambdas_all.push(l_end);
                    skipped_l_end_w += 1;
                }

                if let Some((duration, k_start)) = s.dur.zip(s.k_start) {
                    let l_start = lambda(s.location, k_start, traj.point, traj.direction);
                    let dist = l_end - l_start;
                    speeds.push(dist / duration);
                }
            }
        }

        eprintln!("skipped {skipped_l_end_w}(weight) + {skipped_l_end_dir}(direction) l_end out of {total_l_end}");

        let l_end = lambdas.median();
        let l_end_all = lambdas_all.median();

        dbg!(l_end - l_end_all);

        let speed = speeds.median();
        (traj.point + traj.direction.into_inner() * l_end, speed)
    }

    pub fn evaluate_traj(&self, sample: &DataSample, traj: Line) -> Evaluation {
        Evaluation {
            start: None,
            end: if !self.params.da_only {
                sample.z0.map(|z0| {
                    let zf = sample.calc_azimuth(traj.point);
                    // TODO: what when wrong?
                    // let Azimuthal { z: zf, h: _ } =
                    //     UnitVec3::new_normalize(traj.point - sample.location)
                    //         .to_local(sample.geo_location)
                    //         .into_inner()
                    //         .into();
                    angle_diff(z0, zf)
                })
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

pub struct Evaluation {
    pub start: Option<f64>,
    pub end: Option<f64>,
    pub da: Option<f64>,
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

fn lerp(a: f64, b: f64, w: f64) -> f64 {
    a + (b - a) * w
}
