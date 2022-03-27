//! The implementation

use rand::{Rng, SeedableRng};
use std::f64::consts::FRAC_PI_2;

// use crate::constants::EARTH_R;
use crate::data::{Data, DataSample};
use crate::maths::*;
use crate::quick_median::quick_median;
use crate::structs::*;

pub mod pair;
use pair::PairTrajectory;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Parameters to tweak the solution algorithm
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Whether to tilt the observation to better match given descent angle
    pub da_correction: f64,
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

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        // TODO explain what's going on here
        if self.params.da_correction > 0.0 {
            for s in &mut self.data.samples {
                if let Some(((da, axis), comp_da)) = s.da.zip(s.axis).zip(s.calculated_da()) {
                    let angle = angle_diff(comp_da, da) * self.params.da_correction;
                    let q = UnitQuaternion::new(axis, angle);
                    s.plane.as_mut().map(|k| *k = q * *k);
                    s.k_start.as_mut().map(|k| *k = q * *k);
                    s.k_end.as_mut().map(|k| *k = q * *k);
                }
            }
        }

        // The initial guess is the weighted sum of pairwise plane crossings
        let traj = {
            let mut traj = Line::default();
            let mut sum = 0.0;
            for (i, s1) in self.data.iter().enumerate() {
                for s2 in self.data.iter().skip(i + 1) {
                    if let Some(line) = PairTrajectory::calculate(s1, s2) {
                        traj += line.line * line.weight;
                        sum += line.weight;
                    }
                }
            }
            traj.point /= sum;
            traj.direction.normalize();
            traj
        };

        self.data.answer.map(|a| a.compare(traj, "Initial guess"));

        // TODO remove (use traj instead)
        let mut point = traj.point;
        let mut vel = traj.direction;

        let mut buf = Vec::new();
        let mut eval_lms = |vel: Vec3, point: Vec3| -> f64 {
            buf.clear();
            for s in &self.data.samples {
                let (eb, ee, eda) = Solver::evaluate_traj(s, point, vel);
                eb.map(|e| buf.push(e));
                ee.map(|e| buf.push(e));
                eda.map(|e| buf.push(e));
            }
            quick_median(&mut buf)
        };

        const ITER_CNT: usize = 500;
        const ITER_INNER_CNT: usize = 100;

        let mut rng = rand::rngs::SmallRng::from_entropy();
        let mut best_err = eval_lms(vel, point);
        let mut dif = 5000.0;
        for _ in 0..ITER_CNT {
            let distr = rand_distr::Normal::new(0.0, dif).unwrap();
            for _ in 0..ITER_INNER_CNT {
                let mut p = point;
                p.x += rng.sample(distr);
                p.y += rng.sample(distr);
                p.z += rng.sample(distr);
                let e = eval_lms(vel, p);
                if e < best_err {
                    // dbg!(e);
                    best_err = e;
                    point = p;
                }
            }
            dif *= 0.995;

            for _ in 0..ITER_INNER_CNT {
                let v = vel.tilt_random(4f64.to_radians(), &mut rng);
                let e = eval_lms(v, point);
                if e < best_err {
                    // dbg!(e);
                    best_err = e;
                    vel = v;
                }
            }
        }

        let best_err = best_err.sqrt() * 1.483 * (1. + 5. / (self.data.samples.len() - 6) as f64);

        // dbg!(best_err);

        let get_weight = |err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const F: f64 = 0.7;
            const O: f64 = 2.0;
            0.5 * (1.0 - (((err / best_err) - O) / F).tanh())
        };

        let mut weights: Vec<(f64, f64, f64)> = Vec::with_capacity(self.data.samples.len());
        let mut sw = Vec::new();
        let mut ew = Vec::new();
        let mut daw = Vec::new();
        for s in &self.data.samples {
            let (eb, ee, eda) = Solver::evaluate_traj(s, point, vel);
            let b1 = eb.map(f64::sqrt).map_or(0.0, get_weight);
            let b2 = ee.map(f64::sqrt).map_or(0.0, get_weight);
            let b3 = eda.map(f64::sqrt).map_or(0.0, get_weight);
            weights.push((b1, b2, b3));
            eb.map(|_| sw.push(b1));
            ee.map(|_| ew.push(b2));
            eda.map(|_| daw.push(b3));
        }

        if !sw.is_empty() {
            println!("Start ({}):", sw.len());
            histogram::draw_hitogram(&sw, 5);
        }
        if !ew.is_empty() {
            println!("End ({}):", ew.len());
            histogram::draw_hitogram(&ew, 5);
        }
        if !daw.is_empty() {
            println!("DA ({}):", daw.len());
            histogram::draw_hitogram(&daw, 5);
        }

        let eval_ls = |vel: Vec3, point: Vec3| -> f64 {
            let mut sum = 0.0;
            let mut cnt = 0.0;
            let mut push = |x, w| {
                sum += x * w;
                cnt += w;
            };
            for (s, &(w1, w2, w3)) in self.data.iter().zip(&weights) {
                let (eb, ee, eda) = Solver::evaluate_traj(s, point, vel);
                eb.map(|x| push(x, w1));
                ee.map(|x| push(x, w2));
                eda.map(|x| push(x, w3));
            }
            sum / cnt
        };

        let mut best_err = eval_ls(vel, point);
        let mut dif = 5000.0;
        for _ in 0..ITER_CNT {
            let distr = rand_distr::Normal::new(0.0, dif).unwrap();
            for _ in 0..ITER_INNER_CNT {
                let mut p = point;
                p.x += rng.sample(distr);
                p.y += rng.sample(distr);
                p.z += rng.sample(distr);
                let e = eval_ls(vel, p);
                if e < best_err {
                    // dbg!(e);
                    best_err = e;
                    point = p;
                }
            }
            dif *= 0.999;

            for _ in 0..ITER_INNER_CNT {
                let v = vel.tilt_random(5f64.to_radians(), &mut rng);
                let e = eval_ls(v, point);
                if e < best_err {
                    // dbg!(e);
                    best_err = e;
                    vel = v;
                }
            }
        }

        self.data.answer.map(|a| {
            a.compare(
                Line {
                    point,
                    direction: vel,
                },
                "Final answer",
            )
        });

        if let Some(answer) = self.data.answer {
            use scad_gen::constructors::*;
            use scad_gen::scope;
            let mut doc = scad_gen::Document::default();

            let scale = 1e5;
            macro_rules! cube {
                () => {
                    cube(Vec3::new(1., 1., 1.) / 2.0)
                };
            }
            let ans = point;
            let real_vel = vel;
            let vel = vel * 30.;

            doc.color(
                "red",
                Some(0.5),
                scope!(ctx, {
                    ctx.hull(|ctx| {
                        ctx.translate(ans / scale - vel * 10.0, cube!());
                        ctx.translate(ans / scale, cube!());
                    })
                }),
            );

            doc.color(
                "green",
                None,
                scope!(ctx, {
                    ctx.hull(|ctx| {
                        ctx.translate(answer.0.point / scale - answer.0.direction * 10.0, cube!());
                        ctx.translate(answer.0.point / scale, cube!());
                    })
                }),
            );

            doc.color(
                "blue",
                Some(0.25),
                sphere((crate::constants::EARTH_R / scale) as f32),
            );

            {
                for s in &self.data.samples {
                    if !s.observation_matches(point, real_vel) {
                        continue;
                    }
                    if let Some(k_end) = s.k_end {
                        let l_end = lambda(s.location, k_end, point, real_vel);
                        let p = point + real_vel * l_end;
                        doc.translate(p / scale, cube!());
                    }
                }
            }

            use std::fs::File;
            use std::io::prelude::*;
            let mut file = File::create("visual.scad").unwrap();
            write!(file, "{doc}").unwrap();
        }

        //         vec.sort_unstable_by(|a, b| a.0.x.partial_cmp(&b.0.x).unwrap());
        //         let x = vec[vec.len() / 2].0.x;
        //         vec.sort_unstable_by(|a, b| a.0.y.partial_cmp(&b.0.y).unwrap());
        //         let y = vec[vec.len() / 2].0.y;
        //         vec.sort_unstable_by(|a, b| a.0.z.partial_cmp(&b.0.z).unwrap());
        //         let z = vec[vec.len() / 2].0.z;
        //         let point = Vec3 { x, y, z };
        //
        //         vec.sort_unstable_by(|a, b| a.1.x.partial_cmp(&b.1.x).unwrap());
        //         let x = vec[vec.len() / 2].1.x;
        //         vec.sort_unstable_by(|a, b| a.1.y.partial_cmp(&b.1.y).unwrap());
        //         let y = vec[vec.len() / 2].1.y;
        //         vec.sort_unstable_by(|a, b| a.1.z.partial_cmp(&b.1.z).unwrap());
        //         let z = vec[vec.len() / 2].1.z;
        //         let vel = Vec3 { x, y, z }.normalized();

        // dbg!(vec.len());

        // let mut s_point = 0.;
        // let mut s_vel = 0.;
        // for (p, v) in &vec {
        //     s_point += (point - *p).len().powi(2);
        //     s_vel += vel.dot(*v).abs().min(1.).acos().powi(2);
        // }
        // s_point = (s_point / vec.len() as f64).sqrt();
        // s_vel = (s_vel / vec.len() as f64).sqrt();

        // dbg!(s_point);
        // dbg!(s_vel.to_degrees());

        let (flash, speed) = self.calc_flash_and_speed(point, vel);
        let velocity = vel * speed;

        Solution {
            // error: self.evaluate_traj(point, vel),
            flash,
            velocity,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, point: Vec3, v: Vec3) -> (Vec3, f64) {
        let mut lambdas = Vec::new();
        let mut speeds = Vec::new();

        for s in &self.data.samples {
            if !s.observation_matches(point, v) {
                continue;
            }

            if let Some(k_end) = s.k_end {
                let l_end = lambda(s.location, k_end, point, v);
                // let l_end = lambda_corssing(s.location, k_end, point, v);
                lambdas.push(l_end);

                if let Some((duration, k_start)) = s.dur.zip(s.k_start) {
                    let l_start = lambda(s.location, k_start, point, v);
                    // let l_start = lambda_corssing(s.location, k_start, point, v);
                    let dist = l_end - l_start;
                    debug_assert!(dist > 0.0);
                    speeds.push(dist / duration);
                }
            }
        }

        let l_end = quick_median(&mut lambdas);
        let speed = quick_median(&mut speeds);
        (point + v * l_end, speed)
    }

    /// Calculate the mean of squared errors (less is better)
    pub fn evaluate_traj(
        sample: &DataSample,
        p: Vec3,
        vel: Vec3,
    ) -> (Option<f64>, Option<f64>, Option<f64>) {
        let plane = (p - sample.location).cross(vel).normalized();
        let perpendic = vel.cross(plane);

        let mut es = None;
        let mut ee = None;
        let mut eda = None;

        if let Some(k) = sample.k_start {
            es = if perpendic.dot(k) > 0. {
                Some(plane.dot(k).asin().powi(2))
            } else {
                Some(FRAC_PI_2 * FRAC_PI_2)
            };
        }

        if let Some(k) = sample.k_end {
            ee = if perpendic.dot(k) > 0. {
                Some(plane.dot(k).asin().powi(2))
            } else {
                Some(FRAC_PI_2 * FRAC_PI_2)
            };
        }

        if let Some(da) = sample.da {
            let k = (p - sample.location).normalized();
            let err = angle_diff(descent_angle(sample.location, k, vel), da);
            eda = Some(err * err);
        }

        let scale_err = |e| e * 4.0;
        (es.map(scale_err), ee.map(scale_err), eda)

        // (None, None, eda)
    }
}
