//! The implementation itself

use crate::data::Data;
use crate::maths::*;
use crate::structs::*;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Contains all necessary information to solve the problem
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub initial_range: f64,
    pub initial_iterations: usize,
    pub main_iterations: usize,
    pub threads: usize,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub error: f64,
    pub flash: Spherical,
    pub velocity: Vec3,
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver { data, params }
    }

    /// Find the solution
    pub fn solve(&self) -> Solution {
        // Get a rough estimation
        let (p1, p2) = self.initial_mc();

        // Construct a tunnel
        let k = (p2 - p1).normalized();
        let (l1, l2) = self.lambdas(p1, k);
        let mut tunnel = Tunnel {
            mid_point: p1 + k * ((l1 + l2) * 0.5),
            k,
            r: 200_000.,
        };

        // Run 10 iterations
        // TODO make configurable
        for _ in 0..10 {
            let (mid_point, k) = self.mc(tunnel);
            //let (l1, l2) = self.lambdas(mid_point, k);
            //tunnel.mid_point = mid_point + k * ((l1 + l2) * 0.5);
            tunnel.mid_point = mid_point;
            tunnel.k = k;
            tunnel.r *= 0.5;
        }

        // calculate flash location as well as speed
        let (flash, speed) = self.calc_flash_and_speed(tunnel.mid_point, tunnel.k);
        let velocity = tunnel.k * speed;

        // return
        Solution {
            error: self.evaluate_traj(tunnel.mid_point, tunnel.k),
            flash,
            velocity,
        }
    }

    fn initial_mc(&self) -> (Vec3, Vec3) {
        let mid_point = self.data.mean_pos;
        let range = self.params.initial_range;

        let answers = crossbeam_utils::thread::scope(|scope| {
            let mut threads = Vec::with_capacity(self.params.threads);
            for _ in 0..self.params.threads {
                threads.push(scope.spawn(|_| {
                    let mut error = std::f64::INFINITY;
                    let mut answer = (mid_point, mid_point);
                    for _ in 0..(self.params.initial_iterations / self.params.threads) {
                        // Generate two points
                        let p1 = mid_point + Vec3::rand_uniform(range);
                        let p2 = mid_point + Vec3::rand_uniform(range);

                        let err = self.evaluate_traj(p1, (p2 - p1).normalized());
                        if err < error {
                            answer = (p1, p2);
                            error = err;
                        }
                    }
                    (error, answer)
                }));
            }
            threads
                .into_iter()
                .map(|t| t.join().unwrap())
                .collect::<Vec<(f64, (Vec3, Vec3))>>()
        })
        .unwrap();

        answers
            .iter()
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .unwrap()
            .1
    }

    fn mc(&self, tunnel: Tunnel) -> (Vec3, Vec3) {
        let answers = crossbeam_utils::thread::scope(|scope| {
            let mut threads = Vec::with_capacity(self.params.threads);
            for _ in 0..self.params.threads {
                threads.push(scope.spawn(|_| {
                    let mut error = std::f64::INFINITY;
                    let mut answer = tunnel;
                    for _ in 0..(self.params.main_iterations / self.params.threads) {
                        let rand_ans = tunnel.random();
                        let err = self.evaluate_traj(rand_ans.mid_point, rand_ans.k);
                        if err < error {
                            answer = rand_ans;
                            error = err;
                        }
                    }
                    (error, (answer.mid_point, answer.k))
                }));
            }
            threads
                .into_iter()
                .map(|t| t.join().unwrap())
                .collect::<Vec<(f64, (Vec3, Vec3))>>()
        })
        .unwrap();

        answers
            .iter()
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .unwrap()
            .1
    }

    fn lambdas(&self, point: Vec3, vel: Vec3) -> (f64, f64) {
        let mut l_begin_vec = Vec::with_capacity(self.data.samples.len());
        let mut l_end_vec = Vec::with_capacity(self.data.samples.len());

        for sample in &self.data.samples {
            let pip = point - sample.global_pos;
            let plane = pip.cross(vel).normalized();
            let perpendic = vel.cross(plane).normalized();

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let trust_start = sample.trust_start && perpendic.dot(k_start) > 0.;
            let trust_end = sample.trust_end && perpendic.dot(k_end) > 0.;

            if trust_start {
                l_begin_vec
                    .push(pip.dot(perpendic * (k_start.dot(vel) / k_start.dot(perpendic)) - vel));
            }
            if trust_end {
                l_end_vec.push(pip.dot(perpendic * (k_end.dot(vel) / k_end.dot(perpendic)) - vel));
            }
        }

        let l_begin = quick_median(&mut l_begin_vec);
        let l_end = quick_median(&mut l_end_vec);

        (l_begin, l_end)
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, point: Vec3, v: Vec3) -> (Spherical, f64) {
        let mut speed_vec = Vec::new();
        let mut l_end_vec = Vec::new();

        for sample in &self.data.samples {
            let pip = point - sample.global_pos;
            let plane = pip.cross(v).normalized();
            let perpendic = v.cross(plane).normalized();

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let trust_start = sample.trust_start && perpendic.dot(k_start) > 0.;
            let trust_end = sample.trust_end && perpendic.dot(k_end) > 0.;

            if trust_end {
                let l_end = pip.dot(perpendic * (k_end.dot(v) / k_end.dot(perpendic)) - v);
                l_end_vec.push(l_end);
                if trust_start && sample.duration > 0. {
                    let l_start =
                        pip.dot(perpendic * (k_start.dot(v) / k_start.dot(perpendic)) - v);
                    speed_vec.push((l_end - l_start) / sample.duration);
                }
            }
        }

        let l_end = quick_median(&mut l_end_vec);
        let speed = quick_median(&mut speed_vec);

        //dbg!(l_end_vec.len());
        //dbg!(speed_vec.len());

        ((point + v * l_end).into(), speed)
    }

    /// Calculate the mean of squared errors (less is better)
    pub fn evaluate_traj(&self, p: Vec3, vel: Vec3) -> f64 {
        let mut count = 0.;
        let mut error = 0.;
        for sample in &self.data.samples {
            let plane = (p - sample.global_pos).cross(vel).normalized();
            //let perpendic = vel.cross(plane);

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let e_start = plane.dot(k_start).asin();
            let e_end = plane.dot(k_end).asin();

            //if sample.trust_start && perpendic.dot(k_start) > 0. {
            if sample.trust_start {
                error += e_start * e_start;
                count += 1.;

                let angle = descent_angle(sample.global_pos, k_start, vel);
                let diff = angle_diff(angle, sample.descent_angle);
                error += diff * diff * 0.5;
                //error += diff * diff;
                count += 1.;
            }
            //if sample.trust_end && perpendic.dot(k_end) > 0. {
            if sample.trust_end {
                error += e_end * e_end;
                count += 1.;
            }
        }

        //dbg!(count);

        error / count
    }
}
