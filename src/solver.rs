//! The implementation itself

use crate::constants::EARTH_R;
use crate::data::{Data, DataSample};
use crate::maths::*;
use crate::structs::*;
// use std::f64::consts::{FRAC_PI_2, PI};

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Contains all necessary information to solve the problem
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub threads: usize,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    // pub error: f64,
    pub flash: Vec3,
    pub velocity: Vec3,
}

fn solve_pair(s1: &DataSample, s2: &DataSample) -> Option<(Vec3, Vec3)> {
    let n1 = s1.global_start.cross(s1.global_end).normalized();
    let n2 = s2.global_start.cross(s2.global_end).normalized();

    let mut v = n1.cross(n2).normalized();
    if !v.is_normal() {
        return None;
    }

    let vp1 = s1.global_end - s1.global_start;
    let vp2 = s2.global_end - s2.global_start;
    match (v.dot(vp1) > 0., v.dot(vp2) > 0.) {
        (true, false) | (false, true) => return None,
        (false, false) => v = -v,
        _ => (),
    }

    let l = (s2.global_pos - s1.global_pos).dot(n2) / s1.global_end.dot(n2);
    let point = s1.global_pos + s1.global_end * l;

    let zz: Spherical = point.into();
    if zz.r - EARTH_R > 100_000. {
        return None;
    }

    Some((point, v))
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver { data, params }
    }

    /// Find the solution
    pub fn solve(&self) -> Solution {
        // TODO multithreading

        let mut vec = Vec::new();
        let mut point = Vec3::default();
        let mut vel = Vec3::default();

        for i in 0..(self.data.samples.len() - 1) {
            let s1 = &self.data.samples[i];
            if s1.trust_start && s1.trust_end {
                for j in (i + 1)..self.data.samples.len() {
                    let s2 = &self.data.samples[j];
                    if s2.trust_start && s2.trust_end {
                        if let Some((p, v)) = solve_pair(s1, s2) {
                            point += p;
                            vel += v;
                            vec.push((p, v));
                        }
                    }
                }
            }
        }

        point /= vec.len() as f64;
        vel.normalize();

        dbg!(vec.len());

        let mut s_point = 0.;
        let mut s_vel = 0.;
        for (p, v) in &vec {
            s_point += (point - *p).length().powi(2);
            s_vel += vel.dot(*v).abs().min(1.).acos().powi(2);
        }
        s_point = (s_point / vec.len() as f64).sqrt();
        s_vel = (s_vel / vec.len() as f64).sqrt();

        dbg!(s_point);
        dbg!(s_vel.to_degrees());

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

        (point + v * l_end, speed)
    }

    // /// Calculate the mean of squared errors (less is better)
    //     pub fn evaluate_traj(&self, p: Vec3, vel: Vec3) -> f64 {
    //         let mut error = 0.;
    //         for sample in &self.data.samples {
    //             let plane = (p - sample.global_pos).cross(vel).normalized();
    //             let perpendic = vel.cross(plane);
    //
    //             let k_start = sample.global_start;
    //             let k_end = sample.global_end;
    //
    //             let e_start = plane.dot(k_start).asin();
    //             let e_end = plane.dot(k_end).asin();
    //
    //             if sample.trust_start {
    //                 let c = perpendic.dot(k_start);
    //                 if c > 0. {
    //                     error += e_start * e_start;
    //
    //                     let angle = descent_angle(sample.global_pos, k_start, vel);
    //                     let diff = angle_diff(angle, sample.descent_angle);
    //                     error += diff * diff;
    //                 } else {
    //                     //let c = c.clamp(0., 0.1) * 10.;
    //                     //error += c * e_start * e_start + (1. - c) * FRAC_PI_2 * FRAC_PI_2;
    //                     error += FRAC_PI_2 * FRAC_PI_2;
    //                     error += PI * PI;
    //                 }
    //                 //error += diff * diff;
    //             }
    //             if sample.trust_end {
    //                 let c = perpendic.dot(k_end);
    //                 if c > 0. {
    //                     error += e_end * e_end;
    //                 } else {
    //                     //let c = c.clamp(0., 0.1) * 10.;
    //                     //error += c * e_end * e_end + (1. - c) * FRAC_PI_2 * FRAC_PI_2;
    //                     error += FRAC_PI_2 * FRAC_PI_2;
    //                 }
    //             }
    //         }
    //         error
    //     }
}
