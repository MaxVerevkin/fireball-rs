//! The implementation itself

use std::sync::Arc;

use crate::data::Data;
use crate::maths::*;
use crate::structs::*;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Arc<Data>,
    params: Params,
}

/// Contains all necessary information to solve the problem
#[derive(Clone, Copy)]
pub struct Params {
    pub min_match: f64,
    pub range: f64,
}

/// The answer
pub struct Solution {
    pub error: f64,
    pub flash: Spherical,
    pub velocity: Vec3,
    pub raw: (Vec3, Vec3),
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver {
            data: Arc::new(data),
            params,
        }
    }

    /// Find the solution
    pub fn solve(&self) -> Solution {
        let (p1, p2) = self.monte_carlo(self.data.mean_pos, self.params.range, 5_000., 0.35);

        //#[cfg(debug_assertions)]
        //{
        //// Draw a plot
        //use std::io::Write;
        //let mut file = std::fs::File::create("data_solved.dat").expect("create failed");
        //for dx in -1000..1000 {
        //let x = dx as f64 * 200.;
        //file.write_all(
        //format!(
        //"{} {}\n",
        //x / 1000.,
        //self.evaluate_traj(p1 + Vec3 { x, y: 0., z: 0. }, p2)
        //)
        //.as_bytes(),
        //)
        //.unwrap();
        //}
        //}

        // calculate the velocity
        let velocity = (p2 - p1).normalized();

        // calculate flash location as well as speed
        let (flash, speed) = self.calc_flash_and_speed(p1, velocity);

        // return
        Solution {
            error: (self.evaluate_traj(p1, p2)),
            flash,
            velocity: velocity * speed,
            raw: (p1, p2),
        }
    }

    fn monte_carlo(
        &self,
        mid_point: Vec3,
        range: f64,
        target_range: f64,
        range_mul: f64,
    ) -> (Vec3, Vec3) {
        let fx = || {
            let mut error = std::f64::INFINITY;
            let mut answer = (mid_point, mid_point);
            for _ in 0..20_000 {
                // Generate two points
                let p1 = mid_point + Vec3::rand(range);
                let p2 = mid_point + Vec3::rand(range);

                let err = self.evaluate_traj(p1, p2);
                if err < error {
                    answer = (p1, p2);
                    error = err;
                }
            }
            (answer, error)
        };

        let (a0, a1, a2, a3) = crossbeam_utils::thread::scope(|scope| {
            let a1 = scope.spawn(|_| fx());
            let a2 = scope.spawn(|_| fx());
            let a3 = scope.spawn(|_| fx());

            let a0 = fx();

            let a1 = a1.join().unwrap();
            let a2 = a2.join().unwrap();
            let a3 = a3.join().unwrap();

            (a0, a1, a2, a3)
        })
        .unwrap();

        let answer = if a0.1 < a1.1 { a0 } else { a1 };
        let answer = if answer.1 < a2.1 { answer } else { a2 };
        let answer = if answer.1 < a3.1 { answer } else { a3 };
        let answer = answer.0;

        if range > target_range {
            self.monte_carlo(
                (answer.0 + answer.1) * 0.5,
                range * range_mul,
                target_range,
                range_mul,
            )
        } else {
            answer
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, point: Vec3, v: Vec3) -> (Spherical, f64) {
        let mut l_end_mean = 0.;
        let mut speed = 0.;
        let mut l_count = 0.;
        let mut speed_count = 0.;
        for sample in &self.data.samples {
            let pip = point - sample.global_pos;
            let plane = pip.cross(v).normalized();
            let perpendic = v.cross(plane).normalized();

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            //let trust_start = sample.trust_start && perpendic.dot(k_start) > self.params.min_match;
            //let trust_end = sample.trust_end && perpendic.dot(k_end) > self.params.min_match;

            //if trust_end {
            let l_end = pip.dot(perpendic * (k_end.dot(v) / k_end.dot(perpendic)) - v);
            l_end_mean += l_end;
            l_count += 1.;
            //if trust_start {
            let l_start = pip.dot(perpendic * (k_start.dot(v) / k_start.dot(perpendic)) - v);
            speed += (l_end - l_start) / sample.duration;
            speed_count += 1.;
            //}
            //}
        }
        l_end_mean /= l_count;
        speed /= speed_count;
        //dbg!(l_count);
        //dbg!(speed_count);
        ((point + v * l_end_mean).into(), speed)
    }

    /// Calculate the mean of squared errors (less is better)
    pub fn evaluate_traj(&self, p1: Vec3, p2: Vec3) -> f64 {
        // Get normalized velocity
        let mut vel = p1 - p2;
        // Check if two points are too close.
        if vel.length() < 10. {
            return f64::INFINITY;
        }
        vel.normalize();

        let mut count = 0.;
        let mut error = 0.;
        for sample in &self.data.samples {
            let plane = (p1 - sample.global_pos).cross(vel).normalized();
            //let perpendic = vel.cross(plane);

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let e_start = plane.dot(k_start).asin();
            let e_end = plane.dot(k_end).asin();

            //let trust_start = sample.trust_start && perpendic.dot(k_start) > self.params.min_match;
            //let trust_end = sample.trust_end && perpendic.dot(k_end) > self.params.min_match;

            //if trust_start {
            error += e_start * e_start;
            count += 1.;

            let z = sample.global_pos.normalized();
            let x = k_start.cross(z).normalized();

            let nx = plane.dot(x);
            let nz = plane.dot(z);
            let mut alpha = f64::atan2(nx, nz) + std::f64::consts::FRAC_PI_2;
            if alpha < 0. {
                alpha += std::f64::consts::TAU;
            }

            let diff = angle_diff(alpha, sample.descent_angle);
            error += diff * diff * 0.5;
            count += 1.;
            //}
            //if trust_end {
            error += e_end * e_end;
            count += 1.;
            //}
        }

        //dbg!(count);

        error / count
    }
}
