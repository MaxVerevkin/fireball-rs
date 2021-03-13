//! The implementation itself

use rand::prelude::*;

use crate::data::Data;
use crate::params::Params;

use crate::constants::*;
use crate::maths::*;
use crate::structs::*;

/// Contains all necessary information to solve the problem
pub struct Solver<'a> {
    data: &'a Data,
    params: &'a Params,
}
/// The answer
pub struct Solved {
    pub error: f64,
    pub flash: Spherical,
    pub velocity: Vec3,
    pub raw: [f64; 6],
}

impl<'a> Solver<'a> {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: &'a Data, params: &'a Params) -> Solver<'a> {
        Solver { data, params }
    }

    /// Find the solution
    pub fn solve(&self) -> Solved {
        let range = 10_000.0;
        let approx = self.monte_carlo(
            (self.data.mean_pos, self.data.mean_pos),
            self.params.range,
            range,
            0.6,
        );

        // Default bounds
        let def_min = [
            approx.0.x - range,
            approx.0.y - range,
            approx.0.z - range,
            approx.1.x - range,
            approx.1.y - range,
            approx.1.z - range,
        ];
        let def_max = [
            approx.0.x + range,
            approx.0.y + range,
            approx.0.z + range,
            approx.1.x + range,
            approx.1.y + range,
            approx.1.z + range,
        ];

        // Initial trajectory
        let mut traj = [
            approx.0.x, approx.0.y, approx.0.z, approx.1.x, approx.1.y, approx.1.z,
        ];

        // Binary search over 6 variables that
        // minimizes evaluate_traj function
        for _ in 0..self.params.repeat {
            for i in 0..6 {
                let mut min = def_min[i];
                let mut max = def_max[i];
                for _ in 0..self.params.depth {
                    let tmp = (min + max) * 0.5;
                    let change = (max - min) * 0.15;

                    traj[i] = tmp - change;
                    let e1 = self.evaluate_traj(
                        Vec3 {
                            x: traj[0],
                            y: traj[1],
                            z: traj[2],
                        },
                        Vec3 {
                            x: traj[3],
                            y: traj[4],
                            z: traj[5],
                        },
                    );

                    traj[i] = tmp + change;
                    let e2 = self.evaluate_traj(
                        Vec3 {
                            x: traj[0],
                            y: traj[1],
                            z: traj[2],
                        },
                        Vec3 {
                            x: traj[3],
                            y: traj[4],
                            z: traj[5],
                        },
                    );

                    if e1 < e2 {
                        max = tmp + change;
                    } else {
                        min = tmp - change;
                    }
                }
                traj[i] = (min + max) * 0.5;
            }
        }

        use std::io::Write;
        let mut file = std::fs::File::create("data_solved.dat").expect("create failed");
        for dx in -1000..1000 {
            let x = dx as f64 * 200.;
            file.write(
                format!(
                    "{} {}\n",
                    x / 1000.,
                    self.evaluate_traj(
                        Vec3 {
                            x: traj[0] + x,
                            y: traj[1],
                            z: traj[2],
                        },
                        Vec3 {
                            x: traj[3],
                            y: traj[4],
                            z: traj[5],
                        },
                    )
                )
                .as_bytes(),
            )
            .unwrap();
        }

        // p1 and p2 are now known
        let p1 = Vec3 {
            x: traj[0],
            y: traj[1],
            z: traj[2],
        };
        let p2 = Vec3 {
            x: traj[3],
            y: traj[4],
            z: traj[5],
        };

        // calculate the velocity
        let velocity = (p2 - p1).normalized();

        // calculate flash location as well as speed
        let (flash, speed) = self.calc_flash_and_speed(p1, velocity);

        // return
        Solved {
            error: (self.evaluate_traj(p1, p2)),
            flash,
            velocity: velocity * speed,
            raw: traj,
        }
    }

    fn monte_carlo(
        &self,
        mut answer: (Vec3, Vec3),
        range: f64,
        target_range: f64,
        k: f64,
    ) -> (Vec3, Vec3) {
        let def_min = [
            answer.0.x - range,
            answer.0.y - range,
            answer.0.z - range,
            answer.1.x - range,
            answer.1.z - range,
            answer.1.y - range,
        ];
        let def_max = [
            answer.0.x + range,
            answer.0.y + range,
            answer.0.z + range,
            answer.1.x + range,
            answer.1.z + range,
            answer.1.y + range,
        ];

        let mut error = std::f64::INFINITY;
        for _ in 0..1_000_000 {
            let ans = (
                Vec3 {
                    x: random::<f64>() * (def_max[0] - def_min[0]) + def_min[0],
                    y: random::<f64>() * (def_max[1] - def_min[1]) + def_min[1],
                    z: random::<f64>() * (def_max[2] - def_min[2]) + def_min[2],
                },
                Vec3 {
                    x: random::<f64>() * (def_max[3] - def_min[3]) + def_min[3],
                    y: random::<f64>() * (def_max[4] - def_min[4]) + def_min[4],
                    z: random::<f64>() * (def_max[5] - def_min[5]) + def_min[5],
                },
            );
            let err = self.evaluate_traj(ans.0, ans.1);
            if err < error {
                answer = ans;
                error = err;
            }
        }

        if range > target_range {
            self.monte_carlo(answer, range * k, target_range, k)
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

            if sample.trust_end && perpendic.dot(k_end) > self.params.min_match {
                let l_end = pip.dot(perpendic * (k_end.dot(v) / k_end.dot(perpendic)) - v);
                l_end_mean += l_end;
                l_count += 1.;
                if sample.trust_start && perpendic.dot(k_start) > self.params.min_match {
                    let l_start =
                        pip.dot(perpendic * (k_start.dot(v) / k_start.dot(perpendic)) - v);
                    speed += (l_end - l_start) / sample.duration;
                    speed_count += 1.;
                }
            }
        }
        l_end_mean /= l_count;
        speed /= speed_count;
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
            let perpendic = vel.cross(plane);

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let e_start = plane.dot(k_start).asin();
            let e_end = plane.dot(k_end).asin();

            //let trust_start = sample.trust_start;
            //let trust_end = sample.trust_end;
            let trust_start = sample.trust_start && perpendic.dot(k_start) > self.params.min_match;
            let trust_end = sample.trust_end && perpendic.dot(k_end) > self.params.min_match;

            if trust_start {
                error += e_start * e_start;
                count += 1.;

                if trust_end && sample.trust_da {
                    let start: Azimuthal = (k_start - plane * k_start.dot(plane))
                        .to_local(sample.geo_pos)
                        .into();
                    let end: Azimuthal = (k_end - plane * k_end.dot(plane))
                        .to_local(sample.geo_pos)
                        .into();

                    if let Some(da) = descent_angle(start, end) {
                        let diff = angle_diff(da, sample.descent_angle);
                        error += diff * diff * 0.5;
                        count += 1.;
                    }
                }
            }
            if trust_end {
                error += e_end * e_end;
                count += 1.;
            }
        }

        error / count
    }
}
