//! The implementation itself

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
}

impl<'a> Solver<'a> {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: &'a Data, params: &'a Params) -> Solver<'a> {
        Solver { data, params }
    }

    /// Solution
    pub fn solve(&self) -> Solved {
        let mid = Spherical {
            lat: self.data.mean_lat,
            lon: self.data.mean_lon,
            r: EARTH_R,
        }
        .to_vec3();

        // Default bounds
        let def_min = [
            mid.x - 700_000.,
            mid.y - 700_000.,
            mid.z - 700_000.,
            mid.x - 700_000.,
            mid.y - 700_000.,
            mid.z - 700_000.,
        ];
        let def_max = [
            mid.x + 700_000.,
            mid.y + 700_000.,
            mid.z + 700_000.,
            mid.x + 700_000.,
            mid.y + 700_000.,
            mid.z + 700_000.,
        ];

        // Initial trajectory
        let mut traj = [mid.x, mid.y, mid.z, mid.x, mid.y, mid.z + 1000.];

        // Binary search over 6 variables that
        // minimizes evaluate_traj function
        let mut tmp;
        for _ in 0..self.params.repeat {
            for i in 0..6 {
                let mut min = def_min[i];
                let mut max = def_max[i];
                for _ in 0..self.params.depth {
                    tmp = (min + max) * 0.5;
                    let change = (max - min) * 0.15;

                    traj[i] = tmp - change;
                    let e1 = self.evaluate_traj(
                        &Vec3 {
                            x: traj[0],
                            y: traj[1],
                            z: traj[2],
                        },
                        &Vec3 {
                            x: traj[3],
                            y: traj[4],
                            z: traj[5],
                        },
                    );

                    traj[i] = tmp + change;
                    let e2 = self.evaluate_traj(
                        &Vec3 {
                            x: traj[0],
                            y: traj[1],
                            z: traj[2],
                        },
                        &Vec3 {
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
        let velocity = (p2 - p1).normalize();

        // FIXME
        eprintln!(
            "Distance from the 'true' answer: {} km",
            (p1 - (Spherical {
                lat: 33.1f64.to_radians(),
                lon: 34.3f64.to_radians(),
                r: EARTH_R + 43_300.
            })
            .to_vec3())
            .cross(&velocity)
            .length()
                / 1000.
        );

        // calculate flash location as well as speed
        let (flash, speed) = self.calc_flash_and_speed(&p1, &velocity);

        // return
        Solved {
            error: (self.evaluate_traj(&p1, &p2)),
            flash,
            velocity: velocity * speed,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, point: &Vec3, v: &Vec3) -> (Spherical, f64) {
        let mut l_end_mean = 0.;
        let mut speed = 0.;
        let mut l_count = 0.;
        let mut speed_count = 0.;
        for sample in &self.data.samples {
            let pip = *point - sample.global_pos;
            let plane = pip.cross(&v).normalize();
            let perpendic = v.cross(&plane).normalize();

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            if sample.trust_end && perpendic.dot(&k_end) > self.params.min_match {
                let l_end = pip.dot(&(perpendic * (k_end.dot(&v) / k_end.dot(&perpendic)) - *v));
                l_end_mean += l_end;
                l_count += 1.;
                if sample.trust_start && perpendic.dot(&k_start) > self.params.min_match {
                    let l_start =
                        pip.dot(&(perpendic * (k_start.dot(&v) / k_start.dot(&perpendic)) - *v));
                    speed += (l_end - l_start) / sample.duration;
                    speed_count += 1.;
                }
            }
        }
        l_end_mean /= l_count;
        speed /= speed_count;
        ((*point + *v * l_end_mean).to_spherical(), speed)
    }

    /// Calculate the mean of squared errors (less is better)
    pub fn evaluate_traj(&self, p1: &Vec3, p2: &Vec3) -> f64 {
        let vel = (*p1 - *p2).normalize();
        let mut error = 0.;
        let mut count = 0.;
        for sample in &self.data.samples {
            let plane = (*p1 - sample.global_pos).cross(&vel).normalize();
            let perpendic = vel.cross(&plane);

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let e_start = plane.dot(&k_start).asin();
            let e_end = plane.dot(&k_end).asin();

            let trust_start = sample.trust_start && perpendic.dot(&k_start) > self.params.min_match;
            let trust_end = sample.trust_end && perpendic.dot(&k_end) > self.params.min_match;

            if trust_start {
                error += e_start * e_start;
                count += 1.;

                if trust_end && sample.trust_da {
                    let start = (k_start - plane * k_start.dot(&plane))
                        .to_local(&sample.geo_pos)
                        .to_azimuthal();
                    let end = (k_end - plane * k_end.dot(&plane))
                        .to_local(&sample.geo_pos)
                        .to_azimuthal();

                    if let Some(da) = descent_angle(&start, &end) {
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
