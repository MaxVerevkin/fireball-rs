//! The implementation itself

use crate::data::Data;
use crate::params::Params;

use crate::constants::*;
use crate::structs::*;

use std::f64::consts::{PI, TAU};

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
        // Default bounds
        // 0.08 ~ 4.58 degrees
        let def_min = [
            self.data.mean_lat - 0.08,
            self.data.mean_lon - 0.08,
            EARTH_R + 30_000.,
            self.data.mean_lat - 0.08,
            self.data.mean_lon - 0.08,
            EARTH_R + 10_000.,
        ];
        let def_max = [
            self.data.mean_lat + 0.08,
            self.data.mean_lon + 0.08,
            EARTH_R + 120_000.,
            self.data.mean_lat + 0.08,
            self.data.mean_lon + 0.08,
            EARTH_R + 90_000.,
        ];

        // Initial trajectory
        let mut traj = [
            self.data.mean_lat,
            self.data.mean_lon,
            EARTH_R + 75_000.,
            self.data.mean_lat,
            self.data.mean_lon,
            EARTH_R + 50_000.,
        ];

        // Bounds
        let mut min = def_min.clone();
        let mut max = def_max.clone();

        // Binary search over 6 variables that
        // minimizes evaluate_traj function
        let mut tmp;
        for _ in 0..self.params.repeat {
            for i in 0..6 {
                min[i] = def_min[i];
                max[i] = def_max[i];
                for _ in 0..self.params.depth {
                    tmp = (min[i] + max[i]) * 0.5;
                    let change = (max[i] - min[i]) * 0.15;

                    traj[i] = tmp - change;
                    let e1 = self.evaluate_traj(
                        &Spherical {
                            lat: traj[0],
                            lon: traj[1],
                            r: traj[2],
                        },
                        &Spherical {
                            lat: traj[3],
                            lon: traj[4],
                            r: traj[5],
                        },
                    );

                    traj[i] = tmp + change;
                    let e2 = self.evaluate_traj(
                        &Spherical {
                            lat: traj[0],
                            lon: traj[1],
                            r: traj[2],
                        },
                        &Spherical {
                            lat: traj[3],
                            lon: traj[4],
                            r: traj[5],
                        },
                    );

                    if e1 < e2 {
                        max[i] = tmp + change;
                    } else {
                        min[i] = tmp - change;
                    }
                }
                traj[i] = (min[i] + max[i]) * 0.5;
            }
        }

        // p1 and p2 are now known
        let p1 = Spherical {
            lat: traj[0],
            lon: traj[1],
            r: traj[2],
        };
        let p2 = Spherical {
            lat: traj[3],
            lon: traj[4],
            r: traj[5],
        };

        // calculate the velocity
        let point = p1.to_vec3();
        let velocity = (p2.to_vec3() - point).normalize();

        // FIXME
        eprintln!(
            "Distance from the 'true' answer: {} km",
            (point
                - (Spherical {
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
        let (flash, speed) = self.calc_flash_and_speed(&point, &velocity);

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
    fn evaluate_traj(&self, p1: &Spherical, p2: &Spherical) -> f64 {
        let start = p1.to_vec3();
        let end = p2.to_vec3();
        let velocity = (end - start).normalize();

        let mut error = 0.;
        let mut count = 0.;
        for sample in &self.data.samples {
            let plane = (end - sample.global_pos).cross(&velocity).normalize();
            let perpendic = velocity.cross(&plane);

            let e_start = plane.dot(&sample.global_start).asin();
            let e_end = plane.dot(&sample.global_end).asin();

            let trust_start =
                sample.trust_start && perpendic.dot(&sample.global_start) > self.params.min_match;
            let trust_end =
                sample.trust_end && perpendic.dot(&sample.global_end) > self.params.min_match;

            if trust_start {
                error += e_start * e_start;
                count += 1.;

                if trust_end && sample.trust_da {
                    let start = (sample.global_start - plane * sample.global_start.dot(&plane))
                        .to_local(&sample.geo_pos)
                        .to_azimuthal();
                    let end = (sample.global_end - plane * sample.global_end.dot(&plane))
                        .to_local(&sample.geo_pos)
                        .to_azimuthal();

                    let diff = angle_diff(descent_angle(&start, &end), sample.descent_angle);
                    error += diff * diff * 0.5;
                    count += 1.;
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

fn angle_diff(a1: f64, a2: f64) -> f64 {
    (a2 - a1 + PI).rem_euclid(TAU) - PI
}

fn descent_angle(start: &Azimuthal, end: &Azimuthal) -> f64 {
    let dz = angle_diff(start.z, end.z);

    let cos_l = start.h.sin() * end.h.sin() + start.h.cos() * end.h.cos() * dz.cos();
    let sin_l = (1. - cos_l * cos_l).sqrt();

    let a = ((end.h.sin() - start.h.sin() * cos_l) / (start.h.cos() * sin_l)).acos();
    if a.is_nan() {
        return 0.;
    }

    if dz > 0. {
        a
    } else {
        TAU - a
    }
}

#[test]
fn angle_diff_test() {
    assert_eq!(angle_diff(PI, TAU).abs(), PI);
    assert_eq!(angle_diff(PI, 10. * TAU).abs(), PI);
    assert_eq!(angle_diff(PI, 0.).abs(), PI);

    assert!(angle_diff(PI + 0.1, TAU) - (PI - 0.1) < 1e-10);
    assert!(angle_diff(PI - 0.1, 0.) - (-PI + 0.1) < 1e-10);

    assert!(angle_diff(PI - 1., -0.1) - (-1.1) < 1e-10);
    assert!(angle_diff(PI + 1., TAU - 0.1) - (PI - 1.1) < 1e-10);
}

#[test]
fn descent_angle_test() {
    // basic vertical
    assert_eq!(
        descent_angle(&Azimuthal { z: 1., h: 0.1 }, &Azimuthal { z: 1., h: 0.2 }),
        0.,
    );
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 1., h: 0.2 }) - (PI)) < 1e-5,
    );

    // basic horizontal
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 1.01, h: 0.3 }) - (PI / 2.))
            < 1e-5,
    );
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 0.99, h: 0.3 })
            - (3. * PI / 2.))
            < 1e-2,
    );

    // different quarters
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 1.01, h: 1.01 }) - (PI / 4.))
            < 1e-2,
    );
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 1.01, h: 0.99 })
            - (PI / 4. * 3.))
            < 0.3,
    );
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 0.99, h: 0.99 })
            - (PI / 4. * 5.))
            < 0.3,
    );
    assert!(
        (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 0.99, h: 1.01 })
            - (PI / 4. * 7.))
            < 0.3,
    );
}
