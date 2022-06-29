use rand::prelude::*;

use std::f64::consts::*;
use std::fmt::Write;

use fireball::constants::EARTH_R;
use fireball::data::{Answer, Data};
use fireball::maths::*;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

const IDEAL_DATA_TESTS: usize = 10;
const IDEAL_DATA_SAMPLES: usize = 100;

const NOT_SO_IDEAL_DATA_TESTS: usize = 5;
const NOT_SO_IDEAL_DATA_SAMPLES: usize = 100;

const TEN_PERCENT_RANDOM_TESTS: usize = 20;
const TEN_PERCENT_RANDOM_SAMPLES: usize = 100;

#[test]
fn ideal_data() {
    let mut done = 0usize;
    let mut rng = thread_rng();

    for _ in 0..IDEAL_DATA_TESTS {
        let mean_lat = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
        let mean_lon = rng.gen_range(-PI..PI);

        let (flash, vel) = gen_rand_traj(mean_lat, mean_lon, &mut rng);
        let (data, data_str) = gen_data(
            mean_lat,
            mean_lon,
            flash,
            vel,
            IDEAL_DATA_SAMPLES,
            0,
            None,
            None,
            &mut rng,
        );

        let mut solver = Solver::new(data, Params { da_correction: 0.0 });
        let solution = solver.solve();

        let flash_distance = (flash - solution.flash).len();
        let vel_angle = vel
            .normalized()
            .dot(solution.velocity.normalized())
            .min(1.0)
            .acos();
        let vel_diff = (vel.len() - solution.velocity.len()).abs();

        if flash_distance > 200. || vel_angle > 1f64.to_radians() || vel_diff > vel.len() * 0.03 {
            eprintln!("==== FAILED TESTCASE ==== ");
            eprintln!("Flash error = {}", flash_distance);
            eprintln!("Velocity angle = {}", vel_angle.to_degrees());
            eprintln!("Velocity rel error = {}", vel_diff / vel.len());
        } else {
            done += 1;
        }
    }

    if done < IDEAL_DATA_TESTS {
        panic!("{done}/{IDEAL_DATA_TESTS}");
    }
}

#[test]
fn not_so_ideal_data() {
    let mut done = 0usize;
    let mut rng = thread_rng();

    const K_SIGMA: f64 = 1.0;
    const DA_SIGMA: f64 = 1.0;

    for _ in 0..NOT_SO_IDEAL_DATA_TESTS {
        let mean_lat = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
        let mean_lon = rng.gen_range(-PI..PI);

        let (flash, vel) = gen_rand_traj(mean_lat, mean_lon, &mut rng);
        let (data, data_str) = gen_data(
            mean_lat,
            mean_lon,
            flash,
            vel,
            NOT_SO_IDEAL_DATA_SAMPLES,
            0,
            Some(K_SIGMA.to_radians()),
            Some(DA_SIGMA.to_radians()),
            &mut rng,
        );

        let mut solver = Solver::new(data, Params { da_correction: 0.0 });
        let solution = solver.solve();

        let flash_distance = (flash - solution.flash).len();
        let vel_angle = vel
            .normalized()
            .dot(solution.velocity.normalized())
            .min(1.0)
            .acos();
        let vel_diff = (vel.len() - solution.velocity.len()).abs();

        if flash_distance > 5_000. || vel_angle > 1f64.to_radians() || vel_diff > vel.len() * 0.03 {
            eprintln!("==== FAILED TESTCASE ==== ");
            eprintln!("Flash error = {}", flash_distance);
            eprintln!("Velocity angle = {}", vel_angle.to_degrees());
            eprintln!("Velocity rel error = {}", vel_diff / vel.len());
        } else {
            done += 1;
        }
    }

    if done < NOT_SO_IDEAL_DATA_TESTS {
        panic!("{done}/{NOT_SO_IDEAL_DATA_TESTS}");
    }
}

#[test]
fn ten_percent_random() {
    let mut done = 0usize;
    let mut rng = thread_rng();

    const K_SIGMA: f64 = 1.0;
    const DA_SIGMA: f64 = 1.0;

    for _ in 0..TEN_PERCENT_RANDOM_TESTS {
        let mean_lat = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
        let mean_lon = rng.gen_range(-PI..PI);

        let (flash, vel) = gen_rand_traj(mean_lat, mean_lon, &mut rng);
        let (data, data_str) = gen_data(
            mean_lat,
            mean_lon,
            flash,
            vel,
            TEN_PERCENT_RANDOM_SAMPLES,
            // 0, // 20 out of 20
            TEN_PERCENT_RANDOM_SAMPLES * 10 / 100, // 20 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 20 / 100, // 20 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 30 / 100, // 20 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 40 / 100, // 20 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 45 / 100, // 18 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 49 / 100, // 14 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 50 / 100, // 13 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 51 / 100, // 5 out of 20
            // TEN_PERCENT_RANDOM_SAMPLES * 55 / 100, // 0 out of 20
            Some(K_SIGMA.to_radians()),
            Some(DA_SIGMA.to_radians()),
            &mut rng,
        );

        let mut solver = Solver::new(data, Params { da_correction: 0.0 });
        let solution = solver.solve();

        let flash_distance = (flash - solution.flash).len();
        let vel_angle = vel
            .normalized()
            .dot(solution.velocity.normalized())
            .min(1.0)
            .acos();
        let vel_diff = (vel.len() - solution.velocity.len()).abs();

        if flash_distance > 5_000. || vel_angle > 1f64.to_radians() || vel_diff > vel.len() * 0.05 {
            eprintln!("==== FAILED TESTCASE ==== ");
            eprintln!("Flash error = {}", flash_distance);
            eprintln!("Velocity angle = {}", vel_angle.to_degrees());
            eprintln!("Velocity rel error = {:.1}%", vel_diff / vel.len() * 100.0);

            // Save failed testcase
            {
                use std::io::Write;
                let mut file =
                    std::fs::File::create(format!("failed_tests/test_{flash_distance:.2}.toml"))
                        .unwrap();
                file.write_all(data_str.as_bytes()).unwrap();
            }
        } else {
            done += 1;
        }
    }

    if done < TEN_PERCENT_RANDOM_TESTS {
        panic!("{done}/{TEN_PERCENT_RANDOM_TESTS}");
    }

    panic!("{done}/{TEN_PERCENT_RANDOM_TESTS}");
}

// generate some random data
fn gen_data(
    mean_lat: f64,
    mean_lon: f64,
    p2: Vec3,
    vel: Vec3,
    n_total: usize,
    n_random: usize,
    k_sigma: Option<f64>,
    da_sigma: Option<f64>,
    rng: &mut ThreadRng,
) -> (Data, String) {
    let mut data_str = String::new();

    for i in 0..n_total {
        let duration = rng.gen_range(5f64..30f64);
        let geo_pos = Spherical {
            lat: rng.gen_range((mean_lat - 3f64.to_radians())..(mean_lat + 3f64.to_radians())),
            lon: rng.gen_range((mean_lon - 3f64.to_radians())..(mean_lon + 3f64.to_radians())),
            r: rng.gen_range(EARTH_R..(EARTH_R + 2_000.)),
        };
        let global_pos: Vec3 = geo_pos.into();

        if i < n_random {
            writeln!(
            data_str,
            "[[sample]]\nlat = {}\n lon = {}\n h = {}\na = {}\n zb = {}\n hb = {}\n z0 = {}\n h0 = {}\n t = {}",
            geo_pos.lat.to_degrees(),
            geo_pos.lon.to_degrees(),
            geo_pos.r - EARTH_R,
            rng.gen_range(0f32..360f32),
            rng.gen_range(0f32..360f32),
            rng.gen_range(0f32..90f32),
            rng.gen_range(0f32..360f32),
            rng.gen_range(0f32..90f32),
            duration,
            ).unwrap();
        } else {
            let p1 = p2 - vel * duration;
            let k_start_golbal = p1 - global_pos;
            let k_end_golbal = p2 - global_pos;
            let mut k_start: Azimuthal = k_start_golbal.to_local(geo_pos).into();
            let mut k_end: Azimuthal = k_end_golbal.to_local(geo_pos).into();

            if let Some(k_sigma) = k_sigma {
                let dist = rand_distr::Normal::new(0.0, k_sigma).unwrap();
                k_start.h = (k_start.h + rng.sample(dist)).clamp(0.0, FRAC_PI_2);
                k_end.h = (k_end.h + rng.sample(dist)).clamp(0.0, FRAC_PI_2);
                k_start.z = (k_start.z + rng.sample(dist)).rem_euclid(TAU);
                k_end.z = (k_end.z + rng.sample(dist)).rem_euclid(TAU);
            }

            let mut da = descent_angle(
                global_pos,
                (k_start_golbal + k_end_golbal).normalized(),
                vel,
            );
            if let Some(da_sigma) = da_sigma {
                let dist = rand_distr::Normal::new(0.0, da_sigma).unwrap();
                da = (da + rng.sample(dist)).rem_euclid(TAU);
            }

            writeln!(
            data_str,
            "[[sample]]\nlat = {}\n lon = {}\n h = {}\na = {}\n zb = {}\n hb = {}\n z0 = {}\n h0 = {}\n t = {}",
            geo_pos.lat.to_degrees(),
            geo_pos.lon.to_degrees(),
            geo_pos.r - EARTH_R,
            da.to_degrees(),
            k_start.z.to_degrees(),
            k_start.h.to_degrees(),
            k_end.z.to_degrees(),
            k_end.h.to_degrees(),
            duration,
        ).unwrap();
        }
    }

    let p2_geo: Spherical = p2.into();
    writeln!(
        data_str,
        "[answer]\nlat={}\nlon={}\nh={}\nvx={}\nvy={}\nvz={}",
        p2_geo.lat.to_degrees(),
        p2_geo.lon.to_degrees(),
        (p2_geo.r - EARTH_R),
        vel.x * 1e-3,
        vel.y * 1e-3,
        vel.z * 1e-3,
    )
    .unwrap();

    (toml::from_str(&data_str).unwrap(), data_str)
}

// genrate random trajectory
fn gen_rand_traj(mean_lat: f64, mean_lon: f64, rng: &mut ThreadRng) -> (Vec3, Vec3) {
    let da = 3f64.to_radians();
    let p1 = Spherical {
        lat: rng.gen_range((mean_lat - da)..(mean_lat + da)),
        lon: rng.gen_range((mean_lon - da)..(mean_lon + da)),
        r: rng.gen_range((EARTH_R + 60_000.)..(EARTH_R + 90_000.)),
    };
    let p2 = Spherical {
        lat: rng.gen_range((mean_lat - da)..(mean_lat + da)),
        lon: rng.gen_range((mean_lon - da)..(mean_lon + da)),
        r: rng.gen_range((EARTH_R + 20_000.)..(EARTH_R + 50_000.)),
    };

    let flash: Vec3 = p2.into();
    let dir = (flash - Vec3::from(p1)).normalized() * rng.gen_range(10_000.0..30_000.0);

    (flash, dir)
}
