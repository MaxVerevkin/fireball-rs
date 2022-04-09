use rand::prelude::*;

use std::f64::consts::*;
use std::fmt::Write;

use fireball::constants::EARTH_R;
use fireball::data::{Answer, Data};
use fireball::maths::*;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

const IDEAL_DATA_TESTS: usize = 20;
const IDEAL_DATA_SAMPLES: usize = 20;

const NOT_SO_IDEAL_DATA_TESTS: usize = 5;
const NOT_SO_IDEAL_DATA_SAMPLES: usize = 100;

#[test]
fn ideal_data() {
    let mut done = 0usize;
    let mut rng = thread_rng();

    for _ in 0..IDEAL_DATA_TESTS {
        let mean_lat = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
        let mean_lon = rng.gen_range(-PI..PI);

        let (flash, vel) = gen_rand_traj(mean_lat, mean_lon, &mut rng);
        let data = gen_data(
            mean_lat,
            mean_lon,
            flash,
            vel,
            IDEAL_DATA_SAMPLES,
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
        let data = gen_data(
            mean_lat,
            mean_lon,
            flash,
            vel,
            NOT_SO_IDEAL_DATA_SAMPLES,
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

        if flash_distance > 200. || vel_angle > 1f64.to_radians() || vel_diff > vel.len() * 0.03 {
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

// generate some random data
fn gen_data(
    mean_lat: f64,
    mean_lon: f64,
    p2: Vec3,
    vel: Vec3,
    n: usize,
    k_sigma: Option<f64>,
    da_sigma: Option<f64>,
    rng: &mut ThreadRng,
) -> Data {
    let mut data_str = String::new();

    for _ in 0..n {
        let duration = rng.gen_range(2f64..5f64);
        let geo_pos = Spherical {
            lat: rng.gen_range((mean_lat - 3f64.to_radians())..(mean_lat + 3f64.to_radians())),
            lon: rng.gen_range((mean_lon - 3f64.to_radians())..(mean_lon + 3f64.to_radians())),
            r: rng.gen_range(EARTH_R..(EARTH_R + 2_000.)),
        };
        let global_pos: Vec3 = geo_pos.into();

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
        )
        .unwrap();
    }

    let mut data: Data = toml::from_str(&data_str).unwrap();
    data.answer = Some(Answer(Line {
        point: p2,
        direction: vel,
    }));
    data
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
    let dir = (flash - Vec3::from(p1)).normalized() * rng.gen_range(5_000.0..30_000.0);

    (flash, dir)
}
