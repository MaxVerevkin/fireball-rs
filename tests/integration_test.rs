use rand::prelude::*;

use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt::Write;

use fireball::constants::EARTH_R;
use fireball::data::Data;
use fireball::maths::*;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

#[test]
fn ideal_data() {
    const TESTS: usize = 10;
    let mut done = 0usize;
    let mut rng = thread_rng();

    for _ in 0..TESTS {
        let mean_lat = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
        let mean_lon = rng.gen_range(-PI..PI);

        let (flash, vel) = gen_rand_traj(mean_lat, mean_lon, &mut rng);
        let data = gen_data(mean_lat, mean_lon, flash, vel, 50);

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

    if done < TESTS {
        panic!("{done}/{TESTS}");
    }
}

// generate some random data
fn gen_data(mean_lat: f64, mean_lon: f64, p2: Vec3, vel: Vec3, n: u32) -> Data {
    let mut data_str = String::new();

    for _ in 0..n {
        let duration = 1. + random::<f64>() * 3.;
        let geo_pos = Spherical {
            lat: mean_lat - (3f64).to_radians() + random::<f64>() * (6f64).to_radians(),
            lon: mean_lon - (3f64).to_radians() + random::<f64>() * (6f64).to_radians(),
            r: EARTH_R + random::<f64>() * 2_000.,
        };
        let global_pos: Vec3 = geo_pos.into();

        let p1 = p2 - vel * duration;

        let k_start_golbal = p1 - global_pos;
        let k_end_golbal = p2 - global_pos;

        let k_start: Azimuthal = k_start_golbal.to_local(geo_pos).into();
        let k_end: Azimuthal = k_end_golbal.to_local(geo_pos).into();

        writeln!(
            data_str,
            "[[sample]]\nlat = {}\n lon = {}\n h = {}\na = {}\n zb = {}\n hb = {}\n z0 = {}\n h0 = {}\n t = {}",
            geo_pos.lat.to_degrees(),
            geo_pos.lon.to_degrees(),
            geo_pos.r - EARTH_R,
            descent_angle(global_pos, (k_start_golbal + k_end_golbal).normalized(), vel).to_degrees(),
            k_start.z.to_degrees(),
            k_start.h.to_degrees(),
            k_end.z.to_degrees(),
            k_end.h.to_degrees(),
            duration,
        )
        .unwrap();
    }

    toml::from_str(&data_str).unwrap()
}

// genrate random trajectory
fn gen_rand_traj(mean_lat: f64, mean_lon: f64, rng: &mut ThreadRng) -> (Vec3, Vec3) {
    let da = 5f64.to_radians();
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

    (
        p2.into(),
        (Vec3::from(p2) - Vec3::from(p1)).normalized() * rng.gen_range(5_000.0..30_000.0),
    )
}
