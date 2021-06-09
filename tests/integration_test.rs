use rand::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use fireball::constants::EARTH_R;
use fireball::data::{Data, DataSample};
use fireball::maths::descent_angle;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

#[test]
fn test() {
    let tests: usize = 100;
    let must_pass: usize = 95;
    let mut passed: usize = 0;
    for _ in 0..tests {
        // (-pi/2, pi/2)
        let mean_lat = -FRAC_PI_2 + random::<f64>() * PI;
        // (-pi, pi)
        let mean_lon = -PI + random::<f64>() * TAU;

        // gnenerate trajectory
        let (flash, vel) = gen_traj(mean_lat, mean_lon);
        let data = gen_data(mean_lat, mean_lon, &flash, vel, 200);
        let flash = flash.to_vec3();

        let params = Params {
            initial_range: 500_000.,
            initial_iterations: 10_000,
            main_iterations: 10_000,
            threads: 4,
        };

        let solver = Solver::new(data, params);
        let solution = solver.solve();

        // Draw a plot for generated data.
        //let plot = || {
        //use std::{fs::File, io::Write};
        //let point = flash;
        //let mut offset = Vec3::default();
        //offset.x = -1_000_000.;
        //let mut file = File::create("data_real.dat").unwrap();
        //for _ in 0..2_000 {
        //let point = point + offset;
        //write!(
        //file,
        //"{} {}\n",
        //offset.x / 1_000.,
        //solver.evaluate_traj(point, vel.normalized())
        //)
        //.unwrap();
        //offset.x += 1_000.;
        //}
        //// Draw a plot for solution.
        //let point = solution.flash;
        //let mut offset = Vec3::default();
        //offset.x = -1_000_000.;
        //let mut file = File::create("data_sol.dat").unwrap();
        //for _ in 0..2_000 {
        //let point = point + offset;
        //write!(
        //file,
        //"{} {}\n",
        //offset.x / 1_000.,
        //solver.evaluate_traj(point, solution.velocity.normalized())
        //)
        //.unwrap();
        //offset.x += 1_000.;
        //}
        //};

        let flash_distance = (flash - solution.flash).length();
        let vel_angle = vel.normalized().dot(solution.velocity.normalized()).acos();
        let vel_diff = (vel.length() - solution.velocity.length()).abs();

        let max_vel_diff = vel.length() * 0.1;

        if flash_distance < 20_000. && vel_angle < 0.17 && vel_diff < max_vel_diff {
            passed += 1;
        } else {
            eprintln!("====FAILED TESTCASE====");
            eprintln!("Solution error = {}", solution.error);
            eprintln!("Flash error = {}", flash_distance);
            eprintln!("Velocity angle = {}", vel_angle.to_degrees());
            eprintln!("Velocity rel error = {}", vel_diff / vel.length());
            eprintln!();

            //plot();
            //dbg!(passed);
            //panic!();
        }
    }
    dbg!(passed);
    assert!(passed >= must_pass);
}

// generate some random data
fn gen_data(mean_lat: f64, mean_lon: f64, flash: &Spherical, vel: Vec3, n: u32) -> Data {
    let p2 = flash.to_vec3();
    let mut data = Data {
        samples: Vec::new(),
        mean_pos: Spherical {
            lat: mean_lat,
            lon: mean_lon,
            r: EARTH_R,
        }
        .to_vec3(),
    };

    for _ in 0..n {
        let duration = 1. + random::<f64>() * 3.;
        let geo_pos = Spherical {
            lat: mean_lat - (3f64).to_radians() + random::<f64>() * (6f64).to_radians(),
            lon: mean_lon - (3f64).to_radians() + random::<f64>() * (6f64).to_radians(),
            r: EARTH_R + random::<f64>() * 2_000.,
        };
        let global_pos = geo_pos.to_vec3();

        let p1 = p2 - vel * duration;

        let k_start_golbal = p1 - global_pos;
        let k_end_golbal = p2 - global_pos;

        let k_start: Azimuthal = k_start_golbal.to_local(geo_pos).into();
        let k_end: Azimuthal = k_end_golbal.to_local(geo_pos).into();

        data.samples.push(
            DataSample::from_text(&format!(
                "{} {} {} {} {} {} {} {} {}",
                geo_pos.lat.to_degrees(),
                geo_pos.lon.to_degrees(),
                geo_pos.r - EARTH_R,
                descent_angle(global_pos, k_start_golbal, vel).to_degrees(),
                k_start.z.to_degrees(),
                k_start.h.to_degrees(),
                k_end.z.to_degrees(),
                k_end.h.to_degrees(),
                duration
            ))
            .unwrap(),
        );
    }

    data
}

// genrate some random trajectory
fn gen_traj(mean_lat: f64, mean_lon: f64) -> (Spherical, Vec3) {
    let p1 = Spherical {
        lat: mean_lat - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
        lon: mean_lon - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
        r: EARTH_R + 80_000. + random::<f64>() * 30_000.,
    };
    let p2 = Spherical {
        lat: mean_lat - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
        lon: mean_lon - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
        r: EARTH_R + 20_000. + random::<f64>() * 30_000.,
    };

    (
        p2,
        (p2.to_vec3() - p1.to_vec3()).normalized() * (5_000. + random::<f64>() * 20_000.),
    )
}
