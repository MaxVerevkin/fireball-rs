use rand::prelude::*;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use fireball::constants::EARTH_R;
use fireball::data::{Data, DataSample};
use fireball::maths::descent_angle;
use fireball::params::Params;
use fireball::solver::Solver;
use fireball::structs::*;

#[test]
fn test() {
    let tests: usize = 1;
    let must_pass: usize = 1;
    let mut passed: usize = 0;
    for _ in 0..tests {
        // (-pi/2, pi/2)
        let mean_lat = -FRAC_PI_2 + random::<f64>() * PI;
        // (-pi, pi)
        let mean_lon = -PI + random::<f64>() * TAU;

        // gnenerate trajectory
        let (flash, vel) = gen_traj(mean_lat, mean_lon);
        let data = gen_data(mean_lat, mean_lon, &flash, &vel, 300);

        let params = Params {
            file_name: String::new(),
            range: 500_000.,
            depth: 25,
            repeat: 500,
            min_match: 0.0,
        };

        let solver = Solver::new(&data, &params);
        let solved = solver.solve();

        // Draw a plot for generated data.
        use std::io::Write;
        let xyz1 = flash.to_vec3();
        let xyz2 = xyz1 + vel * 20_000.;
        let xyz1 = [xyz1.x, xyz1.y, xyz1.z];
        let mut file = std::fs::File::create("data_real.dat").expect("create failed");
        for dx in -1000..1000 {
            let x = dx as f64 * 200.;
            file.write(
                format!(
                    "{} {}\n",
                    x / 1000.,
                    solver.evaluate_traj(
                        &Vec3 {
                            x: xyz1[0] + x,
                            y: xyz1[1],
                            z: xyz1[2],
                        },
                        &xyz2
                    )
                )
                .as_bytes(),
            )
            .unwrap();
        }

        dbg!(solved.error);
        dbg!(solved.raw);

        if ((flash.r - solved.flash.r).abs() < 1000.)
            && ((flash.lat - solved.flash.lat).abs() < 0.01)
            && ((flash.lon - solved.flash.lon).abs() < 0.01)
            && ((vel.x - solved.velocity.x).abs() < 500.)
            && ((vel.y - solved.velocity.y).abs() < 500.)
            && ((vel.z - solved.velocity.z).abs() < 500.)
        {
            passed += 1;
        } else {
            dbg!(flash);
            dbg!(solved.flash);
            dbg!(vel);
            dbg!(solved.velocity);

            let flash = flash.to_vec3();
            dbg!(solver.evaluate_traj(&flash, &(flash - vel)));
        }
    }
    dbg!(passed);
    assert!(passed > must_pass);
}

// generate some random data
fn gen_data(mean_lat: f64, mean_lon: f64, flash: &Spherical, vel: &Vec3, n: u32) -> Data {
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
            lat: mean_lat - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
            lon: mean_lon - (5f64).to_radians() + random::<f64>() * (10f64).to_radians(),
            r: EARTH_R + random::<f64>() * 2_000.,
        };
        let global_pos = geo_pos.to_vec3();

        let p1 = p2 - *vel * duration;

        let k_start_golbal = p1 - global_pos;
        let k_end_golbal = p2 - global_pos;

        let k_start = k_start_golbal.to_local(&geo_pos).to_azimuthal();
        let k_end = k_end_golbal.to_local(&geo_pos).to_azimuthal();

        data.samples.push(
            DataSample::from_text(&format!(
                "{} {} {} {} {} {} {} {} {}",
                geo_pos.lat.to_degrees(),
                geo_pos.lon.to_degrees(),
                geo_pos.r - EARTH_R,
                descent_angle(&k_start, &k_end).unwrap().to_degrees(),
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
        (p2.to_vec3() - p1.to_vec3()).normalize() * (5_000. + random::<f64>() * 20_000.),
    )
}

// Compute azimuth given observer location and point which he observes
//fn azimuth(observer: &Spherical, point: &Spherical) -> Option<f64> {
//descent_angle(
//&Azimuthal {
//z: observer.lon,
//h: observer.lat,
//},
//&Azimuthal {
//z: point.lon,
//h: point.lat,
//},
//)
//}
