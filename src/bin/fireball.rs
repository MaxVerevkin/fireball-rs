//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! Run:
//! ```bash
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.

use clap::*;

use fireball::constants::EARTH_R;
use fireball::data::Data;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

fn main() {
    let matches = app_from_crate!()
        .arg(
            Arg::with_name("file")
                .help("Specify the file with data")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("threads")
                .help("The number of threads this program will spawn")
                .default_value("1")
                .takes_value(true)
                .long("threads")
                .short("j"),
        )
        .get_matches();

    let file_name = matches.value_of("file").unwrap();
    let threads: usize = matches
        .value_of("threads")
        .unwrap()
        .parse()
        .expect("failed to parse `threads`");

    // Create new data structure
    let data = match Data::from_file(file_name) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("error: {}", e);
            std::process::exit(1);
        }
    };

    // Solve!
    let mut solver = Solver::new(data, Params { threads });
    let solution = solver.solve();
    let flash: Spherical = solution.flash.into();

    // Print the answer
    // Note: \u{00b0} -- the degree symbol
    // println!("Error: {:.6}", solution.error);
    println!(
        "Flash location:\n  lat: {:.3}\u{00b0}\n  lon: {:.3}\u{00b0}\n  h: {:.3} km",
        flash.lat.to_degrees(),
        flash.lon.to_degrees(),
        (flash.r - EARTH_R) / 1000.
    );
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solution.velocity.len().abs() / 1000.,
        solution.velocity.x / 1000.,
        solution.velocity.y / 1000.,
        solution.velocity.z / 1000.,
    );

    // use fireball::structs::*;
    // let v = Vec3::new(1., 1., 2.);
    // let axis = Vec3::new(0., 1., 0.);
    // let angle = 180f64;
    // dbg!(v);
    // dbg!(axis);
    // dbg!(angle);
    // let q = RotationQuaternion::new(axis, angle.to_radians());
    // dbg!(q);
    // dbg!(q.apply_rotation(v));
}
