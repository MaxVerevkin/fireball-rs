//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! Run:
//! ```bash
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.

use clap::*;

use common::obs_data::Data;
use common::structs::*;
use fireball::solver::{Params, Solver};

fn main() {
    let matches = app_from_crate!()
        .arg(
            Arg::with_name("file")
                .help("The file with data")
                .takes_value(true)
                .required(true)
        )
        .arg(
            Arg::with_name("dac")
                .help("Whether to tilt the observation to better match given descent angle (0 means no, 1 means yes")
                .default_value("0.5")
                .takes_value(true)
                .long("dac")
        )
        .arg(
            Arg::with_name("dac-max")
                .help("Max correction (in degrees) that DAC is allowed to apply")
                .default_value("20")
                .takes_value(true)
                .long("dac-max")
        )
        .arg(
            Arg::with_name("da-only")
                .required(false)
                .long("da-only")
        )
        .get_matches();

    let file_name = matches.value_of("file").unwrap();
    let da_correction: f64 = matches.value_of("dac").unwrap().parse().unwrap();
    let dac_max = matches
        .value_of("dac-max")
        .unwrap()
        .parse::<f64>()
        .unwrap()
        .to_radians();
    let da_only = matches.is_present("da-only");

    // Read data
    let mut data = Data::read_from_toml(file_name);
    data.apply_da_correction(0.064); // ../../da-corrector/log.txt : 265

    // Solve!
    let mut solver = Solver::new(
        data,
        Params {
            da_correction,
            dac_max,
            da_only,
        },
    );
    let solution = solver.solve();
    let flash: Spherical = solution.flash.into();

    // Print the answer
    println!("Flash location: {}", flash);
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solution.velocity.norm() / 1000.,
        solution.velocity.x / 1000.,
        solution.velocity.y / 1000.,
        solution.velocity.z / 1000.,
    );
}
