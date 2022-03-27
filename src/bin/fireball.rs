//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! Run:
//! ```bash
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.

use std::fs::File;
use std::io::{BufReader, Read};

use clap::*;

use fireball::data::Data;
use fireball::solver::{Params, Solver};
use fireball::structs::*;

fn main() {
    let matches = app_from_crate!()
        .arg(
            Arg::with_name("file")
                .help("The file with data")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("dac")
                .help("Whether to tilt the observation to better match given descent angle (0 means no, 1 means yes")
                .default_value("0.5")
                .takes_value(true)
                .long("dac"),
        )
        .get_matches();

    let file_name = matches.value_of("file").unwrap();
    let da_correction: f64 = matches.value_of("dac").unwrap().parse().unwrap();

    // Read data
    let mut file = BufReader::new(File::open(file_name).expect("Failed to open file."));
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)
        .expect("Failed to read from data file");
    let data: Data =
        toml::from_str(&String::from_utf8_lossy(&buf)).expect("Failed to deserialize data");

    // Solve!
    let mut solver = Solver::new(data, Params { da_correction });
    let solution = solver.solve();
    let flash: Spherical = solution.flash.into();

    // Print the answer
    // println!("Error: {:.6}", solution.error);
    println!("Flash location: {}", flash);
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solution.velocity.len().abs() / 1000.,
        solution.velocity.x / 1000.,
        solution.velocity.y / 1000.,
        solution.velocity.z / 1000.,
    );
}
