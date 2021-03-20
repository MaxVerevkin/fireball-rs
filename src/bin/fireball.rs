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

fn main() {
    let matches = app_from_crate!()
        .arg(
            Arg::with_name("file")
                .help("Specify the file with data")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("range")
                .help("The range of values to search in (km)")
                .default_value("700000.0")
                .takes_value(true)
                .long("range"),
        )
        .arg(
            Arg::with_name("min_match")
                .help(
                    "Number that represents how accurate the observation should be to be included",
                )
                .default_value("0.0")
                .takes_value(true)
                .long("min-match"),
        )
        .get_matches();

    let file_name = matches.value_of("file").unwrap();
    let range: f64 = matches
        .value_of("range")
        .unwrap()
        .parse()
        .expect("failed to parse `range`");
    let min_match: f64 = matches
        .value_of("min_match")
        .unwrap()
        .parse()
        .expect("failed to parse `min_match`");

    // Create new data structure
    let data = match Data::from_file(file_name) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("error: {}", e);
            std::process::exit(1);
        }
    };

    // Solve!
    let solver = Solver::new(data, Params { range, min_match });
    let solution = solver.solve();

    // Print the answer
    // Note: \u{00b0} -- the degree symbol
    println!("Error: {:.6}", solution.error);
    println!(
        "Flash location:\n  lat: {:.3}\u{00b0}\n  lon: {:.3}\u{00b0}\n  h: {:.3} km",
        solution.flash.lat.to_degrees(),
        solution.flash.lon.to_degrees(),
        (solution.flash.r - EARTH_R) / 1000.
    );
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solution.velocity.length().abs() / 1000.,
        solution.velocity.x / 1000.,
        solution.velocity.y / 1000.,
        solution.velocity.z / 1000.,
    );
}
