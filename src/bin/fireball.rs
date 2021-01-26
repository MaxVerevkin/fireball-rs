//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! Run:
//! ```bash
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.

use fireball::constants::EARTH_R;
use fireball::data::Data;
use fireball::params::Params;
use fireball::solver::Solver;

fn main() {
    // Parse argruments
    let params = match Params::from_cmd() {
        Some(p) => p,
        None => {
            eprintln!("error: Could not initialize command line arguments");
            std::process::exit(1);
        }
    };

    // Create new data structure
    let data = match Data::from_file(&params.file_name) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("error: {}", e);
            std::process::exit(2);
        }
    };

    // Solve!
    let solved = Solver::new(&data, &params).solve();

    // Print the answer
    // Note: \u{00b0} -- the degree symbol
    println!("Error: {:.6}", solved.error);
    println!(
        "Flash location:\n  lat: {:.3}\u{00b0}\n  lon: {:.3}\u{00b0}\n  h: {:.3} km",
        solved.flash.lat.to_degrees(),
        solved.flash.lon.to_degrees(),
        (solved.flash.r - EARTH_R) / 1000.
    );
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solved.velocity.length().abs() / 1000.,
        solved.velocity.x / 1000.,
        solved.velocity.y / 1000.,
        solved.velocity.z / 1000.,
    );
}
