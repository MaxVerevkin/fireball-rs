//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! The current implementation uses binary search to find two points
//! in the space that represents the trajectory of a fireball. The binary search
//! itself minimizes the mean-square-error of observations given trajectory.
//!
//! Run:
//! ```
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.
//!
//! TODO Write tets and some more documentation

mod constants;
mod data;
mod maths;
mod params;
mod solver;
mod structs;

use constants::EARTH_R;
use data::Data;
use params::Params;
use solver::Solver;

fn main() {
    // Parse argruments
    let params = Params::new();

    // Create new data structure
    let data = match Data::from_file(&params.file_name) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("error: {}", e);
            std::process::exit(1);
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
