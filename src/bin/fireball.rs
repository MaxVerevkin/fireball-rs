//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! Run:
//! ```bash
//! $ fireball -h
//! ```
//! to see which parameters of the implementation can be tweaked.

use std::path::PathBuf;

use clap::Parser;
use common::obs_data::RawData;
use fireball::solver::{Params, Solver};

#[derive(Debug, Parser)]
struct CliArgs {
    /// The file with data
    file: PathBuf,

    /// Use only descent angles in the trajectory evaluation
    #[arg(long)]
    correct_altitudes: bool,

    /// Do not 'flip' the decent angles
    #[arg(long)]
    no_da_flip: bool,

    /// Ignore altitudes
    #[arg(long)]
    no_altitudes: bool,

    /// Ignore azimuths
    #[arg(long)]
    no_azimuths: bool,

    /// Parameter of descent angle correction function da_corrected = da - k * sin(da * 2.0)
    #[arg(long, default_value_t = 0.0)]
    da_k: f64,

    /// Parameter of azimuth correction function az_corrected = az + k * sin(az)
    #[arg(long, default_value_t = 0.0)]
    az_k: f64,
}

fn main() {
    let args = CliArgs::parse();

    // Read data
    let mut data = RawData::read_from_toml_file(&args.file)
        .da_correction(args.da_k)
        .az_correction(args.az_k);
    if args.correct_altitudes {
        data = data.altitudes_correction();
    }
    let data = data.finalize();

    // Solve!
    let mut solver = Solver::new(
        data,
        Params {
            no_da_flip: args.no_da_flip,
            no_altitudes: args.no_altitudes,
            no_azimuths: args.no_azimuths,
            correct_altitudes: args.correct_altitudes,
            da_k: args.da_k,
            az_k: args.az_k,
        },
    );
    let solution = solver.solve();

    // Print the answer
    println!("Flash location: {}", solution.flash);
    println!(
        "Velocity: {:.3} km/s\n  x: {:.3} km/s\n  y: {:.3} km/s\n  z: {:.3} km/s",
        solution.velocity.norm() / 1000.,
        solution.velocity.x / 1000.,
        solution.velocity.y / 1000.,
        solution.velocity.z / 1000.,
    );
}
