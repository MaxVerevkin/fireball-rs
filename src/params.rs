//! Use clap crate to parse command line arguments
//!
//! This program uses command line arguments to tweak some implementation details.

use clap::*;

/// Represents the command line arguments
pub struct Params {
    pub file_name: String,
    pub depth: u32,
    pub repeat: u32,
    pub min_match: f64,
}

impl Params {
    /// Automaticaly parse arguments using clap
    pub fn new() -> Self {
        let matches = app_from_crate!()
            .arg(
                Arg::with_name("file")
                    .help("Specify the file with data")
                    .takes_value(true)
                    .required(true),
            )
            .arg(
                Arg::with_name("depth")
                    .help("The depth of binary search")
                    .default_value("16")
                    .takes_value(true)
                    .long("depth")
                    .short("d"),
            )
            .arg(
                Arg::with_name("repeat")
                    .help("Number of repetitions of binary search")
                    .default_value("50")
                    .takes_value(true)
                    .long("repeat")
                    .short("r"),
            )
            .arg(
                Arg::with_name("min_match")
                    .help("Number that represents how accurate the observation should be to be included")
                    .default_value("0.3")
                    .takes_value(true)
                    .long("min-match")
                    .short("m"),
            )
            .get_matches();

        Params {
            file_name: matches.value_of("file").unwrap().to_owned(),
            depth: matches.value_of("depth").unwrap().parse().unwrap(),
            repeat: matches.value_of("repeat").unwrap().parse().unwrap(),
            min_match: matches.value_of("min_match").unwrap().parse().unwrap(),
        }
    }
}
