//! Use clap crate to parse command line arguments
//!
//! This program uses command line arguments to tweak some implementation details.

use clap::*;

/// Represents the command line arguments
pub struct Params {
    pub file_name: String,
    pub range: f64,
    pub depth: u32,
    pub repeat: u32,
    pub min_match: f64,
}

impl Params {
    /// Automaticaly parse arguments using clap
    pub fn from_cmd() -> Option<Self> {
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
                    .long("range")
            )
            .arg(
                Arg::with_name("depth")
                    .help("The depth of binary search")
                    .default_value("20")
                    .takes_value(true)
                    .long("depth")
            )
            .arg(
                Arg::with_name("repeat")
                    .help("Number of repetitions of binary search")
                    .default_value("50")
                    .takes_value(true)
                    .long("repeat")
            )
            .arg(
                Arg::with_name("min_match")
                    .help("Number that represents how accurate the observation should be to be included")
                    .default_value("0.0")
                    .takes_value(true)
                    .long("min-match")
            )
            .get_matches();

        Some(Params {
            file_name: matches.value_of("file")?.to_owned(),
            range: matches.value_of("range")?.parse().ok()?,
            depth: matches.value_of("depth")?.parse().ok()?,
            repeat: matches.value_of("repeat")?.parse().ok()?,
            min_match: matches.value_of("min_match")?.parse().ok()?,
        })
    }
}
