//! The data structures
//!
//! The data is given in a form:
//! ```text
//!     lat lon h A z_start h_start z_end h_end t
//! ```
//! where:
//! * `lat` - geographical latitude of observer's location,
//! * `lon` - geographical longitude of observer's location,
//! * `h` - the observer's altitude above the sea,
//! * `A` - observed descent angle,
//! * `z_start` - azimuth of the begining of the observation,
//! * `h_start` - altitude of the begining of the observation,
//! * `z_end` - azimuth of the end of the observation,
//! * `h_end` - altitude of the end of the observation,
//! * `t` - the duration of the observation.
//!
//! The file with data should contain one observation per line.

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::{Error, ErrorKind};

use crate::aprox_eq::AproxEq;
use crate::constants::*;
use crate::structs::*;

/// Represents data given by a witness
#[derive(Debug, Clone)]
pub struct DataSample {
    pub geo_pos: Spherical,
    pub global_pos: Vec3,

    pub descent_angle: f64,
    pub start: Azimuthal,
    pub end: Azimuthal,
    pub duration: f64,

    pub global_start: Vec3,
    pub global_end: Vec3,

    pub trust_da: bool,
    pub trust_start: bool,
    pub trust_end: bool,
}

/// A collenction of observations
#[derive(Clone)]
pub struct Data {
    pub samples: Vec<DataSample>,
    pub mean_pos: Vec3,
}

impl DataSample {
    /// Create new instance of `DataSample` from text
    pub fn from_text(line: &str) -> Option<Self> {
        // parse the line into a list of f64
        let mut nums = line
            .split_whitespace()
            .filter_map(|x| x.parse::<f64>().ok());

        let geo_pos = Spherical {
            lat: nums.next()?.to_radians(),
            lon: nums.next()?.to_radians(),
            r: EARTH_R + nums.next()?,
        };
        let descent_angle = nums.next()?.to_radians();
        let start = Azimuthal {
            z: nums.next()?.to_radians(),
            h: nums.next()?.to_radians(),
        };
        let end = Azimuthal {
            z: nums.next()?.to_radians(),
            h: nums.next()?.to_radians(),
        };
        let duration = nums.next()?;

        // do not trust witness if he claims that start == end
        // and check if line is too long
        if nums.next().is_some() || (start.z.aprox_eq(end.z) && start.h.aprox_eq(end.h)) {
            None
        } else {
            Some(Self {
                geo_pos,
                global_pos: geo_pos.to_vec3(),

                descent_angle,
                start,
                end,
                duration,

                global_start: start.to_vec3().to_global(geo_pos),
                global_end: end.to_vec3().to_global(geo_pos),

                trust_da: descent_angle >= 0.,
                trust_start: start.z >= 0. && start.h >= 0.,
                trust_end: end.z >= 0. && end.h >= 0.,
            })
        }
    }
}

impl Data {
    /// Initialize from file
    pub fn from_file(file_name: &str) -> Result<Self, Error> {
        let mut data = Data {
            samples: Vec::new(),
            mean_pos: Vec3::default(),
        };

        // Open data file
        let mut file = BufReader::new(File::open(file_name)?);
        let mut buffer = String::new();

        // Read data file
        while file.read_line(&mut buffer)? != 0 {
            // Fill the data in
            match DataSample::from_text(&buffer) {
                Some(sample) => {
                    data.mean_pos += sample.global_pos;
                    data.samples.push(sample);
                }
                None => {
                    eprintln!("warning: Skip witness {}", data.samples.len());
                }
            };

            // Clear buffer
            buffer.clear();
        }

        // Two or more witnesses are required
        if data.samples.len() < 2 {
            Err(Error::new(
                ErrorKind::InvalidData,
                "File must contain at least 2 lines",
            ))
        } else {
            data.mean_pos /= data.samples.len() as f64;
            Ok(data)
        }
    }
}
