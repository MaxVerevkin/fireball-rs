//! The data structures
//!
//! The data is given in a form:
//! ```
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
//! The file with data should contain only the numbers, one observation per line.

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::{Error, ErrorKind};

use crate::constants::*;
use crate::structs::*;

/// Represents data given by a witness
#[derive(Debug)]
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
pub struct Data {
    pub samples: Vec<DataSample>,
    pub mean_lon: f64,
    pub mean_lat: f64,
}

impl DataSample {
    /// Construct new sample from an array of `f64` given in the same order as in file
    pub fn from_array(arr: &[f64; 9]) -> Self {
        let geo_pos = Spherical {
            lat: arr[0].to_radians(),
            lon: arr[1].to_radians(),
            r: EARTH_R + arr[2],
        };
        let descent_angle = arr[3].to_radians();
        let mut start = Azimuthal {
            z: arr[4].to_radians(),
            h: arr[5].to_radians(),
        };
        let mut end = Azimuthal {
            z: arr[6].to_radians(),
            h: arr[7].to_radians(),
        };

        // do not trust witness if he claims that start == end
        if (start.z - end.z).abs() < 1e-5 && (start.h - end.h).abs() < 1e-4 {
            end.z = -1.;
            start.z = -1.;
        }

        Self {
            geo_pos,
            global_pos: geo_pos.to_vec3(),

            descent_angle,
            start,
            end,
            duration: arr[8],

            global_start: start.to_vec3().to_global(&geo_pos),
            global_end: end.to_vec3().to_global(&geo_pos),

            trust_da: descent_angle >= 0.,
            trust_start: start.z >= 0. && start.h >= 0.,
            trust_end: end.z >= 0. && end.h >= 0.,
        }
    }
}

impl Data {
    /// Initialize from file
    pub fn from_file(file_name: &str) -> Result<Self, Error> {
        let mut data = Data {
            samples: Vec::new(),
            mean_lon: 0.,
            mean_lat: 0.,
        };

        // Use to compute mean longitude
        let mut mean_lon_sin = 0.;
        let mut mean_lon_cos = 0.;

        // Open data file
        let file = File::open(file_name)?;
        let mut file = BufReader::new(file);
        let mut buffer = String::new();

        // Read data file
        while file.read_line(&mut buffer)? != 0 {
            // Check line
            let mut nums = [0f64; 9];
            let mut words = buffer.split_whitespace();
            let bad_line_error = Err(Error::new(
                ErrorKind::InvalidData,
                "Each line in the file must contain 9 numbers",
            ));
            for i in 0..9 {
                let word = match words.next() {
                    Some(w) => w,
                    None => return bad_line_error,
                };
                nums[i] = match word.parse() {
                    Ok(n) => n,
                    Err(_) => return bad_line_error,
                };
            }
            if words.next().is_some() {
                return bad_line_error;
            }

            // Fill the data in
            let sample = DataSample::from_array(&nums);
            data.mean_lat += sample.geo_pos.lat;
            mean_lon_sin += sample.geo_pos.lon.sin();
            mean_lon_cos += sample.geo_pos.lon.cos();
            data.samples.push(sample);

            // Clear buffer
            buffer.clear();
        }

        // Mean lon and mean lat
        data.mean_lat /= data.samples.len() as f64;
        data.mean_lon = f64::atan2(mean_lon_sin, mean_lon_cos);

        // Two or more witnesses are required
        if data.samples.len() < 2 {
            Err(Error::new(
                ErrorKind::InvalidData,
                "File must contain at least 2 lines",
            ))
        } else {
            Ok(data)
        }
    }
}
