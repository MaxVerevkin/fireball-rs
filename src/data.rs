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
use std::io::{self, BufReader};

use crate::constants::*;
use crate::structs::*;

/// Represents data given by a witness
#[derive(Debug, Clone, Copy)]
pub struct DataSample {
    pub pos: Vec3,
    pub geo_pos: Spherical,

    pub descent_angle: f64,
    pub start: Azimuthal,
    pub end: Azimuthal,
    pub duration: f64,

    pub global_start: Vec3,
    pub global_end: Vec3,
}

impl DataSample {
    /// Create new instance of `DataSample` from text
    pub fn from_text(line: &str) -> Option<Self> {
        // Get exactly 9 numbers
        let mut nums = [0f64; 9];
        let mut words = line.split_ascii_whitespace();
        for n in &mut nums {
            *n = words.next()?.parse::<f64>().ok()?;
        }
        if words.next().is_some() {
            return None;
        }

        let geo_pos = Spherical {
            lat: nums[0].to_radians(),
            lon: nums[1].to_radians(),
            r: EARTH_R + nums[2],
        };
        let descent_angle = nums[3].to_radians();
        let start = Azimuthal {
            z: nums[4].to_radians(),
            h: nums[5].to_radians(),
        };
        let end = Azimuthal {
            z: nums[6].to_radians(),
            h: nums[7].to_radians(),
        };
        let duration = nums[8];

        Some(Self {
            geo_pos,
            pos: geo_pos.into(),

            descent_angle,
            start,
            end,
            duration,

            global_start: Vec3::from(start).to_global(geo_pos),
            global_end: Vec3::from(end).to_global(geo_pos),
        })
    }

    // pub fn to_text(&self) -> String {
    //     format!(
    //         "{} {} {} {} {} {} {} {} {}",
    //         self.geo_pos.lat.to_degrees(),
    //         self.geo_pos.lon.to_degrees(),
    //         self.geo_pos.r - EARTH_R,
    //         self.descent_angle.to_degrees(),
    //         self.start.z.to_degrees(),
    //         self.start.h.to_degrees(),
    //         self.end.z.to_degrees(),
    //         self.end.h.to_degrees(),
    //         self.duration
    //     )
    // }
}

/// A collenction of observations
#[derive(Debug, Clone)]
pub struct Data {
    pub samples: Vec<DataSample>,
}

impl Data {
    /// Initialize from file
    pub fn from_file(file_name: &str) -> Result<Self, io::Error> {
        let mut file = BufReader::new(File::open(file_name)?);
        let mut samples = Vec::with_capacity(400);
        let mut buf = String::new();
        let mut skipped = 0usize;

        while {
            buf.clear();
            file.read_line(&mut buf)?
        } != 0
        {
            match DataSample::from_text(&buf) {
                Some(s) => samples.push(s),
                None => skipped += 1,
            }
        }

        if skipped > 0 {
            eprintln!("warning: {} lines skipped", skipped);
        }

        Ok(Self { samples })
    }

    // pub fn to_file(&self, file: &str) -> std::io::Result<()> {
    //     let mut file = File::create(file)?;
    //     for sample in &self.samples {
    //         writeln!(file, "{}", sample.to_text())?;
    //     }
    //     Ok(())
    // }
}
