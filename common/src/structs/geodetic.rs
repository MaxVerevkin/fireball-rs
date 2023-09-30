use std::f64::consts::{PI, TAU};
use std::fmt;

use super::{UnitVec3, Vec3};
use crate::constants::DEGREE_SYM;

// As per GRS80
// https://en.wikipedia.org//wiki/Geodetic_Reference_System_1980
const A: f64 = 6_378_137.0;
const E_SQ: f64 = 0.006_694_380_022_903_416;
const FRAC_B_A: f64 = 0.996_647_189_318_816_363;

#[derive(Debug, Clone, Copy)]
pub struct Geodetic {
    pub lat: f64,
    pub lon: f64,
    pub h: f64,
}

impl Geodetic {
    pub fn new_from_degrees_m(lat: f64, lon: f64, h: f64) -> Self {
        Self {
            lat: lat.to_radians(),
            lon: lon.to_radians(),
            h,
        }
    }

    pub fn local_cartesian_triple(self) -> (UnitVec3, UnitVec3, UnitVec3) {
        let (sin_lon, cos_lon) = self.lon.sin_cos();
        let (sin_lat, cos_lat) = self.lat.sin_cos();
        let east = UnitVec3::new_unchecked(Vec3::new(-sin_lon, cos_lon, 0.0));
        let north =
            UnitVec3::new_unchecked(Vec3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat));
        let zenith =
            UnitVec3::new_unchecked(Vec3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat));
        (east, north, zenith)
    }

    // https://en.wikipedia.org/wiki/Geodetic_coordinates#Conversion
    pub fn into_geocentric_cartesian(self) -> Vec3 {
        let (sin_lat, cos_lat) = self.lat.sin_cos();
        let prime_vertical_radius = A / f64::sqrt(1.0 - E_SQ * sin_lat * sin_lat);
        let x = (prime_vertical_radius + self.h) * cos_lat * self.lon.cos();
        let y = (prime_vertical_radius + self.h) * cos_lat * self.lon.sin();
        let z = (FRAC_B_A * FRAC_B_A * prime_vertical_radius + self.h) * sin_lat;
        Vec3::new(x, y, z)
    }

    // https://en.wikipedia.org/wiki/Geodetic_coordinates#Conversion
    pub fn from_geocentric_cartesian(vec: Vec3, iters: usize) -> Self {
        let p = f64::hypot(vec.x, vec.y);

        let mut retval = Self {
            lon: f64::atan2(vec.y, vec.x),
            lat: f64::atan(vec.z / p / (1.0 - E_SQ)),
            h: 0.0,
        };

        for _ in 0..iters {
            let (sin_lat, cos_lat) = retval.lat.sin_cos();
            let prime_vertical_radius = A / f64::sqrt(1.0 - E_SQ * sin_lat * sin_lat);
            retval.h = p / cos_lat - prime_vertical_radius;
            retval.lat = f64::atan(
                vec.z
                    / p
                    / (1.0 - E_SQ * prime_vertical_radius / (prime_vertical_radius + retval.h)),
            );
        }

        retval
    }
}

impl fmt::Display for Geodetic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut lon = self.lon;
        if lon > PI {
            lon -= TAU;
        }

        let lat_dir = if self.lat >= 0. { 'N' } else { 'S' };
        let lon_dir = if lon >= 0. { 'E' } else { 'W' };
        write!(
            f,
            "({:.3}{}{}, {:.3}{}{}, {:.3} KM)",
            self.lat.to_degrees().abs(),
            DEGREE_SYM,
            lat_dir,
            self.lon.to_degrees().abs(),
            DEGREE_SYM,
            lon_dir,
            self.h * 1e-3,
        )
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::constants::EARTH_R;

//     use super::super::*;

//     #[test]
//     fn basic() {
//         let lat = 33.6;
//         let lon = 35.0;
//         let h = 40_234.0;

//         let lat2 = 34.6;
//         let lon2 = 33.0;
//         let h2 = 20_234.0;

//         let geodetic = Geodetic::new_from_degrees_m(lat, lon, h);
//         let spherical = Spherical {
//             lat: lat.to_radians(),
//             lon: lon.to_radians(),
//             r: EARTH_R + h,
//         };

//         let geodetic2 = Geodetic::new_from_degrees_m(lat2, lon2, h2);
//         let spherical2 = Spherical {
//             lat: lat2.to_radians(),
//             lon: lon2.to_radians(),
//             r: EARTH_R + h2,
//         };

//         let from_geodetic = geodetic.into_geocentric_cartesian();
//         let from_spherical = Vec3::from(spherical);

//         let from_geodetic2 = geodetic2.into_geocentric_cartesian();
//         let from_spherical2 = Vec3::from(spherical2);

//         let dist_geodetic = (from_geodetic - from_geodetic2).norm();
//         let dist_spherical = (from_spherical - from_spherical2).norm();

//         dbg!(dist_geodetic - dist_spherical);

//         for iters in 0..=5 {
//             eprintln!();
//             let recreated = Geodetic::from_geocentric_cartesian(from_geodetic, iters);
//             dbg!(iters);
//             dbg!((geodetic.lat - recreated.lat).to_degrees());
//             dbg!(geodetic.h - recreated.h);
//         }

//         eprintln!();
//         let (e, n, z) = geodetic.local_cartesian_triple();
//         dbg!(e.dot(*n));
//         dbg!(e.dot(*z));
//         dbg!(n.dot(*z));

//         panic!()
//     }
// }
