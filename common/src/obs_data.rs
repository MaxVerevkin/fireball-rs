use std::cmp::Ordering;
use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::fs::File;
use std::io::Read;
use std::ops::{Deref, DerefMut};
use std::path::Path;

use serde::{Deserialize, Serialize};

use nalgebra::Unit;

use crate::constants::*;
use crate::maths::*;
use crate::quick_median::SliceExt;
use crate::structs::*;

#[derive(Deserialize, Serialize)]
struct RawSample {
    lat: f64,
    lon: f64,
    h: f64,
    a: Option<f64>,
    zb: Option<f64>,
    hb: Option<f64>,
    z0: Option<f64>,
    h0: Option<f64>,
    t: Option<f64>,
    name: Option<String>,
    exp: Option<f64>,
}

#[derive(Deserialize, Debug, Clone)]
#[serde(from = "RawSample")]
pub struct DataSample {
    pub geo_location: Spherical,
    pub location: Vec3,
    pub east_dir: UnitVec3,
    pub north_dir: UnitVec3,

    pub da: Option<f64>,
    pub k_start: Option<UnitVec3>,
    pub k_end: Option<UnitVec3>,
    pub axis: Option<UnitVec3>,
    pub plane: Option<UnitVec3>,
    pub dur: Option<f64>,
    pub name: Option<String>,
    pub exp: Option<f64>,

    pub z0: Option<f64>,
}

impl From<RawSample> for DataSample {
    fn from(raw: RawSample) -> Self {
        let geo_location = Spherical {
            lat: raw.lat.to_radians(),
            lon: raw.lon.to_radians(),
            r: EARTH_R + raw.h,
        };

        let location: Vec3 = geo_location.into();

        let da = raw.a.map(f64::to_radians);

        let make_k_vec = move |z: Option<f64>, h: Option<f64>| {
            Some(
                UnitVec3::from(Azimuthal {
                    z: z?.to_radians(),
                    h: h?.to_radians(),
                })
                .to_global(geo_location),
            )
        };

        let mut k_start = make_k_vec(raw.zb, raw.hb);
        let k_end = make_k_vec(raw.z0, raw.h0);

        if let Some(ks) = k_start && let Some(ke) = k_end && ks.dot(&ke).abs() > 0.9999 {
            k_start = None;
        }

        let axis = match (k_start, k_end) {
            (Some(a), Some(b)) => Some(Unit::new_normalize(a.into_inner() + b.into_inner())),
            (Some(k), None) | (None, Some(k)) => Some(k),
            (None, None) => None,
        };

        let plane = if let (Some(v1), Some(v2)) = (k_start, k_end) {
            Some(Unit::new_normalize(v2.cross(&v1)))
        } else if let (Some(axis), Some(da)) = (axis, da) {
            // TODO mention
            if axis.normalize().dot(&location.normalize()) > 0.999 {
                None
            } else {
                let plane = Unit::new_normalize(axis.cross(&location));
                let q = UnitQuaternion::from_axis_angle(&axis, da + PI);
                Some(q * plane)
            }
        } else {
            None
        };

        Self {
            geo_location,
            location,
            east_dir: geo_location.east_direction(),
            north_dir: geo_location.north_direction(),
            // east_dir: Vec3::x_axis().to_local(geo_location),
            // north_dir: Vec3::y_axis().to_local(geo_location),
            da,
            k_start,
            k_end,
            axis,
            plane,
            dur: raw.t,
            name: raw.name,
            exp: raw.exp,

            z0: raw.z0.map(|z| z.to_radians()),
        }
    }
}

// impl From<DataSample> for RawSample {
//     fn from(value: DataSample) -> Self {
//         todo!()
//     }
// }

impl DataSample {
    /// Compute the descent angle given `k_start` and `k_end`
    pub fn calculated_da(&self) -> Option<f64> {
        let axis = self.axis?;
        let plane = self.plane?;

        // if (axis^location) < 5 degrees
        if self.location.normalize().dot(&axis) > 0.996 {
            return None;
        }

        let i = axis.cross(&self.location).normalize();
        let j = i.cross(&axis);

        let x = plane.dot(&i);
        let y = plane.dot(&j);

        let da = f64::atan2(x, y) + FRAC_PI_2;
        Some(if da < 0. { da + TAU } else { da })
    }

    /// Check whether a direction matches the observation
    pub fn direction_matches(&self, direction: UnitVec3) -> bool {
        if let Some(axis) = self.axis && let Some(plane) = self.plane {
            axis.cross(&plane).dot(&direction) > 0.0
        } else {
            false
        }
    }

    /// Returns `false` if this observation conflicts with the proposed trajectory
    pub fn observation_matches(&self, traj: Line) -> bool {
        let Line { point, direction } = traj;

        // W/o the axis we really cannot do anything.
        let axis = match self.axis {
            Some(axis) => axis,
            None => return false,
        };
        //
        // A vertor that describes the hemispace of "allowed" observations
        let perp = direction.cross(&(point - self.location)).cross(&direction);

        if let Some((k_start, k_end)) = self.k_start.zip(self.k_end) {
            // When we have both k_start and k_end, we must ensure that they both lie in the
            // correct hemispace
            if k_start.dot(&perp) <= 0.0 || k_end.dot(&perp) <= 0.0 {
                return false;
            }
            let l_start = lambda(self.location, k_start, point, direction);
            let l_end = lambda(self.location, k_end, point, direction);
            match l_end.partial_cmp(&l_start) {
                Some(Ordering::Less | Ordering::Equal) | None => return false,
                _ => (),
            }
        } else {
            // The best we have is the axis, so let's check it
            if axis.dot(&perp) <= 0.0 {
                return false;
            }
            if let Some(observation_norm) = self.plane {
                let dir = axis.cross(&observation_norm);
                if dir.dot(&direction) <= 0.0 {
                    return false;
                }
            }
        }

        true
    }

    /// Calculate azimuth in range [0, Tau)
    pub fn calc_azimuth(&self, p: Vec3) -> f64 {
        let k = p - self.location;
        let x = k.dot(&self.east_dir);
        let y = k.dot(&self.north_dir);
        let atan = f64::atan2(x, y);
        if atan < 0.0 {
            atan + TAU
        } else {
            atan
        }
    }
}

#[derive(Deserialize)]
pub struct RawAnswer {
    lat: f64,
    lon: f64,
    h: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

#[derive(Deserialize, Debug, Clone, Copy, PartialEq)]
#[serde(from = "RawAnswer")]
pub struct Answer {
    pub traj: Line,
    pub speed: f64,
}

impl From<RawAnswer> for Answer {
    fn from(raw: RawAnswer) -> Self {
        let geo_location = Spherical {
            lat: raw.lat.to_radians(),
            lon: raw.lon.to_radians(),
            r: EARTH_R + raw.h,
        };
        let direction = Vec3::new(raw.vx, raw.vy, raw.vz) * 1_000.0;
        let (direction, speed) = Unit::new_and_get(direction);
        Self {
            traj: Line {
                point: geo_location.into(),
                direction,
            },
            speed,
        }
    }
}

/// A collenction of observations
#[derive(Deserialize, Debug, Clone, Default)]
pub struct Data {
    #[serde(rename = "sample")]
    pub samples: Vec<DataSample>,
    pub answer: Option<Answer>,
}

impl Data {
    pub fn apply_da_correction(&mut self, a: f64) {
        for s in &mut self.samples {
            if let Some(da) = &mut s.da {
                if *da >= FRAC_PI_2 && *da <= (FRAC_PI_2 * 3.0) {
                    *da = *da - a * f64::sin(*da * 2.0);
                }
            }
        }
    }

    // TODO: return errors instead of panicing!
    pub fn read_from_toml(path: impl AsRef<Path>) -> Self {
        let mut in_file = File::open(path).unwrap();
        let mut buf = Vec::new();
        in_file.read_to_end(&mut buf).unwrap();
        toml::from_str(&String::from_utf8_lossy(&buf)).unwrap()
    }

    pub fn compare(&self, other: Line, message: &str) {
        if let Some(answer) = self.answer {
            let vel_error = answer
                .traj
                .direction
                .dot(&other.direction)
                .acos()
                .to_degrees();

            let distance = (answer.traj.point - other.point)
                .cross(&other.direction)
                .norm()
                * 1e-3;

            let mut med_buf = Vec::new();
            med_buf.extend(self.samples.iter().map(|s| s.geo_location.r));
            let med_r = med_buf.median();
            med_buf.clear();
            med_buf.extend(self.samples.iter().map(|s| s.geo_location.lat));
            let med_lat = med_buf.median();
            med_buf.clear();
            med_buf.extend(self.samples.iter().map(|s| s.geo_location.lon));
            let med_lon = med_buf.median();
            let med_observer: Vec3 = Spherical {
                lat: med_lat,
                lon: med_lon,
                r: med_r,
            }
            .into();
            let n = answer
                .traj
                .direction
                .into_inner()
                .cross(&(med_observer - answer.traj.point).normalize());
            let elevation_angle = n.dot(&other.direction.into_inner()).asin().to_degrees();

            let spherical: Spherical = other.point.into();

            eprintln!("--- {message} ---\nVelocity angular error: {vel_error:.1}{DEGREE}\nElevation angle error: {elevation_angle:.1}{DEGREE}\nDistance: {distance:.1}km\nPoint: {spherical}\n");
        }
    }
}

impl Deref for Data {
    type Target = Vec<DataSample>;

    fn deref(&self) -> &Self::Target {
        &self.samples
    }
}

impl DerefMut for Data {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.samples
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn azimuth() {
        for _ in 0..10_000 {
            let raw = RawSample {
                lat: random::<f64>() * 180.0 - 90.0,
                lon: random::<f64>() * 360.0 - 180.0,
                h: random::<f64>() * 2_000.0,
                a: None,
                zb: None,
                hb: None,
                z0: None,
                h0: None,
                t: None,
                name: None,
                exp: None,
            };
            let s = DataSample::from(raw);

            {
                let z_north = s.calc_azimuth(
                    s.location
                        + &*s.north_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.location.normalize() * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(0.0, z_north).abs() < 0.001);
            }

            {
                let z_east = s.calc_azimuth(
                    s.location
                        + &*s.east_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.location.normalize() * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(FRAC_PI_2, z_east).abs() < 0.001);
            }

            {
                let z_south = s.calc_azimuth(
                    s.location - &*s.north_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.location.normalize() * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(PI, z_south).abs() < 0.001);
            }

            {
                let z_west = s.calc_azimuth(
                    s.location - &*s.east_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.location.normalize() * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(FRAC_PI_2 * 3.0, z_west).abs() < 0.001);
            }
        }
    }
}
