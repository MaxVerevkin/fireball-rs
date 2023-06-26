use std::cmp::Ordering;
use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::ops::{Deref, DerefMut};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::constants::*;
use crate::maths::*;
use crate::quick_median::SliceExt;
use crate::structs::*;

/// The sample as it is in the file. Angles are meausred in degrees.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct RawSample {
    pub lat: f64,
    pub lon: f64,
    pub h: f64,
    pub a: Option<f64>,
    pub zb: Option<f64>,
    pub hb: Option<f64>,
    pub z0: Option<f64>,
    pub h0: Option<f64>,
    pub t: Option<f64>,
    pub name: Option<String>,
    pub exp: Option<f64>,
}

/// `RawSample` with degrees converted to radians
#[derive(Deserialize, Serialize, Debug, Clone)]
#[serde(from = "RawSample")]
#[serde(into = "RawSample")]
pub struct NormalizedRawSample(pub RawSample);

impl From<RawSample> for NormalizedRawSample {
    fn from(value: RawSample) -> Self {
        Self(RawSample {
            lat: value.lat.to_radians(),
            lon: value.lon.to_radians(),
            a: value.a.map(f64::to_radians),
            zb: value.zb.map(f64::to_radians),
            hb: value.hb.map(f64::to_radians),
            z0: value.z0.map(f64::to_radians),
            h0: value.h0.map(f64::to_radians),
            ..value
        })
    }
}

impl From<NormalizedRawSample> for RawSample {
    fn from(value: NormalizedRawSample) -> Self {
        Self {
            lat: value.0.lat.to_degrees(),
            lon: value.0.lon.to_degrees(),
            a: value.0.a.map(f64::to_degrees),
            zb: value.0.zb.map(f64::to_degrees),
            hb: value.0.hb.map(f64::to_degrees),
            z0: value.0.z0.map(f64::to_degrees),
            h0: value.0.h0.map(f64::to_degrees),
            ..value.0
        }
    }
}

/// Pre-processed sample. Angles are in radians.
#[derive(Deserialize, Debug, Clone)]
#[serde(from = "NormalizedRawSample")]
pub struct DataSample {
    pub geo_location: Geodetic,
    pub location: Vec3,
    pub east_dir: UnitVec3,
    pub north_dir: UnitVec3,
    pub zenith_dir: UnitVec3,

    pub da: Option<f64>,
    pub k_start: Option<UnitVec3>,
    pub k_end: Option<UnitVec3>,
    pub axis: Option<UnitVec3>,
    pub plane: Option<UnitVec3>,
    pub dur: Option<f64>,
    pub name: Option<String>,
    pub exp: Option<f64>,

    pub z0: Option<f64>,
    pub h0: Option<f64>,
}

impl From<NormalizedRawSample> for DataSample {
    fn from(raw: NormalizedRawSample) -> Self {
        let geo_location = Geodetic {
            lat: raw.0.lat,
            lon: raw.0.lon,
            h: raw.0.h,
        };
        let location = geo_location.into_geocentric_cartesian();
        let (east_dir, north_dir, zenith_dir) = geo_location.local_cartesian_triple();

        let make_k_vec = move |z: Option<f64>, h: Option<f64>| {
            let local = UnitVec3::from(Azimuthal { z: z?, h: h? });
            let v = local.x * east_dir.into_inner()
                + local.y * north_dir.into_inner()
                + local.z * zenith_dir.into_inner();
            Some(UnitVec3::new_unchecked(v))
        };

        let mut k_start = make_k_vec(raw.0.zb, raw.0.hb);
        let k_end = make_k_vec(raw.0.z0, raw.0.h0);

        if let Some(ks) = k_start && let Some(ke) = k_end && ks.dot(*ke).abs() > 0.9999 {
            k_start = None;
        }

        let axis = match (k_start, k_end) {
            (Some(a), Some(b)) => Some(UnitVec3::new_normalize(a.into_inner() + b.into_inner())),
            (Some(k), None) | (None, Some(k)) => Some(k),
            (None, None) => None,
        };

        let plane = if let (Some(v1), Some(v2)) = (k_start, k_end) {
            Some(UnitVec3::new_normalize(v2.cross(*v1)))
        } else if let (Some(axis), Some(da)) = (axis, raw.0.a) {
            // TODO mention
            if axis.normalize().dot(location.normalize()) > 0.999 {
                None
            } else {
                Some(UnitVec3::new_normalize(axis.cross(location)).rotate(axis, da + PI))
            }
        } else {
            None
        };

        Self {
            geo_location,
            location,
            east_dir,
            north_dir,
            zenith_dir,

            da: raw.0.a,
            k_start,
            k_end,
            axis,
            plane,
            dur: raw.0.t,
            name: raw.0.name,
            exp: raw.0.exp,

            z0: raw.0.z0,
            h0: raw.0.h0,
        }
    }
}

impl DataSample {
    /// Compute the descent angle given `k_start` and `k_end`
    pub fn calculated_da(&self) -> Option<f64> {
        let axis = self.axis?;
        let plane = self.plane?;

        // if (axis^location) < 5 degrees
        if self.location.normalize().dot(*axis) > 0.996 {
            return None;
        }

        let i = axis.cross(self.location).normalize();
        let j = i.cross(*axis);

        let x = plane.dot(i);
        let y = plane.dot(j);

        let da = f64::atan2(x, y) + FRAC_PI_2;
        Some(if da < 0. { da + TAU } else { da })
    }

    /// Check whether a direction matches the observation
    pub fn direction_matches(&self, direction: UnitVec3) -> bool {
        if let Some(axis) = self.axis && let Some(plane) = self.plane {
            axis.cross(*plane).dot(*direction) > 0.0
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
        let perp = direction.cross(point - self.location).cross(*direction);

        if let Some((k_start, k_end)) = self.k_start.zip(self.k_end) {
            // When we have both k_start and k_end, we must ensure that they both lie in the
            // correct hemispace
            if k_start.dot(perp) <= 0.0 || k_end.dot(perp) <= 0.0 {
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
            if axis.dot(perp) <= 0.0 {
                return false;
            }
            if let Some(observation_norm) = self.plane {
                let dir = axis.cross(*observation_norm);
                if dir.dot(*direction) <= 0.0 {
                    return false;
                }
            }
        }

        true
    }

    /// Calculate azimuth in range [0, Tau)
    pub fn calc_azimuth(&self, p: Vec3) -> f64 {
        let k = p - self.location;
        let x = k.dot(*self.east_dir);
        let y = k.dot(*self.north_dir);
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
    h: Option<f64>,
    vx: Option<f64>,
    vy: Option<f64>,
    vz: Option<f64>,
}

#[derive(Deserialize, Debug, Clone, Copy, PartialEq)]
#[serde(from = "RawAnswer")]
pub struct Answer {
    pub traj: Option<Line>,
    pub speed: Option<f64>,
    pub point: Option<Vec3>,
    pub point_no_h: Vec3,
}

impl From<RawAnswer> for Answer {
    fn from(raw: RawAnswer) -> Self {
        let geo_location = Geodetic::new_from_degrees_m(raw.lat, raw.lon, raw.h.unwrap_or(0.0));
        let point_no_h = geo_location.into_geocentric_cartesian();
        let point = raw.h.is_some().then_some(point_no_h);

        let (direction, speed) = if let (Some(vx), Some(vy), Some(vz)) = (raw.vx, raw.vy, raw.vz) {
            Some(UnitVec3::new_and_get(Vec3::new(vx, vy, vz) * 1_000.0)).unzip()
        } else {
            (None, None)
        };

        Self {
            traj: point
                .zip(direction)
                .map(|(point, direction)| Line { point, direction }),
            speed,
            point,
            point_no_h,
        }
    }
}

/// A collenction of observations
#[derive(Debug, Clone)]
pub struct Data {
    pub samples: Vec<DataSample>,
    pub name: Option<String>,
    pub answer: Option<Answer>,
    pub meta: toml::Table,
}

#[derive(Deserialize, Debug, Clone)]
pub struct RawData {
    pub sample: Vec<NormalizedRawSample>,
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub answer: Option<Answer>,
    #[serde(default)]
    pub meta: toml::Table,
}

impl RawData {
    pub fn read_from_toml_file(path: impl AsRef<Path>) -> Self {
        let contents = std::fs::read(path.as_ref()).unwrap();
        let contents = String::from_utf8_lossy(&contents);

        let mut retval: Self = toml::from_str(&contents).unwrap();

        if retval.name.is_none() {
            retval.name = path
                .as_ref()
                .file_name()
                .map(|file_name| file_name.to_string_lossy())
                .and_then(|file_name| Some(file_name.split_once('.')?.0.to_owned()));
        }

        retval
    }

    pub fn da_correction(mut self, da_k: f64) -> Self {
        if da_k != 0.0 {
            for s in &mut self.sample {
                if let Some(da) = &mut s.0.a {
                    if (FRAC_PI_2..(FRAC_PI_2 * 3.0)).contains(da) {
                        *da = *da - da_k * f64::sin(*da * 2.0)
                    }
                }
            }
        }
        self
    }

    pub fn az_correction(mut self, az_k: f64) -> Self {
        if az_k != 0.0 {
            for s in &mut self.sample {
                if let Some(az) = &mut s.0.z0 {
                    *az = *az + az_k * f64::sin(*az)
                }
            }
        }
        self
    }

    pub fn finalize(self) -> Data {
        Data {
            samples: self.sample.into_iter().map(Into::into).collect(),
            name: self.name,
            answer: self.answer,
            meta: self.meta,
        }
    }
}

impl Data {
    pub fn compare(&self, other: Line, message: &str) {
        let Some(answer) = self.answer else { return };
        let Some(answer_point) = answer.point else { return };

        let vel_error = answer
            .traj
            .map(|traj| traj.direction.dot(*other.direction).acos().to_degrees());

        let distance = (answer_point - other.point).cross(*other.direction).norm() * 1e-3;

        let mut med_buf = Vec::new();
        med_buf.extend(self.samples.iter().map(|s| s.geo_location.h));
        let med_h = med_buf.median();
        med_buf.clear();
        med_buf.extend(self.samples.iter().map(|s| s.geo_location.lat));
        let med_lat = med_buf.median();
        med_buf.clear();
        med_buf.extend(self.samples.iter().map(|s| s.geo_location.lon));
        let med_lon = med_buf.median();
        let med_observer = Geodetic {
            lat: med_lat,
            lon: med_lon,
            h: med_h,
        }
        .into_geocentric_cartesian();

        let n = answer.traj.map(|traj| {
            traj.direction
                .into_inner()
                .cross((med_observer - traj.point).normalize())
        });
        let elevation_angle = n.map(|n| n.dot(*other.direction).asin().to_degrees());

        let geo = Geodetic::from_geocentric_cartesian(other.point, 10);

        eprintln!("--- {message} ---");
        eprintln!("Point: {geo}");
        eprintln!("Distance: {distance:.1}km");
        if let Some(vel_error) = vel_error {
            eprintln!("Velocity angular error: {vel_error:.1}{DEGREE}");
        }
        if let Some(elevation_angle) = elevation_angle {
            eprintln!("Elevation angular error: {elevation_angle:.1}{DEGREE}");
        }
        eprintln!();
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
                        + s.north_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.zenith_dir * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(0.0, z_north).abs() < 0.001);
            }

            {
                let z_east = s.calc_azimuth(
                    s.location
                        + s.east_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.zenith_dir * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(FRAC_PI_2, z_east).abs() < 0.001);
            }

            {
                let z_south = s.calc_azimuth(
                    s.location - s.north_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.zenith_dir * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(PI, z_south).abs() < 0.001);
            }

            {
                let z_west = s.calc_azimuth(
                    s.location - s.east_dir * (10.0 + random::<f64>() * 1_000.0)
                        + s.zenith_dir * random::<f64>() * 1_000.0,
                );
                assert!(angle_diff(FRAC_PI_2 * 3.0, z_west).abs() < 0.001);
            }
        }
    }
}
