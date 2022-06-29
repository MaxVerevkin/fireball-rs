use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::ops::{Deref, DerefMut};

use crate::constants::*;
use crate::maths::*;
use crate::structs::*;

use serde::Deserialize;

#[derive(Deserialize)]
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
    pub da: Option<f64>,
    pub k_start: Option<Vec3>,
    pub k_end: Option<Vec3>,
    pub axis: Option<Vec3>,
    pub plane: Option<Vec3>,
    pub dur: Option<f64>,
    pub name: Option<String>,
    pub exp: Option<f64>,
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

        let mut k_start = raw.zb.zip(raw.hb).map(|(z, h)| {
            Vec3::from(Azimuthal {
                z: z.to_radians(),
                h: h.to_radians(),
            })
            .to_global(geo_location)
        });
        let k_end = raw.z0.zip(raw.h0).map(|(z, h)| {
            Vec3::from(Azimuthal {
                z: z.to_radians(),
                h: h.to_radians(),
            })
            .to_global(geo_location)
        });

        if let Some((ks, ke)) = k_start.zip(k_end) {
            if ks.dot(ke).abs() > 0.9999 {
                k_start = None;
            }
        }

        let axis = match (k_start, k_end) {
            (Some(a), Some(b)) => Some((a + b).normalized()),
            (Some(k), None) | (None, Some(k)) => Some(k),
            (None, None) => None,
        };

        let plane = if let (Some(v1), Some(v2)) = (k_start, k_end) {
            Some(v2.cross(v1).normalized())
        } else if let (Some(axis), Some(da)) = (axis, da) {
            // TODO mention
            if axis.normalized().dot(location.normalized()) > 0.999 {
                None
            } else {
                let plane = axis.cross(location).normalized();
                let q = UnitQuaternion::new(axis, da + PI);
                Some(q * plane)
            }
        } else {
            None
        };

        Self {
            geo_location,
            location,
            da,
            k_start,
            k_end,
            axis,
            plane,
            dur: raw.t,
            name: raw.name,
            exp: raw.exp,
        }
    }
}

impl DataSample {
    /// Compute the descent angle given `k_start` and `k_end`
    pub fn calculated_da(&self) -> Option<f64> {
        let (axis, plane) = self.axis.zip(self.plane)?;

        // if (axis^location) < 5 degrees
        if self.location.normalized().dot(axis) > 0.996 {
            return None;
        }

        let i = axis.cross(self.location).normalized();
        let j = i.cross(axis);

        let x = plane.dot(i);
        let y = plane.dot(j);

        let da = f64::atan2(x, y) + FRAC_PI_2;
        Some(if da < 0. { da + TAU } else { da })
    }

    /// Check whether a direction matches the observation
    pub fn direction_matches(&self, direction: Vec3) -> bool {
        let dir = self.axis.and_then(|axis| Some(axis.cross(self.plane?)));
        dir.map_or(false, |dir| dir.dot(direction) > 0.0)
    }

    /// Returns `false` if this observation conflicts with the proposed trajectory
    pub fn observation_matches(&self, point: Vec3, direction: Vec3) -> bool {
        // W/o the axis we really cannot do anything.
        let axis = match self.axis {
            Some(axis) => axis,
            None => return false,
        };

        // A normal to a plane passing through the observer and proposed trajectory
        let norm = direction.cross(point - self.location).normalized();

        // A vertor that describes the hemispace of "allowed" observations
        let perp = norm.cross(direction).normalized();

        if let Some((k_start, k_end)) = self.k_start.zip(self.k_end) {
            // When we have both k_start and k_end, we must ensure that they both lie in the
            // correct hemispace
            if k_start.dot(perp) <= 0.0 || k_end.dot(perp) <= 0.0 {
                return false;
            }
            let l_start = lambda(self.location, k_start, point, direction);
            let l_end = lambda(self.location, k_end, point, direction);
            // Using !(..>..) instead of (..<=..) to catch NaNs
            if !(l_end > l_start) {
                return false;
            }
        } else {
            // The best we have is the axis, so let's check it
            if axis.dot(perp) <= 0.0 {
                return false;
            }
            if let Some(observation_norm) = self.plane {
                let dir = axis.cross(observation_norm);
                if dir.dot(direction) <= 0.0 {
                    return false;
                }
            }
        }

        true
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
pub struct Answer(pub Line);

impl From<RawAnswer> for Answer {
    fn from(raw: RawAnswer) -> Self {
        let geo_location = Spherical {
            lat: raw.lat.to_radians(),
            lon: raw.lon.to_radians(),
            r: EARTH_R + raw.h,
        };
        Self(Line {
            point: geo_location.into(),
            direction: Vec3::new(raw.vx, raw.vy, raw.vz) * 1e6,
        })
    }
}

/// A collenction of observations
#[derive(Deserialize, Debug, Clone)]
pub struct Data {
    #[serde(rename = "sample")]
    pub samples: Vec<DataSample>,
    pub answer: Option<Answer>,
}

impl Data {
    pub fn compare(&self, other: Line, message: &str) {
        if let Some(Answer(Line { point, direction })) = self.answer {
            eprintln!(
                "--- {message} ---\nVelocity angular error: {:.1}{DEGREE}\nDistance: {:.1}km\n",
                direction
                    .normalized()
                    .dot(other.direction.normalized())
                    .acos()
                    .to_degrees(),
                (point - other.point)
                    .cross(other.direction.normalized())
                    .len()
                    / 1e3
            );
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
