//! Basic structures such as Vec3

use std::f64::consts::TAU;
use std::fmt::{self, Display};

use crate::constants::DEGREE_SYM;

mod vec3;
pub use vec3::*;

mod line;
pub use line::{Line, LineGrad};

mod geodetic;
pub use geodetic::*;

/// Azimuthal coordinates tuple
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Azimuthal {
    pub z: f64,
    pub h: f64,
}

impl Azimuthal {
    pub fn to_vec3_diff(self, z_diff: f64, h_diff: f64) -> Vec3 {
        let (h_sin, h_cos) = self.h.sin_cos();
        let h_sin_diff = h_cos * h_diff;
        let h_cos_diff = -h_sin * h_diff;

        let (z_sin, z_cos) = self.z.sin_cos();
        let z_sin_diff = z_cos * z_diff;
        let z_cos_diff = -z_sin * z_diff;

        Vec3::new(
            h_cos * z_sin_diff + h_cos_diff * z_sin,
            h_cos * z_cos_diff + h_cos_diff * z_cos,
            h_sin_diff,
        )
    }
}

impl Display for Azimuthal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({:.3}{DEGREE_SYM}, {:.3}{DEGREE_SYM})",
            self.z.to_degrees(),
            self.h.to_degrees()
        )
    }
}

impl From<Vec3> for Azimuthal {
    fn from(val: Vec3) -> Self {
        let xy = f64::hypot(val.x, val.y);
        let z = f64::atan2(val.x, val.y);
        Self {
            z: if z < 0. { z + TAU } else { z },
            h: f64::atan(val.z / xy),
        }
    }
}

impl From<UnitVec3> for Azimuthal {
    fn from(val: UnitVec3) -> Self {
        let xy = f64::hypot(val.x, val.y);
        let z = f64::atan2(val.x, val.y);
        Self {
            z: if z < 0. { z + TAU } else { z },
            h: f64::atan(val.z / xy),
        }
    }
}

impl From<Azimuthal> for UnitVec3 {
    fn from(a: Azimuthal) -> Self {
        // a = (h, z)
        Self::new_unchecked(Vec3::new(
            a.h.cos() * a.z.sin(),
            a.h.cos() * a.z.cos(),
            a.h.sin(),
        ))
    }
}
