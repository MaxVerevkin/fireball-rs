//! Basic structures such as Vec3

use std::f64::consts::TAU;
use std::fmt::{self, Display};

use crate::constants::DEGREE;

mod vec3;
pub use vec3::*;

mod line;
pub use line::Line;

mod geodetic;
pub use geodetic::*;

/// Azimuthal coordinates tuple
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Azimuthal {
    pub z: f64,
    pub h: f64,
}

impl Display for Azimuthal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({:.3}{DEGREE}, {:.3}{DEGREE})",
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
        Self::new_unchecked(Vec3::new(
            a.h.cos() * a.z.sin(),
            a.h.cos() * a.z.cos(),
            a.h.sin(),
        ))
    }
}
