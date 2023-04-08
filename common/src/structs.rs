//! Basic structures such as Vec3

use std::f64::consts::{PI, TAU};
use std::fmt::{self, Display};

use crate::constants::{DEGREE, EARTH_R};

mod vec3;
pub use vec3::*;

mod line;
pub use line::Line;

pub type UnitQuaternion = nalgebra::UnitQuaternion<f64>;

/// Spherical coordinates triple
#[derive(Debug, Clone, Copy)]
pub struct Spherical {
    pub lat: f64,
    pub lon: f64,
    pub r: f64,
}

impl Display for Spherical {
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
            self.lat.to_degrees(),
            DEGREE,
            lat_dir,
            self.lon.to_degrees().abs(),
            DEGREE,
            lon_dir,
            (self.r - EARTH_R) * 1e-3,
        )
    }
}

impl Spherical {
    pub fn north_direction(self) -> UnitVec3 {
        let z = self.lat.cos();
        let x = -self.lat.sin() * self.lon.cos();
        let y = -self.lat.sin() * self.lon.sin();
        UnitVec3::new_unchecked(Vec3::new(x, y, z))
    }

    pub fn east_direction(self) -> UnitVec3 {
        let x = -self.lon.sin();
        let y = self.lon.cos();
        UnitVec3::new_unchecked(Vec3::new(x, y, 0.0))
    }
}

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

pub type Matrix33 = nalgebra::Matrix3<f64>;

pub trait Matrix33Ext {
    /// Generate x-rotation matrix
    fn rx(angle: f64) -> Self;
    /// Generate z-rotation matrix
    fn rz(angle: f64) -> Self;
}

impl Matrix33Ext for Matrix33 {
    /// ```text
    /// | 1  0       0       |
    /// | 0  cos(a)  -sin(a) |
    /// | 0  sin(a)  cos(a)  |
    /// ```
    fn rx(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        nalgebra::matrix![
            1.0, 0.0,  0.0;
            0.0, cos, -sin;
            0.0, sin,  cos;
        ]
    }

    /// ```text
    /// | cos(a)  -sin(a)  0 |
    /// | sin(a)  cos(a)   0 |
    /// | 0       0        1 |
    /// ```
    fn rz(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        nalgebra::matrix![
            cos, -sin,  0.0;
            sin,  cos,  0.0;
            0.0,  0.0,  1.0;
        ]
    }
}

//
// Convertations
//

impl From<Vec3> for Spherical {
    fn from(val: Vec3) -> Self {
        let xy = f64::hypot(val.x, val.y);
        Spherical {
            lat: f64::atan(val.z / xy),
            lon: f64::atan2(val.y, val.x),
            r: val.norm(),
        }
    }
}

impl From<Spherical> for Vec3 {
    fn from(s: Spherical) -> Self {
        Self::new(
            s.lat.cos() * s.lon.cos(),
            s.lat.cos() * s.lon.sin(),
            s.lat.sin(),
        ) * s.r
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

impl From<Azimuthal> for UnitVec3 {
    fn from(a: Azimuthal) -> Self {
        Self::new_unchecked(Vec3::new(
            a.h.cos() * a.z.sin(),
            a.h.cos() * a.z.cos(),
            a.h.sin(),
        ))
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::approx_eq::ApproxEq;
//     use std::f64::consts::*;
//
//     #[test]
//     fn quaternion_rotate() {
//         let x = Vec3::x();
//         let y = Vec3::y();
//         let z = Vec3::z();
//
//         let q = UnitQuaternion::new(z * FRAC_PI_2);
//         let x_rot = q * x;
//
//         assert!(x_rot.x.approx_eq(y.x));
//         assert!(x_rot.y.approx_eq(y.y));
//         assert!(x_rot.z.approx_eq(y.z));
//     }
// }
