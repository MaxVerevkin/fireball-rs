//! Basic structures such as Vec3

use std::f64::consts::{PI, TAU};
use std::fmt::{self, Display};
use std::ops;

use crate::constants::{DEGREE, EARTH_R};

mod vec3;
pub use vec3::Vec3;

mod line;
pub use line::Line;

/// Rotate a vector about any axis
#[derive(Debug, Clone, Copy)]
pub struct UnitQuaternion {
    s: f64,
    v: Vec3,
}

impl UnitQuaternion {
    /// Create a new unit quaternion
    pub fn new(axis: Vec3, angle: f64) -> Self {
        Self {
            s: f64::cos(angle / 2.),
            v: axis * f64::sin(angle / 2.),
        }
    }

    /// Compute `q^-1`
    pub fn inverse(&self) -> Self {
        Self {
            s: self.s,
            v: -self.v,
        }
    }
}

impl ops::Mul for UnitQuaternion {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let a = self.s;
        let b = self.v.x;
        let c = self.v.y;
        let d = self.v.z;

        let ap = rhs.s;
        let bp = rhs.v.x;
        let cp = rhs.v.y;
        let dp = rhs.v.z;

        Self {
            s: a * ap - b * bp - c * cp - d * dp, // 1
            v: Vec3 {
                x: a * bp - d * cp + b * ap + c * dp, // i
                y: a * cp - b * dp + c * ap + d * bp, // j
                z: a * dp - c * bp + b * cp + d * ap, // k
            },
        }

        // Self {
        //     s: self.s * rhs.s - self.v.dot(rhs.v),
        //     v: self.v * rhs.s + rhs.v * self.s + self.v.cross(rhs.v),
        // }
    }
}

impl ops::Mul<Vec3> for UnitQuaternion {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        let p = Self { s: 0., v: rhs };
        let q = self * p * self.inverse();
        q.v
    }
}

impl ops::Mul<&Vec3> for UnitQuaternion {
    type Output = Vec3;

    fn mul(self, rhs: &Vec3) -> Vec3 {
        let p = Self { s: 0., v: *rhs };
        let q = self * p * self.inverse();
        q.v
    }
}

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

/// Azimuthal coordinates tuple
#[derive(Debug, Clone, Copy)]
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

/// 3x3 Matirix (column-major)
#[derive(Debug, Clone, Copy)]
pub struct Matrix33([f64; 9]);

impl Matrix33 {
    /// Generate x-rotation matrix
    ///
    /// ```text
    /// | 1  0       0       |
    /// | 0  cos(a)  -sin(a) |
    /// | 0  sin(a)  cos(a)  |
    /// ```
    pub fn rx(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self([1., 0., 0., 0., cos, sin, 0., -sin, cos])
    }

    /// Generate z-rotation matrix
    ///
    /// ```text
    /// | cos(a)  -sin(a)  0 |
    /// | sin(a)  cos(a)   0 |
    /// | 0       0        1 |
    /// ```
    pub fn rz(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self([cos, sin, 0., -sin, cos, 0., 0., 0., 1.])
    }

    /// Compute `Matrix33`*`Vec3`
    pub fn mul_vec(&self, v: &Vec3) -> Vec3 {
        Vec3 {
            x: v.x * self.0[0] + v.y * self.0[3] + v.z * self.0[6],
            y: v.x * self.0[1] + v.y * self.0[4] + v.z * self.0[7],
            z: v.x * self.0[2] + v.y * self.0[5] + v.z * self.0[8],
        }
    }

    /// Compute `Matrix33`*`Matrix33`
    pub fn mul_mat(&self, other: &Self) -> Self {
        let s = &self.0;
        let o = &other.0;
        Self([
            o[0] * s[0] + o[1] * s[3] + o[2] * s[6],
            o[0] * s[1] + o[1] * s[4] + o[2] * s[7],
            o[0] * s[2] + o[1] * s[5] + o[2] * s[8],
            o[3] * s[0] + o[4] * s[3] + o[5] * s[6],
            o[3] * s[1] + o[4] * s[4] + o[5] * s[7],
            o[3] * s[2] + o[4] * s[5] + o[5] * s[8],
            o[6] * s[0] + o[7] * s[3] + o[8] * s[6],
            o[6] * s[1] + o[7] * s[4] + o[8] * s[7],
            o[6] * s[2] + o[7] * s[5] + o[8] * s[8],
        ])
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
            r: val.len(),
        }
    }
}

impl From<Spherical> for Vec3 {
    fn from(s: Spherical) -> Self {
        Self {
            x: s.lat.cos() * s.lon.cos(),
            y: s.lat.cos() * s.lon.sin(),
            z: s.lat.sin(),
        } * s.r
    }
}

impl From<Vec3> for Azimuthal {
    fn from(val: Vec3) -> Self {
        let xy = f64::hypot(val.x, val.y);
        let z = f64::atan2(val.x, val.y);
        Azimuthal {
            z: if z < 0. { z + TAU } else { z },
            h: f64::atan(val.z / xy),
        }
    }
}

impl From<Azimuthal> for Vec3 {
    fn from(a: Azimuthal) -> Self {
        Self {
            x: a.h.cos() * a.z.sin(),
            y: a.h.cos() * a.z.cos(),
            z: a.h.sin(),
        }
    }
}

//
// Unit
//

//
// impl Length for Vec3 {
//     fn len(&self) -> f64 {
//         f64::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
//     }
// }
//
// pub struct Unit<T: Length>(T);
//
// impl<T: Length + ops::DivAssign<f64>> Unit<T> {
//     pub fn new(mut elem: T) -> Self {
//         let len = elem.len();
//         elem /= len;
//         Self(elem)
//     }
//
//     pub fn new_unchecked(elem: T) -> Self {
//         Self(elem)
//     }
//
//     pub fn to_inner(self) -> T {
//         self.0
//     }
// }
//
// impl<T: Length> ops::Deref for Unit<T> {
//     type Target = T;
//
//     fn deref(&self) -> &Self::Target {
//         &self.0
//     }
// }

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aprox_eq::AproxEq;
    use std::f64::consts::*;

    #[test]
    fn quaternion_rotate() {
        let x = Vec3::x_axis();
        let y = Vec3::y_axis();
        let z = Vec3::z_axis();

        let q = UnitQuaternion::new(z, FRAC_PI_2);
        let x_rot = q * x;

        assert!(x_rot.x.aprox_eq(y.x));
        assert!(x_rot.y.aprox_eq(y.y));
        assert!(x_rot.z.aprox_eq(y.z));
    }
}
