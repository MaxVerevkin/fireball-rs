//! Basic structures such as Vec3

pub type Vec3 = nalgebra::Vector3<f64>;

// use std::f64::consts::{FRAC_PI_2, TAU};
// use std::ops;

/*

/// 3D vector
// #[derive(Debug, Default, Copy, Clone)]
// pub struct Vec3 {
//     pub x: f64,
//     pub y: f64,
//     pub z: f64,
// }
//
// impl Vec3 {
//     /// Create new vector
//     pub fn new(x: f64, y: f64, z: f64) -> Self {
//         Self { x, y, z }
//     }
//
//     /// Create normalied vector
//     pub fn normalized(self) -> Self {
//         self / self.len()
//     }
//
//     /// Normalize this vector in-place
//     pub fn normalize(&mut self) {
//         *self /= self.len();
//     }
//
//     /// Translate local cartesian coordinates (East, North, Zenith) to global (x, y, z)
//     pub fn to_global(&self, pos: Spherical) -> Self {
//         Matrix33::rz(FRAC_PI_2 + pos.lon)
//             .mul_mat(&Matrix33::rx(FRAC_PI_2 - pos.lat))
//             .mul_vec(self)
//     }
//
//     /// Translate global cartesian coordinates (x, y, z) to local (East, North, Zenith)
//     pub fn to_local(&self, pos: Spherical) -> Self {
//         Matrix33::rx(-FRAC_PI_2 + pos.lat)
//             .mul_mat(&Matrix33::rz(-FRAC_PI_2 - pos.lon))
//             .mul_vec(self)
//     }
//
//     /// Dot product of two vectors
//     pub fn dot(&self, other: Vec3) -> f64 {
//         self.x * other.x + self.y * other.y + self.z * other.z
//     }
//
//     /// Cross product of two vectors
//     pub fn cross(&self, other: Vec3) -> Self {
//         Self {
//             x: self.y * other.z - self.z * other.y,
//             y: self.z * other.x - self.x * other.z,
//             z: self.x * other.y - self.y * other.x,
//         }
//     }
//
//     /// Same as `f64::is_normal()` but accepts zero
//     pub fn is_normal(&self) -> bool {
//         (self.x.is_normal() || self.x == 0.0)
//             && (self.y.is_normal() || self.y == 0.0)
//             && (self.z.is_normal() || self.z == 0.0)
//     }
// }
//
// impl rand::distributions::Distribution<Vec3> for rand::distributions::Standard {
//     fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
//         Vec3 {
//             x: rng.gen_range(-1.0..=1.0),
//             y: rng.gen_range(-1.0..=1.0),
//             z: rng.gen_range(-1.0..=1.0),
//         }
//     }
// }
//
// impl ops::Neg for Vec3 {
//     type Output = Self;
//     fn neg(self) -> Self {
//         Self {
//             x: -self.x,
//             y: -self.y,
//             z: -self.z,
//         }
//     }
// }
//
// impl ops::Add for Vec3 {
//     type Output = Self;
//     fn add(self, other: Self) -> Self {
//         Self {
//             x: self.x + other.x,
//             y: self.y + other.y,
//             z: self.z + other.z,
//         }
//     }
// }
//
// impl ops::AddAssign for Vec3 {
//     fn add_assign(&mut self, other: Self) {
//         self.x += other.x;
//         self.y += other.y;
//         self.z += other.z;
//     }
// }
//
// impl ops::Sub for Vec3 {
//     type Output = Self;
//     fn sub(self, other: Self) -> Self {
//         Self {
//             x: self.x - other.x,
//             y: self.y - other.y,
//             z: self.z - other.z,
//         }
//     }
// }
//
// impl ops::SubAssign for Vec3 {
//     fn sub_assign(&mut self, other: Self) {
//         self.x -= other.x;
//         self.y -= other.y;
//         self.z -= other.z;
//     }
// }
//
// impl ops::Mul<f64> for Vec3 {
//     type Output = Self;
//     fn mul(self, other: f64) -> Self {
//         Self {
//             x: self.x * other,
//             y: self.y * other,
//             z: self.z * other,
//         }
//     }
// }
//
// impl ops::MulAssign<f64> for Vec3 {
//     fn mul_assign(&mut self, other: f64) {
//         self.x *= other;
//         self.y *= other;
//         self.z *= other;
//     }
// }
//
// impl ops::Div<f64> for Vec3 {
//     type Output = Self;
//     fn div(self, other: f64) -> Self {
//         self * (1. / other)
//     }
// }
//
// impl ops::DivAssign<f64> for Vec3 {
//     fn div_assign(&mut self, other: f64) {
//         let k = 1. / other;
//         self.x *= k;
//         self.y *= k;
//         self.z *= k;
//     }
// }

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

    /// Aplly rotation to a vector
    pub fn apply_rotation(self, v: Vec3) -> Vec3 {
        let p = Self { s: 0., v };
        let q = self * p * self.inverse();
        q.v
    }
}

impl ops::Mul for UnitQuaternion {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            s: self.s * rhs.s - self.v.dot(&rhs.v),
            v: self.v * rhs.s + rhs.v * self.s + self.v.cross(&rhs.v),
        }
    }
}

/// Spherical coordinates triple
#[derive(Debug, Clone, Copy)]
pub struct Spherical {
    pub lat: f64,
    pub lon: f64,
    pub r: f64,
}

/// Azimuthal coordinates tuple
#[derive(Debug, Clone, Copy)]
pub struct Azimuthal {
    pub z: f64,
    pub h: f64,
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

// pub trait Norm {
//     fn len(&self) -> f64;
// }
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
    use super::{Matrix33, Vec3};
    use crate::aprox_eq::AproxEq;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6};

    #[test]
    fn rz_multiply_rx_test() {
        let m1 = Matrix33::rz(FRAC_PI_2 + FRAC_PI_3);
        let m2 = Matrix33::rx(FRAC_PI_2 - FRAC_PI_6);
        let m3 = m1.mul_mat(&m2);

        let x = m3.mul_vec(&Vec3 {
            x: 1.,
            y: 0.,
            z: 0.,
        });
        let y = m3.mul_vec(&Vec3 {
            x: 0.,
            y: 1.,
            z: 0.,
        });
        let z = m3.mul_vec(&Vec3 {
            x: 0.,
            y: 0.,
            z: 1.,
        });

        assert!(x.x.aprox_eq(-FRAC_PI_3.sin()));
        assert!(x.y.aprox_eq(FRAC_PI_3.cos()));
        assert!(x.z.aprox_eq(0.0));
        assert!((y.x + FRAC_PI_6.sin() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((y.y + FRAC_PI_6.sin() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((y.z - FRAC_PI_6.cos()).abs() < 1e-10);
        assert!((z.x - FRAC_PI_6.cos() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((z.y - FRAC_PI_6.cos() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((z.z - FRAC_PI_6.sin()).abs() < 1e-10);
    }
}

*/
