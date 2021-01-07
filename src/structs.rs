//! Basic structures such as Vec3

use std::f64::consts::FRAC_PI_2;
use std::ops::{Add, Div, Mul, Sub};

/// 3D vector
#[derive(Debug, Copy, Clone)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Spherical coordinates triple
#[derive(Debug, Copy, Clone)]
pub struct Spherical {
    /// Positive latitude means North, negative means South
    pub lat: f64,
    /// Positive longitude means East, negative means West
    pub lon: f64,
    /// Relative to the center of the Earth
    pub r: f64,
}

/// Azimuthal coordinates tuple
#[derive(Debug, Copy, Clone)]
pub struct Azimuthal {
    /// Azimuth is equal to zero when points to North,
    /// 90 degrees when points to West, and so on
    pub z: f64,
    pub h: f64,
}

/// 3x3 Matirix
pub struct Matrix33 {
    val: [f64; 9],
}

impl Vec3 {
    /// Compute the magnitude of the vector
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Return a unit vector (i.e. with the magnitude of one) with the same direction
    pub fn normalize(mut self) -> Self {
        let l = self.length();
        self.x /= l;
        self.y /= l;
        self.z /= l;
        self
    }

    /// Translate cartesian coordinates to spherical
    pub fn to_spherical(&self) -> Spherical {
        let xy = f64::hypot(self.x, self.y);
        Spherical {
            lat: f64::atan(self.z / xy),
            lon: f64::atan2(self.y, self.x),
            r: self.length(),
        }
    }

    /// Translate cartesian coordinates to azimuthal
    pub fn to_azimuthal(&self) -> Azimuthal {
        let xy = f64::hypot(self.x, self.y);
        Azimuthal {
            z: f64::atan2(self.x, self.y),
            h: f64::atan(self.z / xy),
        }
    }

    /// Translate local cartesian coordinates (East, North, Zenith) to global (x, y, z)
    pub fn to_global(&self, pos: &Spherical) -> Self {
        Matrix33::rz(FRAC_PI_2 + pos.lon)
            .mul_mat(&Matrix33::rx(FRAC_PI_2 - pos.lat))
            .mul_vec(self)
    }

    /// Translate global cartesian coordinates (x, y, z) to local (East, North, Zenith)
    pub fn to_local(&self, pos: &Spherical) -> Self {
        Matrix33::rx(-FRAC_PI_2 + pos.lat)
            .mul_mat(&Matrix33::rz(-FRAC_PI_2 - pos.lon))
            .mul_vec(self)
    }

    /// Compute the `dot` product of two vectors
    pub fn dot(&self, other: &Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Compute the `cross` product of two vectors
    pub fn cross(&self, other: &Vec3) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}
impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, other: f64) -> Self {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl Spherical {
    /// Translate spherical coordinates to cartesian
    pub fn to_vec3(&self) -> Vec3 {
        let xy = self.r * self.lat.cos();
        Vec3 {
            x: xy * self.lon.cos(),
            y: xy * self.lon.sin(),
            z: self.r * self.lat.sin(),
        }
    }
}

impl Azimuthal {
    /// Translate azimuthal coordinates to an unit cartesian vector
    pub fn to_vec3(&self) -> Vec3 {
        Vec3 {
            x: self.h.cos() * self.z.sin(),
            y: self.h.cos() * self.z.cos(),
            z: self.h.sin(),
        }
    }
}

impl Matrix33 {
    /// Generate Rx matrix
    ///
    /// Rx matrix represents rotation about X-axis
    ///
    /// | 1  0       0       |
    /// | 0  cos(a)  -sin(a) |
    /// | 0  sin(a)  cos(a)  |
    ///
    pub fn rx(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self {
            val: [1., 0., 0., 0., cos, sin, 0., -sin, cos],
        }
    }

    /// Generate Rz matrix
    ///
    /// Rz matrix represents rotation about Z-axis
    ///
    /// | cos(a)  -sin(a)  0 |
    /// | sin(a)  cos(a)   0 |
    /// | 0       0        1 |
    ///
    pub fn rz(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self {
            val: [cos, sin, 0., -sin, cos, 0., 0., 0., 1.],
        }
    }

    /// Compute M*v
    pub fn mul_vec(&self, v: &Vec3) -> Vec3 {
        Vec3 {
            x: v.x * self.val[0] + v.y * self.val[3] + v.z * self.val[6],
            y: v.x * self.val[1] + v.y * self.val[4] + v.z * self.val[7],
            z: v.x * self.val[2] + v.y * self.val[5] + v.z * self.val[8],
        }
    }

    /// Compute M1*M2
    pub fn mul_mat(&self, other: &Self) -> Self {
        Self {
            val: [
                other.val[0] * self.val[0]
                    + other.val[1] * self.val[3]
                    + other.val[2] * self.val[6],
                other.val[0] * self.val[1]
                    + other.val[1] * self.val[4]
                    + other.val[2] * self.val[7],
                other.val[0] * self.val[2]
                    + other.val[1] * self.val[5]
                    + other.val[2] * self.val[8],
                other.val[3] * self.val[0]
                    + other.val[4] * self.val[3]
                    + other.val[5] * self.val[6],
                other.val[3] * self.val[1]
                    + other.val[4] * self.val[4]
                    + other.val[5] * self.val[7],
                other.val[3] * self.val[2]
                    + other.val[4] * self.val[5]
                    + other.val[5] * self.val[8],
                other.val[6] * self.val[0]
                    + other.val[7] * self.val[3]
                    + other.val[8] * self.val[6],
                other.val[6] * self.val[1]
                    + other.val[7] * self.val[4]
                    + other.val[8] * self.val[7],
                other.val[6] * self.val[2]
                    + other.val[7] * self.val[5]
                    + other.val[8] * self.val[8],
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Matrix33, Vec3};
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

        assert!((x.x + FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((x.y - FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((x.z - 0.).abs() < 1e-10);
        assert!((y.x + FRAC_PI_6.sin() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((y.y + FRAC_PI_6.sin() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((y.z - FRAC_PI_6.cos()).abs() < 1e-10);
        assert!((z.x - FRAC_PI_6.cos() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((z.y - FRAC_PI_6.cos() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((z.z - FRAC_PI_6.sin()).abs() < 1e-10);
    }
}
