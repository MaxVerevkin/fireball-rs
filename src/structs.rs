//! Basic structures such as Vec3

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

impl Vec3 {
    /// Compute the magnitude of the vector
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Return a unit vector (i.e. with the magnitude of one) with the same direction
    pub fn normalize(&self) -> Self {
        let l = self.length();
        Self {
            x: self.x / l,
            y: self.y / l,
            z: self.z / l,
        }
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

    /// Translate local cartesian coordinates (East, North, Zenith) to global (x, y, z)
    pub fn to_global(&self, pos: &Spherical) -> Self {
        let (east, north, normal) = pos.directions();
        east * self.x + north * self.y + normal * self.z
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

    /// Calculate a tuple of three unit vectors: the direction to the East,
    /// North and Zenith for an observer with a specified spherical location
    pub fn directions(&self) -> (Vec3, Vec3, Vec3) {
        let east = Vec3 {
            x: -self.lon.sin(),
            y: self.lon.cos(),
            z: 0.,
        };
        let north = Vec3 {
            x: -self.lat.sin() * self.lon.cos(),
            y: -self.lat.sin() * self.lon.sin(),
            z: self.lat.cos(),
        };
        let zenith = Vec3 {
            x: self.lat.cos() * self.lon.cos(),
            y: self.lat.cos() * self.lon.sin(),
            z: self.lat.sin(),
        };
        (east, north, zenith)
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
