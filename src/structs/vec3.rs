use super::{Matrix33, Spherical, UnitQuaternion};

use rand::distributions::Distribution;
use rand::Rng;

use std::f64::consts::{FRAC_PI_2, TAU};
use std::intrinsics::unlikely;
use std::ops;

/// 3D vector
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    /// Create new vector
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// X-Axis vector
    pub fn x_axis() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }

    /// Y-Axis vector
    pub fn y_axis() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }

    /// Z-Axis vector
    pub fn z_axis() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }

    /// The length of a vector squared
    pub fn len_sq(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// The length of a vector
    pub fn len(&self) -> f64 {
        f64::sqrt(self.len_sq())
    }

    /// Create normalied vector
    pub fn normalized(self) -> Self {
        self / self.len()
    }

    /// Normalize this vector in-place
    pub fn normalize(&mut self) {
        *self /= self.len();
    }

    /// Translate local cartesian coordinates (East, North, Zenith) to global (x, y, z)
    pub fn to_global(&self, pos: Spherical) -> Self {
        Matrix33::rz(FRAC_PI_2 + pos.lon)
            .mul_mat(&Matrix33::rx(FRAC_PI_2 - pos.lat))
            .mul_vec(self)
    }

    /// Translate global cartesian coordinates (x, y, z) to local (East, North, Zenith)
    pub fn to_local(&self, pos: Spherical) -> Self {
        Matrix33::rx(-FRAC_PI_2 + pos.lat)
            .mul_mat(&Matrix33::rz(-FRAC_PI_2 - pos.lon))
            .mul_vec(self)
    }

    /// Dot product of two vectors
    pub fn dot(&self, other: Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Cross product of two vectors
    pub fn cross(&self, other: Vec3) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Same as `f64::is_normal()` but accepts zero
    pub fn is_normal(&self) -> bool {
        (self.x.is_normal() || self.x == 0.0)
            && (self.y.is_normal() || self.y == 0.0)
            && (self.z.is_normal() || self.z == 0.0)
    }

    /// The distance from `self` to a line (`p`, `k`), where `p` is a point on the line and `k` is
    /// the unit vector of direction
    pub fn dist_to_line_sq(self, p: Self, k: Self) -> f64 {
        k.cross(p - self).len_sq()
    }

    /// Tilt `self` by a random angle in normal distribution
    pub fn tilt_random(self, sigma: f64, rng: &mut impl Rng) -> Self {
        let alpha: f64 = rng.sample(rand_distr::Normal::new(0.0, sigma).unwrap());
        let beta = rng.gen::<f64>() * TAU;

        let mut perp = self.cross(Self::x_axis());
        if unlikely(!perp.is_normal()) {
            perp = self.cross(Self::y_axis());
        }
        perp.normalize();

        UnitQuaternion::new(self.normalized(), beta) * UnitQuaternion::new(perp, alpha) * self
    }
}

impl From<Vec3> for scad_gen::Vec3 {
    fn from(vec: Vec3) -> Self {
        scad_gen::Vec3::new(vec.x as _, vec.y as _, vec.z as _)
    }
}

impl Distribution<Vec3> for rand::distributions::Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        Vec3 {
            x: rng.gen_range(-1.0..=1.0),
            y: rng.gen_range(-1.0..=1.0),
            z: rng.gen_range(-1.0..=1.0),
        }
    }
}

impl Distribution<Vec3> for rand_distr::StandardNormal {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        Vec3 {
            x: self.sample(rng),
            y: self.sample(rng),
            z: self.sample(rng),
        }
    }
}

impl ops::Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl ops::Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl ops::AddAssign for Vec3 {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl ops::Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl ops::SubAssign for Vec3 {
    fn sub_assign(&mut self, other: Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl ops::Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl ops::MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl ops::Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, other: f64) -> Self {
        self * (1. / other)
    }
}

impl ops::DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, other: f64) {
        let k = 1. / other;
        self.x *= k;
        self.y *= k;
        self.z *= k;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aprox_eq::AproxEq;
    use rand::{rngs::SmallRng, SeedableRng};
    use std::f64::consts::FRAC_PI_6;

    #[test]
    fn titl_random() {
        let mut rng = SmallRng::from_entropy();
        for _ in 0..10 {
            let sigma = rng.gen::<f64>() * FRAC_PI_6;

            const N: usize = 15_000;
            let mut s = 0.0;
            for _ in 0..N {
                let v1: Vec3 = rng.gen();
                let v2 = v1.tilt_random(sigma, &mut rng);

                assert!(v1.len_sq().aprox_eq(v2.len_sq()));

                let a = (v1.dot(v2) / v1.len_sq()).acos();
                s += a * a;
            }

            let s = (s / N as f64).sqrt();

            assert!((s - sigma).abs() / s.max(sigma) < 0.05);
        }
    }
}
