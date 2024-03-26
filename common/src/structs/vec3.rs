use rand::Rng;
use rand_distr::StandardNormal;
use serde::{Deserialize, Serialize};

use std::f64::consts::TAU;
use std::ops;

#[derive(Serialize, Deserialize, Debug, Default, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn x() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }

    pub fn y() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }

    pub fn z() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }

    pub fn norm_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    pub fn norm_diff(self, self_diff: Self) -> f64 {
        // norm = (x**2 + y**2 + z**2)**0.5
        // norm,i = 0.5 * (x**2 + y**2 + z**2)**(-0.5) * (2x*x,i + 2y*y,i + 2z*z,i)
        //        = dot(self, self_diff) / L
        self.dot(self_diff) / self.norm()
    }

    pub fn normalize_mut(mut self: &mut Self) {
        self /= self.norm();
    }

    pub fn normalize(mut self) -> Self {
        self /= self.norm();
        self
    }

    pub fn normalize_diff(self, self_diff: Self) -> Self {
        // L = (vec.x**2 + vec.y**2 + vec.z**2)**0.5
        // L,i = 0.5 * (vec.x**2 + vec.y**2 + vec.z**2)**(-0.5) * (2*vec.x*vec.x,i + 2.vec.y....)
        //     = dot(vec, vec,i) / L

        // result.x = vec.x / L
        // result.x,i = vec.x,i / L - vec.x / L**2 * L,i
        //            = vec.x,i / L - vec.x / L**3 * dot(vec, vec,i)

        // result,i = vec,i / L - vec / L**3 * dot(vec, vec,i)
        let l = self.norm();
        self_diff / l - self / (l * l * l) * self.dot(self_diff)
    }

    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn cross(self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    // https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    pub fn rotate(self, axis: UnitVec3, angle: f64) -> Self {
        let (sin_a, cos_a) = angle.sin_cos();
        self * cos_a + axis.cross(self) * sin_a + (*axis) * (1.0 - cos_a) * axis.dot(self)
    }

    pub fn rotate_diff(
        self,
        axis: UnitVec3,
        angle: f64,
        self_diff: Vec3,
        axis_diff: Vec3,
        angle_diff: f64,
    ) -> Vec3 {
        let (sin_a, cos_a) = angle.sin_cos(); // sin_a = sin(angle)
        let sin_a_diff = cos_a * angle_diff;
        let cos_a_diff = -sin_a * angle_diff;

        self * cos_a_diff
            + self_diff * cos_a
            + axis.cross(self) * sin_a_diff
            + sin_a * (axis.cross(self_diff) + axis_diff.cross(self))
            + axis_diff * ((1.0 - cos_a) * axis.dot(self))
            + axis
                * ((1.0 - cos_a) * (axis.dot(self_diff) + axis_diff.dot(self))
                    - axis.dot(self) * cos_a_diff)
    }
}

impl ops::Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl ops::Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        rhs * self
    }
}

impl ops::MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl ops::MulAssign<f64> for &mut Vec3 {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl ops::Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl ops::DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl ops::DivAssign<f64> for &mut Vec3 {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl ops::Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl ops::AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl ops::AddAssign for &mut Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl ops::Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl ops::SubAssign for Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl ops::SubAssign for &mut Vec3 {
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

impl ops::Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq)]
#[serde(transparent)]
pub struct UnitVec3 {
    inner: Vec3,
}

impl UnitVec3 {
    pub fn new_unchecked(vec: Vec3) -> Self {
        Self { inner: vec }
    }

    pub fn new_normalize(vec: Vec3) -> Self {
        Self::new_unchecked(vec.normalize())
    }

    pub fn new_and_get(vec: Vec3) -> (Self, f64) {
        let norm = vec.norm();
        (Self::new_unchecked(vec / norm), norm)
    }

    pub fn into_inner(self) -> Vec3 {
        self.inner
    }

    pub fn rotate(self, axis: UnitVec3, angle: f64) -> Self {
        Self::new_unchecked(self.inner.rotate(axis, angle))
    }

    // TODO: do not require unit vector
    pub fn tilt_random(self, sigma: f64, rng: &mut impl Rng) -> Self {
        let alpha: f64 = rng.sample::<f64, _>(StandardNormal) * sigma;
        let beta = rng.gen::<f64>() * TAU;

        let mut perp = self.cross(Vec3::x());
        if self.dot(Vec3::x()) < 0.01 {
            perp = self.cross(Vec3::y());
        }
        let perp = Self::new_normalize(perp);

        self.rotate(perp, alpha).rotate(self, beta)
    }
}

impl ops::Deref for UnitVec3 {
    type Target = Vec3;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl ops::Neg for UnitVec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self { inner: -self.inner }
    }
}

impl ops::Mul<f64> for UnitVec3 {
    type Output = Vec3;

    fn mul(self, rhs: f64) -> Self::Output {
        self.inner * rhs
    }
}

impl ops::Mul<UnitVec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: UnitVec3) -> Self::Output {
        rhs.inner * self
    }
}

impl ops::Add<Vec3> for UnitVec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Self::Output {
        self.inner + rhs
    }
}

impl ops::Add<UnitVec3> for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: UnitVec3) -> Self::Output {
        self + rhs.inner
    }
}
