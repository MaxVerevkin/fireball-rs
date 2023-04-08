use super::{Matrix33, Matrix33Ext, Spherical, UnitQuaternion};
//
// use rand::distributions::Distribution;
use rand::Rng;
use rand_distr::StandardNormal;
//
use std::f64::consts::{FRAC_PI_2, TAU};
use std::intrinsics::unlikely;
// use std::ops;

pub type Vec3 = nalgebra::Vector3<f64>;

pub type UnitVec3 = nalgebra::UnitVector3<f64>;

// TODO: Why do we require unit vectors?????
pub trait UnitVec3Ext {
    /// Translate local cartesian coordinates (East, North, Zenith) to global (x, y, z)
    fn to_global(self, pos: Spherical) -> Self;

    /// Translate global cartesian coordinates (x, y, z) to local (East, North, Zenith)
    fn to_local(self, pos: Spherical) -> Self;

    /// Tilt `self` by a random angle in normal distribution
    fn tilt_random(self, sigma: f64, rng: &mut impl Rng) -> Self;
}

impl UnitVec3Ext for UnitVec3 {
    fn to_global(self, pos: Spherical) -> Self {
        Self::new_unchecked(
            Matrix33::rz(FRAC_PI_2 + pos.lon)
                * Matrix33::rx(FRAC_PI_2 - pos.lat)
                * self.into_inner(),
        )
    }

    fn to_local(self, pos: Spherical) -> Self {
        Self::new_unchecked(
            Matrix33::rx(-FRAC_PI_2 + pos.lat)
                * Matrix33::rz(-FRAC_PI_2 - pos.lon)
                * self.into_inner(),
        )
    }

    fn tilt_random(self, sigma: f64, rng: &mut impl Rng) -> Self {
        let alpha: f64 = rng.sample::<f64, _>(StandardNormal) * sigma;
        let beta = rng.gen::<f64>() * TAU;

        let mut perp = self.cross(&Vec3::x());
        if unlikely(self.dot(&Vec3::x()) < 0.01) {
            perp = self.cross(&Vec3::y());
        }
        perp.normalize_mut();

        UnitQuaternion::new(self.into_inner() * beta)
            * UnitQuaternion::new(perp.normalize() * alpha)
            * self
    }
}

// /// The distance from `self` to a line (`p`, `k`), where `p` is a point on the line and `k` is
// /// the unit vector of direction
// pub fn dist_to_line_sq(self, p: Self, k: Self) -> f64 {
//     k.cross(p - self).len_sq()
// }

// impl From<Vec3> for scad_gen::Vec3 {
//     fn from(vec: Vec3) -> Self {
//         scad_gen::Vec3::new(vec.x as _, vec.y as _, vec.z as _)
//     }
// }

// impl Distribution<Vec3> for rand::distributions::Standard {
//     fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
//         Vec3 {
//             x: rng.gen_range(-1.0..=1.0),
//             y: rng.gen_range(-1.0..=1.0),
//             z: rng.gen_range(-1.0..=1.0),
//         }
//     }
// }
//
// impl Distribution<Vec3> for rand_distr::StandardNormal {
//     fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
//         Vec3 {
//             x: self.sample(rng),
//             y: self.sample(rng),
//             z: self.sample(rng),
//         }
//     }
// }

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
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::approx_eq::ApproxEq;
//     use rand::{rngs::SmallRng, SeedableRng};
//     use std::f64::consts::FRAC_PI_6;
//
//     #[test]
//     fn titl_random() {
//         let mut rng = SmallRng::from_entropy();
//         for _ in 0..10 {
//             let sigma = rng.gen::<f64>() * FRAC_PI_6;
//
//             const N: usize = 15_000;
//             let mut s = 0.0;
//             for _ in 0..N {
//                 let v1: Vec3 = rng.gen();
//                 let v2 = v1.tilt_random(sigma, &mut rng);
//
//                 assert!(v1.len_sq().approx_eq(v2.len_sq()));
//
//                 let a = (v1.dot(v2) / v1.len_sq()).acos();
//                 s += a * a;
//             }
//
//             let s = (s / N as f64).sqrt();
//
//             assert!((s - sigma).abs() / s.max(sigma) < 0.05);
//         }
//     }
// }
