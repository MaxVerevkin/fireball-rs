//! Some useful mathematical functions

use crate::structs::*;
use std::f64::consts::*;

pub const fn radians(degrees: f64) -> f64 {
    degrees * (PI / 180.0)
}

/// TODO: figure out why this function is less precice than others (is this still true?)
pub fn descent_angle(observer: Vec3, k: Vec3, vel: Vec3) -> f64 {
    let x = k.cross(observer).normalize();
    let y = x.cross(k).normalize();
    let angle = f64::atan2(vel.dot(x), vel.dot(y));
    if angle < 0. {
        angle + TAU
    } else {
        angle
    }
}

pub fn angle_diff(a1: f64, a2: f64) -> f64 {
    (a2 - a1 + PI).rem_euclid(TAU) - PI
}

// pub fn lambda_corssing(p2: Vec3, k2: Vec3, p: Vec3, k: Vec3) -> f64 {
//     let k1k2 = k.dot(k2);
//     ((p2 - p).dot(k - k2 * k1k2)) / (1.0 - k1k2 * k1k2)
// }

pub fn lambda(observer: Vec3, k: UnitVec3, point: Vec3, v: UnitVec3) -> f64 {
    let p1p0 = point - observer;

    let a = k.dot(*v);
    let b = p1p0.dot(*v);
    let c = k.dot(p1p0);

    (b * c - a * p1p0.norm_squared()) / (a * b - c)
}

pub fn lambda_old(observer: Vec3, k: UnitVec3, point: Vec3, v: UnitVec3) -> f64 {
    let p1p0 = point - observer;

    let n = v.cross(p1p0);
    let p = n.cross(*k);
    p1p0.dot(p) / v.dot(-p)
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     // use proptest::prelude::*;
//     use rand::prelude::*;
//
//     #[test]
//     fn lambda_test() {
//         let mut rng = thread_rng();
//         for _ in 0..1000 {
//             let point = rng.gen::<Vec3>() * 100.;
//             let dir = rng.gen::<Vec3>().normalized();
//             if !dir.is_normal() {
//                 continue;
//             }
//
//             let dist = rng.gen_range(0.0..100.0);
//             let point_on_line = point + (dir * dist);
//             let calc_dist = lambda(Vec3::default(), point_on_line, point, dir);
//             assert_approx_eq!(dist, calc_dist);
//         }
//     }
//
//     #[test]
//     fn da_test() {
//         let ob = Spherical {
//             lat: 30f64.to_radians(),
//             lon: 90f64.to_radians(),
//             r: 1.0,
//         };
//         let ob: Vec3 = ob.into();
//
//         let v = Vec3::x_axis();
//         let k = Vec3::y_axis();
//
//         let comp = descent_angle(ob, k, v);
//         assert_approx_eq!(comp, FRAC_PI_2);
//     }
// }
