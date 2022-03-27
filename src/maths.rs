//! Some useful mathematical functions

use crate::structs::*;
use std::f64::consts::{PI, TAU};

pub fn descent_angle(observer: Vec3, k: Vec3, vel: Vec3) -> f64 {
    let x = k.cross(observer).normalized();
    let y = x.cross(k).normalized();
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

pub fn lambda_corssing(p2: Vec3, k2: Vec3, p: Vec3, k: Vec3) -> f64 {
    let k1k2 = k.dot(k2);
    ((p2 - p).dot(k - k2 * k1k2)) / (1.0 - k1k2 * k1k2)
}

pub fn lambda(observer: Vec3, observation: Vec3, point: Vec3, vel: Vec3) -> f64 {
    let p1p0 = point - observer;
    let n = vel.cross(p1p0);
    let p = n.cross(observation);
    p1p0.dot(p) / vel.dot(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    // use proptest::prelude::*;
    use rand::prelude::*;

    #[test]
    fn lambda_test() {
        for _ in 0..100 {
            let p = random::<Vec3>() * 10.;
            let k = random::<Vec3>().normalized();
            let p2 = random::<Vec3>() * 10.;
            let k2 = random::<Vec3>().normalized();
            let l = lambda_corssing(p, k, p2, k2);

            let dist = (p + k * l).dist_to_line_sq(p2, k2);
            let dist_1 = (p + k * (l + 1.)).dist_to_line_sq(p2, k2);
            let dist_2 = (p + k * (l - 1.)).dist_to_line_sq(p2, k2);

            assert!(dist <= dist_1);
            assert!(dist <= dist_2);
        }
    }
}
