//! Some useful mathematical functions

use crate::structs::*;
use std::f64::consts::*;

pub const fn to_radians(degrees: f64) -> f64 {
    degrees * (PI / 180.0)
}

/// TODO: figure out why this function is less precice than others (is this still true?)
pub fn descent_angle(zenith: UnitVec3, k: UnitVec3, vel: Vec3) -> f64 {
    // let a = k.dot(*zenith).min(1.0).acos().to_degrees(); -> 3deg
    // if a < 5.0 {
    //     dbg!(a);
    // }
    let x = k.cross(*zenith); // -> +0
    let y = x.cross(*k); // -> +0
    f64::atan2(vel.dot(x), vel.dot(y))
}

pub fn angle_diff(a1: f64, a2: f64) -> f64 {
    (a2 - a1 + PI).rem_euclid(TAU) - PI
}

pub fn lambda(observer: Vec3, k: UnitVec3, point: Vec3, v: UnitVec3) -> f64 {
    let p1p0 = point - observer;

    let a = k.dot(*v);
    let b = p1p0.dot(*v);
    let c = k.dot(p1p0);

    (b * c - a * p1p0.norm_squared()) / (a * b - c)
}
