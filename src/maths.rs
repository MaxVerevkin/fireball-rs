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

/// Compute the smallest difference of two angles
pub fn angle_diff(a1: f64, a2: f64) -> f64 {
    (a2 - a1 + PI).rem_euclid(TAU) - PI
}

fn quick_sort(list: &mut [f64], target_index: usize) {
    let len = list.len();

    if len < 2 {
        return;
    } else if len == 2 {
        if list[0] > list[1] {
            list.swap(0, 1);
        }
        return;
    } else if len == 3 {
        if list[0] > list[1] {
            list.swap(0, 1);
        }
        if list[0] > list[2] {
            list.swap(0, 2);
        }
        if list[1] > list[2] {
            list.swap(1, 2);
        }
        return;
    }

    let piv = list[0];
    let mut left = 1;
    let mut right = len - 1;

    while left <= right {
        if list[left] <= piv {
            left += 1;
        } else {
            list.swap(left, right);
            right -= 1;
        }
    }

    list.swap(0, left - 1);

    if target_index <= left {
        quick_sort(&mut list[..(left - 1)], target_index);
    }
    if target_index >= right {
        quick_sort(
            &mut list[(right + 1)..],
            target_index.saturating_sub(right + 1),
        );
    }
}

pub fn quick_median(list: &mut [f64]) -> f64 {
    let len = list.len();
    if len == 0 {
        f64::INFINITY
    } else if len == 1 {
        list[0]
    } else if len == 2 {
        (list[0] + list[1]) * 0.5
    } else {
        quick_sort(list, len / 2);
        if len % 2 == 0 {
            (list[len / 2] + list[len / 2 - 1]) * 0.5
        } else {
            list[len / 2]
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::maths::*;

    #[test]
    fn angle_diff_test() {
        assert_eq!(angle_diff(PI, TAU).abs(), PI);
        assert_eq!(angle_diff(PI, 10. * TAU).abs(), PI);
        assert_eq!(angle_diff(PI, 0.).abs(), PI);

        assert!(angle_diff(PI + 0.1, TAU) - (PI - 0.1) < 1e-10);
        assert!(angle_diff(PI - 0.1, 0.) - (-PI + 0.1) < 1e-10);

        assert!(angle_diff(PI - 1., -0.1) - (-1.1) < 1e-10);
        assert!(angle_diff(PI + 1., TAU - 0.1) - (PI - 1.1) < 1e-10);
    }

    //#[test]
    //fn descent_angle_test() {
    //assert_eq!(
    //descent_angle(Azimuthal { z: 1., h: 0.1 }, Azimuthal { z: 1., h: 0.2 }),
    //Some(0.)
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 1. }, Azimuthal { z: 1.01, h: 1.01 }).unwrap()
    //- (PI / 4.))
    //< 1e-2
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 0.3 }, Azimuthal { z: 1.01, h: 0.3 }).unwrap()
    //- (PI / 2.))
    //< 1e-5
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 1. }, Azimuthal { z: 1.01, h: 0.99 }).unwrap()
    //- (PI / 4. * 3.))
    //< 0.3
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 0.3 }, Azimuthal { z: 1., h: 0.2 }).unwrap()
    //- (PI))
    //< 1e-5
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 1. }, Azimuthal { z: 0.99, h: 0.99 }).unwrap()
    //- (PI / 4. * 5.))
    //< 0.3
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 0.3 }, Azimuthal { z: 0.99, h: 0.3 }).unwrap()
    //- (3. * PI / 2.))
    //< 1e-2
    //);
    //assert!(
    //(descent_angle(Azimuthal { z: 1., h: 1. }, Azimuthal { z: 0.99, h: 1.01 }).unwrap()
    //- (PI / 4. * 7.))
    //< 0.3
    //);
    //}
}
