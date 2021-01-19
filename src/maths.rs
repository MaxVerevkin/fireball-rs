//! Some useful mathematical functions

use crate::structs::*;
use std::f64::consts::{PI, TAU};

/// Compute the smallest difference of two angles
pub fn angle_diff(a1: f64, a2: f64) -> f64 {
    (a2 - a1 + PI).rem_euclid(TAU) - PI
}

/// Compute descent angle
pub fn descent_angle(start: &Azimuthal, end: &Azimuthal) -> Option<f64> {
    let dz = angle_diff(start.z, end.z);
    if dz.abs() < 1e-5 {
        return if (end.h - start.h).abs() < 1e-5 {
            None
        } else {
            Some(if end.h > start.h { 0. } else { PI })
        };
    }

    let cos_l = start.h.sin() * end.h.sin() + start.h.cos() * end.h.cos() * dz.cos();
    let sin_l = (1. - cos_l * cos_l).sqrt();

    let a = ((end.h.sin() - start.h.sin() * cos_l) / (start.h.cos() * sin_l)).acos();
    if a.is_nan() {
        return None;
    }

    Some(if dz > 0. { a } else { TAU - a })
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

    #[test]
    fn descent_angle_test() {
        assert_eq!(
            descent_angle(&Azimuthal { z: 1., h: 0.1 }, &Azimuthal { z: 1., h: 0.2 }).unwrap(),
            0.
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 1.01, h: 1.01 }).unwrap()
                - (PI / 4.))
                < 1e-2
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 1.01, h: 0.3 }).unwrap()
                - (PI / 2.))
                < 1e-5
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 1.01, h: 0.99 }).unwrap()
                - (PI / 4. * 3.))
                < 0.3
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 1., h: 0.2 }).unwrap()
                - (PI))
                < 1e-5
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 0.99, h: 0.99 }).unwrap()
                - (PI / 4. * 5.))
                < 0.3
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 0.3 }, &Azimuthal { z: 0.99, h: 0.3 }).unwrap()
                - (3. * PI / 2.))
                < 1e-2
        );
        assert!(
            (descent_angle(&Azimuthal { z: 1., h: 1. }, &Azimuthal { z: 0.99, h: 1.01 }).unwrap()
                - (PI / 4. * 7.))
                < 0.3
        );
    }
}
