//! This module is just a collection of functions to manipulate matices

use super::*;

pub fn rx(angle: f64) -> [f64; 9] {
    let (sin, cos) = angle.sin_cos();
    [1., 0., 0., 0., cos, sin, 0., -sin, cos]
}

pub fn rz(angle: f64) -> [f64; 9] {
    let (sin, cos) = angle.sin_cos();
    [cos, sin, 0., -sin, cos, 0., 0., 0., 1.]
}

pub fn mul_vec(mat: &[f64; 9], v: &Vec3) -> Vec3 {
    Vec3 {
        x: v.x * mat[0] + v.y * mat[3] + v.z * mat[6],
        y: v.x * mat[1] + v.y * mat[4] + v.z * mat[7],
        z: v.x * mat[2] + v.y * mat[5] + v.z * mat[8],
    }
}

pub fn mul_mat(mat: &[f64; 9], other: &[f64; 9]) -> [f64; 9] {
    [
        other[0] * mat[0] + other[1] * mat[3] + other[2] * mat[6],
        other[0] * mat[1] + other[1] * mat[4] + other[2] * mat[7],
        other[0] * mat[2] + other[1] * mat[5] + other[2] * mat[8],
        other[3] * mat[0] + other[4] * mat[3] + other[5] * mat[6],
        other[3] * mat[1] + other[4] * mat[4] + other[5] * mat[7],
        other[3] * mat[2] + other[4] * mat[5] + other[5] * mat[8],
        other[6] * mat[0] + other[7] * mat[3] + other[8] * mat[6],
        other[6] * mat[1] + other[7] * mat[4] + other[8] * mat[7],
        other[6] * mat[2] + other[7] * mat[5] + other[8] * mat[8],
    ]
}

#[cfg(test)]
mod tests {
    use crate::structs::matrix33::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6};

    #[test]
    fn identity_test() {
        let m1 = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
        let m2 = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
        let m3 = mul_mat(&m1, &m2);
        assert_eq!(m1, m2);
        assert_eq!(m1, m3);
    }

    #[test]
    fn rz_multiply_rx_test() {
        let m1 = rz(FRAC_PI_2 + FRAC_PI_3);
        let m2 = rx(FRAC_PI_2 - FRAC_PI_6);
        let m3 = mul_mat(&m1, &m2);
        assert!((m3[0] + FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((m3[1] - FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((m3[2] - 0.).abs() < 1e-10);
        assert!((m3[3] + FRAC_PI_6.sin() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((m3[4] + FRAC_PI_6.sin() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((m3[5] - FRAC_PI_6.cos()).abs() < 1e-10);
        assert!((m3[6] - FRAC_PI_6.cos() * FRAC_PI_3.cos()).abs() < 1e-10);
        assert!((m3[7] - FRAC_PI_6.cos() * FRAC_PI_3.sin()).abs() < 1e-10);
        assert!((m3[8] - FRAC_PI_6.sin()).abs() < 1e-10);
    }
}
