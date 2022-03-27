use super::Vec3;
use std::ops;

/// A (directed) line in space
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Line {
    pub point: Vec3,
    pub direction: Vec3,
}

impl ops::Add for Line {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl ops::AddAssign for Line {
    fn add_assign(&mut self, rhs: Self) {
        self.point += rhs.point;
        self.direction += rhs.direction;
    }
}

impl ops::Mul<f64> for Line {
    type Output = Self;
    fn mul(mut self, rhs: f64) -> Self {
        self *= rhs;
        self
    }
}

impl ops::MulAssign<f64> for Line {
    fn mul_assign(&mut self, rhs: f64) {
        self.point *= rhs;
        self.direction *= rhs;
    }
}
