use super::{UnitVec3, Vec3};
use rand::Rng;
use rand_distr::StandardNormal;
use serde::Serialize;
use std::ops;

/// A (directed) line in space
#[derive(Serialize, Debug, Copy, Clone, PartialEq)]
pub struct Line {
    pub point: Vec3,
    pub direction: UnitVec3,
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct LineGrad {
    pub point: Vec3,
    pub direction: Vec3,
}

impl Line {
    pub fn adjust_position(
        &mut self,
        rng: &mut impl Rng,
        sigma: f64,
        current_value: &mut f64,
        heuristic: impl FnOnce(Self) -> f64,
    ) {
        let mut copy = *self;
        copy.point.x += sigma * rng.sample::<f64, _>(StandardNormal);
        copy.point.y += sigma * rng.sample::<f64, _>(StandardNormal);
        copy.point.z += sigma * rng.sample::<f64, _>(StandardNormal);
        let new_value = heuristic(copy);
        if new_value < *current_value {
            *self = copy;
            *current_value = new_value;
        }
    }

    pub fn adjust_direction(
        &mut self,
        rng: &mut impl Rng,
        sigma: f64,
        current_value: &mut f64,
        heuristic: impl FnOnce(Self) -> f64,
    ) {
        let mut copy = *self;
        copy.direction = copy.direction.tilt_random(sigma, rng);
        let new_value = heuristic(copy);
        if new_value < *current_value {
            *self = copy;
            *current_value = new_value;
        }
    }

    pub fn adjust_both(
        &mut self,
        rng: &mut impl Rng,
        sigma_point: f64,
        sigma_vel: f64,
        current_value: &mut f64,
        heuristic: impl FnOnce(Self) -> f64,
    ) {
        let mut copy = *self;
        copy.point.x += sigma_point * rng.sample::<f64, _>(StandardNormal);
        copy.point.y += sigma_point * rng.sample::<f64, _>(StandardNormal);
        copy.point.z += sigma_point * rng.sample::<f64, _>(StandardNormal);
        copy.direction = copy.direction.tilt_random(sigma_vel, rng);
        let new_value = heuristic(copy);
        if new_value < *current_value {
            *self = copy;
            *current_value = new_value;
        }
    }

    pub fn offset_point(self, offset: Vec3) -> Self {
        Self {
            point: self.point + offset,
            direction: self.direction,
        }
    }

    pub fn offset_direction(self, offset: Vec3) -> Self {
        Self {
            point: self.point,
            direction: UnitVec3::new_normalize(self.direction + offset),
        }
    }
}

impl ops::AddAssign for LineGrad {
    fn add_assign(&mut self, rhs: Self) {
        self.point += rhs.point;
        self.direction += rhs.direction;
    }
}

impl ops::Mul<f64> for LineGrad {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            point: self.point * rhs,
            direction: self.direction * rhs,
        }
    }
}

impl ops::Div<f64> for LineGrad {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            point: self.point / rhs,
            direction: self.direction / rhs,
        }
    }
}
