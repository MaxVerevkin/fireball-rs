//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! The current implementation uses binary search to find two points
//! in the space that represents the trajectory of a fireball. The binary search
//! itself minimizes the mean-square-error of observations given trajectory.

#![feature(core_intrinsics)]
#![allow(clippy::option_map_unit_fn)]

#[macro_use]
extern crate assert_approx_eq;

pub mod constants;
pub mod data;
pub mod maths;
pub mod quick_median;
pub mod solver;
pub mod structs;
pub mod util;
