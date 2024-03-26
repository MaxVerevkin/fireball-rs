//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! The current implementation uses binary search to find two points
//! in the space that represents the trajectory of a fireball. The binary search
//! itself minimizes the mean-square-error of observations given trajectory.

#![feature(const_fn_floating_point_arithmetic)]
// #![allow(clippy::option_map_unit_fn)]

#[macro_use]
extern crate assert_approx_eq;

pub mod approx_eq;
pub mod db;
pub mod solver;
