//! This program aims to reconstruct the trajectory of a fireball given
//! only data collected from witnesses.
//!
//! The current implementation uses binary search to find two points
//! in the space that represents the trajectory of a fireball. The binary search
//! itself minimizes the mean-square-error of observations given trajectory.

pub mod constants;
pub mod data;
pub mod maths;
pub mod params;
pub mod solver;
pub mod structs;
