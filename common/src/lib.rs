#![feature(core_intrinsics, let_chains)]

pub mod constants;
pub mod histogram;
pub mod maths;
pub mod obs_data;
pub mod plot;
pub mod quick_median;
pub mod structs;

pub use nalgebra;
pub use rand;
pub use rand_distr;
pub use serde;
pub use toml;
