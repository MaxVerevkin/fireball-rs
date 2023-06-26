#![feature(core_intrinsics, let_chains, const_fn_floating_point_arithmetic)]
#![allow(clippy::excessive_precision)]

pub use rand;
pub use rand_distr;
pub use serde;
pub use toml;

pub mod constants;
pub mod histogram;
pub mod maths;
pub mod obs_data;
pub mod one_dim_search;
pub mod plot;
pub mod quick_median;
pub mod structs;
pub mod utils;
