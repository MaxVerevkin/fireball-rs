//! Compare two objects with some tolerance

/// f64 tolerance.
/// Equals to `f64::EPSILON.sqrt()`
pub const F64_TOLERANCE: f64 = 0.000000014901161193847656_f64;

pub trait ApproxEq {
    /// Compare two objects with some tolerance
    fn approx_eq(self, rhs: Self) -> bool;
}

impl ApproxEq for f64 {
    #[inline]
    fn approx_eq(self, rhs: f64) -> bool {
        (self - rhs).abs() < F64_TOLERANCE
    }
}
