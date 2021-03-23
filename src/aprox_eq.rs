//! Compare two objects with some tolerance

/// f64 tolerance.
/// Equals to `f64::EPSILON.sqrt()`
pub const F64_TOLERANCE: f64 = 0.000000014901161193847656_f64;

pub trait AproxEq {
    /// Compare two objects with some tolerance
    fn aprox_eq(self, rhs: Self) -> bool;
}

impl AproxEq for f64 {
    #[inline]
    fn aprox_eq(self, rhs: f64) -> bool {
        (self - rhs).abs() < F64_TOLERANCE
    }
}
