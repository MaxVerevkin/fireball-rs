use std::fmt;

use crate::constants::PLUS_MINUS;

#[derive(Debug, Clone, Copy)]
pub struct DisplayWithSigmaPercents {
    value: f64,
    sigma: f64,
}

impl DisplayWithSigmaPercents {
    pub fn new(value: f64, sigma: f64) -> Self {
        Self { value, sigma }
    }
}

impl fmt::Display for DisplayWithSigmaPercents {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rel = self.sigma / self.value;
        let precision = self
            .value
            .is_normal()
            .then(|| (2.0 - f64::min(2.0, self.value.abs().log10().floor())) as usize)
            .unwrap_or(0);
        write!(
            f,
            "{:.*}{PLUS_MINUS}{:.0}%",
            precision,
            self.value,
            rel * 100.0
        )
    }
}
