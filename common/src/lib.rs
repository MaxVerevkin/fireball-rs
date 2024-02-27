#![feature(let_chains, const_fn_floating_point_arithmetic, const_mut_refs)]
#![allow(clippy::excessive_precision)]

use std::ops;

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

use plotters::style::RGBColor;
use structs::Vec3;

pub const COLORS: &[RGBColor] = {
    const fn rgb(c: u32) -> RGBColor {
        let [_, r, g, b] = c.to_be_bytes();
        RGBColor(r, g, b)
    }
    &[
        rgb(0x000000),
        rgb(0x00FF00),
        rgb(0x0000FF),
        rgb(0xFF0000),
        rgb(0x01FFFE),
        rgb(0xFFA6FE),
        rgb(0xFFDB66),
        rgb(0x006401),
        rgb(0x010067),
        rgb(0x95003A),
        rgb(0x007DB5),
        rgb(0xFF00F6),
        rgb(0xFFEEE8),
        rgb(0x774D00),
        rgb(0x90FB92),
        rgb(0x0076FF),
        rgb(0xFF937E),
        rgb(0x6A826C),
        rgb(0xFF029D),
        rgb(0xFE8900),
        rgb(0x7A4782),
        rgb(0x7E2DD2),
        rgb(0x85A900),
        rgb(0xFF0056),
        rgb(0xA42400),
        rgb(0x00AE7E),
        rgb(0x683D3B),
        rgb(0xBDC6FF),
        rgb(0x263400),
        rgb(0xBDD393),
        rgb(0x00B917),
        rgb(0x9E008E),
        rgb(0x001544),
        rgb(0xC28C9F),
        rgb(0xFF74A3),
        rgb(0x01D0FF),
        rgb(0x004754),
        rgb(0xE56FFE),
        rgb(0x788231),
        rgb(0x0E4CA1),
        rgb(0x91D0CB),
        rgb(0xBE9970),
        rgb(0x968AE8),
        rgb(0xBB8800),
        rgb(0x43002C),
        rgb(0xDEFF74),
        rgb(0x00FFC6),
        rgb(0xFFE502),
        rgb(0x620E00),
        rgb(0x008F9C),
        rgb(0x98FF52),
        rgb(0x7544B1),
        rgb(0xB500FF),
        rgb(0x00FF78),
        rgb(0xFF6E41),
        rgb(0x005F39),
        rgb(0x6B6882),
        rgb(0x5FAD4E),
        rgb(0xA75740),
        rgb(0xA5FFD2),
        rgb(0xFFB167),
        rgb(0x009BFF),
        rgb(0xE85EBE),
        rgb(0xD5FF00),
    ]
};

// This version makes the program up to 3x faster.
pub const fn hilbert(mut a: f64, mut iters: u32) -> Vec3 {
    let mut cur = (0.5, 0.5);
    let mut x = (0.25, 0.0);
    let mut y = (0.0, 0.25);
    while iters > 0 {
        iters -= 1;
        a *= 4.0;
        if a < 1.0 {
            cur.0 -= x.0;
            cur.1 -= x.1;
            cur.0 += y.0;
            cur.1 += y.1;
            x = (x.1 * 0.5, x.0 * -0.5);
            y = (y.1 * -0.5, y.0 * 0.5);
        } else if a < 2.0 {
            a -= 1.0;
            cur.0 -= x.0;
            cur.1 -= x.1;
            cur.0 -= y.0;
            cur.1 -= y.1;
            x = (x.1 * 0.5, x.0 * 0.5);
            y = (y.1 * 0.5, y.0 * 0.5);
        } else if a < 3.0 {
            a -= 2.0;
            cur.0 += x.0;
            cur.1 += x.1;
            cur.0 -= y.0;
            cur.1 -= y.1;
            x = (x.1 * 0.5, x.0 * 0.5);
            y = (y.1 * 0.5, y.0 * 0.5);
        } else {
            a -= 3.0;
            cur.0 += x.0;
            cur.1 += x.1;
            cur.0 += y.0;
            cur.1 += y.1;
            x = (x.1 * -0.5, x.0 * 0.5);
            y = (y.1 * 0.5, y.0 * -0.5);
        }
    }
    Vec3::new(cur.0, cur.1, 0.0)
}

pub const fn old_hilbert(a: f64, iters: u32) -> Vec3 {
    if iters == 0 {
        Vec3::new(0.5, 0.5, 0.0)
    } else if a < 0.25 {
        let h = hilbert(a * 4.0, iters - 1);
        Vec3::new(0.5 - h.y * 0.5, 1.0 - h.x * 0.5, 0.0)
    } else if a < 0.5 {
        let h = hilbert(a * 4.0 - 1.0, iters - 1);
        Vec3::new(h.x * 0.5, h.y * 0.5, 0.0)
    } else if a < 0.75 {
        let h = hilbert(a * 4.0 - 2.0, iters - 1);
        Vec3::new(0.5 + h.x * 0.5, h.y * 0.5, 0.0)
    } else {
        let h = hilbert(a * 4.0 - 3.0, iters - 1);
        Vec3::new(0.5 + h.y * 0.5, 0.5 + h.x * 0.5, 0.0)
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Sigmas {
    /// Meters
    pub x: f64,
    /// Meters
    pub y: f64,
    /// Meters
    pub z: f64,
    /// Radians
    pub v_angle: f64,
    /// Meters/Seconds
    pub speed: f64,
}

impl Sigmas {
    pub fn sqrt(self) -> Self {
        Self {
            x: self.x.sqrt(),
            y: self.y.sqrt(),
            z: self.z.sqrt(),
            v_angle: self.v_angle.sqrt(),
            speed: self.speed.sqrt(),
        }
    }

    pub fn add_point(&mut self, p_diff: Vec3, error: f64) {
        self.x += (p_diff.x * error).powi(2);
        self.y += (p_diff.y * error).powi(2);
        self.z += (p_diff.z * error).powi(2);
    }
}

impl ops::AddAssign for Sigmas {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
        self.v_angle += rhs.v_angle;
        self.speed += rhs.speed;
    }
}

impl ops::Add for Sigmas {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl std::iter::Sum for Sigmas {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |s, x| s + x)
    }
}
