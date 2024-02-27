use std::env;
use std::f64::consts::{FRAC_PI_2, PI};

use common::hilbert;
use common::plot::plotters::style::RED;
use common::{
    obs_data::{Data, RawData},
    plot::{draw_plot_svg_with_named_axes, plotters::style::BLACK, RGBColor},
    structs::Vec3,
};

fn main() {
    let data_sets: Vec<_> = env::args()
        .skip(1)
        .map(|x| RawData::read_from_toml_file(x).finalize())
        .collect();
    if data_sets.is_empty() {
        eprintln!("provide arguments");
        return;
    }

    let points: Vec<_> = data_sets
        .iter()
        .enumerate()
        .flat_map(|(i, data)| get_altitudes(i, data).into_iter())
        .collect();
    if points.is_empty() {
        eprintln!("no points");
        return;
    }

    let mut plot_points: Vec<_> = points
        .iter()
        .map(|p| {
            let c = BLACK;
            (p.observed.to_degrees(), p.calculated.to_degrees(), c, 1.0)
        })
        .collect();

    let mut buckets = [(0.0, 0); 91];
    for &(observed, calculated, _, _) in &plot_points {
        let i = observed.round() as i32;
        assert!(i >= 0 && i <= 90);
        buckets[i as usize].0 += calculated;
        buckets[i as usize].1 += 1;
    }

    for i in 0..91 {
        if buckets[i].1 > 0 {
            let mean = buckets[i].0 / buckets[i].1 as f64;
            plot_points.push((i as f64, mean, RED, 5.0));
        }
    }

    // const CURVES: &[&dyn Curve] = &[&Line, &Cos, &Bezier, &Power, &FiveOverTwo, &Combined];
    const CURVES: &[&dyn Curve] = &[&Bezier];
    for &curve in CURVES {
        let name = curve.name();
        let k = curve.select_k(&points);
        let err = curve.evaluate_k(k, &points);
        println!("--- {name} ---");
        println!("k: {k}");
        println!("err: {err}");

        draw_plot_svg_with_named_axes(
            &format!("plots/altitudes-{name}.svg"),
            &plot_points,
            &[],
            &[],
            &[
                (&|x| x, RGBColor(0, 0, 255)),
                (
                    &|x| curve.predict_calculated(k, x.to_radians()).to_degrees(),
                    RGBColor(255, 0, 0),
                ),
            ],
            "Observed Altitude (degrees)",
            "Calculated Altitude (degrees)",
        )
        .unwrap();
    }
}

#[derive(Debug, Clone, Copy)]
struct Point {
    observed: f64,
    calculated: f64,
    #[allow(dead_code)]
    left_to_right: Option<bool>,
    #[allow(dead_code)]
    exp: u8,
    #[allow(dead_code)]
    data_set: usize,
}

fn get_altitudes(data_set: usize, data: &Data) -> Vec<Point> {
    let Some(traj) = data.answer.and_then(|a| a.traj) else {
        return vec![];
    };

    data.samples
        .iter()
        .filter_map(|s| {
            let observed = s.h0?;
            let k = (traj.point - s.location).normalize();
            let calculated = k.dot(*s.zenith_dir).asin();
            Some(Point {
                observed,
                calculated,
                left_to_right: None,
                exp: s.exp.unwrap_or(1),
                data_set,
            })
        })
        .collect()
}

trait Curve {
    fn name(&self) -> &str;

    fn predict_calculated(&self, k: f64, observed: f64) -> f64;

    fn k_range(&self) -> (f64, f64);

    fn evaluate_k(&self, k: f64, points: &[Point]) -> f64 {
        points
            .iter()
            .map(|p| self.predict_calculated(k, p.observed) - p.calculated)
            .map(|err| err * err)
            .sum()
    }

    fn select_k(&self, points: &[Point]) -> f64 {
        let (mut k, k_max) = self.k_range();
        let step = (k_max - k) / 2_000.0;
        let mut best_k = k;
        let mut best_eval = f64::INFINITY;
        while k <= k_max {
            let eval = self.evaluate_k(k, &points);
            if eval < best_eval {
                best_k = k;
                best_eval = eval;
            }
            k += step;
        }
        best_k
    }
}

struct Line;
impl Curve for Line {
    fn name(&self) -> &str {
        "line"
    }

    fn k_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        observed * k
    }
}

struct Cos;
impl Curve for Cos {
    fn name(&self) -> &str {
        "cos"
    }

    fn k_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        (1.0 - f64::cos(2.0 * k * observed)) * FRAC_PI_2 / (1.0 - f64::cos(k * PI))
    }
}

struct Bezier;
impl Curve for Bezier {
    fn name(&self) -> &str {
        "bezier"
    }

    fn k_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        let p1 = FRAC_PI_2 * hilbert(k, 16);
        let p2 = FRAC_PI_2 * Vec3::new(1.0, 1.0, 0.0);
        let p = |t: f64| t * (2.0 * (1.0 - t) * p1 + t * p2);

        let d = (p1.x * p1.x + observed * (p2.x - 2.0 * p1.x)).sqrt();
        let t1 = (p1.x + d) / (2.0 * p1.x - p2.x);
        let t2 = (p1.x - d) / (2.0 * p1.x - p2.x);

        let e = 0.001;
        let t = if t1 >= -e && t1 <= 1.0 + e {
            t1
        } else if t2 >= -e && t2 <= 1.0 + e {
            t2
        } else {
            dbg!(t1, t2);
            unreachable!()
        };

        p(t).y
    }
}

struct Power;
impl Curve for Power {
    fn name(&self) -> &str {
        "power"
    }

    fn k_range(&self) -> (f64, f64) {
        (0.0, 10.0)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        FRAC_PI_2 * (observed / FRAC_PI_2).powf(k)
    }
}

struct FiveOverTwo;
impl Curve for FiveOverTwo {
    fn name(&self) -> &str {
        "fiveovertwo"
    }

    fn k_range(&self) -> (f64, f64) {
        (0.0, 1.0)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        let x = observed / FRAC_PI_2;
        FRAC_PI_2 * x.sqrt() * (k + (1.0 - k) * x * x)
    }
}

struct Combined;
impl Curve for Combined {
    fn name(&self) -> &str {
        "combined"
    }

    fn k_range(&self) -> (f64, f64) {
        (1.7, 1.9)
    }

    fn predict_calculated(&self, k: f64, observed: f64) -> f64 {
        const P1: Vec3 = hilbert(0.668, 32);
        let bez = bezier_predict_calculated(FRAC_PI_2 * P1, observed);
        1.8216 * (observed / FRAC_PI_2).sqrt() * (1.0 - (1.6335 * observed + 0.65775).tanh()) + bez
    }
}

fn bezier_predict_calculated(p1: Vec3, observed: f64) -> f64 {
    let p1 = FRAC_PI_2 * p1;
    let p2 = FRAC_PI_2 * Vec3::new(1.0, 1.0, 0.0);
    let p = |t: f64| t * (2.0 * (1.0 - t) * p1 + t * p2);

    let d = (p1.x * p1.x + observed * (p2.x - 2.0 * p1.x)).sqrt();
    let t1 = (p1.x + d) / (2.0 * p1.x - p2.x);
    let t2 = (p1.x - d) / (2.0 * p1.x - p2.x);

    let e = 0.001;
    let t = if t1 >= -e && t1 <= 1.0 + e {
        t1
    } else if t2 >= -e && t2 <= 1.0 + e {
        t2
    } else {
        dbg!(t1, t2);
        unreachable!()
    };

    p(t).y
}
