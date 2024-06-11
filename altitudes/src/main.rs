use std::env;
use std::f64::consts::{FRAC_PI_2, PI};

use common::hilbert;
use common::plot::plotters::style::*;
use common::structs::UnitVec3;
use common::{
    obs_data::{Data, RawData},
    plot::draw_plot_svg_with_named_axes,
    quick_median::*,
    structs::Vec3,
};

fn _main() {
    let data_sets: Vec<_> = env::args()
        .skip(1)
        .map(|x| {
            RawData::read_from_toml_file(x)
                .da_correction(0.12)
                .az_correction(0.21)
                .altitudes_correction()
                .finalize()
        })
        .collect();
    if data_sets.is_empty() {
        eprintln!("provide arguments");
        return;
    }

    let mut points = Vec::new();

    let db = fireball::db::Db::read(Some("../db.json")).unwrap();
    let key = "da_k=0.12;az_k=0.21;";

    for (i, data) in data_sets.iter().enumerate() {
        if let Some(event) = db.events.get(data.name.as_deref().unwrap()) {
            if let Some(run) = event.runs.get(key) {
                let traj = run.last().unwrap().traj;
                points.extend(get_altitudes(i, traj.point, traj.direction, data));
            }
        }
    }

    if points.is_empty() {
        eprintln!("no points");
        return;
    }

    let mut end_plot_points: Vec<_> = points
        .iter()
        .filter_map(|p| {
            let c = BLACK;
            Some((
                p.end_lambda_rad?.to_degrees(),
                p.end_height_rad?.to_degrees(),
                c,
                1.0,
            ))
        })
        .collect();
    let mut start_plot_points: Vec<_> = points
        .iter()
        .filter_map(|p| {
            let c = RED;
            Some((
                p.start_lambda_rad?.to_degrees(),
                p.start_height_rad?.to_degrees(),
                c,
                1.0,
            ))
        })
        .collect();

    let med_end_lambda_deg = end_plot_points.median_by_key(|p| p.0);
    let med_end_height_deg = end_plot_points.median_by_key(|p| p.1);

    let med_start_lambda_deg = start_plot_points.median_by_key(|p| p.0);
    let med_start_height_deg = start_plot_points.median_by_key(|p| p.1);

    // dbg!(med_end_lambda_deg - med_start_lambda_deg);
    // let mut dls: Vec<_> = points
    //     .iter()
    //     .filter_map(|p| Some(p.end_lambda_rad?.to_degrees() - p.start_lambda_rad?.to_degrees()))
    //     .collect();
    // dbg!(dls.median());

    let mut end_distances = end_plot_points
        .iter()
        .map(|p| {
            ((p.0 - med_end_lambda_deg).to_radians().cos()
                * (p.1 - med_end_height_deg).to_radians().cos())
            .acos()
            .to_degrees()
        })
        .collect::<Vec<_>>();
    end_distances.sort_unstable_by(f64::total_cmp);
    let end_dist_sigma = end_distances[end_distances.len() * 68 / 100];

    let mut start_distances = start_plot_points
        .iter()
        .map(|p| {
            ((p.0 - med_start_lambda_deg).to_radians().cos()
                * (p.1 - med_start_height_deg).to_radians().cos())
            .acos()
            .to_degrees()
        })
        .collect::<Vec<_>>();
    start_distances.sort_unstable_by(f64::total_cmp);
    let start_dist_sigma = start_distances[start_distances.len() * 68 / 100];

    let points_to_draw = start_plot_points
        .into_iter()
        .chain(end_plot_points)
        .collect::<Vec<_>>();

    draw_plot_svg_with_named_axes(
        &format!("plots/lambda-height-ours.svg"),
        &points_to_draw,
        &[
            (0.0, GREEN),
            (med_end_lambda_deg, BLACK),
            (med_start_lambda_deg, RED),
        ],
        &[
            (0.0, GREEN),
            (med_end_height_deg, BLACK),
            (med_start_height_deg, RED),
        ],
        &[
            (
                &|x| {
                    let dx = x - med_end_lambda_deg;
                    if dx <= -end_dist_sigma || dx >= end_dist_sigma {
                        med_end_height_deg
                    } else {
                        let dy = (end_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_end_height_deg + dy
                    }
                },
                BLACK,
            ),
            (
                &|x| {
                    let dx = x - med_end_lambda_deg;
                    if dx <= -end_dist_sigma || dx >= end_dist_sigma {
                        med_end_height_deg
                    } else {
                        let dy = (end_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_end_height_deg - dy
                    }
                },
                BLACK,
            ),
            (
                &|x| {
                    let dx = x - med_start_lambda_deg;
                    if dx <= -start_dist_sigma || dx >= start_dist_sigma {
                        med_start_height_deg
                    } else {
                        let dy = (start_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_start_height_deg + dy
                    }
                },
                RED,
            ),
            (
                &|x| {
                    let dx = x - med_start_lambda_deg;
                    if dx <= -start_dist_sigma || dx >= start_dist_sigma {
                        med_start_height_deg
                    } else {
                        let dy = (start_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_start_height_deg - dy
                    }
                },
                RED,
            ),
        ],
        "Lambda (degrees)",
        "Height (degrees)",
    )
    .unwrap();

    // let mut plot_points: Vec<_> = points
    //     .iter()
    //     .map(|p| {
    //         let c = BLACK;
    //         (p.observed.to_degrees(), p.calculated.to_degrees(), c, 1.0)
    //     })
    //     .collect();

    // let mut buckets = [(0.0, 0); 91];
    // for &(observed, calculated, _, _) in &plot_points {
    //     let i = observed.round() as i32;
    //     assert!(i >= 0 && i <= 90);
    //     buckets[i as usize].0 += calculated;
    //     buckets[i as usize].1 += 1;
    // }

    // for i in 0..91 {
    //     if buckets[i].1 > 0 {
    //         let mean = buckets[i].0 / buckets[i].1 as f64;
    //         plot_points.push((i as f64, mean, RED, 5.0));
    //     }
    // }

    // // const CURVES: &[&dyn Curve] = &[&Line, &Cos, &Bezier, &Power, &FiveOverTwo, &Combined];
    // const CURVES: &[&dyn Curve] = &[&Bezier];
    // for &curve in CURVES {
    //     let name = curve.name();
    //     let k = curve.select_k(&points);
    //     let err = curve.evaluate_k(k, &points);
    //     println!("--- {name} ---");
    //     println!("k: {k}");
    //     println!("err: {err}");

    //     draw_plot_svg_with_named_axes(
    //         &format!("plots/altitudes-{name}.svg"),
    //         &plot_points,
    //         &[],
    //         &[],
    //         &[
    //             (&|x| x, RGBColor(0, 0, 255)),
    //             (
    //                 &|x| curve.predict_calculated(k, x.to_radians()).to_degrees(),
    //                 RGBColor(255, 0, 0),
    //             ),
    //         ],
    //         "Observed Altitude (degrees)",
    //         "Calculated Altitude (degrees)",
    //     )
    //     .unwrap();
    // }
}

fn main() {
    let data_sets: Vec<_> = env::args()
        .skip(1)
        .map(|x| {
            RawData::read_from_toml_file(x)
                .da_correction(0.12)
                .az_correction(0.21)
                .altitudes_correction()
                .finalize()
        })
        .collect();
    if data_sets.is_empty() {
        eprintln!("provide arguments");
        return;
    }

    let points: Vec<_> = data_sets
        .iter()
        .filter_map(|data| Some((data, data.answer.and_then(|a| a.traj)?)))
        .enumerate()
        .flat_map(|(i, (data, traj))| {
            get_altitudes(i, traj.point, traj.direction, data).into_iter()
        })
        .collect();
    if points.is_empty() {
        eprintln!("no points");
        return;
    }

    let mut end_plot_points: Vec<_> = points
        .iter()
        .filter_map(|p| {
            let c = BLACK;
            Some((
                p.end_lambda_rad?.to_degrees(),
                p.end_height_rad?.to_degrees(),
                c,
                1.0,
            ))
        })
        .collect();
    let mut start_plot_points: Vec<_> = points
        .iter()
        .filter_map(|p| {
            let c = RED;
            Some((
                p.start_lambda_rad?.to_degrees(),
                p.start_height_rad?.to_degrees(),
                c,
                1.0,
            ))
        })
        .collect();

    let med_end_lambda_deg = end_plot_points.median_by_key(|p| p.0);
    let med_end_height_deg = end_plot_points.median_by_key(|p| p.1);

    let med_start_lambda_deg = start_plot_points.median_by_key(|p| p.0);
    let med_start_height_deg = start_plot_points.median_by_key(|p| p.1);

    let mut end_distances = end_plot_points
        .iter()
        .map(|p| {
            ((p.0 - med_end_lambda_deg).to_radians().cos()
                * (p.1 - med_end_height_deg).to_radians().cos())
            .acos()
            .to_degrees()
        })
        .collect::<Vec<_>>();
    end_distances.sort_unstable_by(f64::total_cmp);
    let end_dist_sigma = end_distances[end_distances.len() * 68 / 100];

    let mut start_distances = start_plot_points
        .iter()
        .map(|p| {
            ((p.0 - med_start_lambda_deg).to_radians().cos()
                * (p.1 - med_start_height_deg).to_radians().cos())
            .acos()
            .to_degrees()
        })
        .collect::<Vec<_>>();
    start_distances.sort_unstable_by(f64::total_cmp);
    let start_dist_sigma = start_distances[start_distances.len() * 68 / 100];

    // let points_to_draw = start_plot_points;
    let points_to_draw = end_plot_points;
    // let points_to_draw = start_plot_points
    //     .into_iter()
    //     .chain(end_plot_points)
    //     .collect::<Vec<_>>();

    draw_plot_svg_with_named_axes(
        &format!("plots/lambda-height.svg"),
        &points_to_draw,
        &[
            (0.0, GREEN),
            (med_end_lambda_deg, BLACK),
            (med_start_lambda_deg, RED),
        ],
        &[
            (0.0, GREEN),
            (med_end_height_deg, BLACK),
            (med_start_height_deg, RED),
        ],
        &[
            (
                &|x| {
                    let dx = x - med_end_lambda_deg;
                    if dx <= -end_dist_sigma || dx >= end_dist_sigma {
                        med_end_height_deg
                    } else {
                        let dy = (end_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_end_height_deg + dy
                    }
                },
                BLACK,
            ),
            (
                &|x| {
                    let dx = x - med_end_lambda_deg;
                    if dx <= -end_dist_sigma || dx >= end_dist_sigma {
                        med_end_height_deg
                    } else {
                        let dy = (end_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_end_height_deg - dy
                    }
                },
                BLACK,
            ),
            (
                &|x| {
                    let dx = x - med_start_lambda_deg;
                    if dx <= -start_dist_sigma || dx >= start_dist_sigma {
                        med_start_height_deg
                    } else {
                        let dy = (start_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_start_height_deg + dy
                    }
                },
                RED,
            ),
            (
                &|x| {
                    let dx = x - med_start_lambda_deg;
                    if dx <= -start_dist_sigma || dx >= start_dist_sigma {
                        med_start_height_deg
                    } else {
                        let dy = (start_dist_sigma.to_radians().cos() / dx.to_radians().cos())
                            .acos()
                            .to_degrees();
                        med_start_height_deg - dy
                    }
                },
                RED,
            ),
        ],
        "Lambda (degrees)",
        "Height (degrees)",
    )
    .unwrap();
}

#[derive(Debug, Clone, Copy)]
struct Point {
    #[allow(dead_code)]
    observed: f64,
    #[allow(dead_code)]
    calculated: f64,

    end_lambda_rad: Option<f64>,
    end_height_rad: Option<f64>,
    start_lambda_rad: Option<f64>,
    start_height_rad: Option<f64>,

    #[allow(dead_code)]
    left_to_right: Option<bool>,
    #[allow(dead_code)]
    exp: u8,
    #[allow(dead_code)]
    data_set: usize,
}

fn get_altitudes(
    data_set: usize,
    traj_point: Vec3,
    traj_direction: UnitVec3,
    data: &Data,
) -> Vec<Point> {
    data.samples
        .iter()
        .filter_map(|s| {
            let observed = s.h0?;
            let k_calculated = UnitVec3::new_normalize(traj_point - s.location);
            let calculated = k_calculated.dot(*s.zenith_dir).asin();

            let p_norm = k_calculated.cross(*traj_direction).normalize();
            let p_norm_up = p_norm * p_norm.dot(*s.zenith_dir).signum();

            let (end_lambda, end_height) = s
                .k_end
                .map(|k| {
                    (
                        p_norm.cross(*k_calculated).normalize().dot(*k).asin(),
                        p_norm_up.dot(*k).asin(),
                    )
                })
                .unzip();
            let (start_lambda, start_height) = s
                .k_start
                .map(|k| {
                    (
                        p_norm.cross(*k_calculated).normalize().dot(*k).asin(),
                        p_norm_up.dot(*k).asin(),
                    )
                })
                .unzip();

            Some(Point {
                observed,
                calculated,

                end_lambda_rad: end_lambda,
                end_height_rad: end_height,
                start_lambda_rad: start_lambda,
                start_height_rad: start_height,

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
