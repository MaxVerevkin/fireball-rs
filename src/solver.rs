//! The implementation

use rand::{Rng, SeedableRng};
use std::f64::consts::FRAC_PI_2;

use crate::constants::*;
use crate::data::{Data, DataSample};
use crate::maths::*;
use crate::quick_median::quick_median;
use crate::structs::*;

pub mod pair;
use pair::PairTrajectory;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Parameters to tweak the algorithm
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Whether to tilt the observation to better match given descent angle
    pub da_correction: f64,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub flash: Vec3,
    pub velocity: Vec3,
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver { data, params }
    }

    /// The weighted sum of pairwise plane crossings
    pub fn pairwise(&self) -> Line {
        let mut traj = Line::default();
        let mut sum = 0.0;
        for (i, s1) in self.data.iter().enumerate() {
            for s2 in self.data.iter().skip(i + 1) {
                if let Some(line) = PairTrajectory::calculate(s1, s2) {
                    traj += line.line * line.weight;
                    sum += line.weight;
                }
            }
        }
        traj.point /= sum;
        traj.direction.normalize();
        traj
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        // TODO explain what's going on here
        if self.params.da_correction > 0.0 {
            for s in &mut self.data.samples {
                if let Some(((da, axis), comp_da)) = s.da.zip(s.axis).zip(s.calculated_da()) {
                    let angle = angle_diff(comp_da, da) * self.params.da_correction;
                    let q = UnitQuaternion::new(axis, angle);
                    s.plane.as_mut().map(|k| *k = q * *k);
                    s.k_start.as_mut().map(|k| *k = q * *k);
                    s.k_end.as_mut().map(|k| *k = q * *k);
                }
            }
        }

        // The initial guess is the weighted sum of pairwise plane crossings
        let mut traj = self.pairwise();
        self.data.compare(traj, "Initial guess");

        traj.direction *= self.calc_flash_and_speed(traj).1;

        let mut buf = Vec::new();
        let mut eval_lms = |traj: Line| -> f64 {
            buf.clear();
            for s in &self.data.samples {
                let (eb, ee, eda) = Self::evaluate_traj(s, traj);
                eb.map(|e| buf.push(e));
                ee.map(|e| buf.push(e));
                eda.map(|e| buf.push(e));
            }
            quick_median(&mut buf)
        };

        let mut rng = rand::rngs::SmallRng::from_entropy();
        let mut best_err = eval_lms(traj);
        let mut dif = 20_000.0;
        for _ in 0..3_000 {
            let distr = rand_distr::Normal::new(0.0, dif).unwrap();
            for _ in 0..100 {
                let mut t = traj;
                t.point.x += rng.sample(distr);
                t.point.y += rng.sample(distr);
                t.point.z += rng.sample(distr);
                let e = eval_lms(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
            dif *= 0.999;

            let distr = rand_distr::Normal::new(0.0, 500.0).unwrap();
            for _ in 0..100 {
                let mut t = traj;
                t.direction.x += rng.sample(distr);
                t.direction.y += rng.sample(distr);
                t.direction.z += rng.sample(distr);
                let e = eval_lms(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
        }

        self.data.compare(traj, "After LMS");

        let best_err = best_err.sqrt() * 1.483 * (1. + 5. / (self.data.samples.len() - 6) as f64);
        let get_weight = |err: f64| -> f64 {
            // Smooth transition from 1 to 0
            // Visualization: https://www.desmos.com/calculator/qxgcyxc3dc
            const O: f64 = 1.5;
            const F: f64 = 0.6;
            0.5 * (1.0 - (((err / best_err) - O) / F).tanh())
        };

        let mut weights: Vec<(f64, f64, f64)> = Vec::with_capacity(self.data.samples.len());
        let mut sw = Vec::new();
        let mut ew = Vec::new();
        let mut daw = Vec::new();
        // let mut errors = Vec::<(Option<f64>, Option<f64>, Option<f64>)>::new();
        for s in &self.data.samples {
            let (eb, ee, eda) = Self::evaluate_traj(s, traj);
            let b1 = eb.map(f64::sqrt).map_or(0.0, get_weight);
            let b2 = ee.map(f64::sqrt).map_or(0.0, get_weight);
            let b3 = eda.map(f64::sqrt).map_or(0.0, get_weight);
            weights.push((b1, b2, b3));
            // errors.push((eb, ee, eda));
            eb.map(|_| sw.push(b1));
            ee.map(|_| ew.push(b2));
            eda.map(|_| daw.push(b3));
        }

        if !sw.is_empty() {
            println!("Start ({}):", sw.len());
            histogram::draw_hitogram(&sw, 5);
        }
        if !ew.is_empty() {
            println!("End ({}):", ew.len());
            histogram::draw_hitogram(&ew, 5);
        }
        if !daw.is_empty() {
            println!("DA ({}):", daw.len());
            histogram::draw_hitogram(&daw, 5);
        }

        let eval_ls = |traj: Line| -> f64 {
            let mut sum = 0.0;
            let mut cnt = 0.0;
            let mut push = |x, w| {
                sum += x * w;
                cnt += w;
            };
            for (s, &(w1, w2, w3)) in self.data.iter().zip(&weights) {
                let (eb, ee, eda) = Self::evaluate_traj(s, traj);
                eb.map(|x| push(x, w1));
                ee.map(|x| push(x, w2));
                eda.map(|x| push(x, w3));
            }
            sum / cnt
        };

        let mut best_err = eval_ls(traj);
        let dif = 10_000.0;
        for _ in 0..2_000 {
            let distr = rand_distr::Normal::new(0.0, dif).unwrap();
            for _ in 0..50 {
                let mut t = traj;
                t.point.x += rng.sample(distr);
                t.point.y += rng.sample(distr);
                t.point.z += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }

            let distr = rand_distr::Normal::new(0.0, 500.0).unwrap();
            for _ in 0..50 {
                let mut t = traj;
                t.direction.x += rng.sample(distr);
                t.direction.y += rng.sample(distr);
                t.direction.z += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
        }

        let mut best_err = eval_ls(traj);
        let dif = 1_000.0;
        for _ in 0..2_000 {
            let distr = rand_distr::Normal::new(0.0, dif).unwrap();
            for _ in 0..50 {
                let mut t = traj;
                t.point.x += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
            for _ in 0..50 {
                let mut t = traj;
                t.point.y += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
            for _ in 0..50 {
                let mut t = traj;
                t.point.z += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }

            let distr = rand_distr::Normal::new(0.0, 100.0).unwrap();
            for _ in 0..50 {
                let mut t = traj;
                t.direction.x += rng.sample(distr);
                t.direction.y += rng.sample(distr);
                t.direction.z += rng.sample(distr);
                let e = eval_ls(t);
                if e < best_err {
                    best_err = e;
                    traj = t;
                }
            }
        }

        // for s in &self.data.samples {
        //     let (eb, ee, eda) = Solver::evaluate_traj(s, point, vel);
        //     eprintln!(
        //         "{} {} {}",
        //         eb.unwrap_or(-1.0),
        //         ee.unwrap_or(-1.0),
        //         eda.unwrap_or(-1.0)
        //     );
        // }

        self.data.compare(traj, "Final answer");

        // if let Some(answer) = self.data.answer {
        //     use scad_gen::constructors::*;
        //     use scad_gen::scope;
        //     let mut doc = scad_gen::Document::default();
        //
        //     let scale = 1e5;
        //     macro_rules! cube {
        //         () => {
        //             cube(Vec3::new(1., 1., 1.) / 2.0)
        //         };
        //     }
        //     let ans = traj.point;
        //     let real_vel = traj.direction;
        //     let vel = traj.direction * 30.;
        //
        //     doc.color(
        //         "red",
        //         Some(0.5),
        //         scope!(ctx, {
        //             ctx.hull(|ctx| {
        //                 ctx.translate(ans / scale - vel * 10.0, cube!());
        //                 ctx.translate(ans / scale, cube!());
        //             })
        //         }),
        //     );
        //
        //     doc.color(
        //         "green",
        //         None,
        //         scope!(ctx, {
        //             ctx.hull(|ctx| {
        //                 ctx.translate(answer.0.point / scale - answer.0.direction * 10.0, cube!());
        //                 ctx.translate(answer.0.point / scale, cube!());
        //             })
        //         }),
        //     );
        //
        //     doc.color("blue", Some(0.25), sphere((EARTH_R / scale) as f32));
        //
        //     {
        //         for s in &self.data.samples {
        //             if !s.observation_matches(point, real_vel) {
        //                 continue;
        //             }
        //             if let Some(k_end) = s.k_end {
        //                 let l_end = lambda(s.location, k_end, point, real_vel);
        //                 let p = point + real_vel * l_end;
        //                 doc.translate(p / scale, cube!());
        //             }
        //         }
        //     }
        //
        //     use std::fs::File;
        //     use std::io::prelude::*;
        //     let mut file = File::create("visual.scad").unwrap();
        //     write!(file, "{doc}").unwrap();
        // }

        let diff = |traj: Line| {
            const D_P: f64 = 50.0;
            const D_V: f64 = 50.0;
            let vel = traj.direction;
            let point = traj.point;
            let evall_ls = |direction, point| eval_ls(Line { point, direction });
            let err0 = evall_ls(vel, point).sqrt();
            let diff_x = (evall_ls(vel, point + Vec3::x_axis() * D_P).sqrt() - err0) / D_P;
            let diff_y = (evall_ls(vel, point + Vec3::y_axis() * D_P).sqrt() - err0) / D_P;
            let diff_z = (evall_ls(vel, point + Vec3::z_axis() * D_P).sqrt() - err0) / D_P;
            let diff_vx = (evall_ls(vel + Vec3::x_axis() * D_V, point).sqrt() - err0) / D_V;
            let diff_vy = (evall_ls(vel + Vec3::y_axis() * D_V, point).sqrt() - err0) / D_V;
            let diff_vz = (evall_ls(vel + Vec3::z_axis() * D_V, point).sqrt() - err0) / D_V;
            (err0, diff_x, diff_y, diff_z, diff_vx, diff_vy, diff_vz)
        };

        let before = diff(traj);
        dbg!(before);

        // loop {
        // for _ in 0..10_000 {
        //     let (_, x, y, z, vx, vy, vz) = diff(traj);
        //     traj.point.x -= x * 2e6;
        //     traj.point.y -= y * 2e6;
        //     traj.point.z -= z * 2e6;
        //     traj.direction.x -= vx * 2e6;
        //     traj.direction.y -= vy * 2e6;
        //     traj.direction.z -= vz * 2e6;
        //     // if x.abs() < 5e-8 && y.abs() < 5e-8 && z.abs() < 5e-8 {
        //     //     break;
        //     // }
        // }

        self.data.compare(traj, "After gradient descent");

        // let after = diff(traj);
        // dbg!(after);

        // eprintln!("delta_x = {} km", after.0 / 6.0 / after.1 * 1e-3);
        // eprintln!("delta_y = {} km", after.0 / 6.0 / after.2 * 1e-3);
        // eprintln!("delta_z = {} km", after.0 / 6.0 / after.3 * 1e-3);
        // eprintln!("delta_v_x = {} km/s", after.0 / 6.0 / after.4 * 1e-3);
        // eprintln!("delta_v_y = {} km/s", after.0 / 6.0 / after.5 * 1e-3);
        // eprintln!("delta_v_z = {} km/s", after.0 / 6.0 / after.6 * 1e-3);

        let (flash, _speed) = self.calc_flash_and_speed(traj);

        let mut points = Vec::with_capacity(1_000);
        for i in 0..1_000 {
            let dx = i as f64 / 1_000.0 - 0.5;
            let dx = dx * 20_000.0;

            let mut t = traj;
            t.point.x += dx;
            points.push((dx, eval_ls(t), false, 1.0));
        }

        draw_plot(&format!("plots/plot-{:.0}.png", traj.point.x), &points).unwrap();

        Solution {
            flash,
            velocity: traj.direction,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, traj: Line) -> (Vec3, f64) {
        let point = traj.point;
        let v = traj.direction.normalized();

        let mut lambdas = Vec::new();
        let mut speeds = Vec::new();

        for s in self.data.samples.iter() {
            if !s.observation_matches(point, v) {
                continue;
            }

            if let Some(k_end) = s.k_end {
                let l_end = lambda(s.location, k_end, point, v);
                lambdas.push(l_end);

                if let Some((duration, k_start)) = s.dur.zip(s.k_start) {
                    let l_start = lambda(s.location, k_start, point, v);
                    let dist = l_end - l_start;
                    debug_assert!(dist > 0.0);
                    speeds.push(dist / duration);
                }
            }
        }

        let l_end = quick_median(&mut lambdas);
        let speed = quick_median(&mut speeds);
        (point + v * l_end, speed)
    }

    /// Calculate the mean of squared errors (less is better)
    pub fn evaluate_traj(
        sample: &DataSample,
        traj: Line,
    ) -> (Option<f64>, Option<f64>, Option<f64>) {
        let p = traj.point;
        let vel = traj.direction;

        let plane = (p - sample.location).cross(vel).normalized();
        let perpendic = vel.cross(plane);

        let es = sample.k_start.map(|k| {
            if perpendic.dot(k) > 0. {
                match sample.dur {
                    Some(dur) => (p - vel * dur - sample.location)
                        .normalized()
                        .dot(k)
                        .acos()
                        .powi(2),
                    None => plane.dot(k).asin().powi(2),
                }
            } else {
                FRAC_PI_2 * FRAC_PI_2
            }
        });

        let ee = sample.k_end.map(|k| {
            if perpendic.dot(k) > 0. {
                (p - sample.location).normalized().dot(k).acos().powi(2)
            } else {
                FRAC_PI_2 * FRAC_PI_2
            }
        });

        let eda = sample.da.map(|da| {
            let k1 = (p - sample.location).normalized();
            let k2 = (p - vel * sample.dur.unwrap_or(0.0) - sample.location).normalized();
            let k = (k1 + k2).normalized();
            angle_diff(descent_angle(sample.location, k, vel.normalized()), da).powi(2)
        });

        let scale_err = |k| move |e| if e < 1e-7 { 1e-7 * k } else { e * k };
        (
            es.map(scale_err(4.0)),
            ee.map(scale_err(4.0)),
            eda.map(scale_err(1.0)),
        )
    }
}

pub fn draw_plot(
    file_name: &str,
    points: &[(f64, f64, bool, f64)],
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let x_min = points
        .iter()
        .map(|x| x.0)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let x_max = points
        .iter()
        .map(|x| x.0)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let y_min = points
        .iter()
        .map(|x| x.1)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let y_max = points
        .iter()
        .map(|x| x.1)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    let root = BitMapBackend::new(file_name, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_by_breakpoints([944], [80]);

    let mut scatter_ctx = ChartBuilder::on(&areas[2])
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    scatter_ctx
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;
    scatter_ctx.draw_series(points.iter().map(|(x, y, is_red, size)| {
        Circle::new(
            (*x, *y),
            *size,
            if *is_red {
                RED.filled()
            } else {
                GREEN.filled()
            },
        )
    }))?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", file_name);

    Ok(())
}
