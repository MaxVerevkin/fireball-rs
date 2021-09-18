//! The implementation

// use crate::constants::EARTH_R;
use crate::data::{Data, DataSample};
use crate::maths::*;
use crate::structs::*;

/// Contains all necessary information to solve the problem
#[derive(Clone)]
pub struct Solver {
    data: Data,
    params: Params,
}

/// Contains all necessary information to solve the problem
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub threads: usize,
}

/// The answer
#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub flash: Vec3,
    pub velocity: Vec3,
}

fn solve_pair(s1: &DataSample, s2: &DataSample) -> Option<(Vec3, Vec3)> {
    let n1 = s1.global_start.cross(s1.global_end).normalized();
    let n2 = s2.global_start.cross(s2.global_end).normalized();

    let v = n1.cross(n2).normalized();
    if !v.is_normal() {
        // Something went wrong
        return None;
    }

    // Check velocity direction
    let vp1 = s1.global_end - s1.global_start;
    let vp2 = s2.global_end - s2.global_start;
    let v = match (v.dot(vp1) > 0., v.dot(vp2) > 0.) {
        // Observers are giving opposite directions
        (true, false) | (false, true) => return None,
        // Both observers argee that velocity should be negated
        (false, false) => -v,
        // Both observers argee that velocity should not be changed
        _ => v,
    };

    if s1.pos.dot(v) > 0. {
        return None;
    }

    let l1_end = (s2.pos - s1.pos).dot(n2) / s1.global_end.dot(n2);
    let l2_end = (s1.pos - s2.pos).dot(n1) / s2.global_end.dot(n1);
    let l1_begin = (s2.pos - s1.pos).dot(n2) / s1.global_start.dot(n2);
    let l2_begin = (s1.pos - s2.pos).dot(n1) / s2.global_start.dot(n1);

    if l1_end < 0. || l2_end < 0. || l1_begin < 0. || l2_begin < 0. {
        // Trajectory is below the ground
        return None;
    }

    let p1 = s1.pos + s1.global_end * l1_end;
    let p2 = s2.pos + s2.global_end * l2_end;
    let p = (p1 + p2) * 0.5;

    Some((p, v))
}

impl Solver {
    /// Construct new solver given dataset and paramesters
    pub fn new(data: Data, params: Params) -> Solver {
        Solver { data, params }
    }

    /// Find the solution
    pub fn solve(&mut self) -> Solution {
        // let set = [
        //     1, 21, 61, 76, 91, 112, 118, 133, 136, 139, 173, 174, 195, 202,
        // ];

        // let mut new_samples = Vec::with_capacity(set.len());
        // for &i in &set {
        // new_samples.push(self.data.samples[i].clone());
        // }

        let mut vec = Vec::new();
        let mut point = Vec3::default();
        let mut vel = Vec3::default();

        for i in 0..(self.data.samples.len() - 1) {
            // for i in 0..(set.len() - 1) {
            let s1 = &self.data.samples[i];
            // let s1 = &new_samples[i];
            for j in (i + 1)..self.data.samples.len() {
                // for j in (i + 1)..set.len() {
                let s2 = &self.data.samples[j];
                // let s2 = &new_samples[j];
                if let Some((p, v)) = solve_pair(s1, s2) {
                    // println!("{} {}", i, j);
                    point += p;
                    vel += v;
                    vec.push((p, v));
                }
            }
        }

        if false {
            //             use crate::scad;
            //
            //             impl scad::ScadValue for Vec3 {
            //                 fn scad_str(&self) -> String {
            //                     format!("[{},{},{}]", self.x, self.y, self.z)
            //                 }
            //             }
            //
            //             let mut doc = scad::Document::default();
            //             doc.var("scale", 1e5);
            //             doc.var("data", &*vec);
            //             doc.for_loop("line = data", |ctx| {
            //                 ctx.hull(|ctx| {
            //                     ctx.translate("line[0] / scale", |ctx| {
            //                         ctx.var("hello", 1);
            //                     });
            //                     ctx.var("bye", 3);
            //                 })
            //             });
            //
            //             dbg!(&doc);
            //             eprintln!("\n# DOC START");
            //             eprint!("{}", doc);
            //             eprintln!("# DOC END\n");

            use std::fs::File;
            use std::io::prelude::*;
            let mut file = File::create("visual.scad").unwrap();

            let scale = 1e5;
            writeln!(file, "scale = {};", scale).unwrap();
            writeln!(file, "data = [").unwrap();
            for (p, v) in &vec {
                writeln!(
                    file,
                    "  [[{},{},{}],[{},{},{}]],",
                    p.x, p.y, p.z, v.x, v.y, v.z
                )
                .unwrap();
            }
            writeln!(file, "];").unwrap();

            writeln!(file, "CUBE = [0.1, 0.1, 0.1];").unwrap();
            writeln!(file, "for (line = data) {{").unwrap();
            writeln!(file, "    hull() {{").unwrap();
            writeln!(
                file,
                "        translate(line[0] / scale - line[1] * 10) cube(CUBE / 2);"
            )
            .unwrap();
            writeln!(
                file,
                "        translate(line[0] / scale + line[1] * 0) cube(CUBE / 2);"
            )
            .unwrap();
            writeln!(file, "    }}").unwrap();
            writeln!(file, "}}").unwrap();

            let true_ans: Vec3 = Spherical {
                lat: 33.1f64.to_radians(),
                lon: 34.3f64.to_radians(),
                r: crate::constants::EARTH_R + 43_300.,
            }
            .into();
            writeln!(
                file,
                "ANS = [{},{},{}]; VEL = [-7.5,-23.5,-11.9];",
                true_ans.x, true_ans.y, true_ans.z
            )
            .unwrap();
            writeln!(file, "color(\"red\") hull() {{").unwrap();
            writeln!(
                file,
                "        translate(ANS / scale - VEL * 10) cube(CUBE / 2);"
            )
            .unwrap();
            writeln!(
                file,
                "        translate(ANS / scale - VEL * 0) cube(CUBE / 2);"
            )
            .unwrap();
            writeln!(file, "}}").unwrap();

            writeln!(
                file,
                "color(\"blue\", 0.25) sphere(r = {} / scale);",
                crate::constants::EARTH_R
            )
            .unwrap();

            writeln!(file, "$vpt = data[0][0] / scale;").unwrap();
            writeln!(file, "$vpd = 150;").unwrap();
        }

        point /= vec.len() as f64;
        vel.normalize();

        //         vec.sort_unstable_by(|a, b| a.0.x.partial_cmp(&b.0.x).unwrap());
        //         let x = vec[vec.len() / 2].0.x;
        //         vec.sort_unstable_by(|a, b| a.0.y.partial_cmp(&b.0.y).unwrap());
        //         let y = vec[vec.len() / 2].0.y;
        //         vec.sort_unstable_by(|a, b| a.0.z.partial_cmp(&b.0.z).unwrap());
        //         let z = vec[vec.len() / 2].0.z;
        //         let point = Vec3 { x, y, z };
        //
        //         vec.sort_unstable_by(|a, b| a.1.x.partial_cmp(&b.1.x).unwrap());
        //         let x = vec[vec.len() / 2].1.x;
        //         vec.sort_unstable_by(|a, b| a.1.y.partial_cmp(&b.1.y).unwrap());
        //         let y = vec[vec.len() / 2].1.y;
        //         vec.sort_unstable_by(|a, b| a.1.z.partial_cmp(&b.1.z).unwrap());
        //         let z = vec[vec.len() / 2].1.z;
        //         let vel = Vec3 { x, y, z }.normalized();

        dbg!(vec.len());

        let mut s_point = 0.;
        let mut s_vel = 0.;
        for (p, v) in &vec {
            s_point += (point - *p).len().powi(2);
            s_vel += vel.dot(*v).abs().min(1.).acos().powi(2);
        }
        s_point = (s_point / vec.len() as f64).sqrt();
        s_vel = (s_vel / vec.len() as f64).sqrt();

        dbg!(s_point);
        dbg!(s_vel.to_degrees());

        let (flash, speed) = self.calc_flash_and_speed(point, vel);
        let velocity = vel * speed;

        Solution {
            // error: self.evaluate_traj(point, vel),
            flash,
            velocity,
        }
    }

    /// Calculate the flash location and the speed of the fireball
    fn calc_flash_and_speed(&self, point: Vec3, v: Vec3) -> (Vec3, f64) {
        let mut speed_vec = Vec::new();
        let mut l_end_vec = Vec::new();

        for sample in &self.data.samples {
            let pip = point - sample.pos;
            let plane = pip.cross(v).normalized();
            let perpendic = v.cross(plane).normalized();

            let k_start = sample.global_start;
            let k_end = sample.global_end;

            let trust_start = perpendic.dot(k_start) > 0.;
            let trust_end = perpendic.dot(k_end) > 0.;

            if trust_end {
                let l_end = pip.dot(perpendic * (k_end.dot(v) / k_end.dot(perpendic)) - v);
                l_end_vec.push(l_end);
                if trust_start && sample.duration > 0. {
                    let l_start =
                        pip.dot(perpendic * (k_start.dot(v) / k_start.dot(perpendic)) - v);
                    speed_vec.push((l_end - l_start) / sample.duration);
                }
            }
        }

        let l_end = quick_median(&mut l_end_vec);
        let speed = quick_median(&mut speed_vec);

        //dbg!(l_end_vec.len());
        //dbg!(speed_vec.len());

        (point + v * l_end, speed)
    }

    // /// Calculate the mean of squared errors (less is better)
    //     pub fn evaluate_traj(&self, p: Vec3, vel: Vec3) -> f64 {
    //         let mut error = 0.;
    //         for sample in &self.data.samples {
    //             let plane = (p - sample.global_pos).cross(vel).normalized();
    //             let perpendic = vel.cross(plane);
    //
    //             let k_start = sample.global_start;
    //             let k_end = sample.global_end;
    //
    //             let e_start = plane.dot(k_start).asin();
    //             let e_end = plane.dot(k_end).asin();
    //
    //             if sample.trust_start {
    //                 let c = perpendic.dot(k_start);
    //                 if c > 0. {
    //                     error += e_start * e_start;
    //
    //                     let angle = descent_angle(sample.global_pos, k_start, vel);
    //                     let diff = angle_diff(angle, sample.descent_angle);
    //                     error += diff * diff;
    //                 } else {
    //                     //let c = c.clamp(0., 0.1) * 10.;
    //                     //error += c * e_start * e_start + (1. - c) * FRAC_PI_2 * FRAC_PI_2;
    //                     error += FRAC_PI_2 * FRAC_PI_2;
    //                     error += PI * PI;
    //                 }
    //                 //error += diff * diff;
    //             }
    //             if sample.trust_end {
    //                 let c = perpendic.dot(k_end);
    //                 if c > 0. {
    //                     error += e_end * e_end;
    //                 } else {
    //                     //let c = c.clamp(0., 0.1) * 10.;
    //                     //error += c * e_end * e_end + (1. - c) * FRAC_PI_2 * FRAC_PI_2;
    //                     error += FRAC_PI_2 * FRAC_PI_2;
    //                 }
    //             }
    //         }
    //         error
    //     }
}
