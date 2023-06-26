use std::env;

use common::{
    obs_data::{Data, RawData},
    plot::{draw_plot_svg_with_named_axes, plotters::style::BLACK, RGBColor},
};

#[allow(dead_code)]
const COLORS: &[RGBColor] = {
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

fn main() {
    let data_sets: Vec<_> = env::args()
        .skip(1)
        .map(|x| RawData::read_from_toml_file(x).finalize())
        .collect();
    if data_sets.is_empty() {
        eprintln!("provide arguments");
        return;
    }

    {
        let points: Vec<_> = data_sets
            .iter()
            .enumerate()
            .flat_map(|(i, data)| get_altitudes(i, data).into_iter())
            .map(Point::to_degrees)
            .map(|p| {
                // let c = COLORS[p.data_set];
                let c = BLACK;
                (p.observed, p.calculated, c, 1.0)
            })
            .collect();

        draw_plot_svg_with_named_axes(
            "plots/altitudes.svg",
            &points,
            &[],
            &[],
            &[(&|x| x, RGBColor(0, 0, 255))],
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
    data_set: usize,
}

impl Point {
    fn to_degrees(self) -> Self {
        Self {
            observed: self.observed.to_degrees(),
            calculated: self.calculated.to_degrees(),
            ..self
        }
    }
}

fn get_altitudes(data_set: usize, data: &Data) -> Vec<Point> {
    let Some(traj) = data.answer.and_then(|a| a.traj) else { return vec![] };

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
                data_set,
            })
        })
        .collect()
}
