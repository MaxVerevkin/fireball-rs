use chrono::{DateTime, Utc};
use common::{
    plot::plotters::style::*,
    structs::{Geodetic, Vec3},
};

fn main() {
    let kml = include_str!("../20220516T204421_UT_trajectory.kml");
    let doc = roxmltree::Document::parse(kml).unwrap();
    let mut placemarks: Vec<_> = doc
        .root()
        .children()
        .find(|x| x.tag_name().name() == "kml")
        .unwrap()
        .children()
        .find(|x| x.tag_name().name() == "Document")
        .unwrap()
        .children()
        .filter(|x| x.tag_name().name() == "Placemark")
        .map(parse_placemark)
        .collect();
    placemarks.sort_unstable_by_key(|p| p.timestamp);

    let time_anchor = placemarks.first().unwrap().timestamp;

    let alt_time: Vec<_> = placemarks
        .iter()
        .map(|p| {
            (
                get_dt_secs(time_anchor, p.timestamp),
                p.geo.h * 0.001,
                BLACK,
                1.0,
            )
        })
        .collect();

    common::plot::draw_plot_svg_with_named_axes(
        "plots/alt-time.svg",
        &alt_time,
        &[],
        &[],
        &[],
        "time (s)",
        "altitude (km)",
    )
    .unwrap();

    let vertical_speed_time: Vec<_> = placemarks
        .iter()
        .zip(placemarks.iter().skip(20))
        .map(|(p1, p2)| {
            let dh = p2.geo.h - p1.geo.h;
            let dt = get_dt_secs(p1.timestamp, p2.timestamp);
            (get_dt_secs(time_anchor, p1.timestamp), dh / dt, BLACK, 1.0)
        })
        .collect();

    common::plot::draw_plot_svg_with_named_axes(
        "plots/vertical-speed-time.svg",
        &vertical_speed_time,
        &[],
        &[],
        &[],
        "time (s)",
        "vertical speed (m/s)",
    )
    .unwrap();

    let speed_time: Vec<_> = placemarks
        .iter()
        .zip(placemarks.iter().skip(20))
        .map(|(p1, p2)| {
            let dp = (p2.point - p1.point).norm();
            let dt = get_dt_secs(p1.timestamp, p2.timestamp);
            (get_dt_secs(time_anchor, p1.timestamp), dp / dt, BLACK, 1.0)
        })
        .collect();

    common::plot::draw_plot_svg_with_named_axes(
        "plots/speed-time.svg",
        &speed_time,
        &[],
        &[],
        &[],
        "time (s)",
        "speed (m/s)",
    )
    .unwrap();

    let velocity =
        (placemarks.last().unwrap().point - placemarks.first().unwrap().point).normalize() * 22.0;
    println!(
        "vx = {:.1} # km/s\nvy = {:.1}\nvz = {:.1}",
        velocity.x, velocity.y, velocity.z
    );

    let last_point = placemarks.last().unwrap().geo;
    dbg!(last_point.lat.to_degrees());
    dbg!(last_point.lon.to_degrees());
    dbg!(last_point.h);

    // for n in [2, 5, 10, 20, 100] {
    //     let enter_dir = (placemarks[n - 1].point - placemarks[0].point).normalize();
    //     let final_dir = (placemarks[placemarks.len() - 1].point
    //         - placemarks[placemarks.len() - n].point)
    //         .normalize();
    //     // dbg!(enter_dir);
    //     // dbg!(final_dir);
    //     // dbg!(enter_dir.dot(final_dir));
    //     let da = enter_dir.dot(final_dir).acos().to_degrees();
    //     println!("n={n} -> da={da}");
    // }
}

#[derive(Debug)]
struct Placemark {
    timestamp: DateTime<Utc>,
    point: Vec3,
    geo: Geodetic,
}

fn get_dt_secs(t1: DateTime<Utc>, t2: DateTime<Utc>) -> f64 {
    t2.signed_duration_since(t1).num_milliseconds() as f64 * 0.001
}

fn parse_placemark<'a>(node: roxmltree::Node<'a, '_>) -> Placemark {
    let timestamp_str = node
        .children()
        .find(|x| x.tag_name().name() == "TimeStamp")
        .unwrap()
        .children()
        .find(|x| x.tag_name().name() == "when")
        .unwrap()
        .text()
        .unwrap();
    let timestamp = format!("{timestamp_str}Z").parse().unwrap();

    let mut point_elements = node
        .children()
        .find(|x| x.tag_name().name() == "Point")
        .unwrap()
        .children()
        .find(|x| x.tag_name().name() == "coordinates")
        .unwrap()
        .text()
        .unwrap()
        .split(',');
    let lon = point_elements
        .next()
        .unwrap()
        .parse::<f64>()
        .unwrap()
        .to_radians();
    let lat = point_elements
        .next()
        .unwrap()
        .parse::<f64>()
        .unwrap()
        .to_radians();
    let h = point_elements.next().unwrap().parse().unwrap();

    let geo = Geodetic { lat, lon, h };
    Placemark {
        timestamp,
        point: geo.into_geocentric_cartesian(),
        geo,
    }
}
