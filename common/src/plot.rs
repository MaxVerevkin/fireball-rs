pub use plotters;
pub use plotters::style::RGBColor;

pub fn draw_plot_svg(
    file_name: &str,
    points: &[(f64, f64, RGBColor, f64)],
    vertical_lines: &[(f64, RGBColor)],
    horizontal_lines: &[(f64, RGBColor)],
    fns: &[(&dyn Fn(f64) -> f64, RGBColor)],
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let x_min = points
        .iter()
        .map(|x| x.0)
        .min_by(|a, b| a.total_cmp(b))
        .unwrap();
    let x_max = points
        .iter()
        .map(|x| x.0)
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();
    let y_min = points
        .iter()
        .map(|x| x.1)
        .min_by(|a, b| a.total_cmp(b))
        .unwrap();
    let y_max = points
        .iter()
        .map(|x| x.1)
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();

    let root = SVGBackend::new(file_name, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut scatter_ctx = ChartBuilder::on(&root)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    for &(x, color) in vertical_lines {
        scatter_ctx.draw_series(LineSeries::new([(x, y_min), (x, y_max)], color))?;
    }
    for &(y, color) in horizontal_lines {
        scatter_ctx.draw_series(LineSeries::new([(x_min, y), (x_max, y)], color))?;
    }
    for &(function, color) in fns {
        let mut points = Vec::new();
        let mut cur_x = x_min;
        while cur_x <= x_max {
            points.push((cur_x, function(cur_x)));
            cur_x += 0.01;
        }
        scatter_ctx.draw_series(LineSeries::new(points, color))?;
    }
    scatter_ctx
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;
    scatter_ctx.draw_series(
        points
            .iter()
            .map(|(x, y, color, size)| Circle::new((*x, *y), *size, color)),
    )?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    eprintln!("Result has been saved to {}", file_name);

    Ok(())
}

pub fn weight_to_rgb(w: f64) -> RGBColor {
    let red = 1.0 - w;
    let green = w;
    RGBColor((red * 255.0) as _, (green * 255.0) as _, 0)
}
