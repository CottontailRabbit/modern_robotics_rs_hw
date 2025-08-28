
use nalgebra::Vector3;
use modern_robotics_rs::vec_to_skew3;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let v = Vector3::new(1.0, 2.0, 3.0);
    let skew = vec_to_skew3(&v);
    println!("Skew-symmetric matrix:\n{}", skew);

    let root = BitMapBackend::new("plot.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y = x^2", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..10f32, 0f32..100f32)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
        (0..=10).map(|x| x as f32).map(|x| (x, x * x)),
        &RED,
    ))?;
    root.present()?;
    println!("Saved plot to plot.png");
    Ok(())

