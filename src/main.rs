use nalgebra::{matrix, vector, Matrix4, Vector6}; 

mod solve_utils;
use solve_utils::{ikin_space_log, ikin_body_log, wrap_to_2pi, IKSolveLog, trans_inv};
use modern_robotics_rs::{fkin_space, matrix_log6, se3_to_vec};
use plotters::prelude::*;
use anyhow::Result;

fn plot_log_png(path: &str, title: &str, log: &IKSolveLog) -> anyhow::Result<()> {
    use plotters::prelude::*;
    let root = BitMapBackend::new(path, (900, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let max_it = *log.iters.last().unwrap_or(&0);
    let max_ang = log.ang_norm.iter().cloned().fold(0./*dummy*/, f64::max).max(1e-9);
    let max_lin = log.lin_norm.iter().cloned().fold(0./*dummy*/, f64::max).max(1e-9);
    let max_y = max_ang.max(max_lin);

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0usize..(max_it.max(1)), 0f64..(max_y*1.1))?;

    chart.configure_mesh()
        .x_desc("iteration")
        .y_desc("norm")
        .draw()?;

    chart.draw_series(LineSeries::new(
        log.iters.iter().zip(log.ang_norm.iter()).map(|(i, y)| (*i, *y)),
        &BLUE,
    ))?.label("‖ω_err‖").legend(|(x, y)| PathElement::new([(x, y), (x+20, y)], &BLUE));

    chart.draw_series(LineSeries::new(
        log.iters.iter().zip(log.lin_norm.iter()).map(|(i, y)| (*i, *y)),
        &RED,
    ))?.label("‖v_err‖").legend(|(x, y)| PathElement::new([(x, y), (x+20, y)], &RED));

    chart.configure_series_labels().border_style(&BLACK).draw()?;
    root.present()?;
    Ok(())
}

fn ur5_space_screws_and_home() -> (Vec<Vector6<f64>>, Matrix4<f64>) {
    let h1 = 0.089_f64;
    let h2 = 0.095_f64;
    let w1 = 0.109_f64;
    let w2 = 0.082_f64;
    let l1 = 0.425_f64;
    let l2 = 0.392_f64;

    // Slist: columns are screw axes in the space frame at home
    // (ω, v) with v = -ω × q for a point q on the axis.
    let s_list = vec![
        vector![ 0.0,  0.0,  1.0,         0.0,          0.0,            0.0],
        vector![ 0.0,  1.0,  0.0,       - h1,          0.0,            0.0],
        vector![ 0.0,  1.0,  0.0,       - h1,          0.0,            l1 ],
        vector![ 0.0,  1.0,  0.0,       - h1,          0.0,        l1 + l2],
        vector![ 0.0,  0.0, -1.0,       - w1,       l1 + l2,            0.0],
        vector![ 0.0,  1.0,  0.0,   (h2 - h1),        0.0,        l1 + l2],
    ];

    let m = matrix![
        -1.0, 0.0, 0.0,  l1 + l2;
         0.0, 0.0, 1.0,  w1 + w2;
         0.0, 1.0, 0.0,  h1 - h2;
         0.0, 0.0, 0.0,  1.0
    ];

    (s_list, m)
}
/// ==== TODO: WAM (Sec. 4.1.3) body screw axes Blist and home M ====
fn wam_body_screws_and_home() -> (Vec<Vector6<f64>>, Matrix4<f64>) {
    // Figure 4.8 캡션: L1=550 mm, L2=300 mm, L3=60 mm, W1=45 mm
    let l1 = 0.550;
    let l2 = 0.300;
    let l3 = 0.060;
    let w1 = 0.045;

    // M: end-effector frame at zero position (body form)
    // M = diag(R=I, p=[0, 0, L1+L2+L3])
    let m = matrix![
        1.0, 0.0, 0.0, 0.0;
        0.0, 1.0, 0.0, 0.0;
        0.0, 0.0, 1.0, l1 + l2 + l3;
        0.0, 0.0, 0.0, 1.0;
    ];

    // Blist: 각 열이 하나의 스크류축 (ω, v) in body frame
    // 표: i, ωi, vi
    let b_list: [Vector6<f64>; 7] = [
        vector![0.0, 0.0, 1.0,          0.0,          0.0, 0.0],            // 1
        vector![0.0, 1.0, 0.0,  (l1+l2+l3),          0.0, 0.0],            // 2
        vector![0.0, 0.0, 1.0,          0.0,          0.0, 0.0],            // 3
        vector![0.0, 1.0, 0.0,      (l2+l3),          0.0, w1],            // 4
        vector![0.0, 0.0, 1.0,          0.0,          0.0, 0.0],            // 5
        vector![0.0, 1.0, 0.0,          l3,           0.0, 0.0],            // 6
        vector![0.0, 0.0, 1.0,          0.0,          0.0, 0.0],            // 7
    ];

    (b_list.to_vec(), m)
}

/// ---- Exercise 6.12: UR5 / IKinSpace ----
fn exercise_6_12() -> anyhow::Result<()> {
    // Tsd (문제에서 주어진 목표)
    let tsd = matrix![
        0.0,  1.0,  0.0, -0.5;
        0.0,  0.0, -1.0,  0.1;
       -1.0,  0.0,  0.0,  0.1;
        0.0,  0.0,  0.0,  1.0;
    ];
    let (slist, m) = ur5_space_screws_and_home();

    // 초기값 θ0 = 0.1 rad (모든 관절)
    let dof = slist.len();
    let theta0 = vec![0.1; dof];
    let (theta, ok, log) = solve_utils::ikin_space_log(&m, &slist, &tsd, &theta0, 1e-3, 1e-4, 1000);
    let mut theta_wrapped = theta.clone();
    solve_utils::wrap_to_2pi(&mut theta_wrapped);

    println!("# Exercise 6.12 (UR5, IKinSpace)");
    println!("success: {}", ok);
    println!("theta (rad): {:?}", theta);
    println!("theta wrapped [0,2π): {:?}", theta_wrapped);

    // 플롯
    plot_log_png("ex6_12_convergence.png", "UR5 IKinSpace convergence", &log)?;
    Ok(())
}

/// ---- Exercise 6.13: WAM / IKinBody ----
fn exercise_6_13() -> anyhow::Result<()> {
    // Tsd (문제에서 주어진 목표)
    let tsd = matrix![
        1.0, 0.0, 0.0, 0.5;
        0.0, 1.0, 0.0, 0.0;
        0.0, 0.0, 1.0, 0.4;
        0.0, 0.0, 0.0, 1.0;
    ];
    let (blist, m) = wam_body_screws_and_home();

    let dof = blist.len();
    let theta0 = vec![0.1; dof];
    let (theta, ok, log) = solve_utils::ikin_body_log(&m, &blist, &tsd, &theta0, 1e-3, 1e-4, 200);

    let mut theta_wrapped = theta.clone();
    wrap_to_2pi(&mut theta_wrapped);

    println!("# Exercise 6.13 (WAM, IKinBody)");
    println!("success: {}", ok);
    println!("theta (rad): {:?}", theta);
    println!("theta wrapped [0,2π): {:?}", theta_wrapped);

    plot_log_png("ex6_13_convergence.png", "WAM IKinBody convergence", &log)?;
    Ok(())
}

// fn main() -> anyhow::Result<()> {
//     // exercise_6_12()?;
//     // exercise_6_13()?;
//     // println!("Saved plots: ex6_12_convergence.png, ex6_13_convergence.png");
//     // Ok(())

    
// }

fn main() {
    // ----- Example Input (MR core.py와 동일) -----
    // Slist.T 의 각 "열"이 하나의 스크류축 (ω, v)
    let slist: Vec<Vector6<f64>> = vec![
        vector![0.0, 0.0,  1.0,  4.0, 0.0,    0.0],
        vector![0.0, 0.0,  0.0,  0.0, 1.0,    0.0],
        vector![0.0, 0.0, -1.0, -6.0, 0.0, -0.1],
    ];

    let m: Matrix4<f64> = matrix![
        -1.0, 0.0,  0.0, 0.0;
         0.0, 1.0,  0.0, 6.0;
         0.0, 0.0, -1.0, 2.0;
         0.0, 0.0,  0.0, 1.0
    ];

    let t_target: Matrix4<f64> = matrix![
        0.0, 1.0,  0.0,     -5.0;
        1.0, 0.0,  0.0,      4.0;
        0.0, 0.0, -1.0,  1.6858;
        0.0, 0.0,  0.0,      1.0
    ];

    let theta0 = vec![1.5, 2.5, 3.0];
    let eomg = 0.01;
    let ev   = 0.001;

    // ----- Solve (space frame) -----
    let (theta, ok, _log) = ikin_space_log(&m, &slist, &t_target, &theta0, eomg, ev, 200);

    println!("success: {}", ok);
    println!("theta   : {:?}", theta);

    // ----- 검증: 최종 포즈 오차 -----
    let t_now = fkin_space(&m, &slist, &theta);
    let v = se3_to_vec(&matrix_log6(&(t_target * trans_inv(&t_now))));
    let omg = nalgebra::Vector3::new(v[0], v[1], v[2]).norm();
    let lin = nalgebra::Vector3::new(v[3], v[4], v[5]).norm();
    println!("final error  |ω|={:.6e}, |v|={:.6e}", omg, lin);

    // ----- 기대값과 비교 (참고) -----
    let expected = [1.57073783, 2.99966384, 3.14153420];
    println!("expected     {:?}", expected);
}


