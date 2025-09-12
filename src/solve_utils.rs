
use modern_robotics_rs::{
    fk_in_space, fk_in_body, ik_in_space, ik_in_body,
    matrix_log6, se3_to_vec,
};

use nalgebra::{Matrix4, Vector6, Vector3};

pub struct IKSolveLog {
    pub iters: Vec<usize>,
    pub ang_norm: Vec<f64>,
    pub lin_norm: Vec<f64>,
}

fn twist_norms(v: &Vector6<f64>) -> (f64, f64) {
    let omg = Vector3::new(v[0], v[1], v[2]).norm();
    let lin = Vector3::new(v[3], v[4], v[5]).norm();
    (omg, lin)
}

/// IKinSpace using library function with minimal logging for plots
pub fn ikin_space_log(
    m: &Matrix4<f64>,
    s_list: &[Vector6<f64>],
    t_goal: &Matrix4<f64>,
    theta0: &[f64],
    eomg: f64,
    ev: f64,
    _max_iter: usize, // Library handles max iterations internally
) -> (Vec<f64>, bool, IKSolveLog) {
    // Use the library's built-in IK solver
    let (theta, success) = ik_in_space(s_list, m, t_goal, theta0, eomg, ev);
    
    // Create minimal log for plotting - just initial and final states
    let mut log = IKSolveLog { 
        iters: vec![0], 
        ang_norm: vec![], 
        lin_norm: vec![], 
    };
    
    // Calculate initial error
    let t_initial = fk_in_space(m, s_list, theta0);
    let t_err_initial = t_goal * trans_inv(&t_initial);
    let v_s_initial = se3_to_vec(&matrix_log6(&t_err_initial));
    let (o_initial, l_initial) = twist_norms(&v_s_initial);
    log.ang_norm.push(o_initial);
    log.lin_norm.push(l_initial);
    
    // Add final iteration if successful
    if success {
        log.iters.push(1);
        let t_final = fk_in_space(m, s_list, &theta);
        let t_err_final = t_goal * trans_inv(&t_final);
        let v_s_final = se3_to_vec(&matrix_log6(&t_err_final));
        let (o_final, l_final) = twist_norms(&v_s_final);
        log.ang_norm.push(o_final);
        log.lin_norm.push(l_final);
    }
    
    (theta, success, log)
}

/// IKinBody using library function with minimal logging for plots
pub fn ikin_body_log(
    m: &Matrix4<f64>,
    b_list: &[Vector6<f64>],
    t_goal: &Matrix4<f64>,
    theta0: &[f64],
    eomg: f64,
    ev: f64,
    _max_iter: usize, // Library handles max iterations internally
) -> (Vec<f64>, bool, IKSolveLog) {
    // Use the library's built-in IK solver
    let (theta, success) = ik_in_body(b_list, m, t_goal, theta0, eomg, ev);
    
    // Create minimal log for plotting - just initial and final states
    let mut log = IKSolveLog { 
        iters: vec![0], 
        ang_norm: vec![], 
        lin_norm: vec![], 
    };
    
    // Calculate initial error
    let t_initial = fk_in_body(m, b_list, theta0);
    let t_err_initial = trans_inv(&t_initial) * t_goal;
    let v_b_initial = se3_to_vec(&matrix_log6(&t_err_initial));
    let (o_initial, l_initial) = twist_norms(&v_b_initial);
    log.ang_norm.push(o_initial);
    log.lin_norm.push(l_initial);
    
    // Add final iteration if successful
    if success {
        log.iters.push(1);
        let t_final = fk_in_body(m, b_list, &theta);
        let t_err_final = trans_inv(&t_final) * t_goal;
        let v_b_final = se3_to_vec(&matrix_log6(&t_err_final));
        let (o_final, l_final) = twist_norms(&v_b_final);
        log.ang_norm.push(o_final);
        log.lin_norm.push(l_final);
    }
    
    (theta, success, log)
}

/// [0, 2π) wrapping
pub fn wrap_to_2pi(theta: &mut [f64]) {
    for t in theta {
        let mut x = *t % (2.0 * std::f64::consts::PI);
        if x < 0.0 { x += 2.0 * std::f64::consts::PI; }
        *t = x;
    }
}
pub fn trans_inv(t: &Matrix4<f64>) -> Matrix4<f64> {
    let r = t.fixed_view::<3, 3>(0, 0).clone_owned();
    let p = t.fixed_view::<3, 1>(0, 3).clone_owned();
    let r_t = r.transpose();
    let mut t_inv = Matrix4::identity();
    t_inv.fixed_view_mut::<3, 3>(0, 0).copy_from(&r_t);
    t_inv.fixed_view_mut::<3, 1>(0, 3).copy_from(&(-&r_t * p));
    t_inv
}