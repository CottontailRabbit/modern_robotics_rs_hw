
use modern_robotics_rs::{
    fkin_space, fkin_body, jacobian_space, jacobian_body,
    matrix_log6, se3_to_vec,
};

use nalgebra::{Matrix4, Vector6, Vector3, Matrix6, OMatrix, U6, Dyn};

// 6×n Jacobian 타입 별칭 (여기서 직접 정의)
pub type Mat6xX = OMatrix<f64, U6, Dyn>;

pub struct IKSolveLog {
    pub iters: Vec<usize>,
    pub ang_norm: Vec<f64>,
    pub lin_norm: Vec<f64>,
    pub thetas: Vec<Vec<f64>>, // optional: 각 반복의 θ 저장(원하면)
}

/// DLS step: Δθ = Jᵀ (J Jᵀ + λ² I)⁻¹ V
fn dls_step(j: &Mat6xX, v: &Vector6<f64>, lambda: f64) -> Vec<f64> {
    use nalgebra::{Matrix6, OMatrix, Dyn, U6};
    let jt = j.transpose();               // n×6
    let jj_t: Matrix6<f64> = j * j.transpose(); // 6×6
    
    // Adaptive damping based on conditioning
    let condition_number = jj_t.norm() / (jj_t.try_inverse().map_or(1e-12, |inv| inv.norm()));
    let adaptive_lambda = if condition_number > 1e6 { lambda * 1000.0 } else { lambda };
    
    let a = jj_t + Matrix6::identity() * (adaptive_lambda * adaptive_lambda);

    let y = if let Some(cho) = a.cholesky() {
        cho.solve(v)
    } else {
        a.lu().solve(v).expect("DLS solve failed")
    };
    let delta = jt * y; // n×1
    
    // Limit step size to prevent large jumps
    let delta_vec: Vec<f64> = delta.iter().copied().collect();
    let max_step = 0.1; // Limit individual joint steps to 0.1 radians
    delta_vec.iter().map(|&x| if x.abs() > max_step { x.signum() * max_step } else { x }).collect()
}

fn twist_norms(v: &Vector6<f64>) -> (f64, f64) {
    let omg = Vector3::new(v[0], v[1], v[2]).norm();
    let lin = Vector3::new(v[3], v[4], v[5]).norm();
    (omg, lin)
}

/// 수렴 기준
fn converged(v: &Vector6<f64>, eomg: f64, ev: f64) -> bool {
    let (o, l) = twist_norms(v);
    o < eomg && l < ev
}

/// IKinSpace with logging
pub fn ikin_space_log(
    m: &Matrix4<f64>,
    s_list: &[Vector6<f64>],
    t_goal: &Matrix4<f64>,
    theta0: &[f64],
    eomg: f64,
    ev: f64,
    max_iter: usize,
) -> (Vec<f64>, bool, IKSolveLog) {
    let mut theta = theta0.to_vec();
    let lambda = 1e-3; // Increased damping for better stability
    let mut log = IKSolveLog { iters: vec![], ang_norm: vec![], lin_norm: vec![], thetas: vec![] };

    for it in 0..max_iter {
        let t_now = fkin_space(m, s_list, &theta);
        let t_err = t_goal * trans_inv(&t_now);
        let v_s = se3_to_vec(&matrix_log6(&t_err));

        let (o, l) = twist_norms(&v_s);
        log.iters.push(it);
        log.ang_norm.push(o);
        log.lin_norm.push(l);
        log.thetas.push(theta.clone());

        if converged(&v_s, eomg, ev) {
            return (theta, true, log);
        }

        let js = jacobian_space(s_list, &theta);
        let dtheta = dls_step(&js, &v_s, lambda);

        // 간단 백트래킹
        let mut accepted = false;
        for alpha in [1.0, 0.5, 0.25, 0.1] {
            let mut th_try = theta.clone();
            for i in 0..th_try.len() { th_try[i] += alpha * dtheta[i]; }
            let t_try = fkin_space(m, s_list, &th_try);
            let v_try = se3_to_vec(&matrix_log6(&(t_goal * trans_inv(&t_try))));
            let (o1, l1) = twist_norms(&v_try);
            if o1 <= o && l1 <= l {
                theta = th_try;
                accepted = true;
                break;
            }
        }
        if !accepted {
            // Reduce step size when backtracking fails
            for i in 0..theta.len() { theta[i] += 0.01 * dtheta[i]; }
        }
    }
    (theta, false, log)
}

/// IKinBody with logging (WAM 용)
pub fn ikin_body_log(
    m: &Matrix4<f64>,
    b_list: &[Vector6<f64>],
    t_goal: &Matrix4<f64>,
    theta0: &[f64],
    eomg: f64,
    ev: f64,
    max_iter: usize,
) -> (Vec<f64>, bool, IKSolveLog) {
    let mut theta = theta0.to_vec();
    let lambda = 1e-3; // Increased damping for better stability
    let mut log = IKSolveLog { iters: vec![], ang_norm: vec![], lin_norm: vec![], thetas: vec![] };

    for it in 0..max_iter {
        let t_now = fkin_body(m, b_list, &theta);
        let t_err = trans_inv(&t_now) * t_goal;
        let v_b = se3_to_vec(&matrix_log6(&t_err));

        let (o, l) = twist_norms(&v_b);
        log.iters.push(it);
        log.ang_norm.push(o);
        log.lin_norm.push(l);
        log.thetas.push(theta.clone());

        if converged(&v_b, eomg, ev) {
            return (theta, true, log);
        }

        let jb = jacobian_body(b_list, &theta);
        let dtheta = dls_step(&jb, &v_b, lambda);

        let mut accepted = false;
        for alpha in [1.0, 0.5, 0.25, 0.1] {
            let mut th_try = theta.clone();
            for i in 0..th_try.len() { th_try[i] += alpha * dtheta[i]; }
            let t_try = fkin_body(m, b_list, &th_try);
            let v_try = se3_to_vec(&matrix_log6(&(trans_inv(&t_try) * t_goal)));
            let (o1, l1) = twist_norms(&v_try);
            if o1 <= o && l1 <= l {
                theta = th_try;
                accepted = true;
                break;
            }
        }
        if !accepted {
            // Reduce step size when backtracking fails
            for i in 0..theta.len() { theta[i] += 0.01 * dtheta[i]; }
        }
    }
    (theta, false, log)
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