use nalgebra::{matrix, vector, Matrix4, Vector6};

fn main() {
    // Test data
    let m = Matrix4::identity();
    let slist = vec![Vector6::zeros()];
    let tsd = Matrix4::identity();
    let theta0 = vec![0.0];
    
    // Try calling the IK functions
    let result1 = modern_robotics_rs::ik_in_space(&m, &slist, &tsd, &theta0, 1e-3, 1e-3, 100);
    println!("ik_in_space result: {:?}", result1);
    
    let result2 = modern_robotics_rs::ik_in_body(&m, &slist, &tsd, &theta0, 1e-3, 1e-3, 100);
    println!("ik_in_body result: {:?}", result2);
}