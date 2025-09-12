use nalgebra::{vector, Matrix4, Vector6};
use modern_robotics_rs::{fkin_space, jacobian_space};

#[test]
fn test_fkin_space_basic() {
    // Simple test with identity configuration
    let m = Matrix4::identity();
    let slist = vec![vector![0.0, 0.0, 1.0, 0.0, 0.0, 0.0]];
    let theta = vec![0.0];
    
    let result = fkin_space(&m, &slist, &theta);
    
    // At zero configuration, should return the home configuration
    assert!((result - m).norm() < 1e-10);
}

#[test]
fn test_jacobian_space_basic() {
    // Simple test with identity configuration
    let slist = vec![vector![0.0, 0.0, 1.0, 0.0, 0.0, 0.0]];
    let theta = vec![0.0];
    
    let result = jacobian_space(&slist, &theta);
    
    // Should return a 6x1 matrix
    assert_eq!(result.nrows(), 6);
    assert_eq!(result.ncols(), 1);
}