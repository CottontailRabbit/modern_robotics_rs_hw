use nalgebra::{Matrix3, Vector3};
use modern_robotics_rs::vec_to_skew3;

#[test]
fn computes_skew() {
    let v = Vector3::new(1.0, 2.0, 3.0);
    let skew = vec_to_skew3(&v);
    let expected = Matrix3::new(
        0.0, -3.0, 2.0,
        3.0, 0.0, -1.0,
       -2.0, 1.0, 0.0,
    );
    assert_eq!(skew, expected);
}
