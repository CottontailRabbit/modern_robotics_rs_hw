fn main() -> anyhow::Result<()> {
    // Let's add some debug output to understand what's happening
    println!("=== Testing UR5 Exercise 6.12 ===");
    exercise_6_12()?;
    println!("\n=== Testing WAM Exercise 6.13 ===");
    exercise_6_13()?;
    println!("Saved plots: ex6_12_convergence.png, ex6_13_convergence.png");
    Ok(())
}