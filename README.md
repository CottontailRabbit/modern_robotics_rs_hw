# modern_robotics_rs_hw

This project contains Rust exercises for the course "Modern Robotics". It
demonstrates using the
[`modern_robotics_rs`](https://github.com/CottontailRabbit/modern_robotics_rs)
library for computing skew-symmetric matrices and the
[`plotters`](https://crates.io/crates/plotters) crate for simple plots.

## Setup

Install Rust using rustup:

```bash
curl https://sh.rustup.rs -sSf | sh -s -- -y
. "$HOME/.cargo/env"
```

## Build and Test

```bash
cargo build
cargo test
```

## Run

```bash
cargo run
```

Running the binary prints the skew-symmetric matrix computed from a sample vector.
It also generates a `plot.png` file with a quadratic curve.

