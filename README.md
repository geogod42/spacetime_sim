# Spacetime Simulation in Rust

A simple 3D spacetime simulation using [Kiss3D](https://github.com/sebcrozet/kiss3d) to visualize metric grid deformations caused by a moving mass.

## Features

- 3D grid initialized with the Minkowski metric, evolving over time.
- Mass moves in a circular orbit, affecting spacetime curvature.
- Real-time rendering of grid deformations and mass.

## Usage

1. Install [Rust](https://www.rust-lang.org/tools/install).
2. Install dependencies:
   ```sh
   cargo install kiss3d
   ```
3. Run the simulation:
   ```sh
   cargo run
   ```

### Controls

- `.`: Increase mass speed
- `,`: Decrease mass speed

#### License

GNU GPL-3.0