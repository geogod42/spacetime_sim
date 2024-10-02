use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::nalgebra::{Vector3, Point3, Translation3, Matrix4};
use kiss3d::window::Window;
use kiss3d::camera::FirstPerson;
use kiss3d::light::Light;
use std::f32::consts::PI;
use kiss3d::camera::Camera;


/// Constants for the simulation
const GRID_SIZE: usize = 20; // Grid size (adjust for performance)
const DX: f32 = 1.0;         // Spatial grid spacing
const DT: f32 = 0.01;        // Time step for evolution
const STEPS_PER_FRAME: usize = 1; // Number of simulation steps per frame
const MASS_VALUE: f32 = 1.0e5;    // Mass value

/// Initialize the 3D grid for the metric tensor
fn initialize_metric_grid() -> Vec<Vec<Vec<Matrix4<f32>>>> {
    let mut grid = vec![
        vec![
            vec![Matrix4::identity(); GRID_SIZE];
            GRID_SIZE
        ];
        GRID_SIZE
    ];

    // Initialize the metric tensor to Minkowski metric at each grid point
    for i in 0..GRID_SIZE {
        for j in 0..GRID_SIZE {
            for k in 0..GRID_SIZE {
                grid[i][j][k] = minkowski_metric();
            }
        }
    }

    grid
}

/// Function to return the Minkowski metric
fn minkowski_metric() -> Matrix4<f32> {
    Matrix4::new(
        -1.0, 0.0, 0.0, 0.0, // Time component
         0.0, 1.0, 0.0, 0.0, // X component
         0.0, 0.0, 1.0, 0.0, // Y component
         0.0, 0.0, 0.0, 1.0, // Z component
    )
}

/// Initialize the stress-energy tensor
fn initialize_stress_energy_tensor() -> Vec<Vec<Vec<f32>>> {
    vec![
        vec![
            vec![0.0; GRID_SIZE];
            GRID_SIZE
        ];
        GRID_SIZE
    ]
}

/// Calculate the Ricci scalar at a grid point (simplified)
fn calculate_ricci_scalar(
    grid: &Vec<Vec<Vec<Matrix4<f32>>>>,
    i: usize,
    j: usize,
    k: usize,
) -> f32 {
    let mut ricci_scalar = 0.0;

    // Ensure we are not at the boundary to prevent index out of bounds
    if i > 0 && i < GRID_SIZE - 1 && j > 0 && j < GRID_SIZE - 1 && k > 0 && k < GRID_SIZE - 1 {
        for mu in 0..4 {
            // Second derivatives in x, y, z directions using central differences
            let d2g_dx2 = (grid[i + 1][j][k][(mu, mu)] - 2.0 * grid[i][j][k][(mu, mu)] + grid[i - 1][j][k][(mu, mu)]) / (DX * DX);
            let d2g_dy2 = (grid[i][j + 1][k][(mu, mu)] - 2.0 * grid[i][j][k][(mu, mu)] + grid[i][j - 1][k][(mu, mu)]) / (DX * DX);
            let d2g_dz2 = (grid[i][j][k + 1][(mu, mu)] - 2.0 * grid[i][j][k][(mu, mu)] + grid[i][j][k - 1][(mu, mu)]) / (DX * DX);
            ricci_scalar += d2g_dx2 + d2g_dy2 + d2g_dz2;
        }
    }

    ricci_scalar
}

/// Evolve the metric tensor at each grid point using a simplified Einstein equation
fn evolve_metric(
    grid: &mut Vec<Vec<Vec<Matrix4<f32>>>>,
    stress_energy: &Vec<Vec<Vec<f32>>>,
    dt: f32,
) {
    let old_grid = grid.clone();

    for i in 0..GRID_SIZE {
        for j in 0..GRID_SIZE {
            for k in 0..GRID_SIZE {
                let ricci_scalar = calculate_ricci_scalar(&old_grid, i, j, k);
                let t = stress_energy[i][j][k];
                let delta_g = dt * (ricci_scalar - 8.0 * PI * t);

                for mu in 0..4 {
                    grid[i][j][k][(mu, mu)] += delta_g;
                }
            }
        }
    }
}

/// Update the stress-energy tensor with the new mass position
fn update_stress_energy_tensor(
    stress_energy: &mut Vec<Vec<Vec<f32>>>,
    mass_pos: (f32, f32),
) {
    // Clear previous mass position
    for i in 0..GRID_SIZE {
        for j in 0..GRID_SIZE {
            for k in 0..GRID_SIZE {
                stress_energy[i][j][k] = 0.0;
            }
        }
    }

    // Place the mass at the new position
    let (x, y) = mass_pos;
    let ix = ((x + GRID_SIZE as f32 * DX / 2.0) / DX).round() as usize;
    let iy = ((y + GRID_SIZE as f32 * DX / 2.0) / DX).round() as usize;
    let center_z = GRID_SIZE / 2;

    if ix < GRID_SIZE && iy < GRID_SIZE {
        stress_energy[ix][iy][center_z] = MASS_VALUE;
    }
}

/// Compute the new position of the mass in a 2D circular orbit
fn update_mass_position(time: f32, orbit_radius: f32, orbit_speed: f32) -> (f32, f32) {
    let x = orbit_radius * (orbit_speed * time).cos();
    let y = orbit_radius * (orbit_speed * time).sin();
    (x, y)
}

/// Render the metric grid using Kiss3D
fn render_metric(
    window: &mut Window,
    grid: &Vec<Vec<Vec<Matrix4<f32>>>>,
    grid_center: f32,
    scaling_factor: f32,
) {
    // Clear previous grid lines
    window.scene_mut().unlink();

    // Compute displaced positions based on curvature
    let mut positions = vec![vec![vec![Vector3::zeros(); GRID_SIZE]; GRID_SIZE]; GRID_SIZE];

    for i in 0..GRID_SIZE {
        for j in 0..GRID_SIZE {
            for k in 0..GRID_SIZE {
                let x = i as f32 * DX as f32 - grid_center;
                let y = j as f32 * DX as f32 - grid_center;
                let z = k as f32 * DX as f32 - grid_center;

                let position = Vector3::new(x, y, z);
                let curvature = (grid[i][j][k] - minkowski_metric()).norm();
                let displacement_magnitude = curvature * scaling_factor;

                let mass_center = Vector3::new(0.0, 0.0, 0.0); // Assuming mass at origin

                // Avoid division by zero
                let direction = if (position - mass_center).magnitude() > 0.0 {
                    (position - mass_center).normalize()
                } else {
                    Vector3::zeros()
                };

                let displacement = direction * displacement_magnitude;
                let new_position = position + displacement;

                positions[i][j][k] = new_position;
            }
        }
    }

    // Draw lines between adjacent grid points to form the grid
    for i in 0..GRID_SIZE {
        for j in 0..GRID_SIZE {
            for k in 0..GRID_SIZE {
                let p = Point3::from(positions[i][j][k]);

                // Draw lines in x-direction
                if i < GRID_SIZE - 1 {
                    let p_x = Point3::from(positions[i + 1][j][k]);
                    window.draw_line(&p, &p_x, &Point3::new(1.0, 1.0, 1.0));
                }

                // Draw lines in y-direction
                if j < GRID_SIZE - 1 {
                    let p_y = Point3::from(positions[i][j + 1][k]);
                    window.draw_line(&p, &p_y, &Point3::new(1.0, 1.0, 1.0));
                }

                // Draw lines in z-direction
                if k < GRID_SIZE - 1 {
                    let p_z = Point3::from(positions[i][j][k + 1]);
                    window.draw_line(&p, &p_z, &Point3::new(1.0, 1.0, 1.0));
                }
            }
        }
    }
}

fn main() {
    let mut window = Window::new("spacetime_sim");
    window.set_light(Light::StickToCamera);
    window.set_background_color(0.0, 0.0, 0.0);

    // Initialize the metric grid and stress-energy tensor
    let mut metric_grid = initialize_metric_grid();
    let mut stress_energy_tensor = initialize_stress_energy_tensor();

    // Initialize the camera with custom clipping planes
    let mut camera = FirstPerson::new_with_frustrum(
        45.0,                    // Field of view
        0.1,                     // Near clipping plane
        10000.0,                 // Far clipping plane
        Point3::new(50.0, 50.0, 50.0), // Eye position
        Point3::origin(),         // Look-at point
    );

    // Add the moving mass as a single sphere
    let mut mass_sphere = window.add_sphere(4.0);
    mass_sphere.set_color(0.0, 0.0, 0.0); 

    let mut time = 0.0;
    let mut mass_speed = 1.0;  // Initial speed for the mass orbit

    // Precompute the center of the grid for displacement calculations
    let grid_center = (GRID_SIZE as f32 * DX as f32) / 2.0;
    let scaling_factor = 1e-2; // Scaling factor for displacement based on curvature

    while window.render_with_camera(&mut camera) {
        // Event handling for key presses
        for event in window.events().iter() {
            match event.value {
                WindowEvent::Key(Key::Period, Action::Press, _) => {
                    mass_speed += 0.25; // Increase speed when `.` key is pressed
                }
                WindowEvent::Key(Key::Comma, Action::Press, _) => {
                    mass_speed -= 0.25; // Decrease speed when `,` key is pressed
                }
                _ => {}
            }
        }

        // Compute the new position of the mass with the updated speed
        let mass_pos = update_mass_position(time, GRID_SIZE as f32 * DX / 4.0, mass_speed);

        // Update the stress-energy tensor with the new mass position
        update_stress_energy_tensor(&mut stress_energy_tensor, mass_pos);

        // Evolve the metric grid based on the updated stress-energy tensor
        for _ in 0..STEPS_PER_FRAME {
            evolve_metric(&mut metric_grid, &stress_energy_tensor, DT);
            time += DT;
        }

        // Calculate the distance from the camera's current position (eye) to the origin
        let eye_position = camera.eye();

        // Set movement step to increase with the square root of the distance
        camera.set_move_step((eye_position.coords.magnitude()).sqrt() * 0.05);


        // Render the deformed grid
        render_metric(&mut window, &metric_grid, grid_center, scaling_factor);

        // Update the position of the mass sphere
        let mass_x = mass_pos.0;
        let mass_y = mass_pos.1;
        let mass_z = 0.0; // Assuming 2D orbit in the XY plane

        // Update the sphere's translation using Translation3
        mass_sphere.set_local_translation(Translation3::new(mass_x, mass_y, mass_z));
    }
}
