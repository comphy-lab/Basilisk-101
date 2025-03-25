/**
 ## 1D Steady Heat Conduction Solver (Basilisk version)
 
 This program solves the heat conduction equation to reach steady state:
 
 $$
 \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 $$
 
 ## Subject to Dirichlet boundary conditions:
 
 - T(0) = 0
 - T(1) = 1
 
 The steady-state analytical solution is the linear profile $T(x) = x$.
*/

#include "grid/cartesian1D.h" // 1D grid
#include "run.h" // utility functions for running the simulation

// Declare scalar field for temperature
scalar T[];
T[left] = dirichlet(0.0); // left boundary
T[right] = dirichlet(1.0); /**
 * @brief Initializes and starts the 1D steady-state heat conduction simulation.
 *
 * Configures the computational domain by setting the domain length (L0 = 1.0), left boundary (X0 = 0.0),
 * and number of cells (N = 200) for spatial discretization. Calculates the time step (DT) using the CFL stability
 * criterion for the explicit scheme, creates an output directory ("intermediate") for simulation results, and 
 * invokes the simulation routine via run().
 */

int main() {
  // Domain setup
  L0 = 1.0;     // Domain length
  X0 = 0.0;    // Left boundary
  N = 200;       // Number of cells
  
  // Set the timestep based on stability criterion (CFL condition)
  // dt = dx^2/2 for explicit scheme
  DT = (L0/N)*(L0/N)/2;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field
 
 Starting with zero temperature everywhere, letting the boundary
 conditions drive the evolution toward steady state.
*/
event init (t = 0) {
  foreach()
    T[] = 0.0; 
}

/**
 ## Time integration using explicit finite volume method
 */
/**
 * @brief Advances the temperature field by one time step.
 *
 * This event implements an explicit finite-difference update for the 1D heat conduction equation.
 * It computes the discrete Laplacian using the values from the neighboring cells and updates the
 * temperature field with the time step DT. The update is applied to every cell in the domain.
 */
event marching (i++) {
  foreach() {
    // Proper heat equation time stepping
    T[] += DT*(T[1] - 2*T[] + T[-1])/(Delta*Delta);
  }
}

/**
 * @brief Monitors the convergence of the temperature field to the steady state.
 *
 * This event computes the normalized error between the simulated temperature field and 
 * the analytical solution T(x) = x by averaging the absolute differences over all grid cells.
 * The error is printed every 1000 iterations for monitoring purposes. When the normalized error 
 * falls below 1e-10, the event prints a convergence message and returns 1 to signal simulation termination.
 *
 * @return int Returns 1 when the convergence criterion is met.
 */
event testingConvergence (i++; i < 100000) {
  double error = 0.0;
  foreach() {
    // Calculate error relative to analytical solution T(x) = x
    error += fabs(T[] - x);
  }
  error = error/N; // Normalize by number of cells
  
  // Only print error every 1000 iterations to avoid flooding console
  if (i % 1000 == 0) {
    printf("i: %d, Error: %g\n", i, error);
  }
  
  // Stop if we've reached steady state
  if (error < 1e-10) {
    printf("Converged at iteration %d with error %g\n", i, error);
    return 1; // Trigger end of simulation
  }
}

/**
 * @brief Saves the final simulation results to a CSV file.
 *
 * This event is triggered at the end of the simulation and writes the computed x-coordinates
 * and corresponding temperature values to "conduction-simple.csv". The output file provides data
 * that can be used to compare the numerical solution with the analytical solution.
 */
event end (t = end) {
  char filename[100];
  sprintf(filename, "conduction-simple.csv");

  FILE *fp = fopen(filename, "w");
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);
}
