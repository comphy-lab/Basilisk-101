/**
 ## 1D Transient Heat Conduction Solver (Basilisk version)
 
 This program solves the transient heat conduction equation in one dimension:
 
 $$
 \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 $$
 
 - Subject to no-flux boundary conditions on both ends of the domain.
 - Initial condition is a "Dirac delta" approximated by a thin rectangle
 centered at $x=0$ with total integral = 1.
 
 The exact self-similar analytical solution is:
 
 $$
 T(x,t) = \frac{1}{2\sqrt{\pi t}}e^{-\frac{x^2}{4t}}
 $$
*/

#include "grid/cartesian1D.h"
#include "run.h"

// Declare scalar field for temperature
scalar T[];

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1
double K;
/**
 * @brief Initializes simulation parameters and starts the transient heat conduction solver.
 *
 * This function sets up the computational domain by specifying the domain length, left boundary,
 * and the number of cells for discretization. It computes a stable timestep based on the CFL
 * condition for an explicit finite volume scheme and creates an output directory ("intermediate")
 * for saving simulation snapshots. Finally, it launches the simulation by calling the run() function.
 *
 * @note The stability constant K must be properly defined (currently marked as "XX") to ensure a stable solution.
 *
 * @return int Exit status.
 */
int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  N = 200;       // Number of cells
  
  // Set the timestep based on stability criterion (CFL condition)
  // dt = dx^2/2 for explicit scheme
  K = XX; // fill in the value of K for a stable solution
  DT = (L0/N)*(L0/N)/K;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field

 - Sets up a `Dirac delta` approximated by a thin rectangle
 - centered at $x=0$ with total integral = 1.
*/
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 

}

/**
 ## Time integration using explicit finite volume method
*/
/**
 * @brief Executes a single explicit Euler time integration step for the 1D heat conduction simulation.
 *
 * This event calculates the next time step using dtnext(DT) and updates the temperature field T.
 * It outlines the process for computing temperature fluxes at cell faces and the corresponding rate of 
 * temperature change, which are then used to update the temperature field via an explicit Euler step.
 *
 * @note The detailed implementations for flux calculation and temperature derivative computation are pending.
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // TODO: Implement the heat equation integration
  // HINT: You need to compute fluxes at cell faces and update temperature field
  
  // 1. Define a scalar field for temperature fluxes
  scalar q[];
  
  // 2. Compute fluxes at cell faces
  // q_i = -(T_i - T_{i-1})/Delta
  // YOUR CODE HERE
  
  // 3. Compute temperature change rate
  // dT/dt = -(q_{i+1} - q_i)/Delta
  scalar dT[];
  // YOUR CODE HERE
  
  // 4. Update temperature field using explicit Euler step
  // T = T + dt * dT
  // YOUR CODE HERE
}

/**
 ## Save snapshots at regular intervals
*/
event writingFiles (t += tsnap; t < tmax+tsnap) {
  char filename[100];
  sprintf(filename, "intermediate/snapshot-%5.4f.csv", t);  

  FILE *fp = fopen(filename, "w");  
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);

}

/**
 * @brief Saves the final simulation temperature distribution to a CSV file.
 *
 * At the conclusion of the simulation, this event writes the x-coordinate and corresponding
 * temperature values for each grid cell to "conduction-transient.csv". The resulting file can
 * be used for post-processing and comparison with analytical solutions.
 */
event end (t = end) {
  char filename[100];
  sprintf(filename, "conduction-transient.csv");

  FILE *fp = fopen(filename, "w");
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);
}
