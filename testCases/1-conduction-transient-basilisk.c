/**
 ## 1D Transient Heat Conduction Solver (Basilisk version)
 
 This program solves the transient heat conduction equation in one dimension:
 
 $$
 \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 $$
 
 - Subject to no-flux boundary conditions on both ends of the domain.
 - Initial condition is a "Dirac delta" approximated by a thin rectangle centered at $x=0$ with total integral = 1.
 
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

/**
 * @brief Entry point for the 1D transient heat conduction simulation.
 *
 * This function initializes the simulation domain by setting the domain length, computing the left boundary,
 * defining the number of cells, and calculating the time step based on the CFL condition. It then creates an
 * output directory for intermediate files and starts the simulation by calling the run() function.
 *
 * Side effects:
 * - Creates an "intermediate" directory for storing simulation outputs.
 *
 * @return int The exit status of the program.
 */
int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
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
 
 Sets up a `Dirac delta` approximated by a thin rectangle centered at $x=0$ with total integral = 1.
*/
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 
}

/**
 ## Time integration using explicit finite volume method
*/
/**
 * @brief Advances the temperature field by one time step.
 *
 * This event performs a time integration step for the 1D transient heat conduction simulation using an explicit finite volume method.
 * It computes the fluxes at cell faces based on a central difference:
 *   $$ q_i = -\\frac{(T_i - T_{i-1})}{\\Delta} $$
 * and then evaluates the time derivative of temperature in each cell:
 *   $$ \\frac{dT}{dt} = -\\frac{(q_{i+1} - q_i)}{\\Delta} $$
 * The temperature field is updated via an explicit Euler scheme:
 *   $$ T^{n+1} = T^n + dt \\cdot \\frac{dT}{dt} $$
 * where the timestep, dt, is determined dynamically using dtnext(DT).
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Compute fluxes at the faces of the cells
  // q_i = -(T_i - T_{i-1})/Delta
  scalar q[];
  foreach()
    q[] = -(T[] - T[-1])/Delta;
  
  // Update temperature
  // dT/dt = -(q_{i+1} - q_i)/Delta
  // this enforces the central difference scheme
  scalar dT[];
  foreach()
    dT[] = -(q[1] - q[])/Delta;
  
  // Explicit Euler step
  foreach()
    T[] += dt*dT[];
  
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
 * @brief Saves the final temperature distribution to a CSV file.
 *
 * This event is triggered at the end of the simulation. It writes the spatial coordinate and the corresponding
 * temperature value to "conduction-transient.csv". Each CSV row contains a pair of values (x, T), allowing for 
 * comparison with the analytical solution.
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
