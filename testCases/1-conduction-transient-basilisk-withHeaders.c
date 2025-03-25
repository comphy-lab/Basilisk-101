/**
 ## 1D Transient Heat Conduction Solver (Using Basilisk diffusion.h)
 
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
 
This version uses the diffusion.h module from Basilisk to handle the diffusion equation implicitly.
*/

// Include necessary headers in the correct order for Basilisk
#include "grid/multigrid1D.h"  /* Multigrid solver is required by diffusion.h */
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];
// Boundary conditions
// The diffusion solver will use homogeneous Neumann conditions by default
T[left] = neumann(0.);
T[right] = neumann(0.);

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1

/**
 * @brief Main simulation driver for the 1D transient heat conduction solver.
 *
 * Initializes the simulation domain and parameters—setting domain length (10.0), left boundary (-L0/2), number of cells (10000),
 * and a timestep (1e-4) suitable for the implicit solver. Creates the output directory for intermediate results and starts
 * the simulation by calling the run() function.
 *
 * @return int Exit status of the program (0 indicates successful execution).
 */
int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  N = 10000;       // Number of cells
  
  // We can use a larger timestep with the implicit solver
  // compared to the explicit version which requires dt = dx^2/2
  DT = 1e-4;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field
 
 Sets up a "Dirac delta" approximated by a thin rectangle
 centered at $x=0$ with total integral = 1.
*/
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 
}

/**
 ## Time integration using implicit diffusion solver
*/
/**
 * @brief Advances the temperature field by one time step using diffusion.
 *
 * This event computes the next time step with `dtnext(DT)` and then updates the temperature field `T`
 * by solving the diffusion (heat conduction) equation:
 * 
 *   $$\frac{\partial T}{\partial t} = \nabla^2 T$$
 * 
 * which corresponds to a unit diffusivity.
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  /**
   ## Solve the diffusion equation
   
   The heat equation is: 
   
   $$
   \frac{\partial T}{\partial t} = \nabla^2 T
   $$
   
   which corresponds to diffusion with D = 1
  */  
  diffusion(T, dt);
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
 * @brief Saves the final simulation results to a CSV file.
 *
 * At the end of the simulation, this event writes the spatial coordinates and corresponding
 * temperature values to "conduction-transient.csv". Each line in the file contains the position and
 * temperature, which can be used to compare the numerical results with the analytical solution.
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