/**
 ## Axisymmetric Transient Heat Conduction Solver
 
 This program solves the transient heat conduction equation in an axisymmetric domain:
 
 $$
 \frac{\partial T}{\partial t} = \nabla^2 T
 $$
 
 In cylindrical coordinates with axisymmetry, the Laplacian becomes:
 
 $$
 \nabla^2 T = \frac{1}{r}\frac{\partial}{\partial r}(r\frac{\partial T}{\partial r}) + \frac{\partial^2 T}{\partial z^2}
 $$
 
 ## Domain: 

  - Length: 10.0 units
  - Centered at x=0 (left boundary at -5.0, right boundary at 5.0)
  - Grid resolution: 128 cells (2^7)
 
 ## Subject to boundary conditions:
 
  - Top: T = 0 (Dirichlet)
  - Left (r=0): Gaussian heat flux profile exp(-y²) (Neumann)
  - Right: T = 0 (Dirichlet)
 
 ## Initial condition:
 
  - T = 0 everywhere
  
 The simulation uses an implicit diffusion solver, allowing for larger timesteps compared to explicit methods.
*/

// Include necessary headers in the correct order for Basilisk
#include "axi.h"
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];

// Simulation parameters
#define tmax 100.0
#define tsnap 1.0

mgstats mgd;
char nameOut[80], dumpFile[80], logFile[80];
/**
 * @brief Initializes and runs the transient heat conduction simulation.
 *
 * This function configures the simulation domain and grid for an axisymmetric transient heat
 * conduction problem using the Basilisk framework. It sets the domain length and grid resolution,
 * specifies an implicit time step, creates a directory for simulation snapshots, and defines file names
 * for restart and log outputs. The boundary conditions are established with Dirichlet conditions on
 * the top and right boundaries and a Neumann condition with a Gaussian flux profile on the left boundary.
 * Finally, the simulation is started by invoking the run() function.
 *
 * @return int Exit status of the program.
 */
int main() {
  // Domain setup
  L0 = 4.0;     // Domain length
  X0 = 0.0;    // Left boundary
  init_grid (1 << 7);      // Number of cells
  
  // We can use a larger timestep with the implicit solver
  // compared to the explicit version which requires dt = dx^2/2
  DT = 0.01;
  
  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");
  // Name of the log file. See logWriting event.
  sprintf (logFile, "logData.dat");
  
  // Boundary conditions
  T[top] = dirichlet(0.0);
  T[left] = neumann(4e0*exp(-y*y)); // Gaussian flux profile where y is r, decays away from r=0
  T[right] = dirichlet(0.);
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field
 
 Sets the initial temperature field to zero throughout the domain.
 The Basilisk event system will automatically call this function at t=0.
 */
event init (t = 0) {
  foreach()
    T[] = 0.0; 
}

/**
 ## Time integration using implicit diffusion solver
 
 This event performs the time integration for each timestep using
 the implicit diffusion solver from diffusion.h. The solver handles
 the axisymmetric formulation automatically through the axi.h module.
 
 The diffusion coefficient D is set to 1.0 across the entire domain,
 which corresponds to solving the standard heat equation:
 
 $$
 \frac{\partial T}{\partial t} = \nabla^2 T
 $$
 
 The solver uses multigrid acceleration for improved convergence. Convergence statistics are stored in the mgd variable for logging.
*/
/**
 * @brief Integrates the heat diffusion equation over a single time step.
 *
 * This function computes the next time step using the explicit helper dtnext(DT) and prepares
 * a constant diffusion coefficient field for solving the heat equation, ∂T/∂t = ∇²T, via the implicit
 * diffusion solver. The diffusion coefficient is set uniformly to 1.0 across all cell faces.
 *
 * The solver's convergence metric is stored in the global variable mgd.
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Use the diffusion() function from diffusion.h to solve the equation
  // The heat equation is: ∂T/∂t = ∇²T which corresponds to diffusion with D = 1
  face vector D[];
  foreach_face()
    D.x[] = fm.x[]; // Constant diffusion coefficient of 1.0
  
  mgd = diffusion(T, dt, D);
}

/**
 ## Save snapshots at regular intervals
 
 This event saves simulation data at regular time intervals (tsnap)
 throughout the simulation until the maximum time (tmax) is reached.
 
 Two types of files are generated:
 1. A single restart file named "restart" that gets overwritten at each interval
 2. Time-specific snapshots with timestamps stored in the "intermediate" directory
 
 These snapshot files contain the complete simulation state and can be used
 for visualization or to restart the simulation from a specific time.
*/
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/** 
 ## Write logs for monitoring simulation progress
 
 This event logs simulation metrics at every timestep to both:
 1. Standard error (console output)
 2. A log file named "logData.dat"
 
 For each timestep, it records:

 - i: The timestep number
 - t: The current simulation time
 - dt: The size of the current timestep
 - mgd.i: Number of iterations required by the diffusion solver
 
 This information is valuable for:
 
 - Monitoring simulation progress
 - Verifying solver convergence
 - Debugging numerical issues
 - Performance analysis
 */
event logWriting (i++) {

  if (i == 0) {
    fprintf (stderr, "i t dt mgd.i\n");
    FILE *fp = fopen (logFile, "w");
    fprintf (fp, "i t dt mgd.i\n");
    fclose (fp);
  }
  fprintf (stderr, "%d %g %g %d\n", i, t, dt, mgd.i);
  FILE *fp = fopen (logFile, "a");
  fprintf (fp, "%d %g %g %d\n", i, t, dt, mgd.i);
  fclose (fp);
}