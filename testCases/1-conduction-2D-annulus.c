/**
 
 # 2D Transient Heat Conduction Solver in an Annulus
 
 This program solves the transient heat conduction equation in an annular domain:
 
 $$
 \frac{\partial T}{\partial t} = \nabla^2 T
 $$
 
 ## Subject to boundary conditions:

  - T = 1 on the inner boundary (r = 1)
  - T = 0 on the outer boundary (r = 4)
 
 We use the embedded boundary method to define the annular domain.
*/

// Include necessary headers in the correct order for Basilisk
#include "run.h"
#include "embed.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];

// Define inner and outer radii
#define INNER_RADIUS 1.0
#define OUTER_RADIUS 4.0

// Simulation parameters
#define tmax 100.0
#define tsnap 1.0

mgstats mgd;
char nameOut[80], dumpFile[80], logFile[80];

/**
 * @brief Main entry point for the transient heat conduction simulation.
 *
 * Configures the computational domain and grid to fully enclose the annular geometry, sets the implicit
 * solver time step, and prepares the output folder and filenames for simulation snapshots and logs.
 * The simulation is started by invoking the run() function.
 *
 * The domain size is set to twice the outer radius and centered at the origin, ensuring the annulus is
 * completely contained within the simulation area.
 *
 * @return int Exit status code.
 */
int main() {
  // Domain setup
  // Make the domain large enough to contain the outer circle
  L0 = 2.0 * OUTER_RADIUS;  
  // Center the domain at origin
  origin (-L0/2, -L0/2);    
  // Initialize grid (adjust resolution as needed)
  init_grid (1 << 7);      
  
  // We use the implicit solver for stability
  DT = 0.01;
  
  // Create a folder for simulation snapshots
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");
  sprintf (logFile, "logData.dat");
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field and boundary conditions
*/
event init (t = 0) {
  if (!restore(file = dumpFile)) {
  /**
  Define the geometry of the embedded boundary as an annulus:
    - outer circle: sq(OUTER_RADIUS) - sq(x) - sq(y) (positive inside circle of radius 4)
    - inner circle: sq(INNER_RADIUS) - sq(x) - sq(y) (positive inside circle of radius 1)
    - difference gives the annular region (positive inside the annulus)
  */
  solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                           sq(INNER_RADIUS) - sq(x) - sq(y)));
    
  // Set boundary conditions using dirichlet()
  // T = 1 on inner boundary (r = 1)
  // T = 0 on outer boundary (r = 4)
  T[embed] = dirichlet (sq(x) + sq(y) < sq(INNER_RADIUS + 0.1) ? 1.0 : 0.0);
  }
}

/**
 ## Time integration using implicit diffusion solver
*/
/**
 * @brief Performs time integration by solving the diffusion equation for the temperature field.
 *
 * This event computes the next timestep using dtnext(DT), sets up a diffusion coefficient field
 * on cell faces (assigning a value of 1.0 in fluid cells, determined by fs.x[] > 0), and then
 * advances the temperature field by solving the diffusion equation using an implicit solver while
 * respecting embedded boundaries. The solver's multigrid data is stored in the global variable
 * `mgd`.
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Define diffusion coefficient field (constant = 1.0)
  face vector D[];
  foreach_face()
    D.x[] = fs.x[] > 0. ? 1.0 : 0.0; // Diffusion only in fluid cells
  
  // Solve the diffusion equation while respecting embedded boundaries
  mgd = diffusion (T, dt, D);
}

/**
 * @brief Adapts the computational grid and updates the embedded boundary conditions.
 *
 * This event performs wavelet-based adaptation on the temperature and cell-fraction fields using
 * error tolerances of 1e-4 and 1e-3, respectively, up to a refinement level of 9. It then redefines the
 * solid region for the annular domain by subtracting the inner circle from the outer circle, and sets a
 * Dirichlet condition on the temperature field where cells within a slightly enlarged inner radius are
 * fixed at 1.0, and others at 0.0.
 */
event adapt(i++){
  adapt_wavelet ((scalar *){T, cs},
    (double[]){1e-4, 1e-3},
    9);
  
  solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                           sq(INNER_RADIUS) - sq(x) - sq(y)));
  T[embed] = dirichlet (sq(x) + sq(y) < sq(INNER_RADIUS + 0.1) ? 1.0 : 0.0);
}

/**
 ## Save snapshots at regular intervals
*/
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
 * @brief Logs convergence details of the diffusion solver at each timestep.
 *
 * This event outputs the current iteration count (i), simulation time (t), timestep size (dt),
 * and the iteration number from the multigrid diffusion solver (mgd.i) to both standard error and
 * a designated log file. On the first timestep (i == 0), it writes the column headers to standard error
 * and resets the log file with the same header. On subsequent timesteps, it appends the log entry to the log file.
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