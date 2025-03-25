/**
## 2D Rotating Cylindrical Annulus Simulation
 
Simulates the flow between two concentric cylinders where the inner cylinder rotates.
This is a classic problem in fluid dynamics known as Taylor-Couette flow.

Contributed by [Paula Magrinya Aguilo](https://github.com/Pmagrinya)
*/

#include "embed.h"         // For embedded boundaries
#include "navier-stokes/centered.h"  // Centered Navier-Stokes solver

// Define inner and outer radii
#define INNER_RADIUS 1.0   // Radius of inner cylinder
#define OUTER_RADIUS 1.5   // Radius of outer cylinder

// Velocidad angular de rotación (rad/s)
#define OMEGA 1.0          // Angular velocity of rotation (rad/s)
#define tsnap 1.0          // Time interval between snapshots
#define tmax 10.0          // Maximum simulation time
const face vector muv[] = {1, 1}; // Viscosity coefficient (e.g., nu = 1e-3)

// File naming variables
char logFile[80];          // Log file name
char dumpFile[80];         // Restart file name
char nameOut[80];          /**
 * @brief Initializes the simulation environment and launches the simulation.
 *
 * Configures the simulation domain size based on the outer cylinder radius, centers the domain,
 * and sets a grid resolution of 128 cells per dimension. It also sets the maximum timestep size for
 * the implicit solver, creates the directory for simulation snapshots, assigns file names for dumping
 * restart data and logging, and then starts the simulation.
 *
 * @return int Exit status.
 */

int main() {
  // Ajustamos el tamaño del dominio para que sea [-L0/2,L0/2] x [-L0/2,L0/2] (por ejemplo)
  // Adjust the domain size to be [-L0/2,L0/2] x [-L0/2,L0/2] (for example)
  L0 = 2.0 * OUTER_RADIUS; // Domain size based on outer radius
  // Center the domain at origin
  origin (-L0/2, -L0/2);
  // Initialize grid (adjust resolution as needed)
  init_grid (1 << 7);      // Grid resolution 2^7 = 128 cells per dimension
  
  // We use the implicit solver for stability
  DT = 0.01;               // Maximum timestep size
  
  // Create a folder for simulation snapshots
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");
  sprintf (logFile, "logData.dat");
  
  // Iniciamos la simulación
  // Start the simulation
  run();
}

/**
 ## Initialization event (t=0)

 This event defines the geometry of the cylindrical annulus using the embedded boundary
 method and sets boundary conditions for the velocity field.
 The initial condition imposes rigid body rotation on the inner cylinder.
 */
event init (t = 0) {
    mu = fm;               // Set viscosity field for the fluid domain
    solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                             sq(INNER_RADIUS) - sq(x) - sq(y)));
    
    // Boundary conditions:
    // Normal velocity component (no-penetration)
    u.n[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. : - y);
    // Tangential velocity component (rotation)
    u.t[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. :   x);    
}

/**
 ## Creates simulation snapshots at regular time intervals
 
 This saves:
 1. A restart file for resuming the simulation
 2. Snapshot files in the 'intermediate' directory for visualization and analysis
 */
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
 * @brief Logs simulation statistics at each timestep.
 *
 * This event records the iteration number, simulation time, and timestep size by writing the data to both
 * standard error and a designated log file. On the first call (when i == 0), it outputs a header line ("i t dt")
 * to both destinations.
 */
event logWriting (i++) {
  if (i == 0) {
    fprintf (stderr, "i t dt\n");
    FILE *fp = fopen (logFile, "w");
    fprintf (fp, "i t dt\n");
    fclose (fp);
  }
  fprintf (stderr, "%d %g %g\n", i, t, dt);
  FILE *fp = fopen (logFile, "a");
  fprintf (fp, "%d %g %g\n", i, t, dt);
  fclose (fp);
}
