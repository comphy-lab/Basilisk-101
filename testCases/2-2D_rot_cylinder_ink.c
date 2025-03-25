/**
 ## 2D Rotating Cylindrical Annulus with Dye Injection
 
 This simulation models flow in an annular region between two concentric cylinders,
 with the inner cylinder rotating at a constant angular velocity. A dye is injected
 to visualize the flow patterns. The simulation also includes a time reversal event
 to study mixing and unmixing phenomena.
 
 Contributed by [Paula Magrinya Aguilo](https://github.com/Pmagrinya)
*/

#include "embed.h"               // Embedded boundary module for complex geometries
#include "navier-stokes/centered.h" // Centered finite volume solver for Navier-Stokes
#include "die-injection.h"       // Module for dye injection into the flow field

// Define inner and outer radii of the annular region
#define INNER_RADIUS 1.0
#define OUTER_RADIUS 1.5

// Angular velocity of rotation (rad/s)
#define OMEGA 1.0
// Time interval between snapshots
#define tsnap 0.01
// Maximum simulation time
#define tmax 4.0
// Viscosity field (uniform in this case)
const face vector muv[] = {1e1, 1e1}; 

// File naming variables
char logFile[80];
char dumpFile[80];
char nameOut[80];

/**
 * @brief Configures and initiates the fluid dynamics simulation.
 *
 * This function sets up the computational domain based on the outer cylinder radius,
 * establishes grid resolution, and configures numerical solver parameters including
 * convergence tolerance, CFL condition, and timestep for implicit integration.
 * It also defines dye injection parameters, creates necessary directories for output,
 * and initializes file names for logging and restart data before starting the simulation
 * by invoking the run() function.
 *
 * @return int Exit status of the simulation application.
 */
int main() {
  // Set domain size based on outer cylinder radius
  L0 = 2.0 * OUTER_RADIUS;
  // Center the domain at origin for symmetric boundary conditions
  origin (-L0/2, -L0/2);
  // Initialize grid with resolution level 7 (2^7 = 128 cells per dimension)
  init_grid (1 << 8);
  
  // Numerical solver parameters
  TOLERANCE = 1e-5;        // Convergence tolerance for pressure solver
  CFL = 0.25;              // Courant-Friedrichs-Lewy condition for stability
  
  // Using implicit time stepping with fixed timestep for stability
  DT = 0.01;
  
  // Die injection parameters
  tInjection = 0.0;        // Inject the dye after flow is established
  xInjection = 0.00;        // X position of injection (center of cavity)
  yInjection = 1.25;        // Y position of injection
  dieRadius = 0.15;  // Radius of the injection region
  
  // Create directory for output files and set file names
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");
  sprintf (logFile, "logData.dat");
  
  //mu = muv;
  
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
    mu = fm;
    // Initialize velocity field (will be overwritten by boundary conditions)
    
    // Define the geometry: difference between outer and inner cylinders
    solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                             sq(INNER_RADIUS) - sq(x) - sq(y)));
    
    // Set normal and tangential velocity components at the boundary
    // For rotating inner cylinder: u = (-y, x) corresponds to counterclockwise rotation
    // Outer cylinder remains stationary
    u.n[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. : - y);
    u.t[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. :   x);
    
}

/**
 ## Output event for writing simulation data

 This event saves the full simulation state for restart capability
 and creates snapshots at regular intervals for visualization.
*/
event writingFiles (t = 0.0; t += tsnap; t < tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
 ## Logging event for monitoring simulation progress

 Records iteration number, simulation time, and timestep size
 both to standard error and to a log file for post-processing.
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

event adapt(i++) {
  adapt_wavelet ((scalar *){cs, u.x, u.y, T},
      (double[]){1e-3, 1e-2, 1e-2, 1e-3},
      10, 4);
}



/**
 ## Time reversal event

 At t=1, this event reverses the direction of rotation of the inner cylinder.
 This can be used to study flow reversibility and mixing properties.
 The velocity boundary conditions are flipped to (y, -x) for clockwise rotation.
*/
event timeReversal(t=0.5*tmax){
    u.n[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. : y);
    u.t[embed] = dirichlet (x*x + y*y > sq(OUTER_RADIUS-0.1) ? 0. :  -x);
}
