#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

/**
 ## 1D Transient Heat Conduction Solver
 
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

// Simulation parameters
#define N       1000     // Number of cells
#define L0      10.0    // Domain length (from -L0/2 to +L0/2)
#define X0      (-L0/2) // Left boundary

/**
 * @brief Creates the "intermediate" directory for storing output files.
 *
 * This function checks whether the "intermediate" directory exists. If it does not exist,
 * the function attempts to create it with owner permissions (read, write, and execute).
 * If directory creation fails, an error message is printed to stderr.
 *
 * @return 0 if the directory exists or is successfully created, 1 if directory creation fails.
 */
int create_output_directory() {
  struct stat st = {0};
  
  if (stat("intermediate", &st) == -1) {
    if (mkdir("intermediate", 0700) == -1) {
      fprintf(stderr, "Error creating intermediate directory: %s\n", strerror(errno));
      return 1;
    }
  }
  
  return 0;
}

/**
 * @brief Initializes the temperature field to approximate a Dirac delta function.
 *
 * This function initializes the temperature array using a thin rectangular pulse centered at x = 0. The pulse,
 * which has a width of 2 * eps, is scaled to ensure that its total integral equals 1. For each cell,
 * the cell-center coordinate is computed as X0 + (i + 0.5) * dx. If the absolute value of this coordinate is less than eps,
 * the cell is set to 1 / (2 * eps); otherwise, it is set to 0.
 *
 * @param temperature Pointer to the temperature array.
 * @param dx Spatial cell size.
 * @param eps Half-width of the rectangular pulse.
 */
void initialize_temperature(double *temperature, double dx, double eps) {
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx; // Cell-center coordinate
    if (fabs(x) < eps) {
      temperature[i] = 1.0 / (2.0 * eps);
    } else {
      temperature[i] = 0.0;
    }
  }
}

/**
 * @brief Prints the current temperature field to the console.
 *
 * This function iterates over the temperature array and prints each cell's center position,
 * temperature, and the current simulation time. The cell center is computed using the global offset
 * and the provided cell size.
 *
 * @param temperature Pointer to the temperature array.
 * @param dx The size of each cell.
 * @param time The current simulation time.
 */
void print_temperature(double *temperature, double dx, double time) {
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;
    printf("%g  %g  %g\n", x, temperature[i], time);
  }
  printf("\n\n");
}

/**
 * @brief Saves a snapshot of the temperature field to a CSV file.
 *
 * This function writes each cell's position and corresponding temperature value to a CSV file
 * in the "intermediate" directory. The filename is generated using the current simulation time.
 * If the file cannot be opened for writing, an error message is printed to stderr and the
 * function returns 1.
 *
 * @param temperature Array of temperature values for each cell.
 * @param dx Width of each cell (used to calculate cell positions).
 * @param time Current simulation time, used in the snapshot filename.
 * @return Returns 0 on success and 1 on failure.
 */
int save_snapshot(double *temperature, double dx, double time) {
  char filename[100];
  sprintf(filename, "intermediate/snapshot-%5.4f.csv", time);
  
  FILE *snap_file = fopen(filename, "w");
  if (snap_file == NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return 1;
  }
  
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;
    fprintf(snap_file, "%g,%g\n", x, temperature[i]);
  }
  
  fclose(snap_file);
  printf("Saved snapshot at time = %5.4f\n", time);
  
  return 0;
}

/**
 * @brief Computes heat fluxes at cell interfaces using a finite difference approximation.
 *
 * This function calculates the heat fluxes based on the temperature differences between adjacent cells.
 * The flux at each internal interface is computed as:
 *
 *     q[i+0.5] = - (T[i+1] - T[i]) / dx
 *
 * No-flux boundary conditions are enforced by setting the flux at both the left (index 0) and right (index N)
 * boundaries to zero.
 *
 * @param temperature Array containing temperature values at cell centers.
 * @param flux Array to store the computed fluxes, with a size of N+1.
 * @param dx The cell size.
 */
void compute_fluxes(double *temperature, double *flux, double dx) {
  flux[0] = 0.0; // Left boundary (no flux)
  flux[N] = 0.0; // Right boundary (no flux)
  
  for (int i = 0; i < N - 1; i++) {
    flux[i+1] = - (temperature[i+1] - temperature[i]) / dx;
  }
}

/**
 * @brief Updates the temperature array using flux differences.
 *
 * This function computes a finite-volume update for the temperature field by subtracting
 * the scaled flux differences from the current temperature. A no-flux condition is applied
 * at the left and right boundaries by assuming zero flux outside the domain.
 *
 * @param t_current Pointer to the current temperature array.
 * @param t_new Pointer to the array where the updated temperature will be stored.
 * @param flux Pointer to the array of fluxes at cell interfaces.
 * @param dx Cell size.
 * @param dt Time step.
 */
void update_temperature(double *t_current, double *t_new, double *flux, 
                        double dx, double dt) {
  for (int i = 0; i < N; i++) {
    double dq;
    
    if (i == 0) {
      // Leftmost cell: flux difference is q[0.5] - q[-0.5], but q[-0.5]=0
      dq = flux[1];
    } else if (i == N-1) {
      // Rightmost cell: flux difference is q[N-0.5] - q[N-1.5], but q[N+0.5]=0
      dq = - flux[N-1];
    } else {
      dq = flux[i+1] - flux[i];
    }
    
    t_new[i] = t_current[i] - (dt/dx) * dq;
  }
}

/**
 * @brief Copies updated temperature values into the current temperature array.
 *
 * This function iterates over the temperature arrays and copies each element from the new temperature array
 * (source) to the current temperature array (destination), updating the entire field.
 *
 * @param t_current Destination array to store the updated temperature values.
 * @param t_new Source array containing the new temperature values.
 */
void swap_temperature(double *t_current, double *t_new) {
  for (int i = 0; i < N; i++) {
    t_current[i] = t_new[i];
  }
}

/**
 * @brief Saves the final simulation temperature distribution to a CSV file.
 *
 * This function writes the cell-center coordinates and the corresponding temperature
 * values into "conduction-transient.csv". The x-coordinate for each cell is computed as
 * X0 + (i + 0.5) * dx.
 *
 * @param temperature Final temperature array.
 * @param dx Cell size used for computing cell-center coordinates.
 * @return int 0 if the file was successfully written, 1 if an error occurred while opening the file.
 */
int save_final_results(double *temperature, double dx) {
  FILE *file = fopen("conduction-transient.csv", "w");
  if (file == NULL) {
    fprintf(stderr, "Error opening final output file\n");
    return 1;
  }
  
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;  // Cell-center coordinate
    fprintf(file, "%g,%g\n", x, temperature[i]);
  }
  
  fclose(file);
  return 0;
}

/**
 * @brief Runs the transient heat conduction simulation.
 *
 * This function sets up and executes a one-dimensional transient heat conduction simulation using an explicit time integration scheme. It allocates memory for temperature and flux arrays, initializes the temperature field with an approximation of a Dirac delta function, and then iteratively computes fluxes, updates temperatures, and advances the simulation time. Periodic console outputs and CSV snapshots capture the simulation state, and the final temperature distribution is saved to a CSV file.
 *
 * @return 0 on success, non-zero on failure.
 */
int run_simulation() {
  // Numerical parameters
  const double dx   = L0/N;          // Cell size
  const double eps  = 0.1;           // Half-width of the initial "rectangle"
  
  // Time-stepping control (explicit)
  // Stability condition for 1D heat eq.: dt <= dx^2/2
  double dt        = 0.5 * dx * dx;  // Time step
  double tmax      = 1.0;            // Final time
  double tprint    = 0.1;            // Print interval for console output
  double tsnap     = 0.1;            // Snapshot interval for file output
  double nextPrint = 0.0;
  double nextSnap  = 0.0;

  // Create output directory
  if (create_output_directory() != 0) {
    return 1;
  }
  
  // Allocate arrays
  double *T  = (double*) calloc(N, sizeof(double)); // Temperature
  double *Tn = (double*) calloc(N, sizeof(double)); // Updated temperature
  double *q  = (double*) calloc(N+1, sizeof(double)); // Flux at interfaces
  
  if (T == NULL || Tn == NULL || q == NULL) {
    fprintf(stderr, "Error allocating memory\n");
    free(T);
    free(Tn);
    free(q);
    return 1;
  }

  // Initialize temperature field
  initialize_temperature(T, dx, eps);

  // Time integration
  double time = 0.0;
  while (time < tmax) {
    // Determine the time step for this iteration
    double current_dt = dt;
    
    // If next time step would cross a snapshot time, adjust dt to hit it exactly
    // !Highlight: this will be done automatically in Basilisk
    if (time < nextSnap && time + dt > nextSnap) {
      current_dt = nextSnap - time;
    }
    
    // Output data to console if print time is reached
    if (time >= nextPrint) {
      print_temperature(T, dx, time);
      nextPrint += tprint;
    }

    // Save snapshot if snapshot time is reached
    if (fabs(time - nextSnap) < 1e-10) {
      save_snapshot(T, dx, time);
      nextSnap += tsnap;
    }

    // Compute fluxes at cell interfaces
    compute_fluxes(T, q, dx);

    // Update temperature field
    update_temperature(T, Tn, q, dx, current_dt);

    // Swap arrays for next time step
    swap_temperature(T, Tn);

    // Advance time
    time += current_dt;
  }

  // Save final results
  int status = save_final_results(T, dx);

  // Free memory
  free(T);
  free(Tn);
  free(q);

  return status;
}

int main() {
  int status = run_simulation();
  return status;
}