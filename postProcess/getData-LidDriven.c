/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025
*/

#include "utils.h"
#include "output.h"
#include "poisson.h"

vector u[];
char filename[1000];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar T[], vel[], psi[], omega[];
scalar * list = NULL;

/**
   * @brief Main entry point for processing fluid dynamics simulation data.
   *
   * This function validates the command-line arguments and initializes the simulation
   * by reading in the file name and domain parameters (xmin, ymin, xmax, ymax, ny). It registers
   * the necessary scalar fields (T, vel, psi), restores simulation data, applies Dirichlet
   * boundary conditions for both the velocity and streamfunction, computes the velocity magnitude
   * and vorticity, and then solves the Poisson equation to update the streamfunction. Finally, it
   * calculates grid spacing and interpolates the scalar fields over the designated grid.
   *
   * Command-line arguments:
   * - arguments[0]: Program name.
   * - arguments[1]: Filename of the simulation snapshot.
   * - arguments[2]: Lower bound (xmin) of the x-domain.
   * - arguments[3]: Lower bound (ymin) of the y-domain.
   * - arguments[4]: Upper bound (xmax) of the x-domain.
   * - arguments[5]: Upper bound (ymax) of the y-domain.
   * - arguments[6]: Number of grid points along the y-direction (ny).
   *
   * @param a Number of command-line arguments (expected to be 7).
   * @param arguments Array of command-line argument strings.
   *
   * @return Returns 1 if the argument validation fails; otherwise, the program proceeds with simulation processing.
   */
  int main(int a, char const *arguments[])
{
  if (a != 7) {
    fprintf(ferr, "Error: Expected 6 arguments\n");
    fprintf(ferr, "Usage: %s <filename> <xmin> <ymin> <xmax> <ymax> <ny>\n", arguments[0]);
    return 1;
  }

  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);

  list = list_add (list, T);
  list = list_add (list, vel);
  list = list_add (list, psi);

  /*
  Actual run and codes!
  */
 restore (file = filename);

  // Top moving wall
  u.t[top] = dirichlet(1);
  /**
  For the other no-slip boundaries this gives */
  u.t[bottom] = dirichlet(0);
  u.t[left]   = dirichlet(0);
  u.t[right]  = dirichlet(0);
  // solve for the streamfunction
  psi[top] = dirichlet(0);
  psi[bottom] = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);

  foreach() {
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));
  }

  foreach() {
    omega[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    psi[] = 0.;
  }
  poisson (psi, omega);

  FILE * fp = ferr;
  Deltay = (double)((ymax-ymin)/(ny));
  // fprintf(ferr, "%g\n", Deltay);
  nx = (int)((xmax - xmin)/Deltay);
  // fprintf(ferr, "%d\n", nx);
  Deltax = (double)((xmax-xmin)/(nx));
  // fprintf(ferr, "%g\n", Deltax);
  len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);

}
