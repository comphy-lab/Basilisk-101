/* Title: getting Data from viscoelastic simulation snapshot
# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025
*/

#include "utils.h"
#include "output.h"

scalar f[];
vector u[];
scalar A11[], A12[], A22[]; // conformation tensor
scalar AThTh[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar D2c[], vel[], trA[];
scalar * list = NULL;

/**
   * @brief Entry point for processing viscoelastic simulation snapshot data.
   *
   * This function reads simulation parameters and a data file specified via command-line
   * arguments, computes derived quantities (such as shear rate and trace of the conformation
   * tensor) from the simulation data, and interpolates the results onto a uniform grid.
   *
   * The function expects exactly 7 command-line arguments:
   *   1. Program name.
   *   2. Input simulation data filename.
   *   3. Minimum x-coordinate.
   *   4. Minimum y-coordinate.
   *   5. Maximum x-coordinate.
   *   6. Maximum y-coordinate.
   *   7. Number of grid points in the y-direction.
   *
   * If the number of arguments is incorrect, an error message is printed to stderr and the
   * program exits with a status of 1.
   *
   * This function utilizes external global variables and helper routines (e.g., restore(),
   * list_add(), matrix_new()) and processes the simulation grid data using a foreach() macro.
   *
   * @return int Returns 0 upon successful execution, or 1 if the command-line arguments are invalid.
   */
  int main(int a, char const *arguments[])
{
  if (a != 7) {
    fprintf(stderr, "Error: Expected 6 arguments\n");
    fprintf(stderr, "Usage: %s <filename> <xmin> <ymin> <xmax> <ymax> <ny>\n", arguments[0]);
    return 1;
  }

  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);

  list = list_add (list, D2c);
  list = list_add (list, vel);
  list = list_add (list, trA);

  /*
  Actual run and codes!
  */
  restore (file = filename);

  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = f[]*D2;
    
    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }

    vel[] = sqrt(sq(u.x[])+sq(u.y[]));

    trA[] = f[]*((A11[] + A22[] + AThTh[])/3.0-1.0);

    // if (trA[] > 0.){
    //   trA[] = log(trA[])/log(10);
    // } else {
    //   trA[] = -10;
    // }

  }

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
