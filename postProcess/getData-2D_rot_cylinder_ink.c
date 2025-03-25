/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025
*/

#include "utils.h"
#include "output.h"

scalar T[];
vector u[];
scalar vel[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar * list = NULL;

/**
   * @brief Entry point for processing simulation snapshot data.
   *
   * This function validates and parses command-line arguments, restores simulation data from a snapshot
   * file, computes velocity magnitudes from the simulation's vector components, calculates grid spacing,
   * interpolates scalar fields over a generated grid based on the specified coordinate bounds and resolution,
   * and writes the interpolated results to an output file.
   *
   * If the number of arguments is not equal to 7, the function prints an error message along with usage
   * instructions to stderr and exits with status 1.
   *
   * @param a Number of command-line arguments.
   * @param arguments Array of command-line argument strings, where:
   *        - arguments[0]: Program name.
   *        - arguments[1]: Filename for the simulation snapshot.
   *        - arguments[2]: Minimum x-coordinate.
   *        - arguments[3]: Minimum y-coordinate.
   *        - arguments[4]: Maximum x-coordinate.
   *        - arguments[5]: Maximum y-coordinate.
   *        - arguments[6]: Number of divisions along the y-axis.
   *
   * @return int Exit status (0 on success, 1 if argument validation fails).
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

  list = list_add (list, T);
  list = list_add (list, vel);

  /*
  Actual run and codes!
  */
  restore (file = filename);

  foreach() {
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));
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
