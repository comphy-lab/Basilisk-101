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
scalar cs[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar * list = NULL;

/**
 * @brief Main entry point for simulation data interpolation.
 *
 * This function validates the command-line arguments, parses input parameters,
 * restores simulation data, adjusts scalar fields with weighting factors, and performs interpolation
 * over a computed grid. It outputs the grid coordinates along with the interpolated values.
 *
 * The function expects exactly 6 additional arguments:
 *   <filename> <xmin> <ymin> <xmax> <ymax> <ny>
 * where:
 *   - <filename> specifies the input data file.
 *   - <xmin> and <ymin> define the lower bounds of the computational domain.
 *   - <xmax> and <ymax> define the upper bounds of the computational domain.
 *   - <ny> is the number of grid divisions along the y-axis.
 *
 * If the argument count is incorrect, it prints an error message with usage instructions to stderr
 * and exits with a status of 1.
 *
 * @param a Number of command-line arguments (including the program name).
 * @param arguments Array of command-line argument strings.
 *
 * @return int Returns 1 if the argument count is incorrect; otherwise, completes the data processing.
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

  /*
  Actual run and codes!
  */
  restore (file = filename);

  // if sum of all cs is 0, then set cs to 1
  double sum = 0;
  foreach() {
    sum += cs[];
  }
  if (sum == 0) {
    foreach() {
      cs[] = 1;
    }
  } else {
    for (scalar s in list){
      foreach(){
        s[] *= cs[];
      }
    }
  }

  list = list_add (list, cs);

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
