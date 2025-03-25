/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

scalar f[], D2[];
vector u[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar D2c[], vel[];
scalar * list = NULL;

/**
   * @brief Main entry point for simulation data processing.
   *
   * This function parses command-line arguments to obtain simulation parameters,
   * including the snapshot filename, spatial boundaries (xmin, ymin, xmax, ymax),
   * and the number of grid points in the y-direction. It expects 7 arguments (the
   * program name plus 6 user-supplied values). If the argument count is incorrect, it
   * prints an error message with usage instructions and exits with an error code.
   *
   * Upon successful argument parsing, the function restores simulation data from the
   * specified file, computes derived scalar fields—such as a transformed derivative field
   * and the velocity magnitude—and interpolates these values over a structured grid.
   * The resulting data is then written to the output file.
   *
   * @return int Returns 1 if the command-line arguments are incorrect; otherwise, an exit
   *             status indicating the result of the simulation data processing.
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

  /*
  Actual run and codes!
  */
  restore (file = filename);

  #define UseD2c 1
  foreach() {
    #if UseD2c
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = sqrt(sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = f[]*(D2)/sqrt(2.);
    
    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }
    #else
    D2c[] = D2[];
    #endif

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
