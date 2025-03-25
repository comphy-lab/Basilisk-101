/* Title: Getting y data from Basilisk file along the x = 0 line
Compile this code using:
qcc -O2 -Wall getDatay.c -o getDatay -lm
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

vector u[];
char filename[80];
face vector shear[], D2p[];


/**
 * @brief Main entry point for extracting y data from a Basilisk simulation file.
 *
 * This function sets up boundary conditions (periodic on the right, slip at the top, and no slip at the bottom),
 * restores simulation data from a file specified by the first command-line argument, and computes shear and
 * deformation rates using finite difference approximations on the velocity field. The computed values are 
 * interpolated along the x = 0 line for y-coordinates in the range [-0.5, 0.5) with a step of 0.01, and the 
 * results are output to the standard error stream.
 *
 * The function assumes that the input file exists, is accessible, and that all restoration operations succeed.
 */
int main(int a, char const *arguments[])
{
  // Boundary condition: periodic right - left
  periodic (right);
  // Slip at the top
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  // No slip at the bottom
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
  sprintf (filename, "%s", arguments[1]);
  
  restore (file = filename);
  boundary((scalar *){u});
  
  FILE * fp = ferr;

  foreach_face(y){
    shear.y[] = (u.x[0,0] - u.x[0,-1])/(Delta);
  }

  foreach_face(y) {
    double D11 = (u.y[] - u.y[0,-1]);
    double D22 = ((u.x[1,0]-u.x[-1,0])+(u.x[1,-1]-u.x[-1,-1]))/4.0;
    double D12 = 0.5*(((u.y[1,0]-u.y[-1,0])+(u.y[1,-1]-u.y[-1,-1]))/4.0 + (u.x[] - u.x[0,-1]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    // D2 /= sqrt(2.0); // note that this is $\|D_{ij}\|$
    D2p.y[] = D2;
  }

  for (double y = -0.5; y < 0.5; y += 1e-2){
    fprintf(ferr, "%g %g %g %g\n", y, interpolate(u.x, 0.0, y),
            interpolate(shear.y, 0.0, y), interpolate(D2p.y, 0.0, y));
  }
  fflush (fp);
  fclose (fp);
}
