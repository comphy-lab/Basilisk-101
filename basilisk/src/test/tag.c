#include "utils.h"
#include "tag.h"

int main()
{
  size (1.[0]);
  init_grid (64);
  X0 = Y0 = -0.5;
  scalar t[];
  foreach()
    t[] = sin(4*pi*x)*cos(4.*pi*y) > 0.5*(x + 1.);
  int n = tag (t);
  foreach() {
    if (t[] == 0)
      t[] = nodata;
    else
      fprintf (qerr, "%g %g %g\n", x, y, t[]);
  }
  output_ppm (t, file = "t.png", checksum = stderr, min = 1, max = n, n = 256);
  
  foreach()
    t[] = sin(4*pi*x)*cos(4.*pi*y) > 0.5*(x + 1.);
  remove_droplets (t, -3);
  output_ppm (t, file = "t1.png", checksum = stderr, min = 0, max = 1, n = 256);
}
