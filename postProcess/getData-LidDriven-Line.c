/* Title: Getting y data from Basilisk file
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

vector u[];
char filename[80];

int main(int a, char const *arguments[])
{
  if (a != 2) {
    fprintf(stderr, "Error: Expected 1 argument\n");
    fprintf(stderr, "Usage: %s <filename>\n", arguments[0]);
    return 1;
  }

  sprintf (filename, "%s", arguments[1]);
  
  restore (file = filename);

  FILE * fp = ferr;
  for (double y = -0.5; y < 0.5; y += 1e-2){
    fprintf(fp, "%g %g\n", y, interpolate(u.x, 0.0, y));
  }
  fflush (fp);
  fclose (fp);

}
