/* Title: Getting y data from Basilisk file
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"

vector u[];
char filename[80];

/**
 * @brief Entry point for extracting and outputting interpolation data from a Basilisk file.
 *
 * This function processes a Basilisk file specified via command-line arguments. It first checks that exactly one argument 
 * (the filename) is provided; if not, it prints an error message with usage instructions to stderr and exits with code 1.
 * When a valid filename is provided, the function copies the filename to a global variable, restores the file data using 
 * the restore() function, and iterates over y-values in the range [-0.5, 0.5) in increments of 0.01. For each y-value, 
 * it computes an interpolated x-value using the interpolate() function and writes the pair (y, interpolated x) to the output stream.
 *
 * @param a Number of command-line arguments.
 * @param arguments Array of command-line arguments where the first element is the program name and the second is the filename.
 *
 * @return int Returns 1 if an incorrect number of arguments is provided; otherwise, returns 0 upon successful completion.
 */
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
