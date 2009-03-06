#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "tools.h"
#include "graph.h"
#include "matrix.h"

#define EPS 1.e-4

int
main()
{
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  gsl_vector *spec;
  int i;
  double rightResult[11];

  /* The expected results */
  rightResult[0] = 0.0000;
  rightResult[1] = 0.0000;
  rightResult[2] = 0.2043;
  rightResult[3] = 0.5405;
  rightResult[4] = 1.0000;
  rightResult[5] = 1.0000;
  rightResult[6] = 1.5989;
  rightResult[7] = 2.0000;
  rightResult[8] = 2.4425;
  rightResult[9] = 4.0170;
  rightResult[10] = 5.1969;

  /* Build the network */
  infile = fopen("test.dat", "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);

  /* Calculate the spectrum of the Laplacian */
  spec = LaplacianSpectrum(net);

  /* Check results */
  for (i=0; i<11; i++) {
    if (fabs(gsl_vector_get(spec, i) - rightResult[i]) < EPS)
      fprintf(stdout, "%lf\tOK\n", gsl_vector_get(spec, i));
    else
      return 1;
  }
  fprintf(stdout, "Synchronizability = %g\n", Synchronizability(net));

  /* Free memory */
  RemoveGraph(net);
  gsl_vector_free(spec);

  /* Done */
  return 0;
}
