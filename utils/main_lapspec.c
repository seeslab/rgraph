#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "tools.h"
#include "graph.h"
#include "matrix.h"

int
main(int argc, char **argv)
{
  FILE *infile=NULL;
  struct node_gra *net=NULL;
  gsl_vector *spec;
  int i, nnod;
  char *netF;

  /* Command line parameters */
  if (argc < 2) {
    printf("\nUse: lapspec net_file\n\n");
    return -1;
  }
  netF = argv[1];

  /* Build the network */
  infile = fopen(netF, "r");
  net = FBuildNetwork(infile, 0, 0, 0, 1);
  fclose(infile);
  nnod = CountNodes(net);

  /* Calculate the spectrum of the Laplacian */
  spec = LaplacianSpectrum(net);

  /* Check results */
  for (i=0; i<nnod; i++)
    fprintf(stdout, "%d %lf\n", i+1, gsl_vector_get(spec, i));

  /* Free memory */
  RemoveGraph(net);
  gsl_vector_free(spec);

  /* Done */
  return 0;
}
