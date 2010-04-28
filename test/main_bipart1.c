/*
  main_bipart1.c
  $LastChangedDate$
  $Revision$
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include "graph.h"
#include "bipartite.h"

int main()
{
  struct binet *binet = NULL;
  gsl_rng *randGen;

  /*
    --------------------------------------------------------------
    Initialize the random number generator
    --------------------------------------------------------------
  */
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, 3333);

  /*
    ------------------------------------------------------------
    Build the network
    ------------------------------------------------------------
  */
  fprintf(stderr, "Creating the network...\n");
  binet = BuildModularBipartiteNetwork(NULL, 32, 4, NULL, 64,
				       7, 7, -1.,
				       0.85, randGen);

  fprintf(stderr, "Inverting the network twice...\n");
  InvertBipart(binet);
  InvertBipart(binet);

  fprintf(stderr, "Randomizing the network...\n");
  RandomizeBipart(binet, 100, randGen);

  fprintf(stderr, "Projecting the network...\n");
  struct node_gra *projection;
  projection = ProjectBipart(binet);
  RemoveGraph(projection);

  FPrintPajekFileBipart ("binet.net", binet, 0, 0);

  /*
    ------------------------------------------------------------
    Network properties
    ------------------------------------------------------------
  */
  printf("L = %d\n", NLinksBipart(binet));
      
  /*
    ------------------------------------------------------------
    Free memory
    ------------------------------------------------------------
  */
  RemoveBipart(binet);
  gsl_rng_free(randGen);

  return 0;
}
